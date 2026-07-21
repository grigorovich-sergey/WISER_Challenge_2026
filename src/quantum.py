from collections.abc import Mapping
from time import perf_counter

import numpy as np
from qiskit.circuit.library import QAOAAnsatz
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import SamplerV2 as AerSampler
from qiskit_ibm_runtime import QiskitRuntimeService
from qiskit_ibm_runtime import SamplerV2 as RuntimeSampler
from qiskit_optimization import QuadraticProgram
from qiskit_optimization.algorithms import MinimumEigenOptimizer
from qiskit_optimization.minimum_eigensolvers import NumPyMinimumEigensolver, QAOA
from qiskit_optimization.optimizers import COBYLA
from qiskit_optimization.utils import algorithm_globals


def solve_exact(quadratic_program: QuadraticProgram) -> dict:
    """Solve a small binary quadratic program exactly."""
    if quadratic_program.get_num_binary_vars() < 1:
        raise ValueError(
            "quadratic_program must contain at least one binary variable."
        )

    optimizer = MinimumEigenOptimizer(NumPyMinimumEigensolver())

    start_time = perf_counter()
    result = optimizer.solve(quadratic_program)
    runtime_seconds = perf_counter() - start_time

    return {
        "best_fval": float(result.fval),
        "runtime_seconds": runtime_seconds,
    }


def make_aer_sampler(
    shots: int = 1024,
    seed: int | None = 666,
    backend_options: dict | None = None,
) -> AerSampler:
    """Create a shot-based Aer SamplerV2."""
    simulator_options = {} if backend_options is None else dict(backend_options)

    return AerSampler(
        default_shots=shots,
        seed=seed,
        options={"backend_options": simulator_options},
    )


def solve_qaoa_aer(
    quadratic_program: QuadraticProgram,
    reps: int = 1,
    shots: int = 1024,
    maxiter: int = 50,
    seed: int | None = 666,
    initial_point: list[float] | np.ndarray | None = None,
) -> dict:
    """Optimize a binary quadratic program with QAOA on Aer."""
    if quadratic_program.get_num_binary_vars() < 1:
        raise ValueError(
            "quadratic_program must contain at least one binary variable."
        )

    if initial_point is not None:
        initial_point = np.asarray(initial_point, dtype=float).reshape(-1)
        expected_parameters = 2 * reps

        if initial_point.size != expected_parameters:
            raise ValueError(
                "initial_point must contain "
                f"{expected_parameters} values for reps={reps}."
            )

    algorithm_globals.random_seed = seed
    aer_backend = AerSimulator()

    pass_manager = generate_preset_pass_manager(
        optimization_level=2,
        backend=aer_backend,
        seed_transpiler=seed,
    )
    sampler = make_aer_sampler(shots=shots, seed=seed)
    qaoa = QAOA(
        sampler=sampler,
        optimizer=COBYLA(maxiter=maxiter),
        reps=reps,
        initial_point=initial_point,
        pass_manager=pass_manager,
    )
    optimizer = MinimumEigenOptimizer(qaoa)

    start_time = perf_counter()
    result = optimizer.solve(quadratic_program)
    runtime_seconds = perf_counter() - start_time

    qaoa_result = result.min_eigen_solver_result
    optimal_point = getattr(qaoa_result, "optimal_point", None)
    optimal_parameters = (
        None
        if optimal_point is None
        else np.asarray(optimal_point, dtype=float).reshape(-1).tolist()
    )
    best_fval = None if result.fval is None else float(result.fval)

    return {
        "optimal_parameters": optimal_parameters,
        "best_fval": best_fval,
        "optimizer_evaluations": getattr(qaoa_result, "cost_function_evals", None),
        "runtime_seconds": runtime_seconds,
    }


def get_ibm_backend(
    min_num_qubits: int,
    backend_name: str | None = None,
):
    """Select an operational IBM Quantum hardware backend."""
    service = QiskitRuntimeService()

    if backend_name is not None:
        backend = service.backend(backend_name.strip())
    else:
        backend = service.least_busy(
            min_num_qubits=min_num_qubits,
            operational=True,
            simulator=False,
        )

    status = backend.status()

    if not status.operational:
        raise ValueError(f"Backend '{backend.name}' is not currently operational.")

    if backend.num_qubits < min_num_qubits:
        raise ValueError(
            f"Backend '{backend.name}' has {backend.num_qubits} qubits, "
            f"but at least {min_num_qubits} are required."
        )

    print(
        f"Selected backend: {backend.name} | "
        f"qubits: {backend.num_qubits} | "
        f"pending jobs: {status.pending_jobs}"
    )

    return backend


def build_qaoa_circuit(
    quadratic_program: QuadraticProgram,
    reps: int = 1,
) -> QAOAAnsatz:
    """Build a reusable, unmeasured QAOA circuit for a binary QUBO."""
    cost_operator, _ = quadratic_program.to_ising()

    return QAOAAnsatz(
        cost_operator=cost_operator,
        reps=reps,
        flatten=True,
    )


def sample_qaoa_parameters(
    quadratic_program: QuadraticProgram,
    parameters: list[float] | np.ndarray,
    backend_mode: str = "aer",
    backend=None,
    reps: int = 1,
    shots: int = 2048,
    seed: int | None = 42,
) -> dict:
    """Sample fixed QAOA parameters on Aer or IBM hardware.

    For IBM jobs, ``runtime_seconds`` is QPU execution time reported by the
    job metrics and does not include queue time.
    """
    num_binary_vars = quadratic_program.get_num_binary_vars()


    circuit = build_qaoa_circuit(
        quadratic_program=quadratic_program,
        reps=reps,
    )
    parameter_values = np.asarray(parameters, dtype=float).reshape(-1)

    bound_circuit = circuit.assign_parameters(
        parameter_values.tolist(),
        inplace=False,
    )
    bound_circuit.measure_all()

    if backend_mode == "aer":
        execution_backend = AerSimulator() if backend is None else backend
        pass_manager = generate_preset_pass_manager(
            optimization_level=2,
            backend=execution_backend,
            seed_transpiler=seed,
        )
        transpiled_circuit = pass_manager.run(bound_circuit)
        sampler = AerSampler.from_backend(
            execution_backend,
            default_shots=shots,
            seed=seed,
        )

        start_time = perf_counter()
        job = sampler.run([transpiled_circuit], shots=shots)
        primitive_result = job.result()
        runtime_seconds = perf_counter() - start_time
        job_id = None

    else:

        status = backend.status()

        if not status.operational:
            raise ValueError(f"Backend '{backend.name}' is not operational.")

        if backend.num_qubits < num_binary_vars:
            raise ValueError(
                f"Backend '{backend.name}' has {backend.num_qubits} qubits, "
                f"but the problem requires {num_binary_vars}."
            )

        execution_backend = backend
        pass_manager = generate_preset_pass_manager(
            optimization_level=2,
            backend=execution_backend,
            seed_transpiler=seed,
        )
        transpiled_circuit = pass_manager.run(bound_circuit)
        sampler = RuntimeSampler(mode=execution_backend)
        job = sampler.run([transpiled_circuit], shots=shots)
        primitive_result = job.result()
        job_id = job.job_id()

        try:
            runtime_seconds = float(job.metrics()["usage"]["qpu_charge_time_seconds"])
        except (AttributeError, KeyError, TypeError, ValueError):
            runtime_seconds = None

    pub_result = primitive_result[0]
    counts = {
        str(bitstring): int(count)
        for bitstring, count in pub_result.data.meas.get_counts().items()
    }
    samples = counts_to_candidate_samples(
        counts=counts,
        quadratic_program=quadratic_program,
    )

    backend_name = getattr(
        execution_backend,
        "name",
        type(execution_backend).__name__,
    )
    if callable(backend_name):
        backend_name = backend_name()

    return {
        "backend_mode": backend_mode,
        "backend_name": str(backend_name),
        "samples": samples,
        "runtime_seconds": runtime_seconds,
        "circuit_depth": transpiled_circuit.depth(),
        "num_qubits": num_binary_vars,
        "job_id": job_id,
    }


def counts_to_candidate_samples(
    counts: Mapping[str, int],
    quadratic_program: QuadraticProgram,
) -> list[dict]:
    """Convert measured bitstring counts into evaluated QUBO candidates."""
    num_binary_vars = quadratic_program.get_num_binary_vars()
    cleaned_counts: dict[str, int] = {}

    for raw_bitstring, count in counts.items():
        bitstring = "".join(raw_bitstring.split())

        cleaned_counts[bitstring] = cleaned_counts.get(bitstring, 0) + int(count)

    total_shots = sum(cleaned_counts.values())
    candidates = []

    for bitstring, count in cleaned_counts.items():
        x = [int(bit) for bit in reversed(bitstring)]
        objective = float(quadratic_program.objective.evaluate(x))
        candidates.append(
            {
                "x": x,
                "bitstring": bitstring,
                "count": count,
                "probability": count / total_shots,
                "objective": objective,
            }
        )

    objective_sense = quadratic_program.objective.sense.name.upper()

    if objective_sense == "MINIMIZE":
        candidates.sort(
            key=lambda candidate: (
                candidate["objective"],
                -candidate["count"],
                candidate["bitstring"],
            )
        )
    elif objective_sense == "MAXIMIZE":
        candidates.sort(
            key=lambda candidate: (
                -candidate["objective"],
                -candidate["count"],
                candidate["bitstring"],
            )
        )
    else:
        raise ValueError(f"Unsupported objective sense: {objective_sense!r}.")

    return candidates


if __name__ == "__main__":
    pass
