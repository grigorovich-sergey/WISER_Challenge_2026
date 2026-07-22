import matplotlib.pyplot as plt
import pandas as pd
import RNA

from model import (
    count_selection_violations,
    pairs_to_dot_bracket,
    repair_stem_selection,
    selected_stems_to_pairs,
)


def dot_bracket_to_pairs(structure):
    """Convert a nonpseudoknotted dot-bracket structure to base-pair indices."""
    stack = []
    pairs = []

    for index, character in enumerate(structure):
        if character == "(":
            stack.append(index)
        elif character == ")":
            if not stack:
                raise ValueError(
                    f"Unmatched closing parenthesis at index {index}."
                )
            pairs.append((stack.pop(), index))
        elif character != ".":
            raise ValueError(
                f"Unsupported character {character!r} at index {index}."
            )

    return sorted(pairs)


def base_pair_metrics(reference_structure, candidate_structure):
    """Calculate base-pair precision, recall, and F1."""
    reference_pairs = set(dot_bracket_to_pairs(reference_structure))
    candidate_pairs = set(dot_bracket_to_pairs(candidate_structure))

    if not reference_pairs and not candidate_pairs:
        precision = recall = f1 = 1.0
    else:
        true_positives = len(reference_pairs & candidate_pairs)
        precision_denominator = len(candidate_pairs)
        recall_denominator = len(reference_pairs)

        precision = (
            true_positives / precision_denominator
            if precision_denominator
            else 0.0
        )
        recall = (
            true_positives / recall_denominator
            if recall_denominator
            else 0.0
        )
        f1 = (
            2.0 * precision * recall / (precision + recall)
            if precision + recall
            else 0.0
        )

    return {
        "pair_precision": precision,
        "pair_recall": recall,
        "pair_f1": f1,
    }


def evaluate_candidate_structure(
    sequence,
    candidate_structure,
    reference_structure,
    reference_mfe,
):
    """Evaluate one candidate RNA structure against a reference structure."""
    fold_compound = RNA.fold_compound(sequence)
    candidate_energy = float(
        fold_compound.eval_structure(candidate_structure)
    )
    metrics = base_pair_metrics(reference_structure, candidate_structure)

    return {
        "candidate_energy": candidate_energy,
        "energy_gap": candidate_energy - float(reference_mfe),
        "pair_f1": metrics["pair_f1"],
    }


def evaluate_solver_sample(sample, sequence_row, metadata, repair=True):
    """Evaluate one binary solver sample and optionally repair conflicts."""
    sequence = str(sequence_row["sequence"])
    sequence_length = int(sequence_row["length"])
    stems = metadata["stems"]
    conflicts = metadata["conflicts"]
    x = [int(value) for value in sample["x"]]

    selected_before = [
        index for index, value in enumerate(x) if value == 1
    ]
    raw_violations = count_selection_violations(
        selected_before,
        stems,
        conflicts,
    )
    raw_overlap = int(raw_violations["overlap"])
    raw_crossing = int(raw_violations["crossing"])

    if repair:
        selected_after = sorted(
            repair_stem_selection(selected_before, stems)
        )
    else:
        selected_after = selected_before.copy()

    final_violations = count_selection_violations(
        selected_after,
        stems,
        conflicts,
    )
    final_valid = (
        int(final_violations["overlap"])
        + int(final_violations["crossing"])
        == 0
    )

    result = {
        "sequence_id": sequence_row["sequence_id"],
        "length": sequence_length,
        "variant": metadata["variant"],
        "sample_probability": sample.get("probability"),
        "bitstring": sample.get(
            "bitstring",
            "".join(str(value) for value in x),
        ),
        "raw_valid": raw_overlap + raw_crossing == 0,
        "raw_overlap_violations": raw_overlap,
        "raw_crossing_violations": raw_crossing,
        "stems_removed_by_repair": (
            len(selected_before) - len(selected_after)
        ),
        "candidate_structure": None,
        "candidate_energy": None,
        "energy_gap": None,
        "pair_f1": None,
    }

    if not final_valid:
        return result

    candidate_pairs = selected_stems_to_pairs(selected_after, stems)
    candidate_structure = pairs_to_dot_bracket(
        sequence_length,
        candidate_pairs,
    )
    candidate_metrics = evaluate_candidate_structure(
        sequence=sequence,
        candidate_structure=candidate_structure,
        reference_structure=sequence_row["reference_structure"],
        reference_mfe=sequence_row["reference_mfe"],
    )

    result["candidate_structure"] = candidate_structure
    result.update(candidate_metrics)
    return result


def evaluate_solver_result(
    solver_result,
    sequence_row,
    metadata,
    max_samples=20,
    repair=True,
):
    """Evaluate retained objective-ranked samples from one solver run."""
    samples = list(solver_result.get("samples", []))

    num_variables = int(
        metadata.get("num_variables", len(metadata["stems"]))
    )
    num_quadratic_terms = int(metadata.get("num_quadratic_terms", 0))
    rows = []

    for sample_rank, sample in enumerate(samples, start=1):
        row = evaluate_solver_sample(
            sample=sample,
            sequence_row=sequence_row,
            metadata=metadata,
            repair=repair,
        )
        row.update(
            {
                "sample_rank": sample_rank,
                "backend_mode": solver_result.get("backend_mode"),
                "backend_name": solver_result.get("backend_name"),
                "sampling_runtime_seconds": solver_result.get(
                    "runtime_seconds"
                ),
                "circuit_depth": solver_result.get("circuit_depth"),
                "num_qubits": solver_result.get("num_qubits"),
                "num_variables": num_variables,
                "num_quadratic_terms": num_quadratic_terms,
            }
        )
        rows.append(row)

    return pd.DataFrame(rows)


def build_per_run_summary(detailed_df):
    """Summarize retained candidates into one row per solver run."""
    group_columns = [
        "sequence_id",
        "length",
        "variant",
        "backend_mode",
        "backend_name",
        "run_id",
        "run_seed",
    ]

    ranked_df = detailed_df.copy()
    ranked_df["_has_candidate"] = ranked_df[
        "candidate_structure"
    ].notna()
    ranked_df = ranked_df.sort_values(
        by=[
            *group_columns,
            "_has_candidate",
            "energy_gap",
            "pair_f1",
            "sample_probability",
            "sample_rank",
        ],
        ascending=[
            *([True] * len(group_columns)),
            False,
            True,
            False,
            False,
            True,
        ],
        na_position="last",
        kind="stable",
    )

    best_df = (
        ranked_df.groupby(group_columns, dropna=False, sort=False)
        .head(1)
        .drop(columns="_has_candidate")
        .reset_index(drop=True)
    )

    behavior_rows = []

    for group_key, group_df in detailed_df.groupby(
        group_columns,
        dropna=False,
        sort=False,
    ):
        probabilities = pd.to_numeric(
            group_df["sample_probability"], errors="coerce"
        ).fillna(0.0)
        probability_mass = float(probabilities.sum())
        raw_valid = group_df["raw_valid"].fillna(False).astype(float)
        removed_stems = pd.to_numeric(
            group_df["stems_removed_by_repair"], errors="coerce"
        ).fillna(0.0)

        if probability_mass > 0:
            weighted_validity = float(
                (probabilities * raw_valid).sum() / probability_mass
            )
            weighted_removed_stems = float(
                (probabilities * removed_stems).sum() / probability_mass
            )
        else:
            weighted_validity = float(raw_valid.mean())
            weighted_removed_stems = float(removed_stems.mean())

        row = dict(zip(group_columns, group_key))
        row.update(
            {
                "retained_sample_count": len(group_df),
                "retained_probability_mass": probability_mass,
                "probability_weighted_raw_validity": weighted_validity,
                "probability_weighted_stems_removed": (
                    weighted_removed_stems
                ),
            }
        )
        behavior_rows.append(row)

    behavior_df = pd.DataFrame(behavior_rows)
    return best_df.merge(
        behavior_df,
        on=group_columns,
        how="left",
        validate="one_to_one",
    )


def aggregate_length_statistics(per_run_df):
    """Aggregate run-level results by length, variant, and backend."""
    group_columns = [
        "length",
        "variant",
        "backend_mode",
        "backend_name",
    ]

    summary_df = (
        per_run_df.groupby(group_columns, dropna=False, sort=False)
        .agg(
            num_sequences=("sequence_id", "nunique"),
            num_sequence_runs=("run_id", "size"),
            mean_energy_gap=("energy_gap", "mean"),
            std_energy_gap=("energy_gap", "std"),
            mean_pair_f1=("pair_f1", "mean"),
            std_pair_f1=("pair_f1", "std"),
            mean_weighted_raw_validity=(
                "probability_weighted_raw_validity",
                "mean",
            ),
            std_weighted_raw_validity=(
                "probability_weighted_raw_validity",
                "std",
            ),
            mean_weighted_stems_removed=(
                "probability_weighted_stems_removed",
                "mean",
            ),
            std_weighted_stems_removed=(
                "probability_weighted_stems_removed",
                "std",
            ),
            mean_retained_probability_mass=(
                "retained_probability_mass",
                "mean",
            ),
            std_retained_probability_mass=(
                "retained_probability_mass",
                "std",
            ),
            mean_num_variables=("num_variables", "mean"),
            std_num_variables=("num_variables", "std"),
            mean_num_quadratic_terms=(
                "num_quadratic_terms",
                "mean",
            ),
            std_num_quadratic_terms=(
                "num_quadratic_terms",
                "std",
            ),
            mean_num_qubits=("num_qubits", "mean"),
            std_num_qubits=("num_qubits", "std"),
            mean_circuit_depth=("circuit_depth", "mean"),
            std_circuit_depth=("circuit_depth", "std"),
            mean_sampling_seconds=("sampling_runtime_seconds", "mean"),
            std_sampling_seconds=("sampling_runtime_seconds", "std"),
        )
        .reset_index()
    )

    return summary_df.sort_values(
        ["backend_mode", "backend_name", "variant", "length"],
        kind="stable",
    ).reset_index(drop=True)


def plot_metric_by_length(
    summary_df,
    mean_column,
    std_column,
    y_label,
    title,
    backend_modes=("aer",),
    y_limits=None,
):
    """Plot one aggregated metric against sequence length."""
    if isinstance(backend_modes, str):
        backend_modes = (backend_modes,)
    else:
        backend_modes = tuple(backend_modes)

    plot_df = summary_df[
        summary_df["backend_mode"].isin(backend_modes)
    ].copy()

    backend_count = len(
        plot_df[["backend_mode", "backend_name"]].drop_duplicates()
    )
    fig, ax = plt.subplots()

    for (mode, backend_name, variant), variant_df in plot_df.groupby(
        ["backend_mode", "backend_name", "variant"],
        dropna=False,
        sort=False,
    ):
        variant_df = variant_df.sort_values("length", kind="stable")
        label = str(variant)

        if backend_count > 1:
            backend_label = (
                str(backend_name)
                if pd.notna(backend_name)
                else str(mode)
            )
            label = f"{variant} — {backend_label}"

        if std_column is None:
            x_values = variant_df["length"]
            y_values = variant_df[mean_column]
            
            ax.plot(
                x_values,
                y_values,
                marker="o",
                label=label,
            )
            continue
        
        else:

            ax.errorbar(
                variant_df["length"],
                variant_df[mean_column],
                yerr=pd.to_numeric(
                    variant_df[std_column], errors="coerce"
                ).fillna(0.0),
                marker="o",
                capsize=3,
                label=label,
            )

    ax.set_xlabel("Sequence length")
    ax.set_ylabel(y_label)
    ax.set_title(title)

    if y_limits is not None:
        ax.set_ylim(*y_limits)

    ax.legend()
    fig.tight_layout()
    return fig, ax


if __name__ == "__main__":
    pass
