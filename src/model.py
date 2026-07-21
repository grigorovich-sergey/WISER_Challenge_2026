from itertools import combinations

from qiskit_optimization import QuadraticProgram


def is_allowed_pair(base_a, base_b):
    """Return True when two nucleotides form an allowed RNA base pair."""
    allowed_pairs = {
        ("A", "U"),
        ("U", "A"),
        ("G", "C"),
        ("C", "G"),
        ("G", "U"),
        ("U", "G"),
    }
    
    return (base_a, base_b) in allowed_pairs


def enumerate_candidate_stems(sequence, min_stem_length=2, min_loop_length=3):
    """Enumerate consecutive antiparallel candidate stems in an RNA sequence."""

    stems = []
    seen_pair_sets = set()
    sequence_length = len(sequence)

    for left in range(sequence_length):
        for right in range(left + 1, sequence_length):
            pairs = []
            inner_left = left
            inner_right = right

            while (
                inner_right - inner_left - 1 >= min_loop_length
                and is_allowed_pair(sequence[inner_left], sequence[inner_right])
            ):
                pairs.append((inner_left, inner_right))
                inner_left += 1
                inner_right -= 1

            for stem_length in range(min_stem_length, len(pairs) + 1):
                stem_pairs = tuple(pairs[:stem_length])

                if stem_pairs in seen_pair_sets:
                    continue

                seen_pair_sets.add(stem_pairs)
                stems.append(
                    {
                        "id": len(stems),
                        "pairs": list(stem_pairs),
                        "length": stem_length,
                        "score": -float(stem_length),
                    }
                )

    return stems


def stems_share_nucleotide(stem_a, stem_b):
    """Return True when two stems use at least one common nucleotide."""
    
    positions_a = {index for pair in stem_a["pairs"] for index in pair}
    positions_b = {index for pair in stem_b["pairs"] for index in pair}
    
    return bool(positions_a & positions_b)


def stems_cross(stem_a, stem_b):
    """Return True when two stems contain a pseudoknot-style crossing."""
    for i, j in stem_a["pairs"]:
        for k, l in stem_b["pairs"]:
            if i < k < j < l or k < i < l < j:
                return True
    
    return False


def build_conflict_lists(stems):
    """Return overlap and crossing conflicts between candidate stems."""
    overlap_conflicts = []
    crossing_conflicts = []

    for stem_a, stem_b in combinations(stems, 2):
        conflict_pair = (stem_a["id"], stem_b["id"])

        if stems_share_nucleotide(stem_a, stem_b):
            overlap_conflicts.append(conflict_pair)
        if stems_cross(stem_a, stem_b):
            crossing_conflicts.append(conflict_pair)

    return {
        "overlap": overlap_conflicts,
        "crossing": crossing_conflicts,
    }


def selected_stems_to_pairs(selected_stem_ids, stems):
    """Combine base pairs from selected stems into a stable sorted list."""
    stems_by_id = {stem["id"]: stem for stem in stems}
    selected_pairs = set()

    for stem_id in selected_stem_ids:
        selected_pairs.update(stems_by_id[stem_id]["pairs"])

    return sorted(selected_pairs)


def pairs_to_dot_bracket(sequence_length, pairs):
    """Convert a valid noncrossing base-pair list to dot-bracket notation."""
    dot_bracket = ["."] * sequence_length
    used_indices = set()
    validated_pairs = []

    for left, right in pairs:
        if not (0 <= left < sequence_length and 0 <= right < sequence_length):
            raise ValueError("Pair index is outside the sequence range")
        if left >= right:
            raise ValueError("Each pair must satisfy left < right")
        if left in used_indices or right in used_indices:
            raise ValueError("A nucleotide cannot appear in more than one pair")

        for other_left, other_right in validated_pairs:
            if (
                left < other_left < right < other_right
                or other_left < left < other_right < right
            ):
                raise ValueError("Crossing base pairs are not supported")

        used_indices.add(left)
        used_indices.add(right)
        validated_pairs.append((left, right))
        dot_bracket[left] = "("
        dot_bracket[right] = ")"

    return "".join(dot_bracket)


def count_selection_violations(selected_stem_ids, stems, conflicts):
    """Count overlap and crossing conflicts in a selected stem set."""
    selected_ids = set(selected_stem_ids)

    overlap_count = sum(
        stem_a in selected_ids and stem_b in selected_ids
        for stem_a, stem_b in conflicts["overlap"]
    )
    crossing_count = sum(
        stem_a in selected_ids and stem_b in selected_ids
        for stem_a, stem_b in conflicts["crossing"]
    )

    return {
        "overlap": overlap_count,
        "crossing": crossing_count,
        "total": overlap_count + crossing_count,
    }


def repair_stem_selection(selected_stem_ids, stems):
    """Greedily repair a stem selection into a valid noncrossing subset."""
    stems_by_id = {stem["id"]: stem for stem in stems}
    unique_selected_ids = list(dict.fromkeys(selected_stem_ids))
    ranked_ids = sorted(
        unique_selected_ids,
        key=lambda stem_id: (
            stems_by_id[stem_id]["score"],
            -stems_by_id[stem_id]["length"],
            stem_id,
        ),
    )

    accepted_ids = []

    for stem_id in ranked_ids:
        candidate = stems_by_id[stem_id]
        is_valid = all(
            not stems_share_nucleotide(candidate, stems_by_id[accepted_id])
            and not stems_cross(candidate, stems_by_id[accepted_id])
            for accepted_id in accepted_ids
        )

        if is_valid:
            accepted_ids.append(stem_id)

    return accepted_ids


def build_qubo_model(
    stems,
    conflicts,
    variant="strict",
    overlap_penalty=None,
    crossing_penalty=None,
):
    """Build a strict, relaxed, or postprocessed RNA stem QUBO model."""
    quadratic_program = QuadraticProgram()
    variable_names = {}
    max_stem_benefit = max(
        (-float(stem["score"]) for stem in stems),
        default=0.0,
    )
    default_penalty = max_stem_benefit + 1.0

    if overlap_penalty is None:
        overlap_penalty = default_penalty
    if crossing_penalty is None:
        crossing_penalty = default_penalty

    for stem in stems:
        stem_id = stem["id"]
        variable_name = f"x_{stem_id}"
        quadratic_program.binary_var(name=variable_name)
        variable_names[stem_id] = variable_name

    linear = {variable_names[stem["id"]]: stem["score"] for stem in stems}
    quadratic = {}

    def add_quadratic_penalty(stem_id_a, stem_id_b, penalty):
        variable_a = variable_names[stem_id_a]
        variable_b = variable_names[stem_id_b]
        variable_pair = tuple(sorted((variable_a, variable_b)))
        quadratic[variable_pair] = quadratic.get(variable_pair, 0.0) + penalty

    if variant in {"strict", "relaxed"}:
        for stem_id_a, stem_id_b in conflicts["overlap"]:
            add_quadratic_penalty(stem_id_a, stem_id_b, overlap_penalty)

    if variant == "strict":
        for stem_id_a, stem_id_b in conflicts["crossing"]:
            add_quadratic_penalty(stem_id_a, stem_id_b, crossing_penalty)

    quadratic_program.minimize(linear=linear, quadratic=quadratic)

    metadata = {
        "variant": variant,
        "num_variables": len(stems),
        "num_quadratic_terms": sum(
            coefficient != 0 for coefficient in quadratic.values()
        ),
        "stems": stems,
        "conflicts": conflicts,
    }

    return quadratic_program, metadata


def build_sequence_qubo(
    sequence_row,
    variant="strict",
    min_stem_length=2,
    min_loop_length=3,
    overlap_penalty=None,
    crossing_penalty=None,
):
    """Build one RNA QUBO directly from a sequence table row."""
    sequence = str(sequence_row["sequence"])
    stems = enumerate_candidate_stems(
        sequence,
        min_stem_length=min_stem_length,
        min_loop_length=min_loop_length,
    )

    if not stems:
        sequence_id = sequence_row.get("sequence_id", "unknown_sequence")
        raise ValueError(f"No candidate stems found for {sequence_id}.")

    conflicts = build_conflict_lists(stems)

    return build_qubo_model(
        stems=stems,
        conflicts=conflicts,
        variant=variant,
        overlap_penalty=overlap_penalty,
        crossing_penalty=crossing_penalty,
    )


if __name__ == "__main__":
    pass
