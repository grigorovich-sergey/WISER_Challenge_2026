import random
from pathlib import Path

import pandas as pd
import RNA


def load_beacon_sequences(path, max_sequences = None, max_length = None):
    """Load RNA sequences from a BEACON file """
    
    raw_df = pd.read_csv(path)

    df = raw_df[["sequence"]].copy()

    df = df.dropna(subset=["sequence"])
    df["sequence"] = (
        df["sequence"]
        .astype(str)
        .str.strip()
        .str.upper()
        .str.replace("T", "U", regex=False)
    )
    df = df[df["sequence"] != ""]

    valid_mask = df["sequence"].str.fullmatch(r"[AUCG]+")
    invalid_count = int((~valid_mask).sum())

    df = df[valid_mask].copy()
    df["length"] = df["sequence"].str.len()

    if max_length is not None:
        df = df[df["length"] <= max_length]

    if max_sequences is not None:
        df = df.head(max_sequences)

    df = df.reset_index(drop=True)

    df["sequence_id"] = [
        f"beacon_{index:06d}" for index in range(1, len(df) + 1)
    ]
    df["source"] = "beacon"

    if df.empty:
        print("WARNING! Sequences loaded: 0")
    else:
        print(f"Sequences loaded: {len(df)}")
        print(
            "Length statistics: "
            f"min={df['length'].min()}, "
            f"max={df['length'].max()}, "
            f"median={df['length'].median():g}"
        )

    return df[["sequence_id", "sequence", "length", "source"]]


def append_synthetic_sequences(df, len_list, n_seq, seed=666):
    """ Generate synthetic RNA sequences from concatenated source"""

    target_lengths = list(len_list)

    source_sequences = (
        df["sequence"]
        .dropna()
        .astype(str)
        .loc[lambda values: values.str.len() > 0]
        .tolist()
    )

    random_generator = random.Random(seed)
    synthetic_rows = []
    synthetic_index = 1

    for target_length in target_lengths:
        required_length = target_length * n_seq
        combined_parts = []
        combined_length = 0

        while combined_length < required_length:
            selected_sequence = random_generator.choice(source_sequences)
            combined_parts.append(selected_sequence)
            combined_length += len(selected_sequence)

        combined_sequence = "".join(combined_parts)[:required_length]

        for sequence_number in range(n_seq):
            start = sequence_number * target_length
            end = start + target_length
            synthetic_sequence = combined_sequence[start:end]

            synthetic_rows.append(
                {
                    "sequence_id": f"synth_{synthetic_index:06d}",
                    "sequence": synthetic_sequence,
                    "length": target_length,
                    "source": "synth",
                }
            )
            synthetic_index += 1

    synthetic_df = pd.DataFrame(
        synthetic_rows,
        columns=[
            "sequence_id",
            "sequence",
            "length",
            "source",
        ],
    )

    return pd.concat(
        [df.copy(), synthetic_df],
        ignore_index=True,
        sort=False,
    )


def add_vienna_references(df):
    """ Add ViennaRNA minimum-free-energy reference structures and energies """
    
    result = df.copy()

    structures = []
    mfe_values = []

    for sequence in result["sequence"]:
        structure, mfe = RNA.fold(sequence)
        structures.append(structure)
        mfe_values.append(float(mfe))

    result["reference_structure"] = structures
    result["reference_mfe"] = mfe_values

    return result


def save_processed_sequences(df, path):
    """ Save RNA sequences + data to a CSV file """
    
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(output_path, index=False)

    print(f"Sequences saved: {len(df)}")
    print(f"Output path: {output_path}")


if __name__ == "__main__":
    beacon_df = load_beacon_sequences(
        "./beacon.csv",
        max_sequences=200,
        max_length=None,
    )

    beacon_with_synthetic = append_synthetic_sequences(
        beacon_df,
        len_list=[15, 20, 25, 30, 35],
        n_seq=5,
        seed=666,
    )

    beacon_with_references = add_vienna_references(beacon_with_synthetic)

    save_processed_sequences(
        beacon_with_references,
        "./data/processed/beacon_sequences.csv",
    )
