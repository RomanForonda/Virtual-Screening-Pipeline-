# Fingerprint-based diversity selection with Tanimoto similarity and Butina clustering.

import os
from pathlib import Path

matplotlib_cache_dir = (
    Path(__file__).resolve().parent
    / ".matplotlib_cache"
)
matplotlib_cache_dir.mkdir(parents=True, exist_ok=True)
os.environ["MPLCONFIGDIR"] = str(matplotlib_cache_dir)

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.decomposition import PCA


def generate_fingerprints(mols):

    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2,   # Each atom considers its environment up to 2 bonds (ECFP4).
        fpSize=2048  # Common fingerprint size; 1024 and 2048 are typical choices.
    )

    return [
        generator.GetFingerprint(mol)
        for mol in mols
    ]


def calculate_neighbor_lists(fps, similarity_threshold):
    neighbor_lists = [
        [idx]
        for idx in range(len(fps))
    ]

    for i in range(1, len(fps)):
        similarities = DataStructs.BulkTanimotoSimilarity(
            fps[i],
            fps[:i]
        )

        for j, similarity in enumerate(similarities):
            if similarity >= similarity_threshold:
                neighbor_lists[i].append(j)
                neighbor_lists[j].append(i)

    return neighbor_lists


def cluster_from_neighbor_lists(neighbor_lists):
    sorted_indices = [
        (len(neighbors), idx)
        for idx, neighbors in enumerate(neighbor_lists)
    ]
    sorted_indices.sort(reverse=True)

    clusters = []
    seen = np.zeros(
        len(neighbor_lists),
        dtype=bool
    )

    while sorted_indices and sorted_indices[0][0] > 1:
        _, idx = sorted_indices.pop(0)

        if seen[idx]:
            continue

        cluster = [idx]
        seen[idx] = True

        for neighbor in neighbor_lists[idx]:
            if not seen[neighbor]:
                cluster.append(neighbor)
                seen[neighbor] = True

        clusters.append(
            tuple(cluster)
        )

    while sorted_indices:
        _, idx = sorted_indices.pop(0)

        if seen[idx]:
            continue

        clusters.append(
            (idx,)
        )

    return tuple(clusters)


def fingerprints_to_array(fps):

    fp_array = np.zeros(
        (len(fps), fps[0].GetNumBits()),
        dtype=np.uint8
    )

    for idx, fp in enumerate(fps):
        DataStructs.ConvertToNumpyArray(
            fp,
            fp_array[idx]
        )

    return fp_array


def calculate_cluster_representatives(clusters, fps):

    representatives = {}

    for cluster_id, cluster in enumerate(clusters, start=1):

        cluster_indices = list(cluster)

        if len(cluster_indices) == 1:
            representatives[cluster_id] = cluster_indices[0]
            continue

        mean_similarities = []

        for mol_idx in cluster_indices:
            other_indices = [
                other_idx
                for other_idx in cluster_indices
                if other_idx != mol_idx
            ]

            similarities = DataStructs.BulkTanimotoSimilarity(
                fps[mol_idx],
                [fps[other_idx] for other_idx in other_indices]
            )

            mean_similarities.append(
                np.mean(similarities)
            )

        best_position = int(
            np.argmax(mean_similarities)
        )

        representatives[cluster_id] = cluster_indices[best_position]

    return representatives


def add_2d_coordinates(df, fps):

    if len(fps) == 0:
        df["Cluster_X"] = []
        df["Cluster_Y"] = []
        return df

    if len(fps) == 1:
        df["Cluster_X"] = 0.0
        df["Cluster_Y"] = 0.0
        return df

    fp_array = fingerprints_to_array(fps)

    pca = PCA(
        n_components=2,
        random_state=42
    )

    coords = pca.fit_transform(
        fp_array
    )

    df["Cluster_X"] = coords[:, 0]
    df["Cluster_Y"] = coords[:, 1]

    return df


def save_cluster_plot(df, plot_path):

    if len(df) == 0:
        return

    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 8))

    scatter = plt.scatter(
        df["Cluster_X"],
        df["Cluster_Y"],
        c=df["Cluster_ID"],
        cmap="tab20",
        s=35,
        alpha=0.8,
        edgecolors="none"
    )

    selected = df[
        df["Selection_Reason"] == "large_cluster_representative"
    ]

    if len(selected) > 0:
        plt.scatter(
            selected["Cluster_X"],
            selected["Cluster_Y"],
            facecolors="none",
            edgecolors="black",
            s=90,
            linewidths=1.2,
            label="Large-cluster representative"
        )

    plt.xlabel("PCA 1")
    plt.ylabel("PCA 2")
    plt.title("Chemical diversity map by Butina cluster")
    plt.colorbar(
        scatter,
        label="Cluster ID"
    )
    if len(selected) > 0:
        plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(
        plot_path,
        dpi=300
    )
    plt.close()


def cluster_molecules(
    input_file,
    output_file,
    selected_output_file="Dataset_selected.xlsx",
    plot_file="clusters_2d.png",
    similarity_threshold=0.7,  # Common similarity threshold for ECFP4 clustering.
    max_cluster_size=10
):

    input_path = os.path.join(
        "results",
        input_file
    )

    output_path = os.path.join(
        "results",
        output_file
    )

    selected_output_path = os.path.join(
        "results",
        selected_output_file
    )

    plot_path = os.path.join(
        "results",
        plot_file
    )

    df = pd.read_excel(input_path)

    print(
        f"\nMolecules read: {len(df)}"
    )

    # -------------------------
    # RDKit molecules
    # -------------------------

    df["mol"] = df["smiles"].apply(
        Chem.MolFromSmiles
    )

    # -------------------------
    # Remove invalid SMILES.
    # This keeps the script robust if it is run directly without
    # first running data_curation.py.
    # -------------------------

    df = df[df["mol"].notna()].reset_index(drop=True)

    # -------------------------
    # Morgan fingerprints
    # ECFP4
    # -------------------------

    fps = generate_fingerprints(
        df["mol"]
    )

    print(
        "Fingerprints generated."
    )

    # -------------------------
    # Tanimoto neighbor lists
    # -------------------------

    neighbor_lists = calculate_neighbor_lists(
        fps,
        similarity_threshold
    )

    print(
        "Similarity neighbors calculated."
    )

    # -------------------------
    # Butina clustering groups molecules based on their neighbors.
    # RDKit's current Python implementation expands the triangular
    # distance list into a full n x n matrix, which is too memory-heavy
    # for larger screening libraries.
    # -------------------------

    clusters = cluster_from_neighbor_lists(
        neighbor_lists
    )

    print(
        f"Clusters found: "
        f"{len(clusters)}"
    )

    # -------------------------
    # Cluster_ID
    # -------------------------

    cluster_ids = [-1] * len(df)

    for cluster_id, cluster in enumerate(
        clusters,
        start=1
    ):
        for mol_idx in cluster:
            cluster_ids[mol_idx] = (
                cluster_id
            )

    df["Cluster_ID"] = cluster_ids

    # -------------------------
    # Cluster_Size
    # -------------------------

    cluster_sizes = {
        cluster_id: len(cluster)
        for cluster_id, cluster in enumerate(
            clusters,
            start=1
        )
    }

    df["Cluster_Size"] = (
        df["Cluster_ID"]
        .map(cluster_sizes)
    )

    # -------------------------
    # Representative molecules
    # -------------------------

    representatives = calculate_cluster_representatives(
        clusters,
        fps
    )

    representative_indices = set(
        representatives.values()
    )

    df["Representative"] = [
        mol_idx in representative_indices
        for mol_idx in range(len(df))
    ]

    # -------------------------
    # Select molecules for docking.
    # Small clusters keep all molecules. Large clusters keep only
    # the cluster medoid, i.e. the molecule most similar to the rest
    # of its cluster.
    # -------------------------

    df["Selected_for_Docking"] = (
        (df["Cluster_Size"] <= max_cluster_size)
        | df["Representative"]
    )

    df["Selection_Reason"] = np.where(
        df["Cluster_Size"] <= max_cluster_size,
        "small_cluster",
        np.where(
            df["Representative"],
            "large_cluster_representative",
            "large_cluster_not_selected"
        )
    )

    # -------------------------
    # 2D chemical diversity map
    # -------------------------

    df = add_2d_coordinates(
        df,
        fps
    )

    # -------------------------
    # Clean up
    # -------------------------

    df = df.drop(
        columns=["mol"]
    )

    # -------------------------
    # Save
    # -------------------------

    df.to_excel(
        output_path,
        index=False
    )

    selected_df = (
        df[df["Selected_for_Docking"]]
        .reset_index(drop=True)
    )

    selected_df.to_excel(
        selected_output_path,
        index=False
    )

    save_cluster_plot(
        df,
        plot_path
    )

    print(
        f"\nResults saved to:\n"
        f"{output_path}"
    )

    print(
        f"\nSelected molecules saved to:\n"
        f"{selected_output_path}"
    )

    print(
        f"\n2D cluster plot saved to:\n"
        f"{plot_path}"
    )

    print(
        f"\nMolecules selected for docking: "
        f"{len(selected_df)} / {len(df)}"
    )

    if len(clusters) > 0:
        print(
            f"Mean cluster size: "
            f"{len(df)/len(clusters):.2f}"
        )
    else:
        print("No clusters found.")

    return df
