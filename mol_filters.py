import os
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski

from rdkit.Chem.FilterCatalog import (
    FilterCatalog,
    FilterCatalogParams
)


# -------------------------
# PAINS catalog
# -------------------------

params_pains = FilterCatalogParams()
params_pains.AddCatalog(
    FilterCatalogParams.FilterCatalogs.PAINS
)

pains_catalog = FilterCatalog(params_pains)

# -------------------------
# Brenk catalog
# -------------------------

params_brenk = FilterCatalogParams()
params_brenk.AddCatalog(
    FilterCatalogParams.FilterCatalogs.BRENK
)

brenk_catalog = FilterCatalog(params_brenk)


def passes_pains(mol):
    return pains_catalog.GetFirstMatch(mol) is None


def passes_brenk(mol):
    return brenk_catalog.GetFirstMatch(mol) is None


def apply_filters(input_file, output_file):

    input_path = os.path.join("results", input_file)
    output_path = os.path.join("results", output_file)

    df = pd.read_excel(input_path)

    print(f"\nInput molecules: {len(df)}")

    # -------------------------
    # RDKit molecules
    # -------------------------

    df["mol"] = df["smiles"].apply(
        lambda s: Chem.MolFromSmiles(s)
    )

    # -------------------------
    # Remove invalid SMILES.
    # This keeps the script robust if it is run directly without
    # first running data_curation.py.
    # -------------------------

    df = df[df["mol"].notna()].reset_index(drop=True)

    # -------------------------
    # Molecular descriptors
    # -------------------------

    df["MW"] = df["mol"].apply(
        Descriptors.MolWt
    )

    df["LogP"] = df["mol"].apply(
        Descriptors.MolLogP
    )

    df["HBD"] = df["mol"].apply(
        Lipinski.NumHDonors
    )

    df["HBA"] = df["mol"].apply(
        Lipinski.NumHAcceptors
    )

    df["TPSA"] = df["mol"].apply(
        Descriptors.TPSA
    )

    df["RotB"] = df["mol"].apply(
        Lipinski.NumRotatableBonds
    )

    df["FractionCSP3"] = df["mol"].apply(
        Lipinski.FractionCSP3
    )

    df["RingCount"] = df["mol"].apply(
        Lipinski.RingCount
    )

    df["HeavyAtomCount"] = df["mol"].apply(
        Lipinski.HeavyAtomCount
    )

    # -------------------------
    # Lipinski
    # Allow one violation.
    # -------------------------

    lipinski_violations = (
        (df["MW"] > 500).astype(int)
        + (df["LogP"] > 5).astype(int)
        + (df["HBD"] > 5).astype(int)
        + (df["HBA"] > 10).astype(int)
    )

    df["Lipinski"] = (
        lipinski_violations <= 1
    )

    # -------------------------
    # Veber
    # -------------------------

    df["Veber"] = (
        (df["TPSA"] <= 140)
        & (df["RotB"] <= 10)
    )

    # -------------------------
    # PAINS
    # -------------------------

    df["PAINS"] = df["mol"].apply(
        passes_pains
    )

    # -------------------------
    # Brenk
    # -------------------------

    df["Brenk"] = df["mol"].apply(
        passes_brenk
    )

    # -------------------------
    # Statistics
    # -------------------------

    print("\nIndividual filter results:")

    print(
        f"Pass Lipinski: "
        f"{df['Lipinski'].sum()} / {len(df)}"
    )

    print(
        f"Pass Veber: "
        f"{df['Veber'].sum()} / {len(df)}"
    )

    print(
        f"PAINS-free: "
        f"{df['PAINS'].sum()} / {len(df)}"
    )

    print(
        f"Brenk-free: "
        f"{df['Brenk'].sum()} / {len(df)}"
    )

    # -------------------------
    # Final filter
    # -------------------------

    final_filter = (
        df["Lipinski"]
        & df["Veber"]
        & df["PAINS"]
        & df["Brenk"]
    )

    filtered_df = df[final_filter].copy()

    print(
        f"\nFinal molecules: "
        f"{len(filtered_df)} / {len(df)}"
    )

    # -------------------------
    # Save the original data plus calculated descriptors.
    # -------------------------

    columns_to_drop = [
        "mol"
    ]

    filtered_df = filtered_df.drop(
        columns=columns_to_drop
    )

    filtered_df = filtered_df.reset_index(drop=True)

    filtered_df.to_excel(
        output_path,
        index=False
    )

    print(
        f"\nFiltered data saved to:\n"
        f"{output_path}"
    )

    return filtered_df
