import os

import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


def choose_main_fragment(mol):
    """
    If a SMILES contains a salt, RDKit reads it as multiple fragments.
    Keep the main compound fragment, prioritizing organic fragments and
    then the highest heavy-atom count.
    """
    fragments = Chem.GetMolFrags(
        mol,
        asMols=True,
        sanitizeFrags=True
    )

    if not fragments:
        return None

    def fragment_score(fragment):
        has_carbon = any(atom.GetAtomicNum() == 6 for atom in fragment.GetAtoms())
        heavy_atoms = fragment.GetNumHeavyAtoms()
        return has_carbon, heavy_atoms

    return max(fragments, key=fragment_score)


def clean_smiles(smiles):
    """
    Convert an input SMILES into a curated version:
    1. remove invalid SMILES,
    2. keep the main compound when salts are present,
    3. neutralize charges when RDKit can do so,
    4. return a canonical SMILES.
    """
    if not isinstance(smiles, str):
        return None

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    mol = choose_main_fragment(mol)

    if mol is None:
        return None

    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)

    return Chem.MolToSmiles(mol, canonical=True)


def curate_data(input_file, output_file):
    input_path = os.path.join("data", input_file)
    output_path = os.path.join("results", output_file)

    os.makedirs("results", exist_ok=True)

    df = pd.read_excel(input_path)

    smiles_column = next(
        (
            column
            for column in df.columns
            if str(column).strip().lower() == "smiles"
        ),
        None
    )

    if smiles_column is None:
        raise ValueError(
            "The input file must contain a SMILES column, for example "
            "'smiles' or 'SMILES'."
        )

    if smiles_column != "smiles":
        df = df.rename(columns={smiles_column: "smiles"})

    id_column = next(
        (
            column
            for column in df.columns
            if str(column).strip().lower() == "id"
        ),
        None
    )

    if id_column is not None and id_column != "id":
        df = df.rename(columns={id_column: "id"})

    original_count = len(df)

    df["original_smiles"] = df["smiles"]
    df["smiles"] = df["smiles"].apply(clean_smiles)

    df = df[df["smiles"].notna()]
    valid_count = len(df)

    df = df.drop_duplicates(subset="smiles")
    unique_count = len(df)

    df = df.reset_index(drop=True)

    df.to_excel(output_path, index=False)

    print(f"Initial molecules: {original_count}")
    print(f"Molecules with valid SMILES: {valid_count}")
    print(f"Molecules after salt removal and deduplication: {unique_count}")
    print(f"Curated data saved to: {output_path}")

    return df
