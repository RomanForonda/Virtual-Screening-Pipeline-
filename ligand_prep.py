import os
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem


def prepare_ligands(input_file, output_file):
    input_path = os.path.join("results", input_file)
    output_path = os.path.join("results", output_file)

    df = pd.read_excel(input_path)

    if "id" not in df.columns:
        raise ValueError(
            "The input file must contain a column named 'id' so the original "
            "compound identifier can be preserved."
        )

    print(f"\nMolecules read: {len(df)}")

    writer = Chem.SDWriter(output_path)

    successful = 0
    failed = 0

    for _, row in df.iterrows():
        smiles = row["smiles"]
        compound_id = str(row["id"])

        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            failed += 1
            continue

        try:
            mol = Chem.AddHs(mol)

            params = AllChem.ETKDGv3()
            params.randomSeed = 42

            embed_result = AllChem.EmbedMolecule(mol, params)
            if embed_result != 0:
                failed += 1
                continue

            AllChem.MMFFOptimizeMolecule(mol)

            # _Name is the main name field that many SDF converters read.
            # compound_id is an explicit copy of the original compound ID.
            mol.SetProp("_Name", compound_id)
            mol.SetProp("compound_id", compound_id)

            writer.write(mol)
            successful += 1

        except Exception:
            failed += 1

    writer.close()

    print(f"\nLigands prepared: {successful}")
    print(f"Failed: {failed}")
    print(f"\nFile generated:\n{output_path}")
