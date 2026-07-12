import os
import re
import subprocess

from rdkit import Chem


def safe_filename(name):
    """
    Convert a compound ID into a safe filename.
    Keep letters, numbers, hyphens, underscores and dots.
    """
    name = str(name).strip()
    name = re.sub(r"[^A-Za-z0-9_.-]", "_", name)
    return name or "ligand"


def sdf_to_pdbqt(
    input_sdf,
    output_dir="results/pdbqt_ligands",
    obabel_path="obabel"
):
    """
    Convert a multi-molecule SDF into separate PDBQT files while
    preserving the original compound ID as the output filename.
    """
    input_path = os.path.join("results", input_sdf)
    temp_dir = os.path.join("results", "temp_sdf_ligands")

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)

    print("\n==============================")
    print("SDF -> PDBQT conversion")
    print("==============================")
    print(f"Input:  {input_path}")
    print(f"Output: {output_dir}")

    supplier = Chem.SDMolSupplier(input_path, removeHs=False)

    converted = 0
    failed = 0

    for i, mol in enumerate(supplier, start=1):
        if mol is None:
            failed += 1
            continue

        if mol.HasProp("compound_id"):
            compound_id = mol.GetProp("compound_id")
        elif mol.HasProp("_Name"):
            compound_id = mol.GetProp("_Name")
        else:
            compound_id = f"ligand_{i}"

        filename_id = safe_filename(compound_id)
        temp_sdf = os.path.join(temp_dir, f"{filename_id}.sdf")
        output_pdbqt = os.path.join(output_dir, f"{filename_id}.pdbqt")

        writer = Chem.SDWriter(temp_sdf)
        writer.write(mol)
        writer.close()

        command = [
            obabel_path,
            temp_sdf,
            "-O",
            output_pdbqt,
        ]

        result = subprocess.run(
            command,
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            print(f"\nError converting {compound_id}:")
            print(result.stderr)
            failed += 1
            continue

        converted += 1

    print("\nConversion complete")
    print(f"PDBQT files generated: {converted}")
    print(f"Failed: {failed}")
    print(f"Saved to: {output_dir}")

    return output_dir
