from pdbfixer import PDBFixer
from openmm.app import PDBFile
from Bio.PDB import PDBParser, PDBIO
import os


def fix_protein(input_file, output_file, hydrogenation_ph=7.4):

    input_path = os.path.join("results", input_file)
    output_path = os.path.join("results", output_file)

    fixer = PDBFixer(filename=input_path)

    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=hydrogenation_ph)

    temp_path = os.path.join("results", "temp.pdb")

    with open(temp_path, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", temp_path)

    # -------------------------------------------------------
    # Rename each histidine (HIS) according to its actual
    # protonation state.
    # -------------------------------------------------------
    # PDBFixer already decides, residue by residue, whether each
    # histidine has a hydrogen on ND1 (HD1), on NE2 (HE2), or on
    # both atoms. Here we translate that decision into the residue
    # names expected by docking preparation tools: HID, HIE or HIP.
    # -------------------------------------------------------

    for model in structure:
        for chain in model:
            for residue in chain:

                if residue.get_resname() != "HIS":
                    continue

                atom_names = {atom.get_name() for atom in residue}
                has_hd1 = "HD1" in atom_names
                has_he2 = "HE2" in atom_names

                if has_hd1 and has_he2:
                    residue.resname = "HIP"
                elif has_hd1:
                    residue.resname = "HID"
                elif has_he2:
                    residue.resname = "HIE"
                else:
                    # Rare case with no titratable hydrogen detected.
                    # HIE is the most common form at physiological pH.
                    residue.resname = "HIE"

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path)

    os.remove(temp_path)

    print("Protein ready for docking:", output_path)

    return output_path


if __name__ == "__main__":
    import argparse

    cli = argparse.ArgumentParser()
    cli.add_argument("--input", required=True, help="Input cleaned PDB filename inside results/")
    cli.add_argument("--output", required=True, help="Output repaired PDB filename inside results/")
    cli.add_argument("--ph", type=float, default=7.4)
    args = cli.parse_args()

    fix_protein(args.input, args.output, hydrogenation_ph=args.ph)
