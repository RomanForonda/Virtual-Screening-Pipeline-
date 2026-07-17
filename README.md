# Virtual Screening Pipeline

A Python pipeline for early-stage computer-aided drug discovery. It takes a compound dataset and a protein structure, prepares ligands and receptor files, runs AutoDock Vina docking, and writes a ranked results table.

## Quick Start

If you are new to this kind of workflow, follow these steps in order:

1. Install Conda or Miniconda.
2. Install OpenBabel, ADFRsuite, and AutoDock Vina.
3. Clone or download this repository.
4. Create the two Conda environments:

   ```bash
   conda env create -f environment-rdkit.yml
   conda env create -f environment-docking.yml
   ```

5. Put your input files in the `data/` folder:

   ```text
   data/Dataset.xlsx
   data/1M17.pdb
   ```

6. Open `configs/config.yaml` and check that the file names, Vina executable, and ligand residue name match your project.
7. Activate the main environment:

   ```bash
   conda activate rdkit_env
   ```

8. Run the pipeline:

   ```bash
   python main.py
   ```

9. Open the final results:

   ```text
   results/docking_results.xlsx
   ```

In most cases, you only need to edit `configs/config.yaml`. You should not need to edit the Python files unless you want to change how the pipeline works internally.

## What It Does

The pipeline runs these steps:

1. Curates the compound dataset and removes invalid or duplicate SMILES.
2. Applies drug-likeness filters, including Lipinski, Veber, PAINS, and Brenk.
3. Clusters molecules by chemical similarity and selects representative compounds.
4. Generates 3D ligand conformers.
5. Converts ligands to PDBQT with OpenBabel.
6. Cleans and repairs the protein structure.
7. Converts the receptor to PDBQT with ADFRsuite.
8. Runs molecular docking with AutoDock Vina.
9. Saves ranked docking results to Excel.

Final results are saved in:

```text
results/docking_results.xlsx
```

## Requirements

You need:

- Conda or Miniconda.
- OpenBabel.
- ADFRsuite.
- AutoDock Vina.

This project uses two Conda environments:

- `rdkit_env` for dataset processing, ligand preparation, ligand conversion, and docking.
- `docking_env` for protein cleaning and repair.

## Installation

Clone the repository:

```bash
git clone <your-repo-url>
cd Virtual_Screening_Pipeline
```

Create the Conda environments:

```bash
conda env create -f environment-rdkit.yml
conda env create -f environment-docking.yml
```

Install the external programs separately and make sure they work from the terminal:

```bash
obabel -V
prepare_receptor -h
vina_1.2.7_win.exe --help
```

On Windows, the Vina executable may have a versioned name such as `vina_1.2.7_win.exe`. If your executable has a different name, update it in `configs/config.yaml`.

## Input Files

Create a `data/` folder in the project root and place your input files there:

```text
data/Dataset.xlsx
data/1M17.pdb
```

The Excel dataset must contain at least these columns:

```text
id
smiles
```

The PDB file should contain the original co-crystallized ligand. The pipeline uses that ligand to define the docking box.

## Configuration

Before running a new project, edit:

```text
configs/config.yaml
```

At minimum, check these values:

```yaml
paths:
  dataset_file: Dataset.xlsx
  protein_file: 1M17.pdb

protein:
  ligand_resname: AQ4

tools:
  vina: vina_1.2.7_win.exe

docking:
  exhaustiveness: 8
  cpu: 8
  max_ligands: null
```

Important parameters:

- `dataset_file`: Excel file inside `data/`.
- `protein_file`: PDB file inside `data/`.
- `ligand_resname`: three-letter residue name of the co-crystallized ligand in the PDB.
- `vina`: exact name of your Vina executable.
- `exhaustiveness`: Vina search thoroughness. Use `8` for a fast first pass, `16` for more careful docking, or `32` for slower final/top-hit docking.
- `cpu`: number of CPU threads used by Vina.
- `max_ligands`: use `null` to dock all ligands, or an integer such as `10` for a short test run.

## Run The Pipeline

Activate the main environment:

```bash
conda activate rdkit_env
```

Run:

```bash
python main.py
```

The pipeline automatically reads:

```text
configs/config.yaml
```

## Outputs

Generated files are saved in `results/`.

Main outputs:

```text
results/docking_results.xlsx
results/pipeline_summary.txt
results/clusters_2d.png
results/pdbqt_ligands/
results/docking_poses/
```

The final table contains:

- `ligand`: original compound ID.
- `affinity_kcal_mol`: best predicted Vina binding affinity. More negative values are better predicted binders.

## Restarting Or Re-running

The pipeline skips steps whose output already exists in `results/`.

If you change important settings and want to regenerate outputs, delete the relevant files or folders in `results/` before running again.

For example, to force docking to run again, delete:

```text
results/docking_poses/
results/docking_results.xlsx
```

## Troubleshooting

- **Module not found**: make sure you activated `rdkit_env` before running `python main.py`.
- **Command not found**: check that OpenBabel, ADFRsuite, or Vina is installed and available on PATH.
- **Docking box missing**: check that `protein.ligand_resname` in `configs/config.yaml` matches the ligand residue name in the PDB file.
- **Vina executable not found**: check the `tools.vina` value in `configs/config.yaml`.

## License

This project is licensed under the MIT License.
