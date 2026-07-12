import os
import glob
import subprocess
import pandas as pd


# -------------------------------------------------------------
# CONFIGURATION
# -------------------------------------------------------------

# Exact name of the Vina executable. Change this value if you
# update Vina or rename the executable.
VINA_EXE = "vina_1.2.7_win.exe"

RECEPTOR = os.path.join("results", "1M17_fixed.pdbqt")
LIGANDS_DIR = os.path.join("results", "pdbqt_ligands")
OUTPUT_DIR = os.path.join("results", "docking_poses")
BINDING_SITE_FILE = os.path.join("results", "binding_site.txt")
RESULTS_EXCEL = os.path.join("results", "docking_results.xlsx")

# How thoroughly Vina searches for each ligand.
# Higher values are slower and may improve docking quality.
# 8 is Vina's default and is a practical first-pass screening value.
EXHAUSTIVENESS = 8

# Number of CPU threads Vina can use for each ligand.
# 4 is conservative: faster than the default without saturating most PCs.
CPU = 8

# Maximum number of new ligands to dock in one run.
# Use an integer for short test runs, or None to dock every remaining ligand.
MAX_LIGANDS = None


def settings_from_config(config):
    paths = config["outputs"]
    tools = config["tools"]
    docking = config["docking"]

    return {
        "vina_exe": tools["vina"],
        "receptor": os.path.join("results", paths["receptor_pdbqt"]),
        "ligands_dir": os.path.join("results", paths["pdbqt_ligands_dir"]),
        "output_dir": os.path.join("results", paths["docking_poses_dir"]),
        "binding_site_file": os.path.join("results", paths["binding_site"]),
        "results_excel": os.path.join("results", paths["docking_results"]),
        "exhaustiveness": docking["exhaustiveness"],
        "cpu": docking["cpu"],
        "max_ligands": docking["max_ligands"],
    }


def default_settings():
    return {
        "vina_exe": VINA_EXE,
        "receptor": RECEPTOR,
        "ligands_dir": LIGANDS_DIR,
        "output_dir": OUTPUT_DIR,
        "binding_site_file": BINDING_SITE_FILE,
        "results_excel": RESULTS_EXCEL,
        "exhaustiveness": EXHAUSTIVENESS,
        "cpu": CPU,
        "max_ligands": MAX_LIGANDS,
    }


def read_docking_box(path):
    """
    Read the docking box center and size from binding_site.txt.
    """
    values = {}

    with open(path) as f:
        for line in f:
            if "=" not in line:
                continue
            key, value = line.split("=")
            values[key.strip()] = float(value.strip())

    center = (values["center_x"], values["center_y"], values["center_z"])

    # If the file does not include box size values, use a reasonable default.
    size = (
        values.get("size_x", 20.0),
        values.get("size_y", 20.0),
        values.get("size_z", 20.0),
    )

    return center, size


def dock_ligand(ligand_path, center, size, settings):
    """
    Run Vina for one ligand against the receptor.
    Return (ligand_name, best_affinity_kcal_mol).
    """

    name = os.path.splitext(os.path.basename(ligand_path))[0]
    output_path = os.path.join(settings["output_dir"], f"{name}_out.pdbqt")

    command = [
        settings["vina_exe"],
        "--receptor", settings["receptor"],
        "--ligand", ligand_path,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--exhaustiveness", str(settings["exhaustiveness"]),
        "--cpu", str(settings["cpu"]),
        "--out", output_path,
    ]

    result = subprocess.run(command, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"Docking failed for {name}:")
        print(result.stderr)
        return name, None

    # -------------------------------------------------------
    # Extract the best affinity, mode 1, from Vina stdout.
    # Vina's result table looks like this:
    #
    #    mode |   affinity | ...
    #    -----+------------+----
    #       1       -7.3     ...
    #       2       -6.9     ...
    # -------------------------------------------------------
    affinity = None
    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "1":
            try:
                affinity = float(parts[1])
            except ValueError:
                pass
            break

    return name, affinity


def run_docking(config=None):

    settings = settings_from_config(config) if config else default_settings()

    os.makedirs(settings["output_dir"], exist_ok=True)

    center, size = read_docking_box(settings["binding_site_file"])
    print(f"Docking box center: {center}")
    print(f"Docking box size: {size}")

    ligand_files = sorted(glob.glob(os.path.join(settings["ligands_dir"], "*.pdbqt")))
    print(f"\nLigands to dock: {len(ligand_files)}")

    if settings["max_ligands"] is not None:
        print(f"Docking limit for this run: {settings['max_ligands']} new ligands")

    results = []
    docked_this_run = 0

    for i, ligand_path in enumerate(ligand_files, start=1):

        if (
            settings["max_ligands"] is not None
            and docked_this_run >= settings["max_ligands"]
        ):
            print(
                f"\nReached max_ligands={settings['max_ligands']}; "
                "stopping this run."
            )
            break

        name = os.path.splitext(os.path.basename(ligand_path))[0]
        output_path = os.path.join(settings["output_dir"], f"{name}_out.pdbqt")

        # Per-ligand checkpoint: skip ligands that were already docked.
        if os.path.exists(output_path):
            print(f"[{i}/{len(ligand_files)}] {name}: already done, skipping.")
            continue

        name, affinity = dock_ligand(ligand_path, center, size, settings)
        print(f"[{i}/{len(ligand_files)}] {name}: {affinity} kcal/mol")

        results.append({"ligand": name, "affinity_kcal_mol": affinity})
        docked_this_run += 1

    # -------------------------------------------------------
    # Save or update the Excel summary.
    # -------------------------------------------------------
    if os.path.exists(settings["results_excel"]):
        previous_df = pd.read_excel(settings["results_excel"])
        new_df = pd.DataFrame(results)
        df = pd.concat([previous_df, new_df], ignore_index=True)
    else:
        df = pd.DataFrame(results)

    df = df.sort_values("affinity_kcal_mol").reset_index(drop=True)
    df.to_excel(settings["results_excel"], index=False)

    print(f"\nResults saved to: {settings['results_excel']}")

    return df


if __name__ == "__main__":
    import argparse

    from config import load_config

    cli = argparse.ArgumentParser()
    cli.add_argument("--config", default="configs/config.yaml")
    args = cli.parse_args()

    run_docking(load_config(args.config))
