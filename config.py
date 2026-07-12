from copy import deepcopy
from pathlib import Path

try:
    import yaml
except ModuleNotFoundError as error:
    raise ModuleNotFoundError(
        "PyYAML is required to read config files. Activate rdkit_env or install "
        "it with: conda install -c conda-forge pyyaml"
    ) from error


DEFAULT_CONFIG = {
    "paths": {
        "dataset_file": "Dataset.xlsx",
        "protein_file": "1M17.pdb",
    },
    "outputs": {
        "curated_dataset": "Dataset_curated.xlsx",
        "filtered_dataset": "Dataset_filtered.xlsx",
        "clustered_dataset": "Dataset_clustered.xlsx",
        "selected_dataset": "Dataset_selected.xlsx",
        "cluster_plot": "clusters_2d.png",
        "prepared_ligands": "Ligands_prep.sdf",
        "pdbqt_ligands_dir": "pdbqt_ligands",
        "clean_protein": "1M17_clean.pdb",
        "fixed_protein": "1M17_fixed.pdb",
        "receptor_pdbqt": "1M17_fixed.pdbqt",
        "binding_site": "binding_site.txt",
        "docking_poses_dir": "docking_poses",
        "docking_results": "docking_results.xlsx",
        "summary": "pipeline_summary.txt",
    },
    "environments": {
        "docking_env_name": "docking_env",
    },
    "protein": {
        "ligand_resname": "AQ4",
        "padding_angstroms": 8.0,
        "hydrogenation_ph": 7.4,
    },
    "clustering": {
        "similarity_threshold": 0.7,
        "max_cluster_size": 10,
    },
    "tools": {
        "obabel": "obabel",
        "prepare_receptor": "prepare_receptor",
        "vina": "vina_1.2.7_win.exe",
    },
    "docking": {
        "exhaustiveness": 8,
        "cpu": 8,
        "max_ligands": None,
    },
}


def deep_update(base, updates):
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_update(base[key], value)
        else:
            base[key] = value
    return base


def load_config(config_path="configs/config.yaml"):
    config = deepcopy(DEFAULT_CONFIG)
    path = Path(config_path)

    if not path.exists():
        if config_path == "configs/config.yaml":
            return config
        raise FileNotFoundError(f"Config file not found: {path}")

    with path.open("r", encoding="utf-8") as f:
        user_config = yaml.safe_load(f) or {}

    return deep_update(config, user_config)
