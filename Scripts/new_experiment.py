"""
new_experiment.py

Sets up a new experiment folder with template input CSVs.

Makes:  experiments/<experiment>/inputs/proteins.csv   (header row only)
        experiments/<experiment>/inputs/ligands.csv    (header row only)
        experiments/<experiment>/jobs/  outputs/  logs/

The name may contain '/' to group experiments in a subfolder, to any depth. Any
missing parent folders are made. The same name is then what every other script
takes, including submit_boltz.sh.

Example:
    python Scripts/new_experiment.py --experiment my_screen
    python Scripts/new_experiment.py --experiment cdk4/005_apo
"""

import argparse
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# 1. Read the command-line arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Set up a new experiment folder.")
parser.add_argument("--experiment", required=True, help="Name for the new folder inside experiments/")
args = parser.parse_args()


# ---------------------------------------------------------------------------
# 2. Check the name is usable as a folder path
# ---------------------------------------------------------------------------
# '/' is allowed so that experiments can be grouped in subfolders. Everything that
# would cause a confusing error later is rejected here. In particular the filepaths
# have to stay inside /experiments/.

name = args.experiment.strip().replace("\\", "/")

if name == "":
    sys.exit("ERROR: the experiment name is empty")
if " " in name:
    sys.exit(f"ERROR: the experiment name must not contain a space: {name}")
if name.startswith("/") or (len(name) > 1 and name[1] == ":"):
    sys.exit(f"ERROR: the experiment name must be relative to experiments/, "
             f"not an absolute path: {name}")
if name.endswith("/"):
    sys.exit(f"ERROR: the experiment name must not end with '/': {name}")
if "//" in name:
    sys.exit(f"ERROR: the experiment name has an empty folder in it: {name}")
if any(part == ".." for part in name.split("/")):
    sys.exit(f"ERROR: the experiment name must stay inside experiments/, so it "
             f"cannot contain '..': {name}")
if any(part == "." for part in name.split("/")):
    sys.exit(f"ERROR: the experiment name must not contain '.' as a folder: {name}")


# ---------------------------------------------------------------------------
# 3. Work out the folder paths
# ---------------------------------------------------------------------------
# This script lives in <project>/Scripts/, so the project folder is one level up.
project_dir = Path(__file__).resolve().parent.parent
experiments_dir = project_dir / "experiments"
experiment_dir = experiments_dir / name
inputs_dir = experiment_dir / "inputs"

# The checks above work on the text of the name. This one works on the path it
# actually resolves to, which is what finally decides where files get written, and
# so also covers a symlinked parent that the text checks cannot see.
if not experiment_dir.resolve().is_relative_to(experiments_dir.resolve()):
    sys.exit(f"ERROR: '{name}' resolves to {experiment_dir.resolve()}, which is "
             f"outside {experiments_dir}")

if experiment_dir.exists():
    print(f"Note: {experiment_dir} already exists. Only the missing parts will be made.")


# ---------------------------------------------------------------------------
# 4. Make the folders
# ---------------------------------------------------------------------------
# msa_jobs/ and msa_outputs/ are left out on purpose: run_pipeline.sh makes them
# when it needs them and deletes them again at step 3, so an experiment that has
# not been run should not have them. exist_ok=True leaves an existing folder alone.
for folder in ["inputs", "jobs", "outputs", "logs"]:
    (experiment_dir / folder).mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# 5. Write the two template CSVs
# ---------------------------------------------------------------------------
# Just the header row. A file that already exists is left untouched, so this
# cannot wipe out input you have already typed in.
templates = {
    "proteins.csv": "ID,SEQUENCE,COFACTORS,COVALENT_SITE,MODIFICATIONS\n",
    "ligands.csv": "ID,SMILES,WARHEAD\n",
}

for file_name, header in templates.items():
    file_path = inputs_dir / file_name
    if file_path.exists():
        print(f"  kept:    {file_name} (already there, left as it is)")
    else:
        with open(file_path, "w", encoding="utf-8", newline="") as f:
            f.write(header)
        print(f"  created: {file_name}")


# ---------------------------------------------------------------------------
# 6. Report, and say what to do next
# ---------------------------------------------------------------------------
print(f"\nExperiment folder ready: {experiment_dir}")
print("\nNext:")
print(f"  1. Fill in {inputs_dir / 'proteins.csv'}")
print("     ID,SEQUENCE,COFACTORS,COVALENT_SITE,MODIFICATIONS  - one row per protein.")
print("     COVALENT_SITE may be left empty if this is not a covalent experiment.")
print("     SEQUENCE separates chains with ':'. COFACTORS holds CCD codes")
print("     separated by ';' (e.g. ZN;MG) and may be left empty.")
print("     MODIFICATIONS names non-standard residues such as phosphosites, as")
print("     <chain>:<residue letter><position>:<CCD code> separated by ';'")
print("     (e.g. A:S42:SEP;B:Y117:PTR). Position counts from 1 within that chain.")
print("     SEP/TPO/PTR are phospho-Ser/Thr/Tyr. It may be left empty.")
print(f"  2. Fill in {inputs_dir / 'ligands.csv'}")
print("     ID,SMILES,WARHEAD  - one row per ligand.")
print("     WARHEAD may be left empty if the ligand is not covalent.")
print("     IDs must be unique and must not contain spaces, slashes or '__'.")
print(f"  3. On the cluster, submit it:  ./Scripts/submit_boltz.sh {name}")
