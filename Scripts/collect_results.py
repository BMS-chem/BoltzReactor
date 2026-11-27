"""
collect_results.py

Collects the Boltz result files for one experiment into a single table.

Reads:  experiments/<experiment>/outputs/boltz_results_jobs/predictions/<id>/
            affinity_<id>.json              (if affinity was on)
            confidence_<id>_model_0.json
Writes: experiments/<experiment>/results.csv

The <id> is the job YAML's file name without '.yaml', i.e. <proteinID>__<ligandID>,
so the two IDs are split back out into their own columns.

A run submitted with --n_seeds N writes outputs/seed_1/ .. outputs/seed_N/ instead
of a single outputs/boltz_results_jobs/. This script writes one row per pair per seed, 
with the seed in a SEED column. The spread across seeds is the noise floor.

Example:
    python Scripts/collect_results.py --experiment <experiment>
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# 1. Read the command-line arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Collect the Boltz results into one table.")
parser.add_argument("--experiment", required=True, help="Name of the folder inside experiments/")
args = parser.parse_args()


# ---------------------------------------------------------------------------
# 2. Work out the folder paths
# ---------------------------------------------------------------------------
# This script lives in <project>/Scripts/, so the project folder is one level up.
project_dir = Path(__file__).resolve().parent.parent
experiment_dir = project_dir / "experiments" / args.experiment
outputs_dir = experiment_dir / "outputs"
# The table sits at the top of the experiment folder, beside reference/ and
# enrichment/, not inside outputs/ where Boltz writes its own files.
results_file = experiment_dir / "results.csv"

if not experiment_dir.exists():
    sys.exit(f"ERROR: experiment folder not found: {experiment_dir}")


# ---------------------------------------------------------------------------
# 3. Find the run folders
# ---------------------------------------------------------------------------
# An unseeded run writes straight into outputs/. A run submitted with --n_seeds N
# writes outputs/seed_1/ .. outputs/seed_N/, one complete set of predictions each.
# Each of those is read as its own run, with its seed carried into every row.
if not outputs_dir.exists():
    sys.exit(f"ERROR: no outputs folder: {outputs_dir}\n"
             "       Has the Boltz prediction step run yet?")

seed_dirs = sorted(d for d in outputs_dir.glob("seed_*") if d.is_dir())

# A seed folder has to be seed_<number>, because the number goes in the table.
for d in seed_dirs:
    suffix = d.name.split("_", 1)[1]
    if not suffix.isdigit():
        sys.exit(f"ERROR: '{d.name}' is not a seed folder name (expected seed_<number>): {d}")

if seed_dirs:
    # A top-level boltz_results_* beside the seed folders is left over from an
    # earlier unseeded run. Mixing the two would put rows from different settings
    # in one table without saying so.
    strays = sorted(d for d in outputs_dir.glob("boltz_results_*") if d.is_dir())
    if strays:
        found = "\n  ".join(str(d) for d in strays)
        sys.exit(f"ERROR: {outputs_dir} has both seed folders and an unseeded run:\n  {found}\n"
                 "       Delete the one you do not want and run this again.")
    runs = [(int(d.name.split("_", 1)[1]), d) for d in seed_dirs]
    runs.sort()
    print(f"Found {len(runs)} seed folder(s): " + ", ".join(str(s) for s, _ in runs))
else:
    runs = [(None, outputs_dir)]


def find_predictions_dir(run_dir):
    """The predictions folder inside one run folder.

    Boltz writes to <run_dir>/boltz_results_<name of the jobs folder>/predictions/,
    where <name> is normally 'jobs'. rglob searches for it rather than assuming the
    name, so this still works if the folder is ever renamed.
    """
    found = sorted(run_dir.rglob("predictions"))
    if not found:
        sys.exit(f"ERROR: no predictions folder found under {run_dir}\n"
                 "       Has the Boltz prediction step run yet?")
    if len(found) > 1:
        listed = "\n  ".join(str(d) for d in found)
        sys.exit(f"ERROR: more than one predictions folder found under {run_dir}:\n  {listed}\n"
                 "       Delete the outputs of the run you do not want and try again.")
    return found[0]


# ---------------------------------------------------------------------------
# 4. Read one row per prediction folder, per run
# ---------------------------------------------------------------------------
# Each pair has its own folder, named after the job YAML (<proteinID>__<ligandID>).
rows = []
n_no_confidence = 0
n_no_affinity = 0

for seed, run_dir in runs:

    predictions_dir = find_predictions_dir(run_dir)
    print(f"Reading results from {predictions_dir}")

    n_pairs_this_run = 0

    for pair_dir in sorted(predictions_dir.iterdir()):

        if not pair_dir.is_dir():
            continue

        pair_id = pair_dir.name
        n_pairs_this_run += 1

        # Split the ID back into the protein and the ligand. Both make_jobs.py and the
        # input checks forbid '__' inside an ID, so there is exactly one separator.
        if "__" in pair_id:
            protein_id, ligand_id = pair_id.split("__", 1)
        else:
            protein_id, ligand_id = pair_id, ""

        row = {"PROTEIN_ID": protein_id, "LIGAND_ID": ligand_id, "PAIR_ID": pair_id}
        if seed is not None:
            row["SEED"] = seed

        # 4a. Affinity. The file is missing when the run used --no_affinity, which is fine.
        affinity_file = pair_dir / f"affinity_{pair_id}.json"
        if affinity_file.exists():
            with open(affinity_file) as f:
                affinity = json.load(f)
            # affinity_pred_value is log10(IC50 in uM): lower means a tighter binder.
            row["AFFINITY_PRED_VALUE"] = affinity.get("affinity_pred_value")
            row["AFFINITY_PROBABILITY_BINARY"] = affinity.get("affinity_probability_binary")
        else:
            n_no_affinity += 1

        # 4b. Confidence, for model 0 (the best-ranked pose).
        confidence_file = pair_dir / f"confidence_{pair_id}_model_0.json"
        if confidence_file.exists():
            with open(confidence_file) as f:
                confidence = json.load(f)
            # Only the overall scores are taken. The per-chain scores (chains_ptm and
            # pair_chains_iptm) are nested, so they are left in the JSON files.
            for key in ["confidence_score", "ptm", "iptm", "ligand_iptm", "protein_iptm",
                        "complex_plddt", "complex_iplddt", "complex_pde", "complex_ipde"]:
                row[key.upper()] = confidence.get(key)
        else:
            # No confidence file means this pair did not finish, so it is worth knowing about.
            where = pair_id if seed is None else f"{pair_id} (seed {seed})"
            print(f"  WARNING: no confidence file for '{where}'")
            n_no_confidence += 1

        # 4c. Where the structure is, so the row can be traced back to a file.
        structure_file = pair_dir / f"{pair_id}_model_0.cif"
        if not structure_file.exists():
            structure_file = pair_dir / f"{pair_id}_model_0.pdb"
        row["STRUCTURE"] = str(structure_file.relative_to(project_dir)) if structure_file.exists() else ""

        rows.append(row)

    if n_pairs_this_run == 0:
        sys.exit(f"ERROR: {predictions_dir} has no prediction folders in it")

if not rows:
    sys.exit(f"ERROR: no prediction folders found under {outputs_dir}")


# ---------------------------------------------------------------------------
# 5. Write the table
# ---------------------------------------------------------------------------
results = pd.DataFrame(rows)

# With seeds, keep each pair's seeds together. Without seeds, sort the tightest 
# predicted binders first (a lower log10 IC50 is tighter); na_position="last" 
# keeps any pair without an affinity value at the bottom.
if "SEED" in results.columns:
    results = results.sort_values(["PROTEIN_ID", "LIGAND_ID", "SEED"])
elif "AFFINITY_PRED_VALUE" in results.columns:
    results = results.sort_values("AFFINITY_PRED_VALUE", na_position="last")
else:
    results = results.sort_values(["PROTEIN_ID", "LIGAND_ID"])

results.to_csv(results_file, index=False, encoding="utf-8")


# ---------------------------------------------------------------------------
# 6. Report
# ---------------------------------------------------------------------------
if "SEED" in results.columns:
    n_pairs = results["PAIR_ID"].nunique()
    print(f"Wrote {len(results)} rows to {results_file} "
          f"({n_pairs} pair(s) x {results['SEED'].nunique()} seed(s))")
    # A pair missing from one seed's folder would otherwise be invisible: the
    # table would just have fewer rows for it, and a spread computed over two
    # seeds instead of three reads exactly like one computed over three.
    incomplete = results.groupby("PAIR_ID")["SEED"].count()
    incomplete = incomplete[incomplete != results["SEED"].nunique()]
    if len(incomplete):
        print(f"WARNING: {len(incomplete)} pair(s) do not appear under every seed, "
              "so their spread is over fewer runs than the rest:")
        for pair_id, n in incomplete.items():
            print(f"  {pair_id}: {n} seed(s)")
else:
    print(f"Wrote {len(results)} rows to {results_file}")
if n_no_affinity:
    print(f"Note: {n_no_affinity} pair(s) had no affinity file (expected if the run used --no_affinity)")
if n_no_confidence:
    print(f"WARNING: {n_no_confidence} pair(s) had no confidence file: these predictions did not finish")
