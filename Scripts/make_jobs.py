"""
make_jobs.py

Makes one Boltz YAML file for every protein x ligand pair in an experiment.

Reads:  experiments/<experiment>/inputs/proteins.csv   (ID, SEQUENCE, COFACTORS[, COVALENT_SITE][, MODIFICATIONS])
        experiments/<experiment>/inputs/ligands.csv    (ID, SMILES[, WARHEAD])
        experiments/<experiment>/inputs/covalent_prepared.csv   (if this is a covalent experiment)
Writes: experiments/<experiment>/jobs/<proteinID>__<ligandID>.yaml
        experiments/<experiment>/jobs/<proteinID>.yaml          (apo experiment, see below)

If covalent_prepared.csv is present, every job gets a bond constraint joining that atom to the protein's COVALENT_SITE. 

An experiment with no ligands - either a ligands.csv with its header and no rows, or
no ligands.csv at all - is an apo experiment: one job per protein, holding only the
protein chains and any cofactors. There is no ligand chain, so there is no affinity
prediction and no covalent bond constraint.

Example:
    python Scripts/make_jobs.py --experiment <experiment>
    python Scripts/make_jobs.py --experiment <experiment> --no_affinity
"""

import argparse
import hashlib
import sys
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem


# ---------------------------------------------------------------------------
# 0. Defining modified residues
# ---------------------------------------------------------------------------
# 
# MODIFICATIONS optionally names non-standard residues (e.g. phospho PTMs) as a
# ';'-separated list of <chain>:<residue letter><position>:<CCD code>, for example
# "A:S42:SEP" or "A:Y117:PTR;B:T5:TPO".
# Position is counted from 1 within that chain, not across the whole SEQUENCE.
# 
# This list is used to parse that modifications have been written correctly before
# committing to processing the experiment.

MODIFICATION_PARENTS = {
    "SEP": "S",   # phosphoserine
    "TPO": "T",   # phosphothreonine
    "PTR": "Y",   # O-phosphotyrosine
    "NEP": "H",   # N1-phosphonohistidine
    "ALY": "K",   # N(6)-acetyllysine
    "MLZ": "K",   # N-methyl-lysine
    "MLY": "K",   # N-dimethyl-lysine
    "M3L": "K",   # N-trimethyllysine
    "KCX": "K",   # lysine NZ-carboxylic acid
    "HYP": "P",   # 4-hydroxyproline
    "CSO": "C",   # S-hydroxycysteine
    "CSD": "C",   # 3-sulfinoalanine
    "OCS": "C",   # cysteinesulfonic acid
    "CME": "C",   # S,S-(2-hydroxyethyl)thiocysteine
    "HIC": "H",   # 4-methyl-histidine
    "PCA": "Q",   # pyroglutamic acid
    "MSE": "M",   # selenomethionine (Boltz maps this back to MET, so it does nothing)
}

# Boltz's parser rewrites MSE to MET before building the residue, so asking for it
# produces an ordinary methionine. Worth saying rather than letting it look applied.
MODIFICATIONS_WITHOUT_EFFECT = {"MSE"}

chain_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


# ---------------------------------------------------------------------------
# 1. Read the command-line arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Make one Boltz YAML per protein x ligand pair.")
parser.add_argument("--experiment", required=True, help="Name of the folder inside experiments/")
parser.add_argument("--no_affinity", action="store_true", help="Leave out the affinity prediction")
args = parser.parse_args()


# ---------------------------------------------------------------------------
# 2. Work out the folder paths
# ---------------------------------------------------------------------------
# This script lives in <project>/Scripts/, so the project folder is one level up.
# Using this (rather than the current folder) means the script works wherever you run it from.
project_dir = Path(__file__).resolve().parent.parent
experiment_dir = project_dir / "experiments" / args.experiment
inputs_dir = experiment_dir / "inputs"
jobs_dir = experiment_dir / "jobs"

if not experiment_dir.exists():
    sys.exit(f"ERROR: experiment folder not found: {experiment_dir}")


# ---------------------------------------------------------------------------
# 3. Read the two CSV files
# ---------------------------------------------------------------------------
# encoding="utf-8-sig" removes the invisible marker some text editors put at the start of a file.
# dtype=str and keep_default_na=False keep every cell as text, with empty cells as "" (not NaN).
proteins = pd.read_csv(inputs_dir / "proteins.csv", encoding="utf-8-sig", dtype=str, keep_default_na=False)

# An apo experiment has no ligands at all, written either as a ligands.csv with its
# header and no rows, or as no ligands.csv at all. Both are read as an empty table.
ligands_file = inputs_dir / "ligands.csv"
if ligands_file.exists():
    ligands = pd.read_csv(ligands_file, encoding="utf-8-sig", dtype=str, keep_default_na=False)
else:
    ligands = pd.DataFrame({"ID": pd.Series(dtype=str), "SMILES": pd.Series(dtype=str)})

# Remove stray spaces around column names and values
proteins.columns = proteins.columns.str.strip()
ligands.columns = ligands.columns.str.strip()
proteins = proteins.apply(lambda column: column.str.strip())
if len(ligands) > 0:
    ligands = ligands.apply(lambda column: column.str.strip())

print(f"Read {len(proteins)} proteins and {len(ligands)} ligands")

# No ligands means an apo experiment: predict each protein on its own.
apo = len(ligands) == 0
if apo:
    if ligands_file.exists():
        print("ligands.csv has no rows: apo experiment, one job per protein.")
    else:
        print("No ligands.csv: apo experiment, one job per protein.")

# If prepare_covalents.py has run and found warheads, it leaves this file behind.
# Its presence is what makes this a covalent experiment.
prepared_file = inputs_dir / "covalent_prepared.csv"
covalent = prepared_file.exists()
prepared_by_id = {}

if covalent:
    prepared = pd.read_csv(prepared_file, encoding="utf-8-sig", dtype=str, keep_default_na=False)
    prepared.columns = prepared.columns.str.strip()
    prepared = prepared.apply(lambda column: column.str.strip())
    prepared_by_id = {row["ID"]: row for _, row in prepared.iterrows()}
    print(f"Covalent experiment: read {len(prepared_by_id)} prepared ligand(s)")


# ---------------------------------------------------------------------------
# 4. Check the inputs
# ---------------------------------------------------------------------------
# Every problem found is added to this list, and the script stops at the end
# of this section if there are any, so you can fix them all in one go.
problems = []

# 4a. The expected columns are there
for column in ["ID", "SEQUENCE", "COFACTORS"]:
    if column not in proteins.columns:
        problems.append(f"proteins.csv has no '{column}' column")
for column in ["ID", "SMILES"]:
    if column not in ligands.columns:
        problems.append(f"ligands.csv has no '{column}' column")
if problems:
    sys.exit("ERROR:\n  " + "\n  ".join(problems))

# 4b. IDs are filled in, unique, and safe to use in a file name
# ("__" is not allowed because it separates the protein and ligand IDs in the job name)
for file_name, table in [("proteins.csv", proteins), ("ligands.csv", ligands)]:
    for id_ in table["ID"]:
        if id_ == "":
            problems.append(f"{file_name}: a row has an empty ID")
        elif " " in id_ or "/" in id_ or "\\" in id_ or "__" in id_:
            problems.append(f"{file_name}: ID '{id_}' contains a space, a slash or '__'")
    duplicated = table["ID"][table["ID"].duplicated()].unique()
    for id_ in duplicated:
        problems.append(f"{file_name}: ID '{id_}' appears more than once")

# 4c. Protein sequences only contain the 20 standard amino acid letters
amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
for _, row in proteins.iterrows():
    if row["SEQUENCE"] == "":
        problems.append(f"proteins.csv: '{row['ID']}' has an empty SEQUENCE")
        continue
    for chain in row["SEQUENCE"].split(":"):
        bad_letters = set(chain.upper()) - amino_acids
        if chain == "":
            problems.append(f"proteins.csv: '{row['ID']}' has an empty chain (check the ':' separators)")
        elif bad_letters:
            problems.append(f"proteins.csv: '{row['ID']}' has letters that are not amino acids: {sorted(bad_letters)}")

# 4d. SMILES can be read by RDKit
for _, row in ligands.iterrows():
    if row["SMILES"] == "":
        problems.append(f"ligands.csv: '{row['ID']}' has an empty SMILES")
    elif Chem.MolFromSmiles(row["SMILES"]) is None:
        problems.append(f"ligands.csv: '{row['ID']}' has a SMILES that RDKit cannot read: {row['SMILES']}")

# 4e. Prevent make_jobs.py from being run on a stale or incomplete covalent_prepared.csv.
if covalent:
    for ligand_id in ligands["ID"]:
        if ligand_id not in prepared_by_id:
            problems.append(f"'{ligand_id}' is in ligands.csv but not in covalent_prepared.csv. "
                            "Re-run prepare_covalents.py.")
    if "COVALENT_SITE" not in proteins.columns:
        problems.append("this is a covalent experiment, but proteins.csv has no "
                        "'COVALENT_SITE' column")
    else:
        for _, row in proteins.iterrows():
            if row["COVALENT_SITE"] == "":
                problems.append(f"proteins.csv: '{row['ID']}' has an empty COVALENT_SITE")

# 4f. Read and check the MODIFICATIONS column.
# Each protein's entry is turned into {chain index: [{"position": int, "ccd": str}, ...]},
# which is what section 6 puts into the YAML.
modifications_by_protein = {}
unknown_codes = set()
inert_codes = set()

for _, row in proteins.iterrows():
    protein_id = row["ID"]
    modifications_by_protein[protein_id] = {}

    if "MODIFICATIONS" not in proteins.columns or row["MODIFICATIONS"] == "":
        continue

    # A protein whose SEQUENCE is empty or malformed has already been reported by 4c.
    # Checking positions against it as well would only add confusing extra errors.
    if row["SEQUENCE"] == "":
        continue

    chains = [chain.upper() for chain in row["SEQUENCE"].split(":")]
    seen_positions = set()

    for entry in row["MODIFICATIONS"].split(";"):
        entry = entry.strip()
        if entry == "":
            continue

        parts = entry.split(":")
        if len(parts) != 3:
            problems.append(f"proteins.csv: '{protein_id}' MODIFICATIONS entry '{entry}' should be "
                            "<chain>:<residue letter><position>:<CCD code>, e.g. A:S42:SEP")
            continue

        chain_letter = parts[0].strip().upper()
        residue_part = parts[1].strip().upper()
        ccd_code = parts[2].strip().upper()

        if len(residue_part) < 2 or not residue_part[0].isalpha() or not residue_part[1:].isdigit():
            problems.append(f"proteins.csv: '{protein_id}' MODIFICATIONS entry '{entry}' should give a "
                            "residue letter then a position, e.g. A:S42:SEP")
            continue

        residue_letter = residue_part[0]
        position = int(residue_part[1:])

        if ccd_code == "" or not ccd_code.isalnum() or len(ccd_code) > 5:
            problems.append(f"proteins.csv: '{protein_id}' MODIFICATIONS entry '{entry}' has an "
                            f"unusable CCD code '{ccd_code}'. Expected something like SEP or PTR.")
            continue

        if chain_letter not in chain_letters[:len(chains)]:
            problems.append(f"proteins.csv: '{protein_id}' MODIFICATIONS entry '{entry}' names chain "
                            f"'{chain_letter}', but this protein has {len(chains)} chain(s) "
                            f"({chain_letters[:len(chains)]})")
            continue

        chain_index = chain_letters.index(chain_letter)
        chain_sequence = chains[chain_index]

        if not 1 <= position <= len(chain_sequence):
            problems.append(f"proteins.csv: '{protein_id}' MODIFICATIONS entry '{entry}' position "
                            f"{position} is outside chain {chain_letter}, which is "
                            f"{len(chain_sequence)} residues long")
            continue

        # Check that the typed residue matches the actual residue in the sequence
        actual = chain_sequence[position - 1]
        if actual != residue_letter:
            problems.append(f"proteins.csv: '{protein_id}' MODIFICATIONS entry '{entry}' says residue "
                            f"{position} of chain {chain_letter} is '{residue_letter}', but the "
                            f"sequence has '{actual}' there")
            continue

        expected_parent = MODIFICATION_PARENTS.get(ccd_code)
        if expected_parent is None:
            unknown_codes.add(ccd_code)
        elif expected_parent != residue_letter:
            problems.append(f"proteins.csv: '{protein_id}' MODIFICATIONS entry '{entry}' puts "
                            f"{ccd_code} on a '{residue_letter}', but {ccd_code} is made from "
                            f"'{expected_parent}'")
            continue

        if ccd_code in MODIFICATIONS_WITHOUT_EFFECT:
            inert_codes.add(ccd_code)

        if (chain_letter, position) in seen_positions:
            problems.append(f"proteins.csv: '{protein_id}' has more than one modification at "
                            f"{chain_letter}:{position}")
            continue
        seen_positions.add((chain_letter, position))

        # A covalent ligand is attached by atom name to the residue at COVALENT_SITE.
        # Replacing that residue changes which atoms it has, so the bond would be made
        # to an atom that may no longer exist.
        if covalent and "COVALENT_SITE" in proteins.columns and row["COVALENT_SITE"] != "":
            site_text = row["COVALENT_SITE"].strip().upper()
            if ":" in site_text:
                site_chain_letter, site_residue_part = site_text.split(":", 1)
                if (site_chain_letter.strip() == chain_letter
                        and site_residue_part.strip()[1:].isdigit()
                        and int(site_residue_part.strip()[1:]) == position):
                    problems.append(f"proteins.csv: '{protein_id}' modifies {chain_letter}:{position}, "
                                    f"which is also its COVALENT_SITE. The covalent bond is made to a "
                                    f"named atom of that residue, so it cannot also be replaced.")
                    continue

        modifications_by_protein[protein_id].setdefault(chain_index, []).append(
            {"position": position, "ccd": ccd_code}
        )

if problems:
    sys.exit("ERROR: please fix these problems in the input files:\n  " + "\n  ".join(problems))

# 4g. An apo job has no ligand to attach, so a COVALENT_SITE cannot be used.
sites_ignored = 0
if apo and "COVALENT_SITE" in proteins.columns:
    sites_ignored = sum(1 for _, row in proteins.iterrows() if row["COVALENT_SITE"] != "")

print("Input checks passed")


# ---------------------------------------------------------------------------
# 5. Empty the jobs folder
# ---------------------------------------------------------------------------
# Boltz runs every .yaml file in this folder, so old files from an earlier run
# (e.g. a ligand you have since removed from the CSV) are deleted first.
jobs_dir.mkdir(exist_ok=True)
for old_file in jobs_dir.glob("*.yaml"):
    old_file.unlink()


# ---------------------------------------------------------------------------
# 6. Write one YAML per protein x ligand pair
# ---------------------------------------------------------------------------
missing_msas = set()
n_jobs = 0
n_modified_chains = 0

for _, protein in proteins.iterrows():

    # The MSA cache folder for this protein is named after a hash of the full SEQUENCE text,
    # and each chain's MSA file is named after a hash of that chain's sequence.
    complex_hash = hashlib.sha256(protein["SEQUENCE"].encode()).hexdigest()[:16]
    chains = [chain.upper() for chain in protein["SEQUENCE"].split(":")]

    cofactors = []
    if protein["COFACTORS"] != "":
        cofactors = [code.strip().upper() for code in protein["COFACTORS"].split(";") if code.strip() != ""]

    # Where a covalent ligand attaches, e.g. "A:C16" -> chain A, residue 16.
    # prepare_covalents.py has already checked the chain exists, the position is
    # inside it, and the residue there really is the letter written down, so this
    # only has to pull the two numbers back out.
    site_chain = None
    site_position = None
    if covalent:
        chain_part, residue_part = protein["COVALENT_SITE"].split(":", 1)
        site_chain = chain_part.strip().upper()
        site_position = int(residue_part.strip()[1:])   # drop the residue letter

    # An apo protein still gets one job; its single "ligand" is None.
    ligand_rows = [None] if apo else [row for _, row in ligands.iterrows()]

    for ligand in ligand_rows:

        sequences = []      # the list of chains that goes into the YAML
        letter_index = 0    # which chain letter to use next

        # Protein chains first: A, B, ...
        for chain_index, chain in enumerate(chains):
            chain_hash = hashlib.sha256(chain.encode()).hexdigest()[:16]
            msa_path = f"cache_msa/{complex_hash}/{chain_hash}.csv"   # relative to the project folder
            if not (project_dir / msa_path).exists():
                missing_msas.add(msa_path)
            chain_entry = {
                "id": chain_letters[letter_index],
                "sequence": chain,
                "msa": msa_path,
            }

            # Any modifications for this chain. The MSA is deliberately still the one
            # for the unmodified sequence: Boltz applies modifications after the
            # sequence has been used for the alignment, so the two stay consistent.
            chain_modifications = modifications_by_protein[protein["ID"]].get(chain_index, [])
            if chain_modifications:
                chain_entry["modifications"] = sorted(
                    chain_modifications, key=lambda mod: mod["position"]
                )
                n_modified_chains += 1

            sequences.append({"protein": chain_entry})
            letter_index += 1

        # Then the cofactors, given as CCD codes (e.g. ZN)
        for code in cofactors:
            sequences.append({"ligand": {
                "id": chain_letters[letter_index],
                "ccd": code,
            }})
            letter_index += 1

        # The screened ligand goes last. In a covalent experiment it goes in as its
        # bound (open-valence) form, which is what the bond constraint attaches to.
        # An apo job stops here, with only the protein chains and any cofactors.
        ligand_letter = None
        if not apo:
            ligand_letter = chain_letters[letter_index]
            if covalent:
                ligand_smiles = prepared_by_id[ligand["ID"]]["BOUND_SMILES"]
            else:
                ligand_smiles = ligand["SMILES"]
            sequences.append({"ligand": {
                "id": ligand_letter,
                "smiles": ligand_smiles,
            }})

        job = {"version": 1, "sequences": sequences}

        # The covalent bond itself. A SMILES ligand is one residue, numbered 1, and
        # its atoms are named by Boltz's own rule - see PLAN.md for both facts.
        if covalent:
            prepared_row = prepared_by_id[ligand["ID"]]
            job["constraints"] = [{"bond": {
                "atom1": [site_chain, site_position, prepared_row["TARGET_ATOM"]],
                "atom2": [ligand_letter, 1, prepared_row["ATTACHMENT_ATOM"]],
            }}]

        # Affinity is a property of a binder, so an apo job never has one.
        if not args.no_affinity and not apo:
            job["properties"] = [{"affinity": {"binder": ligand_letter}}]

        # width=100000 stops long sequences being split across lines
        if apo:
            job_file = jobs_dir / f"{protein['ID']}.yaml"
        else:
            job_file = jobs_dir / f"{protein['ID']}__{ligand['ID']}.yaml"
        with open(job_file, "w") as f:
            yaml.safe_dump(job, f, sort_keys=False, width=100000)
        n_jobs += 1


# ---------------------------------------------------------------------------
# 7. Report
# ---------------------------------------------------------------------------
print(f"Wrote {n_jobs} job files to {jobs_dir}")
if n_modified_chains:
    print(f"{n_modified_chains} protein chain(s) across these jobs carry modified residues")
if unknown_codes:
    print(f"Note: {', '.join(sorted(unknown_codes))} are not in this script's table of known "
          "modifications, so the residue they are made from was not checked. Boltz will fail "
          "the job with 'CCD component ... not found!' if a code does not exist.")
if inert_codes:
    print(f"Note: {', '.join(sorted(inert_codes))} is rewritten to the standard residue by "
          "Boltz's parser, so it has no effect on the prediction.")
if covalent:
    print("These are covalent jobs: each one bonds the ligand to its protein's COVALENT_SITE")
if apo:
    print("These are apo jobs: protein chains and cofactors only, with no affinity prediction")
    if sites_ignored:
        print(f"Note: {sites_ignored} protein(s) have a COVALENT_SITE, which an apo job cannot use. "
              "It is ignored.")
elif args.no_affinity:
    print("Affinity prediction is OFF for these jobs")

if missing_msas:
    # This is expected on the laptop or before the MSA step has run, so it is a note, not an error.
    print(f"Note: {len(missing_msas)} MSA file(s) are not in cache_msa/ yet. "
          "Boltz will fail on these jobs unless the MSA step runs first.")
