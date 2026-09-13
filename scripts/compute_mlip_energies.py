"""Single-point MLIP energies for the generated crystal structures in this repository.

Self-contained: needs only the packages in scripts/requirements-mlip.txt and a
MACE model file. Reads data/gen_crystals_<generator>.csv (a ``cif`` column),
checks that no two atoms are closer than 0.7 A (the LeMat-GenBench validity
cutoff), evaluates every structure with the MACE model, and adds Materials
Project MP2020 corrections with the run type and Hubbard U values that
MPRelaxSet would assign to the composition. Writes one CSV per generator to
data/energies/<model-name>/ plus a manifest.json with checksums.

Usage:
    python scripts/compute_mlip_energies.py --model models/mace-omat-0-medium.model
"""

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
import warnings
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path

import numpy as np
from pymatgen.core import Structure
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp.sets import MPRelaxSet

REPO = Path(__file__).resolve().parents[1]
MIN_DISTANCE_ANGSTROM = 0.7
MODEL_SOURCES = {
    "mace-omat-0-medium": {
        "url": "https://github.com/ACEsuit/mace-foundations/releases/download/mace_omat_0/mace-omat-0-medium.model",
        "sha256": "d4b14be9afa294eebdbe31a0280b26a0fa29715771e978cbfd3ca24e0d90307a",
        "license": "ASL (Academic Software Licence)",
        "training_data": "OMat24 (PBE/PBE+U)",
    }
}
COLUMNS = [
    "index",
    "status",
    "formula",
    "chemsys",
    "nelements",
    "nsites",
    "composition",
    "min_distance_angstrom",
    "energy_total_eV",
    "energy_per_atom_eV",
    "uncorrected_energy_total_eV",
    "mp2020_correction_eV",
    "mp2020_adjustments",
    "run_type",
    "hubbards",
    "error",
]


def relative_to_repo(path):
    path = Path(path).resolve()
    return str(path.relative_to(REPO)) if path.is_relative_to(REPO) else str(path)


def sha256(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def parse_structure(cif):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Structure.from_str(cif, fmt="cif")


def min_distance(structure):
    if len(structure) < 2:
        return math.inf
    dm = structure.distance_matrix
    np.fill_diagonal(dm, np.inf)
    return float(dm.min())


def mp_parameters(structure):
    """run_type and Hubbard U values as Materials Project's MPRelaxSet assigns them."""
    vis = MPRelaxSet(structure)
    incar = vis.incar
    symbols = vis.poscar.site_symbols
    if incar.get("LDAU"):
        hubbards = {s: float(u) for s, u in zip(symbols, incar["LDAUU"])}
        return {"run_type": "GGA+U", "is_hubbard": True, "hubbards": hubbards}
    return {"run_type": "GGA", "is_hubbard": False, "hubbards": {}}


def evaluate(structure, calc, compat, entry_id):
    atoms = AseAtomsAdaptor.get_atoms(structure)
    atoms.calc = calc
    raw = float(atoms.get_potential_energy())
    if not math.isfinite(raw):
        raise ValueError("non-finite MLIP energy")
    parameters = mp_parameters(structure)
    entry = ComputedStructureEntry(structure, raw, parameters=parameters, entry_id=entry_id)
    processed = compat.process_entries([entry], clean=True, inplace=False)
    if not processed:
        raise ValueError("MP2020 rejected the entry")
    entry = processed[0]
    return {
        "energy_total_eV": entry.energy,
        "energy_per_atom_eV": entry.energy_per_atom,
        "uncorrected_energy_total_eV": entry.uncorrected_energy,
        "mp2020_correction_eV": entry.correction,
        "mp2020_adjustments": json.dumps(
            [{"name": a.name, "value": a.value} for a in entry.energy_adjustments]
        ),
        "run_type": parameters["run_type"],
        "hubbards": json.dumps(parameters["hubbards"]) if parameters["hubbards"] else "",
    }


def process_file(path, calc, compat, min_dist, limit):
    with path.open() as f:
        cifs = [row["cif"] for row in csv.DictReader(f)]
    if limit:
        cifs = cifs[:limit]
    rows, counts = [], {}
    for i, cif in enumerate(cifs):
        row = {"index": i}
        try:
            s = parse_structure(cif)
            row.update(
                formula=s.composition.reduced_formula,
                chemsys=s.composition.chemical_system,
                nelements=len(s.composition.elements),
                nsites=len(s),
                composition=json.dumps(s.composition.as_dict()),
                min_distance_angstrom=min_distance(s),
            )
            if row["min_distance_angstrom"] < min_dist:
                raise ValueError("too_close")
            row.update(evaluate(s, calc, compat, f"{path.stem}-{i}"))
            row["status"] = "ok"
        except Exception as exc:  # record the failure kind and continue
            row["status"] = "too_close" if str(exc) == "too_close" else type(exc).__name__
            row["error"] = str(exc)[:200]
        counts[row["status"]] = counts.get(row["status"], 0) + 1
        rows.append(row)
        if (i + 1) % 500 == 0:
            print(f"{path.name}: {i + 1}/{len(cifs)}", flush=True)
    return rows, counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, required=True, help="MACE model file")
    parser.add_argument(
        "--model-name",
        help="Directory name under data/energies/ (default: model file stem)",
    )
    parser.add_argument(
        "--inputs",
        type=Path,
        nargs="+",
        default=sorted((REPO / "data").glob("gen_crystals_*.csv")),
    )
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--dtype", default="float64", choices=["float32", "float64"])
    parser.add_argument("--min-distance", type=float, default=MIN_DISTANCE_ANGSTROM)
    parser.add_argument("--limit", type=int, help="Structures per file (testing)")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    model_name = args.model_name or args.model.stem
    output = args.output_dir or REPO / "data/energies" / model_name
    output.mkdir(parents=True, exist_ok=True)
    model_sha = sha256(args.model)
    source = MODEL_SOURCES.get(model_name, {})
    if source and source["sha256"] != model_sha:
        raise ValueError("Model file does not match the recorded SHA256")

    from mace.calculators import MACECalculator

    calc = MACECalculator(
        model_paths=[str(args.model)], device=args.device, default_dtype=args.dtype
    )
    compat = MaterialsProject2020Compatibility(check_potcar=False)
    manifest = {
        "status": "in_progress",
        "started_at_utc": datetime.now(UTC).isoformat(),
        "model_name": model_name,
        "model_sha256": model_sha,
        "model_source": source,
        "device": args.device,
        "dtype": args.dtype,
        "relaxation": "none; single point on the generated structure",
        "min_distance_angstrom": args.min_distance,
        "energy_convention": (
            "energy_total_eV = MLIP total energy + MP2020 corrections (pymatgen "
            "MaterialsProject2020Compatibility with run type and Hubbard U from "
            "MPRelaxSet); comparable only with references on the same footing"
        ),
        "packages": {n: version(n) for n in ["mace-torch", "torch", "pymatgen", "ase", "numpy"]},
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "script_sha256": sha256(Path(__file__)),
        "files": {},
    }
    manifest_path = output / "manifest.json"
    for path in args.inputs:
        rows, counts = process_file(path, calc, compat, args.min_distance, args.limit)
        target = output / f"{path.stem}_energies.csv"
        with target.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
        manifest["files"][target.name] = {
            "input": relative_to_repo(path),
            "input_sha256": sha256(path),
            "n_rows": len(rows),
            "status_counts": counts,
            "csv_sha256": sha256(target),
        }
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        print(f"{target.name}: {counts}", flush=True)
    manifest["status"] = "complete"
    manifest["finished_at_utc"] = datetime.now(UTC).isoformat()
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
