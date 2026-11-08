#!/usr/bin/env python3
"""
Alkamine isomer resolvability: how often can MS/MS tell same-mass library
members apart, and does spectral similarity track structural similarity?

Usage
-----
    python run_analysis.py --mgf-dir data --out results

Writes
------
    best_partner.csv   one row per distinct spectrum that has a same-formula
                       partner, with its best-matching partner and the scores
    summary.txt        the counts quoted in the manuscript
"""
from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import pandas as pd

import alkamine as ak

LIBRARIES = {
    "ALKYLAMINES": "GNPS-ALKYLAMINES-LIBRARY.mgf",
    "BILE-ACIDS": "GNPS-ALKYLAMINES-BILE-ACIDS-LIBRARY.mgf",
    "DRUG": "Drug_alkamines_JA.mgf",
}


def best_partners(spectra: list[ak.Spectrum], fingerprints: dict) -> pd.DataFrame:
    """
    For each spectrum, find the spectrum of a *different* compound that shares
    its precursor mass and molecular formula and gives the highest cosine.
    """
    prepared = {s.index: ak.preprocess(s) for s in spectra}
    best: dict[int, dict] = {}
    has_same_mass_partner: set[int] = set()

    for group in ak.bin_by_precursor(spectra):
        if len(group) < 2:
            continue
        for first, second in itertools.combinations(group, 2):
            if first.isomeric_smiles == second.isomeric_smiles:
                continue  # replicate spectra of one compound
            has_same_mass_partner.update((first.index, second.index))
            if first.formula != second.formula:
                continue  # same mass, different formula: an adduct coincidence

            score, matched = ak.cosine(prepared[first.index], prepared[second.index])
            similarity = ak.tanimoto(
                fingerprints, first.flat_smiles, second.flat_smiles
            )
            same_flat = first.flat_smiles == second.flat_smiles

            for this, other in ((first, second), (second, first)):
                previous = best.get(this.index)
                better = previous is None or (score, similarity) > (
                    previous["cosine"],
                    previous["tanimoto"],
                )
                if better:
                    best[this.index] = {
                        "spectrum": this.index,
                        "partner": other.index,
                        "cosine": score,
                        "matched_peaks": matched,
                        "tanimoto": similarity,
                        "same_flat_structure": same_flat,
                        "name": this.name,
                        "partner_name": other.name,
                        "smiles": this.smiles,
                        "partner_smiles": other.smiles,
                        "formula": this.formula,
                        "precursor_mz": this.precursor_mz,
                        "library": this.library,
                    }

    table = pd.DataFrame(best.values())
    table.attrs["same_mass_partner"] = has_same_mass_partner
    return table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mgf-dir", default="data", type=Path,
                        help="directory holding the three library MGF files")
    parser.add_argument("--out", default="results", type=Path,
                        help="directory for the output files")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    paths = {name: args.mgf_dir / filename for name, filename in LIBRARIES.items()}
    missing = []
    for name, path in paths.items():
        try:
            paths[name] = ak.resolve(path)
        except FileNotFoundError:
            missing.append(f"{path} (or {path}.gz)")
    if missing:
        raise SystemExit("missing library file(s):\n  " + "\n  ".join(missing))

    print("reading libraries ...")
    spectra = ak.read_libraries(paths)
    ak.annotate_structures(spectra)
    compounds = {s.isomeric_smiles for s in spectra if s.isomeric_smiles}
    print(f"  {len(spectra):,} spectra, {len(compounds):,} compounds")

    print("collapsing identical peak lists ...")
    distinct, collapse = ak.collapse_identical_spectra(spectra)
    distinct_compounds = {s.isomeric_smiles for s in distinct if s.isomeric_smiles}
    print(f"  {collapse['spectra']:,} spectra in {collapse['groups']:,} groups "
          f"({collapse['compounds']:,} compounds) collapsed")
    print(f"  {len(distinct):,} distinct spectra, "
          f"{len(distinct_compounds):,} compounds remain")

    print("scoring same-formula pairs within 5 ppm precursor bins ...")
    fingerprints = ak.fingerprint_table(distinct)
    table = best_partners(distinct, fingerprints)
    table.to_csv(args.out / "best_partner.csv", index=False)

    indistinguishable = table[table.cosine >= ak.COSINE_CUTOFF]
    by_compound = (
        table.assign(compound=table.spectrum.map(
            {s.index: s.isomeric_smiles for s in distinct}))
        .sort_values(["compound", "cosine", "tanimoto"], ascending=[True, False, False])
        .groupby("compound").first()
    )
    compound_hits = by_compound[by_compound.cosine >= ak.COSINE_CUTOFF]

    with_mass_only = sum(
        1 for s in distinct
        if s.index in table.attrs["same_mass_partner"]
        and s.index not in set(table.spectrum)
    )
    no_partner = len(distinct) - len(table) - with_mass_only

    lines = [
        "Alkamine isomer resolvability",
        "=" * 60,
        f"library spectra                               {len(spectra):>8,}",
        f"compounds (unique isomeric SMILES)            {len(compounds):>8,}",
        "",
        f"spectra sharing an identical peak list        {collapse['spectra']:>8,}"
        f"  ({100 * collapse['spectra'] / len(spectra):.1f}%)",
        f"  in groups                                   {collapse['groups']:>8,}",
        f"distinct spectra after collapsing             {len(distinct):>8,}",
        f"distinct compounds                            {len(distinct_compounds):>8,}",
        "",
        f"spectra with a same-formula partner           {len(table):>8,}",
        f"  indistinguishable (cosine >= "
        f"{ak.COSINE_CUTOFF:.2f})          {len(indistinguishable):>8,}"
        f"  ({100 * len(indistinguishable) / len(distinct):.1f}% of distinct spectra)",
        f"    same flat structure (stereoisomers)       "
        f"{int(indistinguishable.same_flat_structure.sum()):>8,}",
        f"    constitutional isomers                    "
        f"{int((~indistinguishable.same_flat_structure).sum()):>8,}",
        f"  resolved (cosine < {ak.COSINE_CUTOFF:.2f})                     "
        f"{int((table.cosine < ak.COSINE_CUTOFF).sum()):>8,}",
        "",
        f"spectra with a same-mass but different-formula partner  "
        f"{with_mass_only:>8,}",
        f"spectra with no partner of a different compound         {no_partner:>8,}",
        "",
        f"compounds with a same-formula partner         {len(by_compound):>8,}",
        f"  indistinguishable                           {len(compound_hits):>8,}"
        f"  ({100 * len(compound_hits) / len(compounds):.1f}% of the library)",
        "",
        "structural similarity of the indistinguishable constitutional isomers",
        f"  median Tanimoto                             "
        f"{indistinguishable[~indistinguishable.same_flat_structure].tanimoto.median():>8.2f}",
        f"  fraction at Tanimoto >= 0.4                 "
        f"{(indistinguishable[~indistinguishable.same_flat_structure].tanimoto >= 0.4).mean():>8.1%}",
    ]
    report = "\n".join(lines)
    print("\n" + report)
    (args.out / "summary.txt").write_text(report + "\n")
    print(f"\nwritten to {args.out}/")


if __name__ == "__main__":
    main()
