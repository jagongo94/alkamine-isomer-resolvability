"""
Shared routines for the alkamine isomer-resolvability analysis.

Reproduces the comparison of every MS/MS spectrum in the alkamine spectral
libraries against every other spectrum of the same precursor mass, and relates
spectral similarity (cosine) to structural similarity (Tanimoto).

Reference: Agongo et al., "Combinatorial Steroid and Drug Metabolism Hidden in
Plain Sight", J. Am. Chem. Soc.
"""
from __future__ import annotations

import collections
import gzip
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")

# ---------------------------------------------------------------- parameters
PRECURSOR_PPM = 5.0      # precursor grouping tolerance
PRECURSOR_WINDOW = 17.0  # Da around the precursor removed before scoring
FRAGMENT_TOL = 0.05      # Da, fragment matching tolerance
COSINE_CUTOFF = 0.80     # at or above this, a pair is called indistinguishable
MORGAN_RADIUS = 1
MORGAN_BITS = 2048


@dataclass
class Spectrum:
    """One library entry: a peak list plus its structural annotation."""

    index: int
    library: str
    name: str
    smiles: str
    precursor_mz: float
    mz: np.ndarray
    intensity: np.ndarray
    isomeric_smiles: str | None = None
    flat_smiles: str | None = None
    formula: str | None = None
    meta: dict = field(default_factory=dict)


# -------------------------------------------------------------------- input
def _open_text(path: Path):
    """Open a plain or gzip-compressed text file."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, encoding="utf-8", errors="ignore")


def resolve(path: str | Path) -> Path:
    """Accept either `name.mgf` or `name.mgf.gz`, whichever is present."""
    path = Path(path)
    if path.exists():
        return path
    gzipped = path.with_suffix(path.suffix + ".gz")
    if gzipped.exists():
        return gzipped
    raise FileNotFoundError(f"neither {path} nor {gzipped} exists")


def read_mgf(path: str | Path, library: str, start_index: int = 0) -> list[Spectrum]:
    """Parse one MGF file (optionally gzip-compressed) into Spectrum objects."""
    spectra: list[Spectrum] = []
    header: dict[str, str] = {}
    mz: list[float] = []
    inten: list[float] = []
    inside = False

    with _open_text(resolve(path)) as handle:
        for raw in handle:
            line = raw.strip()
            if line == "BEGIN IONS":
                inside, header, mz, inten = True, {}, [], []
            elif line == "END IONS":
                if inside:
                    spectra.append(
                        Spectrum(
                            index=start_index + len(spectra),
                            library=library,
                            name=header.get("NAME", ""),
                            smiles=header.get("SMILES", "").strip(),
                            precursor_mz=float(header.get("PEPMASS", "0").split()[0]),
                            mz=np.asarray(mz, dtype=float),
                            intensity=np.asarray(inten, dtype=float),
                            meta=header,
                        )
                    )
                inside = False
            elif inside:
                if "=" in line and not line[0].isdigit():
                    key, value = line.split("=", 1)
                    header[key] = value
                elif line and (line[0].isdigit() or line[0] == "."):
                    parts = line.split()
                    if len(parts) >= 2:
                        mz.append(float(parts[0]))
                        inten.append(float(parts[1]))
    return spectra


def read_libraries(paths: dict[str, str | Path]) -> list[Spectrum]:
    """Parse several MGF files, keeping a single continuous index."""
    spectra: list[Spectrum] = []
    for library, path in paths.items():
        spectra.extend(read_mgf(path, library, start_index=len(spectra)))
    return spectra


# --------------------------------------------------------------- structures
def annotate_structures(spectra: list[Spectrum]) -> None:
    """Attach canonical SMILES (with and without stereochemistry) and formula."""
    cache: dict[str, tuple[str, str, str] | None] = {}
    for spectrum in spectra:
        smiles = spectrum.smiles
        if smiles not in cache:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                cache[smiles] = None
            else:
                isomeric = Chem.MolToSmiles(mol, isomericSmiles=True)
                flat = Chem.MolToSmiles(
                    Chem.MolFromSmiles(Chem.MolToSmiles(mol, isomericSmiles=False))
                )
                cache[smiles] = (isomeric, flat, CalcMolFormula(mol))
        annotation = cache[smiles]
        if annotation is not None:
            spectrum.isomeric_smiles, spectrum.flat_smiles, spectrum.formula = annotation


def fingerprint_table(spectra: list[Spectrum]) -> dict[str, object]:
    """Morgan fingerprints keyed by stereochemistry-collapsed SMILES."""
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=MORGAN_RADIUS, fpSize=MORGAN_BITS
    )
    fingerprints: dict[str, object] = {}
    for spectrum in spectra:
        key = spectrum.flat_smiles
        if key and key not in fingerprints:
            mol = Chem.MolFromSmiles(key)
            fingerprints[key] = generator.GetFingerprint(mol) if mol else None
    return fingerprints


def tanimoto(fingerprints: dict, a: str, b: str) -> float:
    """Tanimoto coefficient between two stereochemistry-collapsed structures."""
    fp_a, fp_b = fingerprints.get(a), fingerprints.get(b)
    if fp_a is None or fp_b is None:
        return float("nan")
    return DataStructs.TanimotoSimilarity(fp_a, fp_b)


# ------------------------------------------------------- duplicate collapse
def collapse_identical_spectra(
    spectra: list[Spectrum], ppm: float = PRECURSOR_PPM
) -> tuple[list[Spectrum], dict[str, int]]:
    """
    Keep one spectrum per group of byte-identical peak lists.

    Entries whose complete peak list matches that of a different compound at the
    same precursor mass report a single acquisition assigned to two or more
    candidate structures, not two independent measurements. Scoring them against
    each other would compare a spectrum with itself, so each group is reduced to
    its first member in library order.
    """
    groups: dict[tuple, list[Spectrum]] = collections.defaultdict(list)
    for spectrum in spectra:
        key = (tuple(np.round(spectrum.mz, 6)), tuple(np.round(spectrum.intensity, 6)))
        groups[key].append(spectrum)

    kept: list[Spectrum] = []
    stats = {"groups": 0, "spectra": 0, "compounds": 0}
    collapsed_compounds: set[str] = set()

    for members in groups.values():
        structures = {m.isomeric_smiles for m in members if m.isomeric_smiles}
        if len(members) > 1 and len(structures) > 1:
            masses = np.array([m.precursor_mz for m in members])
            spread_ppm = (masses.max() - masses.min()) / masses.mean() * 1e6
            if spread_ppm <= ppm:
                stats["groups"] += 1
                stats["spectra"] += len(members)
                collapsed_compounds |= structures
                kept.append(min(members, key=lambda m: m.index))
                continue
        kept.extend(members)

    stats["compounds"] = len(collapsed_compounds)
    return kept, stats


# ------------------------------------------------------------------ binning
def bin_by_precursor(
    spectra: list[Spectrum], ppm: float = PRECURSOR_PPM
) -> list[list[Spectrum]]:
    """
    Group spectra by precursor m/z.

    Spectra are sorted by precursor mass and each bin is anchored at its lowest
    member, taking everything within `ppm` of that member.
    """
    ordered = sorted(spectra, key=lambda s: s.precursor_mz)
    bins: list[list[Spectrum]] = []
    current: list[Spectrum] = [ordered[0]]
    anchor = ordered[0].precursor_mz

    for spectrum in ordered[1:]:
        if spectrum.precursor_mz - anchor <= anchor * ppm * 1e-6:
            current.append(spectrum)
        else:
            bins.append(current)
            current = [spectrum]
            anchor = spectrum.precursor_mz
    bins.append(current)
    return bins


# ------------------------------------------------------------------- cosine
def preprocess(spectrum: Spectrum) -> tuple[np.ndarray, np.ndarray]:
    """
    Prepare a peak list for scoring.

    Peaks within `PRECURSOR_WINDOW` of the precursor are discarded so that the
    precursor and its immediate neutral losses cannot inflate the score.
    Intensities are square-root transformed and the vector is L2-normalised.
    """
    keep = np.abs(spectrum.mz - spectrum.precursor_mz) > PRECURSOR_WINDOW
    mz = spectrum.mz[keep]
    intensity = np.sqrt(spectrum.intensity[keep])
    norm = np.linalg.norm(intensity)
    order = np.argsort(mz)
    if norm == 0:
        return mz[order], intensity[order]
    return mz[order], (intensity / norm)[order]


def cosine(
    peaks_a: tuple[np.ndarray, np.ndarray],
    peaks_b: tuple[np.ndarray, np.ndarray],
    tolerance: float = FRAGMENT_TOL,
) -> tuple[float, int]:
    """
    Classical cosine similarity with greedy one-to-one fragment matching.

    Returns the score and the number of matched fragment pairs. There is no
    precursor-mass shift: only fragments at the same m/z are matched.
    """
    mz_a, int_a = peaks_a
    mz_b, int_b = peaks_b
    if len(mz_a) == 0 or len(mz_b) == 0:
        return 0.0, 0

    low = np.searchsorted(mz_b, mz_a - tolerance)
    high = np.searchsorted(mz_b, mz_a + tolerance)
    candidates = sorted(
        (
            (int_a[i] * int_b[j], i, j)
            for i in range(len(mz_a))
            for j in range(low[i], high[i])
        ),
        reverse=True,
    )

    used_a: set[int] = set()
    used_b: set[int] = set()
    score = 0.0
    matched = 0
    for product, i, j in candidates:
        if i in used_a or j in used_b:
            continue
        used_a.add(i)
        used_b.add(j)
        score += product
        matched += 1
    return score, matched
