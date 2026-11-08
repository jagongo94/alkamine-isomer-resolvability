# Alkamine isomer resolvability

Code accompanying **"Combinatorial Steroid and Drug Metabolism Hidden in Plain
Sight"** (*J. Am. Chem. Soc.*).

How often can tandem mass spectrometry tell apart two different alkamine
conjugates that reach the mass spectrometer at the same precursor mass, and
does spectral similarity track structural similarity? This repository answers
both questions directly from the deposited spectral libraries. No molecular
networking is involved: every spectrum is compared with every other spectrum in
its precursor-mass bin.

## Data

The three spectral libraries are included in `data/`, gzip-compressed (9.4 MB
in total; 33 MB uncompressed). The code reads them compressed, so nothing needs
unpacking.

| file | spectra | compounds |
| --- | ---: | ---: |
| `GNPS-ALKYLAMINES-LIBRARY.mgf.gz` | 13,854 | 6,529 |
| `GNPS-ALKYLAMINES-BILE-ACIDS-LIBRARY.mgf.gz` | 2,148 | 969 |
| `Drug_alkamines_JA.mgf.gz` | 2,029 | 1,076 |
| **total** | **18,031** | **8,511** |

Compound counts are unique SMILES within each file; the total is lower than
their sum because 63 structures appear in two files. Every entry carries a
`SMILES` field, which is the only structural input used. `data/SHA256SUMS.txt`
holds checksums; verify with `sha256sum -c SHA256SUMS.txt` from inside `data/`.

The same libraries are public on GNPS2 as `GNPS-ALKYLAMINES-LIBRARY` and
`GNPS-ALKYLAMINES-BILE-ACIDS-LIBRARY`.

## Install and run

```bash
pip install -r requirements.txt

python run_analysis.py --mgf-dir data --out results      # ~40 s
python make_figures.py --mgf-dir data --results results --out figures
```

`run_analysis.py` writes `results/best_partner.csv` (one row per distinct
spectrum with a same-formula partner) and `results/summary.txt` (the counts
quoted in the paper). `make_figures.py` writes panels a–g2 as SVG, PDF and PNG.

## Method

1. **Parse** the three MGF files; canonicalise each `SMILES` with RDKit and
   record the molecular formula, the isomeric SMILES (the compound identity)
   and the stereochemistry-collapsed SMILES (used for fingerprints).
2. **Collapse identical spectra.** Entries whose complete peak list matches
   that of a *different* compound at the same precursor mass within 5 ppm
   report one acquisition assigned to several candidate structures, not two
   independent measurements. Each such group is reduced to its first member.
3. **Bin by precursor mass.** Spectra are sorted by precursor *m/z*; each bin is
   anchored at its lowest member and takes everything within 5 ppm of it.
4. **Score each pair** of spectra from different compounds in a bin that share a
   molecular formula. Peaks within 17 Da of the precursor are removed so that
   only fragment ions contribute; intensities are square-root transformed and
   L2-normalised; fragments are matched one-to-one within 0.05 Da by greedy
   assignment, highest product first. This is the classical cosine, with no
   precursor-mass shift.
5. **Keep the best partner** for each spectrum (highest cosine, ties broken by
   structural similarity) and record the Tanimoto coefficient of their Morgan
   fingerprints (radius 1, 2048 bits, stereochemistry collapsed).
6. **Call a pair indistinguishable** at a cosine of 0.80 or above.

Parameters live at the top of `alkamine.py` and can be changed in one place.

## Expected output

```
library spectra                                 18,031
compounds (unique isomeric SMILES)               8,511

spectra sharing an identical peak list           4,640  (25.7%)
  in groups                                      2,250
distinct spectra after collapsing               15,641
distinct compounds                               7,411

spectra with a same-formula partner              9,655
  indistinguishable (cosine >= 0.80)             3,067  (19.6% of distinct spectra)
    same flat structure (stereoisomers)             20
    constitutional isomers                       3,047
  resolved (cosine < 0.80)                        6,588

compounds with a same-formula partner            4,657
  indistinguishable                              1,793  (21.1% of the library)

structural similarity of the indistinguishable constitutional isomers
  median Tanimoto                                 0.67
  fraction at Tanimoto >= 0.4                    95.3%
```

## Figures

| file | content |
| --- | --- |
| `Panel_a_cosine_vs_tanimoto` | cosine against Tanimoto, one point per distinct spectrum, coloured by the structural bin of its best partner |
| `Panel_b` … `Panel_g2` | mirror plots of representative pairs, one group per structural bin |

The example pairs are listed in `examples.tsv`; edit that file to change which
pairs are drawn. A blank `target_cosine` selects the most fragment-rich
comparison for that compound pair, and a value selects the acquisition closest
to it, which matters because several compounds have more than one spectrum.

## Notes for reviewers

- **Why identical peak lists are collapsed.** 4,640 of 18,031 spectra (25.7%)
  share a peak list with a different compound. These are single acquisitions
  deposited under two or more candidate structures, with the alternative
  recorded in the library's `INCHIAUX` field; in most cases the two candidates
  name the same amine and carbonyl and differ only in where the bond formed.
  Scoring them would compare a spectrum with itself and inflate the number of
  "indistinguishable" pairs, so they are removed. The reported figure is
  therefore conservative.
- **Why the fingerprint radius is 1.** Morgan radius 2 (ECFP4) and radius 3
  score constitutional isomers lower, moving the same pairs into less similar
  bins, while radius 0 sees only atom composition and calls every isomer pair
  identical. Set `MORGAN_RADIUS` in `alkamine.py` to repeat the analysis at any
  radius; the conclusions about the cosine are unaffected, only the bin labels
  move.
- **Fragment tolerance.** The analysis uses 0.05 Da, matching the setting used
  for library matching and molecular networking elsewhere in the study. The
  result is insensitive to the choice: at 0.01, 0.02 and 0.05 Da the number of
  indistinguishable compounds is 1,760, 1,765 and 1,793 (20.7%, 20.7% and
  21.1% of the library), so the conclusion does not depend on it. Change
  `FRAGMENT_TOL` in `alkamine.py` to repeat the test.
- **Reproducibility.** The analysis is deterministic. Greedy fragment matching
  resolves ties by descending intensity product, and compound identity is the
  RDKit canonical isomeric SMILES rather than a name, so duplicate or ambiguous
  compound names do not affect the counts.

## Citation

> Agongo, J. *et al.* Combinatorial Steroid and Drug Metabolism Hidden in Plain
> Sight. *J. Am. Chem. Soc.* (2026).

## License

MIT, see `LICENSE`.
