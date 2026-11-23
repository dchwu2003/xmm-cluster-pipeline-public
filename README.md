# XMM-Newton cluster spectroscopy: reproducing the literature, learning initial models

X-ray spectral fitting of galaxy clusters is still largely done by hand. For every cluster an expert
chooses the region, the background treatment, the energy band and the starting XSPEC model, and two
careful analyses of the same observation often differ by 10-15 %. This project writes that judgement
down as explicit rules, checks the rules by reproducing published temperatures with one code base and
no tuning toward the published numbers, and logs every fit as training data for a model that proposes
the starting point of the next fit.

*Work in progress.* <!-- numbers:snapshot -->
Snapshot of the data: 2026-11-22.
<!-- /numbers -->

```mermaid
flowchart LR
    A[ingest<br/>ODF, SAS chains,<br/>flare filtering] --> B[extract<br/>source + sky spectra,<br/>particle background]
    B --> C[model<br/>XSPEC model from<br/>rules + initial kT]
    C --> D[fit<br/>XSPEC]
    D --> E[(fit log<br/>every fit, every<br/>parameter)]
    E --> F[reproduce<br/>compare with<br/>published values]
    E --> G[learn<br/>initial-model<br/>prediction]
    G --> C
    F -. offsets reported,<br/>rules revised .-> C
```

## Rules for the agent

Most of the code was written by an AI coding agent (see [Roles](#roles)). It worked under these rules:

- **Published values are targets only.** They are stored next to each fit for the comparison; no fit,
  rule or model input reads them, and no setting is changed to move a result toward them.
- **Quality rejection never sees the published value.** A cluster is left out only by data-quality rules
  (`xmmcluster/rules.py`: flares, untrusted fits, temperatures at the fit limit), applied whatever its ratio.
- **Every server command is logged**, and the expensive parameter-search commands (`error`,
  `steppar`, `chain`) are refused unless authorised.
- **Offsets are reported, not tuned away.** Where a paper's method cannot be followed, the departure is
  written down and the offset is kept in the statistics.

## Results so far

<!-- numbers:headline -->
- **9 papers reproduced**: 598 temperatures of 234 paper-cluster pairs (160 clusters). Median ratio ours / published **0.98**; **70%** within 10 % and **90%** within 20 %. 39 further values are left out by the data-quality rules (rejected data, or kT at the fit limit), whatever their ratio.
- **Initial-model prototype**: from band fluxes measured before any fit, extra trees predict the starting kT of a cluster they never saw (4,957 spectra of 189 clusters, leave-one-cluster-out) with a median error of **7.2%**; 85% are within 20 %. On the 3,537 spectra (114 clusters) that have a catalogue temperature, the prediction's median error is 6.3% against 13.6% for the catalogue kT; the median kT of the training clusters is off by 32.3%.
<!-- /numbers -->

![Reproduced temperatures per paper](docs/figures/reproductions.png)

<!-- numbers:papers -->
| paper | compared | clusters | values | median ratio | 16-84 % | within 10 % | within 20 % |
|---|---|---|---|---|---|---|---|
| Zhang et al. 2008 (LoCuSS) | kT in 0.2-0.5 r500 | 24 | 24 | 0.94 | 0.89-1.01 | 75% | 100% |
| Planck Collaboration 2011 XI | T_X in 0.15-0.75 R500 | 32 | 32 | 0.97 | 0.89-1.08 | 69% | 91% |
| Planck Collaboration 2011 IX, 2012 I, 2013 IV | T_X in 0.15-0.75 R500 | 27 | 27 | 0.95 | 0.87-1.02 | 59% | 89% |
| Lovisari et al. 2015 | group kT | 11 | 22 | 1.14 | 0.97-1.31 | 45% | 68% |
| Lovisari et al. 2020 | kT in 0-1 and 0.15-1 R500 | 54 | 108 | 1.03 | 0.95-1.12 | 66% | 88% |
| Rossetti et al. 2024 (CHEX-MATE) | kT profiles (annuli) | 26 | 277 | 0.96 | 0.89-1.03 | 73% | 93% |
| Johnson et al. 2009 (2dXGS) | group kT in 0.1-0.3 r500 | 23 | 23 | 1.13 | 0.99-1.41 | 43% | 65% |
| Pratt et al. 2009 (REXCESS) | kT in 0.15-1 R500 | 26 | 74 | 0.96 | 0.91-1.03 | 84% | 99% |
| O'Sullivan et al. 2017 (CLoGS) | group system kT | 11 | 11 | 0.98 | 0.93-1.03 | 73% | 82% |
<!-- /numbers -->

![Predicted against fitted kT](docs/figures/initmodel.png)

Median error of the starting kT by fitted temperature (median |ln(kT_start / kT_fit)|, cluster held out):

<!-- numbers:bins -->
| fitted kT | spectra | clusters | prediction | catalogue kT |
|---|---|---|---|---|
| 0-1 keV | 289 | 27 | 10.1% | 19.3% |
| 1-2 keV | 412 | 36 | 8.7% | 21.8% |
| 2-4 keV | 619 | 50 | 8.1% | 19.0% |
| 4-8 keV | 2,641 | 116 | 5.1% | 10.6% |
| 8-12 keV | 858 | 59 | 10.9% | 16.1% |
| 12 keV and above | 138 | 17 | 36.3% | 28.8% |
<!-- /numbers -->

The prediction needs no catalogue, so it also works for clusters that are in none. Its weak end is the
hottest clusters, which are few in the training set. The catalogue baseline is the MCXC mass (from
ROSAT luminosities) turned into kT with the M-T relation of Arnaud et al. 2005, so it is a rough
starting point rather than a spectroscopic temperature.

## Quick start

```bash
pip install -e .
python -m demo.reproductions   # per-paper table and docs/figures/reproductions.png
python -m demo.initmodel       # leave-one-cluster-out training and docs/figures/initmodel.png
```

The initial-model demo takes a minute or two on a laptop; `--folds 10` (ten folds of whole clusters) is faster.
`pytest` runs the tests. `python -m demo.readme_numbers` rewrites the numbers above from the bundled CSVs.

## What is in this repository

- `xmmcluster/compare.py`: the per-paper comparison and the left-out rules.
- `xmmcluster/rules.py`: the decision rules of the automatic fit (bands, screening, regions, quality flags, verdict).
- `xmmcluster/fitlog_schema.py`: the fit-log record and its validation.
- `xmmcluster/initmodel/`: features, extra trees with leave-one-cluster-out evaluation, the replay
  summary, and the band-flux measurement (reference; it needs the spectra).
- `data/summary/reproduction_points.csv`: one row per compared value (ours, published, ratio, data-quality verdict).
- `data/initmodel/spectra.csv`: one row per fitted spectrum (band fluxes, z, N_H, region, fitted kT and Z,
  catalogue kT). Column definitions are in `xmmcluster/initmodel/data.py`.
- `pipeline/`: interfaces of the server-side stages and `config.example.yaml`, with placeholders only.

## What is not in this repository, and why

- **Server-side code** (ingest, extraction, XSPEC execution, job control). It is written for one analysis
  server and its file layout. `pipeline/interfaces.py` describes what each stage does.
- **Raw and intermediate data** (observation data files, event files and spectra, far too large for a
  repository) and the full fit log, which holds server paths. The bundled tables are derived from them.
- **Data outside the scope of this snapshot.** It includes only clusters whose observations were
  downloaded from the XMM-Newton Science Archive and processed end to end by this pipeline. Results on
  data processed separately before this pipeline existed are left out, together with the papers that
  were reproduced only on those data.
- **Unpublished results**, including work in preparation on individual clusters.

## Roles

I (Daocheng Wu) set the design and the direction: what to reproduce and in which order, the rules for
the agent above, and the acceptance checks a reproduction has to pass before it counts. The code was
written largely with Claude Code (an AI coding agent) under my direction. I reviewed the results,
decided how to treat offsets and failures, and chose what is published here.

The project is work in progress: more papers are being reproduced, the coolest groups and the hottest
clusters, the thin ends of the training set, are being added, and the server-side code is not yet in a shareable form.

## License

MIT (see `LICENSE`).
