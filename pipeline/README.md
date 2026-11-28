# Server-side stages (not included)

The stages that touch raw data run on an analysis server with XMM-SAS, HEASoft/XSPEC and the ESAS
calibration files:

- **ingest**: ODF download from the archive, SAS chains, soft-proton flare filtering;
- **automatic fit**: regions from r500, recentring on the X-ray centroid, screening of every exposure,
  source and sky-annulus spectra with particle background from filter-wheel-closed data, the XSPEC
  model from `xmmcluster.rules`, joint fits of all kept cameras, r500 iterated from the fitted kT, and
  a verdict (ok / check / rejected) from the quality flags;
- **sky spectrum**: a ROSAT All-Sky Survey spectrum around the cluster, fetched from HEASARC and
  fitted with the XMM spectra;
- **new sample**: one command per step from a CSV of clusters (name, position, redshift) to the
  verdicts: pick the pointings, ingest, automatic fit, sky spectrum, refit, report. Each step skips
  what is already done, so an interrupted run is restarted as it is.

That code is written for one server and its file layout, so it is not part of this repository.
`interfaces.py` lists what each stage takes and returns. `config.example.yaml` at the top level lists
the settings a deployment needs, with placeholders only. The pieces that do not need the server are
in `xmmcluster/`, and the demos run on their output (`data/`).

## First run on a new sample

The whole chain was run once on three clusters that are in none of the reproduced samples, given only
their names, positions and redshifts (2026-11-26). It took 3 h 20 min of wall-clock time on a shared
server, archive download and SAS included, with no manual step:

- A1068 (z 0.139) and ZwCl 3146 (z 0.291): **ok**. For A1068 one MOS camera was dropped for soft
  protons.
- A2626 (z 0.055): **rejected**, as designed for nearby clusters: its r500 (about 16') does not fit in
  the field, so the sky annulus contains cluster emission. Such clusters need an explicitly set region.

An "ok" means that no quality flag was raised, not that the temperature was checked against another
measurement. These clusters are held out of the initial-model training and are not in the bundled data.
