# Server-side stages (not included)

The stages that touch raw data run on an analysis server with XMM-SAS, HEASoft/XSPEC and the ESAS
calibration files: ingest (ODF download, SAS chains, flare filtering), extraction (source and
sky-annulus spectra, particle background from filter-wheel-closed data, responses), model building
(XSPEC model files from `xmmcluster.rules` and the initial model) and fitting. That code is written
for one server and its file layout, so it is not part of this repository.

`interfaces.py` lists what each stage takes and returns. `config.example.yaml` at the top level lists
the settings a deployment needs, with placeholders only. The pieces that do not need the server are
in `xmmcluster/`, and the demos run on their output (`data/`).
