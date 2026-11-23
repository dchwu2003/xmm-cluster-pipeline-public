"""Interfaces of the server-side stages (not runnable here; the demos do not need them).

The full pipeline runs on an analysis server with XMM-SAS, HEASoft/XSPEC and the ESAS calibration
files, and processes raw observation data files (ODF, ~0.5-1 GB per observation). That code is
specific to one server and stays private; these protocols describe what each stage takes and
returns, so the parts in this repository (rules, fit-log schema, comparison, initial model) can be
read in context. Settings a real deployment needs are listed in ``config.example.yaml``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol


@dataclass
class Observation:
    obsid: str                      # 10-digit XMM ObsID
    clean_events: dict[str, str]    # camera -> filtered event file (flares removed)
    exposure_s: dict[str, float]


@dataclass
class Spectrum:
    camera: str                     # mos1 / mos2 / pn (or rass for the ROSAT sky spectrum)
    region: str                     # e.g. r015-r500, or 'outer' for the sky annulus
    files: dict[str, str]           # spectrum, QPB background, ARF, RMF
    counts: int
    area_arcmin2: float


class Ingest(Protocol):
    def process(self, obsid: str) -> Observation:
        """Download the ODF, run the SAS chains, filter soft-proton flares, check the exposure."""


class Extract(Protocol):
    def spectra(self, obs: list[Observation], ra: float, dec: float, regions: list[dict],
                sky_annulus: tuple[float, float]) -> list[Spectrum]:
        """Source and sky-annulus spectra with particle background (filter-wheel-closed, scaled by the
        unexposed corners) and responses; point sources masked."""


class Model(Protocol):
    def xcm(self, spectra: list[Spectrum], z: float, nh: float, kT0: float, band: dict) -> str:
        """The complete XSPEC model file: ICM (APEC) + sky (local bubble, Galactic halo, CXB) + soft
        protons + instrumental lines, with initial values (kT0 from the initial model)."""


class Fit(Protocol):
    def run(self, xcm_text: str) -> dict:
        """Run XSPEC on the model file and return every final parameter and the statistic. Every
        command sent to the server is logged; error/steppar/chain are refused unless authorised."""


class FitLog(Protocol):
    def log(self, record: dict) -> None:
        """Append one record following ``xmmcluster.fitlog_schema``."""
