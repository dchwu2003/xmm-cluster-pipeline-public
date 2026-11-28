"""Interfaces of the server-side stages (not runnable here; the demos do not need them).

The full pipeline runs on an analysis server with XMM-SAS, HEASoft/XSPEC and the ESAS calibration
files, and processes raw observation data files (ODF, ~0.5-1 GB per observation). That code is
specific to one server and stays private; these protocols describe what each stage takes and
returns, so the parts in this repository (rules, fit-log schema, comparison, initial model) can be
read in context. Settings a real deployment needs are listed in ``config.example.yaml``.

The stages run in this order for a new list of clusters (``NewSample``): pick the observations,
ingest, the automatic fit (regions, screening, model, fit, r500 iteration, verdict), a ROSAT sky
spectrum where the field lacks one, a refit with it, and a report. Every fit goes to the fit log.
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
class ClusterInput:
    name: str
    ra: float                       # degrees; the X-ray centroid is searched within 4' of it
    dec: float
    z: float
    obsids: list[str]               # empty: chosen from the archive
    kT0: float | None = None        # starting kT (keV); without it 5 keV, then r500 follows the fit


@dataclass
class Spectrum:
    camera: str                     # mos1 / mos2 / pn (or rass for the ROSAT sky spectrum)
    region: str                     # e.g. r015-r500, or 'outer' for the sky annulus
    files: dict[str, str]           # spectrum, QPB background, ARF, RMF
    counts: int
    area_arcmin2: float


class Archive(Protocol):
    def pointings(self, ra: float, dec: float, max_obs: int = 1) -> list[str]:
        """Public, non-mosaic XMM pointings: the deepest within 3' of the position, else within 10'."""


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


@dataclass
class AutofitResult:
    name: str
    verdict: str                    # ok / check / rejected (``xmmcluster.rules.verdict``)
    kT: float | None                # joint fit in the primary region, 0.15-1 r500 (capped at the FOV)
    e_kT: float | None
    Z: float | None
    cameras: list[str]              # exposures kept by the screening
    flags: list[str]                # every plan and fit flag, with the reason recorded per choice
    r500_kpc: float                 # from the fitted kT (M-T relation), iterated


class Autofit(Protocol):
    def run(self, cluster: ClusterInput, observations: list[Observation]) -> AutofitResult:
        """Regions from r500 (``rules.regions``), recentring on the X-ray centroid, screening of every
        exposure (``rules.screen_spectrum``), bands and N_H from the rules, joint fits of all kept
        cameras with the sky annulus, r500 moved and the fits redone while it changes by more than
        10 % (only from a trusted fit), then the quality flags and the verdict."""


class SkySpectrum(Protocol):
    def rass(self, ra: float, dec: float) -> Spectrum:
        """ROSAT All-Sky Survey spectrum of an annulus around the cluster (HEASARC X-ray background
        tool), fitted with the XMM spectra so the sky components are constrained beyond the FOV."""


class NewSample(Protocol):
    """One command per step from a CSV of clusters to verdicts; every step skips what is done, so an
    interrupted run is restarted as it is. The clusters are held out of the initial-model training."""

    def prepare(self, rows: list[dict]) -> tuple[list[ClusterInput], list[str]]:
        """Check each row (name, position, z), pick the pointings; return the clusters and the problems."""

    def run(self, clusters: list[ClusterInput]) -> None:
        """Ingest -> autofit -> RASS sky spectrum -> refit with it, as job queues on the server."""

    def report(self, clusters: list[ClusterInput]) -> list[AutofitResult]:
        """One line per cluster: verdict, kT and Z, cameras, whether RASS constrained the sky, flags."""
