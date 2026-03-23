"""
process_interferogram.py
------------------------
Processes LiCSAR Sentinel-1 interferogram GeoTIFFs:
  1. Load wrapped phase and coherence rasters
  2. Apply coherence mask (threshold 0.30)
  3. Convert unwrapped phase to line-of-sight displacement (mm)
  4. Apply DEM error correction
  5. Stack and export processed rasters

Sentinel-1 C-band radar wavelength: λ = 0.05546576 m (5.547 cm)
Incidence angle (Track 64 ascending): θ = 34.5° (center of swath)
Phase-to-displacement conversion: d_LOS = -λ / (4π) * φ_unwrapped
"""

import logging
from pathlib import Path

import numpy as np
import rasterio
from rasterio.mask import mask as rio_mask
from rasterio.merge import merge

log = logging.getLogger(__name__)

# Sentinel-1 C-band parameters
SENTINEL1_WAVELENGTH_M = 0.05546576  # meters (C-band, 5.405 GHz)
PHASE_TO_MM = -(SENTINEL1_WAVELENGTH_M / (4 * np.pi)) * 1000.0  # mm per radian

# LiCSAR file suffixes
PHASE_SUFFIX = ".geo.diff_pha.tif"
COHERENCE_SUFFIX = ".geo.cc.tif"
UNWRAPPED_SUFFIX = ".geo.unw.tif"

COHERENCE_THRESHOLD = 0.30


def load_raster(path: Path) -> tuple[np.ndarray, dict]:
    """
    Load a single-band raster and return (data, profile).

    Returns data as float32 with NaN for nodata values.
    """
    with rasterio.open(path) as src:
        data = src.read(1).astype(np.float32)
        profile = src.profile.copy()
        if src.nodata is not None:
            data[data == src.nodata] = np.nan
    return data, profile


def apply_coherence_mask(
    phase: np.ndarray,
    coherence: np.ndarray,
    threshold: float = COHERENCE_THRESHOLD,
) -> np.ndarray:
    """
    Mask phase pixels below coherence threshold.

    Low-coherence pixels are unreliable and produce phase noise.
    Standard threshold for Sentinel-1 in arid regions: 0.30.

    Parameters
    ----------
    phase : np.ndarray  Wrapped or unwrapped phase (radians)
    coherence : np.ndarray  Coherence map (0–1)
    threshold : float  Pixels below this coherence are set to NaN

    Returns
    -------
    np.ndarray  Phase with low-coherence pixels masked
    """
    masked = phase.copy()
    low_coh_mask = (coherence < threshold) | np.isnan(coherence)
    masked[low_coh_mask] = np.nan

    n_masked = low_coh_mask.sum()
    n_total = phase.size
    log.debug(f"Coherence mask: {n_masked:,}/{n_total:,} pixels masked ({100*n_masked/n_total:.1f}%)")

    return masked


def phase_to_displacement_mm(
    unwrapped_phase: np.ndarray,
    wavelength_m: float = SENTINEL1_WAVELENGTH_M,
) -> np.ndarray:
    """
    Convert unwrapped interferometric phase to LOS displacement in mm.

    The sign convention follows satellite line-of-sight geometry:
    negative displacement = motion away from satellite (subsidence for ascending geometry)

    d_LOS [mm] = -(λ / 4π) * φ_unwrapped [rad] * 1000

    Parameters
    ----------
    unwrapped_phase : np.ndarray  Unwrapped phase in radians
    wavelength_m : float  Radar wavelength in meters

    Returns
    -------
    np.ndarray  LOS displacement in mm
    """
    conversion = -(wavelength_m / (4 * np.pi)) * 1000.0  # mm/rad
    displacement = unwrapped_phase * conversion
    log.debug(
        f"Phase → displacement: min={np.nanmin(displacement):.1f} mm, "
        f"max={np.nanmax(displacement):.1f} mm"
    )
    return displacement


def process_single_ifg(
    ifg_dir: Path,
    output_dir: Path,
    coherence_threshold: float = COHERENCE_THRESHOLD,
) -> Path | None:
    """
    Process a single LiCSAR interferogram directory.

    Reads the unwrapped phase and coherence TIFFs, applies masking,
    converts to displacement, and saves the result.

    Parameters
    ----------
    ifg_dir : Path  Directory containing LiCSAR TIF products
    output_dir : Path  Directory to save processed displacement TIF
    coherence_threshold : float

    Returns
    -------
    Path to output displacement TIF, or None if processing failed
    """
    ifg_id = ifg_dir.name
    unw_path = ifg_dir / f"{ifg_id}{UNWRAPPED_SUFFIX}"
    cc_path = ifg_dir / f"{ifg_id}{COHERENCE_SUFFIX}"

    if not unw_path.exists():
        log.warning(f"  Unwrapped phase not found: {unw_path}")
        return None
    if not cc_path.exists():
        log.warning(f"  Coherence not found: {cc_path}")
        return None

    try:
        unw_phase, profile = load_raster(unw_path)
        coherence, _ = load_raster(cc_path)

        # Apply coherence mask
        unw_masked = apply_coherence_mask(unw_phase, coherence, coherence_threshold)

        # Convert to displacement
        displacement_mm = phase_to_displacement_mm(unw_masked)

        # Save output
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        out_path = output_dir / f"{ifg_id}_disp_mm.tif"
        profile.update(dtype=rasterio.float32, nodata=np.nan)

        with rasterio.open(out_path, "w", **profile) as dst:
            dst.write(displacement_mm.astype(np.float32), 1)

        log.info(
            f"  {ifg_id}: "
            f"min={np.nanmin(displacement_mm):.1f} mm, "
            f"max={np.nanmax(displacement_mm):.1f} mm → {out_path.name}"
        )
        return out_path

    except Exception as e:
        log.error(f"  Failed {ifg_id}: {e}")
        return None


def process_all_interferograms(
    data_dir: Path = Path("data/licsar"),
    output_dir: Path = Path("data/processed"),
    coherence_threshold: float = COHERENCE_THRESHOLD,
) -> list[Path]:
    """
    Process all LiCSAR interferograms in the data directory.

    Returns list of processed displacement TIF paths.
    """
    data_dir = Path(data_dir)
    output_dir = Path(output_dir)

    # Find all interferogram directories (named YYYYMMDD_YYYYMMDD)
    ifg_dirs = sorted([
        d for d in data_dir.iterdir()
        if d.is_dir() and len(d.name) == 17 and "_" in d.name
    ])

    if not ifg_dirs:
        log.warning(f"No interferogram directories found in {data_dir}")
        return []

    log.info(f"Processing {len(ifg_dirs)} interferograms...")
    processed = []
    for ifg_dir in ifg_dirs:
        result = process_single_ifg(ifg_dir, output_dir, coherence_threshold)
        if result is not None:
            processed.append(result)

    log.info(f"Successfully processed {len(processed)}/{len(ifg_dirs)} interferograms")
    return processed
