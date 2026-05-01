#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Orbit visualizer and pass predictor for satellites.
Reads orbital elements from a text file and plots ground track,
finds passes over a given ground station.

Usage examples:
    python tle_visualizer.py orbital_elements_0804.txt --plot-orbits 3 --plot-step-s 15
    python tle_visualizer.py orbital_elements_2503.txt --plot-orbits 1 --plot-step-s 15 --station-lat 69.29 --station-lon 16.02
    python tle_visualizer.py orbital_elements_test.txt --plot-orbits 2 --plot-step-s 30 --station-lat 63.41866 --station-lon 10.3951 --pass-orbits 3
"""
#cd C:\Users\hemin\OneDrive\Dokumenter\Verv\Orbit\ISAR\TLE
import argparse
import datetime as dt
import math
import sys
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt

# Optional dependencies
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAVE_CARTOPY = True
except ImportError:
    HAVE_CARTOPY = False

try:
    from zoneinfo import ZoneInfo
    OSLO_TZ = ZoneInfo("Europe/Oslo")
except ImportError:
    OSLO_TZ = None

# ============================================================
# KONSTANTER
# ============================================================
MU_EARTH_KM3_S2 = 398600.4418          # Gravitational parameter, km^3/s^2
WGS84_A_KM = 6378.137                 # Earth semi-major axis, km
WGS84_F = 1.0 / 298.257223563         # Flattening
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)  # First eccentricity squared


# ============================================================
# GEODETISK TRANSFORMASJON
# ============================================================
def ecef_to_geodetic(x: float, y: float, z: float) -> Tuple[float, float]:
    """
    Convert ECEF coordinates (km) to geodetic latitude (deg) and longitude (deg).
    Uses iterative algorithm (Bowring).
    """
    lon = math.degrees(math.atan2(y, x))
    p = math.sqrt(x * x + y * y)
    if p < 1e-12:
        lat = 90.0 if z >= 0.0 else -90.0
        return lat, lon

    lat = math.atan2(z, p * (1.0 - WGS84_E2))
    for _ in range(10):
        sin_lat = math.sin(lat)
        N = WGS84_A_KM / math.sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat)
        lat_new = math.atan2(z + WGS84_E2 * N * sin_lat, p)
        if abs(lat_new - lat) < 1e-13:
            lat = lat_new
            break
        lat = lat_new
    return math.degrees(lat), lon


# ============================================================
# HJELPEFUNKSJONER FOR TID OG ROTASJON
# ============================================================
def normalize_angle_deg(angle_deg: float) -> float:
    """Normalize angle to [0, 360) degrees."""
    return angle_deg % 360.0


def true_to_mean_anomaly(true_anomaly_deg: float, eccentricity: float) -> float:
    """Convert true anomaly (deg) to mean anomaly (deg)."""
    nu = math.radians(true_anomaly_deg)
    e = eccentricity
    E = 2.0 * math.atan2(
        math.sqrt(1.0 - e) * math.sin(nu / 2.0),
        math.sqrt(1.0 + e) * math.cos(nu / 2.0)
    )
    M = E - e * math.sin(E)
    return normalize_angle_deg(math.degrees(M))


def solve_kepler(M: np.ndarray, e: float, tol: float = 1e-12, max_iter: int = 50) -> np.ndarray:
    """
    Solve Kepler's equation M = E - e*sin(E) for eccentric anomaly E (rad).
    M is array of mean anomalies in radians.
    """
    M = np.asarray(M, dtype=float)
    # Initial guess
    if e < 0.8:
        E = M.copy()
    else:
        E = np.full_like(M, math.pi)
    
    for iteration in range(max_iter):
        f = E - e * np.sin(E) - M
        fp = 1.0 - e * np.cos(E)
        dE = -f / fp
        E = E + dE
        if np.max(np.abs(dE)) < tol:
            break
    else:
        # Warn if not converged
        max_error = np.max(np.abs(E - e * np.sin(E) - M))
        if max_error > 1e-6:
            print(f"Warning: Kepler solver did not converge. Max error: {max_error}")
    
    return E


def r1(theta: float) -> np.ndarray:
    """Rotation matrix about x-axis."""
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([
        [1.0, 0.0, 0.0],
        [0.0, c,   -s],
        [0.0, s,  c]
    ], dtype=float)


def r3(theta: float) -> np.ndarray:
    """Rotation matrix about z-axis."""
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([
        [c,   -s,  0.0],
        [s,  c,  0.0],
        [0.0, 0.0, 1.0]
    ], dtype=float)


def coe_to_eci_matrix(raan_rad: float, inc_rad: float, argp_rad: float) -> np.ndarray:
    """Transformation from perifocal to ECI coordinates."""
    return r3(raan_rad) @ r1(inc_rad) @ r3(argp_rad)



def datetime_to_julian_date(d: dt.datetime) -> float:
    """Convert UTC datetime to Julian Date."""
    if d.tzinfo is None:
        d = d.replace(tzinfo=dt.timezone.utc)
    else:
        d = d.astimezone(dt.timezone.utc)

    year = d.year
    month = d.month
    day = d.day
    hour = d.hour
    minute = d.minute
    second = d.second + d.microsecond / 1e6

    if month <= 2:
        year -= 1
        month += 12

    A = year // 100
    B = 2 - A + (A // 4)

    day_fraction = (hour + minute / 60.0 + second / 3600.0) / 24.0

    jd = (math.floor(365.25 * (year + 4716)) +
          math.floor(30.6001 * (month + 1)) +
          day + day_fraction + B - 1524.5)
    return jd


def gmst_rad(d: dt.datetime) -> float:
    """Greenwich Mean Sidereal Time in radians at given UTC datetime."""
    jd = datetime_to_julian_date(d)
    jd_ut1 = jd  # Simplification: assuming UT1 ≈ UTC
    T = (jd_ut1 - 2451545.0) / 36525.0
    
    # GMST at 0h UT
    gmst_0h = 280.46061837 + 360.98564736629 * (jd_ut1 - 2451545.0) + \
              0.000387933 * T**2 - T**3 / 38710000.0
    
    # Normalize to [0, 360)
    gmst_deg = gmst_0h % 360.0
    return math.radians(gmst_deg)


def parse_utc_datetime(value: str) -> dt.datetime:
    """Parse ISO datetime string to UTC datetime."""
    value = value.strip()
    if value.endswith('Z'):
        value = value[:-1] + '+00:00'
    parsed = dt.datetime.fromisoformat(value)
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=dt.timezone.utc)
    return parsed.astimezone(dt.timezone.utc)


def split_on_dateline(lons: np.ndarray, lats: np.ndarray, threshold: float = 180.0) -> Tuple[np.ndarray, np.ndarray]:
    """Insert NaN where longitude jumps across ±180°."""
    lons = np.asarray(lons, dtype=float)
    lats = np.asarray(lats, dtype=float)
    if len(lons) == 0:
        return lons, lats

    out_lons = [lons[0]]
    out_lats = [lats[0]]
    for i in range(1, len(lons)):
        diff = abs(lons[i] - lons[i-1])
        if diff > threshold:  # Using parameter instead of hardcoded 180
            out_lons.append(np.nan)
            out_lats.append(np.nan)
        out_lons.append(lons[i])
        out_lats.append(lats[i])
    return np.array(out_lons), np.array(out_lats)


def format_dt(d: dt.datetime) -> str:
    """Format datetime to UTC and local (Oslo) time."""
    utc = d.astimezone(dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S.%f")[:-3] + " UTC"
    if OSLO_TZ is None:
        return utc
    local = d.astimezone(OSLO_TZ).strftime("%Y-%m-%d %H:%M:%S.%f")[:-3] + " Europe/Oslo"
    return f"{utc} / {local}"


# ============================================================
# DATACLASSER
# ============================================================
@dataclass
class GroundStation:
    name: str
    lat_deg: float
    lon_deg: float
    alt_m: float = 0.0

    def ecef_km(self) -> np.ndarray:
        """ECEF coordinates (km) of the station."""
        lat = math.radians(self.lat_deg)
        lon = math.radians(self.lon_deg)
        alt_km = self.alt_m / 1000.0
        sin_lat = math.sin(lat)
        cos_lat = math.cos(lat)
        N = WGS84_A_KM / math.sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat)
        x = (N + alt_km) * cos_lat * math.cos(lon)
        y = (N + alt_km) * cos_lat * math.sin(lon)
        z = (N * (1.0 - WGS84_E2) + alt_km) * sin_lat
        return np.array([x, y, z], dtype=float)

    def enu_basis(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return east, north, up unit vectors in ECEF coordinates."""
        lat = math.radians(self.lat_deg)
        lon = math.radians(self.lon_deg)
        east = np.array([-math.sin(lon), math.cos(lon), 0.0], dtype=float)
        north = np.array([
            -math.sin(lat) * math.cos(lon),
            -math.sin(lat) * math.sin(lon),
            math.cos(lat)
        ], dtype=float)
        up = np.array([
            math.cos(lat) * math.cos(lon),
            math.cos(lat) * math.sin(lon),
            math.sin(lat)
        ], dtype=float)
        return east, north, up


@dataclass
class VisibilityPass:
    rise_utc: dt.datetime
    set_utc: dt.datetime
    tca_utc: dt.datetime
    max_elevation_deg: float
    slant_range_km: float

    @property
    def duration_s(self) -> float:
        return (self.set_utc - self.rise_utc).total_seconds()


@dataclass
class OrbitalElements:
    epoch: dt.datetime
    semi_major_axis_km: float
    eccentricity: float
    inclination_deg: float
    raan_deg: float
    argument_of_perigee_deg: float
    true_anomaly_deg: float
    period_s: float
    name: Optional[str] = None
    satid: int = 0

    @property
    def mean_motion_rev_per_day(self) -> float:
        return 86400.0 / self.period_s

    @property
    def mean_anomaly_deg(self) -> float:
        return true_to_mean_anomaly(self.true_anomaly_deg, self.eccentricity)

    def propagate_eci(self, times_s: np.ndarray) -> np.ndarray:
        """
        Propagate satellite in ECI frame (km) at given times (seconds from epoch).
        Returns array shape (N, 3).
        """
        a = self.semi_major_axis_km
        e = self.eccentricity
        inc = math.radians(self.inclination_deg)
        raan = math.radians(self.raan_deg)
        argp = math.radians(self.argument_of_perigee_deg)

        M0 = math.radians(self.mean_anomaly_deg)
        n = math.sqrt(MU_EARTH_KM3_S2 / (a ** 3))  # rad/s
        M = np.mod(M0 + n * times_s, 2.0 * math.pi)
        E = solve_kepler(M, e)

        x_pf = a * (np.cos(E) - e)
        y_pf = a * (np.sqrt(1.0 - e * e) * np.sin(E))
        z_pf = np.zeros_like(x_pf)

        r_pf = np.vstack((x_pf, y_pf, z_pf))  # 3 x N
        Q = coe_to_eci_matrix(raan, inc, argp)
        r_eci = (Q @ r_pf).T  # N x 3
        return r_eci

    def propagate_ecef(self, times_s: np.ndarray) -> np.ndarray:
        """Propagate in ECEF frame (km) at given times (seconds from epoch)."""
        r_eci = self.propagate_eci(times_s)
        r_ecef = np.empty_like(r_eci)
        for i, t in enumerate(times_s):
            theta = gmst_rad(self.epoch + dt.timedelta(seconds=float(t)))
            r_ecef[i] = r3(-theta) @ r_eci[i]
        return r_ecef

    def ground_track(self, duration_s: Optional[float] = None,
                     orbits: Optional[float] = None,
                     step_s: float = 20.0,
                     start_s: float = 0.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute ground track.
        Returns (times_s, lats_deg, lons_deg) arrays.
        """
        if duration_s is None:
            if orbits is not None:
                duration_s = orbits * self.period_s
            else:
                duration_s = self.period_s
        if duration_s <= 0:
            raise ValueError("duration_s must be > 0")
        if step_s <= 0:
            raise ValueError("step_s must be > 0")

        times = np.arange(start_s, start_s + duration_s + step_s, step_s)
        r_ecef = self.propagate_ecef(times)

        lats = np.empty(len(times), dtype=float)
        lons = np.empty(len(times), dtype=float)
        for i, (x, y, z) in enumerate(r_ecef):
            lat, lon = ecef_to_geodetic(float(x), float(y), float(z))
            lats[i] = lat
            lons[i] = lon
        return times, lats, lons

    def topocentric_metrics(self, t_s: float, station: GroundStation) -> Tuple[float, float]:
        """Return elevation (deg) and slant range (km) at time t_s."""
        r_sat = self.propagate_ecef(np.array([t_s], dtype=float))[0]
        r_sta = station.ecef_km()
        rho = r_sat - r_sta
        east, north, up = station.enu_basis()
        e = float(np.dot(rho, east))
        n = float(np.dot(rho, north))
        u = float(np.dot(rho, up))
        slant = float(np.linalg.norm(rho))
        elev = math.degrees(math.atan2(u, math.sqrt(e*e + n*n)))
        return elev, slant

    def elevation_deg_at_time(self, t_s: float, station: GroundStation) -> float:
        elev, _ = self.topocentric_metrics(t_s, station)
        return elev

    def _refine_crossing(self, t1_s: float, t2_s: float, station: GroundStation,
                         target_elev_deg: float, tol_s: float = 0.01) -> float:
        """Bisection to find time when elevation = target."""
        f1 = self.elevation_deg_at_time(t1_s, station) - target_elev_deg
        f2 = self.elevation_deg_at_time(t2_s, station) - target_elev_deg
        if f1 == 0.0:
            return t1_s
        if f2 == 0.0:
            return t2_s
        lo, hi = t1_s, t2_s
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            fm = self.elevation_deg_at_time(mid, station) - target_elev_deg
            if abs(hi - lo) < tol_s:
                return mid
            if (f1 <= 0.0 and fm >= 0.0) or (f1 >= 0.0 and fm <= 0.0):
                hi = mid
                f2 = fm
            else:
                lo = mid
                f1 = fm
        return 0.5 * (lo + hi)

    def find_passes(self, station: GroundStation,
                    duration_s: Optional[float] = None,
                    orbits: Optional[float] = None,
                    step_s: float = 10.0,
                    min_elev_deg: float = 5.0,
                    start_s: float = 0.0) -> List[VisibilityPass]:
        """Find all passes above min_elev_deg."""
        if duration_s is None:
            if orbits is not None:
                duration_s = orbits * self.period_s
            else:
                duration_s = 24.0 * 3600.0
        if duration_s <= 0:
            raise ValueError("duration_s must be > 0")
        if step_s <= 0:
            raise ValueError("step_s must be > 0")

        times = np.arange(start_s, start_s + duration_s + step_s, step_s)
        elevs = np.empty(len(times), dtype=float)
        ranges = np.empty(len(times), dtype=float)
        for i, t in enumerate(times):
            elevs[i], ranges[i] = self.topocentric_metrics(float(t), station)

        above = elevs >= min_elev_deg
        passes = []
        i = 0
        while i < len(times):
            if not above[i]:
                i += 1
                continue
            seg_start = i
            while i+1 < len(times) and above[i+1]:
                i += 1
            seg_end = i

            # refine rise time
            if seg_start == 0:
                rise_s = times[0]
            else:
                rise_s = self._refine_crossing(times[seg_start-1], times[seg_start],
                                               station, min_elev_deg)
            # refine set time
            if seg_end == len(times)-1:
                set_s = times[-1]
            else:
                set_s = self._refine_crossing(times[seg_end], times[seg_end+1],
                                              station, min_elev_deg)
            # TCA at max elevation within segment
            seg_slice = slice(seg_start, seg_end+1)
            local_idx = seg_start + int(np.argmax(elevs[seg_slice]))
            tca_s = times[local_idx]
            max_el = elevs[local_idx]
            slant = ranges[local_idx]

            passes.append(VisibilityPass(
                rise_utc=self.epoch + dt.timedelta(seconds=rise_s),
                set_utc=self.epoch + dt.timedelta(seconds=set_s),
                tca_utc=self.epoch + dt.timedelta(seconds=tca_s),
                max_elevation_deg=max_el,
                slant_range_km=slant
            ))
            i = seg_end + 1
        return passes

    def plot_ground_track_map(self, duration_s: Optional[float] = None,
                              orbits: Optional[float] = None,
                              step_s: float = 20.0,
                              annotate_every_n: int = 0,
                              save_path: Optional[str] = None,
                              show: bool = True,
                              station: Optional[GroundStation] = None) -> None:
        """Plot ground track on a map."""
        _, lats, lons = self.ground_track(duration_s, orbits, step_s)
        lons_plot, lats_plot = split_on_dateline(lons, lats)

        title_name = self.name if self.name else "Satellite"
        title = f"Ground track for {title_name}"

        if HAVE_CARTOPY:
            proj = ccrs.PlateCarree()
            fig = plt.figure(figsize=(16, 8))
            ax = plt.axes(projection=proj)
            ax.set_global()
            ax.add_feature(cfeature.OCEAN, facecolor="#dceaf7", alpha=0.95)
            ax.add_feature(cfeature.LAND, facecolor="#e6e6e6", alpha=0.85)
            ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
            ax.add_feature(cfeature.BORDERS, linewidth=0.4, alpha=0.6)
            gl = ax.gridlines(draw_labels=True, linewidth=0.4, color="gray", alpha=0.5, linestyle="--")
            gl.top_labels = False
            gl.right_labels = False
            ax.plot(lons_plot, lats_plot, transform=proj, linewidth=1.8, label="Ground track")
            ax.scatter([lons[0]], [lats[0]], transform=proj, marker="o", s=50, label="Start")
            ax.scatter([lons[-1]], [lats[-1]], transform=proj, marker="x", s=55, label="End")
            if station:
                ax.scatter([station.lon_deg], [station.lat_deg], transform=proj,
                           marker="^", s=90, label=station.name)
            if annotate_every_n and annotate_every_n > 0:
                times_s, lats_raw, lons_raw = self.ground_track(duration_s, orbits, step_s)
                for i in range(0, len(times_s), annotate_every_n):
                    t_utc = self.epoch + dt.timedelta(seconds=float(times_s[i]))
                    label = t_utc.strftime("%H:%M")
                    ax.annotate(label, (lons_raw[i], lats_raw[i]), transform=proj,
                                fontsize=8, xytext=(4,4), textcoords="offset points")
            ax.set_title(title)
            ax.legend(loc="lower left")
        else:
            print("Warning: Cartopy not installed. Plotting without map projection.")
            fig, ax = plt.subplots(figsize=(12, 6))
            ax.plot(lons_plot, lats_plot, linewidth=1.8, label="Ground track")
            ax.scatter(lons[0], lats[0], marker="o", s=50, label="Start")
            ax.scatter(lons[-1], lats[-1], marker="x", s=55, label="End")
            if station:
                ax.scatter(station.lon_deg, station.lat_deg, marker="^", s=90, label=station.name)
            if annotate_every_n and annotate_every_n > 0:
                times_s, lats_raw, lons_raw = self.ground_track(duration_s, orbits, step_s)
                for i in range(0, len(times_s), annotate_every_n):
                    t_utc = self.epoch + dt.timedelta(seconds=float(times_s[i]))
                    label = t_utc.strftime("%H:%M")
                    ax.annotate(label, (lons_raw[i], lats_raw[i]), xytext=(4,4),
                                textcoords="offset points", fontsize=8)
            ax.set_xlabel("Longitude (deg)")
            ax.set_ylabel("Latitude (deg)")
            ax.grid(True, linestyle="--", alpha=0.7)
            ax.set_title(title)
            ax.legend()
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=220, bbox_inches="tight")
        if show:
            plt.show()
        else:
            plt.close(fig)


# ============================================================
# PARSING AV ELEMENTFIL
# ============================================================
FIELDS = {
    "Spacecraft": ("name", str),
    "UTC time at deployment": ("epoch", parse_utc_datetime),
    "Semi-major axis (km)": ("semi_major_axis_km", float),
    "Eccentricity (-)": ("eccentricity", float),
    "Inclination (deg)": ("inclination_deg", float),
    "RAAN (deg)": ("raan_deg", float),
    "Argument of perigee (deg)": ("argument_of_perigee_deg", float),
    "True anomaly (deg)": ("true_anomaly_deg", float),
    "Period (s)": ("period_s", float),
}

def parse_file(file_obj) -> OrbitalElements:
    """Parse orbital elements from text file."""
    result = {}
    for line in file_obj:
        line = line.strip()
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip()
        if key not in FIELDS:
            continue
        field_name, parser = FIELDS[key]
        result[field_name] = parser(value)
    required = {"epoch", "semi_major_axis_km", "eccentricity", "inclination_deg",
                "raan_deg", "argument_of_perigee_deg", "true_anomaly_deg", "period_s"}
    missing = required - set(result.keys())
    if missing:
        raise ValueError(f"Missing required fields: {sorted(missing)}")
    # Basic validation
    if result["semi_major_axis_km"] <= 0:
        raise ValueError("Semi-major axis must be positive")
    if result["eccentricity"] < 0 or result["eccentricity"] >= 1:
        raise ValueError("Eccentricity must be in [0,1)")
    if result["period_s"] <= 0:
        raise ValueError("Period must be positive")
    return OrbitalElements(**result)


# ============================================================
# UTSKRIFT AV PASS
# ============================================================
def print_passes(passes: List[VisibilityPass], station: GroundStation, min_elev_deg: float) -> None:
    print(f"\nPass over {station.name} with elevation mask {min_elev_deg:.1f}°")
    if not passes:
        print("No passes found in the search window.")
        return
    for idx, p in enumerate(passes, 1):
        print(f"\nPass {idx}")
        print(f"  Rise : {format_dt(p.rise_utc)}")
        print(f"  TCA  : {format_dt(p.tca_utc)}")
        print(f"  Set  : {format_dt(p.set_utc)}")
        print(f"  Max elevation: {p.max_elevation_deg:.2f}°")
        print(f"  Slant range : {p.slant_range_km:.1f} km")
        print(f"  Duration     : {p.duration_s:.1f} s")


# ============================================================
# MAIN
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description="Satellite orbit visualizer and pass predictor",
        epilog="Example: python tle_visualizer.py data.txt --plot-orbits 3 --station-lat 69.29 --station-lon 16.02"
    )
    parser.add_argument("input_file", help="Path to orbital elements file")
    parser.add_argument("--plot-orbits", type=float, default=3.0,
                        help="Number of orbits to plot (default: 3)")
    parser.add_argument("--plot-duration-s", type=float, default=None,
                        help="Alternative to --plot-orbits: duration in seconds")
    parser.add_argument("--plot-step-s", type=float, default=20.0,
                        help="Time step for ground track (seconds, default: 20)")
    parser.add_argument("--annotate-every-n", type=int, default=0,
                        help="Annotate every N points on map (0 = off)")
    parser.add_argument("--save-plot", type=str, default=None,
                        help="Save plot to file (e.g., track.png)")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip ground track plot")
    # Pass search arguments
    parser.add_argument("--pass-orbits", type=float, default=16.0,
                        help="Number of orbits to search for passes (default: 16)")
    parser.add_argument("--pass-duration-s", type=float, default=None,
                        help="Alternative to --pass-orbits: duration in seconds")
    parser.add_argument("--pass-step-s", type=float, default=10.0,
                        help="Time step for pass search (seconds, default: 10)")
    parser.add_argument("--elevation-mask-deg", type=float, default=5.0,
                        help="Minimum elevation for passes (deg, default: 5)")
    parser.add_argument("--start-time-offset", type=float, default=0.0,
                        help="Start time offset from epoch (seconds, default: 0)")
    # Station arguments
    parser.add_argument("--station-name", type=str, default="Trondheim",
                        help="Ground station name (default: Trondheim)")
    parser.add_argument("--station-lat", type=float, default=63.41866,
                        help="Station latitude (deg, default: 63.4305)")
    parser.add_argument("--station-lon", type=float, default=10.3951,
                        help="Station longitude (deg, default: 10.3951)")
    parser.add_argument("--station-alt", type=float, default=51.3,
                        help="Station altitude (m, default: 10)")
    parser.add_argument("--name-override", type=str, default=None,
                        help="Override satellite name in plot")

    args = parser.parse_args()

    # Read file
    try:
        with open(args.input_file, "r", encoding="utf-8") as f:
            elements = parse_file(f)
    except FileNotFoundError:
        print(f"Error: File '{args.input_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing file: {e}")
        sys.exit(1)

    if args.name_override:
        elements.name = args.name_override

    # Create ground station
    station = GroundStation(
        name=args.station_name,
        lat_deg=args.station_lat,
        lon_deg=args.station_lon,
        alt_m=args.station_alt
    )

    # Print summary
    print(f"Satellite: {elements.name if elements.name else 'Unnamed'}")
    print(f"Epoch: {format_dt(elements.epoch)}")
    print(f"Semi-major axis: {elements.semi_major_axis_km:.4f} km")
    print(f"Eccentricity: {elements.eccentricity:.7f}")
    print(f"Inclination: {elements.inclination_deg:.4f} deg")
    print(f"RAAN: {elements.raan_deg:.4f} deg")
    print(f"Argument of perigee: {elements.argument_of_perigee_deg:.4f} deg")
    print(f"True anomaly: {elements.true_anomaly_deg:.4f} deg")
    print(f"Period: {elements.period_s:.3f} s")
    print(f"Mean motion: {elements.mean_motion_rev_per_day:.6f} rev/day")
    print(f"Mean anomaly: {elements.mean_anomaly_deg:.6f} deg")
    print(f"Ground station: {station.name} ({station.lat_deg:.4f}°, {station.lon_deg:.4f}°, {station.alt_m:.1f} m)")

    # Plot ground track
    if not args.no_plot:
        try:
            elements.plot_ground_track_map(
                duration_s=args.plot_duration_s,
                orbits=None if args.plot_duration_s is not None else args.plot_orbits,
                step_s=args.plot_step_s,
                annotate_every_n=args.annotate_every_n,
                save_path=args.save_plot,
                show=True,
                station=station
            )
        except Exception as e:
            print(f"Plotting failed: {e}")
            if not HAVE_CARTOPY:
                print("Consider installing cartopy for better maps: conda install -c conda-forge cartopy")

    # Find and print passes
    pass_duration = args.pass_duration_s
    if pass_duration is None:
        pass_duration = args.pass_orbits * elements.period_s

    try:
        passes = elements.find_passes(
            station=station,
            duration_s=pass_duration,
            step_s=args.pass_step_s,
            min_elev_deg=args.elevation_mask_deg,
            start_s=args.start_time_offset
        )
        print_passes(passes, station, args.elevation_mask_deg)
    except Exception as e:
        print(f"Pass search failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()