from dataclasses import dataclass 
import datetime
import math
from typing import Optional

def true_to_mean_anomaly(true_anomaly, eccentricity):
    true_anomaly_rad = math.radians(true_anomaly)
    mean_eccentricity = 2 * math.atan2(
        math.sqrt(1 - eccentricity) * math.sin(true_anomaly_rad / 2),
        math.sqrt(1 + eccentricity) * math.cos(true_anomaly_rad / 2)
    )
    mean_anomaly_rad = mean_eccentricity - eccentricity * math.sin(mean_eccentricity)
    return math.degrees(mean_anomaly_rad)

def line_checksum(line: str):
    total = 0
    for c in line[:68]:
        if c.isdigit():
            total += int(c)
        elif c == "-":
            total += 1
    return total % 10

@dataclass
class OrbitalElements:
    epoch: datetime.datetime
    semi_major_axis: float
    eccentricity: float
    inclination: float
    raan: float
    argument_of_perigee: float
    true_anomaly: float
    period: float
    name: Optional[str] = None 
    satid: int = 0

    def to_3le(self):
        if self.name == None:
            raise ValueError("Missing name")
        return f"{self.name.upper()}\n{self.to_tle()}"

    def to_tle(self):
        jd_year = str(self.epoch.year)[-2:]
        jd_day = self.epoch.timetuple().tm_yday + self.epoch.hour/24 + self.epoch.minute / (60 * 24) + self.epoch.second / (24 * 3600)
        line1 = f"1 {self.satid:05}U 00000A   {jd_year}{jd_day:012.8f}  .00000000  00000-0 -00000-0 0    0"
        
        eccentricity = f"{self.eccentricity:.7f}"[2:]
        mean_anomaly = true_to_mean_anomaly(self.true_anomaly, self.eccentricity)
        mean_motion = 24 * 3600 / self.period
        line2 = f"2 {self.satid:05} {self.inclination:08.4f} {self.raan:08.4f} {eccentricity} {self.argument_of_perigee:08.4f} {mean_anomaly:08.4f} {mean_motion:011.8f}00000"

        checksum_1 = line_checksum(line1)
        checksum_2 = line_checksum(line2)
        return f"{line1}{checksum_1}\n{line2}{checksum_2}"

fields = {
    "Spacecraft": ("name", str),
    "UTC time at deployment": ("epoch", lambda d: datetime.datetime.strptime(d, "%Y-%m-%d %H:%M:%S.%fZ")),
    "Semi-major axis (km)": ("semi_major_axis", float),
    "Eccentricity (-)": ("eccentricity", float),
    "Inclination (deg)": ("inclination", float),
    "RAAN (deg)": ("raan", float),
    "Argument of perigee (deg)": ("argument_of_perigee", float),
    "True anomaly (deg)": ("true_anomaly", float),
    "Period (s)": ("period", float),
}

def parse_file(file):
    result = {}
    for line in file.readlines():
        if ":" not in line:
            continue
        key, value = line.split(":", maxsplit=1)
        if key not in fields:
            continue

        field_name, field_parser = fields[key]
        result[field_name] = field_parser(value.strip())

    return OrbitalElements(**result)

with open("orbital_elements_0804.txt") as f:
    elements = parse_file(f)

elements.satid = 99524
elements.name = "FRAMSAT-1"
print(elements.to_3le())