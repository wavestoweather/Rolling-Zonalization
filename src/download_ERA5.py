import argparse
import os
import pprint
import cdsapi


request_base = {
    "product_type": "reanalysis",
    "grid": "1.5/1.5",
    "area": "90/-180/-90/180", # central longitude: 0°
    "month": ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"],
    "day": [
        "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12",
        "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24",
        "25", "26", "27", "28", "29", "30", "31"
    ],
    "time": [ "00:00", "06:00", "12:00", "18:00" ],
    "format": "netcdf",
}


parser = argparse.ArgumentParser()
parser.add_argument("year", type=int)
parser.add_argument("--dry-run", action="store_true")
parser.add_argument("--out-dir", default="data/ERA5")

if __name__ == "__main__":
    args = parser.parse_args()
    assert args.year >= 1979

    origin = "reanalysis-era5-pressure-levels"

    request = request_base.copy()
    request["variable"] = ["u_component_of_wind", "v_component_of_wind", "temperature"]
    request["pressure_level"] = [
         "50",  "70", "100", "150", "200", "250", "300", "350", "400", "450",
        "500", "550", "600", "650", "700", "750", "800", "850"
    ]
    request["year"] = str(args.year)

    name = os.path.join(args.out_dir, f"ERA5-{args.year}-tuv-1.5.nc")

    if args.dry_run:
        pprint.pprint(request, compact=True)
    else:
        os.makedirs(args.out_dir, exist_ok=True)
        c = cdsapi.Client()
        c.retrieve(origin, request, name)

