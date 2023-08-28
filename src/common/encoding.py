def encoding_int16(scale_factor, add_offset, zlib=False):
    return {
        "dtype": "int16",
        "scale_factor": scale_factor,
        "add_offset": add_offset,
        "_FillValue": -32768,
        "zlib": zlib
    }

def encoding_int16_range(lower, upper, **kwargs):
    scale = (upper - lower) / (2**16 - 2)
    # translate the range to be symmetric about zero
    offset = lower + (2**15 - 1) * scale
    return encoding_int16(scale, offset)


class preset:
    """Encodings with ranges appropriate for most tropospheric conditions"""

    # Zonal wind [m/s]
    u = encoding_int16_range(-80, 130)
    # Meridional wind [m/s]
    v = encoding_int16_range(-110, 110)
    # Temperature [K]
    t = encoding_int16_range(150, 320)
    # Pressure [hPa]
    pres = encoding_int16_range(0, 1100)
    # Potential vorticity [PVU]
    pv = encoding_int16_range(-25, 25)
    # norm grad log PV gradient
    pv_grad_log = encoding_int16_range(0, 5.0e-6)


_keys = ("dtype", "_FillValue", "scale_factor", "add_offset", "zlib")

def extract_encoding(ds, keys=_keys, skip_coords=True):
    """Copy encoding information from an xarray dataset"""
    encoding = dict()
    for x in ds.variables:
        if skip_coords and x in ds.coords:
            continue
        if hasattr(ds[x], "encoding"):
            encoding[x] = { key: ds[x].encoding[key] for key in keys if key in ds[x].encoding }
    return encoding
