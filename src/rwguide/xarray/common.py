# TODO: actually use this as a look-up table, at the moment only the first name
# from the list is recognized
NAMES = {
    "lat": ["latitude", "lat", "ylat"],
    "lon": ["longitude", "lon", "xlon"],
    "isob": ["level", "lev", "plev", "isobaricInhPa", "isob"],
    "isen": ["isentrope", "isen"],
    # https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation
    "u": ["u", "U"],
    "v": ["v", "V"],
    "t": ["t", "T"],
    "pres": ["pres", "pressure"],
    "pt": ["pt", "potential_temperature", "th", "theta"],
    "sg": ["sg", "isentropic_density", "sigma"],
    "vo": ["vo", "vorticity", "relative_vorticity"],
    "av": ["av", "avo", "absvort", "absolute_vorticity"],
    "pv": ["pv", "potential_vorticity"],
    "ks": ["ks", "statwaven", "Ks"]
}


def get_names(dct, *names):
    """Standard variable names or an override name if provided.

    Use to support user-supplied variable name overrides in xarray wrappers.

    Parameters
    ----------
    dct : Dict | None
        Variable name overrides, mapping variable shorthands to returned names.
    names : str
        Shorthands of the variables to process.

    Returns
    -------
    Tuple[str, ...]
        Variables names.
    """
    if dct is None:
        dct = dict()
    return tuple(dct.get(name, standard_name(name)) for name in names)


def standard_name(identifier):
    """Standard name corresponding to a variable shorthand of rwguide."""
    return NAMES[identifier][0]


def require_lon(da, lon_coord, axis=-1):
    """Blow meridional profiles up to full lat-lon fields

    .. versionadded:: 1.2
    """
    if lon_coord.name not in da.dims:
        da = da.expand_dims({ lon_coord.name: lon_coord }, axis=axis)
    return da

