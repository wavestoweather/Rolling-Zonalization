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
    "pv": ["pv", "potential_vorticity"]
}

def get_names(dct, *names):
    if dct is None:
        dct = dict()
    return tuple(dct.get(name, standard_name(name)) for name in names)

def standard_name(identifier):
    return NAMES[identifier][0]

