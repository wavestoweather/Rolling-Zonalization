import xarray as xr
 

def merge_hemispheres(ds, nh_season, sh_season):
    nh = ds.sel({ "latitude": slice(90,  0), "season": nh_season }, drop=True)
    sh = ds.sel({ "latitude": slice(0, -90), "season": sh_season }, drop=True)
    return xr.concat([nh, sh], dim="latitude")

def select_season(ds, season):
    # Stitch together winter/summer hemispheres if needed
    if season == "winter":
        return merge_hemispheres(ds, "DJF", "JJA")
    elif season == "summer":
        return merge_hemispheres(ds, "JJA", "DJF")
    else:
        return ds.sel({ "season": season })

