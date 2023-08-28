import xarray as xr


def group_count(grouped, grouper):
    # This is simpler and faster than DataArrayGroupBy.count()
    counts = { k: len(g) for k, g in grouped.groups.items() }
    return xr.DataArray(
        data=list(counts.values()),
        dims=[grouper],
        coords={ grouper: list(counts.keys()) },
        name="count"
    )

