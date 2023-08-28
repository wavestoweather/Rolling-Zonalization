import numpy as np
from matplotlib.path import Path
from matplotlib.patches import FancyArrowPatch, PathPatch
from matplotlib.lines import Line2D
import mpl_toolkits.mplot3d.art3d as art3d
import cartopy.crs as ccrs


PLATE = ccrs.PlateCarree()
ROBIN = ccrs.Robinson()


def add_box(fig, rect, label=None, title=None, **kwargs):
    ax = fig.add_axes(rect, **kwargs)
    if title is None and label is not None:
        fig.text(rect[0], rect[1]+rect[3], f"{label}", ha="left", va="top", weight="bold", fontsize="medium")
    elif title is not None:
        ax.set_title(title, loc="center", fontsize="medium")
        if label is not None:
            ax.set_title(f"{label}", loc="left", weight="bold", fontsize="medium")
    return ax

def add_map(fig, rect, label=None, title=None, coastlines=True, anchor="C", projection=ROBIN):
    ax = add_box(fig, rect, label=label, title=title, anchor=anchor, projection=projection)
    if coastlines:
        ax.coastlines(linewidth=0.5, alpha=0.4)
    return ax


def finalize_map(ax, region=None, aspect=None):
    if region is None:
        ax.set_extent([-179.9, 179.9, -83, 83], crs=PLATE)
    elif region == "NH":
        ax.set_extent([-179.9, 179.9, 9, 83], crs=PLATE)
    elif region == "NH*":
        ax.set_extent([-179.9, 179.9, 9, 73], crs=PLATE)
    elif region == "SH":
        ax.set_extent([-179.9, 179.9, -83, -9], crs=PLATE)
    else:
        raise NotImplementedError(f"region {region}")
    if aspect is not None:
        ax.set_aspect(aspect)


def draw_arrow(fig, start, end):
    # https://stackoverflow.com/a/67531807
    arrow = FancyArrowPatch(start, end, shrinkA=0, shrinkB=0, transform=fig.transFigure,
            color="k", arrowstyle="-|>", mutation_scale=15, linewidth=1.5)
    fig.patches.append(arrow)


def draw_compass_3d(ax, x, y, z, scale=1, color="k", linewidth=0.8):
    ew = ax.add_patch(FancyArrowPatch([x-20*scale, y], [x+20*scale, y], arrowstyle="<|-|>",
            color=color, linewidth=linewidth, mutation_scale=50000, mutation_aspect=0.3))
    art3d.pathpatch_2d_to_3d(ew, z=z)
    ns = ax.add_patch(FancyArrowPatch([x, y-9*scale], [x, y+9*scale], arrowstyle="<|-|>",
            color=color, linewidth=linewidth, mutation_scale=50000, mutation_aspect=0.3))
    art3d.pathpatch_2d_to_3d(ns, z=z)
    ax.text(x-28*scale, y, z, "W", va="center", ha="center", color=color, fontsize="x-small")
    ax.text(x, y+12*scale, z, "N", va="center", ha="center", color=color, fontsize="x-small")
    ax.text(x+28*scale, y, z, "E", va="center", ha="center", color=color, fontsize="x-small")
    ax.text(x, y-12*scale, z, "S", va="center", ha="center", color=color, fontsize="x-small")


def feature_to_patches(feature, intersect_with=None, **kwargs):
    """Plot cartopy features in non-GeoAxes"""
    out = []
    for g in feature.geometries():
        if intersect_with is not None:
            g = g.intersection(intersect_with)
        out.extend(polygon_to_patches(g, **kwargs))
    return out

def polygon_to_patches(polygon, **kwargs):
    """Convert a shapely (multi-)polygon to matplotlib patches"""
    # Recursive descent if multi-polygon
    if polygon.geom_type == "MultiPolygon":
        patches = []
        for p in polygon.geoms:
            patches.extend(polygon_to_patches(p, **kwargs))
        return patches
    # Single polygon:
    if polygon.is_empty:
        return []
    # https://github.com/geopandas/geopandas/issues/1039#issuecomment-748625852
    return [PathPatch(
        Path.make_compound_path(
            Path(np.asarray(polygon.exterior.coords)[:, :2]),
            *[Path(np.asarray(ring.coords)[:, :2]) for ring in polygon.interiors]
        ), **kwargs
    )]


cf_style_baro = {
    "colors": [
        # start............: 2.0
        # rotations........: -1.5
        # hue..............: 1.4
        # gamma............: 1.0
        # number of levels.: 13
        "#FFFFFF", "#FFDFE3", "#EED0AA", "#B5D17F", "#70CF83",
        "#48B7AB", "#5288CC", "#7854BD", "#8F307C", "#7B252F"
    ],
    "extend": "both"
}

# Potential Vorticity
cf_style_pv = {
    "levels": [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5],
    "cmap": "PuBu",
    "extend": "both"
}
# "Waveguidability": grad(log(PV)
cf_style_gl = {
    "levels": np.asarray([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0]) * 1.0e-6,
    "colors": [
        # start............: 2.5
        # rotations........: -2
        # hue..............: 1.0
        # gamma............: 1.0
        # number of levels.: 16
        "#FFFFFF", "#EEEEEE", "#DDDDDD", "#CCCCCC", "#DFAF9A",
        "#AEB272", "#74B471", "#4FA48F", "#5280AA", "#6D57A1"
    ],
    "extend": "both"
}
# Area Weights
cf_style_aw = {
    "levels": np.linspace(0, 1, 11),
    "cmap": "binary"
}
# Waveguide occurrence frequency
cf_style_freq = {
    "levels": [10, 20, 30, 40, 50, 60, 70, 80, 90],
    "colors": [
        # start............: 2.0
        # rotations........: -1.5
        # hue..............: 1.4
        # gamma............: 1.0
        # number of levels.: 13
        "#FFFFFF", "#FFDFE3", "#EED0AA", "#B5D17F", "#70CF83",
        "#48B7AB", "#5288CC", "#7854BD", "#8F307C", "#7B252F"
    ],
    "extend": "both"
}
# Meridional wind
ct_style_v = {
    "levels": [-90, -80, -70, -60, -50, -40, -30, -20, 20, 30, 40, 50, 60, 70, 80, 90],
    "colors": ["#000"]*8 + ["#600"]*8,
    "linestyles": ["dashed"]*8 + ["solid"]*8,
    "linewidths": 0.6
}

def add_legend_v(fig, **kwargs):
    pos = Line2D([], [], linewidth=1, color="#600", linestyle="solid")
    neg = Line2D([], [], linewidth=1, color="#003", linestyle="dashed")
    return fig.legend(
        [neg, pos],
        ["$-20$, $-30$, ...", "$+20$, $+30$, ..."],
        title="mer. wind $v$ [$\mathrm{m}\;\mathrm{s}^{-1}$]",
        **kwargs
    )

