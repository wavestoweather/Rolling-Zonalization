from numbers import Number

import numpy as np 
import scipy.signal as sps


def area_weights(lat, radius=1., nlon=1):
    """Area weights for a lat-lon-grid on the sphere."""
    phi = np.full(len(lat)+1, np.pi/180.) # degree to radians factor
    phi[0   ] *= lat[0]
    phi[1:-1] *= 0.5 * (lat[1:] + lat[:-1])
    phi[  -1] *= lat[-1]
    return (2. * np.pi / nlon) * (radius * radius) * (np.sin(phi[:-1]) - np.sin(phi[1:]))


def earth_circumference(lat):
    """Circumference of the Earth (in m) at the given latitude(s)."""
    return 2. * np.pi * 6371.2e3 * np.cos(np.deg2rad(lat))

def wavelength_to_wavenumber(wavelength, lat):
    """Convert wavelength to wavenumber (latitude-dependent)."""
    return np.abs(earth_circumference(lat) / wavelength)

def wavenumber_to_width(wavenumber, n):
    """Convert wavenumber to width (length or # of gridpoints)."""
    return n // wavenumber

def wavelength_to_width(wavelength, lat, n):
    """Convert wavelength to width (latitude-dependent)."""
    return wavenumber_to_width(wavelength_to_wavenumber(wavelength, lat), n)

def oddify(width):
    """Return the next larger or equal odd number."""
    width = np.asarray(width).astype(int)
    return width + (width + 1) % 2 # next larger odd number for even

def clip_width(width, nmin=3, nmax=3600):
    """Convert to integer (# of gridpoints) and clip to the given range."""
    width = np.asarray(width).astype(int)
    return np.clip(width, nmin, nmax)


def normalize(window):
    """Normalize such that window weights add to 1."""
    window = np.asarray(window)
    # In the context of windowed zonalization, negative values should not
    # appear in a weighting function. Fail visibly after a check here.
    if np.any(window < 0):
        raise ValueError("negative value(s) in weighting function")
    # Weights have to add up to 1
    return window / window.sum()

def fold_periodic(window, n, agg="max"):
    """Fold the window into the periodic domain.

    If the window fits into the domain, it is returned unchanged. Otherwise the
    window is centered in the domain (at ``n//2``). Example::

        n = 14, window.size = 18

        |<-window.size-->|        |<-window.size-->|
        |<----n----->|   |        | |<----n----->| |
        |     ..··.. |   |     →  | |   ..··..   | |
        ...···      ···...        ...···   |  ···...
        |------------|----        --|------|-----|--
        012345678901234567          01234567890123
                                           ^ n//2

    Parameters
    ----------
    window : np.ndarray
        The window to fold.
    n : number
        The width of the periodic domain (# of gridpoints).
    agg : "sum" | "clip" | "max", optional
        How to aggregate values in the folded regions. `clip` uses sum but
        restricts values to the maximum of the original window after summation.
    Returns
    -------
    np.ndarray
    """
    # No folding needed when window is smaller than periodic domain
    if n <= 0 or n > window.size:
        return window
    # Extract the maximum value of the window
    maxval = np.max(window)
    # Cut into parts that are of the same size as the domain and collect
    parts = []
    i = 0
    while i < window.size:
        part = window[i:(i+n)]
        part = np.append(part, np.zeros(n - part.size))
        parts.append(part)
        i += n
    # Aggregate the weights of the overlapping parts with the chosen method:
    # clip (sum and clip at previous maximum), sum or max (equivalent to
    # cropping the window for all windows with monotonic lobes, default).
    if agg == "clip":
        folded = np.clip(np.sum(parts, axis=0), 0., maxval)
    elif agg == "sum":
        folded = np.sum(parts, axis=0)
    elif agg == "max":
        folded = np.max(parts, axis=0)
    else:
        raise NotImplementedError(f"unknown aggregation option {agg}")
    # Move the weights such that the original center is at n//2.
    distance = (n // 2) - ((window.size // 2) % n)
    return np.roll(folded, distance)


def get_window(name, width, fold=None):
    """Create a weighting window for the sphere.

    Parameters
    ----------
    name : string | tuple
        Which window function to use. Choose a name from
        :py:mod:`scipy.signal.windows` or specify a tuple containing the name
        and additional arguments for the window creation function. Examples:
        ``"boxcar"``, ``("tukey", 0.5)``.
    width : number | iterable
        Window width in # of gridpoints. Constant or latitude-dependent.
    fold : None | Tuple | number, optional
        Settings for folding the window into the periodic domain. Provide the
        `n` and `agg` arguments for :py:func:`fold_periodic`.

    Returns
    -------
    numpy.ndarray
    """
    if isinstance(width, Number):
        name, *args = (name,) if isinstance(name, str) else name
        out = getattr(sps.windows, name)(width, *args)
        if fold is not None:
            out = fold_periodic(out, *fold)
    else:
        # Stack multiple windows if multiple widths are given
        windows = [get_window(name, w, fold=fold) for w in width]
        ny = len(windows)
        nx = max(map(len, windows))
        out = np.zeros((ny, nx))
        for i, window in enumerate(windows):
            l = (nx - window.size) // 2
            r = l + window.size
            out[i,l:r] = window / np.max(window)
    return out

def product(window_lon, window_lat):
    """2D product of a (1D/2D) longitude and a (1D/2D) latitude window."""
    if window_lon.ndim == 1:
        window_lon = window_lon[None,:]
    if window_lat.ndim == 1:
        window_lat = window_lat[:,None]
    out = window_lon * window_lat
    return out

