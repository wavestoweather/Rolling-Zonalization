import cffi
import numpy as np


ffi = cffi.FFI()

DTYPE = np.float64
PTYPE = "double *"

def as_pointer(arr):
    """Return a CFFI-compatible double pointer to the array."""
    if arr is None:
        return ffi.NULL
    assert arr.dtype == DTYPE
    return ffi.cast(PTYPE, ffi.from_buffer(arr))


def empty(*shape, dtype=np.double, order="C"):
    """Wraps `numpy.empty` with default settings."""
    return np.empty(shape, dtype=dtype, order=order)


def _are_compatible_shapes(xshape, yshape):
    return all((xsize == ysize or ysize is None or xsize is None) for xsize, ysize in zip(xshape, yshape))

def require(arr, dtype=np.double, order="C", shape=None, fold=False):
    """Enforce a specification of shape and type for a given array.

    Parameters
    ----------
    arr : array-like
        The input array.
    dtype : data type, optional
        The required data-type. If None preserve the current dtype.
    order : "C" | "F"
        Ensure a C- or Fortran-contiguous array is returned.
    shape : tuple, optional
        The desired output shape. Use `None` to mark a dimension of unknown
        size and integer to specify a dimension of known size (these will be
        enforced/verified). If the input array has fewer dimensions than
        specified in `shape` it is padded with size-1 dimensions from both
        sides until a compatible shape is found. If no compatible shape can be
        found, an exception is thrown.
    fold : boolean, optional
        When folding is enabled and the input array has more dimensions than
        specified in `shape`, additional dimensions of the array are flattened
        into the first dimension.

    Returns
    -------
    out : numpy.array
        The output array matching the specification.
    out_shape : tuple
        The shape of the output array.
    orig_shape : tuple
        The shape of the input array.

    Raises
    ------
    ValueError
        The specification cannot be enforced.
    """
    # Turn the input array into an appropriately-typed numpy array
    out = np.require(arr, dtype=dtype, requirements=order)
    # Save the original shape of the array before reshaping
    orig_shape = out.shape
    # A shape is specified in addition to an already existing array
    if shape is not None:
        ndim = len(shape)
        # If the given array has fewer dimensions than required, pad size-1
        # dimensions. Padding is applied on both sides, from only-right padding
        # to only-left padding until a compatible shape is found. If no
        # compatible shape can be generated through padding, raise an error.
        if out.ndim < ndim:
            mdim = ndim - out.ndim # number of missing dimensions
            for i in range(mdim+1):
                out_shape = (1,)*(mdim-i) + out.shape + (1,)*i
                if _are_compatible_shapes(out_shape, shape):
                    break
            else:
                raise ValueError(f"array of shape {out.shape} incompatible with required shape of {shape}")
            out = out.reshape(out_shape)
        # Array already has the required number of dimensions, just make sure
        # they are compatible
        elif out.ndim == ndim:
            if not _are_compatible_shapes(out.shape, shape):
                raise ValueError(f"array of shape {out.shape} incompatible with required shape of {shape}")
        # Array has more dimensions that required and folding dimensions is
        # not allowed here, so raise an error
        elif fold:
            fix_shape = out.shape[-ndim+1:]
            if not _are_compatible_shapes(fix_shape, shape[1:]):
                raise ValueError(f"array of shape {out.shape} incompatible with required shape of {shape} even when folding")
            fld_shape = out.shape[:-ndim+1]
            out = out.reshape((np.prod(fld_shape), *fix_shape))
        else:
            raise ValueError(f"the array has more dimensions ({out.ndim}) than it should ({ndim})")
    return out, out.shape, orig_shape

