from __future__ import annotations

from typing import List

from numpy.typing import NDArray

def clusters_from_sparse(indices: NDArray[int]) -> List[NDArray[int]]:
    """
    Find clusters of connected voxel indices from a sparse array of voxel coordinates. The input array should be a Nx3
    array of integer voxel indices. The output is a list of arrays, where each array contains the indices of a single
    cluster of connected voxels.

    The connectivity is defined as 26-connectivity, i.e. each voxel is connected to all neighbors with which it shares
    at least one corner.

    :param indices: Nx3 array of voxel indices
    :return: List of arrays, each containing a Mx3 numpy array of the indices of a single cluster of connected voxels
    """
    ...
