from typing import List

from skan import csr
from skan.csr import (skeleton_to_csgraph,
                      csr_to_nbgraph,
                        NBGraphBool,
                        sparse,
                        _compute_distances)

import numpy as np
import numpy.typing as npt


def _walk_path_lab(
        jgraph, node, neighbor, visited, degrees, indices, path_data, startj
        ):
    """Edited version of the_build_paths function from skan.
    A neighbor is redefined as an adjecent voxel with the 
    same label. This ensures that in infer-subc, branches can only come from
    one object, even if the organelle objects are in contact.

    To view original function (as of 2-20-26),
    ctrl + click https://skeleton-analysis.org/stable/_modules/skan/csr.html
    """
    # walk path troublshoot
    indices[startj] = node
    start_node = node
    path_data[startj] = jgraph.node_properties[node]
    j = startj + 1
    while not visited.edge(node, neighbor):
        visited.set_edge(node, neighbor, True)
        visited.set_edge(neighbor, node, True)
        indices[j] = neighbor
        path_data[j] = jgraph.node_properties[neighbor]
        if degrees[neighbor] != 2 or neighbor == start_node:
            break
        # Iterate to find the next valid neighbor in the same-value path
        nextneighbor = -1
        for n in jgraph.neighbors(neighbor):
            if n != node and jgraph.node_properties[n] == jgraph.node_properties[neighbor]:
                nextneighbor = n
                break
        node, neighbor = neighbor, nextneighbor
        j += 1
    return j - startj + 1

def _build_paths_lab(jgraph, indptr, indices, path_data, visited, degrees):
    """Edited version of the_build_paths function from skan.
    A neighbor is redefined as an adjecent voxel with the 
    same label. This ensures that in infer-subc, branches can only come from
    one object, even if the organelle objects are in contact.
    
    To view original function (as of 2-20-26),
    ctrl + click https://skeleton-analysis.org/stable/_modules/skan/csr.html
    """
    indptr_i = 0
    indices_j = 0
    # first, process all nodes in a path to an endpoint or junction
    for node in range(jgraph.shape[0]):
        if degrees[node] > 2 or degrees[node] == 1:
            for neighbor in jgraph.neighbors(node):
                # Only initiate path if the neighbor shares the same value
                if jgraph.node_properties[node] == jgraph.node_properties[neighbor]:
                    if not visited.edge(node, neighbor):
                        n_steps = _walk_path_lab(
                                jgraph,
                                node,
                                neighbor,
                                visited,
                                degrees,
                                indices,
                                path_data,
                                indices_j,
                                )
                        indptr[indptr_i + 1] = indptr[indptr_i] + n_steps
                        indptr_i += 1
                        indices_j += n_steps
    # everything else is by definition in isolated cycles
    # Find the first unvisited neighbor that shares the same value
    neighbor = -1
    for n in jgraph.neighbors(node):
        if jgraph.node_properties[n] == jgraph.node_properties[node] and not visited.edge(node, n):
            neighbor = n
            break
    if neighbor != -1:
        if degrees[node] > 0:
            neighbor = jgraph.neighbors(node)[0]
            if not visited.edge(node, neighbor):
                n_steps = _walk_path_lab(
                        jgraph,
                        node,
                        neighbor,
                        visited,
                        degrees,
                        indices,
                        path_data,
                        indices_j,
                        )
                indptr[indptr_i + 1] = indptr[indptr_i] + n_steps
                indptr_i += 1
                indices_j += n_steps
    return indptr_i + 1, indices_j

def _build_skeleton_path_graph_lab(graph):
    """Edited version of the _build_skeleton_path_graph function from skan.
    Connectivity is redefined as the amount of neighboring voxels with the 
    same label. This ensures that in infer-subc, branches can only come from
    one object, even if the organelle objects are in contact.
    
    To view original function (as of 2-20-26),
    ctrl + click https://skeleton-analysis.org/stable/_modules/skan/csr.html
    """
    
    max_num_cycles = graph.indices.size // 4
    buffer_size_offset = max_num_cycles
    # Calculate degree based only on label-matching neighbors
    degrees = np.zeros(graph.shape[0], dtype=np.int32)
    for i in range(graph.shape[0]):
        for j in range(graph.indptr[i], graph.indptr[i+1]):
            neighbor = graph.indices[j]
            if graph.node_properties[i] == graph.node_properties[neighbor]:
                degrees[i] += 1
    visited_data = np.zeros(graph.data.shape, dtype=bool)
    visited = NBGraphBool(
            graph.indptr.astype(np.int32, copy=False),
            graph.indices.astype(np.int32, copy=False), visited_data,
            graph.shape, np.broadcast_to(1.0, graph.shape[0])
            )
    endpoints = (degrees != 2)
    endpoint_degrees = degrees[endpoints]
    num_paths = np.sum(endpoint_degrees)
    path_indptr = np.zeros(num_paths + buffer_size_offset, dtype=int)
    # the number of points that we need to save to store all skeleton
    # paths is equal to the number of pixels plus the sum of endpoint
    # degrees minus one (since the endpoints will have been counted once
    # already in the number of pixels) *plus* the number of isolated
    # cycles (since each cycle has one index repeated). We don't know
    # the number of cycles ahead of time, but it is bounded by one quarter
    # of the number of points.
    n_points = (
            graph.indices.size + np.sum(np.maximum(0, endpoint_degrees - 1))
            + buffer_size_offset
            )
    path_indices = np.zeros(n_points, dtype=int)
    path_data = np.zeros(path_indices.shape, dtype=float)
    m, n = _build_paths_lab(
            graph, path_indptr, path_indices, path_data, visited, degrees
            )
    paths = sparse.csr_matrix(
            (path_data[:n], path_indices[:n], path_indptr[:m]),
            shape=(m - 1, n)
            )
    return paths

class Skeleton:
    """
    Modified Skeleton Class from skan to fit new assumptions
    for infer-subc quantification. Branches will only belong
    to one object, i.e., branches spanning multiple objects
    will be broken up by those labels.
    
    To view original class (as of 2-20-26),
    ctrl + click https://skeleton-analysis.org/stable/_modules/skan/csr.html#Skeleton
    
    Object to group together all the properties of a skeleton.

    In the text below, we use the following notation:

    - N: the number of points in the pixel skeleton,
    - ndim: the dimensionality of the skeleton
    - P: the number of paths in the skeleton (also the number of links in the
      junction graph).
    - J: the number of junction nodes
    - Sd: the sum of the degrees of all the junction nodes
    - [Nt], [Np], Nr, Nc: the dimensions of the source image

    Parameters
    ----------
    skeleton_image : array
        The input skeleton (1-pixel/voxel thick skeleton, all other values 0).

    Other Parameters
    ----------------
    spacing : float or array of float, shape ``(ndim,)``
        The scale of the pixel spacing along each axis.
    source_image : array of float, same shape as `skeleton_image`
        The image that `skeleton_image` represents / summarizes / was generated
        from. This is used to produce visualizations as well as statistical
        properties of paths.
    keep_images : bool
        Whether or not to keep the original input images. These can be useful
        for visualization, but they may take up a lot of memory.
    value_is_height : bool
        Whether to consider the value of a float skeleton to be the "height"
        of the image. This can be useful e.g. when measuring lengths along
        ridges in AFM images.

    Attributes
    ----------
    graph : scipy.sparse.csr_matrix, shape (N + 1, N + 1)
        The skeleton pixel graph, where each node is a non-zero pixel in the
        input image, and each edge connects adjacent pixels. The graph is
        represented as an adjacency matrix in SciPy sparse matrix format. For
        more information see the ``scipy.sparse`` documentation as well as
        ``scipy.sparse.csgraph``. Note: pixel numbering starts at 1, so the
        shape of this matrix is ``(N + 1, N + 1)`` instead of ``(N, N)``.
    nbgraph : NBGraph
        A thin Numba wrapper around the ``csr_matrix`` format, this provides
        faster graph methods. For example, it is much faster to get a list of
        neighbors, or test for the presence of a specific edge.
    coordinates : array, shape (N, ndim)
        skeleton_pixel_id i -> coordinates[i]
        The image coordinates of each pixel in the skeleton.
        Some values in this matrix are non-sensical — you should only access
        them from node ids.
    paths : scipy.sparse.csr_matrix, shape (P, N + 1)
        A csr_matrix where element [i, j] is on if node j is in path i. This
        includes path endpoints. The number of nonzero elements is N - J + Sd.
    n_paths : int
        The number of paths, P. This is redundant information given `n_paths`,
        but it is used often enough that it is worth keeping around.
    distances : array of float, shape (P,)
        The distance of each path. Note: not initialized until `path_lengths()`
        is called on the skeleton; use path_lengths() instead
    skeleton_image : array or None
        The input skeleton image. Only present if `keep_images` is True. Set to
        False to preserve memory.
    source_image : array or None
        The image from which the skeleton was derived. Only present if
        `keep_images` is True. This is useful for visualization.
    """
    def __init__(
            self,
            skeleton_image,
            *,
            spacing=1,
            source_image=None,
            keep_images=True,
            value_is_height=False,
            ):
        graph, coords = skeleton_to_csgraph(
                skeleton_image,
                spacing=spacing,
                value_is_height=value_is_height,
                )
        if np.issubdtype(skeleton_image.dtype, np.floating):
            self.pixel_values = skeleton_image[coords]
        elif np.issubdtype(skeleton_image.dtype, np.integer):
            self.pixel_values = skeleton_image.astype(np.float64)[coords]
        else:
            self.pixel_values = None
        self.graph = graph
        self.nbgraph = csr_to_nbgraph(graph, self.pixel_values)
        self.coordinates = np.transpose(coords)
        self.paths = _build_skeleton_path_graph_lab(self.nbgraph)
        self.n_paths = self.paths.shape[0]
        self.distances = np.empty(self.n_paths, dtype=float)
        self._distances_initialized = False
        self.skeleton_image = None
        self.skeleton_shape = skeleton_image.shape
        self.skeleton_dtype = skeleton_image.dtype
        self.source_image = None
        self.degrees = np.diff(self.graph.indptr)
        # added for infer_subc quantification
        # number of neighboring pixels that contain the same label
        self.correct_degrees = np.zeros(self.nbgraph.shape[0], dtype=np.int32)
        for i in range(self.nbgraph.shape[0]):
            for j in range(self.nbgraph.indptr[i], self.nbgraph.indptr[i+1]):
                neighbor = self.nbgraph.indices[j]
                if self.nbgraph.node_properties[i] == self.nbgraph.node_properties[neighbor]:
                    self.correct_degrees[i] += 1
        self.spacing = (
                np.asarray(spacing) if not np.isscalar(spacing) else
                np.full(skeleton_image.ndim, spacing)
                )
        if keep_images:
            self.keep_images = keep_images
            self.skeleton_image = skeleton_image
            self.source_image = source_image
            
    def path(self, index):
            """Return the pixel indices of path number `index`.

            Parameters
            ----------
            index : int
                The desired path.

            Returns
            -------
            path : array of int
                The indices of the pixels belonging to the path, including
                endpoints.
            """
            # The below is equivalent to `self.paths[index].indices`, which is much
            # more elegant. However the below version is about 25x faster!
            # In [14]: %timeit mat[1].indices
            # 128 µs ± 421 ns per loop (mean ± std. dev. of 7 runs, 10000 loops each)
            # In [16]: %%timeit
            # ...: start, stop = mat.indptr[1:3]
            # ...: mat.indices[start:stop]
            # ...:
            # 5.05 µs ± 77.2 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
            start, stop = self.paths.indptr[index:index + 2]
            return self.paths.indices[start:stop]

    def path_coordinates(self, index: int):
        """Return the image coordinates of the pixels in the path.

        Parameters
        ----------
        index : int
            The desired path.

        Returns
        -------
        path_coords : array of float
            The (image) coordinates of points on the path, including endpoints.
        """
        path_indices = self.path(index)
        return self.coordinates[path_indices]

    def path_with_data(self, index: int):
        """Return pixel indices and corresponding pixel values on a path.

        Parameters
        ----------
        index : int
            The desired path.

        Returns
        -------
        path : array of int
            The indices of pixels on the path, including endpoints.
        data : array of float
            The values of pixels on the path.
        """
        start, stop = self.paths.indptr[index:index + 2]
        return self.paths.indices[start:stop], self.paths.data[start:stop]

    def path_lengths(self):
        """Return the length of each path on the skeleton.

        Returns
        -------
        lengths : array of float
            The length of all the paths in the skeleton.
        """
        if not self._distances_initialized:
            _compute_distances(
                    self.nbgraph, self.paths.indptr, self.paths.indices,
                    self.distances
                    )
            self._distances_initialized = True
        return self.distances

    def paths_list(self):
        """List all the paths in the skeleton, including endpoints.

        Returns
        -------
        paths : list of array of int
            The list containing all the paths in the skeleton.
        """
        return [list(self.path(i)) for i in range(self.n_paths)]

    def path_label_image(self):
        """Image like self.skeleton_image with path_ids as values.

        Returns
        -------
        label_image : array of ints
            Image of the same shape as self.skeleton_image where each pixel
            has the value of its branch id + 1.
        """
        image_out = np.zeros(self.skeleton_shape, dtype=int)
        for i in range(self.n_paths):
            coords_to_wipe = self.path_coordinates(i)
            coords_idxs = tuple(np.round(coords_to_wipe).astype(int).T)
            image_out[coords_idxs] = i + 1
        return image_out

    def path_means(self):
        """Compute the mean pixel value along each path.

        Returns
        -------
        means : array of float
            The average pixel value along each path in the skeleton.
        """
        sums = np.add.reduceat(self.paths.data, self.paths.indptr[:-1])
        lengths = np.diff(self.paths.indptr)
        return sums / lengths

    def path_stdev(self):
        """Compute the standard deviation of values along each path.

        Returns
        -------
        stdevs : array of float
            The standard deviation of pixel values along each path.
        """
        data = self.paths.data
        sumsq = np.add.reduceat(data * data, self.paths.indptr[:-1])
        lengths = np.diff(self.paths.indptr)
        means = self.path_means()
        return np.sqrt(np.clip(sumsq/lengths - means*means, 0, None))

    def prune_paths(self, indices: npt.ArrayLike) -> 'eSkeleton':
        """Prune nodes from the skeleton.

        Parameters
        ----------
        indices: List[int]
            List of indices to be removed.

        Retruns
        -------
        Skeleton
            A new Skeleton object pruned.
        """
        # warning: slow
        image_cp = np.copy(self.skeleton_image)
        if not np.all(np.array(indices) < self.n_paths):
            raise ValueError(
                    f'The path index {np.max(indices)} does not exist in this '
                    f'skeleton. (The highest path index is {self.n_paths}.)\n'
                    'If you obtained the index from a summary table, you '
                    'probably need to resummarize the skeleton.'
                    )
        for i in indices:
            pixel_ids_to_wipe = self.path(i)
            junctions = self.degrees[pixel_ids_to_wipe] > 2
            pixel_ids_to_wipe = pixel_ids_to_wipe[~junctions]
            coords_to_wipe = self.coordinates[pixel_ids_to_wipe]
            coords_idxs = tuple(np.round(coords_to_wipe).astype(int).T)
            image_cp[coords_idxs] = 0
        # optional cleanup:
        new_skeleton = morphology.skeletonize(image_cp.astype(bool)) * image_cp
        return eSkeleton(
                new_skeleton,
                spacing=self.spacing,
                source_image=self.source_image,
                keep_images=self.keep_images
                )

    def __array__(self, dtype=None):
        """Array representation of the skeleton path labels."""
        return self.path_label_image()

def skel_euclid_dist(skel: Skeleton) -> np.ndarray:
    """
    A function that returns the euclidean distances for each branch in the skeleton object
    This code is borrowed from skan's csr.summarize() function which can be found here https://skeleton-analysis.org/stable/_modules/skan/csr.html#summarize

    Parameters
    ------------
    skel:
        the Skeleton() object that contains all information about the branches

    Returns
    -------------
    An array of euclidean distances for every branch in the skeleton 
    """
    # this code was borrowed from skan summarize function
    ndim = skel.coordinates.shape[1]
    endpoints_src = skel.paths.indices[skel.paths.indptr[:-1]]
    endpoints_dst = skel.paths.indices[skel.paths.indptr[1:] - 1]
    coords_real_src = skel.coordinates[endpoints_src] * skel.spacing
    coords_real_dst = skel.coordinates[endpoints_dst] * skel.spacing
    euclid_dists = np.sqrt((coords_real_dst - coords_real_src)**2
                    @ np.ones(ndim))
            
    return(euclid_dists)

def branch_type(skel: Skeleton, branch_list: List):
    """
    A function that returns the branch type given a list of branch ids
    This code is borrowed from skan's csr.summarize() function which can be found here https://skeleton-analysis.org/stable/_modules/skan/csr.html#summarize

    Parameters
    ------------
    skel:
        the Skeleton() object that contains all information about the branches
    branch_list:
        a list of branch ids
    Returns
    -------------
    An array of branch types for every branch in the branch_list np.ndarray
    """
    # this is a piece of the summarize function from skan
    # This should be identical logic from the branch table function
    endpoints_src = skel.paths.indices[skel.paths.indptr[:-1]]
    endpoints_dst = skel.paths.indices[skel.paths.indptr[1:] - 1]
    deg_src = skel.correct_degrees[endpoints_src]
    deg_dst = skel.correct_degrees[endpoints_dst]
    kind = np.full(deg_src.shape, 2)  # default: junction-to-junction
    kind[(deg_src == 1) | (deg_dst == 1)] = 1  # tip-junction
    kind[(deg_src == 1) & (deg_dst == 1)] = 0  # tip-tip
    kind[endpoints_src == endpoints_dst] = 3  # cycle
    return kind[branch_list]