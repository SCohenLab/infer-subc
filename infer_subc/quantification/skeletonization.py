from typing import Union, Tuple, List

from skan import csr
from skan.csr import (skeleton_to_csgraph,
                      csr_to_nbgraph,
                        NBGraphBool,
                        sparse,
                        _compute_distances)

from skimage.morphology import skeletonize
from skimage.measure import label
from scipy.sparse import csgraph
from scipy.spatial import cKDTree
from scipy import ndimage
from collections import defaultdict

import numpy as np
import numpy.typing as npt
import pandas as pd

from infer_subc.core.img import apply_mask
from infer_subc.quantification.stats import *

def create_skel(segmentation: np.ndarray) -> np.ndarray:
    """
    A function that skeletonizes the organelle segmentation. This function also generates punctate objects 
    for the round organelle objects and declumped organelles only if they lack a skeleton. This function is
    based off of skimage's skeletonize function. More information about said function can be found here 
    https://scikit-image.org/docs/0.25.x/api/skimage.morphology.html#skimage.morphology.skeletonize

    Parameters
    ------------
    segmentation : array
        the segmentated organelle image as a numpy array. It is assumed that the segmentation has already
        been masked and only contains the organelle of interest.

    Returns
    -------------
    skel_arr : array
        A properly skeletonized np.ndarray with float labels due to skan requirements. Labels correspond to the original
        infer-subc segmentation labels.
    """

    # where the organelles exist
    omask = segmentation > 0

    # This is the raw organelle skeleton, some fixing and relabeling has to be done before we can use the skeleton for computation
    skeleton = skeletonize(omask).astype(bool)

    # relabel segmentation (makes sure that all disconnected components have a skeleton, even in the case of the ER)
    comp_seg = label(segmentation)

    # get maximum value in component segmentation labels
    max_val = np.max(segmentation)

    # calculate the shift needed to label via orginal infer-subc labeling
    if max_val > 1:
        shift = 10**(int(np.floor(np.log10(max_val))) + 1)
        # Use in-place math to save memory
        new_seg = np.zeros_like(segmentation, dtype=np.uint64)
        new_seg[omask] = comp_seg[omask].astype(np.uint64) * shift
        new_seg[omask] += segmentation[omask].astype(np.uint64)
    else:
        new_seg = comp_seg

    # All of disconnected object labels
    all_lab = np.unique(new_seg[omask])

    # Applying the segmentation labels to the raw skeleton
    lab_skel = skeleton * new_seg

    # Labels present in the skeleton
    skel_lab = np.unique(lab_skel[lab_skel > 0])

    # collect the missing labels
    mis_lab = np.setdiff1d(all_lab, skel_lab, assume_unique=True)

    # check if there are missing labels
    if len(mis_lab) > 0:

        # get regionprops for the new segmentation
        # this helps us quickly find the coordinates of the missing labels and assign them to the skeleton
        props = regionprops(new_seg)
    
        # For each missing label, find a coordinate to place a point in the skeleton (centroid)
        for prop in props:
            if prop.label in mis_lab:
                # get the centroid coordinate
                cent_coord = tuple(np.round(prop.centroid).astype(int))
                
                # check if centroid is within the label
                if new_seg[cent_coord] == prop.label:
                    lab_skel[cent_coord] = prop.label
                else:
                    # pick a coordinate listed within the label
                    middle_coord = tuple(prop.coords[len(prop.coords) // 2])
                    lab_skel[middle_coord] = prop.label

    # Reapply original labels and return
    return ((lab_skel > 0) * segmentation).astype(float)

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
    n_branches : int
        The number of paths (branches in infer-subc), P. This is redundant information given `n_branches`,
        but it is used often enough that it is worth keeping around.
    distances : array of float, shape (P,)
        The distance of each path. Note: not initialized until `branch_lengths()`
        is called on the skeleton; use branch_lengths() instead
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
        self.n_branches = self.paths.shape[0]
        self.distances = np.empty(self.n_branches, dtype=float)
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

    def branch_lengths(self):
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
        return [list(self.path(i)) for i in range(self.n_branches)]

    def path_label_image(self):
        """Image like self.skeleton_image with path_ids as values.

        Returns
        -------
        label_image : array of ints
            Image of the same shape as self.skeleton_image where each pixel
            has the value of its branch id + 1.
        """
        image_out = np.zeros(self.skeleton_shape, dtype=int)
        for i in range(self.n_branches):
            coords_to_wipe = self.path_coordinates(i)
            coords_idxs = tuple(np.round(coords_to_wipe).astype(int).T)
            image_out[coords_idxs] = i + 1
        return image_out

    def branch_means(self):
        """Compute the mean pixel value along each path.

        Returns
        -------
        means : array of float
            The average pixel value along each path in the skeleton.
        """
        sums = np.add.reduceat(self.paths.data, self.paths.indptr[:-1])
        lengths = np.diff(self.paths.indptr)
        return sums / lengths

    def branch_stdev(self):
        """Compute the standard deviation of values along each path.

        Returns
        -------
        stdevs : array of float
            The standard deviation of pixel values along each path.
        """
        data = self.paths.data
        sumsq = np.add.reduceat(data * data, self.paths.indptr[:-1])
        lengths = np.diff(self.paths.indptr)
        means = self.branch_means()
        return np.sqrt(np.clip(sumsq/lengths - means*means, 0, None))

    def prune_paths(self, indices: npt.ArrayLike) -> 'Skeleton':
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
        if not np.all(np.array(indices) < self.n_branches):
            raise ValueError(
                    f'The path index {np.max(indices)} does not exist in this '
                    f'skeleton. (The highest path index is {self.n_branches}.)\n'
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
        return Skeleton(
                new_skeleton,
                spacing=self.spacing,
                source_image=self.source_image,
                keep_images=self.keep_images
                )

    def __array__(self, dtype=None):
        """Array representation of the skeleton path labels."""
        return self.path_label_image()

def get_obj_ids(skel: Skeleton) -> np.ndarray:
    """
    A function that returns a np.ndarray (int) of object IDs for each branch in the skeleton object.

    Parameters
    ------------
    skel:
        the Skeleton() graph that contains all information about the skeleton

    Returns
    -------------
    Array of object IDs per branch in the skeleton graph. This is used to assign the original segmentation labels to the skeleton branches.
    """
    # checker to see if all path points and nodes come from the same object (per branch)
    if not np.any(skel.branch_stdev()):
        return skel.branch_means().astype(int)
    else:
        raise ValueError("at least one branch spans across multiple different organelle objects")

def get_skel_branch(skel: Skeleton) -> pd.DataFrame:
    """
    A function that produces branch level quantification for Skeleton() graph objects in table format. 
    This code is based off of skan's csr.summarize() function which can be found here https://skeleton-analysis.org/stable/_modules/skan/csr.html#summarize

    Parameters
    ------------
    skel:
        the Skeleton() graph from which quantification will be based on


    Branch table measurements:
    ------------------------
    'skel_obj_id',
    'point_id_src',
    'point_id_dst',
    'deg_src',
    'deg_dst',
    'branch_length',
    'branch_type',
    'image_coord_src_0',
    'image_coord_src_1',
    'image_coord_src_2',
    'image_coord_dst_0',
    'image_coord_dst_1',
    'image_coord_dst_2',
    'coord_src_0',
    'coord_src_1',
    'coord_src_2',
    'coord_dst_0',
    'coord_dst_1',
    'coord_dst_2',
    'euclidean_distance',
    'str_prop'

    Returns
    -------------
    pandas dataframe of containing measurements (columns) for each branch (rows) in the skeleton 
    """
    summary = {}
    summary['skel_obj_id'] = get_obj_ids(skel)
    ndim = skel.coordinates.shape[1]
    
    endpoints_src = skel.paths.indices[skel.paths.indptr[:-1]]
    endpoints_dst = skel.paths.indices[skel.paths.indptr[1:] - 1]

    summary['point_id_src'] = endpoints_src
    summary['point_id_dst'] = endpoints_dst
    deg_src = skel.correct_degrees[endpoints_src]
    deg_dst = skel.correct_degrees[endpoints_dst]
    summary['deg_src'] = deg_src
    summary['deg_dst'] = deg_dst
    
    kind = np.full(deg_src.shape, 2)  # default: junction-to-junction
    kind[(deg_src == 1) | (deg_dst == 1)] = 1  # tip-junction
    kind[(deg_src == 1) & (deg_dst == 1)] = 0  # tip-tip
    kind[endpoints_src == endpoints_dst] = 3  # cycle
    summary['branch_type'] = kind
    for i in range(ndim):  # keep loops separate for best insertion order
        summary[f'image_coord_src_{i}'] = skel.coordinates[endpoints_src, i]
    for i in range(ndim):
        summary[f'image_coord_dst_{i}'] = skel.coordinates[endpoints_dst, i]
    coords_real_src = skel.coordinates[endpoints_src] * skel.spacing
    for i in range(ndim):
        summary[f'coord_src_{i}'] = coords_real_src[:, i]
    coords_real_dst = skel.coordinates[endpoints_dst] * skel.spacing
    for i in range(ndim):
        summary[f'coord_dst_{i}'] = coords_real_dst[:, i]
        
    summary['branch_length'] = skel.branch_lengths()
    summary['euclidean_distance'] = (
            np.sqrt((coords_real_dst - coords_real_src)**2
                    @ np.ones(ndim))
            )

    summary['str_prop'] = summary['euclidean_distance'] / summary['branch_length']
    return pd.DataFrame(summary).rename_axis('branch_id')

def get_skel_node(skel: Skeleton) -> pd.DataFrame:
    """
    A function that produces node level quantification for Skeleton() graph objects in table format. 

    Parameters
    ------------
    skel:
        the Skeleton() graph from which quantification will be based on


    Branch table measurements:
    ------------------------
    'point_id',
    'node_type',
    'connectivity',
    'image_coord_0',
    'image_coord_1',
    'image_coord_2',
    'coord_0',
    'coord_1',
    'coord_2',
    'branch_ids',
    'obj_id'

    Returns
    -------------
    pandas dataframe of containing measurements (columns) for each node (rows) in the skeleton 
    """
    # gets a list of all non path points (nodes)
    node_list = np.nonzero(skel.correct_degrees != 2)[0]

    # A dictionary that reveals the branches a node is referenced in
    node2branches = dict((point_id,[]) for point_id in node_list)

    # identify the branches each node is involved in
    for branch in range(skel.n_branches):
        for point in skel.path(branch):
            if skel.correct_degrees[point] != 2:
                node2branches[point] += [branch]

    node_lab = []
    # in the case of the ER absolute punctates are not their own skeleton object
    base_n = ["Abs Punctate",
            "Endpoint",
            "Path Point"]

    node_lab = base_n
    if np.max(skel.correct_degrees) > 2:
        for i in np.arange(3, np.max(skel.correct_degrees) + 1):
            node_lab += [f"{i}-way"]
        
    node_table_data = {
        "point_id": node_list,
        "node_type": np.array(node_lab)[skel.correct_degrees[node_list]],
        "connectivity": skel.correct_degrees[node_list],
        "image_coord_0": skel.coordinates[node_list,0],
        "image_coord_1": skel.coordinates[node_list,1],
        "image_coord_2": skel.coordinates[node_list,2],
        "coord_0": skel.coordinates[node_list,0] * skel.spacing[0],
        "coord_1": skel.coordinates[node_list,1] * skel.spacing[1],
        "coord_2": skel.coordinates[node_list,2] * skel.spacing[2],
        "branch_ids": [node2branches[point] for point in node_list],
        'obj_id': [int(skel.pixel_values[point]) for point in node_list]
    }

    return pd.DataFrame(node_table_data)

def skel_euclid_dist(skel: Skeleton) -> np.ndarray:
    """
    A function that returns the euclidean distances for each branch in the skeleton object
    This code is borrowed from skan's csr.summarize() function which can be found here https://skeleton-analysis.org/stable/_modules/skan/csr.html#summarize
    This is used in the get_skel_obj function.

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
    This is used in the get_skel_obj function.
    
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

def skel_width(skel: Skeleton, segmentation: np.ndarray, obj_list: np.ndarray) -> list:
    """
    A function that takes in a list of skeleton object ids and calculates the width of each skeleton object.
    The width is calculated by finding the distance of the nearest boundary point (for each point in the skeleton object)
    and taking the average of those distances (multiplied by 2 to get the full width). This function uses similar logic to 
    the mitograph method for calculating width, which can be found here https://github.com/vianamp/MitoGraph/blob/master/MitoGraph.cxx#L1154
    While this measurement is informative, the function is limited by the resolution of the image. Thus, width values of zero can result 
    from segmented objects that are a few voxels wide.

    Parameters
    ----------
    skel: Skeleton
        the Skeleton() graph from which quantification will be based on
    segmentation: np.ndarray
        the binary segmentation image
    obj_list: np.ndarray
        a list of object ids for which to calculate skeleton widths
        Ideally this list should be ordered in the same fashion as the rows in the skeleton object table

    Returns
    -------
    skel_width
        a dictionary mapping each object id to its calculated skeleton width
    """
    ##############################################
    # CREATE HOLLOWED OUT (SURFACE) SEGMENTATION
    ##############################################

    # erode the organelle segmentation
    eroded = ndimage.binary_erosion(segmentation)
    # perform logical xor to get the hollowed out (surface) organelle segmentation and store in eroded to save memory
    np.logical_xor(segmentation, eroded, out=eroded)
    # apply original labels to the hollowed out segmentation
    hollow_labels = (eroded * segmentation)
    # scale the skeleton coordinates
    coords_scaled = skel.coordinates * skel.spacing
    
    ##########################################################
    # COLLECT INDICES AND LABELS OF HOLLOWED OUT SEGMENTATION
    ##########################################################
    
    # collect the coordinates of the hollowed out segmentation points
    hp_indices = np.argwhere(hollow_labels > 0)
    # collect the labels of the hollowed out segmentation points (.T transposes the array in correct formatting for indexing)
    hp_labels = hollow_labels[tuple((hp_indices).T)]

    # dictionary that contains the coordinates of the hollowed out segmentation points for each label (object)
    hp_by_label = defaultdict(list)
    for idx, label in zip(hp_indices, hp_labels):
        hp_by_label[label].append(idx)

    # dictionary that contains the coordinates of the skeleton points for each label (object)
    sp_by_label = defaultdict(list)
    for idx, label in enumerate(skel.pixel_values):
        sp_by_label[label].append(coords_scaled[idx])

    ##########################################################
    # CALCULATE SKELETAL WIDTHS FOR EACH SKELETON OBJECT
    ##########################################################

    # intialize skeletal width list (will be ordered in the same fashion as the obj_list)
    skel_width = []

    for obj in obj_list:
        # collect the coordinates of the hollowed out (surface)segmentation points that belong to the object (scaled)
        hp = np.array(hp_by_label.get(obj, [])) * skel.spacing
        # collect the coordinates of the skeleton points that belong to the object (scaled)
        sp = np.array(sp_by_label.get(obj, []))
        
        if hp.size > 0 and sp.size > 0:
            # create a KD-tree from the hollowed out segmentation points (KD-tree is a data structure that allows for rapid nearest neighbor searches)
            tree = cKDTree(hp)
            # calculate the distances from each skeleton point (this functions allows all of the skeleton points to be queried at once)
            distances, _ = tree.query(sp, k=1, p=2)
            # average the distances and multiply by 2 to get the average width of the skeleton
            skel_width.append(2 * np.mean(distances))
        else:
            skel_width.append(0)
    
    return skel_width


def get_skel_obj(skel: Skeleton, segmentation: np.ndarray) -> pd.DataFrame:
    """
    A function that produces skeleton object level quantification for Skeleton() graph objects in table format. 

    Parameters
    ------------
    skel:
        the Skeleton() graph from which quantification will be based on
    segmentation:
        the original segmentation image with int labels that correspond to the skeleton object labels


    Branch table measurements:
    ------------------------
    'obj_id',
    'skel_type',
    'skel_type_num',
    'brh_count',
    'branch_ids',
    'min_brh_length',
    'max_brh_length',
    'ave_brh_length',
    'sd_brh_length',
    'med_brh_length',
    'total_length',
    'brh_type_0_tot',
    'brh_type_0_id',
    'brh_type_1_tot',
    'brh_type_1_ids',
    'brh_type_2_tot',
    'brh_type_2_ids',
    'brh_type_3_tot',
    'brh_type_3_ids',
    'comp_count',
    'node_count',
    'abs_punc_count',
    'ep_count',
    'jn_count',
    'ave_jn_deg',
    'max_deg',
    'point_ids',
    'mean_brh_str'
    'med_brh_str'
    'sd_brh_str'
    'width'


    Returns
    -------------
    pandas dataframe of containing measurements (columns) for each skeleton object (rows) in the skeleton 
    """

    #### DICTIONARIES AND ARRAYS ####
    
    # ordered list of all skeleton object ids (including punctates)
    # numpy sorts the unique values automatically
    obj_list = np.unique(skel.pixel_values).astype(int)

    brh_tlist = dict()
    for i in range(4):
        brh_tlist[i] = dict((obj,[]) for obj in obj_list)

    # lengths of the branches in the skeleton object
    obj_bl = dict((obj,[]) for obj in obj_list)

    # branch ids for the branches within a skeleton object
    obj2branch = dict((obj,[]) for obj in obj_list)

    # straightness ratios between distance and branch length for each branch in the skeleton object
    obj_str = dict((obj,[]) for obj in obj_list)

    # returns array of euclidean distances
    euclid_dist_arr = skel_euclid_dist(skel)

    # returns array of branch types
    brh_type = branch_type(skel, range(skel.n_branches))

    #### BRANCH MEASUREMENTS ####

    for br in range(skel.n_branches):
        obj = skel.path_with_data(br)[1][0].astype(int)

        #### OBJECT BRANCHES ####
        obj2branch[obj] += [br]
        
        #### BRANCH TYPE LIST ####
        brh_tlist[brh_type[br]][obj] += [br]

        #### BRANCH LENGTH CALCULATIONS ####
        obj_bl[obj] += [skel.branch_lengths()[br]]
        obj_str[obj] += [euclid_dist_arr[br] / skel.branch_lengths()[br]]
    
    #### SKELETON TYPE ####

    # all referenced dictionaries are in the same order as the obj_list

    # creating a sum branch length array to avoid recalling the dictionary over and over
    sum_bl = [np.sum(bl) for bl in obj_bl.values()]

    # creating a branch count for the same reason
    count_b = [len(bl) for bl in obj_bl.values()]

    # creating a cycle count for the same reason
    count_cyc = [len(ids) for ids in brh_tlist[3].values()]

    # the highest value of length where a type 0 branch will still be considered punctate
    scale = skel.spacing
    p_threshold = min(scale) * 2

    obj_type_n = np.full(obj_list.shape, 1) # Rod as default
    obj_type_n[np.array(sum_bl) <= p_threshold] = 0 # If total length below punctate threshold, then punctate
    obj_type_n[np.array(count_cyc) == 1] = 2 # If object has one cycle, then isolated cycle
    obj_type_n[np.array(count_b) > 1] = 3 # If object has multiple branches, then network

    ### this line is not needed as none of the other conditions would be met in the case of an absolute punctate
    # obj_class_n[np.array(count_b) < 1] = 0 # If object has no branches, then punctate (absolute punctates)

    obj_type = np.array(('Punctate','Rod','Isolated Cycle','Network'))[obj_type_n]
            
    #### CONNECTED COMPONENTS ####
    # using connected component id, ordered by obj_list

    comp_arr = []
    for obj in obj_list:
        comp_ids = csgraph.connected_components(skel.graph, directed=False)[1][skel.pixel_values == obj]
        comp_arr += [len(pd.unique(comp_ids))]

    #### SKELETAL WIDTHS ####
    widths = skel_width(skel, segmentation, obj_list)
    
    #### GENERATE SKELETON OBJECT TABLE ####

    skel_table_data = {
            "obj_id": obj_list,
            "skel_type": obj_type,
            "skel_type_num": obj_type_n,
            "brh_count": count_b,
            "branch_ids": [obj2branch[obj] for obj in obj_list],
            "min_brh_length": [np.min(obj_bl[obj]) if len(obj_bl[obj]) != 0 else np.nan for obj in obj_list],
            "max_brh_length": [np.max(obj_bl[obj]) if len(obj_bl[obj]) != 0 else np.nan for obj in obj_list],
            "ave_brh_length": [np.mean(obj_bl[obj])for obj in obj_list],
            "sd_brh_length": [np.std(obj_bl[obj]) for obj in obj_list],
            "med_brh_length": [np.median(obj_bl[obj])for obj in obj_list],
            "total_length": sum_bl,
            "brh_type_0_tot": [len(brh_tlist[0][obj]) for obj in obj_list],
            "brh_type_0_id": [brh_tlist[0][obj] for obj in obj_list],
            "brh_type_1_tot": [len(brh_tlist[1][obj]) for obj in obj_list],
            "brh_type_1_ids": [brh_tlist[1][obj] for obj in obj_list],
            "brh_type_2_tot": [len(brh_tlist[2][obj]) for obj in obj_list],
            "brh_type_2_ids": [brh_tlist[2][obj] for obj in obj_list],
            "brh_type_3_tot": count_cyc,
            "brh_type_3_ids": [brh_tlist[3][obj] for obj in obj_list],
            "comp_count": comp_arr,
            "node_count": [np.count_nonzero((skel.correct_degrees[skel.pixel_values == i]) != 2) for i in obj_list],
            'abs_punc_count': [np.count_nonzero((skel.correct_degrees[skel.pixel_values == i]) == 0) for i in obj_list],
            'ep_count': [np.count_nonzero((skel.correct_degrees[skel.pixel_values == i]) == 1) for i in obj_list],
            'jn_count' : [np.count_nonzero((skel.correct_degrees[skel.pixel_values == i]) > 2) for i in obj_list],
            'ave_jn_deg' : [np.mean(skel.correct_degrees[(skel.pixel_values == i) & (skel.correct_degrees > 2)]) for i in obj_list],
            'max_deg' : [np.max(skel.correct_degrees[skel.pixel_values == i]) for i in obj_list],
            "point_ids": [np.arange(skel.graph.shape[0])[skel.pixel_values == i] for i in obj_list],
            "mean_brh_str": [np.mean(obj_str[obj]) for obj in obj_list],
            "med_brh_str": [np.median(obj_str[obj]) for obj in obj_list],
            "sd_brh_str": [np.std(obj_str[obj]) for obj in obj_list],
            "width": widths}
    
    return pd.DataFrame(skel_table_data)

def get_skeleton_metrics(org_skel_arr: np.ndarray,
                          seg_name: str,
                          segmentation: np.ndarray,
                          mask_name: str,
                          mask: np.ndarray,
                     scale: Union[Tuple, None] = None,
                     output_all_tables: bool = False):
    
    '''
    A wrapper function that returns quantification describing the skeletonized np.ndarray segmentation.
    The default output is one pandas table where each row corresponds to one organelle object (skeleton object table).
    The branch and node table can also be generated as optional output

    Parameters
    ------------
    org_skel_arr:
        the labeled organelle skeleton np.ndarray (output of create_skel function)
    seg_name:
        the name of the organelle segmentation (string)
    segmentation:
        the original segmentation image with int labels that correspond to the skeleton object labels
    mask_name : str, optional
        The name of the cell mask, by default "cell".
        If None, "whole_image" will be used.
    mask : np.ndarray
        A binary image array representing the cell (or other) mask. 
        If None is provided, a whole image mask will be used.
    scale:
        the real world dimensions of the image (Tuple, ndarray or list of 3 floats)
    output_all_tables:
        False - return skeleton object table
        True - return branch, node and skeleton object table

    Returns
    -------------

    if output_all_table == True:
        pandas dataframe of containing measurements (columns) for each skeleton object (rows) in the skeleton
    if output_all_table == False:
        3 pandas dataframes of containing measurements for each branch, node and skeleton object respectively

    '''

    ##############################################################################
    # GENERATE SKELETON GRAPH
    ##############################################################################
    
    # mask skeleton
    if mask_name is None and mask is not None:
        raise ValueError("The mask_name parameter must be provided if mask is not None")
    elif mask is None and mask_name is None:
        mask_name = "whole_image"
    else:
        org_skel_arr = apply_mask(org_skel_arr, mask)

    # create skeleton
    org_skel = Skeleton(
        skeleton_image = org_skel_arr,
        spacing = scale,
        value_is_height = False)


    ##############################################################################
    # GENERATE AND RETURN TABLES
    ##############################################################################
    
    skel_table = get_skel_obj(org_skel, segmentation)
    skel_table.insert(0, "object", seg_name)
    skel_table.insert(0, "scale", str(scale))
    skel_table.insert(0, column="mask_name", value=mask_name)

    if output_all_tables:
        branch_table = get_skel_branch(org_skel)
        node_table = get_skel_node(org_skel)
        return branch_table, node_table, skel_table.rename(columns={"obj_id": "label"})
    else:
        return skel_table.rename(columns={"obj_id": "label"})
    
def fission_score(skel: Skeleton) -> float:

    """
    This measurement is derived from Spurlock, Xie, Song, Ricketts et al. paper "Mitochondrial fusion and cristae 
    reorganization facilitate acquisition of cardiomyocyte identity during reprogramming of murine fibroblasts"
    This paper can be found here https://pmc.ncbi.nlm.nih.gov/articles/PMC11973714/

    Parameters
    ------------
    skel:
        the Skeleton() graph containing all information about the skeleton

    Returns
    -------------
    fission score for the skeleton
    """

    fis_score = (skel.n_objects + skel.n_nodes + skel.n_branches)/(skel.total_length) if skel.total_length > 0 else 0
    return fis_score

