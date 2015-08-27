
#=======================================================================================================================
# Module implements edge detection and edge visualization.
#=======================================================================================================================

#cython: wraparound=False
#cython: boundscheck=False
#cython: cdivision=True


from libc.stdlib cimport abs, abort, malloc, free
import numpy as np
cimport numpy as cnp



cdef inline cnp.uint8_t _get_edge(cnp.uint8_t[:,::1] image, int z_row, int z_col, \
                                  int v_row, int v_col, int w_row, int w_col, int min_contrast) nogil:
    """
    Determines whether an oriented edge is present. 
    
    Parameters
    ----------
    image : array of shape (nrow, ncol)
        Intensity array.
    (z_row, z_col) : (int, int)
        Location of potential edge.
    (v_row, v_col) : (int, int)
        Direction of potential edge.
    (w_row, w_col) : (int, int)
        Orthogonal direction.
    min_contrast : int
        Minimum intensity difference. 
    
    Returns
    -------
    edge_type : {0, 1, 2} uint8
        0 --- no edge
        1 --- edge with orientation v (higher intensity at z then z+v)
        2 --- edge with orientation -v
    """
    
    # calculate intensity difference for potential edge 
    cdef int y_row = z_row + v_row
    cdef int y_col = z_col + v_col
    cdef cnp.uint8_t image_z = image[z_row, z_col]
    cdef cnp.uint8_t image_y = image[y_row, y_col]
    cdef int d = abs(image_z - image_y)
    
    # no edge if the difference is below the minimum contrast threshold
    if d < min_contrast: return 0
    
    # no edge if the difference is smaller than (or equal to) the difference between other neighbors
    if d <= abs(image_z - image[z_row-v_row, z_col-v_col]): return 0
    if d <= abs(image_z - image[z_row+w_row, z_col+w_col]): return 0
    if d <= abs(image_z - image[z_row-w_row, z_col-w_col]): return 0
    if d <= abs(image_y - image[y_row+v_row, y_col+v_col]): return 0
    if d <= abs(image_y - image[y_row+w_row, y_col+w_col]): return 0
    if d <= abs(image_y - image[y_row-w_row, y_col-w_col]): return 0
    
    # return polarity of edge
    if image_z > image_y: return 1
    else: return 2
    
    
def get_edges(images, int min_contrast=20, int spreading=0, binary=False):
    """
    Creates edge maps corresponding to the four edge types. 
    
    Parameters
    ----------
    images : array of shape ([nimg,] nrow, ncol)
        Intensity array(s).
    min_contrast : int, optional
        Minimum intensity difference.
    spreading : int, optional
        Spreads edges to a (spreading+1)^2 neighborhood. 
    binary = bool, optional
        True returns 8 binary (0=no edge, 1=edge) features,
        False returns 4 ternary (0=no edge, 1=edge of polarity1, 2=edge of polarity2) features.
    
    Returns
    -------
    edge_maps : array of shape ([nimg,] nrow, ncol, ntype)
        Detected edges.
    """
    
    # setup image collection
    cdef cnp.uint8_t[:,:,::1] my_images
    if images.ndim == 2:
        my_images = images.reshape(1, images.shape[0], images.shape[1]) 
    else:
        my_images = images
   
    # create storage for edge maps
    cdef int nimg = my_images.shape[0]
    cdef int nrow = my_images.shape[1]
    cdef int ncol = my_images.shape[2]
    cdef cnp.uint8_t[:,:,:,::1] edges = np.zeros((nimg, nrow, ncol, 4), dtype=np.uint8)    
    
    ### detect edges ###
    
    cdef int i, z_row, z_col, r, c
    cdef cnp.uint8_t[:,::1] my_image
    
    if spreading == 0:
        # - no spreading
        
        for i in range(nimg):
            my_image = my_images[i,]
            
            for z_row in range(1, nrow-2):
                for z_col in range(2, ncol-2):
                    edges[i, z_row, z_col, 0] = _get_edge(my_image, z_row, z_col, 0, 1, 1, 0, min_contrast)
                    edges[i, z_row, z_col, 1] = _get_edge(my_image, z_row, z_col, 1, 1, 1, -1, min_contrast)
                    edges[i, z_row, z_col, 2] = _get_edge(my_image, z_row, z_col, 1, 0, 0, -1, min_contrast)
                    edges[i, z_row, z_col, 3] = _get_edge(my_image, z_row, z_col, 1, -1, -1, -1, min_contrast)
                    
    else:
        # - spreading
        
        for i in range(nimg):
            my_image = my_images[i,]
            
            for z_row in range(1, nrow-2):
                for z_col in range(2, ncol-2):
                    edges[i, z_row, z_col, 0] = _get_edge(my_image, z_row, z_col, 0, 1, 1, 0, min_contrast)
                    if edges[i, z_row, z_col, 0]:                        
                        for r in range(max(0,z_row-spreading), min(z_row+spreading+1,nrow)):
                            for c in range(max(0,z_col-spreading), min(z_col+spreading+1,ncol)):
                                edges[i, r, c, 0] = edges[i, z_row, z_col, 0]
                    edges[i, z_row, z_col, 1] = _get_edge(my_image, z_row, z_col, 1, 1, 1, -1, min_contrast)
                    if edges[i, z_row, z_col, 1]:
                        for r in range(max(0,z_row-spreading), min(z_row+spreading+1,nrow)):
                            for c in range(max(0,z_col-spreading), min(z_col+spreading+1,ncol)):
                                edges[i, r, c, 1] = edges[i, z_row, z_col, 1]
                    edges[i, z_row, z_col, 2] = _get_edge(my_image, z_row, z_col, 1, 0, 0, -1, min_contrast)
                    if edges[i, z_row, z_col, 2]:
                        for r in range(max(0,z_row-spreading), min(z_row+spreading+1,nrow)):
                            for c in range(max(0,z_col-spreading), min(z_col+spreading+1,ncol)):
                                edges[i, r, c, 2] = edges[i, z_row, z_col, 2]
                    edges[i, z_row, z_col, 3] = _get_edge(my_image, z_row, z_col, 1, -1, -1, -1, min_contrast)
                    if edges[i, z_row, z_col, 3]:
                        for r in range(max(0,z_row-spreading), min(z_row+spreading+1,nrow)):
                            for c in range(max(0,z_col-spreading), min(z_col+spreading+1,ncol)):
                                edges[i, r, c, 3] = edges[i, z_row, z_col, 3]
    
    if binary:
        output = np.concatenate((np.array(edges)==1, np.array(edges)==2), 3)
    else:
        output = np.array(edges)
        
    return output if images.ndim==3 else output[0,]