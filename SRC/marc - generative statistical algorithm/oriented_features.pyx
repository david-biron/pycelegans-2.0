
#=======================================================================================================================
# Oriented features for the body and head of a worm.
#=======================================================================================================================

#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True


import numpy as np
cimport numpy as cnp

from libc.stdlib cimport abs
from libc.math cimport floor, sqrt, fabs
from libcpp cimport bool



def feature_detection(cnp.ndarray[cnp.uint8_t, ndim=3, mode='c'] edges, int width, \
                      cnp.ndarray[cnp.uint8_t, ndim=2, mode='c'] coarse, int max_body_intensity=128, \
                      int min_head_intensity=128, int max_head_intensity=230):
    """
    Finds body and head features in the image by counting the number of edges in the appropriate regions.
    For the body features these regions are squares with side length 2*block_size = width/4.
    
    Parameters
    ----------
    edges : array of shape (nrow, ncol, 4)
        Oriented edges.
    width : int
        Worm width.
    coarse : array of shape (nrow//width, ncol//block_size)
        Coarse gray-scale image.
    
    Returns
    -------
    head_features : array of shape (nrow//width, ncol//width, 4)
        Locations of features like "+-+-".
    body_features : array of shape (nrow//width, ncol//width, 4)
        Locations of features like "+  -". 
    doublebody_features : array of shape (nrow//width, ncol//width, 4)
        Locations of features like "+  0 -" (only if gray-scale is rather dark).
    """
    
    cdef int block_size = int(floor(.25 * width + .5))
    cdef int three_half_block_size = block_size + block_size // 2
    cdef int diag = int(floor(sqrt(0.03125)*width + .5))
    cdef int half_diag = diag // 2
    cdef int nrow = edges.shape[0] // block_size
    cdef int ncol = edges.shape[1] // block_size
    
    # create border
    cdef cnp.ndarray[cnp.uint8_t, ndim=3] my_edges = \
        np.zeros((edges.shape[0]+16*block_size, edges.shape[1]+16*block_size, 4), dtype=np.uint8)
    my_edges[8*block_size:edges.shape[0]+8*block_size, 8*block_size:edges.shape[1]+8*block_size,] = edges
    
    # set thresholds
    cdef int threshold = <int>(1.75 * block_size)
    cdef int half_threshold = threshold // 2
    cdef int sum_threshold_outer = <int>(3.75 * block_size)
    cdef int threshold_inner = <int>(2 * block_size)
    cdef int sum_threshold_total = width # 8 block_size

    # create output buffer    
    cdef cnp.ndarray[cnp.uint8_t, ndim=3] head_features = np.zeros((nrow, ncol, 4), dtype=np.uint8)
    cdef cnp.ndarray[cnp.uint8_t, ndim=3] body_features = np.zeros((nrow, ncol, 4), dtype=np.uint8)
    cdef cnp.ndarray[cnp.uint8_t, ndim=3] doublebody_features = np.zeros((nrow, ncol, 4), dtype=np.uint8)

    ### loop through image and detect features ###
    
    cdef int my_row, row, my_r, r, my_col, col, my_c, c, edge_type, other_row, other_col
    cdef int nr1, nr2, nr3, nr4, nr6
    cdef cnp.uint8_t result    
    cdef bool check_body, check_head

    for my_row in range(nrow):
        for my_col in range(ncol):
            check_body = coarse[my_row, my_col] <= max_body_intensity
            check_head = (coarse[my_row, my_col] >= min_head_intensity and coarse[my_row, my_col] <= max_head_intensity)
            for my_r in range(block_size):
                for my_c in range(block_size):
                    
                    row = (my_row+8)*block_size + my_r
                    col = (my_col+8)*block_size + my_c
                    
                    ### vertical ###

                    nr1 = 0; nr2 = 0; nr3 = 0; nr4 = 0; nr6 = 0
                    
                    if check_body:
                   
                        # BODY
                        # left
                        for r in range(row-block_size, row+block_size+1):
                            for c in range(col-3*block_size, col-block_size):
                                nr1 += (my_edges[r, c, 0] == 1)
                        # right
                        for r in range(row-block_size, row+block_size+1):
                            for c in range(col+block_size+1, col+3*block_size+1):
                                nr2 += (my_edges[r, c, 0] == 2)
                        if (nr1 >= threshold and nr2 >= threshold and nr1+nr2 >= sum_threshold_outer):
                            body_features[my_row, my_col, 0] = True
                          
                        # DOUBLE BODY
                        #NOTE: 4-6-8 corresponds to the full worm width.
                        #However, when parts of the worm meet the total width is typically smaller, hence 3-5-7 is used.
                        #(If modified, placement of the adjoint point has to be modified as well.)
                        other_row = my_row
                        other_col = (my_col*block_size + 2*block_size + three_half_block_size + my_c) // block_size
                        if coarse[other_row, other_col] <= max_body_intensity:
                            # right
                            for r in range(row-block_size, row+block_size+1):
                                for c in range(col+5*block_size+1, col+7*block_size+1):
                                    nr6 += (my_edges[r, c, 0] == 2)
                            if (nr1 >= threshold and nr2 < half_threshold and nr6 >= threshold and nr1+nr6 >= sum_threshold_outer): 
                                doublebody_features[my_row, my_col, 0] = True
                                doublebody_features[other_row, other_col, 0] = True
                                
                    if check_head:
                        
                        # HEAD
                        # left
                        for r in range(row-block_size, row+block_size+1):
                            for c in range(col-three_half_block_size, col):
                                nr4 += (my_edges[r, c, 0] == 2)
                        # right
                        for r in range(row-block_size, row+block_size+1):
                            for c in range(col+1, col+three_half_block_size+1):
                                nr3 += (my_edges[r, c, 0] == 1)
                        if (nr3 >= threshold_inner and nr4 >= threshold_inner and nr1+nr2+nr3+nr4 >= sum_threshold_total):
                            head_features[my_row, my_col, 0] = True
                    
                    ### horizontal ###
                    
                    nr1 = 0; nr2 = 0; nr3 = 0; nr4 = 0; nr6 = 0
                    
                    if check_body:
                        
                        # BODY
                        # top
                        for r in range(row-3*block_size, row-block_size):
                            for c in range(col-block_size, col+block_size+1):
                                nr1 += (my_edges[r, c, 2] == 1)
                        # bottom
                        for r in range(row+block_size+1, row+3*block_size+1):
                            for c in range(col-block_size, col+block_size+1):
                                nr2 += (my_edges[r, c, 2] == 2)
                        if (nr1 >= threshold and nr2 >= threshold and nr1+nr2 >= sum_threshold_outer):
                            body_features[my_row, my_col, 2] = True
    
                        # DOUBLE BODY
                        other_row = (my_row*block_size + 2*block_size + three_half_block_size + my_r) // block_size
                        other_col = my_col
                        if coarse[other_row, other_col] <= max_body_intensity:  
                            # bottom
                            for r in range(row+5*block_size+1, row+7*block_size+1):
                                for c in range(col-block_size, col+block_size+1):
                                    nr6 += (my_edges[r, c, 2] == 2)
                            if (nr1 >= threshold and nr2 < half_threshold and nr6 >= threshold and nr1+nr6 >= sum_threshold_outer): 
                                doublebody_features[my_row, my_col, 2] = True
                                doublebody_features[other_row, other_col, 2] = True
                                
                    if check_head:
                        
                        # HEAD
                        # top
                        for r in range(row-three_half_block_size, row):
                            for c in range(col-block_size, col+block_size+1):
                                nr4 += (my_edges[r, c, 2] == 2)
                        # bottom
                        for r in range(row+1, row+three_half_block_size+1):
                            for c in range(col-block_size, col+block_size+1):
                                nr3 += (my_edges[r, c, 2] == 1)
                        if (nr3 >= threshold_inner and nr4 >= threshold_inner and nr1+nr2+nr3+nr4 >= sum_threshold_total):
                            head_features[my_row, my_col, 2] = True
                                
                    ### diagonal (NW-SE) ###
                    
                    nr1 = 0; nr2 = 0; nr3 = 0; nr4 = 0; nr6 = 0
                    
                    if check_body:
                    
                        # BODY
                        # top left
                        for r in range(row-4*diag, row):
                            for c in range(col-4*diag+abs(row-2*diag-r), <int>(col-fabs(row-2*diag-.5-r)+.5)):
                                nr1 += (my_edges[r, c, 1] == 1)
                        # bottom right
                        for r in range(row, row+4*diag+1):
                            for c in range(<int>(col+fabs(row+2*diag-r+.5)+.5), col+4*diag-abs(row+2*diag-r)+1):
                                nr2 += (my_edges[r, c, 1] == 2)
                        if (nr1 >= threshold and nr2 >= threshold and nr1+nr2 >= sum_threshold_outer):
                            body_features[my_row, my_col, 1] = True
                
                        # DOUBLE BODY
                        other_row = (my_row*block_size + 3*diag + half_diag + my_r) // block_size
                        other_col = (my_col*block_size + 3*diag + half_diag + my_c) // block_size
                        if coarse[other_row, other_col] <= max_body_intensity:
                            # bottom right
                            for r in range(row+3*diag, row+7*diag+1):
                                for c in range(<int>(col+3*diag+fabs(row+5*diag-r+.5)+.5), col+7*diag-abs(row+5*diag-r)+1):
                                    nr6 += (my_edges[r, c, 1] == 2)
                            if (nr1 >= threshold and nr2 < half_threshold and nr6 >= threshold and nr1+nr6 >= sum_threshold_outer): 
                                doublebody_features[my_row, my_col, 1] = True
                                doublebody_features[other_row, other_col, 1] = True
                                
                    if check_head:
                        
                        # HEAD
                        # top left
                        for r in range(row-3*diag, row+diag):
                            for c in range(col-2*diag-half_diag+abs(row-half_diag-r), <int>(col+diag-fabs(row-diag-.5-r)+.5)):
                                nr4 += (my_edges[r, c, 1] == 2)
                        # bottom right
                        for r in range(row-diag, row+3*diag+1):
                            for c in range(<int>(col-diag+fabs(row+diag-r+.5)+.5), col+2*diag+half_diag-abs(row+half_diag-r)+1):
                                nr3 += (my_edges[r, c, 1] == 1)
                        if (nr3 >= threshold_inner and nr4 >= threshold_inner and nr1+nr2+nr3+nr4 >= sum_threshold_total):
                            head_features[my_row, my_col, 1] = True
                
                    ### diagonal (NE-SW) ###
                    
                    nr1 = 0; nr2 = 0; nr3 = 0; nr4 = 0; nr6 = 0
                    
                    if check_body:
                    
                        # BODY
                        # top right
                        for r in range(row-4*diag, row):
                            for c in range(<int>(col+fabs(row-2*diag-r-.5)+.5), col+4*diag-abs(row-2*diag-r)+1):
                                nr1 += (my_edges[r, c, 3] == 1)
                        # bottom left
                        for r in range(row, row+4*diag+1):
                            for c in range(col-4*diag+abs(row+2*diag-r), <int>(col-fabs(row+2*diag-r+.5)+.5)):
                                nr2 += (my_edges[r, c, 3] == 2)
                        if (nr1 >= threshold and nr2 >= threshold and nr1+nr2 >= sum_threshold_outer):
                            body_features[my_row, my_col, 3] = True
                                
                        # DOUBLE BODY
                        other_row = (my_row*block_size + 3*diag + half_diag + my_r) // block_size
                        other_col = (my_col*block_size - 3*diag - half_diag + my_c) // block_size
                        if coarse[other_row, other_col] < max_body_intensity: 
                            # bottom left
                            for r in range(row+3*diag, row+7*diag+1):
                                for c in range(col-7*diag+abs(row+5*diag-r), <int>(col-3*diag-fabs(row+5*diag-r+.5)+.5)):
                                    nr6 += (my_edges[r, c, 3] == 2)
                            if (nr1 >= threshold and nr2 < half_threshold and nr6 >= threshold and nr1+nr6 >= sum_threshold_outer): 
                                doublebody_features[my_row, my_col, 3] = True
                                doublebody_features[other_row, other_col, 3] = True
                                
                    if check_head:
                        
                        # HEAD
                        # top right
                        for r in range(row-3*diag, row+diag):
                            for c in range(<int>(col-diag+fabs(row-diag-r-.5)+.5), col+2*diag+half_diag-abs(row-half_diag-r)+1):
                                nr4 += (my_edges[r, c, 3] == 2)
                        # bottom left
                        for r in range(row-diag, row+3*diag+1):
                            for c in range(col-2*diag-half_diag+abs(row+half_diag-r), <int>(col+diag-fabs(row+diag-r+.5)+.5)):
                                nr3 += (my_edges[r, c, 3] == 1)
                        if (nr3 >= threshold_inner and nr4 >= threshold_inner and nr1+nr2+nr3+nr4 >= sum_threshold_total):
                            head_features[my_row, my_col, 3] = True

    return head_features, body_features, doublebody_features