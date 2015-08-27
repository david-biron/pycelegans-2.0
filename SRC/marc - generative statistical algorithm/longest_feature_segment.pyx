
#=======================================================================================================================
# Search for longest continuous feature segment. 
#=======================================================================================================================

#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True


cimport numpy as np



def longest_segment(np.ndarray[np.uint8_t, ndim=3, mode='c'] features):
    """
    Find the best_length continuous segment of blocks where each one has the feature for the appropriate orientation on.
    The actual score is the length plus the number of neighbors which also have the feature on.  
    
    Parameters
    ----------
    features : array of shape (nbest_row, nbest_col, 4)
        Detected oriented edges.
        
    Returns
    -------
    best_row, best_col : int
        Center of optimal segment. 
    best_orientation : int
        Orientation of optimal segment.
    best_length : int
        Length of optimal segment.
    """

    cdef int nrow = features.shape[0], ncol = features.shape[1]
    cdef int best_score = 0, best_length = 0, best_row, best_col, best_orientation
    cdef int score, length
    cdef int r, c, k, l
                
    # horizontal
    for r in range(1, nrow-1):
        length = 0
        score = 0
        for c in range(ncol):
            if features[r, c, 2]:
                length += 1
                score += 1 + features[r-1, c, 2] + features[r+1, c, 2]
                if score > best_score:
                    best_score = score
                    best_length = length
                    best_row = r
                    best_col = c - length // 2
                    best_orientation = 0
            else:
                length = 0
                score = 0
                
    # vertical
    for c in range(1, ncol-1):
        length = 0
        score = 0
        for r in range(nrow):
            if features[r, c, 0]:
                length += 1
                score += 1 + features[r, c-1, 2] + features[r, c+1, 2]
                if score > best_score:
                    best_score = score
                    best_length = length
                    best_row = r - length // 2
                    best_col = c
                    best_orientation = 2
            else:
                length = 0
                score = 0
                
    # diagonal (NW-SE)
    for k in range(nrow-1):  # bottom left
        length = 0
        score = 0
        for l in range(min(nrow-k-1, ncol-1)):
            if features[k+l, l, 3]:
                length += 1
                score += 1 + features[k+l+1, l, 3] + features[k+l, l+1, 3]
                if score > best_score:
                    best_score = score
                    best_length = length
                    best_row = k + l - length // 2
                    best_col = l - length // 2
                    best_orientation = 1
                else:
                    length = 0
                    score = 0        
    for k in range(1, ncol-1):  # top right
        length = 0
        score = 0
        for l in range(min(ncol-k-1, nrow-1)):
            if features[l, k+l, 3]:
                length += 1
                score += 1 + features[l+1, k+l, 3] + features[l, k+l+1, 3]
                if score > best_score:
                    best_score = score
                    best_length = length
                    best_row = l - length // 2
                    best_col = k + l - length // 2
                    best_orientation = 1
            else:
                length = 0
                score = 0         
                
    # diagonal (NE-SW)
    for k in range(1, ncol):   # top left
        length = 0
        score = 0
        for l in range(0, min(k,nrow-1)):
            if features[l, k-l, 1]:
                length += 1
                score += 1 + features[l+1, k-l, 1] + features[l, k-l-1, 1]
                if score > best_score:
                    best_score = score
                    best_length = length
                    best_row = l - length // 2
                    best_col = k-l + length // 2
                    best_orientation = 3
            else:
                length = 0
                score = 0
    for k in range(ncol, ncol+nrow-2):  # bottom right
        length = 0
        score = 0
        for l in range(k-ncol+1, min(k,nrow-1)):
            if features[l, k-l, 1]:
                length += 1
                score += 1 + features[l+1, k-l, 1] + features[l, k-l-1, 1]
                if score > best_score:
                    best_score = score
                    best_length = length
                    best_row = l - length // 2
                    best_col = k-l + length // 2
                    best_orientation = 3
            else:
                length = 0
                score = 0            
    
    return best_row, best_col, best_orientation, best_length