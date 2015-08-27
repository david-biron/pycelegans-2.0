
#=======================================================================================================================
# Dynamic programming for finding sidelines.
#=======================================================================================================================

#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True


from libc.math cimport atan2
from libc.math cimport M_PI as pi
cdef extern from "math.h":
    double INFINITY
    bint isnan(double)
    double sqrt(double)
    double sin(double)
    double pow(double, int)

import numpy as np
cimport numpy as cnp



cdef int _best_ll_ratios(cnp.ndarray[cnp.uint8_t, ndim=3, mode='c'] edges, double ll_ratio_feature, \
                         double ll_ratio_no_feature, cnp.ndarray[cnp.int32_t, ndim=2] region1, \
                         cnp.ndarray[cnp.float64_t, ndim=1] previous_opt_vals, cnp.ndarray[cnp.int32_t, ndim=2] region2, \
                         cnp.ndarray[cnp.int32_t, ndim=1] opt_choices, cnp.ndarray[cnp.float64_t, ndim=1] opt_vals, \
                         double angle_mod, double B, double C):
    """
    Chooses points in the specified regions, determines the appropriate edge types and orientations,
    and stores pairs with optimal log-likelihood ratios in the provided arrays.
    
    Parameters
    ----------
    edges : array
        Oriented edges maps.
    ll_ratio_feature, ll_ratio_no_feature : double
        Log-likelihood ratios for object vs background.
    region1, region2 : array
        Regions for choosing point1 and point2.
    previous_opt_vals : array
        Optimal values for choices of point1.
    opt_choices : array
        Storage for optimal choices of point1 for each point2.
    opt_vals : array
        Storage for optimal values for choices of point2.
    angle_mod : double
        Model angle.
    B : double
        Penalty for deviation from model angle.
    C : double
        Penalty for deviation from model position.
    """
    
    cdef Py_ssize_t row1, col1, row2, col2
    cdef int drow, dcol
    
    cdef double ll_ratio, posterior, angle, delta_angle
    
    cdef int p1, p2
    
    # find optimal point in region1 for each point in region2
    for p2 in range(region2.shape[0]):
        opt_vals[p2] = -INFINITY
        
        for p1 in range(region1.shape[0]):
            
            # calculate angle between point 1 and 2
            row1 = region1[p1, 0]; col1 = region1[p1, 1]
            row2 = region2[p2, 0]; col2 = region2[p2, 1]
            drow = row2 - row1
            dcol = col2 - col1
            if drow == dcol == 0: continue
            angle = atan2(drow, dcol)
            if angle < 0:
                orientation = 2
                angle = angle + pi
            else:
                orientation = 1
            
            # choose appropriate edge types, orientations and determine ll-ratio between points 
            if (angle <= pi/12):
                ll_ratio = _ll_ratio_along_line(edges, 2, orientation, -1, 0, ll_ratio_feature, ll_ratio_no_feature, \
                                                row1, col1, row2, col2)        
            elif (pi/12 < angle < pi/6):
                ll_ratio = _ll_ratio_along_line(edges, 2, orientation, 3, orientation, ll_ratio_feature, \
                                                ll_ratio_no_feature, row1, col1, row2, col2)
            elif (pi/6 <= angle <= pi/3):
                ll_ratio = _ll_ratio_along_line(edges, 3, orientation, -1, 0, ll_ratio_feature, ll_ratio_no_feature, \
                                                row1, col1, row2, col2)
            elif (pi/3 < angle < 5*pi/12):
                ll_ratio = _ll_ratio_along_line(edges, 3, orientation, 0, 2 if orientation == 1 else 1, \
                                                ll_ratio_feature, ll_ratio_no_feature, row1, col1, row2, col2)
            elif (5*pi/12 <= angle <= 7*pi/12):
                ll_ratio = _ll_ratio_along_line(edges, 0, 2 if orientation == 1 else 1, -1, 0, ll_ratio_feature, \
                                                ll_ratio_no_feature, row1, col1, row2, col2)
            elif (7*pi/12 < angle < 2*pi/3):
                ll_ratio = _ll_ratio_along_line(edges, 0, 2 if orientation == 1 else 1, 1, \
                                                2 if orientation == 1 else 1, ll_ratio_feature, ll_ratio_no_feature, \
                                                row1, col1, row2, col2)
            elif (2*pi/3 <= angle <= 5*pi/6):
                ll_ratio = _ll_ratio_along_line(edges, 1, 2 if orientation == 1 else 1, -1, 0, ll_ratio_feature, \
                                                ll_ratio_no_feature, row1, col1, row2, col2)
            elif (5*pi/6 < angle < 11*pi/12):
                ll_ratio = _ll_ratio_along_line(edges, 1, 2 if orientation == 1 else 1, 2, \
                                                2 if orientation == 1 else 1, ll_ratio_feature, ll_ratio_no_feature, \
                                                row1, col1, row2, col2)
            else: # i.e., elif (11*pi/12 <= angle):
                ll_ratio = _ll_ratio_along_line(edges, 2, 2 if orientation == 1 else 1, -1, 0, ll_ratio_feature, \
                                                ll_ratio_no_feature, row1, col1, row2, col2)
            
            # calculate deviation penalty and determine 'posterior'
            if isnan(angle_mod):
                delta_angle = 0
            else:
                delta_angle = angle_mod - angle
            posterior = previous_opt_vals[p1] + ll_ratio - B * abs(sin(delta_angle)) * sqrt(drow**2 + dcol**2) \
                        - C * abs(<double> p1 / region1.shape[0] - .5)            
            
            # check if current choice of p1 is optimal
            if posterior > opt_vals[p2]:
                opt_vals[p2] = posterior
                opt_choices[p2] = p1

    return 0


cdef double _ll_ratio_along_line(cnp.ndarray[cnp.uint8_t, ndim=3, mode='c'] edges, int edge_type1, \
                                 unsigned char orientation1, int edge_type2, unsigned char orientation2, \
                                 double ll_ratio_feature, double ll_ratio_no_feature, \
                                 int row1, int col1, int row2, int col2):
    """
    Calculates the log-likelihood ratio along a line segment for given edge types and orientations.
    
    Parameters
    ----------
    edges : array
        Oriented edges maps.
    edge_type1, edge_type2 : int
        Edge types to be used; -1 for edge_type2 indicates that only one edge type is used.
    orientation1, orientation2 : int
        Orientations of the two edge types.
    ll_ratio_feature, ll_ratio_no_feature : double
        Log-likelihood ratios for object vs background.
    row1, col1, row2, col2 : int
        Coordinates of point1 and point2.
        
    Returns
    -------
    ll_sum : double
        Log-likelihood ratio along the segment.
    """

    cdef int abs_drow = abs(row2 - row1)
    cdef int abs_dcol = abs(col2 - col1)
    cdef Py_ssize_t tmp
    cdef int offset
    cdef double error = 0
    cdef double derror
    cdef double ll_sum = 0
    cdef int row, col
    
    # go along row / column depending on which distance is larger
    if abs_drow > abs_dcol:
        # make sure row end points are in order
        if row1 > row2:
            tmp = row1
            row1 = row2
            row2 = tmp
            tmp = col1
            col1 = col2
            col2 = tmp
        offset = 1 if col1 < col2 else -1
        derror = <double> abs_dcol / abs_drow
        col = col1
        # loop over all rows
        for row in range(row1, row2+1):
            ll_sum += ll_ratio_feature if edges[row, col, edge_type1] == orientation1 else ll_ratio_no_feature
            if edge_type2 != -1:
                ll_sum += ll_ratio_feature if edges[row, col, edge_type2] == orientation2 else ll_ratio_no_feature
            error = error + derror
            if error >= 0.5:
                col += offset
                error -= 1.0
    else:
        # make sure column end points are in order
        if col1 > col2:
            tmp = col1
            col1 = col2
            col2 = tmp
            tmp = row1
            row1 = row2
            row2 = tmp
        offset = 1 if row1 < row2 else -1
        derror = <double> abs_drow / abs_dcol
        row = row1
        # loop over all columns
        for col in range(col1, col2+1):
            ll_sum += ll_ratio_feature if edges[row, col, edge_type1] == orientation1 else ll_ratio_no_feature
            if edge_type2 != -1:
                ll_sum += ll_ratio_feature if edges[row, col, edge_type2] == orientation2 else ll_ratio_no_feature
            error = error + derror
            if error >= 0.5:
                row += offset
                error -= 1.0
                
    return ll_sum


def dynamic_programming(edges, ll_ratio_feature, ll_ratio_no_feature, regions, angles, B=0, C=0):
    """
    Runs dynamic programming for fine detection (from tail to head to tail, background on the left).
    
    Parameters
    ----------
    edges : array
        Oriented edges maps.
    ll_ratio_feature, ll_ratio_no_feature : double
        Log-likelihood ratios for object vs background.
    regions : list
        Admissible regions.
    angles : list
        Model angles.
    B : double, optional
        Penalty for deviation from model angle.
    C : double, optional
        Penalty for deviation from model position.

        
    Returns
    -------
    opt_vals : array
        Optimal values.
    opt_choices : list
        Optimal point choices.
    """ 

    m = len(regions)
    opt_choices = []
    opt_vals = np.zeros(regions[0].shape[0], dtype=np.float64)
    
    for i in range(1,m):
        opt_choices.append(np.zeros(regions[i].shape[0], dtype=np.int32))
        previous_opt_vals = opt_vals
        opt_vals = np.zeros(regions[i].shape[0], dtype=np.float64)
        _best_ll_ratios(edges, ll_ratio_feature, ll_ratio_no_feature, regions[i-1], previous_opt_vals, \
                        regions[i], opt_choices[i-1], opt_vals, angles[i-1], B, C)
        
    return opt_vals, opt_choices