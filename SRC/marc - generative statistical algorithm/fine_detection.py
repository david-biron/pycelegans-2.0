
#=======================================================================================================================
# Fine detection of the worm posture.
#=======================================================================================================================

import numpy as np
from scipy import interpolate

import pyximport; pyximport.install()
import dynamic_programming as DP



PI = np.pi
PI_HALF = np.pi / 2
PI_FOURTH = np.pi / 4
PI_EIGHTS = np.pi / 8



def _add_region(regions, point, angle, im_shape, width, start=.2, end=.8):
    """
    Adds a region, starting from the given point going into the specified direction, to the provided list.
    The region is truncated such that the image is not left.
    
    Parameters
    ----------
    regions : list
        List the region is appended to.
    point : integer pair
        Starting point.
    angle: float
        Direction.
    width : float
        Worm width.
    im_shape : integer tuple
        Shape of the image.
    start, end : double, optional
        Distance into direction.
    """
    
    start = start * width
    end = end * width
    
    region = np.vstack((point[0]+np.sin(angle)*np.arange(start,end), point[1]+np.cos(angle)*np.arange(start,end)))
    region = np.asarray(np.transpose(region), dtype = np.int32)
    #NOTE: some points might be contained multiple times
    
    # truncate region to image size
    while region.shape[0] > 0 and (region[-1,0] < 0 or region[-1,0] >= im_shape[0] \
                                   or region[-1,1] < 0 or region[-1,1] >= im_shape[1]): 
        region = region[:-1,]
    if region.shape[0] == 0:
        region = np.asarray(point, dtype = np.int32).reshape(1,2)
    
    regions.append(region)


def _spline(points, s=None, nr_of_points=None, do_angles=False):
    """
    Calculates a cubic spline for the given sequence of points.
    
    Parameters
    ----------
    points : list of pairs; or pair of lists
        Points to be smoothed.
    s : float, optional
        Smoothing parameter.
    nr_of_points : integer, optional
        Number of points at which the spline is evaluated.
    do_angles : boolean, optional
        Whether angles for the smoothed points and line segments should be calculated.
        
    Returns
    -------
    points : list of pairs
        Smoothed points.
    angles_points : array
        Angles at smoothed points.
    angles_segments : array
        Angles between smoothed points.
    """
    
    if isinstance(points, tuple):
        x = points[0]
        y = points[1]
    else:
        x, y = zip(*points)
        
    n = len(x) if nr_of_points is None else nr_of_points    
        
    # calculate parametric cubic spline
    if s is not None: s = s * n
    tck, __ = interpolate.splprep([x,y], s=s)    
    u = np.linspace(0, 1, n)
    smooth_points = interpolate.splev(u, tck)
    
    if do_angles:
        # calculate angles
        u_der = np.linspace(0, 1, 2*n-1)
        derivatives = interpolate.splev(u_der, tck, der=1)        
        return list(zip(np.array(smooth_points[0],dtype=np.int), np.array(smooth_points[1],dtype=np.int))), \
               np.arctan2(derivatives[0][::2],derivatives[1][::2]), \
               np.arctan2(derivatives[0][1::2],derivatives[1][1::2])
    else:
        return list(zip(np.array(smooth_points[0], dtype=np.int), np.array(smooth_points[1], dtype=np.int)))
    
    
def detect(edges, points, density, width, B, pobj = .7, pbg = .2):
    """
    Finds the sidelines and an improved midline.
    
    Parameters
    ----------
    edges : array
        Oriented edge maps.
    points : list of pairs
        Points from coarse instantiation (tail to head).
    density : integer
        Density of features.
    width : integer
        Width of worm.
    B : double
        Penalty for deviation from model angle (smoothed coarse detection).
    pobj, pbg : double, optional
        Edge probability on the sidelines and on background.
    
    Returns
    -------
    points : list of pairs
        Smoothed worm points.
    regions : list
        Admissible regions for sideline points.
    angles_segments : array
        Angles between smoothed points.
    side : list
        Sideline.
    smooth_side : list
        Smoothed sideline.
    smooth_midline : list
        Smoothed midline.
    ll_edges : float
        Log-likelihood for edges.
    """

    n = len(points)
    
    # rescale and smooth worm points
    x, y = zip(*points)
    x = [density * __ + density//2 for __ in x]
    y = [density * __ + density//2 for __ in y]
    #points, angles_points, angles_segments = _spline((x,y), s=10, do_angles=True)
    #points, angles_points, angles_segments = _spline((x,y), s=5, do_angles=True)
    points, angles_points, angles_segments = _spline((x,y), s=3, do_angles=True)
    
    # choose admissible regions for dynamic programming
    im_shape = edges.shape[:2]
    regions = []

    # - tail
    _add_region(regions, points[0], angles_points[0] - PI_HALF - PI_FOURTH, im_shape, width, start=0, end=1)
    
    # - tail to head
    for i in range(1, n):
        _add_region(regions, points[i], angles_points[i] - PI_HALF, im_shape, width)
        
    # - head
    _add_region(regions, points[n-1], angles_points[n-1] - PI_FOURTH, im_shape, width)
    _add_region(regions, points[n-1], angles_points[n-1] - PI_EIGHTS, im_shape, width, start=0, end=1)
    _add_region(regions, points[n-1], angles_points[n-1], im_shape, width, start=0, end=1)
    _add_region(regions, points[n-1], angles_points[n-1] + PI_EIGHTS, im_shape, width, start=0, end=1)
    _add_region(regions, points[n-1], angles_points[n-1] + PI_FOURTH, im_shape, width)
    
    # - head to tail
    for i in range(n-1, 0, -1):
        _add_region(regions, points[i], angles_points[i] + PI_HALF, im_shape, width)
    
    # - tail
    _add_region(regions, points[0], angles_points[0] + PI_HALF + PI_FOURTH, im_shape, width, start=0, end=1)
    
    # determine model angles for line segments (no prior for head and tail)
    angles_segments = np.hstack((angles_segments, np.ones(6)*np.NAN, \
                                np.where(angles_segments[::-1] >= 0, angles_segments[::-1] - PI, angles_segments[::-1] + PI)))
    angles_segments[:1] = angles_segments[-1:] = np.NAN
    
    # set edge probabilities
    log_feat = np.log(pobj) - np.log(pbg)
    log_no_feat = np.log(1-pobj) - np.log(1-pbg)
    
    # run dynamic programming and store side lines
    m = len(regions)
    opt_vals, opt_choices = DP.dynamic_programming(edges, log_feat, log_no_feat, regions, angles_segments, B)
    opt_choice = np.argmax(opt_vals)
    ll_edges = opt_vals[opt_choice]
    p2 = regions[m-1][opt_choice, :]
    side = [p2]
    for j in range(m-2, -1, -1):
        opt_choice = opt_choices[j][opt_choice]
        p1 = regions[j][opt_choice, :]
        side.append(p1)
        p2 = p1
    
    # smooth side lines
    smooth_side = _spline(side, nr_of_points=2*150)
    
    # create midline
    midline = []
    for j in range(n+1):
        left = side[j]
        right = side[-(j+1)]
        mid = ((left[0]+right[0])//2, (left[1]+right[1])//2)
        midline.append(mid)
    if midline[-1] == midline[-2]: del midline[-1] # otherwise, spline error
        
    # smooth midline (except tip of tail and head)
    smooth_midline = _spline(midline, nr_of_points=150)
    smooth_midline[0] = midline[0]
    smooth_midline[-1] = midline[-1]
        
    return points, regions, side, smooth_side, smooth_midline, ll_edges