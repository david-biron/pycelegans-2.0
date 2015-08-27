
#=======================================================================================================================
# Worm detection (Marc).
#=======================================================================================================================


import os
import scipy.ndimage
from skimage import io
import numpy as np

from marc import fine_detection

import pyximport; pyximport.install()
from marc import edge_detection
from marc import coarse_image
from marc import oriented_features
from marc import longest_feature_segment
from marc import coarse_detection




def main(IP, max_image, worm_width, worm_length):

    # hack to get image over to marc's code
    WORM_AREA_DILATION_ITERATIONS=35
    '''image = IP.OrigImage * scipy.ndimage.morphology.binary_dilation(
            IP.WormBinArr,iterations=WORM_AREA_DILATION_ITERATIONS)'''
    worm_mask = scipy.ndimage.morphology.binary_dilation(IP.WormBinArr,iterations=WORM_AREA_DILATION_ITERATIONS)

    # get a box
    temp = np.argwhere(worm_mask) 
    (ystart, xstart), (ystop, xstop) = temp.min(0), temp.max(0) + 1 
    offsets = {'ystart':ystart,'xstart':xstart,'ystop':ystop,'xstop':xstop}
    trimimage = IP.OrigImage[ystart:ystop, xstart:xstop]
    
    scipy.misc.imsave((os.path.join(IP.scratchdir,IP.NamePrefix+'.handoff.png')),trimimage)
    image = io.imread(os.path.join(IP.scratchdir,IP.NamePrefix+'.handoff.png'))

#---------------------------------------------------------------------------------------------------------------------- 

    bg_image = np.ascontiguousarray(max_image[ystart:ystop, xstart:xstop])

    # detection parameters 
    worm_width = 20  # 20 for young worms, 25 for adults
    worm_length = 1  # 1 or 1.25 or 0.75
    head_lam = 12.5 * worm_length  # expected head length (number of features)
    head_A = .15  # penalty for length violation
    worm_lam = 70 * worm_length  # expected worm length (number of features), for search mode 2
    worm_A = .075  # penalty for length violation
    B = 1.5  # penalty for angle violation, fine search
    
    # detect worm features    
    edges = edge_detection.get_edges(image, min_contrast=25)
    subtract_image = np.where(image > bg_image, 255, 255-bg_image+image)
    bg_edges = edge_detection.get_edges(bg_image, min_contrast=25)
    nobg_edges = np.where(edges == bg_edges, 0, edges)
    density = np.int(np.round(.25*worm_width))
    coarse = coarse_image.coarse_image(subtract_image, density)
    head_features, body_features, doublebody_features = \
        oriented_features.feature_detection(nobg_edges, worm_width, coarse)
    worm_features = body_features + head_features + doublebody_features
    
    ### run coarse detection ###
    
    res = coarse_detection.detect_worm1(np.ascontiguousarray(head_features), head_lam, head_A, \
                                        np.ascontiguousarray(worm_features), worm_lam, worm_A)
    
    if res is not None:
        coarse_points, log_posterior, __ = res
    else:
        print('No head features found!')

    # run fine detection
    worm_points, regions, side, smooth_side, smooth_midline, ll_edges = \
        fine_detection.detect(edges, coarse_points, density, worm_width, B)
    
#---------------------------------------------------------------------------------------------------------------------- 

    # smooth_side = list of (x,y) 
    # smooth_midline = list of (x,y)

    return (smooth_side,smooth_midline,offsets,log_posterior)
