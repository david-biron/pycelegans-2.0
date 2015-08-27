'''This module performs the segmentation and identification of worm.'''

import libcelegans 

def doseg(max_image,OrigImage,THRESH,THRESH_DIFF): 
#  need to make THRESH_DIFF a command line option
    DiffImage = ( (max_image - OrigImage) > THRESH_DIFF ) * OrigImage
    WormBinArr = libcelegans.fill_object(
                     libcelegans.get_biggest_object(
                      libcelegans.get_thresh(DiffImage,thresh=THRESH)))
    WormArr = WormBinArr * OrigImage

    return (DiffImage,WormBinArr,WormArr)
