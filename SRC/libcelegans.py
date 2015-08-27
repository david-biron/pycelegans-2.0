########################################################################
#  Author: Nicholas Labello | labello@uchicago.edu                             
#  The University of Chicago Research Computing Center                         
#                                                                              
#  Library of functions and associated calls for operating on images
#  captured of C. elegans in Biron lab.  
########################################################################
'''
=================
Module pycelegans 
=================

pycelegans is a collection of tools built from gray-scale and binary
morphological functions intended to address image segmentation and 
image analysis problems involving the study of C. elegans.
'''

from collections import namedtuple
import math
import os

import numpy 
import scipy.ndimage 
import scipy.misc    
import skimage.graph
import pymorph

Point = namedtuple('Point', 'x y')
WormBorderPoint = namedtuple('WormBorderPoint', 'x y score')
X = 0 
Y = 1
HIGH_COST = 999999

def get_image(file_pathname):
    '''
* Usage 
    NewArr = get_image(file_pathname)
* Input
    file_pathname: Full path and name to the image on your file system
* Output
    NewArr: the intensity values of the image loaded as integers in a Numpy array 
* Description
    Thin wrapper around scipy.misc.imread.
    '''
    image = scipy.misc.imread(file_pathname,flatten=True)
    return image

def get_thresh(Arr,thresh):
    '''
* Usage
    NewArr = get_thresh(Arr,thresh)
* Input 
    Arr: greyscale image stored as Numpy Array
    thresh: integer that should be between 0-255 is Arr is generated from 8-bit grayscale
* Output
    NewArr: image with original Arr values where Arr > thresh.   
* Description
  thresh refers to the pixel intensity.  Assuming 8-bit grayscale values between
  0 - 255 are logical.  
    '''
    thresh_image = Arr * (Arr < thresh)
    return thresh_image

def get_border(Arr):
    '''
* Usage
    NewArr = get_border(Arr)
* Input 
    Arr: Binary Image stored as Numpy Array
* Output
    NewArr: Binary Image
* Description
    Intended for use with a binary array that contains only one object.
    Takes a binary array, returns a binary array of the same dimensions.  
    The returned array will contain the one-pixel thick border of the
    object in the array passed in.  Array passed in should have one object
    only:  the worm. '''
    return (Arr - scipy.ndimage.binary_erosion(Arr,iterations=1))
    
def fill_object(Arr,i=2): 
    '''
* Usage
    NewArr = get_border(Arr,i=[INTEGER])
* Input 
    Arr: Binary Image stored as Numpy Array
* Output
    NewArr: Binary Image
* Description
    Intended for use with a binary array that contains only one object.
    Takes a binary array, returns a binary array of the same dimensions.  
    The returned array has been dilated by *i* iterations, had a 
    binary_fill_holes morphological operation applied, and then been 
    eroded by *i* iterations.  As a result any holes that may have resulted
    from initial segmentation are filled.  Holes are most common in the 
    lighter regions of the interior of the worm (head/neck region).  
    '''
    filled_object = scipy.ndimage.morphology.binary_dilation(Arr,iterations=i)
    filled_object = scipy.ndimage.morphology.binary_fill_holes(filled_object)
    filled_object = scipy.ndimage.morphology.binary_erosion(filled_object,iterations=i)
    return filled_object

def get_biggest_object(Arr):
    '''
* Usage
    NewArr = get_border(Arr)
* Input 
    Arr: Binary Image stored as Numpy Array
* Output
    NewArr: Binary Image of largest object
* Description
    Takes a binary array, returns a binary array of the same dimensions.  
    The returned array will contain only the largest single object that 
    could be identified in the image.  This is ALWAYS the worm if the
    background image has been correctly substracted. 
    '''
    labeled_array,number_of_features = scipy.ndimage.label(Arr)
# get the size of each object.  biggest *should* be worm, so get the 
# index of the largest object.
    array_of_object_sizes = scipy.ndimage.sum(Arr,labeled_array,xrange(1,number_of_features+1))
    largest_object_index = numpy.argmax(array_of_object_sizes)
    largest_object_label = largest_object_index + 1
    largest_object_binary_array = (labeled_array == largest_object_label)
    return largest_object_binary_array

def get_point_list(Arr):
    '''
* Usage
    ListOfPoints = get_point_list(Arr)
* Input
    Arr: Binary Image stored as Numpy Array
* Output
    ListOfPoints
* Description
    Takes a binary array, returns a Python list of the True or Nonzero
    elements in the array.  Each item in the list is an (x,y) tuple'''
    temp = list(numpy.transpose(numpy.nonzero(Arr)))
    PointList = [Point(item[0],item[1]) for item in temp]
    return PointList

def get_midline(SideOnePath,SideTwoPath,jrange=12):
    '''
* Usage
    ListOfPoints = get_midline(SideOnePath,SideTwoPath)
* Input
    SideOnePath and SideTwoPath: Both are lists of ordered pairs that represent the 
    a path through the perimeter pixels (worm border) of the array, from head to tail.  
    len(SideOnePath) must == len(SideTwoPath).
    jrange: integer
* Output
    ListOfPoints:  a list of :Point:s
* Description
    This function finds the midline of the worm based on the two ordered paths
    passed in.  Each point in SideOnePath is paired to a point in SideTwoPath.  
    Point *i* in SideOnePath is compared to point [i-jrange, i-jrange+1, ... i+jrange]
    points on SideTwoPath.  The closest points in two-dimensional space is used as the
    match.  Midline-Point[i] is taken as the average of these two points.  

    **Note**:  A larger jrange will sample more points and generally be more accurate,
    though exceptions are possible.  Generally this should be based on the number of points
    that is chosen to represent the length of the of worm, len(SideOnePath) or 
    len(SideTwoPath).  jrange of 5-20% of this length is probably reasonable though experimentation
    may be necessary.   
    '''
    LenSideOnePath = len(SideOnePath)
    LenSideTwoPath = len(SideTwoPath)
    offset = 0

    MidPointList = []
    for i in xrange(1,LenSideOnePath-1):     # 1 to -1 excludes the head and tail
        dscore = 9999
        S1Point = Point(SideOnePath[i][X],SideOnePath[i][Y])

        jmin = -offset + i - jrange
        if jmin < 1: jmin = 1
        jmax = offset + i + jrange
        if jmax > LenSideTwoPath-1: jmax = LenSideTwoPath-1
        for j in xrange(jmin,jmax): # 1 to -1 excludes the head and tail
            S2Point = Point(SideTwoPath[j][X],SideTwoPath[j][Y])
            distance = get_distance(S1Point,S2Point)
            if distance < dscore:
                offset = abs(i-j)
                dscore = distance
                MidPoint = Point( (S1Point.x + S2Point.x)/2.0 , (S1Point.y + S2Point.y)/2.0)
                # if midpoints are being identified outside of worm body, add
                # check for whether MidPoint falls in or out here, by bool comparison
                # of MidPoint to WormBinArr.  If fail, simply more forward as if
                # distance was > dscore 
        MidPointList.append(MidPoint)
    return MidPointList 

def get_side_paths(HeadArr,TailArr,BorderRoute):

    headpoint = get_point_list(HeadArr)[0]
    tailpoint = get_point_list(TailArr)[0]
    CostArray = HeadArr + TailArr + HIGH_COST #everything is very costly 
    for point in BorderRoute:
        CostArray[point[X],point[Y]] = 1  # border points are cheap 

    SideOne = skimage.graph.route_through_array(CostArray,headpoint,tailpoint,fully_connected=True)[0]
    SideOneArr = pointlist_to_array2(SideOne,HeadArr)
    ExtraCost = SideOneArr*HIGH_COST # make the first route more costly
    CostArray = CostArray + ExtraCost
    SideTwo = skimage.graph.route_through_array(CostArray,headpoint,tailpoint,fully_connected=True)[0]
    SideTwoArr = pointlist_to_array2(SideTwo,HeadArr)
    return (SideOne,SideTwo,SideOneArr,SideTwoArr) 

def get_border_endpoint(FileName,BorderArr,WormBinArr,SimpleBorderRoute,step_size):
    PointsAndScore = get_dist_step_points(SimpleBorderRoute,step_size)
    LenPointsAndScore = len(PointsAndScore)
#  Sort by score
    PointsAndScoreSorted = sorted(PointsAndScore, key=lambda PointsAndScore: PointsAndScore.score)
    tail = PointsAndScoreSorted[0] # [0] = lowest score = sharpest point (tail)
    tailindex = PointsAndScore.index(tail)
    Tail = pointlist_to_array( (tail,),numpy.zeros(numpy.shape(BorderArr)))

#  Delete a few points around border, near tail, since head will be quite a ways from tail
    delrange = 15 # make commandline option?
    delstart = tailindex - delrange
    delend = tailindex + delrange
    if delstart < 0: 
        delstart = LenPointsAndScore + delstart
        PointsAndScore = PointsAndScore[delend:delstart]

    if delend >= LenPointsAndScore:
        delend = delend - LenPointsAndScore
        PointsAndScore = PointsAndScore[delend:delstart]
    else:
        del PointsAndScore[delstart:delend]

#  Sort again.  This time the lowest score will be the head.
    PointsAndScoreSorted = sorted(PointsAndScore, key=lambda PointsAndScore: PointsAndScore.score) # sort by score 
    head = PointsAndScoreSorted.pop(0) # next lowest score = head
    Head = pointlist_to_array( (head,),numpy.zeros(numpy.shape(BorderArr)))
    return (Head,Tail,head,tail)

def get_dist_step_points(OrderedPoints,step_size):
#  get physical distance between the two points to either side
#  of the point being studied.  this distance will always correlate to the 
#  angle, and is a less expensive calculation.
    NOrderedPoints = len(OrderedPoints)
    PointsAndScore = []
    for i in xrange(0,len(OrderedPoints)):
        p0 = Point(OrderedPoints[i][0],OrderedPoints[i][1])
        p1 = Point(OrderedPoints[i-step_size][0],OrderedPoints[i-step_size][1])
# Handle the case where points at the end of the list must wrap around to the beginning to 
# find their forward point for angle/distance calculation
        if (i + step_size >= NOrderedPoints): 
            pstep_size = step_size - NOrderedPoints   
        else:
            pstep_size = step_size
        p2 = Point(OrderedPoints[i+pstep_size][X],OrderedPoints[i+pstep_size][Y])
        d12 = get_distance(p1,p2)
        PointsAndScore.append(WormBorderPoint(p0.x,p0.y,d12))
    return PointsAndScore 

def get_simple_border_route(BorderRoute,Nworm_divisions):
    Npoints_in_border = len(BorderRoute)
    length_of_worm = int(Npoints_in_border/2)
    step_size = int(length_of_worm / Nworm_divisions)

    if step_size < 1:   # watch for a possible divide by 0 error
        step_size = 1

    SimpleBorderRoute = BorderRoute[0::step_size]
    return SimpleBorderRoute

def get_simple_side_route(SideRoute,Nworm_divisions):
    Npoints_in_side = len(SideRoute)
    length_of_worm = Npoints_in_side
    step_size = int(length_of_worm / Nworm_divisions)
    try:
        SimpleSideRoute = SideRoute[0::step_size]
    except ValueError:
        SimpleSideRoute = SideRoute[:]
    
    return SimpleSideRoute

def pointlist_to_array(PointList,Arr):
#  Simply convert a list of points to be True on the array that is passed in
    NewArr = numpy.copy(Arr)
    for point in PointList:
        NewArr[point.x,point.y] = 1
    return NewArr 

def pointlist_to_array2(PointList,Arr):
    NewArr = numpy.copy(Arr)
    for point in PointList:
        NewArr[point[X],point[Y]] = 1
    return NewArr

def get_border_route(BorderArr):
    BorderRoute = get_cheapest_route(BorderArr,seed=0)
#  In 5-10% of pixels a short route can be found around the high-cost first generations neighbors of the pixel
#  seed, defeating the purpose of the routine which should find the long path, all the way around the border.
#  Simply Use a different pixel seed and try again.
    s = 0
    while len(BorderRoute) < 150:  #150 is specific to worm/resolution.  Need to make option.
        s += 2 # just moves around the border in increments s to find a good starting point 
        BorderRoute = get_cheapest_route(BorderArr,seed=s)
    return BorderRoute

def get_cheapest_route(BorderArr,seed):
    PointList = get_point_list(BorderArr)
    P0 = PointList[seed]
    P0_1Neighbors = get_absolute_pixel_neighbors(P0,k=1)  # first gen neighbors (8 connected)
    P0_1_2Neighbors = get_absolute_pixel_neighbors(P0,k=2) # 1st and 2nd gen neighbors
    P0_2Neighbors = list( set(P0_1_2Neighbors) - set(P0_1Neighbors) ) # 2nd gen neighbors only
    P0_1TrueNeighbors = remove_false_points(BorderArr,P0_1Neighbors) #1st gen True neighbors
    P0_2TrueNeighbors = remove_false_points(BorderArr,P0_2Neighbors) #2nd gen True neighbors

# grab a couple of nearby points on opposite sides of the divide.
    start = P0_2TrueNeighbors[0]
    end = P0_2TrueNeighbors[-1]

# Make the border discontinuous, so that a "start" and "end"  can be fed to
# scikits-route_through_arry
# to do the path finding, (MUCH faster than a naieve implementation).
    BorderArrayDisCont = numpy.copy(BorderArr)
    BorderArrayDisCont[PointList[seed].x,PointList[seed].y] = False
    for point in P0_1TrueNeighbors:
        BorderArrayDisCont[point[X],point[Y]] == False

# make all False points very costly for the path finding
    CostArray = (BorderArrayDisCont == 0) 
    CostArray = CostArray * HIGH_COST

# remove any duplicate points
# get the longpath through the high cost array
    BorderRoute = skimage.graph.route_through_array(CostArray,start,end,fully_connected=True)
    return BorderRoute[0]

def get_distance(p1,p2):
    '''
    * Input p1,p2: Point class
    * Output floating point value, two dimensional distance between p1 and p2
    * Description: Returns the distance between two pixels  
    '''
    return math.sqrt( (p1.x - p2.x)**2 + (p1.y - p2.y)**2)
 
def remove_false_points(Arr,PointList):
    TruePointsList = []
    for p in PointList:
        if Arr[p[X],p[Y]] == True:
            TruePointsList.append(p)
    return TruePointsList
    
def find_true_neighbors(Arr,point):
    '''
* Usage 
    ListOfPoints = find_true_neighbors(Arr,point)
* Input
    Arr: Numpy Array, point: Point class of coordinates in the array
* Output 
    a list of **Point**s
* Description
    Return a list of points that contains all of the True neighbors that exist in the
    binary array for the point passed in.
    '''
    neighbors = get_relative_pixel_neighbors(1) 
    TrueNeighbors = []
    for neigh in neighbors:
        if Arr[point.x+neigh.x,point.y+neigh.y] == True:
            TrueNeighbors.append( Point(point.x+neigh.x,point.y+neigh.y) )
    return TrueNeighbors

def get_relative_pixel_neighbors(k):
    '''
* Usage
    ListOfPoints = get_relative_pixel_neighbors(Arr)
* Input
    k: Any positive integer 
* Output
    ListOfPoints:  a list of **Point**s
* Description
    Returns the transforms necessary to generate the *k*th generation of 
    neighbors around a pixel.  For example, if k = 1, ListOfPoints returns
    the transforms [(-1,-1), (-1,0), ... (1,1)] necessary to get the 8 
    neighboring connected pixels.  (0,0) is excluded. 
    '''
    x = []
    for i in xrange(-k,k+1):
        for j in xrange(-k,k+1):
            if (i == j == 0 ):
                pass
            else:
                x.append(Point(i,j))
    return x

def get_absolute_pixel_neighbors(P,k):
    '''
* Usage
    ListOfPoints = get_absolute_pixel_neighbors(Arr)
* Input
    k: Any positive integer 
* Output
    ListOfPoints:  a list of :Point:s
* Description
    Returns the the *k*th generation of neighbors around a pixel. k = 1 
    returns the 8-connected neighbors.  k = 2 returns the union of k = 1
    and the 8-connected neighbors of all points in the k = 1 set.  
    '''
    relative_neighbors = get_relative_pixel_neighbors(k)
    absolute_neighbors = []
    for rn in relative_neighbors:
        absolute_neighbors.append( (P.x+rn.x,P.y+rn.y) )
    return absolute_neighbors

def perform_head_tail_intensity_test(HeadArr,TailArr,WormArr):
    '''
    * Usage
        bool_result = perform_head_tail_intensity_test(HeadArr,TailArr,WormArr)
    * Input
        Head, Tail, and Worm Arrays
    * Output
        True or False. If Head region is more intense, True, test passed.
    * Description 
        Check regions in worm around what has been identified as Head and Tail.
        This test dilates the head and tail points and averages the intensity of the
        points that fall within the worm body. Since the Head is brighter 
        than the tail, it will both have an average higher intensity. 
        Thus, the value computed for the Head should be greater.  If not, 
        the identification fails this test. 
    '''
    DILATIONS_ITERATIONS = 20
    HeadDilate = scipy.ndimage.morphology.binary_dilation(HeadArr,iterations=DILATIONS_ITERATIONS)
    TailDilate = scipy.ndimage.morphology.binary_dilation(TailArr,iterations=DILATIONS_ITERATIONS)
    HeadRegion = HeadDilate * WormArr
    TailRegion = TailDilate * WormArr
    headmean = scipy.ndimage.mean(HeadRegion)
    tailmean = scipy.ndimage.mean(TailRegion)
    if headmean > tailmean: return True
    else: return False

def perform_head_tail_volume_test(HeadArr,TailArr,WormArr):
    '''
    * Usage
        bool_result = perform_head_tail_volume_test(HeadArr,TailArr,WormArr)
    * Input
        Head, Tail, and Worm Arrays
    * Output
        True or False. If Head region has a greater volume, test passed.
    * Description 
        Check regions in worm around what has been identified as Head and Tail.
        This test dilates the head and tail points and counts the number of 
        points that fall within the worm body. Since the Head is thicker 
        than the tail, it will it should have more pixels. 
        Thus, the value computed for the Head should be greater.  If not, 
        the identification fails this test. 
    '''
    DILATIONS_ITERATIONS = 20
    HeadDilate = scipy.ndimage.morphology.binary_dilation(HeadArr,iterations=DILATIONS_ITERATIONS)
    TailDilate = scipy.ndimage.morphology.binary_dilation(TailArr,iterations=DILATIONS_ITERATIONS)
    HeadRegion = HeadDilate * WormArr
    TailRegion = TailDilate * WormArr
    headvol = numpy.count_nonzero(HeadRegion)
    tailvol = numpy.count_nonzero(TailRegion)
    if headvol > tailvol: return True
    else: return False

def perform_loop_test(SideOnePath,SideTwoPath):
    #  NEED TO MAKE RATIO VALUE TUNABLE AND ADD DOCUMENTATION TO THIS
    #  FUNCTION
    cutoff = 0.6
    L1 = len(SideOnePath)
    L2 = len(SideTwoPath)
    a = max(L1,L2)
    b = min(L1,L2)
    ratio = b/(float(a))
    if ratio <= cutoff: 
        is_loop = True
    elif ratio > cutoff:
        is_loop = False     
    return is_loop

def save_wormviz_image(vizdir,imname,Image,red,green=None,blue=None,
 magneto=None,yellow=None,cyan=None):
    '''
    Save an image with color overlays.  
    vizdir: directory where image will be saved
    imname: name of saved image, include ".jpg"
    Image:  The base grayscale image.
    red,green,blue,magneto,yellow,cyan:  binary arrays that should be
     overlayed. Pass "None" to select a color out of order.  e.g.,
     with one binary array only that should be blue, 
     output = save_wormviz_image(vizdir,imname,Image,None,None,BlueArr) 
    '''
    im = pymorph.overlay(Image,red,green,blue,magneto,yellow,cyan)
    scipy.misc.imsave(os.path.join(vizdir,imname) ,im)
