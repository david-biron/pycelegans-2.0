import os
import sys

##########################################################################

POV_HEADER = '''
#include "colors.inc"    // The include files contain
#include "stones.inc"    // pre-defined scene elements
#include "textures.inc"    // pre-defined scene elements
#include "shapes.inc"
#include "glass.inc"
#include "metals.inc"
#include "woods.inc"
'''

POV_CAMERA = '''camera {
    location <%f,%f,%f>
    look_at  <%f,%f,%f>
    angle %f }
'''

POV_AMBIENT_LIGHT = 'global_settings { ambient_light rgb <%f,%f,%f> }'

POV_LIGHT = '''light_source {
    <%f,%f,%f>
    rgb<%f,%f,%f> }
'''

##########################################################################

def main():
  
    allframes = get_midline(sys.argv[1])
    for frame in allframes:
        process_one_frame(frame)


def process_one_frame( (framenumber,PointList) ):
    for i in xrange(0,len(PointList)):
        try:
            x = PointList[i][0] ; y = PointList[i][1] ; nextx = PointList[i+1][0] ; nexty = PointList[i+1][1]
            # Reads the (x, y) points from a point set and prints out to stdout formatted as a cylinder
            # with the current and next x-y pairs acting as the end points
            # This output can be easily modified in a text editor later to the desired pov-ray specifications
            print "cylinder { < %s, %s, 0>, <%s, %s,  0 >, 10, 1 pigment {BODY_COLOR}  normal {bumps 0.5 scale 0.2 } }" % \
             (x,y,nextx,nexty)
        except IndexError:
            # a pov-ray formatted commented line at the end of each frame; simplifies finding frames
            print '//END FRAME ', framenumber
            break
            

#Opens a file from the output data set of pycelegans and extracts the (x, y) coordinates returning as a
#comma separated list of coordinates with a leading frame number
def get_midline(inputfile):
    filelines = open(inputfile,'r').readlines() 
    allframes = []
    for line in filelines:
        words = line.split(',')
        framenumber = words.pop(0)
        framedata = []
        for pt in words:
            try:
                x,y = pt.split() 
                framedata.append( ( float(x),float(y) ) )
            except ValueError:
                break # reached end of line
        allframes.append( (framenumber,framedata) )
    return allframes
 
if __name__ == '__main__':
    main()
