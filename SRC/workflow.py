'''This module scripts the step-by-step processing of one image.'''

import os

import libcelegans 

import celegansmetadata
import output
import segment
import visualization

import marcdev

def process_frame(file_pathname,max_image,THRESH,THRESH_DIFF,NPOINTS,ANGLE_STEP_SIZE,SCRATCHDIR,DOVIZ,VIZDIR,DETMETHOD,WORMWIDTH,WORMLENGTH): 

    IP = Store()
    IP.scratchdir = SCRATCHDIR
#  get the file meta data
    IP.file_pathname = file_pathname
    IP.FileName,IP.NamePrefix,IP.SeqNumber = celegansmetadata.getjpg_metadata(IP.file_pathname)

#  read the image file
    IP.OrigImage = libcelegans.get_image(file_pathname)

#  perform the image segmentation
    IP.DiffImage,IP.WormBinArr,IP.WormArr = segment.doseg(max_image,IP.OrigImage,THRESH,THRESH_DIFF)


    if DETMETHOD == 0 or DETMETHOD == 1:
    #  get the 1 pixel wide border of the worm
        IP.OrigBorder = libcelegans.get_border(IP.WormBinArr)
    
    #  find the lowest-cost path route around the worm
        IP.BorderRoute = libcelegans.get_border_route(IP.OrigBorder)
    
    #  reduce the number of points in the path
        IP.SimpleBorderRoute = libcelegans.get_simple_border_route(IP.BorderRoute,NPOINTS)
    
    #  get a head and tail assignment (by sharpness)
        IP.Head,IP.Tail,IP.head,IP.tail = libcelegans.get_border_endpoint(
          IP.FileName,IP.OrigBorder,IP.WormBinArr,
          IP.SimpleBorderRoute,step_size=ANGLE_STEP_SIZE)
    
    #  run the intensity test
        IP.PassIntensityTest = libcelegans.perform_head_tail_intensity_test(IP.Head,IP.Tail,IP.WormArr)
    
    #  run the volume test
        IP.PassVolumeTest = libcelegans.perform_head_tail_volume_test(IP.Head,IP.Tail,IP.WormArr)
    
    #  head and tail in hand, use a path path finding algorithm to find the 
    #  lowest-cost path from head to tail, through the border of the worm
        IP.SideOnePath,IP.SideTwoPath,IP.SideOneArr,IP.SideTwoArr = \
         libcelegans.get_side_paths(IP.Head,IP.Tail,IP.BorderRoute)
    
    #  reduce number of points 
        IP.SimpleSideOnePath =  libcelegans.get_simple_side_route(IP.SideOnePath,Nworm_divisions=NPOINTS)
    
    #  reduce number of points
        IP.SimpleSideTwoPath =  libcelegans.get_simple_side_route(IP.SideTwoPath,Nworm_divisions=NPOINTS)
    
    #  compute midline
        IP.MidLine = libcelegans.get_midline(IP.SimpleSideOnePath,IP.SimpleSideTwoPath)
    
    #  test to see if the sides and midline are reasonable
        IP.Is_Loop = libcelegans.perform_loop_test(IP.SideOnePath,IP.SideTwoPath)
        
    #  do the visualization
        if DOVIZ == 'yes':
            visualization.viz(VIZDIR,IP.NamePrefix,IP.DiffImage,IP.OrigImage,IP.SideOnePath,
             IP.SideTwoPath,IP.Tail,IP.Head,IP.MidLine)
        
        
    if DETMETHOD == 0 or DETMETHOD == 2:
    #  begin a section for Marc's development code.  Pass the IP object with all properties that have already been
    #  computed for the image.
        IP.MarcSides,IP.MarcMidLine,IP.MarcOffsets,IP.MarcScore = marcdev.main(IP,max_image,WORMWIDTH,WORMLENGTH)
        
    #  convert MarcMidLine from the TrimImage coordinate system to the OrigImage coordinate system
        IP.OffsetMarcMidLine = [(pt[0]+IP.MarcOffsets['ystart'],pt[1]+IP.MarcOffsets['xstart']) for pt in IP.MarcMidLine]
        IP.OffsetMarcSides = [(pt[0]+IP.MarcOffsets['ystart'],pt[1]+IP.MarcOffsets['xstart']) for pt in IP.MarcSides]

    #  do the visualization
        if DOVIZ == 'yes':
            visualization.mviz(VIZDIR,IP.NamePrefix,IP.OrigImage,IP.OffsetMarcMidLine[-1],IP.OffsetMarcMidLine[0],\
                               IP.OffsetMarcMidLine,IP.OffsetMarcSides)
            
#  write the descriptor file
    DescriptorFile = open(SCRATCHDIR+"/"+IP.NamePrefix+".pickle",'wb')
    if DETMETHOD == 0:
        output.serialize_output(IP.SideOnePath,IP.SideTwoPath,IP.SimpleSideOnePath,IP.SimpleSideTwoPath,
                                 IP.MidLine,IP.head,IP.tail,IP.PassIntensityTest,IP.PassVolumeTest,IP.Is_Loop,
                                 IP.OffsetMarcMidLine,IP.MarcScore,DescriptorFile)
    elif DETMETHOD == 1:
        output.serialize_output(IP.SideOnePath,IP.SideTwoPath,IP.SimpleSideOnePath,IP.SimpleSideTwoPath,
                                 IP.MidLine,IP.head,IP.tail,IP.PassIntensityTest,IP.PassVolumeTest,IP.Is_Loop,
                                 None,None,DescriptorFile)
    elif DETMETHOD == 2:
        output.serialize_output(None,None,None,None,
                                 None,None,None,None,None,None,
                                 IP.OffsetMarcMidLine,IP.MarcScore,DescriptorFile)

class Store:
    pass
