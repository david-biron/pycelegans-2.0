#!/usr/bin/env python
########################################################################
#  Author: Nicholas Labello, labello@uchicago.edu
#  The University of Chicago Research Computing Center
########################################################################
'''
==============
collectoutput
==============
This executable compiles output from the image processing runs. 
'''

import cPickle
import os
import string
import sys

import libutil

def main():

    pickledir=sys.argv[1]
    FileNames = libutil.get_file_names(pickledir,filter='.pickle')
    FileNames.sort()

    outdir = './properties'
    libutil.makedir(outdir)

    s1path_out = open(os.path.join(outdir,"sideonepath.txt"),'w')
    s2path_out = open(os.path.join(outdir,"sidetwopath.txt"),'w')
    simples1path_out = open(os.path.join(outdir,"simplesideonepath.txt"),'w')
    simples2path_out = open(os.path.join(outdir,"simplesidetwopath.txt"),'w')
    midline_out = open(os.path.join(outdir,"midline.txt"),'w')
    nosetip_out = open(os.path.join(outdir,"nosetip.txt"),'w')
    tailtip_out = open(os.path.join(outdir,"tailtip.txt"),'w')
    intensity_test = open(os.path.join(outdir,"intensitystest.txt"),'w')
    volume_test = open(os.path.join(outdir,"volumetest.txt"),'w')
    int_vol_test = open(os.path.join(outdir,"int-vol-test.txt"),'w')   
    is_loop_test = open(os.path.join(outdir,"is_loop.txt"),'w')
    marcmidline_out = open(os.path.join(outdir,"marcmidline.txt"),'w')
    marcscore_out = open(os.path.join(outdir,"marcscore.txt"),'w')

    for file in FileNames:
        NamePrefix = file.strip('.pickle')
        SeqNumber = int(NamePrefix.split('_')[-1])
        try:
            l = cPickle.load(open(file,'rb')) 
            nosetiptxt = "%6i, %4i, %4i \n" % (SeqNumber,l["head"][0],l["head"][1]) 
            tailtiptxt = "%6i, %4i, %4i \n" % (SeqNumber,l["tail"][0],l["tail"][1]) 
            intensitytxt = "%6i, %s \n" %(SeqNumber,l["intensitytest"])
            volumetxt = "%6i, %s \n" %(SeqNumber,l["volumetest"])
            int_vol_txt = "%6i, %s %s \n" %(SeqNumber,l["intensitytest"],l["volumetest"])
            is_loop_txt = "%6i, %s \n" % (SeqNumber,l["is_loop"])
            midlinetxt = format_list_of_points(l["midline"],SeqNumber)
            s1pathtxt = format_list_of_points2(l["s1path"],SeqNumber)
            s2pathtxt = format_list_of_points2(l["s2path"],SeqNumber)
        except:
            nosetiptxt   = "%6i, -1, -1 \n" % SeqNumber
            tailtiptxt   = "%6i, -1, -1 \n" % SeqNumber
            intensitytxt = "%6i, -1 \n" % SeqNumber
            volumetxt = "%6i, -1 \n" % SeqNumber
            int_vol_txt = "%6i -1, -1\n" % SeqNumber
            is_loop_txt = "%6i -1 \n" %SeqNumber
            midlinetxt   = "%6i, -1 \n" % SeqNumber
            s1pathtxt    = "%6i, -1 \n" % SeqNumber 
            s2pathtxt    = "%6i, -1 \n" % SeqNumber
            
        try:
            l = cPickle.load(open(file,'rb')) 
            marcmidlinetxt = format_list_of_points2(l["marcmidline"],SeqNumber)
            marcscoretxt = "%6i, %8.2f \n" % (SeqNumber,l["marcscore"])
        except:
            marcmidlinetxt   = "%6i, -1 \n" % SeqNumber
            marcscoretxt = "%6i, -1 \n" % SeqNumber

        nosetip_out.write(nosetiptxt)
        tailtip_out.write(tailtiptxt)
        intensity_test.write(intensitytxt)
        volume_test.write(volumetxt)
        int_vol_test.write(int_vol_txt)
        is_loop_test.write(is_loop_txt)
        midline_out.write(midlinetxt)
        s1path_out.write(s1pathtxt)
        s2path_out.write(s2pathtxt)
        marcmidline_out.write(marcmidlinetxt)
        marcscore_out.write(marcscoretxt)

def format_list_of_points(L,N):
    FL = ("%6i " % N) 
    for element in L:
        FL =  FL + string.join( (str(element.x),' ',str(element.y),','))
    FL = FL + '\n'
    return FL 

def format_list_of_points2(L,N):
    FL = ("%6i " % N) 
    for element in L:
        FL = FL + string.join( (str(element[0]),' ',str(element[1]),','))
    FL = FL + '\n'
    return FL  
   
if __name__ == '__main__':
    main()


