'''Dump memory structure to a binary file.'''

import cPickle

def serialize_output(SideOnePath,SideTwoPath,SimpleSideOnePath,SimpleSideTwoPath,MidLine,
                     head,tail,PassIntensityTest,PassVolumeTest,Is_Loop,
                     OffsetMarcMidLine,MarcScore,DescriptorFile):

    cPickle.dump(({"s1path":SideOnePath,"s2path":SideTwoPath,
                  "simples1path":SimpleSideOnePath,"simples2path":SimpleSideTwoPath,
                  "midline":MidLine, "head":head, "tail":tail,"intensitytest":PassIntensityTest,
                  "volumetest":PassVolumeTest,"is_loop":Is_Loop,
                  'marcmidline':OffsetMarcMidLine,'marcscore':MarcScore}),
                   DescriptorFile,cPickle.HIGHEST_PROTOCOL)