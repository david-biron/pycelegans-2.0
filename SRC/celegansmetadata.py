'''Define routines to parse file names to retreive meta data.'''


def getjpg_metadata(file_pathname):
    FileName = file_pathname.split('/')[-1].strip()
    NamePrefix = FileName.strip('.jpg')
    SeqNumber = int(NamePrefix.split('_')[-1])
    return (FileName,NamePrefix,SeqNumber)
