def file_filter(file_path, filter_str='#'):
    '''
    Function: Filter the line that starts with filter_str('#') in the file_path
    and then return file handle.
    e.g: 
    file(afile) content:
    #1
    #2
    #3
    123
    ----------------------------------
    >>>afile = file_filter(afile_path)
    >>>afile.readline()
    123
    >>>afile.close()
    Attention:  Close the afile
    '''
    infile = open(file_path)
    for i, line in enumerate(infile):
        if not line.startswith(filter_str):
            break
    infile.close()
    infile = open(file_path)
    for j in xrange(i):
        infile.readline()
    return infile
