try:
    import cPickle as pickle
except:
    import pickle

def data_dump(pathToFile, data):
    '''
    Save a data to a pathToFile according the pickle moudle.
    '''
    print pathToFile
    toFile = open(pathToFile, 'w')
    pickle.dump(data, toFile)
    toFile.close()

def data_load(pathToFile):
    '''
    Load a data from a pathToFile. The dataFile is pickle file(.piso)
    '''
    toFile = open(pathToFile)
    data = pickle.load(toFile)
    toFile.close()
    return data
