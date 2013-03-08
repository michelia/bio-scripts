try:
    import cPickle as pickle
except:
    import pickle

def dump(data, pathToFile):
    '''
    Save a data to a pathToFile according the pickle moudle.
    '''
    toFile = open(pathToFile, 'w')
    pickle.dump(data, toFile)
    toFile.close()

def load(pathToFile):
    '''
    Load a data from a pathToFile. The dataFile is pickle file(.piso)
    '''
    toFile = open(pathToFile)
    data = pickle.load(toFile)
    toFile.close()
    return data
