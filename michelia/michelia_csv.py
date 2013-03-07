import csv
def csvreader(fileHandel):
    return csv.reader(fileHandel, 'excel-tab')
def csvwriter(fileHandel):
    return csv.writer(fileHandel, 'excel-tab')


class CsvReader():
    """docstring for CsvReader
    preLinesNum: pre read lines num"""
    def __init__(self, filePath ,preLinesNum=0, format='excel-tab'):
        self.preLineNum=preLinesNum 
        self.fileHandel = open(filePath)
        self.format = format
        for i in xrange(self.preLineNum):
            self.fileHandel.readline()
    def reader(self):
        for line in csv.reader(self.fileHandel, self.format):
            yield line
    def close(self):
        self.fileHandel.close()

class CsvWriter():
    """docstring for CsvWriter"""
    def __init__(self, filePath, format='excel-tab'):
        self.filePath = filePath
        self.fileHandel = open(filePath)
        self.writer = csv.writer(fileHandel, format)
    def writerow(self, aList):
        self.writer.writerow(aList)
    def close(self):
        self.filePath.close()