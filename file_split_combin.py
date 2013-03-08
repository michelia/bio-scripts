#encoding=utf8
import sys
from os.path import exists
fileCount = 0

'python实现大文件分割与合并'

def splitFile(fileName, partSize=1):
    # 1024 * 1024 = 1048576
    length = partSize * 1048576
    f1 = open(fileName, "rb")
    while True:
        content = f1.read(length)
        if content == "":
            break
        newFile = distFile(fileName)
        f2 = open(newFile, "wb")
        f2.write(content)
        f2.close()
    f1.close()
    print 'split file complete!'
   

def distFile(sourceFile):
    global fileCount
    fileCount = fileCount + 1
    extPos = sourceFile.rfind('.')
    if extPos > 0:
        return sourceFile + '.part' + str(fileCount)
    else:    # extPos == -1
        print 'File type? Can not split!'
        sys.exit(1)

def combine(filename):
    count = 0
    extPos = filename.find('.part')
    if extPos > 0:
        file = filename[:extPos]
        f1 = open(file, "wb")
        while True:
            count = count + 1
            partFile = file + '.part' + str(count)
            if not exists(partFile):
                break
            else:
                f2 = open(partFile, "rb")
                content = f2.read()
                f2.close()
                f1.write(content)
        f1.close()
        print 'combine file complete!'

    else:
        print 'File type? Can not combine!'

def usage():
    print ''''usage is "file.py s[c] filename" or "file.py s[c] filename size(M) \n 's' — 'split', 'c' — 'combine"'''

if __name__ == "__main__":
    if len(sys.argv) !=3 and len(sys.argv) !=4:
        usage()
        sys.exit(1)
    if sys.argv[1] == 's':
        if len(sys.argv) == 3:
            splitFile(sys.argv[2])
        elif len(sys.argv) == 4 and int(sys.argv[3]) > 0:
            splitFile(sys.argv[2], int(sys.argv[3]))
        else:
            usage()
        sys.exit(1)
   
    elif sys.argv[1] =='c':
        if len(sys.argv) == 3:
            combine(sys.argv[2])
        else:
            usage()
    else:
        usage()

  # 使用方法：
  # 分割：file.py s 文件路径 每块的文件大小（不填默认为1M）

  # 合并：file.py c 任意部分文件路径