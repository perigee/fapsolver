#! /usr/bin/env python
import os, glob, sys
import json
import fnmatch


def findFile(pattern,path):
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                return name
             


def fapp(cmd,configFile, nb, outputFile):
    path = '../benchmark/roadef2001/all/'    
    jsonFile = open(configFile)
    data = json.load(jsonFile)
    nbA = len(data["A"])
    nbB = len(data["B"])
    nbX = len(data["X"])
    for i in range(nbA):
        level = data["A"][i]["level"]
        filename = data["A"][i]["name"]
        filenamex = filename + "*.in"
        fileData = os.path.join(path,findFile(filenamex,path))
        for n in range(int(level)):
            fileBegin = "echo f:\t" + filename + "\t"  + str(n) + " >> " + outputFile
            os.system(fileBegin)
            for t in range(int(nb)):
                cmdIn = cmd + fileData + " " + str(n) + " >> " + outputFile
                print(str(t) +  ":\t" + fileBegin)
                os.system(cmdIn)
        

def xx(cmd,configFile, nb, outputFile):
    path = '../benchmark/roadef2001/all/'    
    for infile in glob.glob( os.path.join(path,'fapp*.in') ):
        for cnt in range(11):
            fileBegin = "echo f:\t" + infile + "\t"  + str(cnt) + " >> " + outputFile
            os.system(fileBegin)
            for i in range(int(nb)):
                print("working on " + infile + ' '  + str(cnt))
                cmdIn = cmd + infile + " " +str(cnt) + " >> " + outputFile
                os.system(cmdIn)




def delegator():
    args = sys.argv
    cmd = args[1]
    configFile = args[2]
    nb = args[3]
    outputFile = args[4]
    cmd = "./" + cmd + " "
    fapp(cmd,configFile,nb,outputFile)

if __name__ == '__main__':
    delegator()
