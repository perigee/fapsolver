#! /usr/bin/env python
import os, glob, sys
import json
import fnmatch


# read config file
def findFile(pattern,path):
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                return name
             



def fapp(cmd,configFile, nb, outputFile):
    path = '../benchmark/roadef2001/all/'    
    jsonFile = open(configFile)
    data = json.load(jsonFile)
    for instances in ["A","B","X"]:
        nbIns = len(data[instances])
        for i in range(nbIns):
            level = data[instances][i]["level"]
            filename = data[instances][i]["name"]
            filenamex = filename + "*.in"
            fileData = os.path.join(path,findFile(filenamex,path))
            for n in range(int(level)):
                fileBegin = "echo f:\t" + filename + "\t"  + str(n) + " >> " + outputFile
                os.system(fileBegin)
                for t in range(int(nb)):
                    cmdIn = cmd + fileData + " " + str(n) + " >> " + outputFile
                    print(str(t) +  ":\t" + fileBegin)
                    os.system(cmdIn)
        

# execute the instances without level
def nolevelfapp(cmd,configFile, nb, outputFile):
    nb = 1
    path = '../benchmark/roadef2001/all/'    
    jsonFile = open(configFile)
    data = json.load(jsonFile)
    for instances in ["A","B","X"]:
        nbIns = len(data[instances])
        for i in range(nbIns):
            level = data[instances][i]["level"]
            filename = data[instances][i]["name"]
            filenamex = filename + "*.in"
            fileData = os.path.join(path,findFile(filenamex,path))
            n = 0
            fileBegin = "echo f:\t" + filename + "\t"  + str(n) + " >> " + outputFile
            os.system(fileBegin)
            for t in range(int(nb)):
                cmdIn = cmd + fileData + " " + str(n) + " >> " + outputFile
                print(str(t) +  ":\t" + fileBegin)
                os.system(cmdIn)
        


def delegator():
    args = sys.argv
    cmd = args[1]
    configFile = args[2]
    nb = args[3]
    outputFile = args[4]
    cmd = "./" + cmd + " "
    nolevelfapp(cmd,configFile,nb, outputFile)
    #fapp(cmd,configFile,nb, outputFile)

if __name__ == '__main__':
    delegator()
