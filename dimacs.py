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
             



def dimacs(cmd,configFile, nb, outputFile):
    path = '../benchmark/dimacs/'    
    jsonFile = open(configFile)
    data = json.load(jsonFile)
    for instances in ["Easy"]:
        nbIns = len(data[instances])
        for i in range(nbIns):
            level = data[instances][i]["colors"]
            filename = data[instances][i]["name"]
            filenamex = filename + "*.col"
            fileData = os.path.join(path,findFile(filenamex,path))
            fileBegin = "echo f:\t" + filename + "\t"  + str(level) + " >> " + outputFile
            os.system(fileBegin)
            for t in range(int(nb)):
                cmdIn = cmd + fileData + " " + str(level) + " >> " + outputFile
                print(str(t) +  ":\t" + fileBegin)
                os.system(cmdIn)
        


def delegator():
    args = sys.argv
    cmd = args[1]
    configFile = args[2]
    nb = args[3]
    outputFile = args[4]
    cmd = "./" + cmd + " "
    dimacs(cmd,configFile,nb, outputFile)


if __name__ == '__main__':
    delegator()
