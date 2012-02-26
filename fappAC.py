#! /usr/bin/env python
import os, glob, sys
import json
import fnmatch


def findFile(pattern,path):
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                return name
             


def fapp(cmd,configFile, outputFile):
    path = '../benchmark/roadef2001/all/'    
    jsonFile = open(configFile)
    data = json.load(jsonFile)
    for instance in ["A","B","X"]:
        nbIns = len(data[instance])
        for i in range(nbIns):
            filename = data[instance][i]["name"]
            filenamex = filename + "*.in"
            fileData = os.path.join(path,findFile(filenamex,path))
            fileBegin = "echo f:\t" + filename + "\t"  +  " >> " + outputFile
            os.system(fileBegin)
            print(str(filename) + "\t")
            cmdIn = cmd + fileData + " 0  "  + " >> " + outputFile
            os.system(cmdIn)
        



def delegator():
    args = sys.argv
    cmd = args[1]
    configFile = args[2]
    outputFile = args[3]
    cmd = "./" + cmd + " "
    fapp(cmd,configFile, outputFile)

if __name__ == '__main__':
    delegator()
