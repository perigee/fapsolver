#! /usr/bin/env python
import os
import sys
#from xml.etree import ElementTree
import xml.etree.ElementTree as ET

def maxcut(fname, outFile):
    tree =  ET.ElementTree(file=fname)
    root = tree.getroot()
    iterCtrs = root.getiterator("constraint")
    iterVars = root.getiterator("variable")
    f = os.open(outFile, 'w')
    firstLine = "p edge %s %s"%(len(iterVars),len(iterCtrs)) 
    f.write(firstLine)
    for e in iterCtrs:
        edge = "e %s %s"%(e.get("name"),e.get("scope"))
        f.write("\n" + edge)
    f.close()

def delegator():
    #maxcut("testBench.xml")
    dirList=os.listdir("benchmark/maxcut-30")
    for fname in dirList: 
        #print fname
        maxcut(fname)

def tryDel():
    for root, dirs, files in os.walk('benchmark/maxcut-30'):
        for f in files:
            maxcut(os.path.join(root,f),os.path.join(root, 'c',f))

if __name__ == '__main__':
    #maxcut("testBench.xml", "testBench.data")
    tryDel()
