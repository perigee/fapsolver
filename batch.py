#! /usr/bin/env python
import os
import sys

def graph(nb, outputFile):
    nameFiles = ["03","04","05","06","07","10","11","12","13"]
    for f in nameFiles:
        path = "celar/FullRLFAP/GRAPH/graph" + f
        fileBegin = "echo f:\t" + path  + " >> " + outputFile
        dom = path + "/dom.txt "
        dom = dom + path + "/varm.txt "
        dom = dom + path + "/ctrm.txt "
        cmd = "./testCELAR " + dom  + " >> " + outputFile
        os.system(fileBegin)
        print fileBegin
        for i in range(int(nb)):
            print i
            os.system(cmd)


def celar(nb, outputFile):
    #nameFiles = ["02","03","04","05","06","07","08", "09","11"]
    nameFiles = ["01", "02","03","04","05","06","07","08", "09","11"]
    for f in nameFiles:
        path = "celar/CELAR/scen" + f
        fileBegin = "echo f:\t" + path  + " >> " + outputFile
        dom = path + "/DOMt.TXT "
        dom = dom + path + "/VARt.TXT "
        dom = dom + path + "/CTRt.TXT s "
        cmd = "./testCELAR " + dom  + " >> " + outputFile
        os.system(fileBegin)
        print fileBegin
        for i in range(int(nb)):
            print i
            os.system(cmd)



def delegator():
    args = sys.argv
    nb = args[1]
    outputFile = args[3]
    #print nb
    #print outputFile
    if(args[2] == "graph"):
        graph(nb,outputFile)
    if(args[2] == "celar"):
        celar(nb,outputFile)
    if(args[2] == "rlfap"):
        celar(nb,outputFile)
        graph(nb,outputFile)


if __name__ == '__main__':
    delegator()
