#! /usr/bin/env python
import os
import sys

def graph(cmd, nb, outputFile):
    nameFiles = ["03","04","05","06","07","10","11","12","13"]
    for f in nameFiles:
        path = "celar_bk/FullRLFAP/GRAPH/graph" + f
        fileBegin = "echo f:\t" + path  + " >> " + outputFile
        dom = path + "/dom.txt "
        dom = dom + path + "/varm.txt "
        dom = dom + path + "/ctrm.txt "
        cmdIn = cmd + dom  + " >> " + outputFile
        os.system(fileBegin)
        print fileBegin
        for i in range(int(nb)):
            print i
            os.system(cmdIn)


def celar(cmd, nb, outputFile):
    nameFiles = ["01", "02","03","04","05","06","07","08", "09","11"]
    for f in nameFiles:
        path = "celar_bk/CELAR/scen" + f
        fileBegin = "echo f:\t" + path  + " >> " + outputFile
        dom = path + "/DOMt.TXT "
        dom = dom + path + "/VARt.TXT "
        dom = dom + path + "/CTRt.TXT s "
        cmdIn = cmd + dom  + " >> " + outputFile
        os.system(fileBegin)
        print fileBegin
        for i in range(int(nb)):
            print i
            os.system(cmdIn)



def delegator():
    args = sys.argv
    cmd = args[1]
    nb = args[2]
    outputFile = args[3]
    #print nb
    #print outputFile
    cmd = "./" + cmd + " "
    celar(cmd, nb,outputFile)
    graph(cmd, nb,outputFile)
    


if __name__ == '__main__':
    delegator()
