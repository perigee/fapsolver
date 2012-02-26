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
    #nameFiles = ["02","03","04","05","06","07","08", "09","11"]
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


def dimacEasy(cmd, nb, outputFile):
    #nameFiles = ["myciel3.col 3", "le450_5a.col 4","le450_5b.col 4","le450_5c.col 4","le450_5d.col 4","anna.col 10","david.col 10","huck.col 10", "jean.col 9","DSJC125.1.col 4","DSJR500.1.col 11", "queen5_5.col 4", "queen7_7.col 6","queen10_10.col 9"]
    #nameFiles = ["le450_15b.col 14","le450_15c.col 14","le450_15d.col 14","le450_25a.col 24","le450_25b.col 24","le450_25c.col 24","le450_25d.col 24","myciel4.col 4","myciel5.col 5","myciel6.col 6","myciel7.col 8","queen5_5.col 4","queen6_6.col 6","queen7_7.col 6","queen8_8.col 8","queen9_9.col 9","queen10_10.col 9","queen11_11.col 10","queen12_12.col 11","queen13_13.col 12","queen14_14.col 13","queen15_15.col 14","queen16_16.col 15"]
    #nameFiles = ["DSJC125.1.col 4","DSJR500.1.col 11", "queen5_5.col 4", "queen7_7.col 6","queen10_10.col 9"]
    nameFiles = ["myciel4.col 4","myciel5.col 5","myciel6.col 6","myciel7.col 8","queen5_5.col 4","queen6_6.col 6","queen7_7.col 6","queen8_8.col 8","queen9_9.col 9","queen10_10.col 9","queen11_11.col 10","queen12_12.col 11","queen13_13.col 12","queen14_14.col 13","queen15_15.col 14","queen16_16.col 15"]
    for f in nameFiles:
        path = "benchmark/dimacs/gcp/" + f
        fileBegin = "echo f:\t" + f  + " >> " + outputFile
        cmdIn = cmd + path + " >> " + outputFile
        os.system(fileBegin)
        print fileBegin
        for i in range(int(nb)):
            print i
            os.system(cmdIn)


def maxcut(cmd, nb, outputFile):
    for root, dirs, files in os.walk('benchmark/maxcut-30/c'):
        for f in files:
            path = os.path.join(root,f)
            fileBegin = "echo f:\t" + path  + " >> " + outputFile
            cmdIn = cmd + path  + " >> " + outputFile
            os.system(fileBegin)
            print fileBegin
            for i in range(int(nb)):
                print i
                os.system(cmdIn)



def dimacs(cmd, nb, outputFile):
    #nameFiles = ["myciel3.col 3", "le450_5a.col 4","le450_5b.col 4","le450_5c.col 4","le450_5d.col 4","anna.col 10","david.col 10","huck.col 10", "jean.col 9","DSJC125.1.col 4","DSJR500.1.col 11", "queen5_5.col 4", "queen7_7.col 6","queen10_10.col 9"]
    #nameFiles = ["le450_15b.col 14","le450_15c.col 14","le450_15d.col 14","le450_25a.col 24","le450_25b.col 24","le450_25c.col 24","le450_25d.col 24","myciel4.col 4","myciel5.col 5","myciel6.col 6","myciel7.col 8","queen5_5.col 4","queen6_6.col 6","queen7_7.col 6","queen8_8.col 8","queen9_9.col 9","queen10_10.col 9","queen11_11.col 10","queen12_12.col 11","queen13_13.col 12","queen14_14.col 13","queen15_15.col 14","queen16_16.col 15"]
    #nameFiles = ["DSJC125.1.col 4","DSJR500.1.col 11", "queen5_5.col 4", "queen7_7.col 6","queen10_10.col 9"]
    nameFiles = ["myciel3.col 3", 
                 "myciel4.col 4",
                 "myciel5.col 5",
                 "myciel6.col 6",
                 "myciel7.col 7",
                 "fpsol2.i.1.col 64", 
                 "fpsol2.i.2.col 29",
                 "fpsol2.i.3.col 29",
                 "inithx.i.1.col 53",
                 "inithx.i.2.col 30",
                 "inithx.i.3.col 30",
                 "le450_5a.col 4",
                 "le450_5b.col 4",
                 "le450_5c.col 4",
                 "le450_5d.col 4",
                 "le450_15a.col 14",
                 "le450_15b.col 14",
                 "le450_15c.col 14",
                 "le450_15d.col 14",
                 "le450_25a.col 24",
                 "le450_25b.col 24",
                 "le450_25c.col 24",
                 "le450_25d.col 24",
                 "school1.col 13",
                 "school1_nsh.col 13",
                 "anna.col 10",
                 "david.col 10",
                 "homer.col 12",
                 "huck.col 10",
                 "jean.col 9",
                 "queen5_5.col 4",
                 "queen6_6.col 6",
                 "queen7_7.col 6",
                 "queen8_8.col 8",
                 "queen8_12.col 11",
                 "queen9_9.col 9",
                 "queen10_10.col 9",
                 "queen11_11.col 10",
                 "queen12_12.col 11",
                 "queen13_13.col 12",
                 "queen14_14.col 13",
                 "queen15_15.col 14",
                 "queen16_16.col 15"]                 
    for f in nameFiles:
        path = "../benchmark/" + f
        fileBegin = "echo f:\t" + f  + " >> " + outputFile
        cmdIn = cmd + path + " >> " + outputFile
        os.system(fileBegin)
        print fileBegin
        for i in range(int(nb)):
            print i
            os.system(cmdIn)



def delegator():
    args = sys.argv
    cmd = args[1]
    nb = args[2]
    outputFile = args[4]
    #print nb
    #print outputFile
    cmd = "./" + cmd + " "
    if(args[3] == "graph"):
        graph(cmd, nb,outputFile)
    if(args[3] == "celar"):
        celar(cmd, nb,outputFile)
    if(args[3] == "rlfap"):
        celar(cmd, nb,outputFile)
        graph(cmd, nb,outputFile)
    if(args[3] == "maxcut"):
        maxcut(cmd,nb,outputFile)
    if(args[3] == "gcp"):
        dimacEasy(cmd,nb,outputFile)
    if(args[3] == "color"):
        dimacs(cmd,nb,outputFile)

if __name__ == '__main__':
    delegator()
