from __future__ import division
import ROOT
import PhysicsTools.NanoAODTools.plotting as plot
import os
import numpy as np
from array import array
import sys
from operator import itemgetter
from tmvaBDT import train
import pandas as pd
from plotsbdt import makecutgraph_fortmva

sys.path.append('../plotting')

from plots import saveratioplot, createTH1F, calc_AMS, savegraph


if __name__=="__main__":
    triggers = ["HLT_PFHT300PT30_QuadPFJet_75_60_45_40", "HLT_PFHT180"]
    var = ["NTrees", "MinNodeSize", "BoostType", "SeparationType"]
    #var = ["BoostType"]
    varoptions = {"NTrees":np.array([100, 250, 400, 600, 1000], dtype = str), "MinNodeSize":["1%", "2.5%", "5%", "7.5%", "10%"],
    "BoostType":["AdaBoost", "RealAdaBoost", "Bagging", "Grad"], 
    "SeparationType":["CrossEntropy", "GiniIndex", "GiniIndexWithLaplace", "MisClassificationError", "SDivSqrtSplusB"]}



    sdf = "NTrees=400:MinNodeSize=3.5%:BoostType=RealAdaBoost:UseBaggedBoost:BaggedSampleFraction=0.6:nCuts=250:MaxDepth=20:VarTransform=G:SeparationType=GiniIndexWithLaplace"

    dataexport = pd.DataFrame(columns=[ "trigger", "options", "AMS"])


    for trigger in triggers:
        for curvar in var:
            for curopt in varoptions[curvar]:
                train(trigger, "!V:{}={}:nCuts=250:VarTransform=G".format(curvar, curopt))
                
                infile = ROOT.TFile.Open("tmvaBDT"+trigger+".root")
        
                tree = infile.Get("dataset/TestTree")
                print(tree.GetEntries())
                tree2 = infile.Get("dataset/TrainTree")

                os.chdir("varyingvar")

                value = makecutgraph_fortmva(tree, tree2, "recon_Hmass", "BDT>", np.linspace(75, 160, 60),np.linspace(0,1,100), trigger, dif = ''.join(e for e in curvar+curopt if e.isalnum()))

                dataexport = dataexport.append({"trigger":trigger, "options":curvar+"="+curopt, "AMS":value}, ignore_index=True)                    

                infile.Close()
                os.chdir("../")
                del(infile, tree, tree2)

    dataexport.to_csv("tracktraining.csv")
