from __future__ import division
import ROOT
import PhysicsTools.NanoAODTools.plotting as plot
import os
import numpy as np
from array import array
import sys
from operator import itemgetter

## insert at 1, 0 is the script path (or '' in REPL)
sys.path.append('../plotting')

from plots import saveratioplot, createTH1F, calc_AMS, savegraph

def tmva_create2hist(tree1, tree2, var, cutvar, cut, ranges,norm=False):

    signal1, stats1 = createTH1F(tree1, var, ranges, constraints= "(classID==0&&"+cutvar + str(cut) +")*weight", dif="signal")
    signal2, stats1 = createTH1F(tree2, var, ranges, constraints= "(classID==0&&"+cutvar + str(cut) +")*weight")
    signal1.Add(signal2)
    bckgrd1, trash = createTH1F(tree1, var, ranges, constraints= "(classID==1&&"+cutvar + str(cut) +")*weight", dif="bckgrd")
    bckgrd2, trash = createTH1F(tree2, var, ranges, constraints= "(classID==1&&"+cutvar + str(cut) +")*weight")
    bckgrd1.Add(bckgrd2)
    #print([signal1.GetBinContent(i+1) for i in range(signal1.GetNbinsX())])
    #print([bckgrd1.GetXaxis().GetBinLowEdge(i+1) for i in range(bckgrd1.GetNbinsX())])
    if norm:
        signal1.Scale(1/signal1.Integral())
        bckgrd1.Scale(1/bckgrd1.Integral())
    del(signal2,bckgrd2)
    return(signal1, bckgrd1)

def makecutgraph_fortmva(tree1, tree2, var, cuttingvar,ranges, ncuts, trigger,  dif=""):
    tmp = []
    sim = {"HLT_PFHT300PT30_QuadPFJet_75_60_45_40":"HLT_Quad", "HLT_PFHT180":"HLT_PF"}
    
    maximum = max(tree1.GetMaximum("BDT"), tree2.GetMaximum("BDT"))
    minimum = min(tree2.GetMinimum("BDT"), tree1.GetMinimum("BDT"))
    ranges_BDT =  np.linspace(minimum, maximum, ncuts)

    for i in ranges_BDT:
        signal, bckgrd = tmva_create2hist(tree1, tree2, var, cuttingvar, i, ranges = ranges)
        #print(signal1.Integral()/bckgrd1.Integral())
        s = calc_AMS(signal, bckgrd)
        if not(np.isinf(s)):
            tmp.append([i,s, signal.Integral()])
        
        if i == ranges_BDT[0]:
            #just for comparison---
            # signal.Scale(30000)
            # saveratioplot(signal, bckgrd, var, "signal", "background", title=trigger, calcArea=False, calcAMS=True, dif =dif+"_"+trigger)
            # del(signal,bckgrd)
            print(signal.Integral())
            print(bckgrd.Integral())
            signal, bckgrd = tmva_create2hist(tree1, tree2, "BDT", "classID>","-1000", ranges = ranges_BDT)
            saveratioplot(signal, bckgrd, "BDT output", "H #rightarrow bb", "QCD multijet", calcAMS= False, log = True,  dif = dif+"_"+trigger, ratio = False, name = "BDT_final_{}".format(sim[trigger]))
            
            del(signal,bckgrd)

            # we need more events for final calculation


            signal, _ = createTH1F(tree1, var, ranges, constraints= "(classID==0)*weight")
            bckgrd, _ = createTH1F(tree1, var, ranges, constraints= "(classID==1)*weight")
            finalAMS = calc_AMS(signal, bckgrd)


        del(signal, bckgrd)

    tmp = np.array(tmp)
    index = 0
    value = 0


    for counter, curval in enumerate(tmp[:,1]):
        if tmp[counter,2]>250:
            if value < curval:
                value = curval
                index = counter
    index1, value1 = max(enumerate(tmp[:,1]), key=itemgetter(1))

    sim_trigger = {"HLT_PFHT300PT30_QuadPFJet_75_60_45_40":"HLT_Quad", "HLT_PFHT180":"HLT_PF"}

    signal, bckgrd = tmva_create2hist(tree1, tree2, var, cuttingvar, tmp[index,0], ranges = ranges)
    saveratioplot(signal, bckgrd,  "m_{H} (GeV)",  "H #rightarrow bb", "QCD multijet", log = True, ratio= False, calcAMS = True, name = "recon_Hmass__BDT_{}".format(sim_trigger[trigger]))
    print("{}: {} signal; {} background".format(trigger, signal.Integral(), bckgrd.Integral()))

    savegraph(tmp[:,0], tmp[:,1],"cutplot_"+sim[trigger]+"_final", "BDT output",\
        "Binned s/#sqrt{b}: ", text = "Max {} at {}".format(value1,round(tmp[index1,0],2)))

    return(finalAMS)

    
if __name__=="__main__":
    triggers = ["HLT_PFHT300PT30_QuadPFJet_75_60_45_40", "HLT_PFHT180"]
    ranges = np.linspace(90, 140, 15)

    for j in triggers:

        infile = ROOT.TFile.Open("tmvaBDT"+j+".root")
        
        
        tree = infile.Get("dataset/TestTree")
        tree2 = infile.Get("dataset/TrainTree")

        os.chdir("AMS_onrecon_Hmass")
        makecutgraph_fortmva(tree, tree2, "recon_Hmass", "BDT>", ranges, 30 , j, trigger = j)
        infile.Close()
        os.chdir("../")
        del(infile, tree, tree2)