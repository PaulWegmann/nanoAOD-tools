import subprocess
import shlex
import ROOT
import sys
import os
import array
import numpy as np
import pandas as pd
sys.path.append('../plotting')
import copy

from plots import saveratioplot, createTH1F, calc_AMS, savegraph, saveplot
from tmvaBDT2 import train, evaluate, analyse

def overtrain_ROC(trigger):
    infile = ROOT.TFile.Open("tmvaBDT{}.root".format(trigger), 'read')

    test = infile.Get("dataset/TestTree")
    traintree = infile.Get("dataset/TrainTree")

    overtrainmetric=0
    for compare in ["classID<1", "classID>0"]:
        testhist = ROOT.TH1F("train", "train", 30, minimum, maximum)                        
        trainhist = ROOT.TH1F("test", "test", 30, minimum, maximum)

        test.Draw("BDT>>test", "({})*weight".format(compare))
        traintree.Draw("BDT>>train", "({})*weight".format(compare))

        testhist.Scale(1/testhist.Integral())
        trainhist.Scale(1/trainhist.Integral())

        final = 0
        for i in range(testhist.GetNbinsX()):
            final+= min(testhist.GetBinContent(i+1), trainhist.GetBinContent(i+1)) 
        overtrainmetric+=final #* metric for checking overtraining, checking for similiar shapes of test and train -> "best value" is close to 2
        del(testhist, trainhist)
    
    roc = infile.Get(path+"MVA_BDT_rejBvsS")
    roc_value = roc.Integral()
    del(infile)
    # print("ROC Integral: {}".format(roc_value))

    return(overtrainmetric, roc_value)
    

triggers = ["HLT_PFHT300PT30_QuadPFJet_75_60_45_40", "HLT_PFHT180"]

path = "dataset/Method_BDT/BDT/"

dataexport = pd.DataFrame(columns=[ "trigger", "leftvar", "AMS", "AMSrecon", "MAD", "MSE", "ROCIntegral", "overtrain"])


for trigger in triggers:
    
    allvar= ["recon_Wmass", 
            "recon_Hmass", 
            "deltaR_reconHW",
            "recon_Hpt",
            "recon_Wpt",
            "deltaR_bb", 
            "deltaR_W", 
            "Jet_HT",
            "bJet_pt1",
            "bJet_pt2",
            "WJet_pt1",
            "WJet_pt2",
            "bJet_btagDeepB1",
            "bJet_btagDeepB2",
            "WJet_btagDeepB1",
            "WJet_btagDeepB2",
            "WJet_btagDeepC1", 
            "WJet_btagDeepC2",
            "bJet_mass1",
            "bJet_mass2",
            "WJet_mass1", 
            "WJet_mass2",
            "WJet_neEmEF1",
            "WJet_neEmEF2",
            "WJet_neHEF1", 
            "WJet_neHEF2",
            "WJet_qgl1", 
            "WJet_qgl2",
            "HJet_qgl1",
            "HJet_qgl2",
            "WJet_nConstituents1", 
            "WJet_nConstituents2"

            #"HJet_btagDeepC1", "HJet_btagDeepC2",
            #["bJet_"]
            ] #curcount = ~22


    leavevars = []
    
    process = subprocess.Popen(['python', 'tmvaBDT2.py',  'train', 'none', '{}'.format(trigger)], stdout=subprocess.PIPE)
    stdout = process.communicate()[0]
    #print('STDOUT:{}'.format(stdout))

    found = False

    
    for index, line in enumerate(stdout.splitlines()):
        if 'method specific' in line:
            found = True
            start = index+5
        if found:
            if "Factory" in line:
                found = False
                end = index-2

    for var in allvar:
        if var in stdout.splitlines()[end]:
            lastvar = var
            
    print("lastvar : {}".format(stdout.splitlines()[end]))#, stdout.splitlines()[start])

    evaluate(trigger, variables = allvar)
    leavevars.append(lastvar)    

    mad_val, mse_val, AMSss, AMSrecon, signalevents, bckgrdevents, minimum, maximum = analyse(trigger, return_range = True) # * average mean distance as further parameter to check how distinctive the bdt is
    overtrainmetric, roc_value = overtrain_ROC(trigger)

    dataexport = dataexport.append({"trigger":trigger, "leftvar":"", "AMS":AMSss,"AMSrecon":AMSrecon, "sigevents":signalevents, "bckevents":bckgrdevents, "MAD":mad_val, "MSE":mse_val, "ROCIntegral":roc_value, "overtrain":overtrainmetric}, ignore_index=True)                    

    while (roc_value>82.0):

        process = subprocess.Popen(['python', 'tmvaBDT2.py',  'train', '{}'.format(leavevars), '{}'.format(trigger)], stdout=subprocess.PIPE)
        stdout = process.communicate()[0]
        # print('STDOUT:{}'.format(stdout))

        found = False
        for index, line in enumerate(stdout.splitlines()):
            if 'method specific' in line:
                found = True
                start = index+5
            if found:
                if "Factory" in line:
                    found = False
                    end = index-2

        for var in allvar:
            if var in stdout.splitlines()[end]:
                lastvar = var
                break
                
        print("lastvar : {}".format(stdout.splitlines()[end]))#, stdout.splitlines()[start])

        allvar.remove(leavevars[-1])

        evaluate(trigger, variables = allvar)
        leavevars.append(lastvar)    

        mad_val, mse_val, AMSss, AMSrecon, signalevents, bckgrdevents, minimum, maximum = analyse(trigger, return_range = True) # * average mean distance as further parameter to check how distinctive the bdt is
        overtrainmetric, roc_value = overtrain_ROC(trigger)
        print("ROC: {}".format(roc_value))

        dataexport = dataexport.append({"trigger":trigger, "leftvar":leavevars[-2], "AMS":AMSss,"AMSrecon":AMSrecon, "sigevents":signalevents, "bckevents":bckgrdevents, "MAD":mad_val, "MSE":mse_val, "ROCIntegral":roc_value, "overtrain":overtrainmetric}, ignore_index=True)                    


dataexport.to_csv("leavevar/n200_max3_min5_seed838.csv")

# process = subprocess.Popen(['echo', '"Hello stdout"'], stdout=subprocess.PIPE)
# stdout = process.communicate()[0]
# print 'STDOUT:{}'.format(stdout)



