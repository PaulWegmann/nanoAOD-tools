import ROOT
import sys
import os
import array
import numpy as np
import pandas as pd
sys.path.append('../plotting')

from plots import saveratioplot, createTH1F, calc_AMS, savegraph, saveplot
from tmvaBDT import train, evaluate, analyse

if __name__ == "__main__":
    dataexport = pd.DataFrame(columns=[ "trigger", "options", "AMS", "MAD", "MSE", "ROCIntegral"])

    ntrees = np.linspace(100, 2000, 20)
    ntrees = [int(x) for x in ntrees]
    #ntrees = [250]
    triggers = ["HLT_PFHT300PT30_QuadPFJet_75_60_45_40", "HLT_PFHT180"]
    for trigger in triggers:
        for ntree in ntrees:
            options = "NTrees={}:MinNodeSize=7.5%:BoostType=AdaBoost:SeparationType=MisClassificationError:VarTransform=G:nCuts=400:MaxDepth=5".format(ntree)
            train(trigger, options)
            evaluate(trigger)
            mad_val, mse_val, AMSval = analyse(trigger)

            infile = ROOT.TFile.Open("tmvaBDT{}.root".format(trigger), 'read')
            
            path = "dataset/Method_BDT/BDT/"
            roc = infile.Get(path+"MVA_BDT_rejBvsS")
            roc_value = roc.Integral()
            infile.Close()

            dataexport = dataexport.append({"trigger":trigger, "options":options, "AMS":AMSval, "MAD":mad_val, "MSE":mse_val, "ROCIntegral":roc_value}, ignore_index=True)                    

    dataexport.to_csv("analyse.csv")
