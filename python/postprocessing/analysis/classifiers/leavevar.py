import ROOT
import sys
import os
import array
import numpy as np
import pandas as pd
sys.path.append('../plotting')
import copy

from plots import saveratioplot, createTH1F, calc_AMS, savegraph, saveplot
from tmvaBDT import train, evaluate, analyse


def flatten(A):
    rt = []
    for i in A:
        if isinstance(i,list): rt.extend(flatten(i))
        else: rt.append(i)
    return rt

#from plotsbdt import makecutgraph_fortmva, tmva_create2hist

if __name__=="__main__":

    
    # *  defining variales -----------------------------------
    infile='atlas-higgs-challenge-2014-v2_part.root'
    triggers = ["HLT_PFHT300PT30_QuadPFJet_75_60_45_40", "HLT_PFHT180"]
    #triggers = ["HLT_PFHT180"]

    var = ["NTrees", "MinNodeSize", "BoostType", "SeparationType", "VarTransform"]
    setvars = {"MinNodeSize":"1%", "NTrees":500, "MaxDepth":10}
    
    varoptions = {"NTrees":np.array([100, 250, 400, 600, 1000], dtype = str), "MinNodeSize":["1%", "2.5%", "5%", "7.5%", "10%"],\
        "BoostType":["AdaBoost", "RealAdaBoost", "Bagging", "Grad"],\
        "SeparationType":["CrossEntropy", "GiniIndex", "GiniIndexWithLaplace", "MisClassificationError", "SDivSqrtSplusB"],\
        "VarTransform":["G","N", "D", "P", "U"], "MaxDepth":[3, 5, 7, 10, 15]
        }

    #varoptions = { "MaxDepth":[3], "NTrees":[500], "MinNodeSize":["2.5%"] }
    #varoptions = {"VarTransform":["P", "U", "G", "N", "D"], "BoostType":["AdaBoost", "RealAdaBoost", "Bagging", "Grad"],\
    #    "SeparationType":["CrossEntropy", "GiniIndex", "GiniIndexWithLaplace", "MisClassificationError", "SDivSqrtSplusB"]}
    dataexport = pd.DataFrame(columns=[ "trigger", "leftvar", "AMS", "AMSrecon", "MAD", "MSE", "ROCIntegral", "overtrain"])
    #data = pd.read_csv("analyse/analyse.csv")
    
    # * train -------------------------------------------------
    index = 0

    # list used for determining and saving the best AMS value for default variables, varying just one parameter
    best_AMS = {}
    for curvar in var:
        best_AMS[curvar] = []

    #? blue stuff used for different analysis

    #* very good options filtered: (consindering ROC Integral>0.84)
    #* HLT_PFHT180
    #* VarTransform=U:BoostType=AdaBoost:SeparationType=MisClassificationError:NTrees=500:MinNodeSize=1%:MaxDepth=3
    #* VarTransform=U:BoostType=AdaBoost:SeparationType=MisClassificationError:MinNodeSize=1%:NTrees=400:MaxDepth=3
    #* VarTransform=U:BoostType=AdaBoost:SeparationType=MisClassificationError:MaxDepth=10:MinNodeSize=5%:NTrees=600

    #* HLT_QUAD
    #* VarTransform=U:BoostType=AdaBoost:SeparationType=CrossEntropy:NTrees=500:MinNodeSize=2.5%:MaxDepth=3
    #* VarTransform=U:BoostType=AdaBoost:SeparationType=CrossEntropy:NTrees=500:MinNodeSize=2.5%:MaxDepth=5
    #* VarTransform=U:BoostType=AdaBoost:SeparationType=CrossEntropy:MinNodeSize=1%:NTrees=1000:MaxDepth=3
    #* VarTransform=U:BoostType=AdaBoost:SeparationType=CrossEntropy:MinNodeSize=1%:NTrees=600:MaxDepth=3
    #* VarTransform=U:BoostType=AdaBoost:SeparationType=CrossEntropy:MinNodeSize=1%:NTrees=400:MaxDepth=3

    optionsdict = {"HLT_PFHT300PT30_QuadPFJet_75_60_45_40": ["VarTransform=U:BoostType=AdaBoost:SeparationType=CrossEntropy:NTrees=500:MinNodeSize=2.5%:MaxDepth=3", \
        "VarTransform=U:BoostType=AdaBoost:SeparationType=CrossEntropy:NTrees=500:MinNodeSize=2.5%:MaxDepth=5",
        "VarTransform=U:BoostType=AdaBoost:SeparationType=CrossEntropy:MinNodeSize=1%:NTrees=1000:MaxDepth=3",
        "VarTransform=U:BoostType=AdaBoost:SeparationType=CrossEntropy:MinNodeSize=1%:NTrees=600:MaxDepth=3",
        "VarTransform=U:BoostType=AdaBoost:SeparationType=CrossEntropy:MinNodeSize=1%:NTrees=400:MaxDepth=3"],
        "HLT_PFHT180":["VarTransform=U:BoostType=AdaBoost:SeparationType=MisClassificationError:NTrees=500:MinNodeSize=1%:MaxDepth=3",
        "VarTransform=U:BoostType=AdaBoost:SeparationType=MisClassificationError:MinNodeSize=1%:NTrees=400:MaxDepth=3",
        "VarTransform=U:BoostType=AdaBoost:SeparationType=MisClassificationError:MaxDepth=10:MinNodeSize=5%:NTrees=600"] }

    optionsdict = {"HLT_PFHT300PT30_QuadPFJet_75_60_45_40": ["VarTransform=U:BoostType=AdaBoost:SeparationType=CrossEntropy:MinNodeSize=1%:NTrees=400:MaxDepth=3"], \
        "HLT_PFHT180": ["VarTransform=U:BoostType=AdaBoost:SeparationType=MisClassificationError:MinNodeSize=1%:NTrees=400:MaxDepth=3"]}

    variables = ["recon_Wmass", "recon_Hmass", "deltaR_reconHW","recon_Hpt","recon_Wpt","deltaR_bb", "deltaR_W", 
            "Jet_HT",
            ["bJet_pt1","bJet_pt2"],
            ["WJet_pt1","WJet_pt2"],
            ["bJet_btagDeepB1","bJet_btagDeepB2"],
            ["WJet_btagDeepB1","WJet_btagDeepB2"],
            ["WJet_btagDeepC1", "WJet_btagDeepC2"],
            #["HJet_btagDeepC1", "HJet_btagDeepC2"],
            ["bJet_mass1", "bJet_mass2"],
            ["WJet_mass1", "WJet_mass2"],
            #["bJet_"]
            ["WJet_neEmEF1", "WJet_neEmEF2"],
            ["WJet_neHEF1", "WJet_neHEF2"],
            ["WJet_qgl1", "WJet_qgl2"],
            ["HJet_qgl1", "HJet_qgl2"],
            ["WJet_nConstituents1", "WJet_nConstituents2"]
            ] #curcount = ~22


    ###? leftvars for Quad: NTrees=200:MaxDepth=3:MinNodeSize=5%, SplitSeed=123
    leavevar = ["WJet_neHEF2", "WJet_mass1", "WJet_btagDeepC2", "WJet_mass2", "HJet_qgl1", "WJet_neEmEF1", 
    "bJet_btagDeepB1", "bJet_mass1", "WJet_neEmEF2", "WJet_btagDeepC1", "HJet_qgl2", 
    "WJet_nConstituents2", "WJet_nConstituents1", "WJet_neHEF1", "WJet_pt2", 
    "bJet_btagDeepB2", "recon_Hpt", "bJet_pt1", "WJet_btagDeepB1", "bJet_pt2", 
    "Jet_HT", "WJet_pt1" ]

    ###? leftvars for PFHT180: NTrees=200:MaxDepth=3:MinNodeSize=5%, SplitSeed=123
    leavevar = ["WJet_neEmEF1", "WJet_neHEF2", "WJet_pt1"]

    final_quad = [leavevar]

    variables = flatten(variables)
    

    for trigger in [triggers[0]]:
        for _ in range(1):
            for j in range(len(leavevar)):
                    #if variables[j] == "recon_Hmass":
                    #    continue
                    variables.remove(leavevar[j])

                    #s = copy.copy(variables)
                    #del(s[j])
                    #s = flatten(s)
                    #print(s)
                    
                    #    setvars = {"MinNodeSize":"1%", "NTrees":500, "MaxDepth":10}
                    train(trigger,options = "NTrees=200:MaxDepth=3:MinNodeSize=5%" , variables = variables)
                    print("!!!!!!!!!!!!!!!!!!!!!!!")
                    evaluate(trigger, variables = variables)
                    
                    mad_val, mse_val, AMSss, AMSrecon, signalevents, bckgrdevents, minimum, maximum = analyse(trigger, options=j, return_range = True) # * average mean distance as further parameter to check how distinctive the bdt is

                    
                    # ---------------------------------------------------------------------
                    # * calculate parameter which should describe training and result -----
                    # TODO ROC integral : for test and train? -----------------------------
                    # ? There is a train but no test ROC? ---------------------------------
                    # TODO ranking of parameters?????? ------------------------------------
                    # TODO AMS value (but that happend before) ----------------------------
                    # ---------------------------------------------------------------------
                    # * reading in outfile for analysis

                    infile = ROOT.TFile.Open("tmvaBDT{}.root".format(trigger), 'read')

                    #print(infile.cd('dataset/Method_BDT/BDT'))

                    #hist = ROOT.TH1F("tmp", "tmp", 100, 0, 1)
                    
                    path = "dataset/Method_BDT/BDT/"

                    # hist = infile.Get(path+"BDT_tr_S")
                    roc = infile.Get(path+"MVA_BDT_rejBvsS")

                    #* ---------------------------------------------
                    #* calcualte value for overtraining check ------
                    #* ---------------------------------------------

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

                    roc_value = roc.Integral()
                    print("ROC Integral: {}".format(roc_value))
                    # saveplot(hist, "BDT")
                    
                    #data["ROCIntegral"][index] = roc_value

                    dataexport = dataexport.append({"trigger":trigger, "leftvar":leavevar[j], "AMS":AMSss,"AMSrecon":AMSrecon, "sigevents":signalevents, "bckevents":bckgrdevents, "MAD":mad_val, "MSE":mse_val, "ROCIntegral":roc_value, "overtrain":overtrainmetric}, ignore_index=True)                    
                    index+=1
                    infile.Close()

                    
                    #del(hist)

                        
    dataexport.to_csv("leavevar/leftvar_200_3_fixed.csv")
  
