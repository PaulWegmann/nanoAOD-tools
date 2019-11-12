import ROOT
import sys
import os
import array
import numpy as np
import pandas as pd
sys.path.append('../plotting')

from plots import saveratioplot, createTH1F, calc_AMS, savegraph, saveplot
from tmvaBDT import train, evaluate, analyse
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
    dataexport = pd.DataFrame(columns=[ "trigger", "options", "AMS", "MAD", "MSE", "ROCIntegral", "overtrain"])
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

    for trigger in triggers:
        #! training with fixed options
        for i in range(10): 
                    #? no observed difference for:
                    #? quad = [WJet_btagDeepC2, WJet_mass1, WJet_mass2, WJet_neEmEF1, WJet_neEmEF2, WJet_neHEF2, WJet_qgl2]
                    #? ht180 = [WJet_pt1, bJet_mass1, WJet_mass2, WJet_neEmEF1, WJet_neEmEF2, WJet_neHEF1, WJet_neHEF2]
                    options = "NTrees=200:MaxDepth=3:MinNodeSize=5%"
                    if trigger == triggers[0]:
                        variables = variables= ["recon_Wmass", 
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
                            "WJet_mass1", #*curresults
                            "WJet_mass2", #*curresults
                            "WJet_neEmEF1", #*curresults
                            "WJet_neEmEF2", #*curresults
                            "WJet_neHEF1", #*curresults
                            "WJet_neHEF2", #*curresults
                            "WJet_qgl1", 
                            "WJet_qgl2",
                            "HJet_qgl1", #*curresults
                            "HJet_qgl2", #*curresults
                            "WJet_nConstituents1", 
                            "WJet_nConstituents2"

                            #"HJet_btagDeepC1", "HJet_btagDeepC2",
                            #["bJet_"]
                            ]
                    if trigger == triggers[1]:
                        variables = variables= ["recon_Wmass", 
                            "recon_Hmass", 
                            "deltaR_reconHW",
                            "recon_Hpt",
                            "recon_Wpt",
                            "deltaR_bb", 
                            "deltaR_W", 
                            "Jet_HT",
                            "bJet_pt1",
                            "bJet_pt2",
                            "WJet_pt1", #*curresults
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
                            "WJet_mass2", #*curresults
                            "WJet_neEmEF1", #*curresults
                            "WJet_neEmEF2", #*curresults
                            "WJet_neHEF1", #*curresults
                            "WJet_neHEF2", #*curresults
                            "WJet_qgl1", 
                            "WJet_qgl2",
                            "HJet_qgl1", #*curresults
                            "HJet_qgl2", #*curresults
                            "WJet_nConstituents1", 
                            "WJet_nConstituents2"

                            #"HJet_btagDeepC1", "HJet_btagDeepC2",
                            #["bJet_"]
                            ]
        #! --------------------------------------

        #! training for finally selected options
        # for i in range(10):
        #     for option in optionsdict[trigger]:
        #             options = option
        #! --------------------------------------

        #for index1, curvar1 in enumerate(setvars.keys()): #? change to setvars.keys()
            #for index2, curvar2 in enumerate(tmp.keys()): #? change to setvars
        #    tmp = {"MinNodeSize":"1%", "NTrees":500, "MaxDepth":10}
        #    del tmp[curvar1]
        #    for curopt1 in varoptions[tmp.keys()[0]]:
        #        for curopt2 in varoptions[tmp.keys()[1]]:
                    # options = "VarTransform=U:BoostType=AdaBoost:SeparationType=CrossEntropy:{}={}:{}={}:{}={}".format(curvar1, setvars[curvar1], tmp.keys()[0], curopt1, tmp.keys()[1], curopt2)

                    ##? used for different analysis
                    #?for i in range(3):
                    #?for curopt in varoptions[curvar1]:
                    #?    del setvars[curvar1]
                    #?    options = "SeparationType=GiniIndexWithLaplace:VarTransform=G:{}={}:{}={}:{}={}".format(setvars.keys()[0],\
                    #?        setvars[setvars.keys()[0]], setvars.keys()[1], setvars[setvars.keys()[1]], curvar1, curopt)
                    #?options = "NTrees=500:{}={}".format(curvar1, curopt)
                    optionssimplified = ''.join(e for e in options if e.isalnum())
                    
                    #    setvars = {"MinNodeSize":"1%", "NTrees":500, "MaxDepth":10}

                    train(trigger, options = options, variables = variables)
                    print("!!!!!!!!!!!!!!!!!!!!!!!")
                    evaluate(trigger, variables = variables)
                    

                    mad_val, mse_val, AMSss, AMSrecon, signalevents, bckgrdevents, minimum, maximum = analyse(trigger, options=optionssimplified, return_range = True) # * average mean distance as further parameter to check how distinctive the bdt is

                    
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

                    dataexport = dataexport.append({"trigger":trigger, "options":options, "AMS":AMSss,"AMSrecon":AMSrecon, "sigevents":signalevents, "bckevents":bckgrdevents, "MAD":mad_val, "MSE":mse_val, "ROCIntegral":roc_value, "overtrain":overtrainmetric}, ignore_index=True)                    
                    index+=1
                    infile.Close()


                #del(hist)

                        
    dataexport.to_csv("leavevar/results/var_ntrain10_varU_n500_max2_min5_bobe02.csv")
  
