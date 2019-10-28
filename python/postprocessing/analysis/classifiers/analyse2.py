import ROOT
import sys
import os
import array
import numpy as np
import pandas as pd
sys.path.append('../plotting')
#sys.path.append('../')

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

    setvarsdict = {"HLT_PFHT300PT30_QuadPFJet_75_60_45_40":{"MinNodeSize":"1%", "NTrees":400, "MaxDepth":3}, "HLT_PFHT180": {"MinNodeSize":"1%", "NTrees":1000, "MaxDepth":3}}
    
    varoptions = {"NTrees":np.array([50, 100, 250, 400, 600, 1000], dtype = str), "MinNodeSize":["0.5", "1%", "2.5%", "5%", "7.5%", "10%"],\
        "BoostType":["AdaBoost", "RealAdaBoost", "Bagging", "Grad"],\
        "SeparationType":["CrossEntropy", "GiniIndex", "GiniIndexWithLaplace", "MisClassificationError", "SDivSqrtSplusB"],\
        "VarTransform":["G","N", "D", "P", "U"], "MaxDepth":[3,4, 5, 7, 10, 15], "NodePurityLimit":[0.1, 0.3, 0.5, 0.7, 0.9], "UseFisherCuts":["True", "False"]
        }

    #varoptions = { "BoostType":AdaBoost, "MaxDepth":[3], "NTrees":[500], "MinNodeSize":["2.5%"] }

    standard = {"HLT_PFHT300PT30_QuadPFJet_75_60_45_40":"BoostType=AdaBoost:VarTransform=U:SeparationType=CrossEntropy", "HLT_PFHT180":"BoostType=AdaBoost:VarTransform=U:SeparationType=GiniIndexWithLaplace"}

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

    optionsdict = {"HLT_PFHT300PT30_QuadPFJet_75_60_45_40": ["SeparationType=CrossEntropy:MaxDepth=3:MinNodeSize=1%:NTrees=100",
        "SeparationType=CrossEntropy:MaxDepth=3:MinNodeSize=1%:NTrees=400",
        "SeparationType=CrossEntropy:MaxDepth=3:MinNodeSize=10%:NTrees=250",
        "SeparationType=CrossEntropy:MinNodeSize=1%:NTrees=100:MaxDepth=5",
        "SeparationType=CrossEntropy:MaxDepth=3:MinNodeSize=5%:NTrees=400"],
        "HLT_PFHT180":["SeparationType=GiniIndexWithLaplace:MaxDepth=3:MinNodeSize=1%:NTrees=1000",
        "SeparationType=GiniIndexWithLaplace:MinNodeSize=1%:NTrees=100:MaxDepth=5",
        "SeparationType=GiniIndexWithLaplace:MinNodeSize=1%:NTrees=400:MaxDepth=5",
        "SeparationType=GiniIndexWithLaplace:MinNodeSize=1%:NTrees=600:MaxDepth=5",
        "SeparationType=GiniIndexWithLaplace:NTrees=1000:MinNodeSize=7.5%:MaxDepth=15"
        ] }

    #optionsdict = {"HLT_PFHT300PT30_QuadPFJet_75_60_45_40": ["VarTransform=U:BoostType=AdaBoost:SeparationType=CrossEntropy:MinNodeSize=1%:NTrees=400:MaxDepth=3"], \
    #   "HLT_PFHT180": ["VarTransform=U:BoostType=AdaBoost:SeparationType=MisClassificationError:MinNodeSize=1%:NTrees=400:MaxDepth=3"]}

    for trigger in triggers:
        for i in range(3):
            for var in varoptions.keys():
                for option in varoptions[var]:
                            options = "{}={}".format(var,option)
                

                            """
            for curopt in optionsdict[trigger]:
                            options = "{}:{}".format(standard[trigger], curopt)
                            """
                            """
            for index1, curvar1 in enumerate(setvars.keys()):
                for index2, curvar2 in enumerate(setvars.keys()):
                    if index2<=index1: continue

                    for index3, curopt1 in enumerate(varoptions[curvar1]):
                        for index4, curopt2 in enumerate(varoptions[curvar2]):
                        #if index4<index2: continue
                            #print(curopt1, curopt2)
                            #print(curvar1, curvar2)
                            #print(index1, index2)
                            tmp = ["MinNodeSize", "NTrees", "MaxDepth"]
                            del(tmp[index1], tmp[index2-1]) # is done recursively

                            options = "{}:{}={}:{}={}:{}={}".format(standard[trigger], tmp[0], setvarsdict[trigger][tmp[0]], curvar1, curopt1, curvar2, curopt2)
                            """    

                            optionssimplified = ''.join(e for e in options if e.isalnum())
                            
                            #    setvars = {"MinNodeSize":"1%", "NTrees":500, "MaxDepth":10}

                            train(trigger, options = options)
                            print("!!!!!!!!!!!!!!!!!!!!!!!")
                            evaluate(trigger)
                            
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

                        
    dataexport.to_csv("run2/varyingall.csv")
