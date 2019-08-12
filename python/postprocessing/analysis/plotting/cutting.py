from __future__ import division
import ROOT
import numpy as np
from plots import createTH1F, saveratioplot, calc_AMS, saveplot, savegraph
import os
from array import array
from operator import itemgetter
import pandas as pd

if __name__=="__main__":
    
    infile = ROOT.TFile.Open("../background/bckgrd.root")
    bckgrd = infile.Get("Events")
    infile2 = ROOT.TFile.Open("../signal/signal.root")
    signal = infile2.Get("Events")
    
    cuts = ["deltaR_reconHW<3.3", "deltaR_bb<2", "deltaR_W<2", "Jet_btagDeepB[recon_WJets_id[0]]<0.8", "Jet_btagDeepB[recon_WJets_id[0]]<0.8 && Jet_btagDeepB[recon_WJets_id[1]]<0.85", "recon_Wmass>50&&recon_Wmass<100", "Jet_HT>260", "Jet_pt[recon_bjet_id[0]]>85" ]
    allgoodcuts = "recon_Wmass>50&&recon_Wmass<100 && Jet_btagDeepB[recon_WJets_id[0]]<0.85 && Jet_btagDeepB[recon_WJets_id[1]]<0.9 && deltaR_reconHW<3.5"
    #allgoodcuts = "recon_Wmass>50&&recon_Wmass<100 && Jet_btagDeepB[recon_WJets_id[0]]<0.8 && deltaR_reconHW<3.3"

    selected_triggers = ["HLT_PFHT300PT30_QuadPFJet_75_60_45_40", "HLT_PFHT180"]
    
    
    weightmultiplier = {"HLT_PFHT300PT30_QuadPFJet_75_60_45_40":"*3*10e3", "HLT_PFHT180":"*3*10e3"}
    
    os.chdir("documentation/applyingcuts")
    var = "recon_Hmass"
    ranges = np.linspace(75, 160, 60)
    

    # ---------------------------------------------------------------------------    
    # finding best cut ----------------------------------------------------------
    # ---------------------------------------------------------------------------
    """
    cuttingvars = {"deltaR_reconHW<":[0,5], "deltaR_W<":[0,5], "Jet_btagDeepB[recon_WJets_id[0]]<":[0.5,1], "Jet_btagDeepB[recon_WJets_id[1]]<":[0.5,1], "Jet_HT>":[100,400], "recon_Wmass>":[30,100],"recon_Wmass<":[60,200], "recon_Hpt>":[50,150], "recon_Wpt>":[0, 150], "Jet_pt[recon_WJets_id[1]]>":[0, 60], "Jet_pt[recon_WJets_id[0]]>":[0, 60], "Jet_pt[recon_bjet_id[1]]>":[50, 150], "Jet_pt[recon_bjet_id[0]]>":[50, 150] }
    
    dataexport = pd.DataFrame(columns=["variable", "trigger", "cut at", "AMS"])
    
    for cuttingvar in cuttingvars.keys():
        bestcut = {}
        for j in selected_triggers:
            tmp = []
            for i in np.linspace(cuttingvars[cuttingvar][0],cuttingvars[cuttingvar][1],100):
                    signal1, stats1 = createTH1F(signal, var, ranges, constraints= "("+cuttingvar + str(i) +"&&"+j+ ">0)*weight"+weightmultiplier[j], dif="signal")
                    bckgrd1, trash = createTH1F(bckgrd, var, ranges, constraints= "("+cuttingvar + str(i) +"&&"+j+ ">0)*weight")
                    s = calc_AMS(signal1, bckgrd1)
                    if not(np.isinf(s)):
                        tmp.append([i,s])
                    del(signal1, bckgrd1)
            bestcut.update({j:np.array(tmp)})
        
        
        for i in bestcut.keys():
            maximum1 = max(enumerate(bestcut[i][:,1]), key=itemgetter(1))
            #savegraph(bestcut[i][:,0], bestcut[i][:,1], "cutplot"+cuttingvar[:-1]+"_"+i, cuttingvar, "AMS", text = "Max {} at {}".format(int(maximum1[1]), round(bestcut[i][maximum1[0],0],2)), heading = i)
            dataexport = dataexport.append({"variable":cuttingvar, "trigger":i, "cut at":round(bestcut[i][maximum1[0],0],2),"AMS":int(maximum1[1])}, ignore_index=True)
    
    dataexport.sort_values(by=["variable", "trigger"])
    dataexport.to_csv("bestcuts.csv")
        
        #print((bestcut[i]))
        #n = len(bestcut[i][:,0])
        #x, y = array( 'd' ), array( 'd' )
        #for ele in range(n):
            #x.append(bestcut[i][ele,0])
            #y.append(bestcut[i][ele,1])            
        #gr = ROOT.TGraph(n,x,y)
        #saveplot(gr, "cut")
    """
    # ---------------------------------------------------------------------------
    #creating plots -------------------------------------------------------------
    # ---------------------------------------------------------------------------
    """
    signal1, stats1 = createTH1F(signal, var, ranges, constraints= "("+selected_triggers[1]+">0)*weight"+weightmultiplier[selected_triggers[1]], dif="signal")
    bckgrd1, trash = createTH1F(bckgrd, var, ranges, constraints= "("+selected_triggers[1]+">0)*weight")
    saveratioplot(signal1,bckgrd1, var, "signal","background", title = "no cut", calcArea = False, calcAMS=True, dif = selected_triggers[1]+"_s3e3")
    del(signal1,bckgrd1)
    
    
    signal1, stats1 = createTH1F(signal, var, ranges, constraints= "("+selected_triggers[0]+">0)*weight"+weightmultiplier[selected_triggers[0]], dif="signal")
    bckgrd1, trash = createTH1F(bckgrd, var, ranges, constraints= "("+selected_triggers[0]+">0)*weight")
    saveratioplot(signal1,bckgrd1, var, "signal","background", title = "no cut", calcArea = False, calcAMS=True, dif = selected_triggers[0]+"_s3e3")
    del(signal1,bckgrd1)
    
    
     
    for j in selected_triggers:
        for i in cuts:
            signal1, stats1 = createTH1F(signal, var, ranges, constraints= "("+i+"&&"+j+ ">0)*weight"+weightmultiplier[j], dif="signal")
            bckgrd1, trash = createTH1F(bckgrd, var, ranges, constraints= "("+i+"&&"+j+">0)*weight")
            saveratioplot(signal1,bckgrd1, var, "signal","background", "signal"+weightmultiplier[j] , dif = i+"_"+j+"_s3e3", calcArea = False, calcAMS=True)
            del(signal1,bckgrd1)
        
    # ---------------------------------------------------------------------------
    # combine triggers ----------------------------------------------------------
    # ---------------------------------------------------------------------------
    """
    for j in selected_triggers:
        signal1, stats1 = createTH1F(signal, var, ranges, constraints= "("+allgoodcuts+"&&"+j+ ">0)*weight", dif="signal")
        bckgrd1, trash = createTH1F(bckgrd, var, ranges, constraints= "("+allgoodcuts+"&&"+j+">0)*weight")
        saveratioplot(signal1,bckgrd1, var, "signal","background", title = "allgoodcuts; signal"+weightmultiplier[j] , dif = "allgoodcuts_"+j+"", calcArea = False, calcAMS=True, rdAMS=4)
        del(signal1,bckgrd1)
    
    """
    """
    # ---------------------------------------------------------------------------
    # ---------------------------------------------------------------------------
    # ---------------------------------------------------------------------------
        
        
        
