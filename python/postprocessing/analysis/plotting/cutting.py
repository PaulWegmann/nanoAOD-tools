from __future__ import division
import ROOT
import numpy as np
from plots import createTH1F, saveratioplot, calc_AMS, saveplot, savegraph, createAxisHists
import os
from array import array
from operator import itemgetter
import pandas as pd

if __name__=="__main__":
    
    infile = ROOT.TFile.Open("../classifiers/data/breg/nores/bckgrd.root")
    bckgrd = infile.Get("Events")
    infile2 = ROOT.TFile.Open("../classifiers/data/breg/nores/signal.root")
    signal = infile2.Get("Events")
    
    cuts = ["deltaR_reconHW<3.3", "deltaR_bb<2", "deltaR_W<2", "Jet_btagDeepB[recon_WJets_id[0]]<0.8", "Jet_btagDeepB[recon_WJets_id[0]]<0.8 && Jet_btagDeepB[recon_WJets_id[1]]<0.85", "recon_Wmass>50&&recon_Wmass<100", "Jet_HT>260", "Jet_pt[recon_bjet_id[0]]>85" ]
    allgoodcuts = "recon_Wmass>50&&recon_Wmass<100 && Jet_btagDeepB[recon_WJets_id[0]]<0.85 && Jet_btagDeepB[recon_WJets_id[1]]<0.9 && deltaR_reconHW<3.5 && recon_Wpt>175"#&& deltaR_W<1.16"
    #allgoodcuts = "recon_Wmass>50&&recon_Wmass<100 && Jet_btagDeepB[recon_WJets_id[0]]<0.8 && deltaR_reconHW<3.3"

    selected_triggers = ["HLT_PFHT300PT30_QuadPFJet_75_60_45_40", "HLT_PFHT180"]
    
    
    weightmultiplier = {"HLT_PFHT300PT30_QuadPFJet_75_60_45_40":"*30000", "HLT_PFHT180":"*30000"}

    #os.chdir("documentation/applyingcuts")
    var = "recon_Hmass"
    #ranges = np.linspace(75, 160, 60)
    ranges = np.linspace(90, 160, 16)
    

    #* ---------------------------------------------------------------------------    
    #! finding best cut ----------------------------------------------------------
    #* ---------------------------------------------------------------------------
    
    #* original code -------------------------------------------------------------
    
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



    #cuttingvars = {"deltaR_reconHW<":[0.5,4.5], "deltaR_W<":[0,5], "Jet_btagDeepB[recon_WJets_id[0]]<":[0.5,1], "Jet_btagDeepB[recon_WJets_id[1]]<":[0.5,1], 
    #"Jet_HT>":[100,600], "recon_Wmass>":[30,100],"recon_Wmass<":[60,200], "recon_Hpt>":[50,250], "recon_Wpt>":[0, 250],"recon_Wpt<":[150, 1000], "Jet_pt[recon_WJets_id[1]]>":[0, 60], 
    #"Jet_pt[recon_WJets_id[0]]>":[0, 60], "Jet_pt[recon_bjet_id[1]]>":[50, 150], "Jet_pt[recon_bjet_id[0]]>":[0, 200] }
    
    #cuttingvars = {"Jet_pt[recon_WJets_id[0]]>":[0, 200]}
    cuttingvars = {"deltaR_bb<":[0, 6], "deltaR_reconHW<":[0.5,4.5], "deltaR_W<": [0,7], "WJet_qgl1>":[0,1], "WJet_qgl2>":[0,1], "WJet_neHEF1<":[0,1],
    "WJet_neHEF2<":[0,1], "WJet_pt1>":[20,400], "WJet_pt2>":[20,400], "WJet_pt2>":[20,400], "bJet_pt2>":[20,400], "bJet_pt1>":[20,400], "recon_Hpt>":[20,500],
     "recon_Wmass>":[30,100],"recon_Wmass<":[60,200], "recon_Wpt>":[0, 250], "recon_Wpt<":[150, 500], "Jet_HT>":[100,1000],
     "bJet_btagDeepB1>":[0.5,1],"bJet_btagDeepB2>":[0.5,1], "WJet_btagDeepB2<":[0,1], "WJet_btagDeepB1<":[0,1], "WJet_btagDeepC1>":[0,1],
     "WJet_btagDeepC2>":[0,1], "bJet_mass1>":[0, 50], "bJet_mass2>":[0, 50], "WJet_mass1>":[0,50], "WJet_mass2>":[0,50],
     "WJet_neEmEF1>":[0,1], "WJet_neEmEF2>":[0,1], "HJet_qgl1>":[0,1], "HJet_qgl2>":[0,1], "WJet_nConstituents1>":[0,40], "WJet_nConstituents2>":[0,40]}
    
    sim = {"deltaR_W<":"#DeltaR_{W} <", "deltaR_reconHW<":"#DeltaR_{HW} <", "deltaR_bb<":"#DeltaR_{H} <", 
    "WJet_qgl2>":"Output of qgl classifier for second W-jet >", "WJet_qgl1>":"Output of qgl classifier for first W-jet >",
    "WJet_neHEF1<":"neutral HEF of first W-jet <", "WJet_neHEF2<":"neutral HEF of second W-jet <", "WJet_pt1>":"p_{t} of first W-jet >",
     "bJet_pt2>":"p_{t} of second H-jet >","bJet_pt1>":"p_{t} of first H-jet >", "WJet_pt2>":"p_{t} of second W-jet >", 
    "recon_Hpt>":"p_{t} of Higgs boson >", "recon_Wmass>":"m_{W} >","recon_Wmass<":"m_{W} <", "recon_Wpt>":"p_{t} of W boson >", 
    "recon_Wpt<":"p_{t} of W boson <", "Jet_HT>":"Summed transverse energy >","bJet_btagDeepB1>":"b-tag output for first H-jet >",
    "bJet_btagDeepB2>":"b-tag output for second H-jet >", "WJet_btagDeepB2<":"b-tag output for second W-jet <", 
    "WJet_btagDeepB1<":"b-tag output for first W-jet <", "WJet_btagDeepC1>":"c-tag output for first W-jet >",
    "WJet_btagDeepC2>":"c-tag output for first W-jet >",
    "bJet_mass1>":"m_{1}^{H} >", "bJet_mass2>":"m_{2}^{H} >", "WJet_mass1>":"m_{1}^{W} >", "WJet_mass2>":"m_{2}^{W} >",
     "WJet_neEmEF1>": "neEmEF_{1}^{W} >", "WJet_neEmEF2>":"neEmEF_{2}^{W} >", "HJet_qgl1>":"qgl_{1}^{H} >", "HJet_qgl2>":"qgl_{2}^{H} >",
     "WJet_nConstituents1>":"nConstituents_{1}^{W} >", "WJet_nConstituents2>":"nConstituents_{2}^{W} >"
     }#"Jet_pt[recon_WJets_id[0]]>":"p_{T} of high energy W jet (GeV) <"}
    #cuttingvars = {"Jet_HT>":[100,600],"recon_Wpt>":[0, 250], "recon_Hpt>":[50,250]}
    
    
    os.chdir("pics/cuts")

    original_AMS = {selected_triggers[0]:0, selected_triggers[1]:0}

    for j in selected_triggers:
        signal1, stats1 = createTH1F(signal, var, ranges, constraints= "({}>0)*weight".format(j), dif="signal")
        bckgrd1, trash = createTH1F(bckgrd, var, ranges, constraints= "({}>0)*weight".format(j))

        original_AMS[j] = calc_AMS(signal1, bckgrd1)
        del(signal1,bckgrd1)
    
    print(original_AMS)

    
    

    dataexport = pd.DataFrame(columns=["variable", "trigger", "cut at", "AMS", "improved"])
    
    for cuttingvar in cuttingvars.keys():
        bestcut = {}

        # signal1, stats1 = createTH1F(signal, cuttingvar[:-1], np.linspace(cuttingvars[cuttingvar][0], cuttingvars[cuttingvar][1], 30), constraints= "({}>0)*weight".format(selected_triggers[0]), dif="signal", normalise=True)
        # bckgrd1, trash = createTH1F(bckgrd, cuttingvar[:-1], np.linspace(cuttingvars[cuttingvar][0], cuttingvars[cuttingvar][1], 30), constraints= "({}>0)*weight".format(selected_triggers[0]), normalise=True)

        # signal2, stats1 = createTH1F(signal, cuttingvar[:-1], np.linspace(cuttingvars[cuttingvar][0], cuttingvars[cuttingvar][1], 30), constraints= "({}>0)*weight".format(selected_triggers[1]), dif="signal", normalise=True)
        # bckgrd2, trash = createTH1F(bckgrd, cuttingvar[:-1], np.linspace(cuttingvars[cuttingvar][0], cuttingvars[cuttingvar][1], 30), constraints= "({}>0)*weight".format(selected_triggers[1]), normalise=True)

        

        # canvas = ROOT.TCanvas("c1","c1")


        # signal1.SetLineColor(ROOT.kBlue)
        # signal1.SetMarkerStyle(20)
        # signal1.SetMarkerColor(ROOT.kBlue+3)
        # signal2.SetMarkerStyle(21)
        # signal2.SetLineColor(ROOT.kGreen)
        # signal2.SetMarkerColor(ROOT.kGreen+3)

        # bckgrd1.SetLineColor(ROOT.kRed)
        # bckgrd1.SetMarkerStyle(22)
        # bckgrd1.SetMarkerColor(ROOT.kRed+3)
        # bckgrd2.SetMarkerStyle(23)
        # bckgrd2.SetLineColor(ROOT.kYellow)
        # bckgrd2.SetMarkerColor(ROOT.kYellow+3)


        # signal1.Draw("SAME")
        # signal2.Draw("SAME")
        # bckgrd1.Draw("SAME")
        # bckgrd2.Draw("SAME")
        # canvas.SaveAs("test.pdf")

        # axish = createAxisHists(2,signal1,signal1.GetXaxis().GetXmin(),signal1.GetXaxis().GetXmax()-0.01)
        # axish[0].GetYaxis().SetRangeUser(0,1.5*max([bckgrd1.GetMaximum(), bckgrd1.GetMaximum()]))
        # axish[0].GetYaxis().SetTitle("Normed entries")



        # saveratioplot(signal1, bckgrd1, sim[cuttingvar][:-1], "signal", "background", ratio = False, name="distributions/"+cuttingvar+"_QUAD", calcArea=True)
        # saveratioplot(signal2, bckgrd2, sim[cuttingvar][:-1], "signal", "background", ratio = False, name="distributions/"+cuttingvar+"_PF180", calcArea=True)

        # del(signal1,bckgrd1, signal2, bckgrd2)


        for j in selected_triggers:
            tmp = []
            for i in np.linspace(cuttingvars[cuttingvar][0],cuttingvars[cuttingvar][1],30):
                    signal1, stats1 = createTH1F(signal, var, ranges, constraints= "("+cuttingvar + str(i) +"&&"+j+ ">0)*weight", dif="signal")
                    bckgrd1, trash = createTH1F(bckgrd, var, ranges, constraints= "("+cuttingvar + str(i) +"&&"+j+ ">0)*weight")
                    s = calc_AMS(signal1, bckgrd1)
                    if not(np.isinf(s)):
                        tmp.append([i,s])
                    del(signal1, bckgrd1)
            bestcut.update({j:np.array(tmp)})
        
        
        os.chdir("cutplots")
        for i in bestcut.keys():
            maximum1 = max(enumerate(bestcut[i][:,1]), key=itemgetter(1))
            if maximum1[1]>original_AMS[i]+0.001:
                dataexport = dataexport.append({"variable":cuttingvar, "trigger":i, "cut at":round(bestcut[i][maximum1[0],0],2),"AMS":(maximum1[1]), "improved":True}, ignore_index=True)
            else:
                dataexport = dataexport.append({"variable":cuttingvar, "trigger":i, "cut at":round(bestcut[i][maximum1[0],0],2),"AMS":(maximum1[1]), "improved":False}, ignore_index=True)
            savegraph(bestcut[i][:,0], bestcut[i][:,1], "cutplot"+cuttingvar[:-1]+"_"+i+"notincsig", sim[cuttingvar], 
            "Binned s/#sqrt{b}: ", text = "Max {} at {}".format(maximum1[1], round(bestcut[i][maximum1[0],0],2)))
        os.chdir("../")


    dataexport.sort_values(by=["variable", "trigger"])
    dataexport.to_csv("bestcuts_test.csv")
        
        #print((bestcut[i]))
        #n = len(bestcut[i][:,0])
        #x, y = array( 'd' ), array( 'd' )
        #for ele in range(n):
            #x.append(bestcut[i][ele,0])
            #y.append(bestcut[i][ele,1])            
        #gr = ROOT.TGraph(n,x,y)
        #saveplot(gr, "cut")
    
    #* copied code for normalised distributions but weighted AMS -----------------
    """
    sim = {"HLT_PFHT300PT30_QuadPFJet_75_60_45_40":"HLT_Quad", "HLT_PFHT180":"HLT_PF"}
    cuttingvars = {"Jet_pt[recon_WJets_id[0]]>":[0, 60], "Jet_pt[recon_bjet_id[0]]>":[0, 200], "deltaR_reconHW<":[0.5,4.5], "deltaR_W<": [0,7]}

    #givencuts = [""]

    var = "deltaR_reconHW<"
    for j in selected_triggers:
        tmp = []
        for i in np.linspace(cuttingvars[var][0],cuttingvars[var][1],100):
                signal1, stats1 = createTH1F(signal, "recon_Hmass", ranges, constraints= "("+var + str(i) +"&&"+j+ ">0)*weight", dif="signal")
                bckgrd1, trash = createTH1F(bckgrd, "recon_Hmass", ranges, constraints= "("+var + str(i) +"&&"+j+ ">0)*weight")
                s = calc_AMS(signal1, bckgrd1)
                tmp.append([i,s])
                del(signal1, bckgrd1)

        tmp = np.array(tmp)
        maximum1 = max(enumerate(tmp[:,1]), key=itemgetter(1))
        f(tmp[:,0], tmp[:,1], "cutplot"+var+"_"+j, "\Delta R of reconstructed H and W", 
        "AMS", text = "Max {} at {}".format(maximum1[1], round(tmp[maximum1[0],0], 2)))
    """
    #* here going for calculating weighted AMS values for given cuts
    
    #* ---------------------------------------------------------------------------
    #! creating plots ------------------------------------------------------------
    #* ---------------------------------------------------------------------------
    """
    #* no cuts -------------------------------------------------------------------
    #* without increasing signal as reference ------------------------------------
    
    signal1, stats1 = createTH1F(signal, var, ranges, constraints= "("+selected_triggers[1]+">0)*weight", dif="signal")
    bckgrd1, trash = createTH1F(bckgrd, var, ranges, constraints= "("+selected_triggers[1]+">0)*weight")
    saveratioplot(signal1,bckgrd1, var, "signal","background", title = "no cut", calcArea = False, calcAMS=True, dif = selected_triggers[1]+"", rdAMS=4)
    del(signal1,bckgrd1)
    root2pandas
    
    signal1, stats1 = createTH1F(signal, var, ranges, constraints= "("+selected_triggers[0]+">0)*weight", dif="signal")
    bckgrd1, trash = createTH1F(bckgrd, var, ranges, constraints= "("+selected_triggers[0]+">0)*weight")
    saveratioplot(signal1,bckgrd1, var, "signal","background", title = "no cut", calcArea = False, calcAMS=True, dif = selected_triggers[0]+"", rdAMS=4)
    del(signal1,bckgrd1)
    
    #* increasing signal ---------------------------------------------------------------------------
    signal1, stats1 = createTH1F(signal, var, ranges, constraints= "("+selected_triggers[1]+">0)*weight"+weightmultiplier[selected_triggers[1]], dif="signal")
    bckgrd1, trash = createTH1F(bckgrd, var, ranges, constraints= "("+selected_triggers[1]+">0)*weight")
    saveratioplot(signal1,bckgrd1, var, "signal","background", title = selected_triggers[1]+"; no cut", calcArea = False, calcAMS=True, dif = selected_triggers[1]+"_s3e3")
    del(signal1,bckgrd1)
    
    
    signal1, stats1 = createTH1F(signal, var, ranges, constraints= "("+selected_triggers[0]+">0)*weight"+weightmultiplier[selected_triggers[0]], dif="signal")
    bckgrd1, trash = createTH1F(bckgrd, var, ranges, constraints= "("+selected_triggers[0]+">0)*weight")
    saveratioplot(signal1,bckgrd1, var, "signal","background", title = selected_triggers[0]+"; no cut", calcArea = False, calcAMS=True, dif = selected_triggers[0]+"_s3e3")
    del(signal1,bckgrd1)
    
    
    #* applying cuts ---------------------------------------------------------------------------
    for j in selected_triggers:
        for i in cuts:
            signal1, stats1 = createTH1F(signal, var, ranges, constraints= "("+i+"&&"+j+ ">0)*weight"+weightmultiplier[j], dif="signal")
            bckgrd1, trash = createTH1F(bckgrd, var, ranges, constraints= "("+i+"&&"+j+">0)*weight")
            saveratioplot(signal1,bckgrd1, var, "signal","background", "signal"+weightmultiplier[j] , dif = i+"_"+j+"_s3e3", calcArea = False, calcAMS=True)
            del(signal1,bckgrd1)
    
    #* ---------------------------------------------------------------------------
    #! combine cuts --------------------------------------------------------------
    #* ---------------------------------------------------------------------------
    """
    #get all them good cuts

    for j in selected_triggers:
        allgoodcuts = ""
        allgoodvars = dataexport.variable.values[dataexport.improved.values == True & (dataexport.trigger.values == j)]
        allgoodvals = dataexport["cut at"].values[dataexport.improved.values == True & (dataexport.trigger.values == j)]
        for s in range(len(allgoodvars)):
            allgoodcuts+= "{}{}&&".format(allgoodvars[s], allgoodvals[s])
        allgoodcuts = allgoodcuts[:-2]


    
        signal1, stats1 = createTH1F(signal, var, ranges, constraints= "("+allgoodcuts+"&&"+j+ ">0)*weight", dif="signal")
        bckgrd1, trash = createTH1F(bckgrd, var, ranges, constraints= "("+allgoodcuts+"&&"+j+">0)*weight")
        AMS = calc_AMS(signal1, bckgrd1)
        del(signal1,bckgrd1)

        signal1, stats1 = createTH1F(signal, var, ranges, constraints= "("+allgoodcuts+"&&"+j+ ">0)*weight", dif="signal")
        bckgrd1, trash = createTH1F(bckgrd, var, ranges, constraints= "("+allgoodcuts+"&&"+j+">0)*weight")
        saveratioplot(signal1,bckgrd1, "m_{H}", "H #rightarrow bb","QCD multijet", dif = "allgoodcuts_{}".format(j), calcArea = False, calcAMS=True, ratio = False, log = True)
        del(signal1,bckgrd1)
    """
    """
    #* ---------------------------------------------------------------------------
    #* ---------------------------------------------------------------------------
    #* ---------------------------------------------------------------------------
        
        
        
