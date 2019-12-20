import ROOT
import os
import numpy as np
import sys
from operator import itemgetter
import pandas as pd
import PhysicsTools.NanoAODTools.plotting as plot



from plots import saveplot, createTH1F, saveratioplot, calc_AMS, savegraph, createAxisHists

sys.path.append('../classifiers')

from plotsbdt import makecutgraph_fortmva, tmva_create2hist


if __name__=="__main__":
    # ------------------------------------------------------------------------------------
    # reading in data and definitions

    # var = [0"Jet_pt", 1"Jet_phi", 2"Jet_eta", 3"deltaR", 4"recon_Hmass", 5"nJets", 6"MET_pt", 7"fatjetsused", 8"original_Hmass", 9"genjets_Hmass", 10"nrightjets", 11"Jet_btagDeepB", 12"nmatch", 13"candidate_score_1", 14"wnrightjets", 15"genjets_deltaR_W", 16"recon_WMass", 17"genjets_WMass", 18"wnrightjets_bad"] "original_Hmass" , "nrightjets", "candidate_score_", "wnrightjets", "genjets_deltaR_W", "genjets_WMass", "wnrightjets_bad", "genjets_Hmass", "fatjetsused"
    var = ["Jet_pt",  "Jet_eta", "deltaR_bb", "recon_Hmass", "nJets", "MET_pt" ,  "Jet_btagDeepB", "nmatch",    
    "recon_Wmass",  "Jet_HT", "deltaR_W", "recon_Hpt", "recon_Wpt", "recon_Weta", "deltaR_reconHW", "MET_pt", "nSoftActivityJet"]

    #var = ["Jet_btagDeepC", "Jet_mass", "Jet_nElectrons", "Jet_hadronFlavour"]
    #    ranges = np.linspace(90, 140, 15)


    ranges = {"Jet_pt":np.linspace(0,400, 30), "Jet_phi":np.linspace(-4,4,20), "Jet_eta":np.linspace(-4,4,20), 
    "deltaR_bb":np.linspace(0,6,30), "recon_Hmass":np.linspace(90, 160, 16), "nJets":np.linspace(3.5,12.5,10), 
    "MET_pt":np.linspace(0,250,70), "fatjetsused":np.linspace(0,2,3), "original_Hmass":np.linspace(0,200,50), 
    "genjets_Hmass:":np.linspace(50, 200, 50), "nrightjets":np.linspace(-0.5,2.5,4), "Jet_btagDeepB":np.linspace(0.5,1,30), 
    "nmatch":np.linspace(1.5,6.5,6), "candidate_score_1":np.linspace(1.6,8,20), "wnrightjets":np.linspace(-0.5,2.5,4), 
    "genjets_deltaR_W":np.linspace(0,10,20), "bckgrdrecon_Wmass":np.linspace(0,150, 15), "genjets_WMass":np.linspace(0,150, 40), 
    "wnrightjets_bad":np.linspace(-0.5,2.5,4), "Jet_HT":np.linspace(100, 1000, 30), "weight":np.linspace(0, 10, 50), 
    "deltaPT_W":np.linspace(0, 300, 50), "deltaR_W":np.linspace(0,7,30), "recon_Wpt":np.linspace(0,500,30), 
    "recon_Hpt":np.linspace(0,400,30), "recon_Heta":np.linspace(-3,3,20), "recon_Weta":np.linspace(.3,3), 
    "deltaR_reconHW":np.linspace(0.5,4.5,30), "Jet_btagDeepC":np.linspace(0,1,30), "Jet_mass":np.linspace(0,50,30), 
    "Jet_nElectrons":np.linspace(-0.5,5.5,7), "Jet_hadronFlavour":np.linspace(-20,20,42), "nSoftActivityJet":np.linspace(-0.5,10.5,12),
    "recon_Wmass":np.linspace(20, 140, 15), "Jet_btagDeepB":np.linspace(0,1,20)}

    triggers = ["HLT_AK8PFJet80",  "HLT_PFJet80", "HLT_PFHT300PT30_QuadPFJet_75_60_45_40", "HLT_DiPFJetAve35_HFJEC",  "HLT_PFHT180", "HLT_DoublePFJets40_CaloBTagCSV_p33",  "HLT_AK4PFJet80", "HLT_AK4CaloJet80"]
    selected_triggers = ["HLT_PFHT300PT30_QuadPFJet_75_60_45_40", "HLT_PFHT180" ]

    standardconstraintsH = "deltaR_bjet1 < .3 && deltaR_bjet2 < .3" #first selection is to used to ensure that i can connect the genpart level with the genjet level properly; 
    testing_constraints = "deltaR_bb<1. && deltaR_bb>.0"
    
    standardconstraintsW = "deltaR_W_jet1< .5 && deltaR_W_jet2 < .5"

    infile = ROOT.TFile.Open("../classifiers/data/bckgrd.root")
    bckgrd_nobreg = infile.Get("Events")
    
    infile2 = ROOT.TFile.Open("../signal/signal.root")
    signal_nobreg = infile2.Get("Events")

    os.chdir("pics")

    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! BDT example theory part ------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    
    # path = "dataset/Method_BDT/BDT/"


    # infile = ROOT.TFile.Open("../../classifiers/tmvaBDT{}.root".format(selected_triggers[1]), 'read')
    # tree = infile.Get("dataset/TestTree")

    # roc = infile.Get(path+"MVA_BDT_rejBvsS")

    # canvas = ROOT.TCanvas("c1","c1")
    # roc.SetXTitle("Background Rejection")
    # roc.SetYTitle("Signal Efficency")
    # roc.SetLineColor(ROOT.kRed)
    # roc.Draw("L")
    # canvas.SaveAs("ROC_example.pdf")
    # canvas.Close()

    # signal1, bckgrd1 = tmva_create2hist(tree, tree, "BDT", "BDT>", 0, np.linspace(0,1, 30))

    # saveratioplot(signal1, bckgrd1, "BDT output", "signal", "background", ratio = False, name = "BDT_ROC_example")
    
    

    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! Chi-squared distro -----------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    def chisquare2(x, par):
        if 5.05 >x[0]:
            return(0)
        elif round(x[0],1) == 8:
            return(0)
        else:
            return(1/np.sqrt(np.pi*2*x[0])*np.exp(-x[0]/2))

    chisquare = ROOT.TF1("chisquare", "1/sqrt(TMath::Pi())*sqrt(x)*sqrt(2)*exp(-x/2)", 0, 8)
    canvas = ROOT.TCanvas("c1","c1")
    chisquare_test = ROOT.TF1("chisquare", chisquare2, 5, 8, 10000)
    


    #axish = createAxisHists(2,chisquare,chisquare.GetXaxis().GetXmin(),chisquare.GetXaxis().GetXmax()-0.01)
    #axish[0].GetYaxis().SetRangeUser(0,1.5*max([chisquare.GetMaximum(), chisquare.GetMaximum()]))
    #axish[0].GetYaxis().SetTitle()


    #chisquare.Draw()
    chisquare_test.Draw("F")
    canvas.SaveAs("chi-square.pdf")
    canvas.Close()
    
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! Comparison of with breg and without  -----------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------

    infile3 = ROOT.TFile.Open("../../classifiers/data/breg/nores/bckgrd.root")
    # infile3 = ROOT.TFile.Open("../../background/breg2/bckgrd.root")

    bckgrd = infile3.Get("Events")
    
    infile4 = ROOT.TFile.Open("../../classifiers/data/breg/nores/signal.root")
    # infile4 = ROOT.TFile.Open("../../signal/breg2/signal.root")
    signal = infile4.Get("Events")

    var = "recon_Hmass"
    signal1, _ = createTH1F(signal_nobreg, var, ranges[var], "weight", dif = "signal")
    bckgrd1, _  = createTH1F(signal, var, ranges[var], "weight")


    function1 = ROOT.TF1("m1","gaus",min(ranges["recon_Hmass"]), max(ranges["recon_Hmass"]))
    function2 = ROOT.TF1("m2","gaus",min(ranges["recon_Hmass"]), max(ranges["recon_Hmass"]))
    function1.SetLineColor(ROOT.kBlue)
    function1.SetParameters( np.array([1000, 120, 60]))
    function2.SetParameters( np.array([1000, 120, 60]))

    signal1.Fit("m1")
    bckgrd1.Fit("m2")

    # signal1.Fit("gaus")
    # bckgrd1.Fit("gaus")


    saveratioplot(bckgrd1, signal1, "m_{H} (GeV)", "With b-regression", "Without b-regression", name = "Compare_breg", ratio = False)

    #function1 = ROOT.TF1("m1","gaus",min(ranges["recon_Hmass"]), max(ranges["recon_Hmass"]))
    #function2 = ROOT.TF1("m2","gaus",min(ranges["recon_Hmass"]), max(ranges["recon_Hmass"]))
    #function1.SetParameters( np.array([1000, 120, 60]))
    #function2.SetParameters( np.array([1000, 120, 60]))

    #function1 = signal1.Fit("m1")
    #signal1.SetParameters()
    #function2 = bckgrd1.Fit("m2")
    



    # print("ITS HAPPENING --------------------------")
    # print("without breg:")
    # print("mitte " + function1.GetParameter(1))
    # print("width "+ function1.GetParameter(2))
    # print("with breg:")
    # print("mitte " + function2.GetParameter(1))
    # print("width "+ function2.GetParameter(2))


    del(signal1, bckgrd1)
    signal1, _ = createTH1F(signal, var, ranges[var], "weight", dif = "signal")
    bckgrd1, _  = createTH1F(bckgrd, var, ranges[var], "weight")

    saveratioplot(signal1, bckgrd1, "m_{H} (GeV)", "H #rightarrow bb", "QCD multijet", name = "breg_AMS", ratio=False, log=True, calcAMS=True)
    del(signal1, bckgrd1)

    # signal1, _ = createTH1F(signal_nobreg, var, ranges[var], "weight", dif = "signal")
    # bckgrd1, _  = createTH1F(bckgrd_nobreg, var, ranges[var], "weight")

    # saveratioplot(signal1, bckgrd1, "m_{H} (GeV)", "signal", "bckgrd", name = "nobreg_AMS", ratio=False, log=True, calcAMS=True)
    # del(signal1, bckgrd1)

    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! Reconstruciton plots for Higgs and W  ----------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    
    #* Higgs

    signal1, _ = createTH1F(signal, "nrightjets", ranges["nrightjets"], "weight", dif = "signal")
    print(signal1.Integral())
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

    values = []
    bckgrd1, _  = createTH1F(signal, "nrightjets", ranges["nrightjets"], "((Jet_genJetIdx[0] != asso_bjet_H_id1 && Jet_genJetIdx[0] != asso_bjet_H_id2)&&(Jet_genJetIdx[1] != asso_bjet_H_id1 && Jet_genJetIdx[1] != asso_bjet_H_id2))*weight", dif = "signal")
    values.append(bckgrd1.Integral())

    bckgrd1, _  = createTH1F(signal, "nrightjets", ranges["nrightjets"], "((Jet_genJetIdx[0] == asso_bjet_H_id1 || Jet_genJetIdx[0] == asso_bjet_H_id2)^(Jet_genJetIdx[1] == asso_bjet_H_id1 || Jet_genJetIdx[1] == asso_bjet_H_id2))*weight", dif = "signal")
    values.append(bckgrd1.Integral())

    bckgrd1, _  = createTH1F(signal, "nrightjets", ranges["nrightjets"], "((Jet_genJetIdx[0] == asso_bjet_H_id1 || Jet_genJetIdx[0] == asso_bjet_H_id2)&&(Jet_genJetIdx[1] == asso_bjet_H_id1 || Jet_genJetIdx[1] == asso_bjet_H_id2))*weight", dif = "signal")
    values.append(bckgrd1.Integral())

    bckgrd1 = ROOT.TH1F("tmp", "tmp", len(ranges["nrightjets"])-1, ranges["nrightjets"])
    for i in range(len(ranges["nrightjets"])-1):
        bckgrd1.SetBinContent(i+1, values[i])

    
    
    bckgrd1.Sumw2()
    
    
    saveratioplot(signal1, bckgrd1, "Number of jets assigned correctly", "Max btag", "Max pt", ratio=False, name = "nrightjets_Higgs")
    del(signal1, bckgrd1)

    #* W

    signal1, _ = createTH1F(signal, "wnrightjets", ranges["wnrightjets"], "weight", dif = "signal")

    values = []
    bckgrd1, _  = createTH1F(signal, "wnrightjets", ranges["wnrightjets"], "((Jet_genJetIdx[0] != asso_jet_W_id1 && Jet_genJetIdx[0] != asso_jet_W_id2)&&(Jet_genJetIdx[1] != asso_jet_W_id1 && Jet_genJetIdx[1] != asso_jet_W_id2))*weight", dif = "signal")
    values.append(bckgrd1.Integral())

    bckgrd1, _  = createTH1F(signal, "wnrightjets", ranges["wnrightjets"], "((Jet_genJetIdx[0] == asso_jet_W_id1 || Jet_genJetIdx[0] == asso_jet_W_id2)^(Jet_genJetIdx[1] == asso_jet_W_id1 || Jet_genJetIdx[1] == asso_jet_W_id2))*weight", dif = "signal")
    values.append(bckgrd1.Integral())

    bckgrd1, _  = createTH1F(signal, "wnrightjets", ranges["wnrightjets"], "((Jet_genJetIdx[0] == asso_jet_W_id1 || Jet_genJetIdx[0] == asso_jet_W_id2)&&(Jet_genJetIdx[1] == asso_jet_W_id1 || Jet_genJetIdx[1] == asso_jet_W_id2))*weight", dif = "signal")
    values.append(bckgrd1.Integral())

    bckgrd1 = ROOT.TH1F("tmp", "tmp", len(ranges["nrightjets"])-1, ranges["nrightjets"])
    for i in range(len(ranges["nrightjets"])-1):
        bckgrd1.SetBinContent(i+1, values[i])
    
    bckgrd1.Sumw2()
    
    
    saveratioplot(signal1, bckgrd1, "Number of jets assigned correctly", "Closest to 80 GeV", "Max pt", ratio=False, name = "nrightjets_W")
    del(signal1, bckgrd1)


    #* try to find alternative procedures with closest angle
    infileasdf = ROOT.TFile.Open("../../signal/breg/Wrecons/signal.root")
    tmp = infileasdf.Get("Events")
    
    signal1, _ = createTH1F(signal, "wnrightjets", ranges["wnrightjets"], "weight", dif = "signal")
    bckgrd1, _  = createTH1F(tmp, "wnrightjetsdelta13", ranges["wnrightjets"], "weight")

    signal1, _ = createTH1F(tmp, "genjets_Hpt", ranges["Jet_pt"], "weight")
    saveplot(signal1, "p_{t} (GeV)", "Reconstructed Higgs boson", name = "app/genjets_Hpt")

    
    
    saveratioplot(signal1, bckgrd1, "Number of jets assigned correctly", "Invariant mass closest to 80 GeV", "Angular distance closest to 1.3", ratio=False, name = "nrightjets_W_delta")
    del(signal1, bckgrd1, tmp, infileasdf)

    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! counting Higgs events --------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------

    signal1, _ = createTH1F(signal, "nrightjets", ranges["nrightjets"], "((Sum$(Jet_btagDeepB>0.4941))>2)", dif = "signal")
    print("Number of events which have more than two bjets: {}".format(signal1.GetEntries()))
    del(signal1)
    signal1, _ = createTH1F(signal, "nrightjets", ranges["nrightjets"], dif = "signal")
    print("Out of {} unweighted events".format(signal1.GetEntries()))
    del(signal1)
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! distros utilized for reconstruction ------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------


    distros = ["genjets_deltaR_W", "genjets_WMass", "genjets_Hmass", "Jet_pt[asso_jet_W_id1]", "Jet_pt[asso_jet_W_id2]",
                "Jet_btagDeepB[asso_bjet_H_id1]", "Jet_btagDeepB[asso_bjet_H_id2]", "genjets_deltaR_H"]


    signal1, _ = createTH1F(signal, distros[1], ranges["recon_Wmass"], "weight", dif = "signal")
    saveplot(signal1, "m_{W} (GeV)", "Generated jets", name = "app/genjets_mW")
    del(signal1)

    
    signal1, _ = createTH1F(signal, distros[3], ranges["Jet_pt"], "weight", dif = "signal")
    signal2, _ = createTH1F(signal, distros[4], ranges["Jet_pt"], "weight", dif = "signal")
    signal1.Add(signal2)
    saveplot(signal1, "p_{t} (GeV)", "Generated jets", name = "app/genjets_ptW")
    del(signal1, signal2)

    signal1, _ = createTH1F(signal, distros[5], ranges["Jet_btagDeepB"], "weight", dif = "signal")
    signal2, _ = createTH1F(signal, distros[6], ranges["Jet_btagDeepB"], "weight", dif = "signal")
    signal1.Add(signal2)
    saveplot(signal1, "b-tagging discriminator output", "Generated jets", name = "app/genjets_btagH")
    del(signal1, signal2)


    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! counting alternative reconstructions -----------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------

    infile5 = ROOT.TFile.Open("../../signal/breg/Wrecons/signal.root")
    signalrecons = infile5.Get("Events")
    
    recons = ["wnrightjets", "wnrightjetsdelta13", "wnrightjets_deltaandmass", "wnrightjets_maxpt", "wnrightjets_maxptandmass",
                "nrightjets_deltaandbtag", "nrightjets_delta", "nrightjets_maxptandbtag", "nrightjets_maxpt", "nrightjets"]
    
    columns = ["W normal", "W delta 1.3", "W delta and normal", "W max pt", "W max pt and normal", "H delta3 and max btag", "H delta3", 
            "H max pt and btag", "H max pt", "H max btag"]

    signal1, _ = createTH1F(signalrecons, distros[0], ranges["deltaR_bb"], "weight", dif = "signal")
    saveplot(signal1, "#DeltaR_{W}", "Generated jets", name = "app/genjets_deltarW")
    del(signal1)



    for index, recon in enumerate(recons):
        finaldict = {}
        tmp = []
        for i in range(3):
            signal1, _ = createTH1F(signalrecons, recon, ranges["wnrightjets"], "({}>{} && {}<{})*weight".format(recon, i-0.5, recon, i+0.5), dif = "signal")
            tmp.append(int(signal1.Integral()))
            del(signal1)
        print(columns[index], np.array(tmp)/float(np.sum(tmp)))        
        #finaldict.update({columns[index]:tmp})

    #export = pd.DataFrame(finaldict)
    #export.to_csv("ReconstructionTable.csv", index=True)


    signal1, _ = createTH1F(signalrecons, distros[7], ranges["deltaR_bb"], "weight", dif = "signal")
    saveplot(signal1, "#DeltaR_{H}", "Generated jets", name = "app/genjets_deltarH")
    del(signal1)



    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! Number of rightly assigned jets after selection and reconstruction -----------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    

    signal1, stats1 = createTH1F(signal ,"nrightjets", ranges["nrightjets"], "weight")
    # saveplot(signal1, "number of rightly assigned jets", title = "H Jets")
    bckgrd1, stats1 = createTH1F(signal ,"wnrightjets", ranges["nrightjets"], "weight")
    saveratioplot(signal1, bckgrd1, "Number of jets assigned correctly", "Higgs", "W-Boson", ratio=False, name = "nrightjets")

    del(signal1, bckgrd1)
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! creating plots of no selection applied and calculate AMS value ---------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    
    var = "recon_Hmass"
    
    signal1, _ = createTH1F(signal_nobreg, var, ranges[var], "weight" , dif = "signal")
    bckgrd1, _ = createTH1F(bckgrd_nobreg, var, ranges[var], "weight")
    saveratioplot(signal1, bckgrd1,  "m_{H} (GeV)",  "H #rightarrow bb", "QCD multijet", log = True, ratio= False, calcAMS = True, name = "recon_Hmass_signal_background")
    signal1.Scale(1/signal1.Integral())
    bckgrd1.Scale(1/bckgrd1.Integral())
    saveratioplot(signal1, bckgrd1, "m_{H} (GeV)", "H #rightarrow bb", "QCD multijet", y_label = "Normalized entries", name = "recon_Hmass_signal_background_normed", ratio = False)
    del(signal1, bckgrd1)

    var = "recon_Wmass"
    
    signal1, _ = createTH1F(signal, var, ranges[var], "weight" , dif = "signal")
    bckgrd1, _ = createTH1F(bckgrd, var, ranges[var], "weight")
    saveratioplot(signal1, bckgrd1,  "m_{W} (GeV)",  "H #rightarrow bb", "QCD multijet", log = True, ratio= False, calcAMS = True, name = "recon_Wmass_signal_background")
    signal1.Scale(1/signal1.Integral())
    bckgrd1.Scale(1/bckgrd1.Integral())
    saveratioplot(signal1, bckgrd1, "m_{W} (GeV)", "H #rightarrow bb", "QCD multijet", y_label = "Normalized entries", name = "recon_Wmass_signal_background_normed", ratio = False)
    del(signal1, bckgrd1)

    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! create plots of genlevel data and reconlevel data ----------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------

    var = "recon_Wmass"
    signal1, _ = createTH1F(signal, var, ranges[var], "weight", dif = "signal")
    bckgrd1, _ = createTH1F(signal, "genjets_WMass", ranges[var], "weight")
    saveratioplot(signal1, bckgrd1, "m_{W}  (GeV)", "Reconstructed level", "Generated level", calcArea = False, ratio = False, name = "Wmass_gen-level_recon-level__ratio")
    del(signal1, bckgrd1)

    var = "recon_Hmass"
    signal1, _ = createTH1F(signal, var, ranges[var], "weight", dif = "signal")
    bckgrd1, _ = createTH1F(signal, "genjets_Hmass", ranges[var], "weight")
    helper1, _ = createTH1F(signal, "original_Hmass", ranges[var], "weight")


    """
    canvas = ROOT.TCanvas("c1","c1")

    axish = createAxisHists(2,hist1,hist1.GetXaxis().GetXmin(),hist1.GetXaxis().GetXmax()-0.01)
    axish[0].GetYaxis().SetRangeUser(0,1.5*max([hist1.GetMaximum(), hist2.GetMaximum(), helper1.GetMaximum()]))
    axish[0].GetYaxis().SetTitle("Entries")

    axish[0].Draw()


    hist1.SetLineColor(ROOT.kRed)
    hist1.SetMarkerStyle(20)
    hist1.SetMarkerColor(ROOT.kRed+3)
    hist2.SetMarkerStyle(22)
    hist2.SetLineColor(ROOT.kBlue)
    hist2.SetMarkerColor(ROOT.kBlue+3)
    helper1.SetLineColor(ROOT.kGreen)
    helper1.SetMarkerColor(ROOT.kGreen+3)
    helper1.SetMarkerStyle(24)

    hist1.Draw("SAME")#both histograms are drawn into the same pad
    hist2.Draw("SAME")
    #helper1.Draw("Same")        

    legend = plot.PositionedLegend(0.4, 0.15, 3, 0.015)
    legend.AddEntry(hist1, "reconstructed level")
    legend.AddEntry(hist2, "generated jets")
    #legend.AddEntry(helper1, "generated particles")

    legend.Draw("SAME")
    
    canvas.SaveAs("Hmass_gen-level_recon-level__ratio.pdf")
    
    canvas.Close()
    del(legend, canvas)
    """





    saveratioplot(signal1, bckgrd1,"m_{H} (GeV)", "Reconstructed level", "Generated level", calcArea = False, ratio = False, name = "Hmass_gen-level_recon-level__ratio")
    del(signal1,bckgrd1, helper1)

    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! trigger table ----------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------


    triggers_considered = ["HLT_AK8PFJet80",  "HLT_PFJet80", "HLT_PFHT300PT30_QuadPFJet_75_60_45_40", 
    "HLT_DiPFJetAve35_HFJEC",  "HLT_PFHT180", "HLT_DoublePFJets40_CaloBTagCSV_p33",  "HLT_AK4PFJet80", "HLT_AK4CaloJet80"]
    trigger_table = pd.DataFrame(columns = ["Trigger", "Signal", "Background"])

    for j in triggers_considered:
        signal1, stats1 = createTH1F(signal,j, np.linspace(0.5,1.5,2),constraints= "("+j+">0)*weight", dif="signal")
        ints1 = signal1.Integral()
        del(signal1)
        signal2, stats1 = createTH1F(signal,j, np.linspace(-0.5,1.5,2),constraints= "weight")
        ints2 = signal2.Integral()
        del(signal2)
        bckgrd1, trash = createTH1F(bckgrd, j, np.linspace(0.5,1.5,2), constraints= "("+j+">0)*weight")
        intb1 = bckgrd1.Integral()
        del(bckgrd1)
        bckgrd2, trash = createTH1F(bckgrd,j, np.linspace(-0.5,1.5,2), constraints= "weight")
        intb2 = bckgrd2.Integral()
        print(ints1/ints2, intb1/intb2)
        trigger_table = trigger_table.append({"Trigger": j, "Signal": ints1/ints2, "Background":intb1/intb2}, ignore_index=True)
    
    trigger_table.to_latex("HLT_selection.txt", index = False)


    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! shape comparison trigger -----------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------

    sim = {"HLT_PFHT300PT30_QuadPFJet_75_60_45_40":"HLT_Quad", "HLT_PFHT180":"HLT_PF"}
    sim2 = {"HLT_PFHT300PT30_QuadPFJet_75_60_45_40":"Quad trigger", "HLT_PFHT180":"HT180 trigger"}


    i = "recon_Hmass"
    for j in selected_triggers:
        tmp1, stats1 = createTH1F(signal, i, ranges[i],constraints= "weight", normalise = True, dif="signal")
        tmp2, trash = createTH1F(signal, i, ranges[i], constraints= "("+j+">0)*weight", normalise = True)
        saveratioplot(tmp2, tmp1, "m_{H} (GeV)",  sim2[j], "no trigger",  dif = j, calcArea=False, name = "reconHmass{}shape".format(sim[j]), y_label = "Normalized entries")
        del(tmp1,tmp2)


    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! calcualate AMS value for triggers --------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------

    var = "recon_Hmass"

    for j in selected_triggers:
        signal1, _ = createTH1F(signal, var, ranges[var], "({}>0)*weight".format(j) , dif = "signal")
        bckgrd1, _ = createTH1F(bckgrd, var, ranges[var], "({}>0)*weight".format(j))
        print("{}: Signal {}, Backround {}".format(j, signal1.Integral(), bckgrd1.Integral()) )

        saveratioplot(signal1, bckgrd1,  "m_{H} (GeV)",  "H #rightarrow bb", "QCD multijet", log = True, ratio= False, calcAMS = True, name = "recon_Hmass_{}".format(sim[j]))
        del(signal1, bckgrd1)

    # var = "recon_Wmass"

    # for j in selected_triggers:
    #     signal1, _ = createTH1F(signal, var, ranges[var], "({}>0)*weight".format(j) , dif = "signal")
    #     bckgrd1, _ = createTH1F(bckgrd, var, ranges[var], "({}>0)*weight".format(j))
    #     saveratioplot(signal1, bckgrd1,  "reconstructed m_{W} (GeV)",  "signal", "background", log = True, ratio= False, calcAMS = True, name = "recon_Wmass_{}".format(sim[j]))
    #     del(signal1, bckgrd1)

    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! single cut plots -------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------

    var = ["deltaR_W", "recon_Wpt", "deltaR_reconHW"]
    sims = {"deltaR_W":"#DeltaR_{W}", "recon_Wpt":"#p_{T} of W (GeV)", "deltaR_reconHW":"#DeltaR_{HW}"}
    #ranges.update({"Jet_pt[recon_WJets_id[0]]":np.linspace(0, 200, 30)})
    names = {"deltaR_W":"deltaR_W", "recon_Wpt":"recon_Wpt", "deltaR_reconHW":"deltaR_reconHW"}

    for j in selected_triggers:
        for va in var:
            signal1, _ = createTH1F(signal, va, ranges[va], "({}>0)*weight".format(j) , dif = "signal", normalise = True)
            bckgrd1, _ = createTH1F(bckgrd, va, ranges[va], "({}>0)*weight".format(j), normalise = True)
            saveratioplot(signal1, bckgrd1,  sims[va],  "H #rightarrow bb", "QCD multijet", log = False, ratio= False, calcArea=True, calcAMS = False, name = "{}_{}".format(names[va], sim[j]), y_label = "Normalized entries")
            del(signal1, bckgrd1)


    cuttingvars = {"Jet_pt[recon_WJets_id[0]]>":[0, 200], "deltaR_reconHW<":[0.5,4.5], "deltaR_W<": [0,7]}
    sim = {"deltaR_W<":"#DeltaR_{W}", "deltaR_reconHW<":"#DeltaR_{HW}", "Jet_pt[recon_WJets_id[0]]>":"p_{T} of high energy W jet (GeV)"}
    #cuttingvars = {"Jet_HT>":[100,600],"recon_Wpt>":[0, 250], "recon_Hpt>":[50,250]}
    
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! cutplots ---------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    
    cuttingvars = {"deltaR_W<":[0, 6], "deltaR_reconHW<":[0.5,4.5], "recon_Wpt<": [0,500]}
    sim = {"deltaR_W<":"#DeltaR_{W}", "deltaR_reconHW<":"#DeltaR_{HW}", "recon_Wpt<":"p_{T} of W boson"}#"Jet_pt[recon_WJets_id[0]]>":"p_{T} of high energy W jet (GeV) <"}
    sim_trigger = {"HLT_PFHT300PT30_QuadPFJet_75_60_45_40":"HLT_Quad", "HLT_PFHT180":"HLT_PF"}
    


    for cuttingvar in cuttingvars.keys():
        bestcut = {}
        for j in selected_triggers:
            tmp = []
            for i in np.linspace(cuttingvars[cuttingvar][0],cuttingvars[cuttingvar][1],30):
                    signal1, stats1 = createTH1F(signal, "recon_Hmass", ranges["recon_Hmass"], constraints= "("+cuttingvar + str(i) +"&&"+j+ ">0)*weight", dif="signal")
                    bckgrd1, trash = createTH1F(bckgrd, "recon_Hmass", ranges["recon_Hmass"], constraints= "("+cuttingvar + str(i) +"&&"+j+ ">0)*weight")
                    s = calc_AMS(signal1, bckgrd1)
                    if not(np.isinf(s)):
                        tmp.append([i,s])
                    del(signal1, bckgrd1)
            bestcut.update({j:np.array(tmp)})
        
        
        for i in bestcut.keys():
            maximum1 = max(enumerate(bestcut[i][:,1]), key=itemgetter(1))
            
            savegraph(bestcut[i][:,0], bestcut[i][:,1], "cutplot_{}_{}".format(cuttingvar, sim_trigger[i]), sim[cuttingvar], 
            "Binned s/#sqrt{b}: ", text = "Max {} at {}".format(maximum1[1], round(bestcut[i][maximum1[0],0],2)))
    


    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! all cuts applied plot --------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------

    allgoodcuts_quad = "recon_Wmass>50&&recon_Wmass<100 && Jet_btagDeepB[recon_WJets_id[0]]<0.85 && Jet_btagDeepB[recon_WJets_id[1]]<0.9 && deltaR_reconHW<3.5 && recon_Wpt>175&& deltaR_W<1.16"
    allgoodcuts_pf = "recon_Wmass>50&&recon_Wmass<100 && Jet_btagDeepB[recon_WJets_id[0]]<0.85 && Jet_btagDeepB[recon_WJets_id[1]]<0.9 && deltaR_reconHW<3.5 && recon_Wpt>175&& deltaR_W<1.16"

    va = "recon_Hmass"

    
    
    signal1, _ = createTH1F(signal, va, ranges[va], "({} && {}>0)*weight".format(allgoodcuts_quad, selected_triggers[0]) , dif = "signal")
    bckgrd1, _ = createTH1F(bckgrd, va, ranges[va], "({} && {}>0)*weight".format(allgoodcuts_quad, selected_triggers[0]))

    #print(signal1.Integral(), bckgrd1.Integral())
    saveratioplot(signal1, bckgrd1,  "m_{H} (GeV)",  "H #rightarrow bb", "QCD multijet", log = True, ratio= False, calcArea=False, calcAMS = True, name = "allgoodcuts_{}".format(sim_trigger[selected_triggers[0]]))
    del(signal1, bckgrd1)


    signal1, _ = createTH1F(signal, va, ranges[va], "({} && {}>0)*weight".format(allgoodcuts_pf, selected_triggers[1]) , dif = "signal")
    bckgrd1, _ = createTH1F(bckgrd, va, ranges[va], "({} && {}>0)*weight".format(allgoodcuts_pf, selected_triggers[1]))

    #print(signal1.Integral(), bckgrd1.Integral())
    saveratioplot(signal1, bckgrd1,  "m_{H} (GeV)",  "H #rightarrow bb", "QCD multijet", log = True, ratio= False, calcArea=False, calcAMS = True, name = "allgoodcuts_{}".format(sim_trigger[selected_triggers[1]]))
    del(signal1, bckgrd1)


    # triggers = ["HLT_PFHT300PT30_QuadPFJet_75_60_45_40", "HLT_PFHT180"]
    # ranges = np.linspace(75, 160, 60)

    # for j in triggers:

    #     infile = ROOT.TFile.Open("tmvaBDT"+j+".root")
        
        
    #     tree = infile.Get("dataset/TestTree")
    #     tree2 = infile.Get("dataset/TrainTree")

    #     os.chdir("AMS_onrecon_Hmass")
    #     makecutgraph_fortmva(tree, tree2, "recon_Hmass", "BDT>", ranges,np.linspace(0,1,100), j)
    #     infile.Close()
    #     os.chdir("../")
    #     del(infile, tree, tree2)

    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! get them areas ---------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------

    variables = ["recon_Wmass", "recon_Hmass", "deltaR_reconHW","recon_Hpt","recon_Wpt","deltaR_bb", "deltaR_W", 
    "Jet_HT","bJet_pt1","bJet_pt2","WJet_pt1","WJet_pt2","bJet_btagDeepB1","bJet_btagDeepB2","WJet_btagDeepB1","WJet_btagDeepB2"]

    ranges.update({"bJet_pt1":np.linspace(0,400, 30), "bJet_pt2":np.linspace(0,400, 50), "WJet_pt1":np.linspace(0,400, 50), "WJet_pt2":np.linspace(0,400, 50),
    "bJet_btagDeepB1":np.linspace(0, 1, 30), "bJet_btagDeepB2":np.linspace(0, 1, 30), "WJet_btagDeepB1":np.linspace(0, 1, 30), "WJet_btagDeepB2":np.linspace(0, 1, 30)
    })
    
    for j in selected_triggers:
        for va in variables:
            signal1, _ = createTH1F(signal, va, ranges[va], "({}>0)*weight".format(j) , dif = "signal", normalise = True)
            bckgrd1, _ = createTH1F(bckgrd, va, ranges[va], "({}>0)*weight".format(j), normalise = True)
            value = saveratioplot(signal1, bckgrd1, "", "", "", calcArea=True, return_area=True)
            print(j, va, value)

    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! overlapping areas tables -----------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------

    #for futher input W Jet_neEmEF, W Jet_neHEF, W Jet_qgl, ?H Jet_qgl, ? W Jet_nConstituents ?

    

    varsused = ["recon_Wmass", "recon_Hmass", "deltaR_reconHW","recon_Hpt","recon_Wpt","deltaR_bb", "deltaR_W", 
    "Jet_HT","!HJet_pt", "!WJet_pt", "!WJet_neEmEF", "!WJet_neHEF", "!WJet_qgl", "!HJet_qgl",
    "!WJet_nConstituents", "!WJet_btagDeepB", "!HJet_btagDeepB", "!WJet_btagDeepC", "!WJet_mass", "!HJet_mass"]

    varsconsid = ["!HJet_nMuons", "!WJet_nMuons", "!HJet_nElectrons", "!WJet_nElectrons", "!HJet_chHEF", "!WJet_chHEF", "!HJet_btagDeepC", "!HJet_neHEF", "!HJet_neEmEF", "!HJet_nConstituents"]

    

    ranges.update({"Jet_chEmEF":np.linspace(0,1,20), "Jet_chHEF":np.linspace(0,1,20), "Jet_nConstituents":np.linspace(1, 60, 31),
    "Jet_neEmEF":np.linspace(0,1,20), "Jet_neHEF":np.linspace(0,1,20),"Jet_qgl":np.linspace(-1,1,30),
    "Jet_nElectrons":np.linspace(0,5,6), "Jet_nMuons":np.linspace(0,5,6)})

    nice_plot_dict = {"!WJet_btagDeepC":"btagDeepC of W-jets", 
    "!HJet_btagDeepC":"btagDeepC of H-jets", "!HJet_mass":"Mass of H-jets (GeV)", 
    "!WJet_mass":"Mass of W-jets (GeV)", "!HJet_btagDeepB":"btagDeepB for H-jets",
    "!WJet_btagDeepB":"btagDeepB for W-jets", "recon_Wmass":"m_{W} (GeV)", "recon_Hmass":"m_{H} (GeV)",
    "deltaR_reconHW":"#DeltaR_{HW}", "recon_Hpt":"p_{T} of H (GeV)", "recon_Wpt":"p_{T} of W boson (GeV)",
    "deltaR_bb":"#DeltaR_{H}", "deltaR_W":"#DeltaR_{W}","Jet_HT":"Summed E_{T} (Gev)",
    "!WJet_pt":"p_{T} of W-jets (GeV)", "!HJet_pt":"p_{T} of W-jets (GeV)", 
    "bJet_pt1":"p_{T} of H-jets higher", "bJet_pt2":"p_{T} of H-jets lower",
    "!WJet_neEmEF":"Neutral Electromagnetic Energy Fraction for W-jets",
    "!WJet_neHEF":"Neutral Hadron Energy Fraction for W-jets", 
    "!WJet_qgl":"Quark vs Gluon discriminator for W-jets",
    "!HJet_qgl":"Quark vs Gluon discriminator for H-jets",
    "!WJet_nConstituents":"Number of jet constituents for W-jets", "!WJet_chHEF":"Charged hadronic energy fraction for W-jets"}


    # comparison
    
    for curtab in [varsused, varsconsid]:
        dataframe = pd.DataFrame(columns = ["var", "Quad", "HT180"])

        for va in curtab:
            area = []
            for trigger in selected_triggers:
                if va[0]=="!":
                    #print(va)
                    if va[1]=="W":
                        helper = "recon_WJets_id"
                    elif va[1]=="H":
                        helper = "recon_bjet_id"
                    signal1, _ = createTH1F(signal, va[2:]+"[{}[0]]".format(helper), ranges[va[2:]], "({}>0)*weight".format(trigger) , dif = "signal")
                    bckgrd1, _ = createTH1F(bckgrd, va[2:]+"[{}[0]]".format(helper), ranges[va[2:]], "({}>0)*weight".format(trigger))
                    signal1.Scale(1/signal1.Integral())
                    bckgrd1.Scale(1/bckgrd1.Integral())
                    
                    signal2, _ = createTH1F(signal, va[2:]+"[{}[1]]".format(helper), ranges[va[2:]], "({}>0)*weight".format(trigger))
                    bckgrd2, _ = createTH1F(bckgrd, va[2:]+"[{}[1]]".format(helper), ranges[va[2:]], "({}>0)*weight".format(trigger))
                    signal2.Scale(1/signal2.Integral())
                    bckgrd2.Scale(1/bckgrd2.Integral())
                    

                    area.append([saveratioplot(signal1, bckgrd1, va,  "signal", "background", return_area=True, calcArea=True), 
                    saveratioplot(signal2, bckgrd2, va,  "signal", "background", return_area=True, calcArea=True)])
                    if area[-1][-1]<0.9:
                        saveratioplot(signal2, bckgrd2, nice_plot_dict[va],  "H #rightarrow bb", "QCD multijet", ratio=False, return_area=False, calcArea=True, name = "cuts/distributions/"+va+sim_trigger[trigger])
                    

                    del(signal1, bckgrd1 , signal2, bckgrd2)
                else:
                    signal1, _ = createTH1F(signal, va, ranges[va], "({}>0)*weight".format(trigger) , dif = "signal", normalise=True)
                    bckgrd1, _ = createTH1F(bckgrd, va, ranges[va], "({}>0)*weight".format(trigger), normalise=True)

                    area.append(saveratioplot(signal1, bckgrd1, va,  "signal", "background", return_area=True, calcArea=True))
                    if area[-1]<0.9:
                        saveratioplot(signal1, bckgrd1, nice_plot_dict[va],  "H #rightarrow bb", "QCD multijet", ratio=False, return_area=False, calcArea=True, name = "cuts/distributions/"+va+sim_trigger[trigger])

                    del(signal1, bckgrd1)

            dataframe = dataframe.append({"var":va, "Quad":area[0], "HT180":area[1]}, ignore_index =True)        
        if curtab == varsconsid:
            dataframe.to_latex("notusedvars_latex.txt")
        if curtab == varsused:
            dataframe.to_latex("varsused_latex.txt")
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! variable distributions for Appendix ------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    
    #! ->see overlappingarea table
    """

    os.chdir('app')

    further_variables = {#"!HJet_hadronFlavour":np.linspace(-0.5, 6.5, 8),
    "!WJet_btagDeepC":np.linspace(0,1,20), "!HJet_mass":np.linspace(0,300, 30), "!WJet_mass":np.linspace(0,300, 30),
    "!HJet_btagDeepB":np.linspace(0,1,20), "!WJet_btagDeepB":np.linspace(0,1,20)}
    ranges.update(further_variables)

    nice_plot_dict = {"!WJet_btagDeepC":"btagDeepC of W-jets", 
    "!HJet_btagDeepC":"btagDeepC of H-jets", "!HJet_mass":"Mass of H-jets (GeV)", 
    "!WJet_mass":"Mass of W-jets (GeV)", "!HJet_btagDeepB":"btagDeepB for H-jets",
    "!WJet_btagDeepB":"btagDeepB for W-jets", "recon_Wmass":"m_{W} (GeV)", "recon_Hmass":"m_{H} (GeV)",
    "deltaR_reconHW":"#DeltaR_{HW}", "recon_Hpt":"p_{T} of H (GeV)", "recon_Wpt":"p_{T} of W (GeV)",
    "deltaR_bb":"#DeltaR_{H}", "deltaR_W":"#DeltaR_{W}","Jet_HT":"Summed E_{T} (Gev)",
    "!WJet_pt":"p_{T} of W-jets (GeV)", "!HJet_pt":"p_{T} of W-jets (GeV)", 
    "bJet_pt1":"p_{T} of H-jets higher", "bJet_pt2":"p_{T} of H-jets lower",
    "!WJet_neEmEF":"Neutral Electromagnetic Energy Fraction for W-jets",
    "!WJet_neHEF":"Neutral Hadron Energy Fraction for W-jets", 
    "!WJet_qgl":"Quark vs Gluon likelihood discriminator for W-jets",
    "!HJet_qgl":"Quark vs Gluon likelihood discriminator for H-jets",
    "!WJet_nConstituents":"Number of jet constituents for W-jets"}




    all_vars_consid = np.array(further_variables.keys())
    
    all_vars_consid = np.append(all_vars_consid, ["recon_Wmass", "recon_Hmass", "deltaR_reconHW","recon_Hpt","recon_Wpt","deltaR_bb", "deltaR_W", 
            "Jet_HT","bJet_pt1", "bJet_pt2","!WJet_pt"])
    
    ranges.update({"bJet_pt1":ranges["Jet_pt"], "bJet_pt2":ranges["Jet_pt"]})

    #print([item for sublist in all_vars_consid for item in sublist])
    all_vars_consid = np.array(all_vars_consid)
    all_vars_consid = all_vars_consid.flatten()
    #print(all_vars_consid)

    for va in varsused:
        if va in ["deltaR_HW", "deltaR_W", "deltaR_bb"]:
            continue
        if va[0]=="!":
            #print(va)
            if va[1]=="W":
                helper = "recon_WJets_id"
            elif va[1]=="H":
                helper = "recon_bjet_id"
            

            signal1, _ = createTH1F(signal, va[2:]+"[{}[0]]".format(helper), ranges[va[2:]], "weight" , dif = "signal")
            signal2, _ = createTH1F(signal, va[2:]+"[{}[1]]".format(helper), ranges[va[2:]], "weight")
            signal1.Add(signal2)
            signal1.Scale(1/signal1.Integral())
            del(signal2)

            bckgrd1, _ = createTH1F(bckgrd, va[2:]+"[{}[0]]".format(helper), ranges[va[2:]], "weight")
            bckgrd2, _ = createTH1F(bckgrd, va[2:]+"[{}[1]]".format(helper), ranges[va[2:]], "weight")
            bckgrd1.Add(bckgrd2)
            bckgrd1.Scale(1/bckgrd1.Integral())
            del(bckgrd2)
        else:
            #print(va)
            signal1, _ = createTH1F(signal, va, ranges[va], "weight" , dif = "signal", normalise=True)
            bckgrd1, _ = createTH1F(bckgrd, va, ranges[va], "weight", normalise=True)

        saveratioplot(signal1, bckgrd1, nice_plot_dict[va],  "signal", "background", log = False, ratio= False, calcArea=True, calcAMS = False, name = "ratio_{}".format(va))
        del(signal1, bckgrd1)
    
    """
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! BDT final --------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------

    os.chdir("../../classifiers")

    #triggers = ["HLT_PFHT300PT30_QuadPFJet_75_60_45_40", "HLT_PFHT180"]
    ranges = np.linspace(90, 140, 15)

    for j in selected_triggers:

        infile = ROOT.TFile.Open("tmvaBDT"+j+".root")
        
        
        tree = infile.Get("dataset/TestTree")
        #print(tree.GetEntries("BDT"))
        tree2 = infile.Get("dataset/TrainTree")


        os.chdir("../plotting/pics/BDT")
        makecutgraph_fortmva(tree, tree2, "recon_Hmass", "BDT>", ranges, 30 , j)


        infile.Close()
        os.chdir("../../../classifiers")
        del(infile, tree, tree2)


   