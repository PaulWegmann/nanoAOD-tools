import ROOT
import os
import numpy as np
# import sys
# sys.path.append('plotting')
# from plots import creatTH1F

os.chdir("/pnfs/desy.de/cms/tier2/store/user/adewit/SummerProjectSamples/")

files = ["QCDbEnr_HT100to200_file1.root", "QCDbEnr_HT100to200_file2.root", "QCDbEnr_HT100to200_file3.root", 
"QCDbEnr_HT200to300_file1.root", "QCDbEnr_HT200to300_file2.root", "QCDbEnr_HT200to300_file3.root", "QCDbEnr_HT300to500_file1.root", 
"QCDbEnr_HT300to500_file2.root", "QCDbEnr_HT500to700.root", "QCDbEnr_HT700to1000.root", "QCDbEnr_HT1000to1500_file1.root",
 "QCDbEnr_HT1000to1500_file2.root", "QCDbEnr_HT1500to2000.root", "QCDbEnr_HT2000toInf.root"]
#files = ["WminusH_HToBB_WToQQ_M125.root", "WplusH_HToBB_WToQQ_M125.root"]

var = "recon_Hmass"
ranges = np.linspace(75, 160, 60)

weights = {"HT100to200":1.127e+06,"HT200to300":8.073e+04, "HT300to500":1.668e+04, "HT500to700":1.485e+03*1.9, "HT700to1000":2.976e+02*1.9, 
"HT1000to1500":4.640e+01, "HT1500to2000":3.720e+00, "HT2000toInf":6.438e-01, "Wminus":5.850e-01, "Wplus":5.850e-01}

runcount = {"QCDbEnr_HT100to200_file1.root":1567593, "QCDbEnr_HT100to200_file2.root":1371694, 
"QCDbEnr_HT100to200_file3.root":1385847, "QCDbEnr_HT200to300_file1.root":1259684, "QCDbEnr_HT200to300_file2.root":1349109, 
"QCDbEnr_HT200to300_file3.root":1160070, "QCDbEnr_HT300to500_file1.root":1100312, "QCDbEnr_HT300to500_file2.root":843531, 
"QCDbEnr_HT500to700.root":968568, "QCDbEnr_HT700to1000.root":510681, "QCDbEnr_HT1000to1500_file1.root":49559,
 "QCDbEnr_HT1000to1500_file2.root":46801, "QCDbEnr_HT1500to2000.root":187594, "QCDbEnr_HT2000toInf.root":151292,
  "WminusH_HToBB_WToQQ_M125":484662, "WplusH_HToBB_WToQQ_M125":504997}
luminosity = 58830

multiplier = 0

totalevents = 0
totalweighted = 0
totalpreselect = 0
for f in files:
    for s in weights.keys():
        if s in f:
            multiplier = weights[s]

    infile = ROOT.TFile.Open(f)
    data = infile.Get("Events")
    
    check = ROOT.TH1F("tmp", "tmp", 10, 0, 10000)

    #data.Draw("LHE_HT>>tmp", "(((Sum$((abs(Jet_eta)<2.5 && Jet_pt > 20 && Jet_jetId>0)) >= 2))&&(Sum$(Electron_pt > 20 && Electron_mvaFall17V2Iso_WP90)<1)&&(Sum$(Muon_pt>20 && Muon_tightId)<1) &  Sum$(Jet_pt>20 & Jet_eta<4.8 & Jet_jetId>0)>3)*{}".format(multiplier*luminosity/data.GetEntries())) #Sum$(Jet_pt>20 & abs(Jet_eta)<2.5 & Jet_btagDeepB>0.4941 & Jet_jetId>0)>1 &
    #data.Draw("LHE_HT>>tmp", "(Sum$(Jet_btagDeepB>0.4941) >= 2)*{}".format(multiplier*luminosity/data.GetEntries())) #Sum$(Jet_pt>20 & abs(Jet_eta)<2.5 & Jet_btagDeepB>0.4941 & Jet_jetId>0)>1 &
    preselection_cuts = check.Integral()
    totalpreselect+=preselection_cuts
    del(check)

    totalevents += data.GetEntries()
    totalweighted += multiplier*luminosity

    print("{} with {} Entries and weighted {} and preselection results to {}".format(f, data.GetEntries(), multiplier*luminosity, preselection_cuts))

print("{} entries and {} weighted in total and preselection {}".format(totalevents, totalweighted, totalpreselect))