#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

#from PhysicsTools.NanoAODTools.postprocessing.analysis.higgs.vhbb.VHbbProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
from  PhysicsTools.NanoAODTools.postprocessing.modules.jme.mht import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.BaselineSelection import *

#for importing file in different folder
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.append('./plotting')

from plots import createTH1F, saveplot

# parameter setting ---------------------------------------------
path = "file:///pnfs/desy.de/cms/tier2/store/user/adewit/SummerProjectSamples/"

filesS = ["WminusH_HToBB_WToQQ_M125.root", "WplusH_HToBB_WToQQ_M125.root"]
files = ["QCDbEnr_HT100to200_file2", "QCDbEnr_HT100to200_file3", "QCDbEnr_HT200to300_file2", "QCDbEnr_HT300to500_file2"]
#files = ["QCDbEnr_HT100to200_file1", "QCDbEnr_HT100to200_file2", "QCDbEnr_HT100to200_file3", "QCDbEnr_HT200to300_file1", "QCDbEnr_HT200to300_file2", "QCDbEnr_HT200to300_file3", "QCDbEnr_HT300to500_file1", "QCDbEnr_HT300to500_file2", "QCDbEnr_HT500to700", "QCDbEnr_HT700to1000", "QCDbEnr_HT1000to1500_file1", "QCDbEnr_HT1000to1500_file2", "QCDbEnr_HT1500to2000", "QCDbEnr_HT2000toInf"]


run =  ["HT100to200","HT200to300", "HT300to500", "HT500to700", "HT700to1000", "HT1000to1500", "HT1500to2000", "HT2000toInf"]

weights = {"HT100to200":1.127e+06,"HT200to300":8.073e+04, "HT300to500":1.668e+04, "HT500to700":1.485e+03, "HT700to1000":2.976e+02, "HT1000to1500":4.640e+01, "HT1500to2000":3.720e+00, "HT2000toInf":6.438e-01}

signal = []

luminosity = 58830

selection='''((Sum$((abs(Jet_eta)<2.5 && Jet_pt > 20 && Jet_jetId>0)) >= 2)||(Sum$((abs(FatJet_eta)<2.5 && FatJet_pt > 200 && FatJet_jetId)) >= 1))&&(Sum$(Electron_pt > 20 && Electron_mvaFall17V2Iso_WP90)<1)&&(Sum$(Muon_pt>20 && Muon_tightId)<1)'''



# analysing and saving proberly --------------------------------
for i in files:
    #p=PostProcessor(".",[path+i+".root"],selection.replace('\n',' '),branchsel="keep_and_drop.txt",modules=[basesel()],provenance=True,outputbranchsel="keep_and_drop.txt")
    #p.run()
    
    infile = ROOT.TFile.Open(i+"_Skim.root")
    for j in run:
        if j in i:
            current = j
    
    weight = weights[current]
    
    fout = ROOT.TFile("bckgrd/"+i+".root","recreate")
    fout.cd()
    
    getcount = infile.Get("Runs") # it would be nicer to get number of entries out of this variable getcount -> genEventCount
    getcount.GetEntry()
    entries = getcount.genEventCount 
    
    tree = infile.Get("Events")
    #tree.AddBranch(luminosity*weight/entries)
    tree.Write()
    
    tmp1, trash = createTH1F(tree, "recon_Hmass", np.linspace(80,160, 50))
    tmp2, trash = createTH1F(tree, "recon_WMass", np.linspace(30,160, 70))
    
    saveplot(tmp1, "recon_Hmass", current)
    saveplot(tmp2, "recon_WMass",  current)
    fout.Close()
    
    del(tmp1,tmp2, tree, infile, fout)
    
#fout.Close()
