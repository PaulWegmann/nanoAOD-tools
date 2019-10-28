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


path = "file:///pnfs/desy.de/cms/tier2/store/user/adewit/SummerProjectSamples/"
files = ["QCDbEnr_HT100to200_file1.root", "QCDbEnr_HT100to200_file2.root", "QCDbEnr_HT100to200_file3.root", "QCDbEnr_HT200to300_file1.root", "QCDbEnr_HT200to300_file2.root", "QCDbEnr_HT200to300_file3.root", "QCDbEnr_HT300to500_file1.root", "QCDbEnr_HT300to500_file2.root", "QCDbEnr_HT500to700.root", "QCDbEnr_HT700to1000.root", "QCDbEnr_HT1000to1500_file1.root", "QCDbEnr_HT1000to1500_file2.root", "QCDbEnr_HT1500to2000.root", "QCDbEnr_HT2000toInf.root"]
#files = ["WminusH_HToBB_WToQQ_M125.root", "WplusH_HToBB_WToQQ_M125.root"]
for i in range(len(files)):
    files[i] = path+files[i]

# ||(Sum$((abs(FatJet_eta)<2.5 && FatJet_pt > 200 && FatJet_jetId)) >= 1)
selection='''((Sum$((abs(Jet_eta)<2.5 && Jet_pt > 20 && Jet_jetId>0)) >= 2))&&(Sum$(Electron_pt > 20 && Electron_mvaFall17V2Iso_WP90)<1)&&(Sum$(Muon_pt>20 && Muon_tightId)<1)'''

#selection = ''''''
#what is the last bit?

p=PostProcessor(".",files,selection.replace('\n',' '),branchsel="keep_and_drop.txt",modules=[basesel()],provenance=True,outputbranchsel="keep_and_drop.txt")

p.run()
