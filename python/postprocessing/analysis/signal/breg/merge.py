import ROOT
import numpy as np
import os

files = os.listdir(os.getcwd())
print(files)
files.remove('merge.py')
if 'signal.root' in files:
	files.remove('signal.root')
liste = ROOT.TList()

ch = ROOT.TChain("Events")

final = 0
for fl in files:
    ch.Add(fl)
    
    #infile = ROOT.TFile.Open(fl, "read")
    #tree = infile.Get("Events")
    #final+=tree.GetEntries()
    ##liste.Add(infile.Get("Events"))

#print(final)
ch.Merge("signal.root")

#fout = ROOT.TFile("bckgrd.root", "recreate")

#final = ROOT.TTree.MergeTrees(liste)
#print(final.GetEntries())
##final.SetName("bckgrd")
#fout.cd()
#final.Write()
#fout.Close()
    
    
