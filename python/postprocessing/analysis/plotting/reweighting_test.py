import ROOT
import os
import numpy as np

##for importing file in different folder
#import sys
## insert at 1, 0 is the script path (or '' in REPL)
#sys.path.append('./plotting')

from plots import saveplot, createTH1F, saveratioplot

#changing location to file directory
pathfiles = "/pnfs/desy.de/cms/tier2/store/user/adewit/SummerProjectSamples/"
pathsave = "/nfs/dust/cms/user/wegmannp/CMSSW_10_4_0/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/plotting"
os.chdir(pathfiles)

files = ["QCDbEnr_HT100to200_file1.root", "QCDbEnr_HT100to200_file2.root", "QCDbEnr_HT100to200_file3.root", "QCDbEnr_HT200to300_file1.root", "QCDbEnr_HT200to300_file2.root", "QCDbEnr_HT200to300_file3.root", "QCDbEnr_HT300to500_file1.root", "QCDbEnr_HT300to500_file2.root", "QCDbEnr_HT500to700.root", "QCDbEnr_HT700to1000.root", "QCDbEnr_HT1000to1500_file1.root", "QCDbEnr_HT1000to1500_file2.root", "QCDbEnr_HT1500to2000.root", "QCDbEnr_HT2000toInf.root"]

luminosity = 58830

weights = {"HT100to200":1.127e+06,"HT200to300":8.073e+04, "HT300to500":1.668e+04, "HT500to700":1.485e+03, "HT700to1000":2.976e+02, "HT1000to1500":4.640e+01, "HT1500to2000":3.720e+00, "HT2000toInf":6.438e-01}
#hist = createTH1F("hist", "hist", 100, 0, 3000)


#for fl in files:
infile = ROOT.TFile.Open(files[0])
tree = infile.Get("Events")
histo, trash = createTH1F(tree, "LHE_HT", np.linspace(200,1500,30))
histo2, trash = createTH1F(tree, "LHE_HT", np.linspace(200,1500,30), dif ="weight")


getcount = infile.Get("Runs")
getcount.GetEntry()
entries = getcount.genEventCount
print(entries==tree.GetEntries())

#histo.Scale(weights["HT100to200"]/entries*luminosity) # to make the pdf accordingly for every section

for fl in files[1:]:
    area = ""
    for i in weights.keys():
        if i in fl:
            area=i
            break
    infile1 = ROOT.TFile.Open(fl)
    tree = infile1.Get("Events")
    
    getcount = infile1.Get("Runs") # it would be nicer to get number of entries out of this variable getcount -> genEventCount
    getcount.GetEntry()
    entries = getcount.genEventCount 
    print(entries== tree.GetEntries())
    #print(getcount.GetLeaf("genEventCount").GetValue(0))
    tmp2, trash = createTH1F(tree, "LHE_HT", np.linspace(200,1500,30), dif=fl)
    tmp1, trash = createTH1F(tree, "LHE_HT", np.linspace(200,1500,30))

    if "HT500to700" in fl or "HT700to1000" in fl :
        tmp1.Scale(weights[area]*1.9/entries*luminosity)
    else:
        tmp1.Scale(weights[area]/entries*luminosity)
    histo.Add(tmp2)
    histo2.Add(tmp1)

#tree.Draw("LHE_HT>>hist")
    
#hist.Draw("")
os.chdir(pathsave)
print(histo.GetEntries())
saveplot(histo, "Summed E_{T} (GeV)", log = True)
saveratioplot(histo2, histo, "Summed E_{T}", "weighted", "unweighted", log=True, ratio=False, name="weighting")
