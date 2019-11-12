from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Event
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import InputTree
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
import ROOT
import os
import numpy as np
from plots import saveplot
# import sys
# sys.path.append('plotting')
# from plots import creatTH1F


def deltaR(b, jets): #return min deltaR?
    delta = []
    for jet in jets:
        delta.append(np.sqrt((b.eta-jet.eta)**2+(b.phi-jet.phi)**2))
        
    val, idx = min((val, idx) for (idx, val) in enumerate(delta))
    return(val, idx)

def analyze(event, hist):
    
    bs = []
    wp = [] #"w+ particles"

    
    genparticles = Collection(event, "GenPart")
    genjets = Collection(event, "GenJet")
    
           
            
    for genpart in genparticles:
        #selection of b particles originating from Higgs
        if(abs(genpart.pdgId)==5 and genpart.genPartIdxMother>-1):
            if genparticles[genpart.genPartIdxMother].pdgId==25 and genparticles[genpart.genPartIdxMother].statusFlags&8192==8192: #statusflag is used too make sure it is the "last copy of the particle" due to the nature of the generation algorhythm which often copies particles multiple times to simulate decay
                bs.append(genpart)
        #selection of b particles originating from W+
        #if genpart.genPartIdxMother>-1:
        #    if genparticles[genpart.genPartIdxMother].pdgId==particleid and genparticles[genpart.genPartIdxMother].statusFlags&8192==8192:
        #        wp.append(genpart)

    if len(bs)!=2:
                if prterrmsg:
                    print("not both particles have been found")
                return False

    try:
        delta1, index1 = deltaR(bs[0], genjets) 
        delta2, index2 = deltaR(bs[1], genjets)
        #delta3, index3 = deltaR(wp[0], genjets) 
        #delta4, index4 = deltaR(wp[1], genjets)
        if delta1 >0.4 or delta2 > 0.4:
            return hist
    except:
        return hist

    if index1 == index2:# or index3==index4:
        return hist
    else:        
        hist.Fill((genjets[index1].p4() + genjets[index2].p4()).M())
    return(hist)
    

path = "/pnfs/desy.de/cms/tier2/store/user/adewit/SummerProjectSamples/"

files = ["WminusH_HToBB_WToQQ_M125.root", "WplusH_HToBB_WToQQ_M125.root"]

gen_level_higgs = ROOT.TH1F("tmp", "tmp", 30, 50, 200)

global particleid

for f in files:
    infile = ROOT.TFile.Open(path+f)
    data = infile.Get("Events")
    inTree = InputTree(data)
    #fevent = array.array( 'f',  [0.0])
    #data.SetBranchAddress("fEvent",&fEvent);
    if "Wminus" in f:
        particleid = -24
    elif "Wplus" in f:
        particleid = 24
        

    nentries = data.GetEntries()
    entries = inTree.entries
    for ie,i in enumerate(xrange(entries)):
        e = Event(inTree, i)
        gen_level_higgs = analyze(e, gen_level_higgs)
        
        if i%10000 == 0:
            print("{} of {}".format(i, nentries))
    
#os.chdir()
saveplot(gen_level_higgs, "m_{H}")

        
