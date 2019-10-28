import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import sys
import numpy as np

from operator import itemgetter

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import * 

elements = 10

def deltaR(b, jets):
    delta = []
    for jet in jets:
        delta.append(np.sqrt((b.eta-jet.eta)**2+(b.phi-jet.phi)**2))
        
    val, idx = min((val, idx) for (idx, val) in enumerate(delta))
    return(val, idx)

def second_smallest(numbers):
    m1, m2 = float('inf'), float('inf')
    ind1, ind2 = 0,0
    for index, x in enumerate(numbers):
        if x <= m1:
            m1, m2 = x, m1
            ind1, ind2 = index, ind1
        elif x < m2:
            m2 = x
            ind2 = index
    return ind2, m2



class BaselineSelection(Module):
    def __init__(self):

        pass
    def beginJob(self):
        self.numevents = 0
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        self.out.branch("original_Hmass", "F") #H_mass calculated from gen bquarks 
        self.out.branch("genjets_Hmass", "F") #H_mass calculated from genjets
        
        self.out.branch("asso_bjet_H_id1", "I") # index of genjet associated with genpart b
        self.out.branch("asso_bjet_H_id2", "I") # --
        self.out.branch("deltaR_bjet1", "F") # deltar between genpart and genjet b
        self.out.branch("deltaR_bjet2", "F") # --
        
        
        self.out.branch("asso_jet_W_id1", "I") # same as before but for W+ boson and no specifc genpart quark
        self.out.branch("asso_jet_W_id2", "I") # it is hadronic tough
        self.out.branch("deltaR_W_jet1", "F") # --
        self.out.branch("deltaR_W_jet2", "F")# --
        
        self.out.branch("wnrightjets", "I") # rightly assigned reconJets regarding W Boson
        self.out.branch("genjets_deltaR_W", "F") # deltaR between them
        self.out.branch("genjets_WMass", "F")
        self.out.branch("recon_WMass", "F")
        self.out.branch("wnrightjets_bad", "I") # using maximal pt jets
        self.out.branch("recon_WJets_id", "I", lenVar=str(2)) #original ids of Jets picked for w reconstruction

        


        self.out.branch("Jet_HT","F"); # Sum of transverse P of all Jets
        self.out.branch("nJets", "F"); # number of Jets in this event
        self.out.branch("fatjetsused", "I") # wheter or not a bjet fat jet was used for calculation of H_mass
        
        
        self.out.branch("nmatch", "I") # number of candidates, min 1
        

        self.out.branch("deltaR_bb", "F", lenVar=str(elements)); # deltaR of reconstructed btagged Jets used for H mass calculation 
        self.out.branch("recon_Hmass", "F", lenVar=str(elements)); #H_mass calculated from reconstructed jets (with some algorithm)
        self.out.branch("nrightjets", "I", lenVar=str(elements)) # count of right connection between genjets and reconjets (max 2)
        self.out.branch("asso_bjet_H_id_recon1","I", lenVar=str(elements)) # same as before but on recon level
        self.out.branch("asso_bjet_H_id_recon2", "I", lenVar=str(elements)) # --
        
        
        #variable for selecting hmass if there are multiple bjets
        self.out.branch("best_candidate_id_1", "I") #min(deltaR) ,
        self.out.branch("candidate_score_1", "F")
        
        self.out.branch("best_candidate_id_2", "I") #hmass-110 
        self.out.branch("candidate_score_2", "F")#
        
        self.out.branch("best_hmass110", "I")
        self.out.branch("best_delta2", "I")
        self.out.branch("best_btagdeep", "I") #this turned out to be the best one --------------------
        self.out.branch("best_jetpt65", "I")
        self.out.branch("isright", "O", lenVar=str(2))
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    #def analyze(self, event):
    #    return True
    
    
    def analyze(self, event):
        # ---------------------------------------------------
        # just for easier debugging when needed
        prterrmsg = False
        #self.numevents +=1
        #if self.numevents %10 !=0:
            #return False

        
        # ------------------------------------------------
        # here the GenPart bquarks are tried to connect to GenJets (no selection is done) (added selection of particles originating from W+)
        genparticles = Collection(event, "GenPart")
        genjets = Collection(event, "GenJet")
        
        bs = []
        wp = [] #"w+ particles"
        for genpart in genparticles:
            #selection of b particles originating from Higgs
            if(abs(genpart.pdgId)==5 and genpart.genPartIdxMother>-1):
                if genparticles[genpart.genPartIdxMother].pdgId==25 and genparticles[genpart.genPartIdxMother].statusFlags&8192==8192: #statusflag is used too make sure it is the "last copy of the particle" due to the nature of the generation algorhythm which often copies particles multiple times to simulate decay
                    bs.append(genpart)
            #selection of b particles originating from W+
            if genpart.genPartIdxMother>-1:
                if genparticles[genpart.genPartIdxMother].pdgId==24 and genparticles[genpart.genPartIdxMother].statusFlags&8192==8192:
                    wp.append(genpart)
        
        if len(bs)!=2 or len(wp)!=2:
            if prterrmsg:
                print("not both particles have been found")
            return False
        #print(bs[0].mass, bs[1].mass)
        
        try:
            delta1, index1 = deltaR(bs[0], genjets) 
            delta2, index2 = deltaR(bs[1], genjets)
            delta3, index3 = deltaR(wp[0], genjets) 
            delta4, index4 = deltaR(wp[1], genjets)
            if index1 == index2:# or index3==index4:
                if prterrmsg:
                    print("it is assigned to the same jet")
                return False
            
            # to reconstruct H_mass proberly the Mass of the bottom quarks needs to be added
            
            v1 = ROOT.TLorentzVector()
            v1.SetPtEtaPhiM(bs[0].pt, bs[0].eta, bs[0].phi, 4.18)
            v2 = ROOT.TLorentzVector()
            v2.SetPtEtaPhiM(bs[1].pt, bs[1].eta, bs[1].phi, 4.18)
            self.out.fillBranch("original_Hmass", (v1+v2).M())
            
            self.out.fillBranch("genjets_Hmass", (genjets[index1].p4() + genjets[index2].p4()).M())
            #self.out.fillBranch("original_Hmass", (bs[0].p4()+bs[1].p4()).M())

            #print((bs[0].p4()+bs[1].p4()).M())
            self.out.fillBranch("asso_bjet_H_id1", int(index1))
            self.out.fillBranch("asso_bjet_H_id2", int(index2)) #important!!! these indices only refer to my implementation. They are not equal to the JetId which the particle is constructed with
            self.out.fillBranch("deltaR_bjet1", delta1)
            self.out.fillBranch("deltaR_bjet2", delta2)
            self.out.fillBranch("asso_jet_W_id1", int(index3))
            self.out.fillBranch("asso_jet_W_id2", int(index4))
            self.out.fillBranch("deltaR_W_jet1", delta3)
            self.out.fillBranch("deltaR_W_jet2", delta4)
            self.out.fillBranch("genjets_WMass", (genjets[index3].p4() + genjets[index4].p4()).M())
            self.out.fillBranch("genjets_deltaR_W", np.sqrt((genjets[index3].eta-genjets[index4].eta)**2+(genjets[index3].phi-genjets[index4].phi)**2))
        except:
            if prterrmsg:
                print("The length of genjets is {}".format(len(genjets)))
            return False


        # -------------------------------------------------
        # here the algorhythm which tries to reconstruct the h_mass out of the "detected Jets" is applied
        # description of current algorithm: currently I am saving all events within a certain range
        
        
        jets = list(Collection(event,"Jet"))
        fatJets = list(Collection(event,"FatJet"))
        
        #apply some selection:
        selBJets = [x for x in jets if x.pt>20 and abs(x.eta)<2.5 and x.jetId>0 and x.btagDeepB>0.4941]
        selJets = [x for x in jets if x.pt>20 and abs(x.eta)<4.8 and x.jetId>0]
        
        
        #selfatJets =  [x for x in fatJets if x.pt>20 and abs(x.eta)<2.5 and x.jetId>0]
        #selBfatJets = [x for x in fatJets if x.pt>40 and abs(x.eta)<2.5 and x.jetId>0 and x.btagDeepB>0.4941]
        
        if len(selBJets) < 2:# and len(selBfatJets)<1: #Require at least two b-tagged or at least one fatjet
            return False
        
        Jet_HT=0.
        reco_H = ROOT.TLorentzVector()
        
        
        #if event.MET_pt <50 and event.MET_pt > 150: # seletion for MET
        #    return False
            
        
        #trying to identify the right btagged jets ---------------------------
        
        delta = 0. # this step is important: variables need to be initialised with the proper type otherwise it might won't work and there wont be errors -> actually this is not true, thought so because of another error
        njets = 0.
        
        for j in selJets: #calculating 
            Jet_HT+=(j.pt)
            njets+=1
        #for j in selfatJets:
            #Jet_HT+=(j.pt)
            #njets+=1
        #
        #found = False
        
        
        # gathering all possible Hmasses
        hmasses = np.full(elements, -999, dtype = float)
        deltas = np.full(elements, -999, dtype = float)
        indexes1 = np.full(elements, -999, dtype = int)
        indexes2 = np.full(elements, -999, dtype = int)
        poss = 0

        
        for i1,curjet1 in enumerate(selBJets):
            if poss == elements-1:
                break
            for i2,curjet2 in enumerate(selBJets):
                if poss == elements-1:
                    break
                #p4 is a Lorentzvector so here I take two arbitrary b-jets and treat them as my bottom quarks 
                v1 = curjet1.p4()
                v2 = curjet2.p4()
                hmass = (v1+v2).M()
                delta = np.sqrt((curjet1.eta-curjet2.eta)**2+(curjet1.phi-curjet2.phi)**2)
                if 90 < hmass and hmass < 155 and i2>i1:
                    hmasses[poss] = hmass
                    deltas[poss] = delta
                    indexes1[poss] = i1
                    indexes2[poss] = i2
                    poss+=1
                    
        
        if poss==0:
            return False
        else:
            nrights = np.full(elements, -999, dtype = int)

            for i in range(len(indexes1)):
                if i == elements-1 or indexes1[i] == -999:
                    break
                nright = 0
                if selBJets[indexes1[i]].genJetIdx == index1 or selBJets[indexes1[i]].genJetIdx == index2:
                    nright+=1
                if selBJets[indexes2[i]].genJetIdx == index1 or selBJets[indexes2[i]].genJetIdx == index2:
                    nright+=1
                nrights[i]=nright
            
            
            # selecting favorite candidates: --------------------------------------------------
            ncriteria = np.zeros(poss)
            ncriteria2 = np.zeros(poss)
            if poss > 1:
                helper = np.array(np.array(selBJets)[indexes1[:poss]])
                helper1 = np.array(np.array(selBJets)[indexes2[:poss]])
                helper2 = []
                helper3 = []
                
                for j, i in zip(helper, helper1):
                    helper2.append(j.btagDeepB)
                    helper3.append(i.btagDeepB)
                
                helper2=np.array(helper2)
                helper3 = np.array(helper3)
                
                helper = helper2+helper3
                s = max(enumerate(helper), key=itemgetter(1))[0]
                self.out.fillBranch("best_btagdeep", s)
                
                """ # further analysis for picking favorite candidates when selection isn't unique
                helper = np.array(np.array(selBJets)[indexes1[:poss]])
                helper1 = np.array(np.array(selBJets)[indexes2[:poss]])
                helper2 = []
                helper3 = []
                
                s = min(enumerate(abs(hmasses-110)), key=itemgetter(1))[0]
                self.out.fillBranch("best_hmass110", s)
                ncriteria[s]+=1.6
                ncriteria2[s]+=2
                
                
                
                for j, i in zip(helper, helper1):
                    helper2.append(j.pt-65.)
                    helper3.append(i.pt-65.)
                
                helper2=np.array(helper2)
                helper3 = np.array(helper3)
                
                helper = helper2+helper3
                helper2 = abs(helper2)+abs(helper3)
                
                s = min(enumerate(helper), key=itemgetter(1))[0]
                ncriteria[s]+=1.9
                
                s = min(enumerate(helper2), key=itemgetter(1))[0]
                self.out.fillBranch("best_jetpt65", s)
                ncriteria2[s]+=1.9
                
                #min1 = min(helper)
                #trash, min2 = second_smallest(helper)
                
                #found = False
                
                #for i in range(len(indexes1[:poss])):
                    #if selBJets[indexes1[i]].pt == min1 and selBJets[indexes2[i]].pt == min2 or selBJets[indexes2[i]].pt == min1 and selBJets[indexes1[i]].pt == min2:
                        #ncriteria[i]+=3
                        #found=True
                    #break
                
                #self.noptmin = 0
                #if not(found):
                    #self.noptmin+=1
                    #return False
                
                s = min(enumerate(np.array(deltas[:poss])), key=itemgetter(1))[0]
                ncriteria[s]+=2
                s = min(enumerate(abs(np.array(deltas[:poss])-2)), key=itemgetter(1))[0]
                self.out.fillBranch("best_delta2", s)
                ncriteria2[s]+=2
                
                
                helper = np.array(np.array(selBJets)[indexes1[:poss]])
                helper1 = np.array(np.array(selBJets)[indexes2[:poss]])
                helper2 = []
                helper3 = []
                
                for j, i in zip(helper, helper1):
                    helper2.append(j.btagDeepB)
                    helper3.append(i.btagDeepB)
                
                helper2=np.array(helper2)
                helper3 = np.array(helper3)
                
                helper = helper2+helper3
                

                ncriteria[s]+=1.7
                ncriteria2[s]+=1.7
                
                indextmp = max(enumerate(ncriteria), key=itemgetter(1))[0]
                
                indextmp2 = max(enumerate(ncriteria2), key=itemgetter(1))[0]
                
                
                
                #indextmp = min(enumerate(np.array(deltas[:poss-1])-2), key=itemgetter(1))[0]
                #ncriteria[indextmp]+=2
                #helper = []
                #for j in selBJets[indexes1[:poss-1]]:
                    #helper.append(j.btagDeepB)
                    
                #indextmp = max(enumerate(np.array(helper)), key=itemgetter(1))[0]
                #ncriteria[indextmp]+=1.5
                
                #indextmp = max(enumerate(ncriteria), key=itemgetter(1))[0]
            
            
                    
                self.out.fillBranch("best_candidate_id_1", int(indextmp))
                self.out.fillBranch("candidate_score_1", max(ncriteria))
                
                self.out.fillBranch("best_candidate_id_2", int(indextmp2))
                self.out.fillBranch("candidate_score_2", max(ncriteria2))
                """
            else:
                s=0
                self.out.fillBranch("best_btagdeep", s)
                """ # variables used for analysis of best selection criteria
                self.out.fillBranch("best_candidate_id_1", int(0))
                self.out.fillBranch("candidate_score_1", -999)
                self.out.fillBranch("best_candidate_id_2", int(0))
                self.out.fillBranch("candidate_score_2", -999)
                self.out.fillBranch("best_delta2", 0)
                self.out.fillBranch("best_hmass110", 0)
                self.out.fillBranch("best_jetpt65", 0)
                """
            
            
            # trying to connect Jets with W+ Boson -----------------------------------------------
            
            tmp = []
            for i in selJets:
                if i == selBJets[indexes1[s]] or i == selBJets[indexes2[s]]:
                    continue
                else:
                    tmp.append(i)
            
            tmp2 = [i.pt for i in tmp]
            
            if len(tmp)<2: #skipping events with too few Jets
                return False
            
            # gathering biggest pts
            
            tmp2 = [i.pt for i in tmp]
            windex1,value1 =  min(enumerate(tmp2), key=itemgetter(1))
            windex2, value2 = second_smallest(tmp2)
            
            wnright = 0
            if tmp[windex1].genJetIdx == index3 or tmp[windex1].genJetIdx == index4:
                wnright+=1
            if tmp[windex2].genJetIdx == index3 or tmp[windex2].genJetIdx == index4:
                wnright+=1
            
            self.out.fillBranch("wnrightjets_bad", wnright)
            
            
            # figuring out which one is closest to WMass (~80GeV)
            Wmasses = []
            Wind1 = []
            Wind2 = []
            for j, elj in enumerate(tmp):
                for i, eli in enumerate(tmp):
                    if i>j:
                        Wmass = (elj.p4()+eli.p4()).M()
                        Wmasses.append(Wmass)
                        Wind1.append(j)
                        Wind2.append(i)
                                    
            Wind = min(enumerate(abs(np.array(Wmasses)-80)), key=itemgetter(1))[0]
            
            
            wnright = 0
            isrights = [False, False]
            if tmp[Wind1[Wind]].genJetIdx == index3 or tmp[Wind1[Wind]].genJetIdx == index4:
                wnright+=1
                isrights[0]=True
            if tmp[Wind2[Wind]].genJetIdx == index3 or tmp[Wind2[Wind]].genJetIdx == index4:
                wnright+=1
                isrights[1]=True
            winds = np.array([jets.index(tmp[Wind1[Wind]]), jets.index(tmp[Wind2[Wind]])], dtype=int)
            self.out.fillBranch("wnrightjets", wnright)
            self.out.fillBranch("recon_WMass", Wmasses[Wind])
            self.out.fillBranch("recon_WJets_id", winds)
            self.out.fillBranch("isright", isrights)
            
            # ---------------------- normal filling
            self.out.fillBranch("deltaR_bb", deltas)
            self.out.fillBranch("nmatch", poss)
            self.out.fillBranch("nrightjets", nrights)            
            self.out.fillBranch("fatjetsused", 0)
            self.out.fillBranch("recon_Hmass", hmasses)   
            
            
            
            #transforming indexes to orignal ones
            indexes1 =  [jets.index(selBJets[i]) for i in indexes1[:poss]]
            indexes2 =  [jets.index(selBJets[i]) for i in indexes2[:poss]]
            self.out.fillBranch("asso_bjet_H_id_recon1", indexes1)
            self.out.fillBranch("asso_bjet_H_id_recon2", indexes2)
        
        self.out.fillBranch("Jet_HT",Jet_HT) # what is this? sum of all transverse momenta
        self.out.fillBranch("nJets", njets)
        

        return True

basesel = lambda: BaselineSelection()
