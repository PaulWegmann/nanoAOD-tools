import ROOT
import sys
import os
import array
import numpy as np
from operator import itemgetter
sys.path.append('../plotting')

from plots import saveratioplot, createTH1F, calc_AMS, savegraph, saveplot

def get_variables_2():
    return(["recon_Wmass", "recon_Hmass", "deltaR_reconHW","recon_Hpt","recon_Wpt","deltaR_bb", "deltaR_W", 
            "Jet_HT","bJet_pt1","bJet_pt2","WJet_pt1","WJet_pt2","bJet_btagDeepB1","bJet_btagDeepB2","WJet_btagDeepB1","WJet_btagDeepB2"])



#NTrees=300:MinNodeSize=7.5%:BoostType=AdaBoost:SeparationType=MisClassificationError:VarTransform=G:nCuts=400:MaxDepth=5 #old one
#VarTransform=U:BoostType=AdaBoost:SeparationType=MisClassificationError:MinNodeSize=1%:NTrees=400:MaxDepth=3 #new one which is more shit
def train( trigger, 
options = "VarTransform=U:NTrees=200:MaxDepth=3:MinNodeSize=5%",
#options = "BoostType=AdaBoost:VarTransform=U:SeparationType=CrossEntropy:SeparationType=CrossEntropy:MinNodeSize=1%:NTrees=50:MaxDepth=5:AdaBoostBeta=0.4",
signal = "data/breg/nores/signal.root",bckgrd = "data/breg/nores/bckgrd.root", TMVAoutfile = 'tmvaBDT', variables= ["recon_Wmass", 
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
            "WJet_mass1", #*curresults
            "WJet_mass2", #*curresults
            "WJet_neEmEF1", #*curresults
            "WJet_neEmEF2", #*curresults
            "WJet_neHEF1", #*curresults
            "WJet_neHEF2", #*curresults
            "WJet_qgl1", 
            "WJet_qgl2",
            "HJet_qgl1", #*curresults
            "HJet_qgl2", #*curresults
            "WJet_nConstituents1", 
            "WJet_nConstituents2"

            #"HJet_btagDeepC1", "HJet_btagDeepC2",
            #["bJet_"]
            ] #curcount = ~22
): #,"bJet_btagDeepB1","bJet_btagDeepB2","WJet_btagDeepB1","WJet_btagDeepB2"]): #! ncuts can set to -1 for ?better? !longer! performance
    

    signal_file = ROOT.TFile(signal, 'read')
    bckgrd_file = ROOT.TFile(bckgrd, 'read')
    output_file = ROOT.TFile(TMVAoutfile+trigger+".root", 'recreate')
    
    tmva_factory = ROOT.TMVA.Factory("BDT_"+trigger, output_file, "Color:DrawProgressBar:AnalysisType=Classification")
    tmva_dataloader = ROOT.TMVA.DataLoader("dataset")


    for variable in variables:
        tmva_dataloader.AddVariable(variable, "F")

    tmva_dataloader.AddSignalTree(signal_file.Get("Events" ), 1.)
    tmva_dataloader.AddBackgroundTree(bckgrd_file.Get("Events"), 1.)
    tmva_dataloader.SetWeightExpression("weight")
    
    print("*==* training with {} signal and {} background entries".format(
        signal_file.Get('Events').GetEntries(),
        bckgrd_file.Get('Events').GetEntries()))


    tmva_dataloader.PrepareTrainingAndTestTree(ROOT.TCut(trigger+">0"), ROOT.TCut(trigger+">0"), "SplitMode=Random:MixMode=Random:NormMode=None:SplitSeed=0")
    global_options = 'CreateMVAPdfs:NbinsMVAPdf=100:'

    #vi = ROOT.TMVA.VariableImportance(tmva_dataloader)
        
    #tmva_factory.BookMethod(tmva_dataloader, ROOT.TMVA.Types.kLikelihood, "Likelihood", global_options +\
    #    "TransformOutput:PDFInterpol=Splintriggere2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate") 
    #vi.BookMethod( ROOT.TMVA.Types.kBDT, "BDT" , global_options+ options)
    tmva_factory.BookMethod(tmva_dataloader, ROOT.TMVA.Types.kBDT, "BDT", global_options+options)    
    #tmva_factory.BookMethod(tmva_dataloader, ROOT.TMVA.Types.kBDT, "BDT", global_options +\
    #    "NTrees=1000:VarTransform=G:nCuts=250:MinNodeSize=3%:MaxDepth=5:BoostType=RealAdaBoost:SeparationType=GiniIndexWithLaplace")

    #print(dir(vi))
    tmva_factory.TrainAllMethods() # trains classifiers
    tmva_factory.TestAllMethods() # tests classifiers on test data and writes output into TMVAout file
    tmva_factory.EvaluateAllMethods() # evaluates
    #print(dir(tmva_factory))
    #importance = tmva_factory.EvaluateImportance(tmva_dataloader, 1, ROOT.TMVA.Types.kBDT, "BDT", "NTrees=200:MaxDepth=3")
    #saveplot(importance, "var")

    #results = vi.GetResults()
    #iv = results.GetImportanceValues()
    #for i in variables:
    #    print(i, iv.GetValue(i))
    output_file.Close()

def get_SignalFraction(signal_file, bckgrd_file):
    # Determine signal fraction in training sample from event weights
    label_value = array.array( 'b',  [ord('x')])
    weight_value = array.array( 'f',  [0.0])
    signal_tree = signal_file.Get("Events")
    background_tree = bckgrd_file.Get("Events")
    signal_tree.SetBranchAddress('weight',  weight_value)
    background_tree.SetBranchAddress('weight',  weight_value)

    signal_sum = 0.0
    for ievt in range(signal_tree.GetEntries()):
        signal_tree.GetEntry(ievt)
        signal_sum += weight_value[0]
    background_sum = 0.0
    for ievt in range(background_tree.GetEntries()):
        background_tree.GetEntry(ievt)
        background_sum += weight_value[0]
    # resulting signal fraction
    return signal_sum / (signal_sum + background_sum) 

# columns in .csv file written by function evaluate(), used in function analyse() 
TRUTH, WEIGHT, BDT, \
               BDT_prb,\
               BDT_rar =\
    0, 1, 2, 3, 4
nResult=5

def evaluate(trigger,signal = "data/breg/nores/signal.root",bckgrd = "data/breg/nores/bckgrd.root",variables= ["recon_Wmass", 
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
            "WJet_mass1", #*curresults
            "WJet_mass2", #*curresults
            "WJet_neEmEF1", #*curresults
            "WJet_neEmEF2", #*curresults
            "WJet_neHEF1", #*curresults
            "WJet_neHEF2", #*curresults
            "WJet_qgl1", 
            "WJet_qgl2",
            "HJet_qgl1", #*curresults
            "HJet_qgl2", #*curresults
            "WJet_nConstituents1", 
            "WJet_nConstituents2"

            #"HJet_btagDeepC1", "HJet_btagDeepC2",
            #["bJet_"]
            ]):#,"bJet_btagDeepB1","bJet_btagDeepB2","WJet_btagDeepB1","WJet_btagDeepB2"]):
#
# evaluate the trained MVA methods on the validation sample
#  input: 
#    feature_filename : root file with input data;  
#  output: 
#    results_filename:  important results in csv format 
# 
    file = ROOT.TFile("tmvaBDT"+trigger+".root", 'read')
    tree = file.Get("dataset/TestTree")

    max = 0
    min = 0


    result_filename = "result_"+trigger+".csv"
    
    #feature_tree = signal_tree.Add(bckgrd_tree) #feature_file.Get("validation")
    
    tmva_reader = ROOT.TMVA.Reader("!Color:!Silent")
    
    # Add the variables to TMVA reader
    feature_values = []
    for variable in variables:
        feature_values.append(array.array( 'f',  [0.]))
        tree.SetBranchAddress(variable,  feature_values[-1])
        tmva_reader.AddVariable(variable, feature_values[-1])
     
    label_value = array.array( 'b',  [ord('x')])
    
    weight_value = array.array( 'f',  [0.0])
    # Add label and weight
    tree.SetBranchAddress('classID',  label_value)
    
    tree.SetBranchAddress('weight',  weight_value)
    
    # tmva_reader.BookMVA("Likelihood", "dataset/weights/HiggsClassification_Likelihood.weights.xml")
    tmva_reader.BookMVA("BDT", "dataset/weights/BDT_"+trigger+"_BDT.weights.xml")

# some statistics over events read from validation sample


    result = np.zeros((tree.GetEntries(), nResult))

    # save results
    for ievt in range(tree.GetEntries()):
        tree.GetEntry(ievt)
        result[ievt, TRUTH] = int(label_value[0] == 0)
        result[ievt, WEIGHT] = weight_value[0]
        # save classification results ...
        #   1. as raw classifier:
        result[ievt, BDT]    = tmva_reader.EvaluateMVA('BDT') 
        #   2. transformed as "rarity" (i.e. background distribution flat in [0, 1.]
        result[ievt, BDT_rar]    = tmva_reader.GetRarity('BDT') 
        #   3. as signal probability (using mean signal fraction) 
        #result[ievt, BDT_prb]    = tmva_reader.GetProba('BDT', mean_fsig) 
    np.savetxt(result_filename, result, delimiter=',')
    print("*==* ... stored in file {}".format(result_filename))


def evaluate_all(trigger,signal = "data/breg/nores/signal.root",bckgrd = "data/breg/nores/bckgrd.root", variables = ["recon_Wmass", "recon_Hmass", "deltaR_reconHW","recon_Hpt","recon_Wpt","deltaR_bb", "deltaR_W", 
            "Jet_HT","bJet_pt1","bJet_pt2","WJet_pt1","WJet_pt2"]):
#
# evaluate the trained MVA methods on the validation sample
#  input: 
#    feature_filename : root file with input data;  
#  output: 
#    results_filename:  important results in csv format 
# 
    # Open ROOT files and create TMVA reader

    signal_file = ROOT.TFile(signal, 'read')
    bckgrd_file = ROOT.TFile(bckgrd, 'read')
    signal_tree = signal_file.Get("Events")
    bckgrd_tree = bckgrd_file.Get("Events")

    result_filename = "result_"+trigger+".csv"
    
    #feature_tree = signal_tree.Add(bckgrd_tree) #feature_file.Get("validation")
    
    tmva_reader = ROOT.TMVA.Reader("!Color:!Silent")
    
    # Add the variables to TMVA reader
    feature_values = []
    for variable in variables:
        feature_values.append(array.array( 'f',  [0.]))
        signal_tree.SetBranchAddress(variable,  feature_values[-1])
        bckgrd_tree.SetBranchAddress(variable,  feature_values[-1])
        tmva_reader.AddVariable(variable, feature_values[-1])
     
    #label_value = array.array( 'b',  [ord('x')])
    
    weight_values = array.array( 'f',  [0.0])
    weight_valueb = array.array( 'f' , [0.0])
    # Add label and weight
    #feature_tree.SetBranchAddress('Label',  label_value) #so label_value gets overwritten?
    
    signal_tree.SetBranchAddress('weight',  weight_values)
    bckgrd_tree.SetBranchAddress('weight',  weight_valueb)
    
    # tmva_reader.BookMVA("Likelihood", "dataset/weights/HiggsClassification_Likelihood.weights.xml")
    tmva_reader.BookMVA("BDT", "dataset/weights/BDT_"+trigger+"_BDT.weights.xml")

# some statistics over events read from validation sample
    signal_sum = 0.0
    background_sum = 0.0
    ns=0
    nb=0

    """
    for ievt in range(signal_tree.GetEntries()):
        signal_tree.GetEntry(ievt)
        signal_sum += weight_values[0]
        ns += 1
    for ievt in range(bckgrd_tree.GetEntries()):     
        bckgrd_tree.GetEntry(ievt)
        background_sum += weight_valueb[0]
        nb += 1
    

    print("*==* evaluating validation sample")
    print("{} signal and {} background entries read from validation sample".format(ns,nb))
    print("corresponding to  {} signal and {} background events".format(signal_sum,background_sum))

    # Collect results and write to numpy array
    
    mean_fsig=get_SignalFraction(signal_file, bckgrd_file)
    print("*==* evaluating MVA results ...")
    print("signal fraction in training sample is {}".format(mean_fsig))
    """
    result = np.zeros((signal_tree.GetEntries()+bckgrd_tree.GetEntries(), nResult))

    # save results
    for ievt in range(signal_tree.GetEntries()):
        signal_tree.GetEntry(ievt)
        result[ievt, TRUTH] = 1
        result[ievt, WEIGHT] = weight_values[0]
        # save classification results ...
        #   1. as raw classifier:

        #result[ievt, Likelihood] = tmva_reader.EvaluateMVA('Likelihood')
        result[ievt, BDT]    = tmva_reader.EvaluateMVA('BDT')
        #print("yas")
        #   2. transformed as "rarity" (i.e. background distribution flat in [0, 1.]
        #result[ievt, Likelihood_rar] = tmva_reader.GetRarity('Likelihood')
        result[ievt, BDT_rar]    = tmva_reader.GetRarity('BDT') 
        #   3. as signal probability (using mean signal fraction) 
        #result[ievt, Likelihood_prb] = tmva_reader.GetProba('Likelihood', mean_fsig)
        #result[ievt, BDT_prb]    = tmva_reader.GetProba('BDT', mean_fsig) 
    
    for ievt in np.arange(bckgrd_tree.GetEntries())+signal_tree.GetEntries():
        bckgrd_tree.GetEntry(ievt-signal_tree.GetEntries())
        result[ievt, TRUTH] = 0
        result[ievt, WEIGHT] = weight_valueb[0]
        # save classification results ...
        #   1. as raw classifier:
        #result[ievt, Likelihood] = tmva_reader.EvaluateMVA('Likelihood')
        result[ievt, BDT]    = tmva_reader.EvaluateMVA('BDT') 
        #   2. transformed as "rarity" (i.e. background distribution flat in [0, 1.]
        #result[ievt, Likelihood_rar] = tmva_reader.GetRarity('Likelihood')
        result[ievt, BDT_rar]    = tmva_reader.GetRarity('BDT') 
        #   3. as signal probability (using mean signal fraction) 
        #result[ievt, Likelihood_prb] = tmva_reader.GetProba('Likelihood', mean_fsig)
        #result[ievt, BDT_prb]    = tmva_reader.GetProba('BDT', mean_fsig) 
    
    np.savetxt(result_filename, result, delimiter=',')
    print("*==* ... stored in file {}".format(result_filename))

def AsimovSignificance(s,b):
    #print(2*(s+b)*np.log(s/b)-s)
    return np.sqrt(2*((s+b)*np.log(1+s/b)-s))


def analyse(trigger, options="", return_range=False, create_plots=False):
    # Load results
    result_filename = "result_"+trigger+".csv"



    result = np.loadtxt(result_filename, delimiter=',')
    
    # Calculate Mean of Average Distance (MAD)
    mad = lambda x, y, w: np.sum(w*np.abs(x - y))/np.sum(w)
    print("Mean of Average Distance (MAD):")
    mad1 =  mad(result[:,TRUTH], result[:, BDT], result[:, WEIGHT])
    print("  MAD for BDT", mad1)
    
    # Calculate Mean Squared Error, MSE
    mse = lambda x, y, w: np.sum(w*(x-y)**2)/np.sum(w)
    print("Mean Squared Error (MSE):")
    mse1 = mse(result[:,TRUTH], result[:, BDT], result[:, WEIGHT])
    print("  MSE for BDT", mse1)
    
    def ams(x, y, w, cut):
    # Calculate Average Mean Significane as defined in ATLAS paper
    #    -  approximative formula for large statistics with regularisation
    # x: array of truth values (1 if signal)
    # y: array of classifier result
    # w: array of event weights
    # cut
        t = y > cut 
        s = np.sum((x[t] == 1)*w[t])
        b = np.sum((x[t] == 0)*w[t])
        return s/np.sqrt(b+3.0)

    def find_best_ams(x, y, w):
    # find best value of AMS by scanning cut values; 
    # x: array of truth values (1 if signal)
    # y: array of classifier results
    # w: array of event weights
    #  returns 
    #   ntuple of best value of AMS and the corresponding cut value
    #   list with corresponding pairs (ams, cut) 
    # ----------------------------------------------------------
        ymin=min(y) # classifiers may not be in range [0.,1.]
        ymax=max(y)
        nprobe=200    # number of (equally spaced) scan points to probe classifier 
        amsvec= [(ams(x, y, w, cut), cut) for cut in np.linspace(ymin, ymax, nprobe)] 
        maxams=sorted(amsvec, key=lambda lst: lst[0] )[-1]
        return maxams, amsvec

    maxams_BDT, amsvec_BDT=find_best_ams(result[:,TRUTH], result[:, BDT], result[:, WEIGHT])

    print("Average Mean Sensitivity (AMS) and cut value:")
    print("  - raw classifier")
    print("  AMS for BDT", maxams_BDT)

    maximum, minimum = max(result[:,BDT]), min(result[:,BDT])
    
    
    signal = ROOT.TH1F("signal", "signal", 20, minimum, maximum) # * BDT value is not always between 0 and 1
    bckgrd = ROOT.TH1F("bckgrd", "bckgrd", 20, minimum, maximum)
    for index in range(len(result)):
        if result[index,TRUTH] == 1:
            signal.Fill(result[index,BDT], result[index, WEIGHT])
        else:
            bckgrd.Fill(result[index,BDT], result[index, WEIGHT])
    
    print(bckgrd.Integral())

    AMSval = calc_AMS(signal, bckgrd)
    print("AMS calculated with {} signal and {} background events accounts to {}".format(signal.Integral(), bckgrd.Integral(), AMSval))

    #* calculating final AMS on recon_Hmass for comparison
    #* 
    cuts = np.linspace(maximum, minimum, 15)
    AMS = []
    eventsleft = []
    if True:
        file = ROOT.TFile("tmvaBDT"+trigger+".root", 'read')
    
        test = file.Get("dataset/TestTree")
        train = file.Get("dataset/TrainTree")

        #* find maximum AMS value
        for i in range(len(cuts)):
            #two trees needs to get merched, test and signal
            signal1, _ = createTH1F(test, "recon_Hmass", np.linspace(90, 160, 16), "(BDT>{} &&classID==0)*weight".format(cuts[i]), dif = "test")
            bckgrd1, _ = createTH1F(test, "recon_Hmass", np.linspace(90, 160, 16), "(BDT>{} &&classID>0)*weight".format(cuts[i]), dif = "train")

            #tmpsignal, _ = createTH1F(train, "recon_Hmass", np.linspace(90, 155, 16), "(BDT>{} &&classID==0)*weight".format(cuts[i]), dif = "test1")
            #tmpbckgrd, _ = createTH1F(train, "recon_Hmass", np.linspace(90, 155, 16), "(BDT>{} &&classID>0)*weight".format(cuts[i]), dif = "train1")

            # merge
            #signal1.Add(tmpsignal)
            #bckgrd1.Add(tmpbckgrd)
            #del(tmpsignal, tmpbckgrd)
            
            #if calc_AMS(signal, bckgrd)>AMS:
            AMS.append(calc_AMS(signal1, bckgrd1))            
            eventsleft.append([signal1.Integral(), bckgrd1.Integral()])
            del(signal1, bckgrd1)

        value = 0
        while(value == 0):
            index, helper =  max(enumerate(AMS), key=itemgetter(1))
            if eventsleft[index][0] < 150:
                del(AMS[index], eventsleft[index])
            else:
                value = helper
                
        
        if create_plots:
            savegraph(cuts, AMS, "AMS_graph_{}".format(trigger), "BDT>", "AMS", heading=trigger,  text =  "Max {} at {}".format(value, round(cuts[index],2)))

            
            #* make plot with maximum AMS
            signal1, _ = createTH1F(test, "recon_Hmass", np.linspace(90, 160, 16), "(BDT>{} &&classID==0)*weight".format(cuts[index]), dif = "test")
            bckgrd1, _ = createTH1F(test, "recon_Hmass", np.linspace(90, 160, 16), "(BDT>{} &&classID>0)*weight".format(cuts[index]), dif = "train")

            tmpsignal, _ = createTH1F(train, "recon_Hmass", np.linspace(90, 160, 16), "(BDT>{} &&classID==0)*weight".format(cuts[index]), dif = "test1")
            tmpbckgrd, _ = createTH1F(train, "recon_Hmass", np.linspace(90, 160, 16), "(BDT>{} &&classID>0)*weight".format(cuts[index]), dif = "train1")
            signal1.Add(tmpsignal)
            bckgrd1.Add(tmpbckgrd)
            signal1.Scale(1/signal1.Integral())
            bckgrd1.Scale(1/bckgrd1.Integral())

            del(tmpsignal, tmpbckgrd)

            saveratioplot(signal1, bckgrd1, "recon_Hmass", "signal", "background", calcArea=False, text = "AMS {}".format(AMS[index]), dif = trigger, title = trigger)

            signal1.Scale(1)
            bckgrd1.Scale(1)
            signal1.Add(bckgrd1)
            saveplot(signal1, "recon_Hmass", "data", trigger)
            del(signal1, bckgrd1)



    #* getting that array -------------------------------------
    #* (debugging), and now for calculating the ams value without changing the image
    #* not necessary anymore
    #* it was used as a check for the AMS value, but now it works
    

    signalarray = []
    bckgrdarray = []
    for i in range(signal.GetNbinsX()):
        signalarray.append(signal.GetBinContent(i+1))
        bckgrdarray.append(bckgrd.GetBinContent(i+1))





    signalarray = np.array(signalarray)
    bckgrdarray = np.array(bckgrdarray)

    def AsimovSignificance(s,b):
        return np.sqrt(2*((s+b)*np.log(1+s/b)-s))

    asmvdif = []
    rest = 0

    signalarray = signalarray[::-1]
    bckgrdarray = bckgrdarray[::-1]

    # signalarray = np.array([ 17.210014,  61.499397,  87.940155, 120.51441,  143.79169,  169.26929,
    # 203.37524,  230.97742,  271.31696,  308.65457,  342.6233,   340.9767,
    # 320.40305,  297.67145,  285.51694,  260.3653,   200.59111,  152.03879,
    # 95.85461,   38.200207])
    # bckgrdarray = np.array([1.1238852e+04, 3.5990004e+04, 1.3388028e+05, 4.9230538e+05, 1.0183296e+06,
    # 1.5634670e+06, 2.5431190e+06, 4.1391840e+06, 6.5147880e+06, 8.5776930e+06,
    # 1.2227979e+07, 1.4140604e+07, 1.5295926e+07, 1.6685774e+07, 1.8673644e+07,
    # 2.1107508e+07, 2.2592300e+07, 2.2869652e+07, 2.1650750e+07, 1.1851990e+07])

    # for i in range(len(signalarray)):
    #     signal.SetBinContent(i+1, signalarray[i])
    #     bckgrd.SetBinContent(i+1, bckgrdarray[i])

    # AMSval = calc_AMS(signal, bckgrd)
    # #signal.scale(1/signal)


    for sig, bck in zip(signalarray, bckgrdarray):
        if bck!=0:
            asmvdif.append(AsimovSignificance(sig+rest, bck))
            rest = 0
        else:
            rest+=sig

    asmvdif = np.array(asmvdif)
    amsvalue = np.sqrt(np.sum(asmvdif**2))

    print(amsvalue, AMSval , "!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print("AMS recon : {}".format(value))
    
    #* in princible unnecessary


    signal.Scale(1/signal.Integral())
    bckgrd.Scale(1/bckgrd.Integral())

    
    if create_plots:
        saveratioplot(signal, bckgrd, "BDT", "signal", "bckgrd", calcArea=False, title=trigger, text = "weighted AMS "+str(round(amsvalue,3)), dif = trigger+ options)
    if return_range:
        return(mad1, mse1, AMSval, value, eventsleft[index][0], eventsleft[index][1], minimum, maximum )
    else:
        return(mad1, mse1, AMSval, value, eventsleft[index][0], eventsleft[index][1])
    

    # here, one could test whether using transformed classifiers makes a difference
    #print " - classifier transformed to rarity"
    #     ...
    #print " - classifier expressed as signal probability"
    #     ...
# 
# some plots ...
    # show performance score as a funtion of the cut value
    #    remark: could be matplotlib graphs as well
    # some options first:


    """
    ROOT.gStyle.SetLineColor(38)     
    ROOT.gStyle.SetLineWidth(2) 

    c1=ROOT.TCanvas("Classifier Performance")

    amsvec_BDT=np.asarray(amsvec_BDT) 
    xarr=np.array(amsvec_BDT[:,1])
    yarr=np.array(amsvec_BDT[:,0])
    g_amsBDT= ROOT.TGraph(len(xarr), xarr, yarr)
    g_amsBDT.SetTitle("Significance (AMS) vs. cut on classifier") 
    g_amsBDT.GetXaxis().SetTitle("Cut on BDT classifier")
    g_amsBDT.GetYaxis().SetTitle("Significance (ams)")
    g_amsBDT.Draw("ACLP")

    c1.Update()    

    """


if __name__=="__main__":
    if len(sys.argv) != 2 or (len(sys.argv) == 2 and sys.argv[1] not in ['train', 'evaluate', 'analyse', 'TMVAgui', 'evaluate_all']) :
        print("*** Usage: python {} [train|evaluate|TMVAgui|analyse]".format(sys.argv[0]))
        sys.exit(1)

    if ROOT.gROOT.GetVersionCode() >= 332288 and ROOT.gROOT.GetVersionCode() < 332544:
        print("*** You are running ROOT version 5.18, which has problems in PyROOT such that TMVA")
        print("*** does not run properly (function calls with enums in the argument are ignored).")
        print("*** Solution: either use CINT or a C++ compiled version (see TMVA/macros or TMVA/examples),")
        print("*** or use another ROOT version (e.g., ROOT 5.19).")
        sys.exit(1)

#    ROOT.gROOT.SetBatch() # switch on root batch mode
#                      useful if plots are to ge generated without displaying them

    infile='atlas-higgs-challenge-2014-v2_part.root'
    TMVAoutfile='tmvaBDTHLT_PFHT300PT30_QuadPFJet_75_60_45_40.root'
    resultfile='result.csv'
    triggers = ["HLT_PFHT300PT30_QuadPFJet_75_60_45_40", "HLT_PFHT180"]

    if sys.argv[1] == 'train':
        for i in triggers:
            train(trigger=i)
        #train(triggers[0])
    elif sys.argv[1] == 'evaluate':
        for i in triggers:
            evaluate(i)
    elif sys.argv[1] == 'analyse':
        for i in triggers:
            _, _,_, ams  ,_ ,_ = analyse(i, create_plots=False)
            print("The calculated AMS value accounts to {}".format(ams))
    elif sys.argv[1] == 'evaluate_all':
        for i in triggers:
            evaluate_all(i)
    elif sys.argv[1] == 'TMVAgui':
        ROOT.TMVA.TMVAGui(TMVAoutfile)
        input('Press <ret> to contunue -->')

        #result[ievt, Likelihood_rar] = tmva_reader.GetRarity('Likelihood')