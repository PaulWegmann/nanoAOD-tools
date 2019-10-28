from __future__ import division
import ROOT
import PhysicsTools.NanoAODTools.plotting as plot
import os
import numpy as np
from array import array

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
#ROOT.TH1.AddDirectory(False)
plot.ModTDRStyle()

def createAxisHists(n,src,xmin=0,xmax=499):
  result = []
  for i in range(0,n):
    res = src.Clone()
    res.Reset()
    res.SetTitle("")
    res.SetName("axis%(i)d"%vars())
    res.SetAxisRange(xmin,xmax)
    res.SetStats(0)
    result.append(res)
  return result

def createTH1F(tree, varname,  bins, constraints="", dif = "", normalise = False):
    #! returns th1f files as well as number of entries before and after applying the suggested entries
    
    #* get entries only using elements with reasonable value
    tmp1 = ROOT.TH1F("tmp", "tmp",len(bins)-1,bins)
    tmp1.Sumw2()
    tree.Draw(varname+">>tmp")#, varname+"!=-999")
    before = int(tmp1.GetEntries())
    del(tmp1)
    
    #* normal creating TH1F
    tmp = ROOT.TH1F(varname.lower()+"_"+''.join(e for e in constraints if e.isalnum())+dif, varname.lower()+"_"+''.join(e for e in constraints if e.isalnum())+dif,len(bins)-1,bins)
    tmp.Sumw2()
    after = tmp.GetEntries()
    tree.Draw(varname+">>"+varname.lower()+"_"+''.join(e for e in constraints if e.isalnum())+dif, constraints)
    after = int(tmp.GetEntries())
    if normalise:
        tmp.Scale(1/tmp.Integral())
    return(tmp,[before, after])

# ''.join(e for e in constraints if e.isalnum()) is used to remove special characters out of string, otherwise it doesnt work porperly

def createTH2F(tree, varnamex, varnamey,  binsx, binsy, dif = ""): #currently kind of hard coded to work for entries which are arrays
    #normal creating TH2F
    
    tmp = ROOT.TH2F(varnamex.lower()+"over"+varnamey.lower()+dif, 
                    varnamex.lower()+"over"+varnamey.lower()+dif,
                    len(binsx)-1,binsx, len(binsy)-1, binsy)
    tmp.Sumw2()
    run = 1

    for s in range(tree.GetEntries()):
        tree.GetEntry(s)
        list1 = getattr(tree, varnamex)
        list2 = getattr(tree, varnamey)
        
        #if len(list1)!=len(list2):
            #print("lists must be of same length")
            #return(False, [False, False])
        if list1[1]==-999 or list2[2]==-999:
            continue
        else:
            for index, (elementx, elementy) in enumerate(zip(list1, list2)):
                if elementx == -999 or elementy == -999:
                    break
                else:
                    tmp.Fill(elementx, elementy)

    return(tmp)

def createTH2F_notarrays(tree, varnamex, varnamey,  binsx, binsy, dif = "", ind = ""): #currently kind of hard coded to work for entries which are not arrays
    #normal creating TH2F
    
    tmp = ROOT.TH2F(varnamex.lower()+"over"+varnamey.lower()+dif, 
                    varnamex.lower()+"over"+varnamey.lower()+dif,
                    len(binsx)-1,binsx, len(binsy)-1, binsy)
    tmp.Sumw2()
    run = 1

    for s in range(tree.GetEntries()):
        tree.GetEntry(s)
        elementx = getattr(tree, varnamex)
        elementy = getattr(tree, varnamey)
        if ind != "":
            ind1 = getattr(tree, ind)[0]
            tmp.Fill(elementx[ind1], elementy[ind1])
        #if len(list1)!=len(list2):
            #print("lists must be of same length")
            #return(False, [False, False])
        else:
            tmp.Fill(elementx, elementy)

    return(tmp)
# ''.join(e for e in constraints if e.isalnum()) is used to remove special characters out of string, otherwise it doesnt work porperly

def saveTH2Fplot(hist, xtitle, ytitle):
    canvas = ROOT.TCanvas("c1","c1", 800, 600)
    hist.SetXTitle(xtitle)
    hist.SetYTitle(ytitle)
    canvas.SetRightMargin(0.2)
    ##pads = plot.OnePad()
    #axish = createAxisHists(1,hist,hist.GetXaxis().GetXmin(),hist.GetXaxis().GetXmax())
    ##axish[0].GetYaxis().SetRangeUser(hist.GetMinimum(),hist.GetMaximum())
    #axish[0].GetXaxis().SetTitle(xtitle)
    #axish[0].GetYaxis().SetTitle(ytitle)
    ##pads[0].cd()
    #axish[0].Draw()
    hist.Draw("COLZ")
    canvas.SaveAs(ytitle+"OVER"+xtitle+".pdf")
    canvas.Close()
    #del(axish)
    

def saveplot(hist,  xtitle, channel="", title="", linecolor=ROOT.kRed, stats = [], log = False, ytitle = "Entries"):    
    
    canvas = ROOT.TCanvas("c1","c1")
    pads = plot.OnePad()
    axish = createAxisHists(1,hist,hist.GetXaxis().GetXmin(),hist.GetXaxis().GetXmax()-0.01)
    axish[0].GetYaxis().SetRangeUser(0,1.5*hist.GetMaximum())
    axish[0].GetXaxis().SetTitle(xtitle)
    axish[0].GetYaxis().SetTitle(ytitle)
    pads[0].cd()
    
    if log:
        count1 = 0
        count2 = 0
        hist1bins = []
        hist1y = []
        hist2bins = []
        hist2y = []
        for i in range(hist.GetNbinsX()):
            # create arrays which describe new histogram having no zero entries for log scale
            if hist.GetBinContent(i+1)==0:
                if count1==0:
                    hist1bins.append(hist.GetXaxis().GetBinLowEdge(i+1))
                count1+=1
            else:
                if count1!=0:
                    hist1y.append(hist.GetBinContent(i+1)/count1)
                    count1 = 0
                else:
                    hist1y.append(hist.GetBinContent(i+1))
                    hist1bins.append(hist.GetXaxis().GetBinLowEdge(i+1))
            
            
        if count1!=0:
            hist1y[-1] = hist1y[-1]/count1
            del(hist1bins[-1])

        hist1bins.append(hist.GetXaxis().GetXmax())
        histtitle = hist.GetTitle()
        del(hist)
        hist = ROOT.TH1F(histtitle, histtitle, len(hist1y), np.array(hist1bins))
        for i in range(len(hist1y)):
            hist.SetBinContent(i+1, hist1y[i])


        pads[0].SetLogy()
        hist.GetYaxis().SetTitle(ytitle)
        hist.GetYaxis().SetRangeUser(1,100*hist.GetMaximum())
        hist.GetXaxis().SetTitle(xtitle)
        
        hist.Sumw2()
    else:
        axish[0].Draw()


    hist.SetLineColor(linecolor)
    hist.SetMarkerStyle(20)
    hist.SetMarkerColor(linecolor+3)

    hist.Draw("LSAME")
    if len(stats)==2:
        latex = ROOT.TText(.5,.5, "selected {} of {}".format(int(stats[1]), stats[0]))
        latex.SetX(hist1.GetXaxis().GetXmin()+0.05*abs(hist.GetXaxis().GetXmax()-hist.GetXaxis().GetXmin()))
        latex.SetY(1.35*hist.GetMaximum())
        latex.SetTextSize(20)
        latex.SetTextFont(43)
        latex.Draw("SAME")
        text
    if title!="":
        text2 = ROOT.TText(.5,.5, title)
        text2.SetNDC() #so it orientates on the canvas
        text2.SetX(0.25)
        text2.SetY(0.95)
        text2.SetTextSize(20)
        text2.SetTextFont(43)
        text2.Draw("SAME")
    
        
    

    legend = plot.PositionedLegend(0.2, 0.15,3,0.015)
    legend.AddEntry(hist,channel)
    if channel!="":
        legend.Draw("SAME")    
    canvas.SaveAs(''.join(e for e in channel if e.isalnum())+hist.GetTitle()+".pdf")
    canvas.Close()
    del(axish,pads, legend)


def calc_AMS(hist1,hist2, error=0, rd=-1, AMSmod = 1): # * calculate ams value of 2 histograms, first hist is signal
    final = 0
    rest = 0
    for i in range(hist1.GetNbinsX()):
        if hist2.GetBinContent(i+1)!=0:
            final+= ROOT.RooStats.AsimovSignificance(hist1.GetBinContent(i+1)+rest, hist2.GetBinContent(i+1), error)**2 # last entry is error, maybe add statistical error of histogramms
            rest = 0
        else:
            rest+=hist1.GetBinContent(i+1)
    final = np.sqrt(final*AMSmod)
    if rd == -1:
        if final<10 and final>1:
            final = round(final, 2)
        elif final<1:
            final = round(final,3)
        elif final>10 and final<100:
            final = round(final, 1)
        elif final>100:
            final = int(final)
    else:
        final = round(final, rd)
    return(final)

def savegraph( datax, datay, title, xtitle, ytitle,  text = "", heading = "" ): # get x and y data and plot it
    c1 = ROOT.TCanvas( 'c1', 'c1' )

    
    n = len(datax)
    x, y = array( 'd' ), array( 'd' )
    
    for i in range( n ):
        x.append(datax[i])
        y.append(datay[i])
        #print(' i %i %f %f ' % (i,x[i],y[i]))
    
    gr = ROOT.TGraph( n, x, y )
    gr.SetLineColor( 2 )
    gr.SetLineWidth( 4 )
    gr.SetMarkerColor( 2 )
    gr.SetMarkerStyle( 1 )
    gr.SetTitle( 'a simple graph' )
    gr.GetXaxis().SetTitle(xtitle )
    gr.GetYaxis().SetTitle(ytitle )
    gr.Draw('ACP' )
    
    if text!="":
        text1 = ROOT.TText(.5,.5, text)
        text1.SetNDC() #so it orientates on the canvas
        text1.SetX(0.18)
        text1.SetY(0.88)
        text1.SetTextSize(20)
        text1.SetTextFont(43)
        text1.Draw("SAME")
        
    if heading!="":
        text2 = ROOT.TText(.5,.5, heading)
        text2.SetNDC() #so it orientates on the canvas
        text2.SetX(0.2)
        text2.SetY(0.95)
        text2.SetTextSize(20)
        text2.SetTextFont(43)
        text2.Draw("SAME")
    
    # TCanvas.Update() draws the frame, after which one can change it
    c1.Modified()
    c1.Update()
    c1.SaveAs(title+".pdf")
    c1.Close()
    del(gr)

def saveratioplot(hist1, hist2, xtitle, title1, title2, title = "", linecolor1=ROOT.kRed, return_area=False,
linecolor2=ROOT.kBlue, stats = [], dif="", calcArea=False, calcAMS=False, AMSmod=1, text="", y_label = "Entries", ratio = True, log = False, name = ""):

    canvas = ROOT.TCanvas("c1","c1")

    axish = createAxisHists(2,hist1,hist1.GetXaxis().GetXmin(),hist1.GetXaxis().GetXmax()-0.01)
    axish[0].GetYaxis().SetRangeUser(0,1.5*max([hist1.GetMaximum(), hist2.GetMaximum()]))
    axish[0].GetYaxis().SetTitle(y_label)


    if ratio == True:
        ratio_central = plot.MakeRatioHist(hist1,hist2,True,False) # creates the histogram; takes two histograms; dont know what the boolean values mean tough
        pads = plot.TwoPadSplit(0.29,0.01,0.01)
        axish[0].GetXaxis().SetTitleSize(0)
        axish[0].GetXaxis().SetLabelSize(0)
        axish[1].GetXaxis().SetTitle(xtitle)
        axish[1].GetYaxis().SetTitle("Ratio")
        axish[1].GetYaxis().SetRangeUser(0.2,1.8)
        pads[1].cd()
        axish[1].Draw()
        ratio_central.SetLineColor(linecolor1)
        ratio_central.SetMarkerColor(linecolor2)
        ratio_central.Draw("SAME")
    else:
        axish[0].GetXaxis().SetTitle(xtitle)
        pads = plot.OnePad()

    pads[0].cd()

    if calcAMS:
        AMS = calc_AMS(hist1, hist2, AMSmod=AMSmod)

    if log:
        count1 = 0
        count2 = 0
        hist1bins = []
        hist1y = []
        hist2bins = []
        hist2y = []
        for i in range(hist1.GetNbinsX()):
            # create arrays which describe new histogram having no zero entries for log scale
            if hist1.GetBinContent(i+1)==0:
                if count1==0:
                    hist1bins.append(hist1.GetXaxis().GetBinLowEdge(i+1))
                count1+=1
            else:
                if count1!=0:
                    hist1y.append(hist1.GetBinContent(i+1)/count1)
                    count1 = 0
                else:
                    hist1y.append(hist1.GetBinContent(i+1))
                    hist1bins.append(hist1.GetXaxis().GetBinLowEdge(i+1))
            
            if hist2.GetBinContent(i+1)==0:
                if count2==0:
                    hist2bins.append(hist2.GetXaxis().GetBinLowEdge(i+1))
                count2+=1
            else:
                if count2!=0:
                    hist2y.append(hist2.GetBinContent(i+1)/count2)
                    count2 = 0
                else:
                    hist2y.append(hist2.GetBinContent(i+1))
                    hist2bins.append(hist2.GetXaxis().GetBinLowEdge(i+1))
        
        if count1!=0:
            hist1y[-1] = hist1y[-1]/count1
            del(hist1bins[-1])
        if count2!=0:
            hist2y[-1] = hist2y[-1]/count2
            del(hist2bins[-1])
        hist1bins.append(hist1.GetXaxis().GetXmax())
        hist2bins.append(hist2.GetXaxis().GetXmax())
        
        histtitle1 = hist1.GetTitle()
        histtitle2 = hist2.GetTitle()

        
        del(hist1, hist2)
        hist1 = ROOT.TH1F(histtitle1, histtitle1, len(hist1y), np.array(hist1bins))
        hist2 = ROOT.TH1F(histtitle2, histtitle2, len(hist2y), np.array(hist2bins))
        for i in range(len(hist1y)):
            hist1.SetBinContent(i+1, hist1y[i])
        for i in range(len(hist2y)):
            hist2.SetBinContent(i+1, hist2y[i])

        pads[0].SetLogy()
        hist1.GetYaxis().SetTitle(y_label)
        if log:
            hist1.GetYaxis().SetRangeUser(1,max([hist1.GetMaximum(), hist2.GetMaximum()])**(3/2))
        else:
            hist1.GetYaxis().SetRangeUser(0,1.5*max([hist1.GetMaximum(), hist2.GetMaximum()]))
        hist1.GetXaxis().SetTitle(xtitle)
        
        hist1.Sumw2()
        hist2.Sumw2()
    else:
        axish[0].Draw()


    hist1.SetLineColor(linecolor1)
    hist1.SetMarkerStyle(20)
    hist1.SetMarkerColor(linecolor1+3)
    hist2.SetMarkerStyle(22)
    hist2.SetLineColor(linecolor2)
    hist2.SetMarkerColor(linecolor2+3)




    hist1.Draw("SAME")#both histograms are drawn into the same pad
    hist2.Draw("SAME")
        

    legend = plot.PositionedLegend(0.4, 0.15, 3, 0.015)
    legend.AddEntry(hist1, title1)
    legend.AddEntry(hist2, title2)
    legend.Draw("SAME")
    
    
    #calculating enclosed areas
    if calcArea:
        final = 0
        for i in range(hist1.GetNbinsX()):
            final+= min(hist1.GetBinContent(i+1), hist2.GetBinContent(i+1)) 
            #final+= hist1.GetBinContent(i)
            #if "Jet_btagDeepB" == xtitle:
                #print(hist2.Integral())
        final = round(final, 2)
        text1 = ROOT.TText(.5,.5, "Overlapping Area: "+str(final))
        text1.SetNDC() #so it orientates on the canvas
        text1.SetX(0.18)
        text1.SetY(0.88)
        text1.SetTextSize(20)
        text1.SetTextFont(43)
        text1.Draw("SAME")

    if return_area:
        return final

    # calculating ams value
    if calcAMS:
        # if rdAMS<=0:
        #     text1 = ROOT.TText(.5,.5, "AMS: "+str(int(final)))    
        text1 = ROOT.TText(.5,.5, "AMS: "+str(AMS))
        text1.SetNDC() #so it orientates on the canvas
        text1.SetX(0.18)
        text1.SetY(0.88)
        text1.SetTextSize(20)
        text1.SetTextFont(43)
        text1.Draw("SAME")

    if text!="":
        text1 = ROOT.TText(.5,.5, text)
        text1.SetNDC() #so it orientates on the canvas
        text1.SetX(0.18)
        text1.SetY(0.88)
        text1.SetTextSize(20)
        text1.SetTextFont(43)
        text1.Draw("SAME")
    

        
    

    # print("signal integral {}, background integral {}".format(hist1.Integral(), hist2.Integral()))

    pads[0].cd()
    pads[0].GetFrame().Draw() #
    pads[0].RedrawAxis()# doesnt do anything I realised, so maybe it helps to minimize bugs.

    if len(stats)==2:
        latex = ROOT.TText(.5,.5, "selected {} of {}".format(int(stats[1]), stats[0]))
        latex.SetX(hist1.GetXaxis().GetXmin()+0.2)
        latex.SetY(1.1*hist1.GetMaximum())
        latex.SetTextSize(20)
        latex.SetTextFont(43)
        latex.Draw("SAME")
        
    if title!="":
        text = ROOT.TText(.5,.5, title)
        text.SetNDC() #so it orientates on the canvas
        text.SetX(0.25)
        text.SetY(0.95)
        text.SetTextSize(20)
        text.SetTextFont(43)
        text.Draw("SAME")
    
    diff = ''.join(e for e in dif if e.isalnum())
    if name == "":
        canvas.SaveAs(xtitle+"_"+title1+"_"+title2+"_"+diff+"_ratio.pdf")
    else:
        canvas.SaveAs(name+".pdf")
    canvas.Close()
    del(pads,legend, canvas)

if __name__=="__main__":
    # ------------------------------------------------------------------------------------
    # reading in data and definitions

    # var = [0"Jet_pt", 1"Jet_phi", 2"Jet_eta", 3"deltaR", 4"recon_Hmass", 5"nJets", 6"MET_pt", 7"fatjetsused", 8"original_Hmass", 9"genjets_Hmass", 10"nrightjets", 11"Jet_btagDeepB", 12"nmatch", 13"candidate_score_1", 14"wnrightjets", 15"genjets_deltaR_W", 16"recon_WMass", 17"genjets_WMass", 18"wnrightjets_bad"] "original_Hmass" , "nrightjets", "candidate_score_", "wnrightjets", "genjets_deltaR_W", "genjets_WMass", "wnrightjets_bad", "genjets_Hmass", "fatjetsused"
    var = ["Jet_pt",  "Jet_eta", "deltaR_bb", "recon_Hmass", "nJets", "MET_pt" ,  "Jet_btagDeepB", "nmatch",    
    "recon_Wmass",  "Jet_HT", "deltaR_W", "recon_Hpt", "recon_Wpt", "recon_Weta", "deltaR_reconHW", "MET_pt", "nSoftActivityJet"]

    #var = ["Jet_btagDeepC", "Jet_mass", "Jet_nElectrons", "Jet_hadronFlavour"]
    
    ranges = {"Jet_pt":np.linspace( 0,200, 40), "Jet_phi":np.linspace(-4,4,20), "Jet_eta":np.linspace(-4,4,20), 
    "deltaR_bb":np.linspace(0,6,50), "recon_Hmass":np.linspace(90, 140, 15), "nJets":np.linspace(3.5,12.5,10), 
    "MET_pt":np.linspace(0,250,70), "fatjetsused":np.linspace(0,2,3), "original_Hmass":np.linspace(0,200,50), 
    "genjets_Hmass:":np.linspace(50, 200, 50), "nrightjets":np.linspace(-0.5,2.5,4), "Jet_btagDeepB":np.linspace(0.5,1,50), 
    "nmatch":np.linspace(1.5,6.5,6), "candidate_score_1":np.linspace(1.6,8,20), "wnrightjets":np.linspace(-0.5,2.5,4), 
    "genjets_deltaR_W":np.linspace(0,10,20), "bckgrdrecon_Wmass":np.linspace(0,150, 15), "genjets_WMass":np.linspace(0,150, 40), 
    "wnrightjets_bad":np.linspace(-0.5,2.5,4), "Jet_HT":np.linspace(100, 600, 70), "weight":np.linspace(0, 10, 50), 
    "deltaPT_W":np.linspace(0, 300, 50), "deltaR_W":np.linspace(0,7,50), "recon_Wpt":np.linspace(0,200,50), 
    "recon_Hpt":np.linspace(0,200,50), "recon_Heta":np.linspace(-3,3,20), "recon_Weta":np.linspace(.3,3), 
    "deltaR_reconHW":np.linspace(0.5,4.5,50), "Jet_btagDeepC":np.linspace(0,1,30), "Jet_mass":np.linspace(20,100,70), 
    "Jet_nElectrons":np.linspace(-0.5,5.5,7), "Jet_hadronFlavour":np.linspace(-20,20,42), "nSoftActivityJet":np.linspace(-0.5,10.5,12),
    "recon_Wmass":np.linspace(20, 140, 15)}

    triggers = ["HLT_AK8PFJet80",  "HLT_PFJet80", "HLT_PFHT300PT30_QuadPFJet_75_60_45_40", "HLT_DiPFJetAve35_HFJEC",  "HLT_PFHT180", "HLT_DoublePFJets40_CaloBTagCSV_p33",  "HLT_AK4PFJet80", "HLT_AK4CaloJet80"]
    selected_triggers = ["HLT_PFHT300PT30_QuadPFJet_75_60_45_40"]#, "HLT_PFHT180" ]

    standardconstraintsH = "deltaR_bjet1 < .3 && deltaR_bjet2 < .3" #first selection is to used to ensure that i can connect the genpart level with the genjet level properly; 
    testing_constraints = "deltaR_bb<1. && deltaR_bb>.0"
    
    standardconstraintsW = "deltaR_W_jet1< .5 && deltaR_W_jet2 < .5"

    infile = ROOT.TFile.Open("../classifiers/data/breg/bckgrd.root")
    bckgrd = infile.Get("Events")
    
    infile2 = ROOT.TFile.Open("../classifiers/data/breg/signal.root")
    signal = infile2.Get("Events")
    
    #infile = ROOT.TFile.Open("../WplusH_HToBB_WToQQ_M125_Skim.root")
    #signal = infile.Get("Events")enumerate

    #infile1 = ROOT.TFile.Open("../WminusH_HToBB_WToQQ_M125_Skim.root")
    #signal1 = infile1.Get("Events")

    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! random stuff -----------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    """
    var = "deltaR_bb"
    signal1, _ = createTH1F(signal, var, ranges[var], "weight", dif = "signal", normalise = True)
    bckgrd1, _ = createTH1F(bckgrd, var, ranges[var], "weight", normalise = True)
    saveratioplot(signal1, bckgrd1,  var,  "signal", "background", calcArea = True, ratio = False)
    """
    """
    further_variables = {"Jet_hadronFlavour":np.linspace(-0.5, 6.5, 8),
    "Jet_btagDeepC":np.linspace(0,1,20), "Jet_mass":np.linspace(0,300, 30)}

    ranges.update(further_variables)

    for var in further_variables.keys():        
        signal1, _ = createTH1F(signal, "{}[recon_bjet_id[0]]".format(var), ranges[var], "weight", dif = "signal", normalise = True)
        bckgrd1, _ = createTH1F(bckgrd, "{}[recon_bjet_id[0]]".format(var), ranges[var], "weight", normalise = True)
        saveratioplot(signal1, bckgrd1,  var,  "signal", "background", calcArea = True, ratio = False)

    """
    #### recon_WJets_id, recon_bjet_id
    #for futher input W Jet_neEmEF, W Jet_neHEF, W Jet_qgl, ?H Jet_qgl, ? W Jet_nConstituents ?, ?Jet-HT?

    var = "Jet_btagDeepB"
    
    ranges.update({"Jet_chEmEF":np.linspace(0,1,20), "Jet_chHEF":np.linspace(0,1,20), "Jet_nConstituents":np.linspace(0, 60, 30),
    "Jet_neEmEF":np.linspace(0,1,20), "Jet_neHEF":np.linspace(0,1,20),"Jet_qgl":np.linspace(-1,1,30),
    "Jet_nElectrons":np.linspace(0,5,6), "Jet_nMuons":np.linspace(0,5,6), "Jet_HT":np.linspace(0,700,30)})

    signal1, _ = createTH1F(signal, "{}[recon_WJets_id[0]]".format(var), ranges[var], "weight", dif = "signal", normalise = True)
    bckgrd1, _ = createTH1F(bckgrd, "{}[recon_WJets_id[0]]".format(var), ranges[var], "weight", normalise = True)
    saveratioplot(signal1, bckgrd1,  var,  "signal", "background", calcArea = True, ratio = False, name = "{}_area_W".format(var))
    del(signal1, bckgrd1)
    signal1, _ = createTH1F(signal, "{}[recon_bjet_id[0]]".format(var), ranges[var], "weight", dif = "signal", normalise = True)
    bckgrd1, _ = createTH1F(bckgrd, "{}[recon_bjet_id[0]]".format(var), ranges[var], "weight", normalise = True)
    saveratioplot(signal1, bckgrd1,  var,  "signal", "background", calcArea = True, ratio = False, name = "{}_area_H".format(var))


    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! create plots of genlevel data and reconlevel data ----------------------------------
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    """
    var = "deltaR_bb"
    signal1, _ = createTH1F(signal, "recon_Hmass", np.linspace(0,400, 50), "(nrightjets>1)*weight", dif = "signal")
    bckgrd1, _ = createTH1F(signal, "recon_Wmass", np.linspace(0,400, 50), "(nrightjets>1)*weight")
    saveratioplot(signal1, bckgrd1,"W mass", "gen-level", "recon-level", calcArea = False, ratio = False, log = True)
    
    """

    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #! creating plots of no selection applied and calculate AMS value: (now triggers are applied additionally)
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    """
    var = "recon_Hmass"
    # "({}>0)*weight".format(j)
    # "weight"


    for j in selected_triggers:
        signal1, _ = createTH1F(signal, var, ranges[var], "({}>0)*weight".format(j) , dif = "signal")
        bckgrd1, _ = createTH1F(bckgrd, var, ranges[var], "({}>0)*weight".format(j))
        signal1.Add(bckgrd1)
        saveplot(signal1, var,  "signal + background", title = j)
        #AMS = calc_AMS(signal1, bckgrd1)
        #del(signal1, bckgrd1)
        #signal1, _ = createTH1F(signal, var, ranges[var], "weight", normalise = True, dif = "signal")
        #bckgrd1, _ = createTH1F(bckgrd, var, ranges[var], "weight", normalise = True)
        #saveratioplot(signal1, bckgrd1, var, "signal", "background", calcArea = False, text = "AMS {}".format(AMS), dif = j)
        del(signal1, bckgrd1)
    """
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------
    #* interesting constraints:
    #* ------------------------------------------------------------------------------------
    #* ------------------------------------------------------------------------------------


    # "deltaR_bjet1 < .3 & deltaR_bjet2 < .3" # deltar between genjet and genpart b is not allowed to be too big
    # "nrightjets == 2" # right jets used
    # "MET_pt<70" # Missing Energy is supposed to not be too big
    # "fatjetsused == 0" # no fatjets are used
    
    # calculating fractions after reweighting:
    #tmp1, stats1 = createTH1F(bckgrd,"recon_Hmass",ranges["recon_Hmass"], constraints= "weight", dif="bckgrd")
    #tmp2, stats1 = createTH1F(signal,"recon_Hmass",ranges["recon_Hmass"], constraints= "weight")
    #print("signal {}; background {}".format(tmp2.Integral(), tmp1.Integral() ))
    
    #showing smooth transition of LHE_HT ()

    # tmp1, stats1 = createTH1F(bckgrd,"LHE_HT",np.linspace(300,800,100), constraints= "weight", dif="bckgrd")
    # saveplot(tmp1, "LHE_HT", log = True)
    

    # cur1= 10
    # cur2 =14
    
    #signal1, stats1 = createTH1F(signal ,"nrightjets", ranges["nrightjets"])
    #saveplot(signal1, "nrightjets", title = "H Jets")
    #bckgrd1, stats1 = createTH1F(signal ,"wnrightjets", ranges["nrightjets"])
    #saveplot(bckgrd1, "nrightjets", title = "W Jets")

    
    
    """
    j = "deltaR_reconHW <4.5 && deltaR_reconHW>2."
    signal1, stats1 = createTH1F(bckgrd,triggers[0], np.linspace(-0.5,1.5,2),constraints= "("+j+">0)*weight", dif="signal")
    ints1 = signal1.Integral()
    del(signal1)
    signal2, stats1 = createTH1F(bckgrd,triggers[0], np.linspace(-0.5,1.5,2),constraints= "weight")
    ints2 = signal2.Integral()
    del(signal2)    
    print(ints1/ints2)
    """
    
    #tmp1 = createTH2F_notarrays(signal,"Jet_btagDeepB", "Jet_btagDeepC",np.linspace(0,1,50),np.linspace(0,1,50), ind ="recon_bjet_id")
    #saveTH2Fplot(tmp1, "Jet_btagDeepB", "Jet_btagDeepC")
    
    
    #tmp, trash = createTH1F(signal, "deltaPT_W", np.linspace(0,200,50), "nrightjets==2")
    #tmp2, trash = createTH1F(bckgrd, "deltaPT_W", np.linspace(0,200,50))
    #saveratioplot(tmp, tmp2, "deltaPT_W", "rightly assigned Jets", "background")
    #saveplot(tmp, "deltaPT_W", title="deltaPT_W")
    """# specific TH2F creation, using ids with arrays
    varnamex = "Jet_btagDeepB"
    varnamey = "isright"
    binsx = np.linspace(0,1,50)
    binsy = np.linspace(-0.5,1.5,3)
    
    
    tmp = ROOT.TH2F(varnamex.lower()+"over"+varnamey.lower(), 
                    varnamex.lower()+"over"+varnamey.lower(),
                    len(binsx)-1,binsx, len(binsy)-1, binsy)
    tmp.Sumw2()

    #signal = infile.Get("Events")enumerate

    #infile1 = ROOT.TFile.Open("../WminusH_HToBB_WToQQ_M125_Skim.root")
    #signal1 = infile1.Get("Events")
    
    #* ------------------------------------------------------------------------------------
    #* interesting constraints:



    # "deltaR_bjet1 < .3 & deltaR_bjet2 < .3" # deltar between genjet and genpart b is not allowed to be too big
    # "nrightjets == 2" # right jets used
    # "MET_pt<70" # Missing Energy is supposed to not be too big
    # "fatjetsused == 0" # no fatjets are used
    for s in range(tree.GetEntries()):
        tree.GetEntry(s)
        elementx = getattr(tree, varnamex)
        elementy = getattr(tree, varnamey)
        ids = getattr(tree, "recon_WJets_id")
        tmp.Fill(elementx[ids[0]], elementy[0])
        tmp.Fill(elementx[ids[1]], elementy[1])    
    
    saveTH2Fplot(tmp, "Jet_btagDeepB", "wnrightjets")
    """
    #31.07.19 creating nrights wnrights
    #tmp1, trash1 = createTH1F(tree, var[cur1], ranges[var[cur1]], standardconstraintsH)    
    #tmp2, trash2 = createTH1F(tree, var[cur2], ranges[var[cur2]], standardconstraintsW)
    #saveplot(tmp1, var[cur1], "deltaR_H < 0.3", stats=trash1)
    #saveplot(tmp2, var[cur2], "deltaR_W < 0.3", stats=trash2)
    #del(tmp1, tmp2)
    
    # check wether there is a shape changing regarding triggers:
    
    #i = "recon_Hmass"

    #*--------------------------------------------------------------------------------------------
    #*--------------------------------------------------------------------------------------------
    #! creating plots for presentation -> showing distribution of signal AND background, and both normalised
    #*--------------------------------------------------------------------------------------------
    #*--------------------------------------------------------------------------------------------
    """
    rangeH = np.linspace(90, 140, 15)

    for j in selected_triggers:
        signal1, trash = createTH1F(signal, "recon_Hmass", ranges["recon_Hmass"], constraints= "("+j+">0)*weight", dif = "signal")
        bckgrd1, trash = createTH1F(bckgrd, "recon_Hmass", ranges["recon_Hmass"], constraints= "("+j+">0)*weight", dif = "bckgrd")
        AMSval = calc_AMS(signal1, bckgrd1)
        signal1.Add(bckgrd1)
        saveplot(signal1, "recon_Hmass", "Data", j)

        del(signal1, bckgrd1)
        signal1, trash = createTH1F(signal, "recon_Hmass", ranges["recon_Hmass"], constraints= "("+j+">0)*weight", dif = "signal", normalise=True)
        bckgrd1, trash = createTH1F(bckgrd, "recon_Hmass", ranges["recon_Hmass"], constraints= "("+j+">0)*weight", dif = "bckgrd", normalise=True)
        saveratioplot(signal1, bckgrd1, "recon_Hmass", "Signal", "Background", calcArea=False, text = "AMS: {}".format(AMSval), dif=j)
        del(signal1, bckgrd1)

    """
    #*--------------------------------------------------------------------------------------------
    #*--------------------------------------------------------------------------------------------
    #! checking trigger ratios (again) -----------------------------------------------------------
    #! little bit of an overkill here but one has to be careful since some bins are not displayed in the histogram
    #*--------------------------------------------------------------------------------------------
    #*--------------------------------------------------------------------------------------------
    
    """
    for j in triggers:
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
    """    
    
    #*--------------------------------------------------------------------------------------------
    #*--------------------------------------------------------------------------------------------        
    #! checking PDFS for signal and background for variables for selected triggers including only the selection for selected Jets only -----------------------------------------
    #*--------------------------------------------------------------------------------------------
    #*--------------------------------------------------------------------------------------------
    """
    os.chdir("documentation/selectingparameters/lessbins")
    cuts = ""
    #cuts = "deltaR_reconHW <4.5 && deltaR_reconHW>2. &&"
    for j in selected_triggers:
            for i in var:
                if "Jet" in i and i!="Jet_HT" and i!="nJets" and i!="nSoftActivityJet":
                    signal1, trash = createTH1F(signal, i+"[recon_bjet_id[0]]", ranges[i],constraints= "("+cuts+j+">0)*weight", normalise = True, dif=str(np.random.rand()))
                    signal2, trash = createTH1F(signal, i+"[recon_bjet_id[1]]", ranges[i], constraints= "("+cuts+j+">0)*weight", normalise = True, dif=str(np.random.rand()))

                    signal3, trash = createTH1F(signal, i+"[recon_WJets_id[0]]", ranges[i],constraints=  "("+cuts+j+">0)*weight", normalise = True, dif=str(np.random.rand()))
                    signal4, trash = createTH1F(signal, i+"[recon_WJets_id[1]]", ranges[i], constraints=  "("+cuts+j+">0)*weight", normalise = True, dif=str(np.random.rand()))

                    
                    bckgrd1, trash = createTH1F(bckgrd, i+"[recon_bjet_id[0]]", ranges[i],constraints=  "("+cuts+j+">0)*weight", normalise = True, dif=str(np.random.rand()))
                    bckgrd2, trash = createTH1F(bckgrd, i+"[recon_bjet_id[1]]", ranges[i], constraints=  "("+cuts+j+">0)*weight", normalise = True, dif=str(np.random.rand()))

                    bckgrd3, trash = createTH1F(bckgrd, i+"[recon_WJets_id[0]]", ranges[i],constraints= "("+cuts+j+">0)*weight", normalise = True, dif=str(np.random.rand()))
                    bckgrd4, trash = createTH1F(bckgrd, i+"[recon_WJets_id[1]]", ranges[i], constraints= "("+cuts+j+">0)*weight", normalise = True, dif=str(np.random.rand()))

                    #save first plots
                    saveratioplot(signal1, bckgrd1, i, "signal", "background", title = j+ "; b Jet higher pt", dif = j+"bjhighpt")
                    saveratioplot(signal2, bckgrd2, i, "signal", "background", title = j+ "; b Jet lower pt", dif = j+"bjlowpt")
                    
                    saveratioplot(signal3, bckgrd3, i, "signal", "background", title = j+ "; W Jet higher pt", dif = j+"Wjhighpt")
                    saveratioplot(signal4, bckgrd4, i, "signal", "background", title = j+ "; W Jet lower pt", dif = j+"Wjlowpt")



                    # adding up
                    signal1.Add(signal2)
                    signal1.Scale(0.5)
                    bckgrd1.Add(bckgrd2)
                    bckgrd1.Scale(0.5)
                    
                    signal3.Add(signal4)
                    signal3.Scale(0.5)
                    bckgrd3.Add(bckgrd4)
                    bckgrd3.Scale(0.5)
                    
                    del(signal2, bckgrd2, signal4 ,bckgrd4)
                    
                    
                    saveratioplot(signal1, bckgrd1, i, "signal", "background", title = j+ "; only H Jets", dif = j+"onlyH")
                    saveratioplot(signal3, bckgrd3, i, "signal", "background", title = j+ "; only W Jets", dif = j+"onlyW")
                    signal1.Add(signal3)
                    signal1.Scale(0.5)
                    bckgrd1.Add(bckgrd3)
                    bckgrd1.Scale(0.5)
                    del(signal3,bckgrd3)
                    saveratioplot(signal1, bckgrd1, i, "signal", "background", title = j+ "; all selected Jets", dif = j+"alljets")
                    del(signal1,bckgrd1)
                
                else:
                    signal1, stats1 = createTH1F(signal, i, ranges[i],constraints= "("+cuts+j+">0)*weight", normalise = True, dif="signal")
                    bckgrd1, trash = createTH1F(bckgrd, i, ranges[i], constraints= "("+cuts+j+">0)*weight", normalise = True)
                    #tmp2, trash = createTH1F(tree, var[17],ranges[var[17]])
                    #saveplot(tmp1, i, "background", stats= stats1)
                    #saveplot(tmp2, i, "signal", stats=trash)
                    saveratioplot(signal1,bckgrd1, i, "signal","background", title = j , dif = j)
                    #print("signal: {}; background: {}".format(trash[1]/trash[0], stats1[1]/stats1[0]))
                    del(signal1, bckgrd1)
    """
    """
    #*--------------------------------------------------------------------------------------------
    #*--------------------------------------------------------------------------------------------        
    #! checking shape changing of HLT ------------------------------------------------------------
    #*--------------------------------------------------------------------------------------------
    #*--------------------------------------------------------------------------------------------         
    i = "recon_Hmass"
    for j in selected_triggers:
        tmp1, stats1 = createTH1F(signal, i, ranges[i],constraints= "weight", normalise = True, dif="signal")
        tmp2, trash = createTH1F(signal, i, ranges[i], constraints= "("+j+">0)*weight", normalise = True)
        saveratioplot(tmp1, tmp2, i, "no trigger", j, title = "normed distributions", dif = j, calcArea=False)
        del(tmp1,tmp2)
    """
    #*--------------------------------------------------------------------------------------------
    #*--------------------------------------------------------------------------------------------
    #*--------------------------------------------------------------------------------------------
    #*--------------------------------------------------------------------------------------------

    #saveplot(tmp1, var[16], "100to200")
    #saveplot(tmp2, var[4],  "100to200")
    
    """

    tmp1, trash1 = createTH1F(tree, "Jet_btagDeepB[recon_WJets_id[0]]",np.linspace(0,1,50),"isright[0]>0 &&"+  standardconstraintsH+"&&"+standardconstraintsW)
    tmp2, trash1 = createTH1F(tree, "Jet_btagDeepB[recon_WJets_id[1]]",np.linspace(0,1,50),"isright[1]>0 &&"+  standardconstraintsH+"&&"+standardconstraintsW)
    tmp2.Add(tmp1)
    tmp3, trash1 = createTH1F(tree, "Jet_btagDeepB[recon_WJets_id[0]]",np.linspace(0,1,50),"isright[0]==0 &&"+  standardconstraintsH+"&&"+standardconstraintsW, dif ="s")
    tmp4, trash1 = createTH1F(tree, "Jet_btagDeepB[recon_WJets_id[1]]",np.linspace(0,1,50),"isright[1]==0 &&"+  standardconstraintsH+"&&"+standardconstraintsW, dif = "s")
    tmp4.Add(tmp3)
    
    saveratioplot(tmp2,tmp4, "Jet_btagDeepB", "rightly assigned", "wrongly assigned")
    """
    
    #tmp2, trash = createTH1F(tree, var[cur2], ranges[var[cur2]], ""+standardconstraintsH+"&&deltaR_W_jet1< .3 && deltaR_W_jet2 < .3")    
    #saveratioplot(tmp2, tmp4, var[cur1], "correct assigned jets", "wrongly assigned jets")
    #saveplot(tmp2, var[11], "selected W jets")
    #tmp = createTH2F(tree, var[11], var[10], ranges[var[11]], ranges[var[10]])
    #filein = ROOT.TFile("nrightsOVERbdeep.root", "r")
    #tmp = ROOT.TH2F()
    #filein.GetObject("jet_btagdeepbovernrightjets_", tmp)
    #saveTH2Fplot(tmp, var[11], var[10])

    #tmp, stats1 = createTH1F(tree, var[12], ranges[var[12]], "nmatch>1 &&"+standardconstraints)
    #saveplot(tmp, var[12], "")
    #del(tmp)
    
    
    #!analyzing different picks ---------------------------------
    #*besties: best_hmass110, best_delta2, best_btagdeep, best_jetpt65
    """
    tmp1, trash = createTH1F(tree, var[10]+"[best_candidate_id_1]", ranges[var[10]], "nmatch>1&&"+standardconstraints)
    tmp2, trash = createTH1F(tree, var[10]+"[best_candidate_id_2]", ranges[var[10]], "nmatch>1&&"+standardconstraints)
    tmp3, trash = createTH1F(tree, var[10]+"[0]", ranges[var[10]], "nmatch>1&&"+standardconstraints)
    
    
    tmp4, trash = createTH1F(tree, var[10]+"[best_hmass110]", ranges[var[10]], "nmatch>1&&"+standardconstraints)    
    tmp5, trash = createTH1F(tree, var[10]+"[best_delta2]", ranges[var[10]], "nmatch>1&&"+standardconstraints)     
    tmp6, trash = createTH1F(tree, var[10]+"[best_btagdeep]", ranges[var[10]], "nmatch>1&&"+standardconstraints) 
    tmp7, trash = createTH1F(tree, var[10]+"[best_jetpt65]", ranges[var[10]], "nmatch>1&&"+standardconstraints) 
    
    saveratioplot(tmp4, tmp3, var[10], "best_hmass110", "Max pt")
    saveratioplot(tmp4, tmp2, var[10], "best_hmass110", "Method 2")
    saveratioplot(tmp4, tmp5, var[10], "best_hmass110", "best_delta2")
    saveratioplot(tmp6, tmp3, var[10], "highest btag", "Max pt")
    saveratioplot(tmp4, tmp7, var[10], "best_hmass110", "best_jetpt65")

    
    #tmp1, trash = createTH1F(tree, var[10]+"[0]", ranges[var[10]], "nmatch>1&&"+standardconstraints)
    #tmp2, trash = createTH1F(tree, var[10]+"[best_candidate_id]", ranges[var[10]], "nmatch>1&&"+standardconstraints)
    #saveratioplot(tmp2, tmp1, var[10], "described picking", "max pt")
    

    #fout = ROOT.TFile("nrightsOVERbdeep.root","recreate")
    #fout.cd()
    #tmp.Write()
    #fout.Close()
    """
    #*--------------------------------------------------------------
    #*recon vs gen -----------------------------
    #tmp1, trash = createTH1F(tree, var[4]+"[best_btagdeep]", ranges[var[4]], standardconstraints)
    #tmp2, trash = createTH1F(tree, var[9], ranges[var[4]], standardconstraints)
    #saveratioplot(tmp2, tmp1, "Hmass", "GenJet Hmass", "ReconJet Hmass")
    
    
    
    #* --------------------------------------------------------------------------------
    #* sepecific tasks
    """

    tmp1 = createTH1F(tree, var[9], ranges[var[4]], "deltaR_bjet1 < .3 & deltaR_bjet2 < .3 && nrightjets == 2 && MET_pt<30")
    tmp2 = createTH1F(tree, var[4], ranges[var[4]], "deltaR_bjet1 < .3 & deltaR_bjet2 < .3 && nrightjets == 2 && MET_pt<30")
    saveratioplot(tmp1,tmp2, "H_mass", "genjets", "reconstructed")
    del(tmp1, tmp2)


    tmp1 = createTH1F(tree, var[9], ranges[var[4]], "fatjetsused>0")
    tmp2 = createTH1F(tree, var[4], ranges[var[4]], "fatjetsused>0")
    saveratioplot(tmp1,tmp2, "H_mass", "corresponding H_mass", "only fatjets")
    del(tmp1, tmp2)

    tmp1 = createTH1F(tree, var[4], ranges[var[4]], "")
    tmp2 = createTH1F(tree, var[4], ranges[var[4]], Triggers[-1]+">0")
    saveratioplot(tmp1,tmp2, "H_mass", "selection", "selection + "+Triggers[-1])


    #comparison with actual value, considering if it is a bjet
    tmp1 = createTH1F(tree, "Jet_pt", ranges["Jet_pt"], "GenJet_hadronFlavour==5")
    tmp2 = createTH1F(tree2, "Jet_pt", ranges["Jet_pt"], "GenJet_hadronFlavour==5", "x")

    saveratioplot(tmp2,tmp1, "Jet_pt only actual b jets",  "original data", "skimmed data")


    # checking btag selection
    tmp1 = createTH1F(tree, var[4], ranges[var[4]], "")
    tmp2 = createTH1F(tree, var[4], ranges[var[4]], "Jet_btagDeepC[asso_bjet_H_id_recon1]>0 && Jet_btagDeepC[asso_bjet_H_id_recon2]>0")
    saveratioplot(tmp1,tmp2, "H_mass", "no constraints", "Jet_btagDeepC>0")
    del(tmp1, tmp2)

    #comparing recons_Hmass with genjet_Hmass
    #tmp1, stats1 = createTH1F(tree, var[9], ranges[var[4]], standardconstraints + "&&Sum$("+ testing_constraints+")>0")
    #tmp2, stats2 = createTH1F(tree, var[4], ranges[var[4]], standardconstraints + "&&"+ testing_constraints)

    tmp2, stats2 = createTH1F(tree, var[10], ranges[var[10]], standardconstraints + "&& deltaR_bb > 0")
    #for x in [0.5,1,1.5,2,2.5,3]:
    tmp3, stats3 = createTH1F(tree, var[10], ranges[var[10]], standardconstraints + "&& deltaR_bb > 0 && deltaR_bb>5")
    #saveratioplot(tmp2, tmp3, var[10], "no constraints", r"\delta R < "+ str(x), stats=stats3)
    saveplot(tmp3, var[10], r"\delta R > 5", stats=stats3)
    del(tmp3)
    #saveratioplot(tmp1,tmp2, "H_mass", "gen level level", "recon level")
    #del(tmp1, tmp2)



    """
    #* --------------------------------------------------------------------------------
    #* iterative analysis

    """ # trigger testing on Hmass
    for element in Triggers:
        #tmp = createTH1F(tree, element, ranges[element], "")
        #saveplot(tmp, element, "Z(hh)H(bb)" )
        #del(tmp) #important for not having a memory leak and buggy output
        tmp1  = createTH1F(tree, var[4], ranges[var[4]], standardconstraints)
        tmp2  = createTH1F(tree, var[4], ranges[var[4]], element+">0 && "+ standardconstraints)
        saveratioplot(tmp1,tmp2, var[4], "no fatjetsused", element)
        del(tmp1,tmp2)
    """


    #* ------------------------------------------------------------------------------------
    #* reference stuff

    #for testing stuff
    ##fout = ROOT.TFile("test.root","recreate")
    ##fout.cd()
    ##tmp1.Write()
    ##fout.Close()
    #saveplot(tmp1, var[9], "Z(hh)H(bb)" )
    #-----

    #tmp1 = createTH1F(tree, var[-1], ranges["original_Hmass"], "deltar_bjet1<0.3 && deltar_bjet2<0.3")
    #tmp2 = createTH1F(tree, var[-1], ranges["original_Hmass"], Triggers[2]+">0 && deltar_bjet1<0.3 && deltar_bjet2<0.3")
    #saveratioplot(tmp1,tmp2, var[-1], "no_HLT", Triggers[2])

    """ #definition of variables and Triggers to use
    var = ["Jet_pt", "Jet_phi", "Jet_eta", "deltaR"]
    Triggers = []

    xbins = [0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,210.,220.,240.,260.,280.,300.,350.,400.,450.,500.]

    jet_pt = ROOT.TH1F("jet_pt", "jet_pt", 20, 0, 200)
    jet_pt.Sumw2() #is this necessary tough? yes otherwise the lines will be connected in the plot but i dont know why
    jet_phi = ROOT.TH1F("jet_phi", "jet_phi", 20, -4, 4)
    jet_phi.Sumw2()
    jet_eta = ROOT.TH1F("jet_eta", "jet_eta",20, -4, 4)
    jet_eta.Sumw2()
    jet_dr = ROOT.TH1F("jet_dr", "jet_dr", 20, 0, 10)
    jet_dr.Sumw2()
    #ele_pt = ROOT.TH1F("ele_pt","ele_pt",20,0,100) #For 20 bins between 0 and 100 GeV 
    #ele_pt.Sumw2()
    #h_mass = ROOT.TH1F("h_mass","h_mass",len(xbins)-1,np.array(xbins)) #To use variable bin width based on the array xbins
    #h_mass.Sumw2()

    # maybe try doing it with a loop?
    #for index, element in enumerate(var):
    #    currentvar = ROOT.TH1F(element+str(1), element+str(1))
    #    tree.Draw(element+">>"+element=str(1), "")
    #    currentvar.SetLineColor(ROOT.kBlue)
        
    tree.Draw("Jet_pt>>jet_pt","")
    tree.Draw("Jet_phi>>jet_phi","")
    tree.Draw("Jet_eta>>jet_eta","")
    tree.Draw("deltaR>>jet_dr","")

    #tree.Draw("Electron_pt>>ele_pt","")
    #tree.Draw("H_mass>>h_mass","")

    #ele_pt.SetLineColor(ROOT.kRed)
    #h_mass.SetLineColor(ROOT.kBlue)

    #jet_pt
    canvas = ROOT.TCanvas("c1","c1")
    pads = plot.OnePad()
    axish = createAxisHists(1,jet_pt,jet_pt.GetXaxis().GetXmin(),jet_pt.GetXaxis().GetXmax()-0.01)
    axish[0].GetYaxis().SetRangeUser(0,1.5*jet_pt.GetMaximum())
    axish[0].GetXaxis().SetTitle("Jet p_{T} [GeV]")
    axish[0].GetYaxis().SetTitle("Entries")

    pads[0].cd()
    axish[0].Draw()
    jet_pt.Draw("LSAME")
    legend = plot.PositionedLegend(0.2, 0.15,3,0.015)
    legend.AddEntry(jet_pt,"Z(HH)H(bb)")
    legend.Draw("SAME")

    canvas.SaveAs("Jet_pt.pdf")

    canvas1 = ROOT.TCanvas("c2","c2")
    pads1 = plot.OnePad()
    axish1 = createAxisHists(1,jet_phi,jet_phi.GetXaxis().GetXmin(),jet_phi.GetXaxis().GetXmax()-0.01)
    axish1[0].GetYaxis().SetRangeUser(0,1.5*jet_phi.GetMaximum())
    axish1[0].GetXaxis().SetTitle("Jet phi")
    axish1[0].GetYaxis().SetTitle("Entries")

    pads1[0].cd()
    axish1[0].Draw()
    jet_phi.Draw("LSAME")
    legend1 = plot.PositionedLegend(0.2, 0.15,3,0.015)
    legend1.AddEntry(jet_phi,"Z(ll)H(bb)")
    legend1.Draw("SAME")

    canvas1.SaveAs("jet_phi.pdf")


    canvas2 = ROOT.TCanvas("c3","c3")
    pads2= plot.OnePad()
    axish2= createAxisHists(1,jet_eta,jet_eta.GetXaxis().GetXmin(),jet_eta.GetXaxis().GetXmax()-0.01)
    axish2[0].GetYaxis().SetRangeUser(0,1.5*jet_eta.GetMaximum())
    axish2[0].GetXaxis().SetTitle("Jet eta")
    axish2[0].GetYaxis().SetTitle("Entries")

    pads2[0].cd()
    axish2[0].Draw()
    jet_eta.Draw("LSAME")
    legend2 = plot.PositionedLegend(0.2, 0.15,3,0.015)
    legend2.AddEntry(jet_eta,"Z(ll)H(bb)")
    legend2.Draw("SAME")

    canvas2.SaveAs("jet_eta.pdf")


    canvas3 = ROOT.TCanvas("c4","c4")
    pads3 = plot.OnePad()
    axish3 = createAxisHists(1,jet_dr,jet_dr.GetXaxis().GetXmin(),jet_dr.GetXaxis().GetXmax()-0.01)
    axish3[0].GetYaxis().SetRangeUser(0,1.5*jet_dr.GetMaximum())
    axish3[0].GetXaxis().SetTitle("\delta R [GeV]")
    axish3[0].GetYaxis().SetTitle("Entries")

    pads3[0].cd()
    axish3[0].Draw()
    jet_dr.Draw("LSAME")
    legend3 = plot.PositionedLegend(0.2, 0.15,3,0.015)
    legend3.AddEntry(jet_dr,"Z(ll)H(bb)")
    legend3.Draw("SAME")

    canvas3.SaveAs("Ele_pt.pdf")
    """
