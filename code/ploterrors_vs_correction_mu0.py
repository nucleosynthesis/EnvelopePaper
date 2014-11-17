import sys
import numpy,array
cVal = '0.'
if len(sys.argv)==2:
	cVal = sys.argv[1]

import ROOT
ROOT.gROOT.ProcessLine(".x paperStyle.C");

ROOT.gStyle.SetTextSize(0.056);
ROOT.gStyle.SetTitleSize(0.12,"XYZ");
ROOT.gStyle.SetLabelSize(0.09, "XYZ");
ROOT.gStyle.SetPadBorderMode(0);
ROOT.gStyle.SetTitleYOffset(0.5);

ROOT.gROOT.SetBatch(1)

def trimers(gr):
  np = gr.GetN()
  for n in range(np):
   gr.SetPointError(n,gr.GetErrorX(n)/2,gr.GetErrorY(n));
def getError(gr,x):
  np = gr.GetN()
  arX = array.array('d',[0])
  arY = array.array('d',[0])
  for n in range(np): 
	gr.GetPoint(n,arX,arY)
	if (abs(x-arX[0]) < 0.001) :
		return gr.GetErrorX(n), gr.GetErrorY(n)
def getAsymmError(gr,x):
  np = gr.GetN()
  arX = array.array('d',[0])
  arY = array.array('d',[0])
  for n in range(np): 
	gr.GetPoint(n,arX,arY)
	if (abs(x-arX[0]) < 0.001) :
		return gr.GetErrorYlow(n), gr.GetErrorYhigh(n)
def MattIsBadAtStats(gr):
  np = gr.GetN()
  arX = array.array('d',[0])
  arY = array.array('d',[0])
  for n in range(np):
   gr.GetPoint(n,arX,arY)
   yval = (gr.GetErrorY(n)**2 - arY[0]**2)**0.5
   gr.SetPointError(n,gr.GetErrorX(n),yval);
#fil = ROOT.TFile.Open("biastoys.root")
#fil = ROOT.TFile.Open("BiasNarrowRangeSummary.root")
#fil = ROOT.TFile.Open("BiasFirstOrdersSummary.root")
fil = ROOT.TFile.Open("BiasAllOrdersSummary_morecorrections_errs_v4.root")

gens = [
		#"bestfit_"
		"env_pdf_1_8TeV_exp1_"
		,"env_pdf_1_8TeV_pow1_"
		,"env_pdf_1_8TeV_lau1_"
		,"env_pdf_1_8TeV_bern4_"
		,"env_pdf_1_8TeV_bern5_"
	   ]
gnames = [
 	#"Profiled Function"
	"Exponential (2 pars)"
	,"Power Law (2 pars)"
	,"Laurent (2 pars)"
	,"Polynomial (5 pars)"
 	,"Polynomial (6 pars)"
	]
# these are the same sizes  (add back in poly)
#names = ["Laurent","Power Law","Polynomial","Exponential","Envelope"]
#fits = ["lau1","pow1","bern1","exp1","envelope"]
#styles = [20,24,21,25,23]
#colors = [ROOT.kGreen+2,ROOT.kBlue,ROOT.kMagenta,ROOT.kRed,ROOT.kBlack]
#names = ["Laurent","Power Law","Exponential","Envelope"]
#fits = ["lau1","pow1","exp1","envelope"]
styles = [23,25,28,23]
colors  = [ROOT.kBlue,ROOT.kMagenta+1,ROOT.kGreen+2]
fcolors = [ROOT.kAzure-4,ROOT.kMagenta+1,ROOT.kGreen+2]
fstyles = [1001,1001]
#muvals = ["-0.75","0.","1.","1.5"]
#muvals = ["1.","0."]
muvals = ["0."]
corrections = ["0.25","0.5","0.75","1.","1.25","1.5","1.75","2.","2.25","2.5","2.75","3."]#numpy.arange(0.25,3.25,0.25)
#names  = ["<#hat{#mu} - %s>"%mv for mv in muvals]#["approx. p-value","p-value","Akaike"]
names  = ["#mu = %s"%mv for mv in muvals]#["approx. p-value","p-value","Akaike"]
fits = ["envelope"]
#styles = [23]
#colors = [ROOT.kBlack]

leg = ROOT.TLegend(0.2,0.70,0.54,0.89)
leg.SetNColumns(2)
leg.SetFillColor(0)

c = ROOT.TCanvas("c","c",600,800)
pads = []
np = len(gens)
hists = []
dh = ROOT.TH1F("hd%s","hd",1,0,3.2);
dh.GetYaxis().SetRangeUser(-0.79,1.19);
lineDashed = ROOT.TLine(dh.GetXaxis().GetXmin(),0,dh.GetXaxis().GetXmax(),0)
lineDashed.SetLineColor(1)
lineDashed.SetLineStyle(2)
#dh.GetYaxis().SetTitle("< (#mu - #hat{#mu})/#sigma > or <68.3\% interval>");
#dh.GetXaxis().SetTitle("#mu");
dy = 0.04
width = (1./np) - dy/np

line = ROOT.TLine(dh.GetXaxis().GetXmin(),0,dh.GetXaxis().GetXmax(),0)
line.SetLineColor(1)
line.SetLineStyle(2)
allgrs = []
for m,gen in enumerate(gens):
 latlab = ROOT.TLatex()
 latlab.SetTextSize(0.09)
 latlab.SetNDC()
 latlab.SetTextAngle(90)
 p = ROOT.TPad(gen,gen,0,dy+float(m)*width,1,dy+float(m+1)*width)
 p.SetCanvas(c)
 if not m==np-1   : p.SetTopMargin(0.01)
 if not m==0      : p.SetBottomMargin(0.01)
# p.SetBorderSize(0);
 p.Draw()
 p.cd()
 pads.append(p)
 grs  = []
 grbs = []
 for kk,midtype in enumerate(["mean"]):
  for i,g in enumerate(muvals):
   ng  = ROOT.TGraphErrors()
   ngb = ROOT.TGraphAsymmErrors()
# hists.append(dh)
   f = "envelope"
   #f = "pow1"
   for j,cVal in enumerate(corrections):
   #print g,f
    gr_statL = fil.Get("mean_staterrlow_gen%sfit%s_c%s"%(gen,f,cVal))
    gr_statH = fil.Get("mean_staterrhigh_gen%sfit%s_c%s"%(gen,f,cVal))
    gr      = fil.Get("mean_resid_gen%sfit%s_c%s"%(gen,f,cVal))
    bval = gr.Eval(float(g))
    #if kk==0 : delta = -0.02
    #if kk==1 : delta = +0.02
    delta=0
    ng.SetPoint(j,float(cVal)+delta,bval) 
    ex,ey = getError(gr,float(g))
    #estat_low,estat_high = getAsymmError(gr_stat,float(g))
    print gen, cVal,g,bval,ey
    ng.SetPointError(j,0.05,ey) 
    ngb.SetPoint(j,float(cVal),0)
    print "Errors Vals",gr_statL.Eval(float(g)),gr_statH.Eval(float(g))
    ngb.SetPointError(j,0.125,0.125,gr_statL.Eval(float(g)),gr_statH.Eval(float(g)))
   ng.SetMarkerStyle(styles[i])
   ng.SetMarkerColor(colors[i])
   ng.SetLineColor(colors[i])
   ng.SetMarkerSize(0.9)
   ngb.SetFillColor(fcolors[i])
   ngb.SetFillStyle(fstyles[i])
   #trimers(ng)
   if kk == 0 : 
	ng.SetMarkerStyle(20)
   if kk == 1 : 
	ng.SetMarkerStyle(24)
	#MattIsBadAtStats(ng)
   grs.append(ng.Clone())
   grbs.append(ngb.Clone())
   allgrs.append(grs)
   allgrs.append(grbs)

 print "Drawing", len(pads)
 #dh.GetXaxis().SetTitle("c")
 #dh.GetYaxis().SetTitle("< (#hat{#mu} - #mu)/#sigma >")
 dh.Draw("axis")
 
 #line.Draw()
 for i,gg in enumerate(grbs): 
  gg.Draw("e2same")
 lineDashed.Draw() 
 for i,gg in enumerate(grs): 
  #if m==0:
  # if i< len(names): leg.AddEntry(gg,names[i],"P")
  gg.Draw("ePsame")
 latlab.DrawLatex(0.935,0.15,gnames[m])
 if m==0:
  for i,gg in enumerate(zip(grs,grbs)): 
   if i< len(names): 
   	leg.AddEntry(gg[0],names[i],"P")
	#leg.AddEntry(gg[1],"<68.3% Interval>","F")
 if m==np-1: leg.Draw()
 c.cd()
#raw_input()
#if i==0:
#lat.SetTextAngle(90)
lat = ROOT.TLatex()
lat.SetTextSize(0.05)
lat.SetNDC()
lat.DrawLatex(0.6,0.012,"Correction / par")
#lat.DrawLatex(0.2,0.88,"Correction = %s"%cVal)
lat.SetTextAngle(90)
lat.DrawLatex(0.06,0.45,"<#hat{#mu} - #mu> or <68.3% interval>")
c.SaveAs("../correction/AllOrderFunctions_errors_vs_correction_%s.pdf"%muvals[0])
#raw_input()
