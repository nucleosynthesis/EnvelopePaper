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

def trimers(gr,exp): # also reverses the points (1-coverage)
  np = gr.GetN()
  x = ROOT.Double(0);
  y = ROOT.Double(0);
  for n in range(np):
   gr.GetPoint(n,x,y);
   gr.SetPointError(n,gr.GetErrorX(n)/2,gr.GetErrorY(n)/exp);
   gr.SetPoint(n,x,(1-y)/exp);

def getError(gr,x):
  np = gr.GetN()
  arX = array.array('d',[0])
  arY = array.array('d',[0])
  for n in range(np): 
	gr.GetPoint(n,arX,arY)
	if (abs(x-arX[0]) < 0.001) :
		return gr.GetErrorX(n), gr.GetErrorY(n)
#fil = ROOT.TFile.Open("biastoys.root")
#fil = ROOT.TFile.Open("BiasNarrowRangeSummary.root")
#fil = ROOT.TFile.Open("BiasFirstOrdersSummary.root")
fil = ROOT.TFile.Open("BiasAllOrdersSummary_morecorrections.root")

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
styles = [20,20,20,20]
colors = [ROOT.kBlack,ROOT.kBlue,ROOT.kMagenta+1,ROOT.kGreen+2]
muvals = ["-0.75","0.","1.","1.5"]
corrections = ["0.25","0.5","0.75","1.","1.25","1.5","1.75","2.","2.25","2.5","2.75","3."]#numpy.arange(0.25,3.25,0.25)
names  = ["#mu=%s"%mv for mv in muvals]#["approx. p-value","p-value","Akaike"]
fits = ["envelope"]
#styles = [23]
#colors = [ROOT.kBlack]

leg = ROOT.TLegend(0.2,0.70,0.49,0.89)
leg.SetNColumns(2)
leg.SetFillColor(0)

c = ROOT.TCanvas("c","c",600,800)
pads = []
np = len(gens)
hists = []
dh = ROOT.TH1F("hd%s","hd",1,0,3.2);
#dh.GetYaxis().SetRangeUser(-.99,.99);
dh.GetYaxis().SetRangeUser(0.81,1.19);
#dh.GetYaxis().SetTitle("< (#mu - #hat{#mu})/#sigma >");
#dh.GetXaxis().SetTitle("#mu");
dy = 0.04
width = (1./np) - dy/np

line = ROOT.TLine(dh.GetXaxis().GetXmin(),1.,dh.GetXaxis().GetXmax(),1.)
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
 grs = []
 for i,g in enumerate(muvals):
  ng = ROOT.TGraphErrors()
# hists.append(dh)
  f = "envelope"
  for j,cVal in enumerate(corrections):
   #print g,f
   #gr = fil.Get("pull_mean_gen%sfit%s_c%s"%("bestfit_",f,cVal))
   #gr = fil.Get("pull_mean_gen%sfit%s_c%s"%(gen,f,cVal))
   gr = fil.Get("pull_mean_gen%sfit%s_c%s_cov%s"%(gen,f,cVal,"1."))
   bval = gr.Eval(float(g))
   ng.SetPoint(j,float(cVal),bval) 
   ex,ey = getError(gr,float(g))
   print gen, cVal,g,bval,ey
   ng.SetPointError(j,0.05,ey) 
  ng.SetMarkerStyle(styles[i])
  ng.SetMarkerColor(colors[i])
  ng.SetLineColor(colors[i])
  ng.SetMarkerSize(0.9)
  cv = 1. 
  prob = (1 -  ROOT.Math.chisquared_cdf_c(float(cv)*float(cv),1))
  trimers(ng,prob)
  grs.append(ng.Clone())
  allgrs.append(grs)

 print "Drawing", len(pads)
 #dh.GetXaxis().SetTitle("c")
 #dh.GetYaxis().SetTitle("< (#hat{#mu} - #mu)/#sigma >")
 dh.Draw("axis")
 line.Draw()
 for i,gg in enumerate(grs): 
  if m==0:leg.AddEntry(gg,names[i],"P")
  gg.Draw("ePsame")
 latlab.DrawLatex(0.935,0.15,gnames[m])
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
lat.DrawLatex(0.065,0.22,"Coverage/Expected Coverage (68.3% Interval)")
c.SaveAs("../correction/AllOrderFunctions_coverage_vs_correction.pdf")
#raw_input()
