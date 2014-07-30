import sys
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

def trimers(gr):
  np = gr.GetN()
  for n in range(np):
   gr.SetPointError(n,gr.GetErrorX(n)/2,gr.GetErrorY(n));

#fil = ROOT.TFile.Open("biastoys.root")
#fil = ROOT.TFile.Open("BiasNarrowRangeSummary.root")
#fil = ROOT.TFile.Open("BiasFirstOrdersSummary.root")
fil = ROOT.TFile.Open("BiasAllOrdersSummary.root")

gens = [
		"bestfit_"
		,"env_pdf_1_8TeV_exp1_"
		,"env_pdf_1_8TeV_pow1_"
		,"env_pdf_1_8TeV_lau1_"
		,"env_pdf_1_8TeV_bern4_"
		,"env_pdf_1_8TeV_bern5_"
	   ]
gnames = [
 	"Gen Profiled Function"
	, "Gen Exponential 1"
 	,"Gen Power Law 1"
	,"Gen Laurent 1"
 	,"Gen Polynomial 4"
 	,"Gen Polynomial 5"
	]
# these are the same sizes  (add back in poly)
#names = ["Laurent","Power Law","Polynomial","Exponential","Envelope"]
#fits = ["lau1","pow1","bern1","exp1","envelope"]
#styles = [20,24,21,25,23]
#colors = [ROOT.kGreen+2,ROOT.kBlue,ROOT.kMagenta,ROOT.kRed,ROOT.kBlack]
#names = ["Laurent","Power Law","Exponential","Envelope"]
#fits = ["lau1","pow1","exp1","envelope"]
styles = [20,24,25,23]
colors = [ROOT.kBlack,ROOT.kBlue,ROOT.kRed,ROOT.kGreen+2]
corrections = ["1.","P","2."]
names = ["c=1","c=p(#chi^{2})","c=2"]
fits = ["envelope"]
#styles = [23]
#colors = [ROOT.kBlack]

leg = ROOT.TLegend(0.72,0.52,0.89,0.89)
leg.SetFillColor(0)

c = ROOT.TCanvas("c","c",600,1600)
pads = []
np = len(gens)
hists = []
dh = ROOT.TH1F("hd%s","hd",1,-1.,2.2);
dh.GetYaxis().SetRangeUser(-0.69,0.69);
#dh.GetYaxis().SetTitle("< (#mu - #hat{#mu})/#sigma >");
#dh.GetXaxis().SetTitle("#mu");
dy = 0.04
width = (1./np) - dy/np

line = ROOT.TLine(dh.GetXaxis().GetXmin(),0,dh.GetXaxis().GetXmax(),0)
line.SetLineColor(1)
line.SetLineStyle(2)

for i,g in enumerate(gens):
 #make some pads?
 latlab = ROOT.TLatex()
 latlab.SetTextSize(0.08)
 latlab.SetNDC()
 latlab.SetTextAngle(90)
 p = ROOT.TPad(g,g,0,dy+float(i)*width,1,dy+float(i+1)*width)
 p.SetCanvas(c)
 if not i==np-1   : p.SetTopMargin(0.01)
 if not i==0      : p.SetBottomMargin(0.01)
# p.SetBorderSize(0);
 p.Draw()
 p.cd()
 pads.append(p)
 dh.Draw("axis")
 line.Draw()
 latlab.DrawLatex(0.935,0.15,gnames[i])
# hists.append(dh)
 f = "envelope"
 for j,cVal in enumerate(corrections):
   #print g,f
   p.cd()
   gr = fil.Get("pull_mean_gen%sfit%s_c%s"%(g,f,cVal))
   trimers(gr)
   gr.SetMarkerStyle(styles[j])
   gr.SetMarkerColor(colors[j])
   gr.SetLineColor(colors[j])
   gr.SetMarkerSize(0.9)
   gr.Draw("epsame")
   if i==0:
	leg.AddEntry(gr,names[j],"P")
 if i==np-1: leg.Draw()
 c.cd()

#if i==0:
lat = ROOT.TLatex()
lat.SetNDC()
lat.DrawLatex(0.85,0.012,"#mu")
#lat.DrawLatex(0.2,0.88,"Correction = %s"%cVal)
lat.SetTextAngle(90)
lat.DrawLatex(0.065,0.7,"< (#hat{#mu} - #mu)/#sigma >")
c.SaveAs("../functions/AllOrderFunctions_call.pdf")
#raw_input()
