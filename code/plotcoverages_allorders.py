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

def trimers(gr,exp): # also reverses the points (1-coverage)
  np = gr.GetN()
  x = ROOT.Double(0);
  y = ROOT.Double(0);
  for n in range(np):
   gr.GetPoint(n,x,y);
   gr.SetPointError(n,gr.GetErrorX(n)/2,gr.GetErrorY(n)/exp);
   gr.SetPoint(n,x,(1-y)/exp);

#fil = ROOT.TFile.Open("biastoys.root")
#fil = ROOT.TFile.Open("BiasNarrowRangeSummary.root")
#fil = ROOT.TFile.Open("BiasFirstOrdersSummary.root")
fil = ROOT.TFile.Open("BiasAllOrdersSummary.root")

gens = [
		"bestfit_"
		,"env_pdf_1_8TeV_exp1_"
		,"env_pdf_1_8TeV_pow1_"
		,"env_pdf_1_8TeV_bern1_"
		,"env_pdf_1_8TeV_lau1_"
		,"env_pdf_1_8TeV_bern4_"
		,"env_pdf_1_8TeV_bern5_"
	   ]
gnames = [
 	"Gen Profiled Function"
	, "Gen Exponential 1"
 	,"Gen Power Law 1"
 	,"Gen Polynomial 1"
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
#styles = [20,24,25,23]
#colors = [ROOT.kGreen+2,ROOT.kBlue,ROOT.kRed,ROOT.kBlack]
names = ["Envelope"]
fits = ["envelope"]
styles = [23]
colors = [ROOT.kBlack]

cvals = ["0.5","1.","2.","3."]
covCanvs = []

for i,c in enumerate(cvals):
  	can = ROOT.TCanvas("c%d"%i,"c%d"%i,600,1600)
  	covCanvs.append(can)
  	leg = ROOT.TLegend(0.72,0.52,0.89,0.89)
  	leg.SetFillColor(0)

  	pads = []
  	np = len(gens)
	hists = []
	dh = ROOT.TH1F("hd%s","hd",1,-1.,2.2);
	#dh.GetYaxis().SetRangeUser(0.99,1.01);
	if c=="0.5":
		dh.GetYaxis().SetRangeUser(0.8,1.15);
	elif c=="1.":
		dh.GetYaxis().SetRangeUser(0.9,1.05);
	elif c=="2.":
		dh.GetYaxis().SetRangeUser(0.95,1.05);
	elif c=="3.":
		dh.GetYaxis().SetRangeUser(0.98,1.02);

	#dh.GetYaxis().SetTitle("< (#mu - #hat{#mu})/#sigma >");
	#dh.GetXaxis().SetTitle("#mu");
	dy = 0.04
	width = (1./np) - dy/np


	#cvals = ["1.","2."]
	lines = []

	#for c in cvals:
	prob = (1 -  ROOT.Math.chisquared_cdf_c(float(c)*float(c),1))
	print prob
	line = ROOT.TLine(dh.GetXaxis().GetXmin(),1,dh.GetXaxis().GetXmax(),1)
	line.SetLineColor(1)
	line.SetLineStyle(2)
	lines.append(line)

	for i,g in enumerate(gens):
	 #make some pads?
	 latlab = ROOT.TLatex()
	 latlab.SetTextSize(0.08)
	 latlab.SetNDC()
	 latlab.SetTextAngle(90)
	 p = ROOT.TPad(g,g,0,dy+float(i)*width,1,dy+float(i+1)*width)
	 p.SetCanvas(can)
	 if not i==np-1   : p.SetTopMargin(0.01)
	 if not i==0      : p.SetBottomMargin(0.01)
	# p.SetBorderSize(0);
	 p.Draw()
	 p.cd()
	 pads.append(p)
	 dh.Draw("axis")
	 for line in lines: line.Draw()
	 latlab.DrawLatex(0.935,0.15,gnames[i])
	# hists.append(dh)
	 for j,f in enumerate(fits):
	   #print g,f
	   p.cd()
	   #p.SetLogy()
	   grs = []
	   #for c in cvals:
	   print "pull_mean_gen%sfit%s_c0._cov1%s"%(g,f,c)

	   gr = fil.Get("pull_mean_gen%sfit%s_c0._cov%s"%(g,f,c))
	   trimers(gr,prob)
	   gr.SetMarkerStyle(styles[j])
	   gr.SetMarkerColor(colors[j])
	   gr.SetLineColor(colors[j])
	   gr.SetMarkerSize(0.9)
	   gr.Draw("epsame")
	   grs.append(gr)
	   if i==0:
		leg.AddEntry(grs[0],names[j],"P")
	 if i==np-1: leg.Draw()
	 can.cd()

	#if i==0:
	lat = ROOT.TLatex()
	lat.SetNDC()
	lat.DrawLatex(0.85,0.012,"#mu")
	lat.DrawLatex(0.2,0.96,"%s #sigma Interval"%c)
	lat.DrawLatex(0.2,0.88,"Correction = %s"%cVal)
	lat.SetTextAngle(90)
	lat.DrawLatex(0.065,0.45,"Coverage/Expected Coverage")

	can.SaveAs("../functions/AllOrderFunctions_Coverage_%s_c%s.pdf"%(c,cVal))
	can.Update()
	can.Modified()
	#raw_input()
