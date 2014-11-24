#!/usr/bin/env python

import ROOT as r
import math
import numpy
import time
import sys
import array
#r.gROOT.LoadClass('ProfileMultiplePdfs','lib/libEnvelopeCode.so')
r.gROOT.ProcessLine('.L lib/libEnvelopeCode.so')
r.gROOT.ProcessLine('.x paperStyle.C')
infile = r.TFile("allorderfits.root")

ptypes = ['bern','exp','pow','lau']
pnames = {'bern':'Polynomial', 'exp':'Exponential Sum', 'pow':'Power Law Sum', 'lau':'Laurent Series'}

graphs = []

corr_vals = []
corr_vals = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.25,2.5]
corr_vals = corr_vals + ['P']

r.gROOT.SetBatch()

if len(sys.argv)<=1:
	print 'Hi Nick and Paul. You need to run nllScan_allOrders.C first and make sure the correction value is set to 0 in that script'
	print 'Then you can run this bad boy for various different correction schemes'
	print 'Corrections being run at the moment are:'
	print '\t', corr_vals
	raw_input('Are you happy to go ahead? Press return')
	print '\n\n Segmentation Violation ! \n\n'
	time.sleep(2)
	print 'Just kidding you bus. I\'ll carry on now'
	time.sleep(1)

##
## setup
##
epsilon = 0.01 # bin width offset
bestFitH = r.TH1F("bestFitH","",len(corr_vals),0,len(corr_vals))
bestFitH.SetStats(0)
bestFitH.SetLineWidth(3)
for p in range(len(corr_vals)):
	name = "%s / dof"%str(corr_vals[p])
	if (p==len(corr_vals)-1): name = "p-value"
	if (corr_vals[p]==0.0): name += " = No Correction"
	if (corr_vals[p]==1.0): name += " = Approx. p-value"
	if (corr_vals[p]==2.0): name += " = Akaike"
	bestFitH.GetXaxis().SetBinLabel(p+1,name)

bestFit = r.TGraphAsymmErrors()
bestFit.SetName("bestFit")
bestFit.SetLineColor(r.kRed)
bestFit.SetLineWidth(4)
bestFit.SetMarkerColor(r.kRed)
bestFit.SetMarkerSize(0)

err1sigma = r.TGraphAsymmErrors()
err1sigma.SetName("err1sigma")
err1sigma.SetLineColor(r.kGreen+3);
err1sigma.SetFillColor(r.kGreen+3);

err2sigma = r.TGraphAsymmErrors();
err2sigma.SetName("err2sigma");
err2sigma.SetLineColor(r.kYellow);
err2sigma.SetFillColor(r.kYellow);

def getLabel(corr):
	textStr = "c = %s"%str(corr)
	if corr == "P":
		textStr = "p-value"
	if corr == 0.0:
		textStr = "No Correction"
	if corr == 1.0:
		textStr = "Approx. p-value"
	if corr == 2.0:
		textStr = "Akaike"
	return textStr

def drawLabel(canv,corr):
	text = r.TText()
	text.SetNDC()
	text.SetTextSize(0.045)
	text.SetTextAlign(22) # right align
	textStr = getLabel(corr)
	text.DrawTextNDC(0.42,0.85,textStr)

def getName(pdf,nparams):
	return pnames[pdf]+" (%dpars)"%nparams

def getCorrectedGraph(gr,corr,gr_corr,nparams):
	gr_corrected = r.TGraph()
	if corr=='P':
		assert(gr.GetN()==gr_corr.GetN())
		x = r.Double()
		y = r.Double()
		x_corr = r.Double()
		y_corr = r.Double()
		for p in range(gr.GetN()):
			gr.GetPoint(p,x,y)
			gr_corr.GetPoint(p,x_corr,y_corr)
			assert(r.TMath.Abs(x_corr-x)<1.e-3)
			gr_corrected.SetPoint(p,x,y+y_corr)
	else:
		x = r.Double()
		y = r.Double()
		for p in range(gr.GetN()):
			gr.GetPoint(p,x,y)
			gr_corrected.SetPoint(p,x,y+(corr*nparams))

	return gr_corrected

def getSmoothedGraph(gr):

	# get graph min and assume it isn't a spike
	min_p = -100
	min_y = 1.e6
	x = r.Double()
	y = r.Double()
	for p in range(gr.GetN()):
		gr.GetPoint(p,x,y)
		if float(y) < min_y:
			min_y = float(y)
			min_p = p

	smooth_gr = r.TGraph()
	# now loop again trying to remove spikes
	# but have to do it in two halves
	y_prev = 0.
	new_p = 0

	# lower half (below best fit first)
	for p in range(min_p+1):
		gr.GetPoint(p,x,y)
		if p==0 or float(y)<y_prev:
			smooth_gr.SetPoint(new_p,x,y)
			new_p += 1
			y_prev = float(y)
		else:
				continue

	new_points = []
	# now the upper half (counting backwards)
	for p in range(gr.GetN()-1,min_p,-1):
		gr.GetPoint(p,x,y)
		if p==gr.GetN()-1 or float(y)<y_prev:
			new_points.append([float(x),float(y)])
			y_prev = float(y)

	# now add the upper half points to the graph
	new_points.sort(key=lambda x: x[0])
	for xn, yn in new_points:
		smooth_gr.SetPoint(new_p,xn,yn)
		new_p += 1

	return smooth_gr

def getCorrectedSmoothedGraph(gr,corr,gr_corr,nparams):

	gr_corrected = getCorrectedGraph(gr,corr,gr_corr,nparams)
	gr_smoothed = getSmoothedGraph(gr_corrected)
	return gr_smoothed

def getBestFitGraph(graphs,best_fit_name):
	for name, graph in graphs:
		if graph.GetName()==best_fit_name:
			return graph
	return None

# __main__

dummyHist = r.TH1F("dummy","",1,-1.,2.6)
dummyHist.SetStats(0)
dummyHist.GetXaxis().SetTitle("#mu")
dummyHist.GetYaxis().SetTitle("#Lambda + correction")
color = [r.kBlue, r.kRed, r.kMagenta, r.kGreen+1, r.kOrange+7, r.kAzure+10,r.kBlue]

for p, corr in enumerate(corr_vals):
	pdfit_c=0
	color = [r.kBlue, r.kRed, r.kMagenta, r.kGreen+1, r.kOrange+7, r.kAzure+10,r.kBlue]
	style=1
	profiler = r.ProfileMultiplePdfs()
	leg = r.TLegend(0.55,0.5,0.89,0.91)
	leg.SetFillColor(0)
	graphs = []
	for pdf in ptypes:
		for order in range(1,6):
			nparams = order+1
			if (pdf=='exp' or pdf=='pow') and order%2==0: continue

			gr = infile.Get("env_pdf_1_8TeV_%s%d_gr"%(pdf,order))
			gr_corr = infile.Get("env_pdf_1_8TeV_%s%d_gr_corr"%(pdf,order))

			#gr_corrected = getCorrectedGraph(gr,corr,gr_corr,nparams)
			gr_corrected = getCorrectedSmoothedGraph(gr,corr,gr_corr,nparams)
			gr_corrected.SetName(gr.GetName())

			if pdfit_c!=0 and pdfit_c%6==0:
				style += 1
			gr_corrected.SetLineColor(color[pdfit_c%7])
			gr_corrected.SetLineColor(style)
			pdfit_c += 1

			graphs.append([getName(pdf,nparams),gr_corrected])

			if corr == 'P':
				profiler.addProfile(gr, gr_corr)
			else:
				profiler.addProfile(gr, corr*nparams)

	profiler.constructEnvelope('envelope_c%s'%str(corr))

	# get best fit info
	envelope = profiler.getEnvelopeGraph()
	fit_val = profiler.getEnvelopeBestFitValue()
	fit_pdf = profiler.getEnvelopeBestFitName()
	minnll = profiler.getEnvelopeBestFitNll()
	errHigh1 = profiler.getEnvelopeErrorUp(1.)-fit_val
	errLow1 = fit_val - profiler.getEnvelopeErrorDn(1.)
	errHigh2 = profiler.getEnvelopeErrorUp(2.)-fit_val
	errLow2 = fit_val - profiler.getEnvelopeErrorDn(2.)
	bestFitGraph = getBestFitGraph(graphs,fit_pdf)

	##
	## profiles plot
	##
	canv = r.TCanvas()
	dummyHist.GetYaxis().SetRangeUser(r.TMath.Floor(minnll),r.TMath.Floor(minnll)+16)
	dummyHist.Draw("AXIS")
	style = 1
	for i, graph in enumerate(graphs):
		name = graph[0]
		gr = graph[1]
		if i!=0 and i%6==0: style += 1
		gr.SetLineColor(color[i%7])
		gr.SetLineStyle(style)
		gr.SetLineWidth(2)
		leg.AddEntry(gr,name,"L")
		gr.Draw("Lsame")
	leg.Draw("same")
	pave = r.TPave(0.24,0.84,0.5,0.9,4,"ndc")
	pave.SetShadowColor(0)
	pave.SetLineColor(0)
	pave.SetFillColor(0)
	#pave.Draw()
	drawLabel(canv,corr)
	canv.Update()
	canv.Modified()
	canv.Print("../correction/ProfilesAllOrders%s.pdf"%str(corr))

	##
	## envelope plot
	##
	bestFitGraph.SetLineWidth(2)
	bestFitGraph.SetLineColor(r.kRed)
	bestFitGraph.SetLineStyle(r.kDashed)
	envelope.SetLineColor(r.kBlack)
	envelope.SetLineWidth(2)
	# dashed lines
	line = r.TLine()
	line.SetLineColor(1)
	line.SetLineStyle(2)
	line.SetLineWidth(2)
	# error bands
	gr_err1 = r.TGraphAsymmErrors()
	gr_err2 = r.TGraphAsymmErrors()
	lowBound1 = math.ceil((fit_val-errLow1)/0.01)*0.01
	highBound1 = math.ceil((fit_val+errHigh1)/0.01)*0.01
	lowBound2 = math.ceil((fit_val-errLow2)/0.01)*0.01
	highBound2 = math.ceil((fit_val+errHigh2)/0.01)*0.01

	err1points = numpy.array(numpy.arange(lowBound1,highBound1,0.01))
	err2points = numpy.array(numpy.arange(lowBound2,highBound2,0.01))
	err1points = numpy.insert(err1points,0,fit_val-errLow1)
	err2points = numpy.insert(err2points,0,fit_val-errLow2)
	err1points = numpy.append(err1points,fit_val+errHigh1)
	err2points = numpy.append(err2points,fit_val+errHigh2)

	for np, xp in enumerate(err1points):
		gr_err1.SetPoint(np,xp,minnll)
		gr_err1.SetPointError(np,0.,0.,0.,envelope.Eval(xp)-minnll)
	for np, xp in enumerate(err2points):
		gr_err2.SetPoint(np,xp,minnll)
		gr_err2.SetPointError(np,0.,0.,0.,envelope.Eval(xp)-minnll)

	gr_err1.SetFillColor(r.kGreen+3)
	gr_err1.SetLineColor(r.kGreen+3)
	gr_err1.SetMarkerColor(r.kGreen+3)
	gr_err2.SetFillColor(r.kYellow)
	gr_err2.SetLineColor(r.kYellow)
	gr_err2.SetMarkerColor(r.kYellow)

	canv.Clear()
	dummyHist.Draw("AXIS")
	bestFitGraph.Draw("Lsame")
	envelope.Draw("Lsame")
	gr_err2.Draw("E3same")
	gr_err1.Draw("E3same")
	line.DrawLine(dummyHist.GetBinLowEdge(1),minnll,dummyHist.GetBinLowEdge(dummyHist.GetNbinsX()+1),minnll)
	line.DrawLine(dummyHist.GetBinLowEdge(1),minnll+1,dummyHist.GetBinLowEdge(dummyHist.GetNbinsX()+1),minnll+1)
	line.DrawLine(dummyHist.GetBinLowEdge(1),minnll+4,dummyHist.GetBinLowEdge(dummyHist.GetNbinsX()+1),minnll+4)
	envelope.Draw("Lsame")

	leg = r.TLegend(0.3,0.62,0.89,0.89)
	leg.AddEntry(envelope,"Minimum Envelope (%s)"%getLabel(corr),"L")
	leg.AddEntry(gr_err1,"68.3% Interval")
	leg.AddEntry(gr_err2,"95.4% Interval")
	leg.Draw()
	canv.Update()
	canv.Modified()
	canv.Print("../correction/EnvelopeAllOrders%s.pdf"%str(corr))


	# set values
	bestFitH.SetBinContent(p+1,fit_val)
	bestFit.SetPoint(p,bestFitH.GetBinCenter(p+1),fit_val)
	err1sigma.SetPoint(p,bestFitH.GetBinCenter(p+1),fit_val);
	err2sigma.SetPoint(p,bestFitH.GetBinCenter(p+1),fit_val);
	# set errors
	bestFitH.SetBinError(p+1,0.)
	bestFit.SetPointError(p,bestFitH.GetBinWidth(p+1)/2.-epsilon,bestFitH.GetBinWidth(p+1)/2.-epsilon,0.,0.)
	err1sigma.SetPointError(p,bestFitH.GetBinWidth(p+1)/2.-epsilon,bestFitH.GetBinWidth(p+1)/2.-epsilon,errLow1,errHigh1)
	err2sigma.SetPointError(p,bestFitH.GetBinWidth(p+1)/2.-epsilon,bestFitH.GetBinWidth(p+1)/2.-epsilon,errLow2,errHigh2)

	print 'c =', corr, ' mu =', fit_val, ' +', errHigh1, ' -', errLow1, ' (1sigma) +', errHigh2, ' -', errLow2, ' (2sigma)'

# now do correction plot stuff
canv = r.TCanvas("c","c",100*len(corr_vals),700)
canv.SetBottomMargin(0.25)
canv.SetLeftMargin(0.1)

leg = r.TLegend(0.15,0.7,0.5,0.89)
leg.SetFillColor(0)
leg.SetLineColor(0)
leg.AddEntry(bestFit,"Fit value", "L")
leg.AddEntry(err1sigma,"68.3% Interval", "LF")
leg.AddEntry(err2sigma,"95.4% Interval", "LF")

bestFitH.GetYaxis().SetRangeUser(-1.,4.)
bestFitH.GetYaxis().SetTitle("#mu")
bestFitH.GetYaxis().SetTitleSize(0.08)
bestFitH.GetYaxis().SetTitleOffset(0.4)
bestFitH.GetYaxis().SetLabelSize(0.06)

bestFitH.GetXaxis().LabelsOption("d");
bestFitH.GetXaxis().SetLabelOffset(0.01);
bestFitH.GetXaxis().SetLabelSize(0.07);
bestFitH.GetXaxis().SetTitle("#Lambda correction");
bestFitH.GetXaxis().SetTitleSize(0.08);
bestFitH.GetXaxis().SetTitleOffset(1.4);

line = r.TLine()
line.SetLineWidth(3)
line.SetLineStyle(7)
line.SetLineColor(r.kBlack)

bestFitH.Draw("AXIS")
err2sigma.Draw("2same")
err1sigma.Draw("2same")
bestFit.Draw("Esame")
line.DrawLine(len(corr_vals)-1,bestFitH.GetMinimum(),len(corr_vals)-1,bestFitH.GetMaximum())
leg.Draw()
canv.RedrawAxis()
canv.Update()
canv.Modified()
canv.Print("../correction/correction.pdf")
