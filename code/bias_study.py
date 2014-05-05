#!/usr/bin/env python

# set up stuff - check environment variables are set properly etc
import os
import sys
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i","--input",default="dat/bias_study_config.dat",help="cfg file")
parser.add_option("-t","--toysOnly",default=False,action="store_true")
parser.add_option("-e","--envelopeOnly",default=False,action="store_true")
(options,args) = parser.parse_args()

hostname = os.environ.get('HOST')
rootsys = os.environ.get('ROOTSYS')
if rootsys is None:
	print 'root not setup'
	if 'lxplus' in hostname:
		sys.exit('Run setup_root_lxplus.sh and try again')
	else:
		sys.exit()

import ROOT as r
import array

# load python stuff from python directory 
from python.bias_study_config import *
from python.config import *

cfg = config(options.input)
if cfg.suppressRooFitOutput:
	r.RooMsgService.instance().setGlobalKillBelow(r.RooFit.FATAL)
	r.RooMsgService.instance().setSilentMode(True)

# check c++ compiled library exists
while not os.path.exists('lib/libEnvelopeCode.so'):
	print 'Can\'t find the relevant library: lib/libEnvelopeCode.so'
	print 'Attempting build:'
	if os.system('make')!=0:
		sys.exit('Fail')

r.gROOT.ProcessLine('.L lib/libEnvelopeCode.so')
r.gROOT.SetBatch(cfg.batchmode)
os.system('mkdir -p diagnostics')
r.gROOT.ProcessLine(".x paperStyle.C")

if not options.envelopeOnly:

	infile = r.TFile(cfg.infile_name)
	ws = infile.Get(cfg.ws_name)
	mgg = ws.var(cfg.inv_mass_name)
	dataSet = ws.data(cfg.fit_gen_to_data_first)

	outfile = r.TFile(cfg.outfile_name,"RECREATE")

	# define poi
	poiLow = cfg.poi_scan_range[0] 
	poiHigh = cfg.poi_scan_range[1]
	mu = r.RooRealVar('mu','mu',0.,cfg.poi_range[0],cfg.poi_range[1])
	# build signal model
	mean = r.RooRealVar('mean','mean',125)
	mean.setConstant()
	sigma = r.RooRealVar('sigma','sigma',1.19)
	sigma.setConstant()
	nsignal_const = r.RooRealVar('nsignal_const','nsignal_const',50.8)
	nsignal_const.setConstant()
	sig_pdf = r.RooGaussian('gaus','gaus',mgg,mean,sigma)

	# setup gen pdf
	toySetup = r.PdfModelBuilder()
	toySetup.setObsVar(mgg)
	toySetup.setSignalModifier(mu)
	toySetup.addBkgPdf(cfg.genpdf_name+','+cfg.ws_name+','+cfg.infile_name,False)
	toySetup.setSignalPdf(sig_pdf,nsignal_const)
	toySetup.makeSBPdfs(True)
	if cfg.fit_gen_to_data_first!='' and cfg.fit_gen_to_data_first!='false' and cfg.fit_gen_to_data_first!='False':
		toySetup.setSignalModifierVal(cfg.gen_inj_sig)
		toySetup.setSignalModifierConstant(True)
		toySetup.fitToData(dataSet,False,True,True)
		toySetup.plotPdfsToData(dataSet,80,'diagnostics/genPdfFit',False)
		toySetup.setSignalModifierConstant(False)
	toySetup.setSeed(cfg.seed)

	# setup envelope
	envelopeSetup = r.PdfModelBuilder()
	envelopeSetup.setObsVar(mgg)
	envelopeSetup.setSignalModifier(mu)
	for pdf_name in cfg.envpdf_names:
		envelopeSetup.addBkgPdf(pdf_name+','+cfg.ws_name+','+cfg.infile_name,False)
	envelopeSetup.setSignalPdf(sig_pdf,nsignal_const)
	envelopeSetup.makeSBPdfs(False)

	swToy = r.TStopwatch()
	swToy.Start()
	sw = r.TStopwatch()
	for toy in range(cfg.ntoys):
		sw.Reset()
		sw.Start()
		toySetup.throwToy('truth_toy%d'%toy,int(dataSet.sumEntries()),False,True,True,True)
		toyData = toySetup.getToyDataSingle()
		profiler = r.ProfileMultiplePdfs()
		profiler.addPdfs(envelopeSetup.getSBPdfs())
		profiler.makeProfiles(toyData,mu,poiLow,poiHigh,cfg.points_in_scan,'_toy%d'%toy,cfg.printScanProgress)
		profiler.printProfiles()
		profiler.saveProfiles(outfile)
		del profiler
		sw.Stop()
		print '-----------------------------------------------------'
		print 'Toy took: Real time %3.2f (secs), CPU time %3.2f (secs)'%(sw.RealTime(),sw.CpuTime())
		print '-----------------------------------------------------'

	print 'Output profile curves written to', outfile.GetName()
	outfile.Close()
	swToy.Stop()
	print 'TOTAL FOR ALL TOYS TOOK: Real time %3.2f (secs), CPU time %3.2f (secs)'%(swToy.RealTime(),swToy.CpuTime())
	print '----------------------------------------------------------'

# this can compute bias given toys
if not options.toysOnly:
	
	# setup tree to save results
	infile = r.TFile(cfg.outfile_name)
	outfile = r.TFile(cfg.outfile_name.replace('.root','_hists.root'),'RECREATE')
	tree = r.TTree("BiasResults","BiasResults")
	outfile.mkdir("ProfileCurves")
	outfile.mkdir("ProfileEnvelopes")
	outfile.mkdir("Canvases")
	gen_val = array.array('d', [0])
	fit_val = array.array('d', [0])
	gen_pdf = r.TString() 
	fit_pdf = r.TString() 
	errUp_05sig = array.array('d', [0])
	errDn_05sig = array.array('d', [0])
	errUp_1sig = array.array('d', [0])
	errDn_1sig = array.array('d', [0])
	errUp_2sig = array.array('d', [0])
	errDn_2sig = array.array('d', [0])
	errUp_3sig = array.array('d', [0])
	errDn_3sig = array.array('d', [0])
	
	tree.Branch('gen_val',gen_val,'gen_val/D')
	tree.Branch('fit_val',fit_val,'fit_val/D')
	tree.Branch('gen_pdf','TString',r.AddressOf(gen_pdf))
	tree.Branch('fit_pdf','TString',r.AddressOf(fit_pdf))
	tree.Branch('errUp_05sig',errUp_05sig,'errUp_05sig/D')
	tree.Branch('errDn_05sig',errDn_05sig,'errDn_05sig/D')
	tree.Branch('errUp_1sig',errUp_1sig,'errUp_1sig/D')
	tree.Branch('errDn_1sig',errDn_1sig,'errDn_1sig/D')
	tree.Branch('errUp_2sig',errUp_2sig,'errUp_2sig/D')
	tree.Branch('errDn_2sig',errDn_2sig,'errDn_2sig/D')
	tree.Branch('errUp_3sig',errUp_3sig,'errUp_3sig/D')
	tree.Branch('errDn_3sig',errDn_3sig,'errDn_3sig/D')

	pull_hist = r.TH1F('pull_hist','',50,cfg.gen_inj_sig-3,cfg.gen_inj_sig+3)
	pull_hist.GetXaxis().SetTitle("Pull (#mu)")
	pull_hist.GetYaxis().SetTitle("Number of toy experiments")

	whichpdf_hist = r.TH1F('pdf_choice_hist','',len(cfg.envpdf_names),0,len(cfg.envpdf_names))
	whichpdf_hist.GetYaxis().SetTitle("Number of toy experiments")
	whichpdf_hist.GetXaxis().SetTitle("Best fit pdf")
	whichpdf_dict = {}
	for i, pdf in enumerate(cfg.envpdf_names):
		pdf_name = pdf.split('_')[-1]
		whichpdf_dict[pdf_name] = i
		whichpdf_hist.GetXaxis().SetBinLabel(i+1,pdf_name)
	whichpdf_hist.GetXaxis().SetLabelSize(whichpdf_hist.GetYaxis().GetLabelSize())
	
	# get graphs from current file and compute envelope best fit and error
	list_of_keys = [] 
	for key in infile.GetListOfKeys():
		list_of_keys.append(key.GetName())
	for t in range(cfg.ntoys):
		graphs = []
		for pdf_name in cfg.envpdf_names:
			for key_name in list_of_keys:
				if pdf_name in key_name and 'toy%d'%t in key_name:
					nparams = int(key_name.split('env_pdf_')[1].split('_')[0])
					#print key_name, ' -- ', nparams
					graphs.append([infile.Get(key_name),nparams])
		
		profiler = r.ProfileMultiplePdfs()
		for graph, npar in graphs:
			outfile.cd("ProfileCurves")
			graph.Write()
			profiler.addProfile(graph,npar)

		#profiler.printProfilesWithNPars()
		corrVal = 2.
		profiler.setCorrPdof(corrVal);
		profiler.constructEnvelope('_toy%d'%(t))
		if cfg.plotsForEachToy: profiler.drawEnvelope("diagnostics/envelope_toy%d.pdf"%t,True);
		envelope = profiler.getEnvelopeGraph();
		outfile.cd("ProfileEnvelopes")
		envelope.Write()

		gen_val[0] = cfg.gen_inj_sig
		fit_val[0] = profiler.getEnvelopeBestFitValue()
		gen_pdf.Resize(0)
		gen_pdf.Append(genpdf_name) 
		fit_pdf.Resize(0)
		fit_pdf.Append(profiler.getEnvelopeBestFitName()) 
		errUp_05sig[0] = profiler.getEnvelopeErrorUp(0.5)
		errDn_05sig[0] = profiler.getEnvelopeErrorDn(0.5)
		errUp_1sig[0] = profiler.getEnvelopeErrorUp(1)
		errDn_1sig[0] = profiler.getEnvelopeErrorDn(1)
		errUp_2sig[0] = profiler.getEnvelopeErrorUp(2)
		errDn_2sig[0] = profiler.getEnvelopeErrorDn(2)
		errUp_3sig[0] = profiler.getEnvelopeErrorUp(3)
		errDn_3sig[0] = profiler.getEnvelopeErrorDn(3)
		tree.Fill()
		if (fit_val[0]-gen_val[0]) >= 0.:
			pull_hist.Fill((fit_val[0]-gen_val[0])/(errUp_1sig[0]-fit_val[0]))
		else:
			pull_hist.Fill((fit_val[0]-gen_val[0])/(fit_val[0]-errDn_1sig[0]))
		
		fit_pdf_ind = whichpdf_dict[str(fit_pdf).split('_toy')[0].split('_')[-1]]
		whichpdf_hist.Fill(fit_pdf_ind)

	outfile.cd()
	tree.Write()

	# make pull plot
	canv = r.TCanvas()
	pull_fit = r.TF1('pull_fit','gaus',cfg.gen_inj_sig-3,cfg.gen_inj_sig+3)
	pull_hist.Fit(pull_fit,"QN")
	pull_hist.SetLineColor(r.kBlue+2)
	pull_hist.SetLineWidth(2)
	pull_fit.SetLineColor(r.kRed)
	pull_fit.SetLineWidth(2)
	pull_hist.Draw("LEP")
	# line for hist mean
	mean = r.TLine()
	mean.SetLineWidth(2)
	mean.SetLineColor(r.kBlue+2)
	# box for hist mean +/- 0.14
	box = r.TBox()
	box.SetLineColor(r.kRed)
	box.SetFillColor(r.kRed)
	box.SetFillStyle(3002)
	# line for fit mean
	fit_mean = r.TArrow()
	fit_mean.SetLineColor(r.kRed)
	fit_mean.SetLineWidth(2)
	fit_mean.SetFillColor(r.kRed)
	# line for gen val
	atgen = r.TLine()
	atgen.SetLineWidth(2)
	atgen.SetLineStyle(r.kDashed)
	atgen.SetLineColor(r.kBlack)
	
	box.DrawBox(pull_hist.GetMean()-0.14,0.,pull_hist.GetMean()+0.14,pull_hist.GetBinContent(pull_hist.FindBin(pull_hist.GetMean())))
	fit_mean.DrawArrow(pull_fit.GetParameter(1),0.,pull_fit.GetParameter(1),0.25*pull_hist.GetBinContent(pull_hist.FindBin(pull_fit.GetParameter(1))),0.025,'<|')
	mean.DrawLine(pull_hist.GetMean(),0.,pull_hist.GetMean(),pull_hist.GetBinContent(pull_hist.FindBin(pull_hist.GetMean())))
	atgen.DrawLine(cfg.gen_inj_sig,0.,cfg.gen_inj_sig,pull_hist.GetBinContent(pull_hist.FindBin(pull_hist.GetMean())))
	pull_hist.Draw("LEPsame")
	pull_fit.Draw("Lsame")

	leg = r.TLegend(0.58,0.65,0.89,0.91)
	leg.SetFillColor(0)
	leg.SetLineWidth(1)
	leg.SetHeader("%d Entries"%pull_hist.GetEntries())
	leg.AddEntry(pull_hist,"#splitline{Hist mean  = %5.3f}{Hist sigma = %5.3f}"%(pull_hist.GetMean(),pull_hist.GetRMS()),"LEP")
	leg.AddEntry(pull_fit,"#splitline{Fit mean   = %5.3f}{Fit sigma  = %5.3f}"%(pull_hist.GetMean(),pull_hist.GetRMS()),"L")
	leg.Draw("same")
	
	canv.Print("diagnostics/pull_hist.pdf")
	outfile.cd()
	pull_hist.Write()
	outfile.cd("Canvases")
	canv.SetName("pull_hist")
	canv.Write()

	# make fit pdf plot for number of times refitted
	whichpdf_hist.SetFillColor(r.kBlue+1)
	whichpdf_hist.SetStats(0)
	whichpdf_hist.Draw("HIST")
	lat = r.TLatex()
	lat.SetNDC()
	gen_name = cfg.genpdf_name.split('_')[-1]
	lat.DrawLatex(0.16,0.93,"Generated with %s at #mu=%3.1f"%(gen_name,cfg.gen_inj_sig))
	lat.DrawLatex(0.75,0.93,"%d Entries"%whichpdf_hist.GetEntries())
	canv.Print("diagnostics/pdf_choice_hist.pdf")
	outfile.cd()
	whichpdf_hist.Write()
	outfile.cd("Canvases")
	canv.SetName("pdf_choice_hist")
	canv.Modified()
	canv.Update()
	canv.Write()
	outfile.Close()
