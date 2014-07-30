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
		sys.exit('Run source setup_root_lxplus.sh and try again')
	else:
		sys.exit()

import ROOT as r
import array

# load python stuff from python directory
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
	#mu = r.RooRealVar('mu','mu',0.,cfg.poi_range[0],cfg.poi_range[1])
	# build signal model
	mean = ws.var("mean")#r.RooRealVar('mean','mean',125)
	mean.setConstant()
	sigma = ws.var("sigma")#r.RooRealVar('sigma','sigma',1.19)
	sigma.setConstant()
	#nsignal_const = r.RooRealVar('nsignal_const','nsignal_const',50.8)
	#nsignal_const.setConstant()
	#sig_pdf = r.RooGaussian('gaus','gaus',mgg,mean,sigma)
	mu = ws.var('r')
	mu.setRange(cfg.poi_range[0],cfg.poi_range[1])
	nsignal_const = ws.var("nsignal_const")
	nsignal_const.setConstant()
	sig_pdf = ws.pdf('gaus')

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
		if cfg.plotsForEachToy: toySetup.plotToysWithPdfs('diagnostics/genPdf%s_toy%d'%(cfg.genpdf_name,toy),160,False)
		toyData = toySetup.getToyDataSingle()
		profiler = r.ProfileMultiplePdfs()
		#profiler.wspace = ws
	        profiler.setObsVar(mgg)
		profiler.setSavePVal(cfg.savePVal)
		profiler.addPdfs(envelopeSetup.getSBPdfs())
		profiler.makeProfiles(toyData,mu,poiLow,poiHigh,cfg.points_in_scan,'_toy%d'%toy,cfg.printScanProgress)
		#profiler.printProfiles()
		profiler.saveProfiles(outfile)
		del profiler
		sw.Stop()
		print '-----------------------------------------------------'
		print 'Toy %d/%d took: Real time %3.2f (secs), CPU time %3.2f (secs)'%(toy+1,cfg.ntoys,sw.RealTime(),sw.CpuTime())
		print '-----------------------------------------------------'

	print 'Output profile curves written to', outfile.GetName()
	outfile.Close()
	swToy.Stop()
	print 'TOTAL FOR ALL %d TOYS TOOK: Real time %3.2f (secs), CPU time %3.2f (secs)'%(cfg.ntoys,swToy.RealTime(),swToy.CpuTime())
	print '----------------------------------------------------------'
	print 'Toy throwing done'
	print '----------------------------------------------------------'

# this can compute bias given toys
#if not options.toysOnly:

	# load plugin
	#from python.computeEnvelope import *
	#envelopeComputation(cfg.outfile_name.replace('.root','Summary.root'),[cfg.outfile_name],cfg.corrVals,cfg.genpdf_name,cfg.gen_inj_sig)
