#!/usr/bin/env python
# vim: ts=2 sw=2 noexpandtab

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f","--filename")
parser.add_option("-d","--datfile")
parser.add_option("-D","--outdir")
parser.add_option("-p","--makePlotsOnly",action="store_true",default=False)
parser.add_option("-b","--isBatch",action="store_true",default=False)
parser.add_option("-s","--splitJobs",action="store_true",default=False)
(options,args) = parser.parse_args()

import os
#os.system('bash setup_root_lxplus.sh')

import sys

import ROOT as r
r.gROOT.ProcessLine(".x paperStyle.C")

from python.biasComputer import *
if options.isBatch: r.gROOT.SetBatch()

cfg = {}
pullMeanGraphs={}
pullWidthGraphs={}
coverageGraphs={}
colors = [r.kRed,r.kBlue,r.kGreen+1,r.kMagenta]

def readConfig():
	f = open(options.datfile)
	for line in f.readlines():
		if line.startswith('#'): continue
		if '=' not in line: continue
		line = line.strip('\n')
		line = line.strip()
		option = line.split('=')[0]
		value = line.split('=')[1]
		cfg[option] = value
	f.close()

def remakeDictionaries():
	tf = r.TFile(options.filename)
	for gen_pdf in cfg['gen_pdfs'].split(','):
		pullMeanGraphs[gen_pdf] = {}
		pullWidthGraphs[gen_pdf] = {}
		coverageGraphs[gen_pdf] = {}
		for c in cfg['corrVals'].split(','):
			pullMeanGraphs[gen_pdf][c] = tf.Get('pull_mean_gen%s_c%s'%(gen_pdf,c))
			pullWidthGraphs[gen_pdf][c] = tf.Get('pull_width_gen%s_c%s'%(gen_pdf,c))
			coverageGraphs[gen_pdf][c] = {}
			for cov in [float(x) for x in cfg['coverageValues'].split(',')]:
				coverageGraphs[gen_pdf][c][cov] = tf.Get('cov%3.1f_gen%s_c%s'%(cov,gen_pdf,c))
		if 'savePVal' in cfg.keys() and (cfg['savePVal']=='True' or cfg['savePVal']=='1'):
			pullMeanGraphs[gen_pdf]['P'] = tf.Get('pull_mean_gen%s_cP'%(gen_pdf))
			pullWidthGraphs[gen_pdf]['P'] = tf.Get('pull_width_gen%s_cP'%(gen_pdf))
			coverageGraphs[gen_pdf]['P'] = {}
			for cov in [float(x) for x in cfg['coverageValues'].split(',')]:
				coverageGraphs[gen_pdf]['P'][cov] = tf.Get('cov%3.1f_gen%s_cP'%(cov,gen_pdf))

def makePullPlot(doMean=True):

	histDict = {}
	min = float(cfg['inj_mu_vals'].split(',')[0])-0.5
	max = float(cfg['inj_mu_vals'].split(',')[-1])+1.5
	dummyHist = r.TH1F('dummy','',1,min,max)
	dummyHist.SetStats(0)
	dummyHist.GetXaxis().SetTitle("Generated #mu")
	if doMean:
		histDict = pullMeanGraphs
		dummyHist.GetYaxis().SetTitle("Mean of pull distribution")
		dummyHist.GetYaxis().SetRangeUser(-0.1,0.1)
	else:
		histDict = pullWidthGraphs
		dummyHist.GetYaxis().SetTitle("Width of pull distribution")
		dummyHist.GetYaxis().SetRangeUser(0.9,1.1)

	for c in cfg['corrVals'].split(','):
		canv = r.TCanvas()
		dummyHist.Draw("AXIS")
		leg = r.TLegend(0.7,0.3,0.95,0.8)
		leg.SetFillColor(0)
		leg.SetLineColor(1)
		leg.SetBorderSize(1)
		for i, gen_pdf in enumerate(cfg['gen_pdfs'].split(',')):
			histDict[gen_pdf][c].SetLineColor(colors[i])
			histDict[gen_pdf][c].SetLineWidth(2)
			leg.AddEntry(histDict[gen_pdf][c],gen_pdf,'LEP')
			histDict[gen_pdf][c].Draw("EPsame")
		line = r.TLine()
		line.SetLineColor(r.kBlack)
		line.SetLineStyle(r.kDashed)
		if doMean:
			line.DrawLine(min,0.,max,0.)
			leg.Draw("same")
			canv.Print("results/pull_mean_c%s.pdf"%c)
		else:
			line.DrawLine(min,1.,max,1.)
			leg.Draw("same")
			canv.Print("results/pull_width_c%s.pdf"%c)

def makeCoveragePlot():

	min = float(cfg['inj_mu_vals'].split(',')[0])-0.5
	max = float(cfg['inj_mu_vals'].split(',')[-1])+1.5
	dummyHist = r.TH1F('dummy','',1,min,max)
	dummyHist.SetStats(0)
	dummyHist.GetXaxis().SetTitle("Generated #mu")
	dummyHist.GetYaxis().SetTitle("Coverage")
	dummyHist.GetYaxis().SetRangeUser(1.e-3,1.)

	for c in cfg['corrVals'].split(','):
		canv = r.TCanvas()
		dummyHist.Draw("AXIS")
		leg = r.TLegend(0.7,0.3,0.95,0.8)
		leg.SetFillColor(0)
		leg.SetLineColor(1)
		leg.SetBorderSize(1)
		for i, gen_pdf in enumerate(cfg['gen_pdfs'].split(',')):
			for j, cov in enumerate([float(x) for x in cfg['coverageValues'].split(',')]):
				coverageGraphs[gen_pdf][c][cov].SetLineColor(colors[i])
				coverageGraphs[gen_pdf][c][cov].SetLineWidth(2)
				coverageGraphs[gen_pdf][c][cov].Draw("EPsame")
				if j==0: leg.AddEntry(coverageGraphs[gen_pdf][c][cov],gen_pdf,'LEP')
				if i==0:
					line = r.TLine()
					line.SetLineColor(r.kBlack)
					line.SetLineStyle(r.kDashed)
					line.DrawLine(min,1.-r.TMath.Prob(cov*cov,1),max,1.-r.TMath.Prob(cov*cov,1))
		leg.Draw("same")
		canv.Print("results/coverage_c%s.pdf"%c)

def makeWhichPdf2DPlot():

	env_pdfs = cfg['env_pdfs'].split(',')
	gen_pdfs = cfg['gen_pdfs'].split(',')

	hist = r.TH2F('whichpdf_hist','',len(env_pdfs),0,len(env_pdfs),len(gen_pdfs),0,len(gen_pdfs))

	for x, env_pdf in enumerate(env_pdfs):
		hist.GetXaxis().SetBinLabel(x+1,env_pdf)
	for y, gen_pdf in enumerate(gen_pdfs):
		hist.GetYaxis().SetBinLabel(y+1,gen_pdf)

	for y, gen_pdf in enumerate(gen_pdfs):
		"get on parade"

def makePlots():
	os.system('mkdir -p results')
	makePullPlot(True)
	makePullPlot(False)
	makeCoveragePlot()
	makeWhichPdf2DPlot()

def postProcessToys(outfiledir,gen_pdf,mu_val,env_pdfs):
	# given a file should return pull histogram and which_pdf histogram and save it to a file
	# outfiledir = 'dir/envpdfs/gen_pdf/mu_val/name.root'
	os.system('mkdir -p %s'%outfiledir)
	outf = r.TFile('%s/%s'%(outfiledir,options.filename),'RECREATE')

	file_names = []
	for job in range(int(cfg['njobs'])):
		file_names.append('%s/outfiles/BiasResults_mu%4.2f_gen%s_job%d.root'%(cfg['store_directory'],mu_val,gen_pdf,job))

	biasComp = biasComputer()
	biasComp.setCoverageValues([float(x) for x in cfg['coverageValues'].split(',')])
	biasComp.setOutFile(outf)
	biasComp.addFiles(file_names)
	biasComp.setPullHistBinning(50)
	biasComp.setGenPdf(gen_pdf)
	biasComp.setMuVal(mu_val)
	biasComp.setEnvPdfs(env_pdfs)
	biasComp.doPlots = False

	for c in cfg['corrVals'].split(','):
		biasComp.setCorrection(c)
		biasComp.compute('mu%4.2f_gen%s_c%s'%(mu_val,gen_pdf,c))
		pullHist = biasComp.getPullHist()
		covHist = biasComp.getCoverageHist()
		whichHist = biasComp.getWhichPdfHist()
		outf.cd()
		pullHist.Write()
		covHist.Write()
		whichHist.Write()
	del biasComp
	print 'Computation of bias and coverage done for mu=%4.2f and gen=%s'%(mu_val,gen_pdf)

# MAIN BELOW:

readConfig()
compute_pdf_sets = []
for set in cfg['compute_pdf_sets'].split(':'):
	the_set = set.strip('[').strip(']').split(',')
	compute_pdf_sets.append(the_set)

for pdf_set in compute_pdf_sets:
	location_track = ''
	for pdf in pdf_set: location_track += pdf

	for gen_pdf in cfg['gen_pdfs'].split(','):
		for val in cfg['inj_mu_vals'].split(','):
			mu_val = float(val)
			output_loc = '%s/env_%s/genpdf_%s/mu_%4.2f'%(options.outdir,location_track,gen_pdf,mu_val)
			os.system('mkdir -p %s'%output_loc)
			postProcessToys(output_loc,gen_pdf,mu_val,pdf_set)

sys.exit()

"""
# for plot only option
if options.makePlotsOnly:
	remakeDictionaries()
	makePlots()
	sys.exit()



sys.exit()

# if doing the full garb
outfile = r.TFile(options.filename,'RECREATE')

for gen_pdf in cfg['gen_pdfs'].split(','):
	pullMeanGraphs[gen_pdf] = {}
	pullWidthGraphs[gen_pdf] = {}
	coverageGraphs[gen_pdf] = {}
	for c in cfg['corrVals'].split(','):
		pullMeanGraphs[gen_pdf][c] = r.TGraphErrors()
		pullMeanGraphs[gen_pdf][c].SetName('pull_mean_gen%s_c%s'%(gen_pdf,c))
		pullWidthGraphs[gen_pdf][c] = r.TGraphErrors()
		pullWidthGraphs[gen_pdf][c].SetName('pull_width_gen%s_c%s'%(gen_pdf,c))
		coverageGraphs[gen_pdf][c] = {}
		for cov in [float(x) for x in cfg['coverageValues'].split(',')]:
			coverageGraphs[gen_pdf][c][cov] = r.TGraphErrors()
			coverageGraphs[gen_pdf][c][cov].SetName('cov%3.1f_gen%s_c%s'%(cov,gen_pdf,c))
	if 'savePVal' in cfg.keys() and (cfg['savePVal']=='True' or cfg['savePVal']=='1'):
		pullMeanGraphs[gen_pdf]['P'] = r.TGraphErrors()
		pullMeanGraphs[gen_pdf]['P'].SetName('pull_mean_gen%s_cP'%(gen_pdf))
		pullWidthGraphs[gen_pdf]['P'] = r.TGraphErrors()
		pullWidthGraphs[gen_pdf]['P'].SetName('pull_width_gen%s_cP'%(gen_pdf))
		coverageGraphs[gen_pdf]['P'] = {}
		for cov in [float(x) for x in cfg['coverageValues'].split(',')]:
			coverageGraphs[gen_pdf]['P'][cov] = r.TGraphErrors()
			coverageGraphs[gen_pdf]['P'][cov].SetName('cov%3.1f_gen%s_cP'%(cov,gen_pdf))
	p = 0
	for val in cfg['inj_mu_vals'].split(','):
		mu_val = float(val)

		# get appropriate files
		file_names = []
		for job in range(int(cfg['njobs'])):
			file_names.append('%s/outfiles/BiasResults_mu%4.2f_gen%s_job%d.root'%(cfg['store_directory'],mu_val,gen_pdf,job))

		# do bias and coverage computation
		biasComp = biasComputer()
		biasComp.setCoverageValues([float(x) for x in cfg['coverageValues'].split(',')])
		biasComp.setOutFile(outfile)
		biasComp.addFiles(file_names)
		biasComp.setPullHistBinning(50)
		biasComp.setGenPdf(gen_pdf)
		biasComp.setMuVal(mu_val)
		biasComp.setEnvPdfs(cfg['env_pdfs'].split(','))
		for c in cfg['corrVals'].split(','):
			biasComp.setCorrection(c)
			biasComp.compute('mu%4.2f_gen%s_c%s'%(mu_val,gen_pdf,c))
			pullMeanGraphs[gen_pdf][c].SetPoint(p,mu_val,biasComp.getPullMean())
			pullMeanGraphs[gen_pdf][c].SetPointError(p,0.125,biasComp.getPullMeanError())
			pullWidthGraphs[gen_pdf][c].SetPoint(p,mu_val,biasComp.getPullWidth())
			pullWidthGraphs[gen_pdf][c].SetPointError(p,0.125,biasComp.getPullWidthError())
			for cov in [float(x) for x in cfg['coverageValues'].split(',')]:
				coverageGraphs[gen_pdf][c][cov].SetPoint(p,mu_val,biasComp.getCoverage(cov))
				coverageGraphs[gen_pdf][c][cov].SetPointError(p,0.125,biasComp.getCoverageError(cov))
		if 'savePVal' in cfg.keys():
			if cfg['savePVal']=='True' or cfg['savePVal']=='1':
				biasComp.setCorrection('P')
				biasComp.compute('mu%4.2f_gen%s_cP'%(mu_val,gen_pdf))
				pullMeanGraphs[gen_pdf]['P'].SetPoint(p,mu_val,biasComp.getPullMean())
				pullMeanGraphs[gen_pdf]['P'].SetPointError(p,0.,biasComp.getPullMeanError())
				pullWidthGraphs[gen_pdf]['P'].SetPoint(p,mu_val,biasComp.getPullWidth())
				pullWidthGraphs[gen_pdf]['P'].SetPointError(p,0.,biasComp.getPullWidthError())
				for cov in [float(x) for x in cfg['coverageValues'].split(',')]:
					coverageGraphs[gen_pdf]['P'][cov].SetPoint(p,mu_val,biasComp.getCoverage(cov))
					coverageGraphs[gen_pdf]['P'][cov].SetPointError(p,0.125,biasComp.getCoverageError(cov))
		del biasComp
		print 'Computation of bias and coverage done for mu=%4.2f and gen=%s'%(mu_val,gen_pdf)
		p += 1

for key, item in pullMeanGraphs.items():
	for kkey, graph in item.items():
		outfile.cd()
		graph.Write()

for key, item in pullWidthGraphs.items():
	for kkey, graph in item.items():
		outfile.cd()
		graph.Write()

for key, item in coverageGraphs.items():
	for kkey, iitem in item.items():
		for kkkey, graph in iitem.items():
			outfile.cd()
			graph.Write()

outfile.Close()

makePlots(options.filename)

"""



"""
files = [ x for x in sys.argv[1:] if '.root' in x and options.filename not in x]
if len(files)==0:
	print 'You haven\'t passed any files'
	print parser.print_help()
	sys.exit()
corrs = []
for c in options.list_of_corrs.split(','):
	if c=='P': corrs.append(c)
	else: corrs.append(float(c))
envelopeComputation(options.filename,files,corrs,options.gen_pdf_name,options.gen_inj_sig)
"""
