#!/usr/bin/env python
# vim: ts=2 sw=2 noexpandtab

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-d","--datfile",help="Datfile used to submit jobs")
parser.add_option("-s","--splitJobs",action="store_true",default=False)
parser.add_option("-S","--isSplitJob",action="store_true",default=False)
parser.add_option("-q","--queue")
(options,args) = parser.parse_args()

import os
#os.system('bash setup_root_lxplus.sh')

import sys

import ROOT as r
r.gROOT.ProcessLine(".x paperStyle.C")

from python.biasComputer import *
r.gROOT.SetBatch()

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
	tf = r.TFile(cfg['filename'])
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

def postProcessToys(outfiledir,gen_pdf,mu_val,env_pdfs,file_names):
	# given a file should return pull histogram and which_pdf histogram and save it to a file
	# outfiledir = 'dir/envpdfs/gen_pdf/mu_val/name.root'
	sw = r.TStopwatch()
	sw.Start()
	outfname = ''
	if options.isSplitJob:
		outfname = cfg['filename']
	else:
		outfname = outfiledir+'/'+cfg['filename']
		os.system('mkdir -p %s'%outfiledir)

	outf = r.TFile(outfname,'RECREATE')
	biasTree = r.TTree('BiasTree','BiasTree')

	biasComp = biasComputer()
	biasComp.setCoverageValues([float(x) for x in cfg['coverageValues'].split(',')])
	biasComp.setOutFile(outf)
	biasComp.addFiles(file_names)
	biasComp.setPullHistBinning(50)
	biasComp.setGenPdf(gen_pdf)
	biasComp.setMuVal(mu_val)
	biasComp.setEnvPdfs(env_pdfs)
	biasComp.doPlots = False
	biasComp.setTree(biasTree)

	for c in cfg['corrVals'].split(','):
		biasComp.setCorrection(c)
		biasComp.compute('c%s'%c)
		pullHist = biasComp.getPullHist()
		covHist = biasComp.getCoverageHist()
		whichHist = biasComp.getWhichPdfHist()
		#residHist = biasComp.getResidHist()
		#errorHist = biasComp.getErrorHist()
		outf.cd()
		pullHist.Write()
		covHist.Write()
		whichHist.Write()
		pullCanv = biasComp.pullPlot('',True)
		outf.cd()
		pullCanv.SetName('pull_canv_c%s'%c)
		pullCanv.Write()
		covCanv = biasComp.coverageHistPlot('',True)
		outf.cd()
		covCanv.SetName('coverage_canv_c%s'%c)
		covCanv.Write()
		whichCanv = biasComp.whichPdfPlot('',True)
		outf.cd()
		whichCanv.SetName('pdf_choice_canv_c%s'%c)
		whichCanv.Write()

	outf.cd()
	biasTree.Write()
	outf.Close()
	del biasComp
	sw.Stop()
	print '-----------------------------------------------------'
	print 'Computation of bias and coverage done for mu=%4.2f and gen=%s'%(mu_val,gen_pdf)
	print 'Took: Real time %3.2f (secs), CPU time %3.2f (secs)'%(sw.RealTime(),sw.CpuTime())
	print '-----------------------------------------------------'

def makeNewDatFile(output_loc,pdf_set,gen_pdf,mu_val):
	f = open('%s/cfg.dat'%output_loc,'w')
	for name, value in cfg.items():
		if name=='store_directory':
			f.write('store_directory=.\n')
		elif name=='gen_pdfs':
			f.write('gen_pdfs=%s\n'%gen_pdf)
		elif name=='compute_pdf_sets':
			comp_list = '['
			for pdf in pdf_set: comp_list += pdf+','
			comp_list = comp_list[:-1]+']'
			f.write('compute_pdf_sets=%s\n'%comp_list)
		elif name=='inj_mu_vals':
			f.write('inj_mu_vals=%4.2f\n'%mu_val)
		else:
			f.write('%s=%s\n'%(name,value))
	f.close()

def writeBatchScript(output_loc,pdf_set,gen_pdf,mu_val,file_names):
	makeNewDatFile(output_loc,pdf_set,gen_pdf,mu_val)
	f = open('%s/sub.sh'%output_loc,'w')
	f.write('#!/bin/bash\n')
	f.write('mkdir -p scratch\n')
	f.write('cd scratch\n')
	f.write('cp %s/setup_root_lxplus.sh .\n'%os.getcwd())
	f.write('source setup_root_lxplus.sh\n')
	f.write('mkdir -p outfiles\n')
	for fil in file_names:
		f.write('cp %s outfiles/\n'%fil)
	f.write('cp %s/cfg.dat .\n'%output_loc)
	f.write('cp %s/scripts/extract_bias.py .\n'%os.getcwd())
	f.write('cp -r %s/python .\n'%os.getcwd())
	f.write('cp -r %s/lib .\n'%os.getcwd())
	f.write('cp %s/paperStyle.C .\n'%os.getcwd())
	f.write('if ( ./extract_bias.py -d cfg.dat -S ) then\n')
	f.write('\ttouch %s.done\n'%f.name)
	f.write('\tcp %s %s\n'%(cfg['filename'],output_loc))
	f.write('else\n')
	f.write('\ttouch %s.fail\n'%f.name)
	f.write('fi\n')
	f.close()
	os.system('chmod +x %s'%f.name)
	if options.queue:
		os.system('bsub -q %s -o %s.log %s'%(options.queue,f.name,f.name))


# MAIN BELOW:

readConfig()
outdir = cfg['store_directory']+'/'+'extracted'
if options.isSplitJob: outdir = './'

compute_pdf_sets = []
for set in cfg['compute_pdf_sets'].split(':'):
	the_set = set.strip('[').strip(']').split(',')
	compute_pdf_sets.append(the_set)
compute_pdf_names = cfg['compute_names'].split(',')

print 'Configuring batch scripts for folder %s'%(outdir)
for i, pdf_set in enumerate(compute_pdf_sets):
	compute_name = compute_pdf_names[i]
	for gen_pdf in cfg['gen_pdfs'].split(','):
		for val in cfg['inj_mu_vals'].split(','):
			mu_val = float(val)
			output_loc = '%s/%s/env_%s/genpdf_%s/mu_%4.2f'%(os.getcwd(),outdir,compute_name,gen_pdf,mu_val)
			os.system('mkdir -p %s'%output_loc)
			file_names = []
			for job in range(int(cfg['njobs'])):
				fname = '%s/%s/outfiles/BiasResults_mu%4.2f_gen%s_job%d.root'%(os.getcwd(),cfg['store_directory'],mu_val,gen_pdf,job)
				# only completed jobs
				if os.path.exists(fname):
					file_names.append(fname)

			if not options.splitJobs:
				postProcessToys(output_loc,gen_pdf,mu_val,pdf_set,file_names)
			else:
				writeBatchScript(output_loc,pdf_set,gen_pdf,mu_val,file_names)

os.system('touch %s/envelope_configurations.help'%(outdir))
for i, pdf_set in enumerate(compute_pdf_sets):
	compute_name = compute_pdf_names[i]
	os.system('touch %s/env_%s/gen_pdfs.help'%(outdir,compute_name))
	for gen_pdf in cfg['gen_pdfs'].split(','):
		os.system('touch %s/env_%s/genpdf_%s/mu_vals.help'%(outdir,compute_name,gen_pdf))

