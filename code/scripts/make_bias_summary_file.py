#!/usr/bin/env python
# vim: ts=2 sw=2 noexpandtab

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f","--filename",default='bias_results.root')
parser.add_option("-d","--datfile")
(options,args) = parser.parse_args()

import ROOT as r
r.gROOT.ProcessLine(".x paperStyle.C")
cfg = {}
objects = []

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
		if line.startswith('mu_to_gen_pdf_map'):
			temp_dict={}
			value = value.strip('{').strip('}')
			for entry in value.split(','):
				temp_dict[entry.split(':')[0]] = entry.split(':')[1]
			cfg[option] = temp_dict
	f.close()

def getPullMeanGraph(gen_pdf,env_name,cval,verbose=False):
	# need to know gen_pdf, compute_pdf_set and cval

	objects.append(r.TGraphErrors())
	p=0
	for val in cfg['inj_mu_vals'].split(','):
		mu_val = float(val)
		gpdf = gen_pdf
		if gen_pdf=='bestfit':
			gpdf = cfg['mu_to_gen_pdf_map'][val]

		tf = r.TFile('%s/extracted/env_%s/genpdf_%s/mu_%4.2f/%s'%(cfg['store_directory'],env_name,gpdf,mu_val,cfg['filename']))
		pull_hist = tf.Get('pull_hist_c%s'%cval)
		objects[-1].SetPoint(p,mu_val,pull_hist.GetMean())
		objects[-1].SetPointError(p,0.125,pull_hist.GetMeanError())
		if verbose: print mu_val, ' : ', pull_hist.GetMean()
		tf.Close()
		p+=1

	return objects[-1]

def getPullWidthGraph(gen_pdf,env_name,cval):
	# need to know gen_pdf, compute_pdf_set and cval

	objects.append(r.TGraphErrors())
	p=0
	for val in cfg['inj_mu_vals'].split(','):
		mu_val = float(val)
		gpdf = gen_pdf
		if gen_pdf=='bestfit':
			gpdf = cfg['mu_to_gen_pdf_map'][val]

		tf = r.TFile('%s/extracted/env_%s/genpdf_%s/mu_%4.2f/%s'%(cfg['store_directory'],env_name,gpdf,mu_val,cfg['filename']))
		pull_hist = tf.Get('pull_hist_c%s'%cval)
		objects[-1].SetPoint(p,mu_val,pull_hist.GetRMS())
		objects[-1].SetPointError(p,0.125,pull_hist.GetRMSError())
		tf.Close()
		p+=1

	return objects[-1]

def getCoverageGraph(gen_pdf,env_name,cval,covval):

	objects.append(r.TGraphErrors())
	p=0
	for val in cfg['inj_mu_vals'].split(','):
		mu_val = float(val)
		gpdf = gen_pdf
		if gen_pdf=='bestfit':
			gpdf = cfg['mu_to_gen_pdf_map'][val]

		tf = r.TFile('%s/extracted/env_%s/genpdf_%s/mu_%4.2f/%s'%(cfg['store_directory'],env_name,gpdf,mu_val,cfg['filename']))
		coverage_hist = tf.Get('coverage_hist_c%s'%cval)
		bin = cfg['coverageValues'].split(',').index(covval)+1
		endbin = len(cfg['coverageValues'].split(','))+1
		N = float(coverage_hist.GetBinContent(endbin))
		k = float(coverage_hist.GetBinContent(bin))
		cov = k/N
		covErr = (1./N)*r.TMath.Sqrt(k*(1.-k/N))
		objects[-1].SetPoint(p,mu_val,1.-cov)
		objects[-1].SetPointError(p,0.125,covErr)
		tf.Close()
		p+=1
	return objects[-1]


# MAIN HERE
import os,sys
readConfig()
outf = r.TFile(options.filename,'RECREATE')

names = [x for x in cfg['compute_names'].split(',')]
pdf_sets = []
for set in cfg['compute_pdf_sets'].split(':'):
	set = set.strip('[').strip(']')
	pdf_sets.append(set.split(','))

for cval in cfg['corrVals'].split(','):
	for i, gen_pdf in enumerate(cfg['plot_gen_pdfs'].split(',')):
		for j, pdf_set in enumerate(pdf_sets):
			grMean = getPullMeanGraph(gen_pdf,names[j],cval)
			grMean.SetName('pull_mean_gen%s_fit%s_c%s'%(gen_pdf,names[j],cval))
			grWidth = getPullWidthGraph(gen_pdf,names[j],cval)
			grWidth.SetName('pull_width_gen%s_fit%s_c%s'%(gen_pdf,names[j],cval))
			outf.cd()
			grMean.Write()
			grWidth.Write()
			for covval in cfg['coverageValues'].split(','):
				gr = getCoverageGraph(gen_pdf,names[j],cval,covval)
				gr.SetName('pull_mean_gen%s_fit%s_c%s_cov%s'%(gen_pdf,names[j],cval,covval))
				outf.cd()
				gr.Write()

outf.Close()




