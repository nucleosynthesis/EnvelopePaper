#!/usr/bin/env python
# vim: ts=2 sw=2 noexpandtab

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f","--filename",default='bias_results.root')
parser.add_option("-d","--datfile")
parser.add_option("-D","--outdir",default='./', help='This is where the output from extract_bias went')
parser.add_option("-i","--interactive",default=False,action="store_true")
(options,args) = parser.parse_args()

import ROOT as r
r.gROOT.ProcessLine(".x paperStyle.C")
cfg = {}
objects = []
if not options.interactive: r.gROOT.SetBatch()

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

def getPullMeanGraph(gen_pdf,env_set,cval,verbose=False):
	# need to know gen_pdf, compute_pdf_set and cval
	env_dir_name = ''
	for pdf in env_set: env_dir_name += pdf

	objects.append(r.TGraphErrors())
	p=0
	for val in cfg['inj_mu_vals'].split(','):
		mu_val = float(val)
		gpdf = gen_pdf
		if gen_pdf=='bestfit':
			gpdf = cfg['mu_to_gen_pdf_map'][val]

		tf = r.TFile('%s/env_%s/genpdf_%s/mu_%4.2f/%s'%(options.outdir,env_dir_name,gpdf,mu_val,cfg['filename']))
		pull_hist = tf.Get('pull_hist_c%s'%cval)
		objects[-1].SetPoint(p,mu_val,pull_hist.GetMean())
		objects[-1].SetPointError(p,0.125,pull_hist.GetMeanError())
		if verbose: print mu_val, ' : ', pull_hist.GetMean()
		tf.Close()
		p+=1

	return objects[-1]

def getPullWidthGraph(gen_pdf,env_set,cval):
	# need to know gen_pdf, compute_pdf_set and cval
	env_dir_name = ''
	for pdf in env_set: env_dir_name += pdf

	objects.append(r.TGraphErrors())
	p=0
	for val in cfg['inj_mu_vals'].split(','):
		mu_val = float(val)
		gpdf = gen_pdf
		if gen_pdf=='bestfit':
			gpdf = cfg['mu_to_gen_pdf_map'][val]

		tf = r.TFile('%s/env_%s/genpdf_%s/mu_%4.2f/%s'%(options.outdir,env_dir_name,gpdf,mu_val,cfg['filename']))
		pull_hist = tf.Get('pull_hist_c%s'%cval)
		objects[-1].SetPoint(p,mu_val,pull_hist.GetRMS())
		objects[-1].SetPointError(p,0.125,pull_hist.GetRMSError())
		tf.Close()
		p+=1

	return objects[-1]

def getCoverageGraph(gen_pdf,env_set,cval,covval):
	env_dir_name = ''
	for pdf in env_set: env_dir_name += pdf

	objects.append(r.TGraphErrors())
	p=0
	for val in cfg['inj_mu_vals'].split(','):
		mu_val = float(val)
		gpdf = gen_pdf
		if gen_pdf=='bestfit':
			gpdf = cfg['mu_to_gen_pdf_map'][val]

		tf = r.TFile('%s/env_%s/genpdf_%s/mu_%4.2f/%s'%(options.outdir,env_dir_name,gpdf,mu_val,cfg['filename']))
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
os.system('mkdir -p %s/plots'%options.outdir)
outf = r.TFile(options.filename,'RECREATE')

min = float(cfg['inj_mu_vals'].split(',')[0])-0.25
max = float(cfg['inj_mu_vals'].split(',')[-1])+0.25
dummyHist = r.TH1F('dummy','',1,min,max)
dummyHist.SetStats(0)
dummyHist.SetLineColor(0)
dummyHist.SetFillColor(0)
colors = [int(x) for x in cfg['compute_colors'].split(',')]
names = [x for x in cfg['compute_names'].split(',')]
pdf_sets = []
for set in cfg['compute_pdf_sets'].split(':'):
	set = set.strip('[').strip(']')
	pdf_sets.append(set.split(','))

line = r.TLine()
line.SetLineColor(r.kGray)
line.SetLineWidth(2)
line.SetLineStyle(r.kDashed)

ngen = len(cfg['gen_pdfs'].split(','))

# bias plots
for cval in cfg['corrVals'].split(','):
	#canvs = []
	canv = r.TCanvas('bias_mean_c%s'%cval,'',400,200*ngen)
	#canv.DivideSquare(len(cfg['gen_pdfs'].split(',')))
	for i, gen_pdf in enumerate(cfg['gen_pdfs'].split(',')):
		pad = r.TPad('bias_mean_gen%s_c%s'%(gen_pdf.split('_')[-1],cval),'',0.,float(i)/float(ngen),1.,float(i)/float(ngen)+1./float(ngen))
		#canv.cd(i+1)
		#canv = r.TCanvas('bias_mean_gen%s_c%s'%(gen_pdf.split('_')[-1],cval))
		pad.cd()
		dummyHist.SetTitle('Generated from %s with c=%s'%(gen_pdf.split('_')[-1],cval))
		dummyHist.GetXaxis().SetTitle('Injected signal (#mu)')
		dummyHist.GetYaxis().SetTitle('Pull (#mu)')
		dummyHist.GetYaxis().SetRangeUser(-1.5,1.5)
		dummyHist.Draw()
		line.DrawLine(min,0.,max,0.)
		leg = r.TLegend(0.6,0.6,0.89,0.89)
		leg.SetFillColor(0)
		for j, pdf_set in enumerate(pdf_sets):
			gr = getPullMeanGraph(gen_pdf,pdf_set,cval)
			gr.SetName('pull_mean_gen%s_fit%s_c%s'%(gen_pdf,names[j],cval))
			gr.SetLineWidth(2)
			gr.SetLineColor(colors[j])
			gr.Draw("EPsame")
			leg.AddEntry(gr,names[j],'L')
			outf.cd()
			gr.Write()
		leg.Draw("same")
		canv.cd()
		pad.Draw()
		#canv.Update()
		#canv.Modified()
		#canvs.append(canv.Clone())
		#canv.Print("%s/plots/bias_mean_gen%s_c%s.pdf"%(options.outdir,gen_pdf.split('_')[-1],cval))
		outf.cd()
		pad.Write()
	canv.Update()
	canv.Modified()
	canv.Print("%s/plots/bias_mean_c%s.pdf"%(options.outdir,cval))

	"""
	gCanv = r.TCanvas('bias_mean_c%s'%cval)
	gCanv.DivideSquare(len(cfg['corrVals'].split(',')))
	for i, canv in enumerate(canvs):
		gCanv.cd(i+1)
		canv.Draw()
	gCanv.Print('%s/plots/bias_mean_c%s.pdf'%(options.outdir,cval))
	"""


# coverage plots

for cval in cfg['corrVals'].split(','):
	canvs = []
	#canv = r.TCanvas('coverage_c%s'%cval)
	#canv.DivideSquare(len(cfg['gen_pdfs'].split(',')))
	for i, gen_pdf in enumerate(cfg['gen_pdfs'].split(',')):
		#canv.cd(i)
		canv = r.TCanvas('bias_mean_gen%s_c%s'%(gen_pdf.split('_')[-1],cval))
		dummyHist.SetTitle('Generated from %s with c=%s'%(gen_pdf.split('_')[-1],cval))
		dummyHist.GetYaxis().SetRangeUser(1.e-3,1.)
		dummyHist.GetXaxis().SetTitle('Injected signal (#mu)')
		dummyHist.GetYaxis().SetTitle('1-CL_{s}')
		dummyHist.Draw()
		for j, pdf_set in enumerate(pdf_sets):
			for covval in cfg['coverageValues'].split(','):
				line.DrawLine(min,r.TMath.Prob(float(covval)**2,1),max,r.TMath.Prob(float(covval)**2,1))
				gr = getCoverageGraph(gen_pdf,pdf_set,cval,covval)
				gr.SetName('pull_mean_gen%s_fit%s_c%s_cov%s'%(gen_pdf,names[j],cval,covval))
				gr.SetLineWidth(2)
				gr.SetLineColor(colors[j])
				gr.Draw("EPsame")
				outf.cd()
				gr.Write()
		canv.Update()
		canv.SetLogy()
		leg.Draw("same")
		canv.Print("%s/plots/coverage_gen%s_c%s.pdf"%(options.outdir,gen_pdf,cval))
		outf.cd()
		canv.Write()

outf.Close()




