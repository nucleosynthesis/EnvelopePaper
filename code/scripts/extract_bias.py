#!/usr/bin/env python
# vim: ts=2 sw=2 noexpandtab

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-o","--outfilename")
parser.add_option("-d","--datfile")
parser.add_option("-b","--isBatch",action="store_true",default=False)
(options,args) = parser.parse_args()

import sys
cfg = {}

def readConfig():
	f = open(options.datfile)
  for line in f.readlines():
    if line.startswith('#'): continue
    if '=' not in line: continue
    line = line.strip('\n')
    option = line.split('=')[0]
    value = line.split('=')[1]
    cfg[option] = value
	f.close()

import ROOT as r
from python.computeEnvelope import *
if options.isBatch: r.gROOT.SetBatch()

readConfig()

# dictionary to store results
# given generated value of mu and generated pdf want to know:
# 	-- what is the bias (on pull)
# 	-- what is the coverage at 0.5, 1, 2 and 3 sigma
#		-- what function gives the best fit
# want to pass over:
#		-- generating pdf
#		-- genearting mu value
#		-- which functions to include in envelope
#		-- which correction to use

for val in cfg['inj_mu_vals'].split(','):
	mu_val = float(val)
	for gen_pdf in cfg['gen_pdfs'].split(','):
		file_names = []
		for job in range(int(cfg['njobs'])):
			file_names.append('%s/outfiles/BiasResults_mu%4.2f_gen%s_job%d.root'%(cfg['store_directory'],mu_val,gen_pdf,job))

		envelopeComputer = envelopeComputer()
		envelopeComputer.setOutFile(outfilename)
		envelopeComputer.addFiles(file_names)
		envelopeComputer.setPullHistBinning(50)
		envelopeComputer.setGenPdf(gen_pdf)
		envelopeComputer.setMuVal(mu_val)
		envelopeComputer.setEnvPdfs(cfg['env_pdfs'].split(','))
		for c in cfg['corrVals'].split(','):
			envelopeComputer.setCorrection(c)
			res = envelopeComputer.compute(pass_name_for_plots_here)
			# do something here
		if 'savePVal' in cfg.keys():
			if cfg['savePVal']=='True' or cfg['savePVal']=='1':
				envelopeComputer.setCorrection('P')
				res = envelopeComputer.compute(pass_name_for_plots_here)
				






files = [ x for x in sys.argv[1:] if '.root' in x and options.outfilename not in x]
if len(files)==0:
	print 'You haven\'t passed any files'
	print parser.print_help()
	sys.exit()
corrs = []
for c in options.list_of_corrs.split(','):
	if c=='P': corrs.append(c)
	else: corrs.append(float(c))
envelopeComputation(options.outfilename,files,corrs,options.gen_pdf_name,options.gen_inj_sig)
