#!/usr/bin/env python
# vim: ts=2 sw=2 noexpandtab

import sys
import os
cfg = {}

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-d","--datfile",type="str")
parser.add_option("-q","--queue",type="str",default="8nh")
parser.add_option("--dryRun",default=False,action="store_true")
(options,args) = parser.parse_args()

def readConfig(datfile):

	f = open(datfile)
	for line in f.readlines():
		if line.startswith('#'): continue
		if '=' not in line: continue
		line = line.strip('\n')
		option = line.split('=')[0]
		value = line.split('=')[1]
		cfg[option] = value
	f.close()

def printConfig():
  for key, item in cfg.items():
    print key, ' -- ', item

def writeDatFile(name,gen_pdf,mu_val,outfile):
	f = open(name,'w')
	f.write('batchmode=True\n')
	f.write('suppressRooFitOutput=True\n')
	f.write('printScanProgress=False\n')
	f.write('plotsForEachToy=False\n')
	f.write('infile_name=%s\n'%cfg['infile_name'])
	f.write('ws_name=%s\n'%cfg['ws_name'])
	f.write('inv_mass_name=%s\n'%cfg['inv_mass_name'])
	f.write('genpdf_name=%s\n'%gen_pdf)
	f.write('gen_inj_sig=%4.2f\n'%mu_val)
	f.write('fit_gen_to_data_first=%s\n'%cfg['fit_gen_to_data_first'])
	f.write('envpdf_names=%s\n'%cfg['env_pdfs'])
	f.write('outfile_name=%s\n'%outfile)
	f.write('ntoys=%d\n'%int(cfg['ntoysperjob']))
	f.write('seed=%d\n'%int(cfg['seed']))
	f.write('poi_range=%s\n'%cfg['poi_range'])
	f.write('poi_scan_range=%s\n'%cfg['poi_scan_range'])
	f.write('points_in_scan=%d\n'%int(cfg['points_in_scan']))
	f.write('corrVals=%s\n'%cfg['corrVals'])
	if 'savePVal' in cfg.keys():
		if cfg['savePVal']=='True' or cfg['savePVal']=='1':
			f.write('savePVal=True\n')
		else:
			f.write('savePVal=False\n')
	f.close()

def writeSubScript(scriptname,datfilename,outfilename):
	f = open('%s/%s/scripts/%s'%(os.getcwd(),cfg['store_directory'],scriptname),'w')
	f.write('touch %s.run\n'%f.name)
	f.write('rm -f %s.fail\n'%f.name)
	f.write('rm -f %s.done\n'%f.name)
	f.write('rm -f %s.log\n'%f.name)
	f.write('#!/bin/bash\n')
	f.write('mkdir -p scratch\n')
	f.write('cd scratch\n')
	f.write('mkdir -p lib\n')
	f.write('mkdir -p python\n')
	f.write('cd %s\n'%os.getcwd())
	f.write('source setup_root_lxplus.sh\n')
	f.write('cd -\n')
	f.write('cp %s/lib/libEnvelopeCode.so lib/\n'%os.getcwd())
	f.write('cp %s/bias_study.py .\n'%os.getcwd())
	f.write('cp %s/paperStyle.C .\n'%os.getcwd())
	f.write('cp %s/python/config.py python/\n'%os.getcwd())
	f.write('cp %s/python/__init__.py python/\n'%os.getcwd())
	f.write('cp %s/%s .\n'%(os.getcwd(),cfg['infile_name']))
	f.write('cp %s/%s/dat/%s .\n'%(os.getcwd(),cfg['store_directory'],datfilename))
	f.write('if ( ./bias_study.py -i %s -t ) then\n'%datfilename)
	f.write('\ttouch %s.done\n'%f.name)
	f.write('else\n')
	f.write('\ttouch %s.fail\n'%f.name)
	f.write('fi\n')
	f.write('rm -f %s.run\n'%f.name)
	f.write('cp %s %s/%s/outfiles/\n'%(outfilename,os.getcwd(),cfg['store_directory']))
	f.close()
	os.system('chmod +x %s'%f.name)
	if not options.dryRun:
		os.system('bsub -q %s -o %s.log %s'%(options.queue,f.name,f.name))

readConfig(options.datfile)
corrValsCache = cfg['corrVals']
os.system('mkdir -p %s'%cfg['store_directory'])
os.system('mkdir -p %s/dat'%cfg['store_directory'])
os.system('mkdir -p %s/outfiles'%cfg['store_directory'])
os.system('mkdir -p %s/scripts'%cfg['store_directory'])
for val in cfg['inj_mu_vals'].split(','):
	mu_val = float(val)
	for gen_pdf in cfg['gen_pdfs'].split(','):
		for job in range(int(cfg['njobs'])):
			gen_pdf_title = gen_pdf
			if 'profiled_gen:' in gen_pdf:
				gen_pdf_title = 'profiled_gen'
			outfilename = 'BiasResults_mu%4.2f_gen%s_job%d.root'%(mu_val,gen_pdf_title,job)
			datfilename = 'cfg_mu%4.2f_gen%s_job%d.dat'%(mu_val,gen_pdf_title,job)
			scriptname = 'sub_mu%4.2f_gen%s_job%d.sh'%(mu_val,gen_pdf_title,job)
			if 'profiled_gen:' in gen_pdf: # need to split corr values for this
				for val in corrValsCache.split(','):
					cfg['corrVals'] = val
					ofname = outfilename.replace('.root','_c%s.root'%val)
					dfname = datfilename.replace('.dat','_c%s.dat'%val)
					scname = scriptname.replace('.sh','_c%s.sh'%val)
					writeDatFile(cfg['store_directory']+'/dat/'+dfname,gen_pdf,mu_val,ofname)
					writeSubScript(scname,dfname,ofname)
			else:
				writeDatFile(cfg['store_directory']+'/dat/'+datfilename,gen_pdf,mu_val,outfilename)
				writeSubScript(scriptname,datfilename,outfilename)


