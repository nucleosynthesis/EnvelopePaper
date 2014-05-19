#!/usr/bin/env python

# setup here
store_directory = "batch_jobs"
infile_name = 'envelopews_wsignal_toy1.root'
ws_name = 'multipdf'
inv_mass_name = 'CMS_hgg_mass'
fit_gen_to_data_first='roohist_data_mass_cat1_toy1__CMS_hgg_mass'

#gen_pdfs = ['env_pdf_1_8TeV_exp1','env_pdf_1_8TeV_pow1','env_pdf_1_8TeV_lau1','env_pdf_1_8TeV_bern1']
gen_pdfs = ['env_pdf_1_8TeV_pow1']
#inj_mu_vals = [0.,0.5,1.,1.5,2.]
inj_mu_vals = [1]
env_pdfs = ['env_pdf_1_8TeV_exp1','env_pdf_1_8TeV_exp3','env_pdf_1_8TeV_exp5',
						'env_pdf_1_8TeV_pow1','env_pdf_1_8TeV_pow3','env_pdf_1_8TeV_pow5',
						'env_pdf_1_8TeV_bern1','env_pdf_1_8TeV_bern2','env_pdf_1_8TeV_bern3',
						'env_pdf_1_8TeV_lau1','env_pdf_1_8TeV_lau2','env_pdf_1_8TeV_lau3']

njobs=10
ntoysperjob=100
points_in_scan=100
seed=0
poi_range='-5.,5.'
poi_scan_range='-3.,3.'
corrVals='0.,1.,2.'

import os

def writeDatFile(name,gen_pdf,mu_val,outfile):
	f = open(name,'w')
	f.write('batchmode=True\n')
	f.write('suppressRooFitOutput=True\n')
	f.write('printScanProgress=True\n')
	f.write('plotsForEachToy=False\n')
	f.write('infile_name=%s\n'%infile_name)
	f.write('ws_name=%s\n'%ws_name)
	f.write('inv_mass_name=%s\n'%inv_mass_name)
	f.write('genpdf_name=%s\n'%gen_pdf)
	f.write('gen_inj_sig=%3.1f\n'%mu_val)
	f.write('fit_gen_to_data_first=%s\n'%fit_gen_to_data_first)
	f.write('envpdf_names=')
	for env_pdf in env_pdfs:
		f.write('%s'%env_pdf)
		if env_pdf!=env_pdfs[-1]:
			f.write(',')
	f.write('\n')
	f.write('outfile_name=%s\n'%outfile)
	f.write('ntoys=%d\n'%ntoysperjob)
	f.write('seed=%d\n'%seed)
	f.write('poi_range=%s\n'%poi_range)
	f.write('poi_scan_range=%s\n'%poi_scan_range)
	f.write('points_in_scan=%d\n'%points_in_scan)
	f.write('corrVals=%s\n'%corrVals)
	f.close()

def writeSubScript(scriptname,datfilename,outfilename):
	f = open('%s/%s/scripts/%s'%(os.getcwd(),store_directory,scriptname),'w')
	f.write('#!/bin/bash\n')
	f.write('mkdir -p scratch\n')
	f.write('cd scratch\n')
	f.write('mkdir -p lib\n')
	f.write('mkdir -p python\n')
	f.write('cd %s\n'%os.getcwd())
	f.write('. setup_root_lxplus.sh\n')
	f.write('cd -\n')
	f.write('cp %s/lib/libEnvelopeCode.so lib/\n'%os.getcwd())
	f.write('cp %s/bias_study.py .\n'%os.getcwd())
	f.write('cp %s/paperStyle.C .\n'%os.getcwd())
	f.write('cp %s/python/config.py python/\n'%os.getcwd())
	f.write('cp %s/python/__init__.py python/\n'%os.getcwd())
	f.write('cp %s/%s .\n'%(os.getcwd(),infile_name))
	f.write('cp %s/%s/dat/%s .\n'%(os.getcwd(),store_directory,datfilename))
	f.write('if ( ./bias_study.py -i %s -t ) then\n'%datfilename)
	f.write('\ttouch %s.done\n'%f.name)
	f.write('else\n')
	f.write('\ttouch %s.fail\n'%f.name)
	f.write('fi\n')
	f.write('cp %s %s/%s/outfiles/\n'%(outfilename,os.getcwd(),store_directory))
	f.close()
	os.system('chmod +x %s'%f.name)


os.system('mkdir -p %s'%store_directory)
os.system('mkdir -p %s/dat'%store_directory)
os.system('mkdir -p %s/outfiles'%store_directory)
os.system('mkdir -p %s/scripts'%store_directory)
for mu_val in inj_mu_vals:
	for gen_pdf in gen_pdfs:
		for job in range(njobs):
			outfilename = 'BiasResults_mu%3.1f_gen%s_job%d.root'%(mu_val,gen_pdf,job)
			datfilename = 'cfg_mu%3.1f_gen%s_job%d.dat'%(mu_val,gen_pdf,job)
			scriptname = 'sub_mu%3.1f_gen%s_job%d.sh'%(mu_val,gen_pdf,job)
			writeDatFile(store_directory+'/dat/'+datfilename,gen_pdf,mu_val,outfilename)	
			writeSubScript(scriptname,datfilename,outfilename)
		
