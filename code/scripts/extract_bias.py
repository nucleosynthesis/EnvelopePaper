#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser(usage = "usage: %prog [options] filenames")
parser.add_option("-o","--outfilename")
parser.add_option("-c","--list_of_corrs",type="str",default="1,2")
parser.add_option("-g","--gen_pdf_name",type="str")
parser.add_option("-u","--gen_inj_sig",type="float")
(options,args) = parser.parse_args()

import sys
import ROOT as r
from python.computeEnvelope import *

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
