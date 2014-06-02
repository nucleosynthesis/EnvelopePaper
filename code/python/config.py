import ROOT as r

class config:
	
	def __init__(self,datfilename):
		# ---- set defaults ----
		self.batchmode = True
		self.suppressRooFitOutput = False
		self.infile_name = ''
		self.ws_name = ''
		self.inv_mass_name = ''
		self.genpdf_name = ''
		self.gen_inj_sig = 0.
		self.fit_gen_to_data_first = ''
		self.envpdf_names = ['']
		self.outfile_name = 'BiasStudyProfiles.root'
		self.ntoys = 1
		self.seed = 0
		self.printScanProgress = False
		self.poi_range = [-5.,5.]
		self.poi_scan_range = [-5.,5.]
		self.points_in_scan = 10
		self.plotsForEachToy = False
		self.corrVals = [0.,1.,2.]
		self.savePVal = False
		# ---- end of setting defaults ----

		self.parseDatFile(datfilename)
		print self.batchmode
		r.gROOT.SetBatch(self.batchmode) 

	def parseDatFile(self,datfilename):
		datf = open(datfilename)
		for line in datf.readlines():
			line = line.strip('\n')
			if line==' ': continue
			if line=='\n': continue
			if not line: continue
			if line.startswith('#'): continue
			if '#' in line: line = line.split('#')[0]
			varname = line.split('=')[0]
			varval = line.split('=')[1]
			if ',' in varval: # special case for arrays
				thetype = type(getattr(self,varname)[0])
				setattr(self,varname,[])
				for v in varval.split(','):
					getattr(self,varname).append((thetype)(v))
			elif 'True' in varval:
				setattr(self,varname,True)
			elif 'False' in varval:
				setattr(self,varname,False)
			else:
				setattr(self,varname,type(getattr(self,varname))(varval))
			print varname, type(getattr(self,varname)), getattr(self,varname)

	def plotPdfs(self,mgg,dataSet,pdfs,canv):
		plot = mgg.frame()
		leg = r.TLegend(0.6,0.6,0.89,0.89)
		leg.SetFillColor(0)
		dataSet.plotOn(plot)
		for pdf in pdfs:
			pdf.plotOn(plot)
			leg.AddEntry(plot.getObject(int(plot.numItems())-1),pdf.GetName(),'L')
		plot.Draw()
		leg.Draw("same")
		canv.Modified()
		canv.Update()
