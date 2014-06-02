class envelopeComputer:

	def __init__(self):
		import ROOT as r
		r.gROOT.ProcessLine('.L lib/libEnvelopeCode.so')
		r.gROOT.ProcessLine(".x paperStyle.C")

		self.list_of_files = []
		self.gen_pdf = ''
		self.mu_val = -9999.
		self.env_pdfs = []
		self.correction = ''
		self.pull_hist_binning = 50
	
	def setOutFile(self,outfile):
		self.outfile = outfile

	def addFiles(self,files):
		self.list_of_files = files
	
	def setPullHistBinning(self,binning):
		self.pull_hist_binning = binning

	def setGenPdf(self,name):
		self.gen_pdf = name

	def setMuVal(self,mu_val):
		self.mu_val = mu_val

	def setEnvPdfs(self,env_pdfs):
		self.env_pdfs = env_pdfs
	
	def setCorrection(self,c):
		# check if ok
		if c!='P':
			assert(float(c))
		self.correction = c

	def compute(self,name_ext):
		
		self.pull_hist = r.TH1F('pull_hist_%s'%name,'',50,-3,3)
		self.pull_hist.GetXaxis().SetTitle("Pull (#mu)")
		self.pull_hist.GetYaxis().SetTitle("Number of toy experiments")

		self.whichpdf_hist = r.TH1F('pdf_choice_hist_%s'%name,'',len(self.env_pdfs),0,len(self.env_pdfs))
		self.whichpdf_hist.GetYaxis().SetTitle("Number of toy experiments")
		self.whichpdf_hist.GetXaxis().SetTitle("Best fit pdf")

		self.whichpdf_dict = {} # figure out which bin to put result in
		for i, pdf in enumerate(self.env_pdfs):
			self.whichpdf_dict[pdf] = i
			self.whichpdf_hist.GetXaxis().SetBinLabel(i+1,pdf.split('_')[-1])

		for file in self.list_of_files:
			valid=True
			toyn=0
			while valid:
				profiler = r.ProfileMultiplePdfs()
				for pdf_name in self.env_pdfs:
					graph = infile.Get('%s_toy%d'%(pdf_name,toyn))
					if not graph: valid=False
					if valid:
						print graph.GetName(), ' -- ', graph.GetN()
						nparams = int(graph.GetName().split('_toy')[0][-1])
						if self.correction == 'P':
							corr = infile.Get(

				
			

