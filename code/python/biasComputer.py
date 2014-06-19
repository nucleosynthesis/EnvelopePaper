import ROOT as r

class biasComputer:

	def __init__(self):
		r.gROOT.LoadClass('ProfileMultiplePdfs','lib/libEnvelopeCode.so')
		r.gROOT.ProcessLine(".x paperStyle.C")
		self.canv = r.TCanvas()

		self.list_of_files = []
		self.gen_pdf = ''
		self.mu_val = -9999.
		self.env_pdfs = []
		self.correction = ''
		self.pull_hist_binning = 50
		self.coverageValues = [1.]
		self.outfile = 0
		self.doPlots = False

		# results
		self.ntoys = 0
		self.ntoysInCoverage = { 1. : 0}

	def printCoverageInfo(self):
		print 'nToys Tot: ', self.ntoys
		for c in self.coverageValues:
			print 'nToy cov of', c, ':', self.ntoysInCoverage[c]

	def getCoverage(self,sigma):
		return float(self.ntoysInCoverage[sigma])/float(self.ntoys)

	def getCoverageError(self,sigma):
		N = float(self.ntoys)
		k = float(self.ntoysInCoverage[sigma])
		return (1./N)*r.TMath.Sqrt(k*(1.-k/N))

	def getFracBestFitWith(self,pdf_name):
		return float(self.whichpdf_hist.GetBinContent(self.whichpdf_dict[pdf_name]))/float(self.whichpdf_hist.GetEntries())

	def __del__(self):
		del self.canv

	def setCoverageValues(self,cov_vals):
		self.coverageValues = cov_vals
		for cov in self.coverageValues:
			self.ntoysInCoverage[cov] = 0

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

	def getPullMean(self):
		return self.pull_hist.GetMean()

	def getPullMeanError(self):
		return self.pull_hist.GetMeanError()

	def getPullWidth(self):
		return self.pull_hist.GetRMS()

	def getPullWidthError(self):
		return self.pull_hist.GetRMSError()

	def getPullFitMean(self):
		return self.pull_fit.GetParameter(1)

	def getPullFitMeanError(self):
		return self.pull_fit.GetParError(1)

	def getPullFitWidth(self):
		return self.pull_fit.GetParameter(2)

	def getPullFitWidthError(self):
		return self.pull_fit.GetParError(2)

	def getPullHist(self):
		# is this memory safe?
		return self.pull_hist

	def getPullFit(self):
		# is this memory safe?
		return self.pull_fit

	def getCoverageHist(self):
		# is this memory safe?
		return self.coverage_hist

	def getWhichPdfHist(self):
		# is this memory safe?
		return self.whichpdf_hist

	def setCorrection(self,c):
		# check if ok
		if c!='P':
			try:
				tmp = float(c)
			except:
				sys.exit('%s is not a valid correction value. Must convert to float or be \'P\''%c)
		self.correction = c

	def pullPlot(self,name_ext,ret=False):

		# make self.pull plot
		self.canv.Clear()
		self.canv.SetRightMargin(0.03)
		self.canv.SetLeftMargin(0.15)
		self.canv.Modified()
		self.canv.Update()
		self.pull_hist.GetYaxis().SetTitleOffset(1.4)
		self.pull_fit = r.TF1('self.pull_fit','gaus',-3,3)
		self.pull_hist.Fit(self.pull_fit,"QN")
		self.pull_hist.SetLineColor(r.kBlue+2)
		self.pull_hist.SetLineWidth(2)
		self.pull_fit.SetLineColor(r.kRed)
		self.pull_fit.SetLineWidth(2)
		self.pull_hist.GetXaxis().SetRangeUser(self.pull_hist.GetMean()-3.,self.pull_hist.GetMean()+3.)
		self.pull_hist.Draw("LEP")
		# arrow for hist mean
		hist_mean = r.TArrow()
		hist_mean.SetLineWidth(2)
		hist_mean.SetLineColor(r.kBlue+2)
		# arrow for fit mean
		fit_mean = r.TArrow()
		fit_mean.SetLineColor(r.kRed)
		fit_mean.SetLineWidth(2)
		# line for gen val
		atgen = r.TLine()
		atgen.SetLineWidth(2)
		atgen.SetLineStyle(r.kDashed)
		atgen.SetLineColor(r.kBlack)
		# box for +/- 0.14 band
		box = r.TBox()
		boxCol = r.gROOT.GetColor(r.kGray+1)
		boxCol.SetAlpha(0.5)
		box.SetLineColor(boxCol.GetNumber())
		box.SetFillColor(boxCol.GetNumber())

		box.DrawBox(-0.14,0.,0.14,self.pull_hist.GetBinContent(self.pull_hist.FindBin(0)))
		fit_mean.DrawArrow(self.pull_fit.GetParameter(1),0.,self.pull_fit.GetParameter(1),0.25*self.pull_hist.GetBinContent(self.pull_hist.FindBin(self.pull_fit.GetParameter(1))),0.015,'<')
		hist_mean.DrawArrow(self.pull_hist.GetMean(),0.,self.pull_hist.GetMean(),0.25*self.pull_hist.GetBinContent(self.pull_hist.FindBin(self.pull_fit.GetParameter(1))),0.015,'<')
		atgen.DrawLine(0.,0.,0.,self.pull_hist.GetBinContent(self.pull_hist.FindBin(0)))
		self.pull_hist.Draw("LEPsame")
		self.pull_fit.Draw("Lsame")

		leg = r.TLegend(0.58,0.72,0.89,0.91)
		leg.SetFillColor(0)
		leg.SetLineWidth(1)
		leg.SetTextAlign(32)
		leg.AddEntry(self.pull_hist,"#splitline{Hist mean  = %5.3f #pm %5.3f}{Hist sigma = %5.3f #pm %5.3f}"%(self.pull_hist.GetMean(),self.pull_hist.GetMeanError(),self.pull_hist.GetRMS(),self.pull_hist.GetRMSError()),"LEP")
		leg.AddEntry(self.pull_fit,"#splitline{Fit mean   = %5.3f #pm %5.3f}{Fit sigma  = %5.3f #pm %5.3f}"%(self.pull_fit.GetParameter(1),self.pull_fit.GetParError(1),self.pull_fit.GetParameter(2),self.pull_fit.GetParError(2)),"L")
		leg.Draw("same")

		lat = r.TLatex()
		lat.SetNDC()
		gen_name = self.gen_pdf.split('_')[-1]
		lat.DrawLatex(0.16,0.93,"Generated with %s at #mu=%4.2f and c=%s"%(gen_name,self.mu_val,self.correction))
		lat.DrawLatex(0.75,0.93,"%d Entries"%self.pull_hist.GetEntries())

		if ret:
			return self.canv
		else:
			self.canv.Print("diagnostics/pull_hist_%s.pdf"%name_ext)
			self.outfile.cd()
			self.pull_hist.Write()
		#self.outfile.cd("Canvases")
		#canv.SetName("pull_hist_%s"%name_ext)
		#canv.Write()

	def whichPdfPlot(self,name_ext,ret=False):

		# make fit pdf plot for number of times refitted
		#canv = r.TCanvas()
		self.canv.Clear()
		self.canv.SetRightMargin(0.03)
		self.canv.SetLeftMargin(0.15)
		self.canv.Modified()
		self.canv.Update()
		self.whichpdf_hist.GetYaxis().SetTitleOffset(1.4)
		self.whichpdf_hist.SetFillColor(r.kBlue+1)
		self.whichpdf_hist.SetStats(0)
		self.whichpdf_hist.Draw("HIST")
		self.whichpdf_hist.Draw("TEXTSAME")
		lat = r.TLatex()
		lat.SetNDC()
		gen_name = self.gen_pdf.split('_')[-1]
		lat.DrawLatex(0.16,0.93,"Generated with %s at #mu=%3.1f and c=%s"%(gen_name,self.mu_val,self.correction))
		lat.DrawLatex(0.75,0.93,"%d Entries"%self.whichpdf_hist.GetEntries())

		if ret:
			return self.canv
		else:
			self.canv.Print("diagnostics/pdf_choice_hist_%s.pdf"%name_ext)
			self.outfile.cd()
			self.whichpdf_hist.Write()
		#self.outfile.cd("Canvases")
		#canv.SetName("pdf_choice_hist_%s"%name_ext)
		#canv.Modified()
		#canv.Update()
		#canv.Write()

	def coverageHistPlot(self,name_ext,ret=False):
		self.canv.Clear()
		self.canv.SetRightMargin(0.03)
		self.canv.SetLeftMargin(0.15)
		self.canv.Modified()
		self.canv.Update()
		self.coverage_hist.GetYaxis().SetTitleOffset(1.4)
		self.coverage_hist.SetFillColor(r.kGreen+1)
		self.coverage_hist.SetStats(0)
		self.coverage_hist.Draw("HIST")
		self.coverage_hist.Draw("TEXTSAME")
		lat = r.TLatex()
		lat.SetNDC()
		gen_name = self.gen_pdf.split('_')[-1]
		lat.DrawLatex(0.16,0.93,"Generated with %s at #mu=%3.1f and c=%s"%(gen_name,self.mu_val,self.correction))
		lat.DrawLatex(0.75,0.93,"%d Toys"%self.coverage_hist.GetBinContent(self.coverage_hist.GetNbinsX()))
		if ret:
			return self.canv
		else:
			self.canv.Print("diagnostics/coverage_hist_%s.pdf"%name_ext)
			self.outfile.cd()
			self.coverage_hist.Write()

	def compute(self,name_ext):

		# reset counters
		self.ntoys = 0
		for cov in self.coverageValues:
			self.ntoysInCoverage[cov] = 0
		#if not self.outfile.Get('ProfileEnvelopes'): self.outfile.mkdir('ProfileEnvelopes')
		#if not self.outfile.Get('Canvases'): self.outfile.mkdir('Canvases')

		self.pull_hist = r.TH1F('pull_hist_%s'%name_ext,'',200,-10,10)
		self.pull_hist.GetXaxis().SetTitle("Pull (#mu)")
		self.pull_hist.GetYaxis().SetTitle("Number of toy experiments")

		self.whichpdf_hist = r.TH1F('pdf_choice_hist_%s'%name_ext,'',len(self.env_pdfs),0,len(self.env_pdfs))
		self.whichpdf_hist.GetYaxis().SetTitle("Number of toy experiments")
		self.whichpdf_hist.GetXaxis().SetTitle("Best fit pdf")

		self.whichpdf_dict = {} # figure out which bin to put result in
		for i, pdf in enumerate(self.env_pdfs):
			self.whichpdf_dict['sb_%s'%pdf] = i
			self.whichpdf_hist.GetXaxis().SetBinLabel(i+1,pdf.split('_')[-1])

		self.coverage_hist = r.TH1F('coverage_hist_%s'%name_ext,'',len(self.coverageValues)+1,0,len(self.coverageValues)+1)
		self.coverage_hist.GetYaxis().SetTitle("Number of toy experiments")
		self.coverage_hist.GetXaxis().SetTitle("Coverage (#sigma)")
		self.coverage_hist_dict = {}
		for i, cov in enumerate(self.coverageValues):
			self.coverage_hist_dict[cov] = i
			self.coverage_hist.GetXaxis().SetBinLabel(i+1,'%3.1f'%cov)
		self.coverage_hist.GetXaxis().SetBinLabel(len(self.coverageValues)+1,'Total')

		for file in self.list_of_files:
			infile = r.TFile(file)
			valid=True
			toyn=0
			while valid:
				profiler = r.ProfileMultiplePdfs()
				for pdf_name in self.env_pdfs:
					graph = infile.Get('sb_%s_toy%d'%(pdf_name,toyn))
					if not graph:
						valid=False
						print 'Found %d toys in file %s'%(toyn,infile.GetName())
					if valid:
						#print graph.GetName(), ' -- ', graph.GetN()
						nparams = int(graph.GetName().split('_toy')[0][-1])
						if self.correction == 'P':
							corr = infile.Get('sb_%s_toy%d_correction'%(pdf_name,toyn))
							profiler.addProfile(graph,corr)
						else:
							profiler.addProfile(graph,float(self.correction)*nparams)
				if valid:
					profiler.constructEnvelope('_toy%d_c%s'%(toyn,self.correction))
					#profiler.drawEnvelope('diagnostics/envelope_mu%4.2f_gen%s_c%s_toy%d.pdf'%(self.mu_val,self.gen_pdf,self.correction,toyn),'#mu',True)
					envelope = profiler.getEnvelopeGraph()
					#self.outfile.cd('ProfileEnvelopes')
					#envelope.Write()

					fit_val = profiler.getEnvelopeBestFitValue()
					fit_pdf = profiler.getEnvelopeBestFitName()

					# pull
					if (fit_val-self.mu_val) >= 0:
						self.pull_hist.Fill((fit_val-self.mu_val)/(profiler.getEnvelopeErrorUp(1.)-fit_val))
					else:
						self.pull_hist.Fill((fit_val-self.mu_val)/(fit_val-profiler.getEnvelopeErrorDn(1.)))

					# coverage
					self.coverage_hist.Fill(len(self.coverageValues))
					for cov in self.coverageValues:
						errUp = profiler.getEnvelopeErrorUp(cov)
						errDn = profiler.getEnvelopeErrorDn(cov)
						cov_ind = self.coverage_hist_dict[cov]
						if (self.mu_val >= errDn and self.mu_val <= errUp):
							self.ntoysInCoverage[cov] += 1
							self.coverage_hist.Fill(cov_ind)

					# which pdf
					fit_pdf_ind = self.whichpdf_dict[str(fit_pdf).split('_toy')[0]]
					self.whichpdf_hist.Fill(fit_pdf_ind)
					toyn += 1
					self.ntoys += 1

			infile.Close()

		print '------------------------------------------------'
		print 'Finished bias computation. Found %d toys in total'%self.pull_hist.GetEntries()
		if self.doPlots:
			print 'Writing histograms to file'
			self.pullPlot(name_ext)
			self.whichPdfPlot(name_ext)

