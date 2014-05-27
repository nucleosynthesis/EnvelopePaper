#!/usr/bin/env python
# vim: ts=2 sw=2

import ROOT as r
import array
r.gROOT.ProcessLine('.L lib/libEnvelopeCode.so')
r.gROOT.ProcessLine(".x paperStyle.C")

def pullPlot(outfile,pull_hist,corrName,genpdf_name,gen_inj_sig):
	
	# make pull plot
	canv = r.TCanvas()
	pull_fit = r.TF1('pull_fit','gaus',gen_inj_sig-3,gen_inj_sig+3)
	pull_hist.Fit(pull_fit,"QN")
	pull_hist.SetLineColor(r.kBlue+2)
	pull_hist.SetLineWidth(2)
	pull_fit.SetLineColor(r.kRed)
	pull_fit.SetLineWidth(2)
	pull_hist.Draw("LEP")
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
	
	box.DrawBox(-0.14,0.,0.14,pull_hist.GetBinContent(pull_hist.FindBin(0)))
	fit_mean.DrawArrow(pull_fit.GetParameter(1),0.,pull_fit.GetParameter(1),0.25*pull_hist.GetBinContent(pull_hist.FindBin(pull_fit.GetParameter(1))),0.015,'<')
	hist_mean.DrawArrow(pull_hist.GetMean(),0.,pull_hist.GetMean(),0.25*pull_hist.GetBinContent(pull_hist.FindBin(pull_fit.GetParameter(1))),0.015,'<')
	atgen.DrawLine(gen_inj_sig,0.,gen_inj_sig,pull_hist.GetBinContent(pull_hist.FindBin(0)))
	pull_hist.Draw("LEPsame")
	pull_fit.Draw("Lsame")

	leg = r.TLegend(0.58,0.72,0.89,0.91)
	leg.SetFillColor(0)
	leg.SetLineWidth(1)
	leg.SetTextAlign(32)
	leg.AddEntry(pull_hist,"#splitline{Hist mean  = %5.3f #pm %5.3f}{Hist sigma = %5.3f #pm %5.3f}"%(pull_hist.GetMean(),pull_hist.GetMeanError(),pull_hist.GetRMS(),pull_hist.GetRMSError()),"LEP")
	leg.AddEntry(pull_fit,"#splitline{Fit mean   = %5.3f #pm %5.3f}{Fit sigma  = %5.3f #pm %5.3f}"%(pull_fit.GetParameter(1),pull_fit.GetParError(1),pull_fit.GetParameter(2),pull_fit.GetParError(2)),"L")
	leg.Draw("same")
	
	lat = r.TLatex()
	lat.SetNDC()
	gen_name = genpdf_name.split('_')[-1]
	lat.DrawLatex(0.16,0.93,"Generated with %s at #mu=%3.1f and c=%s"%(gen_name,gen_inj_sig,corrName[1:]))
	lat.DrawLatex(0.75,0.93,"%d Entries"%pull_hist.GetEntries())
	
	canv.Print("diagnostics/pull_hist_%s.pdf"%corrName)
	outfile.cd()
	pull_hist.Write()
	outfile.cd("Canvases")
	canv.SetName("pull_hist_%s"%corrName)
	canv.Write()

def whichPdfPlot(outfile,whichpdf_hist,corrName,genpdf_name,gen_inj_sig):
	
	# make fit pdf plot for number of times refitted
	canv = r.TCanvas()
	whichpdf_hist.SetFillColor(r.kBlue+1)
	whichpdf_hist.SetStats(0)
	whichpdf_hist.Draw("HIST")
	whichpdf_hist.Draw("TEXTSAME")
	lat = r.TLatex()
	lat.SetNDC()
	gen_name = genpdf_name.split('_')[-1]
	lat.DrawLatex(0.16,0.93,"Generated with %s at #mu=%3.1f and c=%s"%(gen_name,gen_inj_sig,corrName[1:]))
	lat.DrawLatex(0.75,0.93,"%d Entries"%whichpdf_hist.GetEntries())
	canv.Print("diagnostics/pdf_choice_hist_%s.pdf"%corrName)
	outfile.cd()
	whichpdf_hist.Write()
	outfile.cd("Canvases")
	canv.SetName("pdf_choice_hist_%s"%corrName)
	canv.Modified()
	canv.Update()
	canv.Write()

def envelopeComputation(outfile_name,list_of_files,list_of_corrs,genpdf_name,gen_inj_sig):

	# correction scheme setup (passed as list_of_corrs)
	# to pass exact p-value type use 'P' in list
	# otherwise should be a float with the corr val per d.o.f
	corrName = {}
	for cval in list_of_corrs:
		if cval=='P': corrName[cval] = 'cpval'
		else:
			assert(type(cval)==float)
			corrName[cval] = 'c%3.1f'%cval

	outfile = r.TFile(outfile_name,"RECREATE")
	tree = r.TTree("BiasResults","BiasResults")
	outfile.mkdir("ProfileCurves")
	outfile.mkdir("ProfileEnvelopes")
	outfile.mkdir("Canvases")
	gen_val = array.array('d', [0])
	fit_val = array.array('d', [0])
	gen_pdf = r.TString() 
	fit_pdf = r.TString()
	corr = r.TString()
	errUp_05sig = array.array('d', [0])
	errDn_05sig = array.array('d', [0])
	errUp_1sig = array.array('d', [0])
	errDn_1sig = array.array('d', [0])
	errUp_2sig = array.array('d', [0])
	errDn_2sig = array.array('d', [0])
	errUp_3sig = array.array('d', [0])
	errDn_3sig = array.array('d', [0])
	
	tree.Branch('gen_val',gen_val,'gen_val/D')
	tree.Branch('fit_val',fit_val,'fit_val/D')
	tree.Branch('gen_pdf','TString',r.AddressOf(gen_pdf))
	tree.Branch('fit_pdf','TString',r.AddressOf(fit_pdf))
	tree.Branch('corr','TString',r.AddressOf(corr))
	tree.Branch('errUp_05sig',errUp_05sig,'errUp_05sig/D')
	tree.Branch('errDn_05sig',errDn_05sig,'errDn_05sig/D')
	tree.Branch('errUp_1sig',errUp_1sig,'errUp_1sig/D')
	tree.Branch('errDn_1sig',errDn_1sig,'errDn_1sig/D')
	tree.Branch('errUp_2sig',errUp_2sig,'errUp_2sig/D')
	tree.Branch('errDn_2sig',errDn_2sig,'errDn_2sig/D')
	tree.Branch('errUp_3sig',errUp_3sig,'errUp_3sig/D')
	tree.Branch('errDn_3sig',errDn_3sig,'errDn_3sig/D')

	env_pdfs = set()
	for file in list_of_files:
		infile = r.TFile(file)
		for key in infile.GetListOfKeys():
			if 'sb_env' not in key.GetName(): continue
			name = key.GetName().split('_toy')[0]
			env_pdfs.add(name)
		infile.Close()

	pull_hist={}
	for cval in list_of_corrs:
		pull_hist[cval] = r.TH1F('pull_hist_%s'%corrName[cval],'',50,gen_inj_sig-3,gen_inj_sig+3)
		pull_hist[cval].GetXaxis().SetTitle("Pull (#mu)")
		pull_hist[cval].GetYaxis().SetTitle("Number of toy experiments")

	whichpdf_hist={}
	for cval in list_of_corrs:
		whichpdf_hist[cval] = r.TH1F('pdf_choice_hist_%s'%corrName[cval],'',len(env_pdfs),0,len(env_pdfs))
		whichpdf_hist[cval].GetYaxis().SetTitle("Number of toy experiments")
		whichpdf_hist[cval].GetXaxis().SetTitle("Best fit pdf")
	
	whichpdf_dict = {} # figure out which bin to put result in

	for i, pdf in enumerate(env_pdfs):
		whichpdf_dict[pdf] = i
		for cval in list_of_corrs:
			whichpdf_hist[cval].GetXaxis().SetBinLabel(i+1,pdf.split('_')[-1])

	toys=set()
	for file in list_of_files:
		infile = r.TFile(file)
		# within each file find the number of toys
		for key in infile.GetListOfKeys():
			if 'sb_env' in key.GetName() and 'correction' not in key.GetName():
				toys.add(int(key.GetName().split('_toy')[-1]))
		print file, len(toys)
		for toy in toys:
			# apply correction to each toy
			for c, cval in enumerate(list_of_corrs):
				profiler = r.ProfileMultiplePdfs()
				for pdf_name in env_pdfs:
					graph = infile.Get('%s_toy%d'%(pdf_name,toy))
					print graph.GetName(), ' -- ', graph.GetN()
					if c==0:
						outfile.cd('ProfileCurves')
						graph.Write()
					nparams = int(graph.GetName().split('_toy')[0][-1])
					# correction specified here
					if cval=='P':
						correction = infile.Get('%s_toy%d_correction'%(pdf_name,toy))
						profiler.addProfile(graph,correction)
					else:
						profiler.addProfile(graph,cval*nparams) 
				profiler.constructEnvelope('_toy%d_%s'%(toy,corrName[cval]))
				# if cfg.plotsForEachToy: profiler.drawEnvelope('diagnostics/envelope_toy%d.pdf'%t,True)
				envelope = profiler.getEnvelopeGraph()
				outfile.cd('ProfileEnvelopes')
				envelope.Write()

				gen_val[0] = gen_inj_sig
				fit_val[0] = profiler.getEnvelopeBestFitValue()
				gen_pdf.Resize(0)
				gen_pdf.Append(genpdf_name) 
				fit_pdf.Resize(0)
				fit_pdf.Append(profiler.getEnvelopeBestFitName())
				corr.Resize(0)
				corr.Append(corrName[cval])
				errUp_05sig[0] = profiler.getEnvelopeErrorUp(0.5)
				errDn_05sig[0] = profiler.getEnvelopeErrorDn(0.5)
				errUp_1sig[0] = profiler.getEnvelopeErrorUp(1)
				errDn_1sig[0] = profiler.getEnvelopeErrorDn(1)
				errUp_2sig[0] = profiler.getEnvelopeErrorUp(2)
				errDn_2sig[0] = profiler.getEnvelopeErrorDn(2)
				errUp_3sig[0] = profiler.getEnvelopeErrorUp(3)
				errDn_3sig[0] = profiler.getEnvelopeErrorDn(3)
				tree.Fill()
				if (fit_val[0]-gen_val[0]) >= 0.:
					pull_hist[cval].Fill((fit_val[0]-gen_val[0])/(errUp_1sig[0]-fit_val[0]))
				else:
					pull_hist[cval].Fill((fit_val[0]-gen_val[0])/(fit_val[0]-errDn_1sig[0]))

				fit_pdf_ind = whichpdf_dict[str(fit_pdf).split('_toy')[0]]
				whichpdf_hist[cval].Fill(fit_pdf_ind)
	
	outfile.cd()
	tree.Write()

	for cval, pull in pull_hist.items():
		pullPlot(outfile,pull,corrName[cval],genpdf_name,gen_inj_sig)
	for cval, whichpdf in whichpdf_hist.items():
		whichPdfPlot(outfile,whichpdf,corrName[cval],genpdf_name,gen_inj_sig)

	outfile.Close()
	"""
	# make pull plot
	canv = r.TCanvas()
	pull_fit = r.TF1('pull_fit','gaus',gen_inj_sig-3,gen_inj_sig+3)
	pull_hist.Fit(pull_fit,"QN")
	pull_hist.SetLineColor(r.kBlue+2)
	pull_hist.SetLineWidth(2)
	pull_fit.SetLineColor(r.kRed)
	pull_fit.SetLineWidth(2)
	pull_hist.Draw("LEP")
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
	
	box.DrawBox(-0.14,0.,0.14,pull_hist.GetBinContent(pull_hist.FindBin(0)))
	fit_mean.DrawArrow(pull_fit.GetParameter(1),0.,pull_fit.GetParameter(1),0.25*pull_hist.GetBinContent(pull_hist.FindBin(pull_fit.GetParameter(1))),0.015,'<')
	hist_mean.DrawArrow(pull_hist.GetMean(),0.,pull_hist.GetMean(),0.25*pull_hist.GetBinContent(pull_hist.FindBin(pull_fit.GetParameter(1))),0.015,'<')
	atgen.DrawLine(gen_inj_sig,0.,gen_inj_sig,pull_hist.GetBinContent(pull_hist.FindBin(0)))
	pull_hist.Draw("LEPsame")
	pull_fit.Draw("Lsame")

	leg = r.TLegend(0.58,0.72,0.89,0.91)
	leg.SetFillColor(0)
	leg.SetLineWidth(1)
	leg.SetTextAlign(32)
	#leg.SetHeader("%d Entries"%pull_hist.GetEntries())
	leg.AddEntry(pull_hist,"#splitline{Hist mean  = %5.3f #pm %5.3f}{Hist sigma = %5.3f #pm %5.3f}"%(pull_hist.GetMean(),pull_hist.GetMeanError(),pull_hist.GetRMS(),pull_hist.GetRMSError()),"LEP")
	leg.AddEntry(pull_fit,"#splitline{Fit mean   = %5.3f #pm %5.3f}{Fit sigma  = %5.3f #pm %5.3f}"%(pull_fit.GetParameter(1),pull_fit.GetParError(1),pull_fit.GetParameter(2),pull_fit.GetParError(2)),"L")
	leg.Draw("same")
	
	lat = r.TLatex()
	lat.SetNDC()
	gen_name = genpdf_name.split('_')[-1]
	lat.DrawLatex(0.16,0.93,"Generated with %s at #mu=%3.1f"%(gen_name,gen_inj_sig))
	lat.DrawLatex(0.75,0.93,"%d Entries"%whichpdf_hist.GetEntries())
	
	canv.Print("diagnostics/pull_hist.pdf")
	outfile.cd()
	pull_hist.Write()
	outfile.cd("Canvases")
	canv.SetName("pull_hist")
	canv.Write()

	# make fit pdf plot for number of times refitted
	whichpdf_hist.SetFillColor(r.kBlue+1)
	whichpdf_hist.SetStats(0)
	whichpdf_hist.Draw("HIST")
	whichpdf_hist.Draw("TEXTSAME")
	lat.DrawLatex(0.16,0.93,"Generated with %s at #mu=%3.1f"%(gen_name,gen_inj_sig))
	lat.DrawLatex(0.75,0.93,"%d Entries"%whichpdf_hist.GetEntries())
	canv.Print("diagnostics/pdf_choice_hist.pdf")
	outfile.cd()
	whichpdf_hist.Write()
	outfile.cd("Canvases")
	canv.SetName("pdf_choice_hist")
	canv.Modified()
	canv.Update()
	canv.Write()
	"""
