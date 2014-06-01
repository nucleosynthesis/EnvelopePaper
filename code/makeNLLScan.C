// Simple Macro to produce the NLL scan of the Toy Data
#include "RooAbsArg.h"
TGraph *nll2scan(double corr=0.5, RooAbsData *dat, RooAbsPdf &pdf, RooRealVar &mu){  // correction to NLL (not 2NLL)
   double xlow = -3.;
   double xhigh= 6.;
   double xstep= 0.1;

   double cfactor = corr*(pdf.getParameters(*dat)->getSize());

   TGraph *graph = new TGraph();
   RooAbsReal *nll = pdf.createNLL(*dat);
   RooMinimizer minim(*nll);
   double minnll = 10;
   double minmu=-1;
   int cpoint=0;
   RooArgSet bfparams;
   RooArgSet *nllparams = pdf.getParameters(*dat);
   nllparams->snapshot(bfparams);
   for (double x=xlow;x<=xhigh;x+=xstep) {
	mu.setVal(x);
	minim.minimize("Minuit","migrad");	
	graph.SetPoint(cpoint,x,2*(nll.getVal()+cfactor));
	if (nll.getVal() < minnll) {
		minnll=nll.getVal();
		minmu = x;
		bfparams.assignValueOnly(*nllparams);
	}
	cpoint++;
   }
   graph->SetLineWidth(2);
   graph->GetXaxis()->SetTitle("#mu");
   mu.setVal(minmu);
   // Set all parameters to best fit ones 
   nllparams->assignValueOnly(bfparams);

   return graph;   
}

void makeNLLScan(){

   TFile *fi = TFile::Open("envelopews.root");
   RooWorkspace *multipdf = fi->Get("multipdf");
   RooRealVar *x    = multipdf->var("CMS_hgg_mass");

   // The data we will use (a Toy dataset);
   RooDataHist *datatoy = multipdf->data("roohist_data_mass_cat1_toy__CMS_hgg_mass");

   // Build signal model
   RooRealVar mean("mean","mean",125); mean.setConstant();
   RooRealVar sigma("sigma","sigma",1.19); sigma.setConstant();
   RooRealVar nsignal_const("nsignal_const","nsignal_const",50.8); nsignal_const.setConstant();

   RooGaussian sig_pdf("gaus","gaus",*x,mean,sigma);
   
   // Should add this to the RooWorkspace
   RooRealVar nbkg("nbkg","nbkg",datatoy->sumEntries(),0,10E6);

   RooRealVar mu("r","r",1,-10,10);
   RooFormulaVar nsig("nsig","nsig","@0*@1",RooArgList(mu,nsignal_const));
   mu.setConstant();

   RooPlot *pl = x->frame();
   datatoy->plotOn(pl,RooFit::Binning(80));

   TCanvas *can = new TCanvas();
   // Loop through and build pdfs!, would have done this from a rooargset but roofit refuses to make one without randomly adding additional crap from the library!?!!?!?!
   int maxpdfs = 6;
   pdfit_c=0;

   int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
   int style=1;

   TList listofpdfs;
   TLegend *leg = new TLegend(0.5,0.5,0.92,0.92);
   leg->SetTextFont(42);
   leg->SetFillColor(0);

   // Bernsteins 
   //while (arg = (RooAbsArg*)pdfit->Next()) {
   for (int typeco=0;typeco<maxpdfs;typeco++){
 	// The pdf made from the new background model
        RooAbsPdf *bkg_pdf = multipdf->pdf(Form("env_pdf_1_8TeV_bern%d",typeco));
	if (!bkg_pdf) continue;
   	RooAddPdf spdf("splusb","splusb",RooArgList(sig_pdf,*bkg_pdf),RooArgList(nsig,nbkg));
	TGraph *gr = nll2scan(0.5,datatoy,spdf,mu);
	if (pdfit_c!=0 && pdfit_c%6==0) style++;
	gr->SetLineColor(color[pdfit_c%7]);
	gr->SetLineStyle(style);
	
	if (pdfit_c == 0) gr->Draw("ALP");
	else gr->Draw("LP");
	spdf.plotOn(pl,RooFit::LineColor(color[pdfit_c%7]),RooFit::LineStyle(style));
	leg->AddEntry(gr,bkg_pdf->GetName(),"L");
	listofpdfs.Add(gr);
        pdfit_c++;
	
   }
   //  
   //while (arg = (RooAbsArg*)pdfit->Next()) {
   for (int typeco=0;typeco<maxpdfs;typeco++){
 	// The pdf made from the new background model
        RooAbsPdf *bkg_pdf = multipdf->pdf(Form("env_pdf_1_8TeV_exp%d",typeco));
	if (!bkg_pdf) continue;
   	RooAddPdf spdf("splusb","splusb",RooArgList(sig_pdf,*bkg_pdf),RooArgList(nsig,nbkg));
	TGraph *gr = nll2scan(0.5,datatoy,spdf,mu);
	if (pdfit_c!=0 && pdfit_c%6==0) style++;
	gr->SetLineColor(color[pdfit_c%7]);
	gr->SetLineStyle(style);
	
	if (pdfit_c == 0) gr->Draw("LP");
	else gr->Draw("LP");
	spdf.plotOn(pl,RooFit::LineColor(color[pdfit_c%7]),RooFit::LineStyle(style));
	leg->AddEntry(gr,bkg_pdf->GetName(),"L");
	listofpdfs.Add(gr);
        pdfit_c++;
	
   }
   //while (arg = (RooAbsArg*)pdfit->Next()) {
   for (int typeco=0;typeco<maxpdfs;typeco++){
 	// The pdf made from the new background model
        RooAbsPdf *bkg_pdf = multipdf->pdf(Form("env_pdf_1_8TeV_pow%d",typeco));
	if (!bkg_pdf) continue;
   	RooAddPdf spdf("splusb","splusb",RooArgList(sig_pdf,*bkg_pdf),RooArgList(nsig,nbkg));
	TGraph *gr = nll2scan(0.5,datatoy,spdf,mu);
	if (pdfit_c!=0 && pdfit_c%6==0) style++;
	gr->SetLineColor(color[pdfit_c%7]);
	gr->SetLineStyle(style);
	
	if (pdfit_c == 0) gr->Draw("LP");
	else gr->Draw("LP");
	spdf.plotOn(pl,RooFit::LineColor(color[pdfit_c%7]),RooFit::LineStyle(style));
	leg->AddEntry(gr,bkg_pdf->GetName(),"L");
	listofpdfs.Add(gr);
        pdfit_c++;
	
   }
   for (int typeco=0;typeco<maxpdfs;typeco++){
 	// The pdf made from the new background model
        RooAbsPdf *bkg_pdf = multipdf->pdf(Form("env_pdf_1_8TeV_lau%d",typeco));
	if (!bkg_pdf) continue;
   	RooAddPdf spdf("splusb","splusb",RooArgList(sig_pdf,*bkg_pdf),RooArgList(nsig,nbkg));
	TGraph *gr = nll2scan(0.5,datatoy,spdf,mu);
	if (pdfit_c!=0 && pdfit_c%6==0) style++;
	gr->SetLineColor(color[pdfit_c%7]);
	gr->SetLineStyle(style);
	
	if (pdfit_c == 0) gr->Draw("LP");
	else gr->Draw("LP");
	spdf.plotOn(pl,RooFit::LineColor(color[pdfit_c%7]),RooFit::LineStyle(style));
	leg->AddEntry(gr,bkg_pdf->GetName(),"L");
	listofpdfs.Add(gr);
        pdfit_c++;
	
   }
   leg->Draw();
   can->Print("Profiles.png");

   TCanvas *can_fits = new TCanvas();
   pl->SetTitle("");
   pl->GetXaxis()->SetTitle("m_{#gamma#gamma}");
   pl->Draw();
   can_fits->Print("BestFits.png");
}
