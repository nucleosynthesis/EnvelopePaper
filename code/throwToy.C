
#include "RooAbsArg.h"
TGraph *nll2scan(double corr=0.5, RooAbsData &dat, RooAbsPdf &pdf, RooRealVar &mu){  // correction to NLL (not 2NLL)
   double xlow = -1.;
   double xhigh= 2.;
   double xstep= 0.005;

   double cfactor = corr*(pdf.getParameters(dat)->getSize());
   mu.setConstant(true);
   TGraph *graph = new TGraph();
   RooAbsReal *nll = pdf.createNLL(dat);
   RooMinimizer minim(*nll);
   double minnll = 10;
   double minmu=-1;
   int cpoint=0;
   RooArgSet bfparams;
   RooArgSet *nllparams = pdf.getParameters(dat);
   nllparams->snapshot(bfparams);
   for (double x=xlow;x<=xhigh;x+=xstep) {
	mu.setVal(x);
	minim.minimize("Minuit","minimize");
	double nllvalue = nll.getVal();	
	std::cout << nllvalue << std::endl;
	graph.SetPoint(cpoint,x,2*(nllvalue+cfactor));
	if (nllvalue < minnll) {
		minnll=nllvalue;
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
   mu.setConstant(false);
   return graph;   
}

void throwToy(){
	// Generate a toy dataset from the equal parts sum of an Exponential, Laurent and power-law (these should already be fit from the real data in the workspace)
	
   gROOT->ProcessLine(".x paperStyle.C");

   TFile *fi = TFile::Open("envelopews_wsignal.root");
   RooWorkspace *multipdf = fi->Get("multipdf");
   RooRealVar *x    = multipdf->var("CMS_hgg_mass");
   RooDataHist *datatoy = multipdf->data("roohist_data_mass_cat1_toy__CMS_hgg_mass");

   double extExp =0.21;
   RooRealVar f1("f1","f1",(1./3)-0.1*extExp,0,1);  
   RooRealVar f2("f2","f2",(1./3)-0.9*extExp,0,1);  
   RooRealVar f3("f3","f3",(1./3)+extExp,0,1);  
 
   // Additional parameters 
   multipdf->var("sigma")->setConstant(true);
   multipdf->var("mean")->setConstant(true);
   multipdf->var("mean")->setVal(125);
   RooRealVar *nbkg = multipdf->var("nbkg");
   RooRealVar *mu = multipdf->var("r");
   RooRealVar *ns = multipdf->var("nsignal_const");
   //RooFormulaVar *nsig = (RooFormulaVar*)multipdf->genobj("nsig");
   RooFormulaVar *nsig = new RooFormulaVar("nsig","nsig","@0*@1",RooArgList(*mu,*ns));

   mu->setConstant(false); 
   // Signal model 
   RooAbsPdf *signal = multipdf->pdf("gaus");
   signal->getParameters(*datatoy)->Print();

   // Get the three pdfs 
   RooAbsPdf *lau =  multipdf->pdf("env_pdf_1_8TeV_lau1");
   RooAbsPdf *pow =  multipdf->pdf("env_pdf_1_8TeV_pow1");
   RooAbsPdf *exp =  multipdf->pdf("env_pdf_1_8TeV_exp1");
   RooRealVar *ha = multipdf->var("env_pdf_1_8TeV_lau1_l1");
   ha->setVal(0.5);

   mu->setVal(1.);
   RooAddPdf toymodel_pre("toygen1","toygen1",RooArgList(*lau,*pow,*exp),RooArgList(f1,f2,f3));
   RooAddPdf toymodel("toygen","toygen",RooArgList(*signal,toymodel_pre),RooArgList(*nsig,*nbkg));

   TCanvas *can = new TCanvas();
   RooPlot *fr = x->frame();
   datatoy->plotOn(fr,RooFit::Binning(80));
   toymodel.plotOn(fr);
   fr->Draw();

   TCanvas *can2 = new TCanvas();
   RooAbsData *newtoyunbinned = (RooAbsData*) toymodel.generate(*x,datatoy->sumEntries());

   x->setBins(320);
   RooDataHist newtoy("newtoy","newtoy",*x,*newtoyunbinned);

   RooAddPdf spdf_lau("splusb_lau","splusb1",RooArgList(*signal,*lau),RooArgList(*nsig,*nbkg));
   RooAddPdf spdf_pow("splusb_pow","splusb1",RooArgList(*signal,*pow),RooArgList(*nsig,*nbkg));
   RooAddPdf spdf_exp("splusb_exp","splusb1",RooArgList(*signal,*exp),RooArgList(*nsig,*nbkg));


   RooPlot *fr2 = x->frame();
   newtoy.plotOn(fr2,RooFit::Binning(80));


   TGraph *gr_lau = nll2scan(0.5,newtoy,spdf_lau,*mu);
   gr_lau->SetLineColor(1);
   spdf_lau.plotOn(fr2,RooFit::LineColor(1));
   double mLau = mu->getVal();

   TGraph *gr_exp = nll2scan(0.5,newtoy,spdf_exp,*mu);
   gr_exp->SetLineColor(2);
   spdf_exp.plotOn(fr2,RooFit::LineColor(2));
   double mExp = mu->getVal();

   TGraph *gr_pow = nll2scan(0.5,newtoy,spdf_pow,*mu);
   gr_pow->SetLineColor(4);
   spdf_pow.plotOn(fr2,RooFit::LineColor(4));
   double mPow = mu->getVal();
   fr2->Draw();
  
   TCanvas *scan = new TCanvas();
   gr_pow->Draw("AL") ;
   gr_exp->Draw("L") ;
   gr_lau->Draw("L") ;
  
   lau->getParameters(newtoy)->Print("v");
   exp->getParameters(newtoy)->Print("v");
   pow->getParameters(newtoy)->Print("v");
   std::cout << "mu (lau,exp,pow) " << mLau << " " << mExp << " " << mPow << std::endl;
  
   // Save new toy in workspace
   TFile *filenew = new TFile("envelopews_wsignal_toy1.root","RECREATE");
   newtoy.SetName("roohist_data_mass_cat1_toy1__CMS_hgg_mass");
   multipdf->import(newtoy);
   filenew->cd();
   multipdf->Write();
   filenew->Close();
   std::cout << f1.getVal() << " " << f2.getVal() << " " << f3.getVal() << std::endl;
}
