
#include "RooAbsArg.h"
#include "getChisq.C"

double YMIN = 204;
double YMAX = 220;

double MULOW = -1;
double MUHIGH = 2.5;
double MUSTEP = 0.05;

void getMinPoint(TGraph *gr, double *xmin, double *ymin){

	double x,y;
	for (int n=0;n<gr->GetN();n++){
		gr->GetPoint(n,x,y);
		if (y < *ymin){
			*ymin = y ;
			*xmin = x ;
		}
	}
}
double getEnvelopeErrorUp(TGraph *envelope, double sigma){
	
	double minimumNll = 10000.;
	double bfval  	  = 0.;
	getMinPoint(envelope,&bfval,&minimumNll);

	TGraph *invertedEnv = new TGraph();
	// best fit value first then the rest
	invertedEnv->SetPoint(0,0.,bfval);
	int newP=1;
	double x,y;
	for (int oldP=0; oldP<envelope->GetN(); oldP++){
		envelope->GetPoint(oldP,x,y);
		if (x>bfval) {
			invertedEnv->SetPoint(newP,y-minimumNll,x);
			newP++;
		}
	}
	double crossing = sigma*sigma;
	double val = invertedEnv->Eval(crossing);
	return val;
}

double getEnvelopeErrorDn(TGraph *envelope, double sigma){

	double minimumNll = 10000.;
	double bfval  	  = 0.;
	getMinPoint(envelope,&bfval,&minimumNll);
	TGraph *invertedEnv = new TGraph();

	int newP=0;
	double x,y;
	for (int oldP=0; oldP<envelope->GetN(); oldP++){
		envelope->GetPoint(oldP,x,y);
		if (x<bfval) {
			invertedEnv->SetPoint(newP,y-minimumNll,x);
			newP++;
		}
	}
	// best fit value last
	invertedEnv->SetPoint(newP,0.,bfval);
	double crossing = sigma*sigma;
	double val = invertedEnv->Eval(crossing);
	return val;
}

void calculateMinAndErrors(TGraph *envelope, double *min,double *max, double sig, int verb){
	
	std::cout << envelope->GetName() << std::endl;
	double minimumNll = 10000.;
	double bfval  	  = 0.;
	getMinPoint(envelope,&bfval,&minimumNll);
        *min = getEnvelopeErrorDn(envelope,sig);
	*max = getEnvelopeErrorUp(envelope,sig);
	double elow =  bfval - *min;
	double ehig =  *max - bfval;
	if (verb) std::cout << " ( @ 2xNll =  " << minimumNll << ") +" << ehig << " -" << elow; 
}

TGraph *nll2scan(double corr=0.5, RooAbsData &dat, RooAbsPdf &pdf, RooRealVar &mu){  // correction to NLL (not 2NLL)
   double xlow = MULOW;
   double xhigh= MUHIGH;
   double xstep= MUSTEP;

   RooRealVar *var=(RooRealVar*)(pdf.getParameters((RooArgSet*)0)->find("CMS_hgg_mass"));
   double cfactor = corr*(pdf.getParameters(dat)->getSize());
   mu.setConstant(true);
   TGraph *graph = new TGraph();
   RooAbsReal *nll = pdf.createNLL(dat);
   RooMinimizer minim(*nll);
   double minnll = 10000000.;
   double minmu=MULOW;
   int cpoint=0;
   RooArgSet bfparams;
   RooArgSet *nllparams = pdf.getParameters(dat);
   nllparams->snapshot(bfparams);
   for (double x=xlow;x<=xhigh;x+=xstep) {
	mu.setVal(x);
	minim.minimize("Minuit","minimize");
	
	double nllvalue = getChisq(dat,pdf,*var,0);
	//std::cout << "nllval " << nllvalue << std::endl;
	graph->SetPoint(cpoint,x,nllvalue+2*cfactor);

	if (nllvalue < minnll) {
		minnll=nllvalue;
		minmu = x;
		bfparams.assignValueOnly(*nllparams);
	}
	cpoint++;
   }
   graph->SetLineWidth(2);
   graph->GetXaxis()->SetTitle("#mu");
   graph->GetYaxis()->SetTitle("-2Log L");
   mu.setVal(minmu);
   // Set all parameters to best fit ones 
   nllparams->assignValueOnly(bfparams);
   mu.setConstant(false);
   // May be best to minimize overall once more to find absolute minimum
   RooMinimizer minim_float(*nll);
   minim_float.minimize("Minuit","minimize");
   graph->SetName(pdf.GetName());
   return graph;   
}

void generateEnvelope(int ng,TGraph **graphs,TGraph *gr_env,double min,double max, double step){
   int p = 0;
   for (double x=min; x<=max;x+=step){
     minnll=1000000.;
     for (int i = 0; i < ng; i++){
	double tnll = graphs[i]->Eval(x);
	if (tnll < minnll){
		minnll = tnll;
	}
     }
     gr_env->SetPoint(p,x,minnll);
     p++;
   }
}

void nllScan_firstOrderFuncs(){
	
   gROOT->ProcessLine(".x paperStyle.C");
   gROOT->ProcessLine(".L getChisq.C");

   TFile *fi = TFile::Open("envelopews_wsignal_toy1_110to150.root");
   RooWorkspace *multipdf = fi->Get("multipdf");
   RooRealVar *x    = multipdf->var("CMS_hgg_mass");
   //RooDataHist *datatoy = multipdf->data("roohist_data_mass_cat1_toy1__CMS_hgg_mass");
  
   RooDataHist *datatoy = multipdf->data("roohist_data_mass_cat1_toy1_cutrange__CMS_hgg_mass");
   //x->setRange("fitrnge",110,150);
   //x->setRange(110,150);
   //RooDataHist *datatoy = datatoy_in->reduce(RooFit::CutRange("fitrnge"));
 
   // Additional parameters 
   multipdf->var("sigma")->setConstant(true);
   multipdf->var("mean")->setConstant(true);
   multipdf->var("mean")->setVal(125);
   RooRealVar *mu = multipdf->var("r");
   RooRealVar *ns = multipdf->var("nsignal_const");
   RooRealVar *nbkg = multipdf->var("nbkg");
   RooFormulaVar *nsig = new RooFormulaVar("nsig","nsig","@0*@1",RooArgList(*mu,*ns));

   RooAbsPdf *lau =  multipdf->pdf("env_pdf_1_8TeV_lau1");
   RooAbsPdf *pow =  multipdf->pdf("env_pdf_1_8TeV_pow1");
   RooAbsPdf *exp =  multipdf->pdf("env_pdf_1_8TeV_exp1");
   RooAbsPdf *pol =  multipdf->pdf("env_pdf_1_8TeV_bern1");
   RooAbsPdf *signal = multipdf->pdf("gaus");
   // This was the function used to throw the toy
   /*
   double extExp =0.21;
   RooRealVar f1("f1","f1",(1./3)-0.1*extExp,0,1);  
   RooRealVar f2("f2","f2",(1./3)-0.9*extExp,0,1);  
   RooRealVar f3("f3","f3",(1./3)+extExp,0,1);  

   //mu->setConstant(false); 
   // Signal model 
   RooAbsPdf *signal = multipdf->pdf("gaus");
   signal->getParameters(*datatoy)->Print();

   // Get the three pdfs 
   RooAbsPdf *lau =  multipdf->pdf("env_pdf_1_8TeV_lau1");
   RooAbsPdf *pow =  multipdf->pdf("env_pdf_1_8TeV_pow1");
   RooAbsPdf *exp =  multipdf->pdf("env_pdf_1_8TeV_exp1");
//   RooRealVar *ha = multipdf->var("env_pdf_1_8TeV_lau1_l1");
//   ha->setVal(0.5);

   //mu->setVal(1.);
   //RooAddPdf toymodel_pre("toygen1","toygen1",RooArgList(*lau,*pow,*exp),RooArgList(f1,f2,f3));
   //RooAddPdf toymodel("toygen","toygen",RooArgList(*signal,toymodel_pre),RooArgList(*nsig,*nbkg));
   */
   /*
   TCanvas *can = new TCanvas();
   RooPlot *fr = x->frame();
   datatoy->plotOn(fr,RooFit::Binning(80));
   toymodel.plotOn(fr);
   fr->Draw();

   TCanvas *can2 = new TCanvas();
   RooAbsData *newtoyunbinned = (RooAbsData*) toymodel.generate(*x,datatoy->sumEntries());

   x->setBins(320);
   RooDataHist newtoy("newtoy","newtoy",*x,*newtoyunbinned);
   */

   RooAddPdf spdf_lau("splusb_lau","splusb1",RooArgList(*signal,*lau),RooArgList(*nsig,*nbkg));
   RooAddPdf spdf_pow("splusb_pow","splusb1",RooArgList(*signal,*pow),RooArgList(*nsig,*nbkg));
   RooAddPdf spdf_exp("splusb_exp","splusb1",RooArgList(*signal,*exp),RooArgList(*nsig,*nbkg));
   RooAddPdf spdf_pol("splusb_pol","splusb1",RooArgList(*signal,*pol),RooArgList(*nsig,*nbkg));
 

   RooPlot *fr2 = x->frame();
   datatoy->plotOn(fr2,RooFit::Binning(x->getMax()-x->getMin()));

   TLegend *leg = new TLegend(0.42,0.62,0.71,0.89);
   leg->SetFillColor(0);
   

   TGraph *gr_lau = nll2scan(0.,*datatoy,spdf_lau,*mu);
   gr_lau->SetLineColor(kGreen+2);
   leg->AddEntry(gr_lau,"Laurent","L");
   spdf_lau.plotOn(fr2,RooFit::LineColor(kGreen+2));
   double mLau = mu->getVal();

   TGraph *gr_exp = nll2scan(0.,*datatoy,spdf_exp,*mu);
   gr_exp->SetLineColor(2);
   leg->AddEntry(gr_exp,"Exponential","L");
   spdf_exp.plotOn(fr2,RooFit::LineColor(2));
   double mExp = mu->getVal();

   TGraph *gr_pow = nll2scan(0.,*datatoy,spdf_pow,*mu);
   gr_pow->SetLineColor(4);
   leg->AddEntry(gr_pow,"Power Law","L");
   spdf_pow.plotOn(fr2,RooFit::LineColor(4));
   double mPow = mu->getVal();

   TLegend *leg2 = (TLegend*)leg->Clone();
   
   TGraph *gr_pol = nll2scan(0.,*datatoy,spdf_pol,*mu);
   gr_pol->SetLineColor(kViolet);
   double mPol = mu->getVal();
   leg->AddEntry(gr_pol,"Polynomial","L");
   spdf_pol.plotOn(fr2,RooFit::LineColor(kViolet));
   leg2->AddEntry(gr_pol,"Polynomial","L");

   TCanvas *fits = new TCanvas();
   fr2->SetTitle("");
   fr2->GetXaxis()->SetTitle("m_{#gamma#gamma}");
   fr2->GetYaxis()->SetTitle("Events / GeV");
   fr2->Draw(); leg2->Draw(); 
   fits->SaveAs("../functions/BestFits.pdf");

   gr_pow->GetXaxis()->SetRangeUser(MULOW,MUHIGH);
   gr_pow->GetYaxis()->SetRangeUser(YMIN,YMAX);
  
   TCanvas *scan = new TCanvas();
   gr_pow->Draw("AL") ; 
   gr_exp->Draw("L") ;
   gr_lau->Draw("L") ;
   gr_pol->Draw("L");
   leg->Draw(); scan->SaveAs("../functions/Profiles.pdf");

   // Construct the envelope ?
   TGraph *gr_env = new TGraph();
   gr_env->SetName("Envelope");
   TGraph *graphs[4] = {gr_pol,gr_exp,gr_lau,gr_pow}; 
   generateEnvelope(4,graphs,gr_env,MULOW,MUHIGH,MUSTEP);  // no correction
   double mlow=0,mhigh=0; 
   double mlow2=0,mhigh2=0; 
   //std::cout << "mu (lau,exp,pow,pol) " << mLau << " " << mExp << " " << mPow << " " << mPol<< std::endl;
   std::cout << " ************* Best Fits and uncertainties **************" << std::endl;
   std::cout << "mu = "<< mPol <<  ", " << calculateMinAndErrors(gr_pol,&mlow,&mhigh,1,1) << std::endl;
   std::cout << "mu = "<< mPow <<  ", " << calculateMinAndErrors(gr_pow,&mlow,&mhigh,1,1) << std::endl;
   std::cout << "mu = "<< mExp <<  ", " << calculateMinAndErrors(gr_exp,&mlow,&mhigh,1,1) << std::endl;
   std::cout << "mu = "<< mLau <<  ", " << calculateMinAndErrors(gr_lau,&mlow,&mhigh,1,1) << std::endl;
   std::cout << calculateMinAndErrors(gr_env,&mlow,&mhigh,1,1) << std::endl;
   std::cout << " ********************************************************" << std::endl;

   std::cout << calculateMinAndErrors(gr_env,&mlow2,&mhigh2,2,0) << std::endl;
   TCanvas *envelope = new TCanvas();
   gr_env->SetLineColor(1); gr_env->SetLineWidth(2);
   double minll = gr_env->Eval( mPow ) ; // Pow is the best fit)
   std::cout << minll << std::endl;
   //gr_env->GetYaxis()->SetRangeUser(minll,minll+6);
   gr_env->GetYaxis()->SetTitle("-2Log L");
   gr_env->GetXaxis()->SetTitle("#mu");
   gr_env->GetXaxis()->SetRangeUser(MULOW,MUHIGH);
   gr_env->GetYaxis()->SetRangeUser(YMIN,YMAX);
   gr_env->Draw("AL");
   TLine *lin = new TLine(MULOW,minll,MUHIGH,minll);
   TLine *lin2 = new TLine(MULOW,minll+1,MUHIGH,minll+1);
   TLine *lin3 = new TLine(MULOW,minll+4,MUHIGH,minll+4);
   lin->SetLineColor(1); 
   lin2->SetLineColor(1); 
   lin2->SetLineColor(1);
   lin->SetLineWidth(2); 
   lin3->SetLineWidth(2); 
   lin2->SetLineWidth(2);
   lin->SetLineStyle(2); 
   lin2->SetLineStyle(2);
   lin3->SetLineStyle(2);
   TGraphAsymmErrors *newg = new TGraphAsymmErrors() ;
   TGraphAsymmErrors *newg2 = new TGraphAsymmErrors() ;
   int np = 0;
    
   for (double xcheck = mlow;xcheck<=mhigh;xcheck+=MUSTEP){
	newg->SetPoint(np,xcheck,minll);
	newg->SetPointError(np,0.01,0,0,gr_env->Eval(xcheck)-minll);
	np++;
   }

   np = 0;
   for (double xcheck = mlow2;xcheck<=mhigh2;xcheck+=MUSTEP*0.1){
	newg2->SetPoint(np,xcheck,minll);
	newg2->SetPointError(np,0.01,0,0,gr_env->Eval(xcheck)-minll);
	np++;
   }

   TLegend *legf = new TLegend(0.30,0.62,0.79,0.89);
   newg->SetFillColor(kGreen+3); newg->SetLineColor(kGreen+3);newg->SetMarkerColor(kGreen+3);
   newg2->SetFillColor(kYellow); newg2->SetLineColor(kYellow);newg2->SetMarkerColor(kYellow);
   newg2->Draw("E3same");
   newg->Draw("E3same");
   gr_env->Draw("same"); 
   lin->Draw(""); 
   lin2->Draw("");
   lin3->Draw("");
    
   legf->AddEntry(gr_env,"Minimum Envelope","L");
   legf->AddEntry(newg,"#pm 1 #sigma Interval");
   legf->AddEntry(newg2,"#pm 2 #sigma Interval");
   legf->Draw();
   envelope->SaveAs("../functions/Envelope.pdf");  
   /* 
   // Save new toy in workspace
   TFile *filenew = new TFile("envelopews_wsignal_toy1.root","RECREATE");
   newtoy.SetName("roohist_data_mass_cat1_toy1__CMS_hgg_mass");
   multipdf->import(newtoy);
   filenew->cd();
   multipdf->Write();
   filenew->Close();
   std::cout << f1.getVal() << " " << f2.getVal() << " " << f3.getVal() << std::endl;
   */
}
