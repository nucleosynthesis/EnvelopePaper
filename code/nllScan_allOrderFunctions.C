#include <algo>
#include <stdlib>
#include "RooAbsArg.h"
#include "getChisq.C"

TRandom3 *rnd = new TRandom3();


bool RUNPVALCORRECTION=true;

double MULOW = -1;
double MUHIGH = 2.5;
double MUSTEP = 0.05;
double CORRECTION = 0.5; // (to NLL (not 2NLL) per parameter)

double YMIN=206;
double YMAX=222;

double MUMIN = -1;
double MUMAX = 2.6;

bool DOEXP=1;
bool DOPOL=1;
bool DOPOW=1;
bool DOLAU=1;

bool DORANDPARS=false;
bool DOREMOVERANGE=false;
bool DOALTPARAM=false;

int BINS = 160;

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

void removeRanges(RooArgSet &pars){
	TIterator *dp = pars.createIterator();
	RooRealVar *arg;
	while (arg = (RooRealVar*)dp->Next()) {
          RooRealVar *cat = dynamic_cast<RooRealVar*>(arg);
          if (cat && !cat->isConstant()) {
		cat->removeRange();
		//cat->setVal(cat->getVal()+rnd->Gaus(0,cat->getError()));
	  }
        }
}

void randomizePars(RooArgSet &pars){

	TIterator *dp = pars.createIterator();
	RooRealVar *arg;
	while (arg = (RooRealVar*)dp->Next()) {
          RooRealVar *cat = dynamic_cast<RooRealVar*>(arg);
          if (cat && !cat->isConstant()) {
		double tval = cat->getVal();
		double dlow = tval-(cat->getMin());
		double dhig = cat->getMax() - tval;
		cat->setVal(rnd->Uniform(tval - dlow/2, tval+dhig/2));
		//cat->setVal(cat->getVal()+rnd->Gaus(0,cat->getError()));
	  }
        }

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
	if (verb) std::cout << "    mu(from-scan) = "<<bfval<<" ( @ 2xNll =  " << minimumNll << ") +" << ehig << " -" << elow << std::endl; 
}

int countNonConstants(RooArgSet *pars){

	int npar =0;
	TIterator *dp = pars->createIterator();
	RooRealVar *arg;
	while (arg = (RooRealVar*)dp->Next()) {
          RooRealVar *cat = dynamic_cast<RooRealVar*>(arg);
          if (cat && !cat->isConstant()) {
		npar++;
		//cat->setVal(cat->getVal()+rnd->Gaus(0,cat->getError()));
	  }
        }
	return npar;
}

TGraph *nll2scan(double corr=0.5, RooAbsData &dat, RooAbsPdf &pdf, RooRealVar &mu){  // correction to NLL (not 2NLL)


   double xlow = MULOW;
   double xhigh= MUHIGH;
   double xstep= MUSTEP;

   RooRealVar *var=(RooRealVar*)(pdf.getParameters((RooArgSet*)0)->find("CMS_hgg_mass"));
   mu.setConstant(true);
   TGraph *graph = new TGraph();
   RooAbsReal *nll = pdf.createNLL(dat);
   RooMinimizer minim(*nll);
   minim.setStrategy(0);

   double cfactor;
   cfactor = 2*corr*(countNonConstants(nll->getParameters(dat))); // remove mu !

//   RooFitResult *res = pdf.fitTo(dat,RooFit::Save(1));
   
   // minimize and them randomize (helps?)


   double minnll = 10000000.;
   double minmu=MULOW;
   int cpoint=0;
   RooArgSet bfparams, prefitparams;
   RooArgSet *nllparams = pdf.getParameters(dat);
   std::cout << " Pre Fit (pre random) params  .... " << std::endl;
   nllparams->Print("v");
   if (DORANDPARS) randomizePars(*nllparams);
   if (DOREMOVERANGE) removeRanges(*nllparams);
   std::cout << "Running scan, " << graph->GetName() <<", npars = " << nllparams->getSize() << std::endl;
   std::cout << " Pre Fit (post random) params  .... " << std::endl;
   nllparams->Print("v");
   nllparams->snapshot(bfparams);
   nllparams->snapshot(prefitparams);

//   nllparams->assignValueOnly((res->randomizePars()));
   for (double x=xlow;x<=xhigh;x+=xstep) {
	nllparams->assignValueOnly(prefitparams);
	mu.setVal(x);
	minim.minimize("Minuit","minimize");
	
	double nllvalue = getChisq(dat,pdf,*var,0);

	graph->SetPoint(cpoint,x,nllvalue+cfactor);

	if (nllvalue < minnll) {
		minnll=nllvalue;
		minmu = x;
		bfparams.assignValueOnly(*nllparams);
	}
	cpoint++;
   }
   graph->SetLineWidth(2);
   graph->GetXaxis()->SetTitle("#mu");
   graph->GetYaxis()->SetTitle("-2Log L + Correction");
   mu.setVal(minmu);
   // Set all parameters to best fit ones 
   nllparams->assignValueOnly(bfparams);
   mu.setConstant(false);
   // May be best to minimize overall once more to find absolute minimum
   RooMinimizer minim_float(*nll);
   minim_float.setStrategy(2);
   minim_float.minimize("Minuit","minimize");
   std::cout << " Best Fit (post-fit params  .... " << std::endl;
   nllparams->Print("v");
//   nllparams->assignValueOnly((res->randomizePars()));
   graph->SetName(pdf.GetName());
   return graph;   
}

const char * GetBkgName(const char* name){

   int npar = 0;
   std::string newname;
   TString str(name);

   if (str.Contains("pow")) {
	newname="Power Law Sum";
   }
   if (str.Contains("bern")) {
	newname="Polynomial";
   }
   if (str.Contains("lau")) {
	newname="Laurent Series";
   }
   if (str.Contains("exp")) {
	newname="Exponential Sum";
   }
   
   for (int typ=0;typ<10;typ++){
	if (str.EndsWith(Form("%d",typ)) ) { npar = typ+1; break;}
   }

//   char np[1];
//   itoa(npar,np,10); 
   newname += std::string(" (")+std::string(Form("%d",npar))+ std::string("pars)");
   return newname.c_str(); 
}

int generateEnvelope(int ng,TGraph *graphs,TGraph *gr_env,double min,double max, double step, double *minimumnll){
   int p = 0;
   double globalMin = 100000.;
   int bestpdf=0;
   for (double x=min; x<=max;x+=step){
     double minnll=1000000.;
     for (int i = 0; i < ng; i++){
	//std::cout << graphs[i]->GetName();
	double tnll = (graphs[i]).Eval(x);
	if (tnll < minnll){
		minnll = tnll;
	}
        if (tnll < globalMin){
		bestpdf = i;
	}
     }
     gr_env->SetPoint(p,x,minnll);
     p++;
     if (minnll < globalMin){
	globalMin = minnll;
     }
   }
   *minimumnll=globalMin;
   return bestpdf; 
}

void nllScan_allOrderFunctions(){
	
   RooFit::RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
   RooFit::RooMsgService::instance().setSilentMode(true);
   gROOT->ProcessLine(".x paperStyle.C");
   gROOT->ProcessLine(".L getChisq.C");

   gROOT->SetBatch(1);

   //TFile *fi = TFile::Open("envelopews_wsignal_toy1.root");
   TFile *fi = TFile::Open("envelopews_wsignal_toy1_110to150.root");
   RooWorkspace *multipdf = fi->Get("multipdf");
   RooRealVar *x    = multipdf->var("CMS_hgg_mass");
   //RooDataHist *datatoy = multipdf->data("roohist_data_mass_cat1_toy1__CMS_hgg_mass");
   RooDataHist *datatoy = multipdf->data("roohist_data_mass_cat1_toy1_cutrange__CMS_hgg_mass");
   //x->setRange("fitrnge",110,150);
   //x->setRange(110,150);
   //RooDataHist *datatoy = datatoy_in->reduce(RooFit::CutRange("fitrnge"));
   // pre-emptive fiddling 
   //multipdf->var("env_pdf_1_8TeV_pow3_f1")->setVal(0.3);
   
   // Make a scan for all of these guys   
 
   // Additional parameters 
   multipdf->var("sigma")->setConstant(true);
   multipdf->var("mean")->setConstant(true);
   RooRealVar *mu = multipdf->var("r");
   RooRealVar *ns = multipdf->var("nsignal_const");
   RooRealVar *nbkg = multipdf->var("nbkg");
   RooFormulaVar *nsig = new RooFormulaVar("nsig","nsig","@0*@1",RooArgList(*mu,*ns));

   RooAbsPdf *sig_pdf = multipdf->pdf("gaus");
   int maxpdfs = 6;
   int pdfit_c=0;

   int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
   int style=1;

   TLegend *leg = new TLegend(0.5,0.5,0.92,0.92);
   leg->SetTextFont(42);
   leg->SetFillColor(0);
   RooPlot *pl = x->frame();
   datatoy->plotOn(pl,RooFit::Binning(x->getMax()-x->getMin())); // main plot
   TCanvas *can = new TCanvas();

   double mlow,mhigh;

   TFile *outfits = new TFile("allorderfits.root","RECREATE");

   TH1F *dummyHist = new TH1F("dummy","dummy;#mu;-2Log L + correction",1,MUMIN,MUMAX);
   dummyHist->SetMinimum(YMIN);
   dummyHist->SetMaximum(YMAX);
   dummyHist->SetTitle("");
   can->cd(); dummyHist->Draw("AXIS");

   TGraph *graphs = new TGraph[50] ;
   int ngCounter=0;
   if (DOPOL){
   // Bernsteins 
   //while (arg = (RooAbsArg*)pdfit->Next()) {
   for (int typeco=0;typeco<maxpdfs;typeco++){
 	// The pdf made from the new background model
        RooAbsPdf *bkg_pdf = multipdf->pdf(Form("env_pdf_1_8TeV_bern%d",typeco));
	if (!bkg_pdf) continue;
   	RooAddPdf spdf(Form("splusb_%s",bkg_pdf->GetName()),"splusb",RooArgList(*sig_pdf,*bkg_pdf),RooArgList(*nsig,*nbkg));
	TGraph *gr = nll2scan(CORRECTION,*datatoy,spdf,*mu);
	if (pdfit_c!=0 && pdfit_c%6==0) style++;
	gr->SetLineColor(color[pdfit_c%7]);
	gr->SetLineStyle(style);
	
	//gr->Draw("AL");
	gr->Draw("L");
	spdf.plotOn(pl,RooFit::LineColor(color[pdfit_c%7]),RooFit::LineStyle(style));
	leg->AddEntry(gr,GetBkgName(bkg_pdf->GetName()),"L");
//	listofpdfs.Add(gr);
        pdfit_c++;
   	calculateMinAndErrors(gr,&mlow,&mhigh,1,1);	
	
	RooPlot *frnew = x->frame(); 
	TCanvas *cnew = new TCanvas(bkg_pdf->GetName(),"canv",800,600);
	cnew->cd();
	datatoy->plotOn(frnew,RooFit::Binning(BINS));	
	spdf.plotOn(frnew);	
	spdf.paramOn(frnew);	
	frnew->Draw();
	outfits->cd();
	cnew->Write();
	
	can->cd();
	graphs[ngCounter]=*(TGraph*)gr->Clone();
	ngCounter++;
   }
   }
   //  
   //while (arg = (RooAbsArg*)pdfit->Next()) {
   if (DOEXP){
   for (int typeco=0;typeco<maxpdfs;typeco++){
 	// The pdf made from the new background model
 	RooAbsPdf *bkg_pdf;
        if (DOALTPARAM) bkg_pdf = multipdf->pdf(Form("env_pdf_1_8TeV_expc%d",typeco));
        else bkg_pdf = multipdf->pdf(Form("env_pdf_1_8TeV_exp%d",typeco));
	if (!bkg_pdf) continue;
   	RooAddPdf spdf(Form("splusb_%s",bkg_pdf->GetName()),"splusb",RooArgList(*sig_pdf,*bkg_pdf),RooArgList(*nsig,*nbkg));
	TGraph *gr = nll2scan(CORRECTION,*datatoy,spdf,*mu);
	if (pdfit_c!=0 && pdfit_c%6==0) style++;
	gr->SetLineColor(color[pdfit_c%7]);
	gr->SetLineStyle(style);
	
	if (pdfit_c == 0) gr->Draw("L");
	else gr->Draw("L");
	spdf.plotOn(pl,RooFit::LineColor(color[pdfit_c%7]),RooFit::LineStyle(style));
	leg->AddEntry(gr,GetBkgName(bkg_pdf->GetName()),"L");
	//listofpdfs.Add(gr);
        pdfit_c++;

   	calculateMinAndErrors(gr,&mlow,&mhigh,1,1);
	RooPlot *frnew = x->frame(); 
	TCanvas *cnew = new TCanvas(bkg_pdf->GetName(),"canv",800,600);
	cnew->cd();
	datatoy->plotOn(frnew,RooFit::Binning(BINS));	
	spdf.plotOn(frnew);	
	spdf.paramOn(frnew);	
	frnew->Draw();
	outfits->cd();
	cnew->Write();
	can->cd();
	
	graphs[ngCounter]=*(TGraph*)gr->Clone();
	ngCounter++;
   }
   }

   if (DOPOW){
   //while (arg = (RooAbsArg*)pdfit->Next()) {
   for (int typeco=0;typeco<maxpdfs;typeco++){
 	// The pdf made from the new background model
 	RooAbsPdf *bkg_pdf;
        if (DOALTPARAM) bkg_pdf = multipdf->pdf(Form("env_pdf_1_8TeV_powc%d",typeco));
        else bkg_pdf = multipdf->pdf(Form("env_pdf_1_8TeV_pow%d",typeco));
	if (!bkg_pdf) continue;
   	RooAddPdf spdf(Form("splusb_%s",bkg_pdf->GetName()),"splusb",RooArgList(*sig_pdf,*bkg_pdf),RooArgList(*nsig,*nbkg));
	TGraph *gr = nll2scan(CORRECTION,*datatoy,spdf,*mu);
	if (pdfit_c!=0 && pdfit_c%6==0) style++;
	gr->SetLineColor(color[pdfit_c%7]);
	gr->SetLineStyle(style);
	
	if (pdfit_c == 0) gr->Draw("L");
	else gr->Draw("L");
	spdf.plotOn(pl,RooFit::LineColor(color[pdfit_c%7]),RooFit::LineStyle(style));
	leg->AddEntry(gr,GetBkgName(bkg_pdf->GetName()),"L");
	//listofpdfs.Add(gr);
        pdfit_c++;
	
   	calculateMinAndErrors(gr,&mlow,&mhigh,1,1) ;
	RooPlot *frnew = x->frame(); 
	TCanvas *cnew = new TCanvas(bkg_pdf->GetName(),"canv",800,600);
	cnew->cd();
	datatoy->plotOn(frnew,RooFit::Binning(BINS));	
	spdf.plotOn(frnew);	
	spdf.paramOn(frnew);	
	frnew->Draw();
	outfits->cd();
	cnew->Write();
	can->cd();

	graphs[ngCounter]=*(TGraph*)gr->Clone();
	ngCounter++;
   }
   }
   
   if (DOLAU){
   for (int typeco=0;typeco<maxpdfs;typeco++){
 	// The pdf made from the new background model
 	RooAbsPdf *bkg_pdf;
        if (DOALTPARAM) bkg_pdf = multipdf->pdf(Form("env_pdf_1_8TeV_lauc%d",typeco));
        else bkg_pdf = multipdf->pdf(Form("env_pdf_1_8TeV_lau%d",typeco));
	if (!bkg_pdf) continue;
   	RooAddPdf spdf(Form("splusb_%s",bkg_pdf->GetName()),"splusb",RooArgList(*sig_pdf,*bkg_pdf),RooArgList(*nsig,*nbkg));
	TGraph *gr = nll2scan(CORRECTION,*datatoy,spdf,*mu);
	if (pdfit_c!=0 && pdfit_c%6==0) style++;
	gr->SetLineColor(color[pdfit_c%7]);
	gr->SetLineStyle(style);
	
	if (pdfit_c == 0) gr->Draw("L");
	else gr->Draw("L");
	spdf.plotOn(pl,RooFit::LineColor(color[pdfit_c%7]),RooFit::LineStyle(style));
	leg->AddEntry(gr,GetBkgName(bkg_pdf->GetName()),"L");
	//listofpdfs.Add(gr);
        pdfit_c++;
   	calculateMinAndErrors(gr,&mlow,&mhigh,1,1);
	RooPlot *frnew = x->frame(); 
	TCanvas *cnew = new TCanvas(bkg_pdf->GetName(),"canv",800,600);
	cnew->cd();
	datatoy->plotOn(frnew,RooFit::Binning(BINS));	
	spdf.plotOn(frnew);	
	spdf.paramOn(frnew);	
	frnew->Draw();
	outfits->cd();
	cnew->Write();
	can->cd();

	graphs[ngCounter]=*(TGraph*)gr->Clone();
	ngCounter++;
   }
   }
   leg->Draw();
   can->Print("../correction/ProfilesAllOrders.pdf");

   TCanvas *can_fits = new TCanvas();
   pl->SetTitle("");
   pl->GetXaxis()->SetTitle("m_{#gamma#gamma}");
   pl->Draw();
   can_fits->Print("../correction/BestFitsAllOrders.pdf");


   // Now make the envelope graph;
   TGraph *gr_env = new TGraph();
   double globalMin;
   //std::cout << (graphs[0]).GetName() << " " << ngCounter << std::endl;;
   std::cout << " ********************************************************" << std::endl;
   int bestPdf = generateEnvelope(ngCounter,graphs,gr_env,MULOW,MUHIGH,MUSTEP,&globalMin);  // no additional correction
   std::cout << "Best pdf " << bestPdf << std::endl;
   std::cout << " ********************************************************" << std::endl;
   double mlow=0,mhigh=0; 
   double mlow2=0,mhigh2=0; 
   std::cout << calculateMinAndErrors(gr_env,&mlow,&mhigh,1,1) << std::endl;
   std::cout << " ********************************************************" << std::endl;

   std::cout << calculateMinAndErrors(gr_env,&mlow2,&mhigh2,2,0) << std::endl;
   TCanvas *envelope = new TCanvas();
   gr_env->SetLineColor(1); gr_env->SetLineWidth(2);
   double minll = globalMin ; // Pow is the best fit)
   //std::cout << minll << std::endl;
   //gr_env->GetYaxis()->SetRangeUser(minll,minll+6);
   gr_env->GetYaxis()->SetTitle("-2Log L + Correction");
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
   graphs[bestPdf].SetLineColor(2);
   graphs[bestPdf].Draw("same");
   gr_env->Draw("same"); 
   lin->Draw(""); 
   lin2->Draw("");
   lin3->Draw("");
    
   //legf->AddEntry(&graphs[bestPdf],"Single function ","L");
   legf->AddEntry(gr_env,"Minimum Envelope","L");
   legf->AddEntry(newg,"#pm 1 #sigma Interval");
   legf->AddEntry(newg2,"#pm 2 #sigma Interval");
   legf->Draw();
   envelope->SaveAs("../correction/EnvelopeAllOrders.pdf");  

   outfits->cd();
   can->Write();
   can_fits->Write();
   envelope->Write();
}
