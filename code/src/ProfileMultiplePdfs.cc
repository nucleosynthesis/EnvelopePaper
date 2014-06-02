#include <map>

#include "TCanvas.h"
#include "TMath.h"
#include "TLatex.h"
#include "RooPlot.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TAxis.h"

#include "../interface/ProfileMultiplePdfs.h"

ProfileMultiplePdfs::ProfileMultiplePdfs(bool verbose):
	verbose_(verbose),
	alreadyChecked_(false),
	savePVal_(false)
{
  listOfPdfs.clear();
  //listOfPdfs = new RooArgList();
}

ProfileMultiplePdfs::~ProfileMultiplePdfs(){
  // I don't own this memory -- delete listOfPdfs;
	listOfPdfs.clear();
	// I do own this memory:
	for (map<string,TGraph*>::iterator m=profileGraphs.begin(); m!=profileGraphs.end(); m++) delete m->second;
	profileGraphs.clear();
}

void ProfileMultiplePdfs::setSavePVal(bool val){
	savePVal_ = val;
}

void ProfileMultiplePdfs::addPdf(RooAbsPdf *pdf, float penaltyTerm){
  pdfsArgSet->add(*pdf);
  listOfPdfs.insert(pair<string,pair<RooAbsPdf*,float> >(pdf->GetName(),make_pair(pdf,penaltyTerm)));
}

void ProfileMultiplePdfs::addPdfs(map<string,RooAbsPdf*> pdfs) {
	for (map<string,RooAbsPdf*>::iterator pdf=pdfs.begin(); pdf!=pdfs.end(); pdf++){
		listOfPdfs.insert(make_pair(pdf->first,make_pair(pdf->second,0.)));
	}
}

void ProfileMultiplePdfs::clearPdfs(){
  listOfPdfs.clear();
}

void ProfileMultiplePdfs::printPdfs(){
  cout << "List of Pdfs: size(" << listOfPdfs.size() << ")" << endl;
  for (map<string,pair<RooAbsPdf*,float> >::iterator m=listOfPdfs.begin(); m!=listOfPdfs.end(); m++) { 
    RooAbsPdf *pdf = m->second.first;
    cout << "\t" << pdf->GetName() << " (penalty=" << m->second.second << ")" << endl;
    cout << "\t -- "; pdf->Print();
  }
}

double ProfileMultiplePdfs::getChisq(RooAbsData &dat, RooAbsPdf &pdf, RooRealVar &var, bool prt) {

  // Find total number of events
  double nEvt;
  double nTot=0.0;

  for(int j=0;j<dat.numEntries();j++) {
    dat.get(j);
    nEvt=dat.weight();
    nTot+=nEvt;    
  }

  // Find chi-squared equivalent 2NLL
  //RooRealVar *var=(RooRealVar*)(pdf.getParameters(*dat)->find("CMS_hgg_mass"));
  double totNLL=0.0;
  double prbSum=0.0;

  for(int j=0;j<dat.numEntries();j++) {
    double m=dat.get(j)->getRealValue(var.GetName());

    // Find probability density and hence probability
    var.setVal(m);
    double prb=0.25*pdf.getVal(var);
    prbSum+=prb;

    dat.get(j);
    nEvt=dat.weight();
	  
    double mubin=nTot*prb;
    double contrib(0.);
    if (nEvt < 1) contrib = mubin;
    else contrib=mubin-nEvt+nEvt*log(nEvt/mubin);
    totNLL+=contrib;
    
    if(prt) cout << "Bin " << j << " prob = " << prb << " nEvt = " << nEvt << ", mu = " << mubin << " contribution " << contrib << endl;    
  }
  
  totNLL*=2.0;
  if(prt) cout << pdf.GetName() << " nTot = " << nTot << " 2NLL constant = " << totNLL << endl;
  
  return totNLL;
}

double ProfileMultiplePdfs::getCorrection(RooAbsData *data, RooAbsPdf *pdf, int add_pars) {
	
	// HARD-CODED!!
	RooRealVar *var = (RooRealVar*)(pdf->getParameters((RooArgSet*)0)->find("CMS_hgg_mass"));
	double minnll = 9999.;
	minnll = getChisq(*data,*pdf,*var);
	int nbins = var->getBins();
	RooArgSet *allParams = (RooArgSet*)pdf->getParameters(*data);
	RooArgSet *constParams = (RooArgSet*)allParams->selectByAttrib("Constant",kTRUE);
	int npars = allParams->getSize()-constParams->getSize()-add_pars; // HACK!
	double pvalEquiv = TMath::Prob(minnll,nbins-npars);
	// protection
	if (pvalEquiv<1.e-16) pvalEquiv=1.e-16;
	double minnllEquiv = TMath::ChisquareQuantile(1.-pvalEquiv,nbins);
	//cout << "HERE: " << nbins << " -- " << npars << " -- " << minnll << " -- " << pvalEquiv << " -- " << minnllEquiv << endl;
	return TMath::Abs(minnll - minnllEquiv);
}

// no penalty
void ProfileMultiplePdfs::makeProfiles(RooAbsData *data, RooRealVar *poi, float poiLow, float poiHigh, int npoints, string name_ext, bool printProgress) {
	
	float step_size = (poiHigh-poiLow)/float(npoints);

	//cout << "pL: " << poiLow << " pH: " << poiHigh << " np: " << npoints << " ss: " << step_size << endl;

	// loop pdfs
	for (map<string,pair<RooAbsPdf*,float> >::iterator pdfIt=listOfPdfs.begin(); pdfIt!=listOfPdfs.end(); pdfIt++){
		string pdf_name = pdfIt->first;
		RooAbsPdf *pdf = pdfIt->second.first;
		TGraph *profileCurve = new TGraph();
		profileCurve->SetName(Form("%s%s",pdf_name.c_str(),name_ext.c_str()));
		profileGraphs.insert(make_pair(pdf_name,profileCurve));
		TGraph *correctionCurve;
		if (savePVal_){
			correctionCurve = new TGraph();
			correctionCurve->SetName(Form("%s%s_correction",pdf_name.c_str(),name_ext.c_str()));
			correctionGraphs.insert(make_pair(pdf_name,correctionCurve));
		}
		// find best fit value
		RooAbsReal *nll = pdf->createNLL(*data);
		poi->setConstant(false);
		RooMinuit(*nll).migrad();
		double bestFitVal=poi->getVal();
		double bestFitNll=nll->getVal();
		double bestFitCorrection=999.;
		if (savePVal_) {
			bestFitCorrection = getCorrection(data,pdf,2); //HACK
		}
		int p=0;
		
		// if best fit val is < low scan point add to graph
		//if (bestFitVal<=poiLow) {
		//		profileGraphs[pdf_name]->SetPoint(p,bestFitVal,2.*bestFitNll);
		//		p++;
		//}
		// scan centered at best fit
		double cen  = 0.5*(poiHigh+poiLow);

		double diff = bestFitVal-cen;

		poiLow  += diff ;// bestFitVal - fabs(poiLow);
		poiHigh += diff ;//+ fabs(poiHigh);
		step_size = (poiHigh-poiLow)/float(npoints);

		cout << "Scanning.... [" << profileCurve->GetName() << " to " << data->GetName() << "]" << endl;
		for (float val=poiLow; val<poiHigh+(step_size/2.); val+=step_size) {
			if (printProgress) printProgressBar(val,poiLow,poiHigh);
			poi->setConstant(false);
			poi->setVal(val);
			poi->setConstant(true);
			RooMinuit(*nll).migrad();
			profileGraphs[pdf_name]->SetPoint(p,val,2.*nll->getVal());
			if (savePVal_) {
				double correction = getCorrection(data,pdf,1);
				correctionGraphs[pdf_name]->SetPoint(p,val,correction);
			}
			p++;
			// check where best fit value is in comparison to scan
			if (bestFitVal>val && bestFitVal<val+step_size) {
				profileGraphs[pdf_name]->SetPoint(p,bestFitVal,2.*bestFitNll);
				if (savePVal_){
					correctionGraphs[pdf_name]->SetPoint(p,bestFitVal,bestFitCorrection);
				}
				p++;
			}
			poi->setConstant(false);
		}
		if (printProgress) cout << endl;
		// if best fit val is > high scan point add to graph
		//if (bestFitVal>=poiHigh){
		//		profileGraphs[pdf_name]->SetPoint(p,bestFitVal,2.*bestFitNll);
		//		p++;
		//}
		cout << "Best fit at: x = " << bestFitVal << " 2nll = " << 2.*bestFitNll << endl;
		if (savePVal_) cout << "\t Saved correction: " << bestFitCorrection << endl;
	}
}

map<string,TGraph*> ProfileMultiplePdfs::returnProfiles(){
	return profileGraphs;
}

void ProfileMultiplePdfs::printProfiles(){
	cout << "Current stored profile TGraphs:" << endl;
	for (map<string,TGraph*>::iterator m=profileGraphs.begin(); m!=profileGraphs.end(); m++){
		cout << "\t" << m->first << " : TGraph = " << m->second->GetName() << endl;
	}
}

void ProfileMultiplePdfs::saveProfiles(TFile *outfile){
	if (verbose_) cout << "Saving stored profile TGraphs" << endl;
	for (map<string,TGraph*>::iterator m=profileGraphs.begin(); m!=profileGraphs.end(); m++){
		outfile->cd();
		m->second->Write();
	}
	if (savePVal_){
		for (map<string,TGraph*>::iterator m=correctionGraphs.begin(); m!=correctionGraphs.end(); m++){
			outfile->cd();
			m->second->Write();
		}
	}
}

void ProfileMultiplePdfs::addProfile(TGraph* graph, float corr_val){
	TGraph *correction = new TGraph();
	double x,y;
	for (int p=0; p<graph->GetN(); p++) {
		graph->GetPoint(p,x,y);
		correction->SetPoint(p,x,corr_val);
	}
	// this will need to be cleaned up later
	cleanUp_.push_back(correction);
	addProfile(graph,correction);
}

void ProfileMultiplePdfs::addProfile(TGraph* graph, TGraph* correction){
	profileGraphsWithCorrectionGraph.insert(make_pair(string(graph->GetName()),make_pair(graph,correction)));
}

/*
void ProfileMultiplePdfs::printProfilesWithNPars(){
	
	cout << "List of profiles with n free params" << endl;
	for (map<string,pair<TGraph*,int> >::iterator it=profileGraphsWithNPars.begin(); it!=profileGraphsWithNPars.end(); it++){
		cout << it->first << " : " << it->second.first->GetName() << " -- " << it->second.second << " -- " << it->second.first->GetN() << endl;
	}
}
*/

void ProfileMultiplePdfs::printProfilesWithCorrectionGraph(){
	
	cout << "List of profiles with corrections" << endl;
	for (map<string,pair<TGraph*,TGraph*> >::iterator it=profileGraphsWithCorrectionGraph.begin(); it!=profileGraphsWithCorrectionGraph.end(); it++){
		cout << it->first << " : " << it->second.first->GetName() << " -- " << it->second.second->GetName() << " -- " << it->second.first->GetN() << " -- " << it->second.second->GetN() << endl;
	}
}

void ProfileMultiplePdfs::checkPointsInGraphs() {
	
	if (alreadyChecked_) return;
	else {
		int npoints = profileGraphsWithCorrectionGraph.begin()->second.first->GetN();
		for (map<string,pair<TGraph*,TGraph*> >::iterator it=profileGraphsWithCorrectionGraph.begin(); it!=profileGraphsWithCorrectionGraph.end(); it++){
			if (it->second.first->GetN()!=npoints){
				cerr << "Unequal numbers of points between graphs to sum. Ahhh!" << endl;
				exit(1);
			}
			if (it->second.second->GetN()!=npoints){
				cerr << "Unequal number of points between graph and correction graph. Ahhh!" << endl;
				exit(1);
			}
		}
	}
}

void ProfileMultiplePdfs::constructEnvelope(string ext){

	if (profileGraphsWithCorrectionGraph.size()==0) {
		cerr << "No graphs to return envelope for! " << endl;
		exit(1);
	}
	else if (profileGraphsWithCorrectionGraph.size()==1) {
		envelope_ = profileGraphsWithCorrectionGraph.begin()->second.first;
	}
	else {
		// check graphs (and correction graphs) have same number of points
		// also find global best fit value for final graph
		envelope_ = new TGraph();
		envelope_->SetName(Form("envelope_%s",ext.c_str()));
		map<string,pair<TGraph*,TGraph*> >::iterator checker = profileGraphsWithCorrectionGraph.begin();
		int np=checker->second.first->GetN();
		envelopeMinNll_ = 9999.;
		envelopeBestFitName_ = "notfound";
		envelopeBestFitVal_ = 999.;

		for (map<string,pair<TGraph*,TGraph*> >::iterator it=profileGraphsWithCorrectionGraph.begin(); it!=profileGraphsWithCorrectionGraph.end(); it++){
			TGraph *graph = it->second.first;
			TGraph *correction = it->second.second;
			if (graph->GetN()!=np) {
				cerr << "Graphs for envelope computation have different number of points. Can't deal with this. Sorry!" << endl;
				exit(1);
			}
			if (correction->GetN()!=np){
				cerr << "Correction graph has different number of points to profile graph. Can't deal with this. Sorry!" << endl;
				exit(1);
			}
			double x,y;
			double x_corr, y_corr;
			for (int p=0; p<np; p++){
				graph->GetPoint(p,x,y);
				correction->GetPoint(p,x_corr,y_corr);
				if (TMath::Abs(x-x_corr)>1.e-3) {
					cerr << "x-values for points between profile graph and correction graph seem to far away. Can't cope. Sorry!" << endl;
					exit(1);
				}
				double corrNll = y+y_corr;
				if (corrNll < envelopeMinNll_) {
					envelopeMinNll_ = corrNll;
					envelopeBestFitName_ = it->first;
					envelopeBestFitVal_ = x;
				}
			}
		}
		// now make envelope graph
		for (int p=0; p<np; p++){
			double x,y;
			// get x value from best fit
			profileGraphsWithCorrectionGraph[envelopeBestFitName_].first->GetPoint(p,x,y);
			double minNll = 9999.;
			for (map<string,pair<TGraph*,TGraph*> >::iterator it=profileGraphsWithCorrectionGraph.begin(); it!=profileGraphsWithCorrectionGraph.end(); it++){
				TGraph *graph = it->second.first;
				TGraph *correction = it->second.second;
				double correctedNll = graph->Eval(x)+correction->Eval(x);
				if (correctedNll<minNll) minNll=correctedNll;
			}
			envelope_->SetPoint(p,x,minNll);
		}
		cout << "Constructed envelope with: " << endl;
		cout << "Best fit: " << envelopeBestFitName_ << endl;
		cout << "At nll = " << envelopeMinNll_ << endl;
		cout << "At val = " << envelopeBestFitVal_ << endl;
	}
}

TGraph *ProfileMultiplePdfs::getEnvelopeGraph(){
	return envelope_;
}

double ProfileMultiplePdfs::getEnvelopeBestFitValue(){
	return envelopeBestFitVal_;
}

string ProfileMultiplePdfs::getEnvelopeBestFitName(){
	return envelopeBestFitName_;
}

double ProfileMultiplePdfs::getEnvelopeErrorUp(double sigma){
	
	TGraph *shiftedEnv = shiftGraph(envelope_,envelopeMinNll_);
	TGraph *invertedEnv = new TGraph();
	// best fit value first then the rest
	invertedEnv->SetPoint(0,0.,envelopeBestFitVal_);
	int newP=1;
	double x,y;
	for (int oldP=0; oldP<shiftedEnv->GetN(); oldP++){
		shiftedEnv->GetPoint(oldP,x,y);
		if (x>envelopeBestFitVal_) {
			invertedEnv->SetPoint(newP,y,x);
			newP++;
		}
	}
	double crossing = sigma*sigma;
	double val = invertedEnv->Eval(crossing);
	delete invertedEnv;
	delete shiftedEnv;
	return val;
}

double ProfileMultiplePdfs::getEnvelopeErrorDn(double sigma){
	
	TGraph *shiftedEnv = shiftGraph(envelope_,envelopeMinNll_);
	TGraph *invertedEnv = new TGraph();
	int newP=0;
	double x,y;
	for (int oldP=0; oldP<shiftedEnv->GetN(); oldP++){
		shiftedEnv->GetPoint(oldP,x,y);
		if (x<envelopeBestFitVal_) {
			invertedEnv->SetPoint(newP,y,x);
			newP++;
		}
	}
	// best fit value last
	invertedEnv->SetPoint(newP,0.,envelopeBestFitVal_);
	double crossing = sigma*sigma;
	double val = invertedEnv->Eval(crossing);
	delete invertedEnv;
	delete shiftedEnv;
	return val;
}

TGraph* ProfileMultiplePdfs::shiftGraph(TGraph *graph, double shiftBy){
	
	TGraph *ret_graph = new TGraph();
	double x,y;
	for (int p=0; p<graph->GetN(); p++){
		graph->GetPoint(p,x,y);
		ret_graph->SetPoint(p,x,y-shiftBy);
	}
	return ret_graph;
}

void ProfileMultiplePdfs::drawEnvelope(string name, string xaxis, bool shiftToZero){
	TCanvas *canv = new TCanvas();
	TGraph *envelope;
	if (shiftToZero) envelope = shiftGraph(envelope_,envelopeMinNll_);
	else envelope = envelope_;
	envelope->SetLineColor(kBlue);
	envelope->SetLineWidth(3);
	envelope->SetLineStyle(kDashed);
	envelope->GetYaxis()->SetTitle("-2#DeltaLL");
	envelope->GetXaxis()->SetTitle(xaxis.c_str());
	envelope->Draw("ALP");
	for (map<string,pair<TGraph*,TGraph*> >::iterator it=profileGraphsWithCorrectionGraph.begin(); it!=profileGraphsWithCorrectionGraph.end(); it++){
		TGraph *graph = it->second.first;
		TGraph *correction = it->second.second;
		TGraph *graphWithCorrection = new TGraph();
		double x,y;
		double x_corr, y_corr;
		for (int p=0; p<graph->GetN(); p++){
			graph->GetPoint(p,x,y);
			correction->GetPoint(p,x_corr,y_corr);
			assert(TMath::Abs(x_corr-x)<1.e-3);
			graphWithCorrection->SetPoint(p,x,y+y_corr);
		}
		TGraph *graphDraw;
		if (shiftToZero) graphDraw = shiftGraph(graphWithCorrection,envelopeMinNll_);
		else graphDraw = graphWithCorrection;
		graphDraw->SetLineColor(kRed);
		graphDraw->SetLineWidth(2);
		graphDraw->Draw("LPsame");
	}
	envelope->Draw("LPsame");
	canv->Print(name.c_str());
}

void ProfileMultiplePdfs::printProgressBar(float val, float low, float high, bool asBar) {
	float percentage = 100.*(val-low)/(high-low);
	string prog = "[";
	for (int i=0; i<=100; i+=2){
		if (percentage>(float(i)-0.001)) prog += "-";
		else prog += " ";
	}
	prog += "]";
	
	if (asBar) cout << "\t" << prog << " " << Form("%3.0f%%\r",100.*(val-low)/(high-low)) << flush;
	else cout << Form("\t %3.0f%%\r",100.*(val-low)/(high-low)) << flush;
}
