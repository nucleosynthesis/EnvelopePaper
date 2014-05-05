#include <map>

#include "TCanvas.h"
#include "TMath.h"
#include "TLatex.h"
#include "RooPlot.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "TF1.h"
#include "TMatrixD.h"

#include "../interface/ProfileMultiplePdfs.h"

ProfileMultiplePdfs::ProfileMultiplePdfs(bool verbose):
	verbose_(verbose),
	corr_pdof_(0.)
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
		// find best fit value
		RooAbsReal *nll = pdf->createNLL(*data);
		poi->setConstant(false);
		RooMinuit(*nll).migrad();
		double bestFitVal=poi->getVal();
		double bestFitNll=nll->getVal();
		int p=0;
		
		// if best fit val is < low scan point add to graph
		if (bestFitVal<=poiLow) {
				profileGraphs[pdf_name]->SetPoint(p,bestFitVal,2.*bestFitNll);
				p++;
		}
		// scan
		cout << "Scanning.... [" << profileCurve->GetName() << " to " << data->GetName() << "]" << endl;
		for (float val=poiLow; val<poiHigh+(step_size/2.); val+=step_size) {
			if (printProgress) printProgressBar(val,poiLow,poiHigh);
			poi->setConstant(false);
			poi->setVal(val);
			poi->setConstant(true);
			RooMinuit(*nll).migrad();
			profileGraphs[pdf_name]->SetPoint(p,val,2.*nll->getVal());
			p++;
			// check where best fit value is in comparison to scan
			if (bestFitVal>val && bestFitVal<val+step_size) {
				profileGraphs[pdf_name]->SetPoint(p,bestFitVal,2.*bestFitNll);
				p++;
			}
			poi->setConstant(false);
		}
		if (printProgress) cout << endl;
		// if best fit val is > high scan point add to graph
		if (bestFitVal>=poiHigh){
				profileGraphs[pdf_name]->SetPoint(p,bestFitVal,2.*bestFitNll);
				p++;
		}
		cout << "Best fit at: x = " << bestFitVal << " 2nll = " << 2.*bestFitNll << endl;
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
}

void ProfileMultiplePdfs::addProfile(TGraph* graph, int npars){
	profileGraphsWithNPars.insert(make_pair(string(graph->GetName()),make_pair(graph,npars)));
}

void ProfileMultiplePdfs::printProfilesWithNPars(){
	
	cout << "List of profiles with n free params" << endl;
	for (map<string,pair<TGraph*,int> >::iterator it=profileGraphsWithNPars.begin(); it!=profileGraphsWithNPars.end(); it++){
		cout << it->first << " : " << it->second.first->GetName() << " -- " << it->second.second << " -- " << it->second.first->GetN() << endl;
	}
}

void ProfileMultiplePdfs::constructEnvelope(string ext){

	if (profileGraphsWithNPars.size()==0) {
		cerr << "No graphs to return envelope for! " << endl;
		exit(1);
	}
	else if (profileGraphsWithNPars.size()==1) {
		envelope_ = profileGraphsWithNPars.begin()->second.first;
	}
	else {
		// check graphs have same number of points
		// also find global best fit value for final graph
		envelope_ = new TGraph();
		envelope_->SetName(Form("envelope_c%3.1f%s",corr_pdof_,ext.c_str()));
		map<string,pair<TGraph*,int> >::iterator checker = profileGraphsWithNPars.begin();
		int np=checker->second.first->GetN();
		envelopeMinNll_ = 9999.;
		envelopeBestFitName_ = "notfound";
		envelopeBestFitVal_ = 999.;

		for (map<string,pair<TGraph*,int> >::iterator it=profileGraphsWithNPars.begin(); it!=profileGraphsWithNPars.end(); it++){
			TGraph *graph = it->second.first;
			int dof = it->second.second;
			if (graph->GetN()!=np) {
				cerr << "Graphs for envelope computation have different number of points. Can't deal with this. Sorry!" << endl;
				exit(1);
			}
			double x,y;
			for (int p=0; p<np; p++){
				graph->GetPoint(p,x,y);
				double corrNll = y+(corr_pdof_*dof);
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
			profileGraphsWithNPars[envelopeBestFitName_].first->GetPoint(p,x,y);
			double minNll = 9999.;
			for (map<string,pair<TGraph*,int> >::iterator it=profileGraphsWithNPars.begin(); it!=profileGraphsWithNPars.end(); it++){
				TGraph *graph = it->second.first;
				int npars = it->second.second;
				double correctedNll = graph->Eval(x)+(npars*corr_pdof_);
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

void ProfileMultiplePdfs::setCorrPdof(double corr){
	corr_pdof_ = corr;
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

void ProfileMultiplePdfs::drawEnvelope(string name, bool shiftToZero){
	TCanvas *canv = new TCanvas();
	TGraph* envelope = shiftGraph(envelope_,envelopeMinNll_);
	envelope->SetLineColor(kBlue);
	envelope->SetLineWidth(3);
	envelope->SetLineStyle(kDashed);
	envelope->Draw("ALP");
	for (map<string,pair<TGraph*,int> >::iterator it=profileGraphsWithNPars.begin(); it!=profileGraphsWithNPars.end(); it++){
		int dof = it->second.second;
		TGraph *graph = shiftGraph(it->second.first,envelopeMinNll_-dof*corr_pdof_);
		graph->SetLineColor(kRed);
		graph->SetLineWidth(2);
		graph->Draw("LPsame");
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
