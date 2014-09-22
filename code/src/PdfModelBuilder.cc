#include "TCanvas.h"

#include "RooPlot.h"
#include "RooBernstein.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
#include "RooExponential.h"
#include "RooKeysPdf.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooFitResult.h"
#include "RooRandom.h"

#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"

#include "../interface/PdfModelBuilder.h"

using namespace std;
using namespace RooFit;
using namespace boost;

PdfModelBuilder::PdfModelBuilder():
  obs_var_set(false),
  signal_modifier_set(false),
  signal_set(false),
  bkgHasFit(false),
  sbHasFit(false),
  verbosity(0)
{
  
  wsCache = new RooWorkspace("PdfModelBuilderCache");

};

PdfModelBuilder::~PdfModelBuilder(){};

void PdfModelBuilder::setObsVar(RooRealVar *var){
  obs_var=var;
  obs_var_set=true;
}

void PdfModelBuilder::setSignalModifier(RooRealVar *var){
  signalModifier=var;
  signal_modifier_set=true;
}

void PdfModelBuilder::setSignalModifierVal(float val){
  signalModifier->setVal(val);
}

void PdfModelBuilder::setSignalModifierConstant(bool val){
  signalModifier->setConstant(val);
}

RooAbsPdf* PdfModelBuilder::getPdfFromFile(string &prefix){
  vector<string> details;
  split(details,prefix,boost::is_any_of(","));

  string fname = details[2];
  string wsname = details[1];
  string pdfname = details[0];

  TFile *tempFile = TFile::Open(fname.c_str());
  if (!tempFile){
    cerr << "PdfModelBuilder::getPdfFromFile -- file not found " << fname << endl;
    assert(0);
  }
  RooWorkspace *tempWS = (RooWorkspace*)tempFile->Get(wsname.c_str());
  if (!tempWS){
    cerr << "PdfModelBuilder::getPdfFromFile -- workspace not found " << wsname << endl;
    assert(0);
  }
  RooAbsPdf *tempPdf = (RooAbsPdf*)tempWS->pdf(pdfname.c_str());
  if (!tempPdf){
    cerr << "PdfModelBuilder::getPdfFromFile -- pdf not found " << pdfname << endl;
    assert(0);
  }
  prefix = pdfname;
  RooAbsPdf *pdf = (RooAbsPdf*)tempPdf->Clone(prefix.c_str());
  tempFile->Close();
  delete tempFile;
  return pdf;
}

void PdfModelBuilder::addBkgPdf(RooAbsPdf *pdf, bool cache){
  assert(pdf);
  if (cache) {
    wsCache->import(*pdf);
    RooAbsPdf *cachePdf = wsCache->pdf(pdf->GetName());
    bkgPdfs.insert(pair<string,RooAbsPdf*>(cachePdf->GetName(),cachePdf));
  }
  else {
    bkgPdfs.insert(pair<string,RooAbsPdf*>(pdf->GetName(),pdf));
  }
	
}

void PdfModelBuilder::addBkgPdf(string name, bool cache){
 
  if (!obs_var_set){
    cerr << "ERROR -- obs Var has not been set!" << endl;
    exit(1);
  }
  RooAbsPdf *pdf = getPdfFromFile(name);

  if (cache) {
    wsCache->import(*pdf);
    RooAbsPdf *cachePdf = wsCache->pdf(pdf->GetName());
    bkgPdfs.insert(pair<string,RooAbsPdf*>(cachePdf->GetName(),cachePdf));
  }
  else {
    bkgPdfs.insert(pair<string,RooAbsPdf*>(pdf->GetName(),pdf));
  }

}

double PdfModelBuilder::getFractionValue(string name){
	return 1.;
}

RooAbsPdf *PdfModelBuilder::returnProfiledBackground(string name){
	
	RooArgList *functions;
	RooArgList *fractions;
	for (map<string,RooAbsPdf*>::iterator it=bkgPdfs.begin(); it!=bkgPdfs.end(); it++){
		functions->add(*(it->second));
		RooRealVar var(Form("f%d",fractions->getSize(),"f",getFractionValue(it->first)));	
		var.setConstant(true);
		fractions->add(var);
	}
	RooAddPdf *theNewPdf = new RooAddPdf(name.c_str(),name.c_str(),*functions,*fractions);
	return theNewPdf;
}

void PdfModelBuilder::setSignalPdf(RooAbsPdf *pdf, RooRealVar *norm){
  sigPdf=pdf;
  sigNorm=norm;
  signal_set=true;
}

void PdfModelBuilder::setSignalPdfFromMC(RooDataSet *data){
  
  RooDataHist *sigMCBinned = new RooDataHist(Form("roohist_%s",data->GetName()),Form("roohist_%s",data->GetName()),RooArgSet(*obs_var),*data);
  sigPdf = new RooHistPdf(Form("pdf_%s",data->GetName()),Form("pdf_%s",data->GetName()),RooArgSet(*obs_var),*sigMCBinned);
  sigNorm = new RooConstVar(Form("sig_events_%s",data->GetName()),Form("sig_events_%s",data->GetName()),data->sumEntries());
  signal_set=true;
}

void PdfModelBuilder::makeSBPdfs(bool cache){
  
  if (!signal_set){
    cerr << "ERROR - no signal model set!" << endl;
    exit(1);
  }
  if (!signal_modifier_set){
    cerr << "ERROR - no signal modifier set!" << endl;
    exit(1);
  }
 
  if (sigNorm) {
    sigYield = new RooProduct("sig_yield","sig_yield",RooArgSet(*signalModifier,*sigNorm));
  }
  else {
    sigYield = signalModifier;
  }

  bkgYield = new RooRealVar("bkg_yield","bkg_yield",7712.,0.,1.e6);

  for (map<string,RooAbsPdf*>::iterator bkg=bkgPdfs.begin(); bkg!=bkgPdfs.end(); bkg++){
    RooAbsPdf *sbMod = new RooAddPdf(Form("sb_%s",bkg->first.c_str()),Form("sb_%s",bkg->first.c_str()),RooArgList(*(bkg->second),*sigPdf),RooArgList(*bkgYield,*sigYield));
    if (cache) {
      wsCache->import(*sbMod,RecycleConflictNodes());
      RooAbsPdf *cachePdf = (RooAbsPdf*)wsCache->pdf(sbMod->GetName());
      signalModifier = (RooRealVar*)wsCache->var(signalModifier->GetName());
      sbPdfs.insert(pair<string,RooAbsPdf*>(cachePdf->GetName(),cachePdf));
    }
    else {
      sbPdfs.insert(pair<string,RooAbsPdf*>(sbMod->GetName(),sbMod));
    }
  }
}

map<string,RooAbsPdf*> PdfModelBuilder::getBkgPdfs(){
  return bkgPdfs;
}

map<string,RooAbsPdf*> PdfModelBuilder::getSBPdfs(){
  return sbPdfs;
}

RooAbsPdf* PdfModelBuilder::getSigPdf(){
  return sigPdf;
}

void PdfModelBuilder::plotPdfsToData(RooAbsData *data, int binning, string name, bool bkgOnly,string specificPdfName){
  
  TCanvas *canv = new TCanvas();
  bool specPdf=false;
  if (specificPdfName!="") specPdf=true;

  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) pdfSet = bkgPdfs;
  else pdfSet = sbPdfs;
  
  for (map<string,RooAbsPdf*>::iterator it=pdfSet.begin(); it!=pdfSet.end(); it++){
    if (specPdf && it->first!=specificPdfName && specificPdfName!="NONE") continue;
    RooPlot *plot = obs_var->frame();
    data->plotOn(plot,Binning(binning));
    if (specificPdfName!="NONE") {
	 it->second->plotOn(plot);
	 it->second->paramOn(plot,RooFit::Layout(0.34,0.96,0.89),RooFit::Format("NEA",AutoPrecision(1)));
    }	
    plot->Draw();
    canv->Print(Form("%s_%s.pdf",name.c_str(),it->first.c_str()));
    canv->Print(Form("%s_%s.png",name.c_str(),it->first.c_str()));
  }
  delete canv;
}

void PdfModelBuilder::fitToData(RooAbsData *data, bool bkgOnly, bool cache, bool print){
  
  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) pdfSet = bkgPdfs;
  else pdfSet = sbPdfs;

  for (map<string,RooAbsPdf*>::iterator it=pdfSet.begin(); it!=pdfSet.end(); it++){
    RooFitResult *fit = (RooFitResult*)it->second->fitTo(*data,Save(true));
    if (print){
      cout << "Fit Res Before: " << endl;
      fit->floatParsInit().Print("v");
      cout << "Fit Res After: " << endl;
      fit->floatParsFinal().Print("v");
    }
    if (cache) {
      RooArgSet *fitargs = (RooArgSet*)it->second->getParameters(*obs_var);
      // remove the signal strength since this will be set AFTER fitting the background 
      fitargs->remove(*signalModifier); 
      wsCache->defineSet(Form("%s_params",it->first.c_str()),*fitargs);
      wsCache->defineSet(Form("%s_observs",it->first.c_str()),*obs_var);
      wsCache->saveSnapshot(it->first.c_str(),*fitargs,true);
      if (print) {
        cout << "Cached values: " << endl;
        fitargs->Print("v");
      }
    }
  }
  if (bkgOnly) bkgHasFit=true;
  else sbHasFit=true;
}

void PdfModelBuilder::setSeed(int seed){
  RooRandom::randomGenerator()->SetSeed(seed);
}

void PdfModelBuilder::throwToy(string postfix, int nEvents, bool bkgOnly, bool binned, bool poisson, bool cache){

  toyData.clear();
  toyDataSet.clear();
  toyDataHist.clear();
  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) {
    pdfSet = bkgPdfs;
    if (!bkgHasFit) cerr << "WARNING -- bkg has not been fit to data. Are you sure this is wise?" << endl; 
  }
  else {
    pdfSet = sbPdfs;
    if (!sbHasFit) cerr << "WARNING -- sb has not been fit to data. Are you sure this is wise?" << endl;
  }
  
  for (map<string,RooAbsPdf*>::iterator it=pdfSet.begin(); it!=pdfSet.end(); it++){
    if (cache) {
      wsCache->loadSnapshot(it->first.c_str());
      cout << "Loaded snapshot, params at.." << endl;
      it->second->getVariables()->Print("v");
    }
    RooAbsData *toy;
    if (binned){
      RooDataHist *toyHist;
      if (poisson) toyHist = it->second->generateBinned(RooArgSet(*obs_var),nEvents,Extended(),Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      else toyHist = it->second->generateBinned(RooArgSet(*obs_var),nEvents,Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      toyDataHist.insert(pair<string,RooDataHist*>(toyHist->GetName(),toyHist));
      toy=toyHist;
    }
    else {
      RooDataSet *toySet;
      if (poisson) toySet = it->second->generate(RooArgSet(*obs_var),nEvents,Extended(),Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      else toySet = it->second->generate(RooArgSet(*obs_var),nEvents,Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      toyDataSet.insert(pair<string,RooDataSet*>(toySet->GetName(),toySet));
      toy=toySet;
    }
    toyData.insert(pair<string,RooAbsData*>(toy->GetName(),toy));
  }
  
}

map<string,RooAbsData*> PdfModelBuilder::getToyData(){
  return toyData;
}

RooAbsData* PdfModelBuilder::getToyDataSingle(){
	if (toyData.size()!=1) {
		cerr << "Terrible way of coding this I know but you can only call PdfModelBuild::getToyDataSingle() when the size of the toyData map has one entry. Sorry!" << endl;
		exit(1);
	}
	return toyData.begin()->second;
}

void PdfModelBuilder::plotToysWithPdfs(string prefix, int binning, bool bkgOnly){
  
  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) {
    pdfSet = bkgPdfs;
  }
  else {
    pdfSet = sbPdfs;
  }
  TCanvas *canv = new TCanvas();
  for (map<string,RooAbsPdf*>::iterator pdfIt = pdfSet.begin(); pdfIt != pdfSet.end(); pdfIt++){
    for (map<string,RooAbsData*>::iterator toyIt = toyData.begin(); toyIt != toyData.end(); toyIt++){
      //cout << "pdf: " << pdfIt->first << " - toy: " << toyIt->first << endl; 
      if (toyIt->first.find(pdfIt->first)!=string::npos){
        RooPlot *plot = obs_var->frame();
        toyIt->second->plotOn(plot,Binning(binning));
        pdfIt->second->plotOn(plot,LineColor(kRed));
        pdfIt->second->paramOn(plot,LineColor(kRed),RooFit::Layout(0.34,0.96,0.89),RooFit::Format("NEA",AutoPrecision(1)));
        plot->Draw();
        canv->Print(Form("%s_%s.pdf",prefix.c_str(),pdfIt->first.c_str()));
        canv->Print(Form("%s_%s.png",prefix.c_str(),pdfIt->first.c_str()));
      }
    }
  }
  delete canv;

}

void PdfModelBuilder::saveWorkspace(string filename){
  
  TFile *outFile = new TFile(filename.c_str(),"RECREATE");
  outFile->cd();
  wsCache->Write();
  outFile->Close();
  delete outFile;
}

void PdfModelBuilder::saveWorkspace(TFile *file){
  file->cd();
  wsCache->Write();
}
