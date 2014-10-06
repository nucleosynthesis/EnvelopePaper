#include "TCanvas.h"
#include "TMath.h"
#include "TIterator.h"
#include "TAxis.h"

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
  correctionValue(0.),
  isPValueCorrection(false),
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

void PdfModelBuilder::setPValueCorrection(bool arg){
  isPValueCorrection=arg;
}

void PdfModelBuilder::setCorrectionValue(double value){
  correctionValue = value;
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
  cachedNllVals.insert(make_pair(pdf->GetName(),999999.));

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
  cachedNllVals.insert(make_pair(pdf->GetName(),999999.));

}

double PdfModelBuilder::getFractionValue(string name){

  // seems a bit silly to re-calculate this on the fly every time
  // but I can't be bothered to setup a caching system for it

  double globalMinOverFunctions=999999.;
  for (map<string,double>::iterator it=cachedNllVals.begin(); it!=cachedNllVals.end(); it++){
    if ( it->second >= 999999. ){
      cout << "ERROR -- it seems the cachedNllVal still has it's initial value. No good." << endl;
      cout << it->first << " " << it->second << endl;
      for (map<string,double>::iterator mapIt=cachedNllVals.begin(); mapIt!=cachedNllVals.end(); mapIt++){
        cout << "\t" << mapIt->first << " -- " << mapIt->second << endl;
      }
      exit(1);
    }

    if ( it->second < globalMinOverFunctions ){
      globalMinOverFunctions = it->second;
    }
  }

  // now get each L value (relative to the difference) and make the sum
  // in other words we are now calculating exp(-0.5*(-2DeltaLL))
  map<string,double> llRelVals;
  double nllSum = 0.;
  for (map<string,double>::iterator it=cachedNllVals.begin(); it!=cachedNllVals.end(); it++){
    double dll = it->second - globalMinOverFunctions;
    double llRel = TMath::Exp(-0.5*dll);
    nllSum += llRel;
    llRelVals.insert(make_pair(it->first,llRel));
    //cout << "ahem: " << it->first << " " << it->second << " " << dll << " " << llRel << " " << nllSum << endl;
  }

  double relval = 0.;
  if (llRelVals.find(name) == llRelVals.end() ) {
    if (llRelVals.find("sb_"+name) == llRelVals.end() ){
      cout << "ERROR - couldn't find a cached value for " << name << endl;
    }
    else {
      relval = llRelVals["sb_"+name];
    }
  }
  else {
    relval = llRelVals[name];
  }

  double frac = relval/nllSum;

  if (frac < 0. || frac > 1.){
    cout << "ERROR -- the fraction for the frequentist toy somehow is not in the range [0,1] " << endl;
    exit(1);
  }

  cout << "Profiled gen frac: " << name << " = " << frac << endl;
  return frac;
}

RooAbsPdf *PdfModelBuilder::returnProfiledBackground(string name){

	RooArgList *functions = new RooArgList();
	RooArgList *fractions = new RooArgList();
  double fracSum=0.;
	for (map<string,RooAbsPdf*>::iterator it=bkgPdfs.begin(); it!=bkgPdfs.end(); it++){
		functions->add(*(it->second));
		RooRealVar *var = new RooRealVar(Form("f%d",fractions->getSize()),"f",getFractionValue(it->first));
		var->setConstant(true);
		fractions->add(*var);
    fracSum += var->getVal();
    // also have to set the background model parameters to those from the SB best fit
    RooArgSet *bkgPars = it->second->getParameters(*obs_var);
    RooArgSet *sbBkgPars = sbPdfs[string("sb_"+it->first)]->getParameters(*obs_var);
    RooRealVar *p;
    TIterator *iter = bkgPars->createIterator();
    while ( (p = (RooRealVar*)iter->Next()) ) {
      p->setVal(sbBkgPars->getRealValue(p->GetName(),-999999.));
      p->setConstant(true);
      if (p->getVal() <= -999999){
        cout << "ERROR -- something went wrong with parameter exchange!" << endl;
        exit(1);
      }
    }
    //cout << "---------- this bit: --------------" << endl;
    //bkgPars->Print("v");
    //sbBkgPars->Print("v");
    //cout << "-----------------------------------" << endl;
  }
  // would hope the fraction sum is bloody close to 1
  if (TMath::Abs(1.-fracSum)>1.e-3) {
    cout << "ERROR - the sum of the fractions for the profiled pdf is not very close to 1! sumF = " << fracSum << endl;
    exit(1);
  }
  functions->Print();
  fractions->Print();
	RooAddPdf *theNewPdf = new RooAddPdf(name.c_str(),name.c_str(),*functions,*fractions);
  theNewPdf->getParameters(*obs_var)->Print("v");
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

  // clear cache nll map
  cachedNllVals.clear();

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
    cachedNllVals.insert(make_pair(sbMod->GetName(),999999.));
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
   RooArgSet *bkgPars;
   if (bkgOnly) {
     bkgPars = it->second->getParameters(*obs_var);
   }
   else {
      bkgPars = bkgPdfs[it->first.substr(3,string::npos)]->getParameters(*obs_var);
      bkgPars->add(*bkgYield);
   }
	 //it->second->paramOn(plot,RooFit::Layout(0.35,1.,1.),RooFit::Format("NEA",FixedPrecision(3)),RooFit::Parameters(*bkgPars),RooFit::ShowConstants(true));
    }
    plot->SetTitle("");
    plot->GetYaxis()->SetRangeUser(0.,plot->GetMaximum()*1.5);
    plot->Draw();
    canv->Print(Form("%s_%s.pdf",name.c_str(),it->first.c_str()));
    canv->Print(Form("%s_%s.png",name.c_str(),it->first.c_str()));
  }
  delete canv;
}

double PdfModelBuilder::getChisq(RooAbsData *dat, RooAbsPdf *pdf, bool prt) {

  //RooRealVar *var = obs_var;
  RooRealVar *var = wsCache->var(obs_var->GetName());
  // Find total number of events
  double nEvt;
  double nTot=0.0;


  for(int j=0;j<dat->numEntries();j++) {
    dat->get(j);
    nEvt=dat->weight();
    nTot+=nEvt;
  }
  //RooAddPdf *pdf = dynamic_cast<RooAddPdf*>(pdfin);
  // Find chi-squared equivalent 2NLL
  double totNLL=0.0;
  double prbSum=0.0;

//  if(prt) cout << "NData " << dat.sumEntries() << endl;
  for(int j=0;j<dat->numEntries();j++) {
    double m=dat->get(j)->getRealValue(var->GetName());

    if ( m < var->getMin() || m > var->getMax())  continue;
    // Find probability density and hence probability
    var->setVal(m);
    double prb = var->getBinWidth(0) * pdf->getVal(*var);
    prbSum+=prb;

    dat->get(j);
    nEvt=dat->weight();

    double mubin=nTot*prb;
    double contrib(0.);
    if (nEvt < 1) contrib = mubin;
    else contrib=mubin-nEvt+nEvt*TMath::Log(nEvt/mubin);
    totNLL+=contrib;

    if(prt) cout << pdf->GetName() << ", Bin " << j << " center = " << m << " prob = " << prb << " nEvt = " << nEvt << ", mu = " << mubin << " contribution " << contrib << endl;
  }

  totNLL*=2.0;
  if(prt) cout << pdf->GetName() << " nTot = " << nTot << " 2NLL constant = " << totNLL << endl;

  return totNLL;
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
      if ( cachedNllVals.find(it->first)==cachedNllVals.end() ) {
        cout << "ERROR -- cachedNllVals map has not been properly initialised. Something went wrong. " << endl;
        cout << it->first << endl;
        for (map<string,double>::iterator mapIt=cachedNllVals.begin(); mapIt!=cachedNllVals.end(); mapIt++){
          cout << "\t" << mapIt->first << " -- " << mapIt->second << endl;
        }
        exit(1);
      }
      RooArgSet *bkgPars = (RooArgSet*)it->second->getParameters(*obs_var);
      if (!bkgOnly) {
        RooArgSet *comps = it->second->getComponents();
        bkgPars = comps->find(it->first.substr(3,string::npos).c_str())->getParameters(*obs_var);
      }

      double nll = getChisq(data,wsCache->pdf(it->second->GetName()));
      // correct
      double correction = correctionValue*bkgPars->getSize();
      if (isPValueCorrection) {
        double pvalue = TMath::Prob(nll,obs_var->getBins()-bkgPars->getSize());
        if (pvalue < 1.e-16) pvalue = 1.e-16;
        double nllEquiv = TMath::ChisquareQuantile(1.-pvalue,obs_var->getBins());
        correction = nllEquiv-nll;
      }
      //cout << "HERE! nll=" << nll << " nbins=" << obs_var->getBins() << " npars=" << bkgPars->getSize() << " corr=" << correction << endl;

      cachedNllVals[it->first] = nll + correction;
      cout << "Cached NLL value: " << it->first << " -- " << cachedNllVals[it->first] << endl;
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
         RooArgSet *bkgPars;
         if (bkgOnly) {
           bkgPars = pdfIt->second->getParameters(*obs_var);
         }
         else {
            bkgPars = bkgPdfs[pdfIt->first.substr(3,string::npos)]->getParameters(*obs_var);
            bkgPars->add(*bkgYield);
         }
         //pdfIt->second->paramOn(plot,RooFit::Layout(0.35,1.,1.),RooFit::Format("NEA",FixedPrecision(3)),RooFit::Parameters(*bkgPars),RooFit::ShowConstants(true));
        plot->SetTitle("");
        plot->GetYaxis()->SetRangeUser(0.,plot->GetMaximum()*1.5);
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
