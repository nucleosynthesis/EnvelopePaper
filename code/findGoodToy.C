#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooRandom.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooMsgService.h"

using namespace std;
using namespace RooFit;

void findGoodToy(){
 
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setSilentMode(true);

  TFile *fi = TFile::Open("envelopews.root");
  RooWorkspace *multipdf = (RooWorkspace*)fi->Get("multipdf");
  RooRealVar *x = (RooRealVar*)multipdf->var("CMS_hgg_mass");
  x->setBins(320);

  TFile *outf = new TFile("new_toys.root","RECREATE");
  RooWorkspace *outws = new RooWorkspace("newtoys");

  float approxMuError = 0.1;

  RooDataHist *datatoy = (RooDataHist*)multipdf->data("roohist_data_mass_cat1_toy__CMS_hgg_mass");

  // get pdfs 
  RooAbsPdf *bern_pdf = multipdf->pdf("env_pdf_1_8TeV_bern1");
  RooAbsPdf *exp_pdf = multipdf->pdf("env_pdf_1_8TeV_exp1");
  RooAbsPdf *pow_pdf = multipdf->pdf("env_pdf_1_8TeV_pow1");
  RooAbsPdf *lau_pdf = multipdf->pdf("env_pdf_1_8TeV_lau1");

  map<string,RooAbsPdf*> bkg_pdfs;
  map<string,RooAbsPdf*> sb_pdfs;

  bkg_pdfs.insert(pair<string,RooAbsPdf*>("bern",bern_pdf));
  bkg_pdfs.insert(pair<string,RooAbsPdf*>("exp",exp_pdf));
  bkg_pdfs.insert(pair<string,RooAbsPdf*>("pow",pow_pdf));
  bkg_pdfs.insert(pair<string,RooAbsPdf*>("lau",lau_pdf));
   
  // Build signal model
  RooRealVar mean("mean","mean",125); 
  mean.setConstant();
  RooRealVar sigma("sigma","sigma",1.19); 
  sigma.setConstant();
  RooRealVar nsignal_const("nsignal_const","nsignal_const",50.8); 
  nsignal_const.setConstant();

  RooGaussian sig_pdf("gaus","gaus",*x,mean,sigma);
  RooRealVar nbkg("nbkg","nbkg",datatoy->sumEntries(),0,10E6);

  RooRealVar mu("r","r",1,-10,10);
  RooFormulaVar nsig("nsig","nsig","@0*@1",RooArgList(mu,nsignal_const));
  mu.setVal(0.);
  mu.setConstant();

  for (map<string,RooAbsPdf*>::iterator pdf=bkg_pdfs.begin(); pdf!=bkg_pdfs.end(); pdf++){
    
    RooAddPdf *sbpdf = new RooAddPdf(Form("sbpdf_%s",pdf->first.c_str()),Form("sbpdf_%s",pdf->first.c_str()),RooArgList(sig_pdf,*(pdf->second)),RooArgList(nsig,nbkg));
    // fit them to our pseduo data thing first
    sbpdf->fitTo(*datatoy);
    sb_pdfs.insert(make_pair(pdf->first,sbpdf));
  }

  // generate exp asimov
  RooDataHist *exp_asimov = (RooDataHist*)bkg_pdfs["exp"]->generateBinned(*x,datatoy->sumEntries(),Asimov());
  // now fit asimov with pow so it looks close
  bkg_pdfs["pow"]->fitTo(*exp_asimov);

  int close_count=0;
  int toyn=0;
  mu.setConstant(false);
  while (close_count<3) {
  //for (int toyn=1; toyn<=1; toyn++){
    toyn++;
    if (toyn%5000==0) cout << "Tried " << toyn << " toys" << endl;

    RooRandom::randomGenerator()->SetSeed(toyn);
    RooDataHist *new_toy = (RooDataHist*)bkg_pdfs["pow"]->generateBinned(*x,datatoy->sumEntries(),Extended());

    map<string,double> fitMuStash;
    map<string,double> fitNllStash;
    double bestFitMu=-999.;
    double globalMinNll=1000.;

    // now fit toy with the four models and try and find one which has what we want
    for (map<string,RooAbsPdf*>::iterator pdf=sb_pdfs.begin(); pdf!=sb_pdfs.end(); pdf++){
      
      RooFitResult *fitRes = pdf->second->fitTo(*new_toy,Save(true));
      double fittedMu = mu.getVal();
      double fittedNll = fitRes->minNll();
      delete fitRes;
      
      fitMuStash.insert(make_pair(pdf->first,fittedMu));
      fitNllStash.insert(make_pair(pdf->first,fittedNll));

      if (fittedNll<globalMinNll){
        globalMinNll = fittedNll;
        bestFitMu = fittedMu;
      }
    }

    // now ascertain if several pdfs are close
    close_count=0;
    for (map<string,double>::iterator muIt=fitMuStash.begin(); muIt!=fitMuStash.end(); muIt++){
      
      double fittedMu = muIt->second;
      double fittedNll = fitNllStash[muIt->first];
      
      if (TMath::Abs(fittedNll-globalMinNll)<0.5 && TMath::Abs(fittedMu-bestFitMu)>approxMuError) { // i.e. close in likelihood but far in mu to give overlapped envelope
        close_count++; // expect at least 1 of course 
      }
    
    }

    if (close_count>=3) {
      cout << "Found with seed " << toyn << endl;
      new_toy->SetName(Form("new_toy_seed%d",toyn));
      outws->import(*new_toy);

      // lets scan the toy
      mu.setConstant(false);
      map<string,TGraph*> nll_graphs;
      for (map<string,RooAbsPdf*>::iterator pdfIt=sb_pdfs.begin(); pdfIt!=sb_pdfs.end(); pdfIt++){
        
        nll_graphs.insert(make_pair(pdfIt->first,new TGraph()));
        nll_graphs[pdfIt->first]->SetName(Form("%s_seed%d",pdfIt->first.c_str(),toyn));
        int p=0;
        for (double val=bestFitMu-2; val<bestFitMu+2.05; val+=0.1){
          
          mu.setVal(val);
          mu.setConstant(true);
          RooFitResult *fR = pdfIt->second->fitTo(*new_toy,Save(true));
          nll_graphs[pdfIt->first]->SetPoint(p,val,fR->minNll()*2.);
          delete fR;
          mu.setConstant(false);
          p++;
        }
        outf->cd();
        nll_graphs[pdfIt->first]->Write();
      }
      TCanvas *canv = new TCanvas();
      nll_graphs.begin()->second->Draw("ALP");
      for (map<string,TGraph*>::iterator graphIt=nll_graphs.begin(); graphIt!=nll_graphs.end(); graphIt++){
        graphIt->second->Draw("LPsame");
      }
      canv->Print("new_toy.pdf");
    }
  }
  
  outf->cd();
  outws->Write();
  outf->Close();
}
