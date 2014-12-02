#include <iostream>
#include <string>
#include <vector>

#include "TH1F.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TF1.h"
#include "TMath.h"

#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"  
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"

using namespace std;
using namespace RooFit;

void envelope_cartoon(){
  gROOT->ProcessLine(".x paperStyle.C");
  // PLOT 1
  TH1F dummy("d","",1,120,130);
  dummy.GetYaxis()->SetRangeUser(0.,5.5);
  dummy.SetStats(0);
  dummy.GetXaxis()->SetTitle("x");
  dummy.GetYaxis()->SetTitle("#Delta#Lambda");

  TCanvas c("c","c",800,600);
  c.SetGrid();

  TLegend *leg = new TLegend(0.2,0.65,0.8,0.89);
  leg->SetFillColor(0);

  dummy.Draw("G");

  TLine l1(120,1,130,1);
  l1.SetLineWidth(1);
  l1.SetLineColor(kBlack);
  l1.SetLineStyle(kDashed);

  TLine l2(120,4,130,4);
  l2.SetLineWidth(1);
  l2.SetLineColor(kBlack);
  l2.SetLineStyle(kDashed);

  
  TF1 f("f1","(0.5*(x-125))*(0.5*(x-125))",120,130);
  f.SetLineColor(kBlack);
  f.SetLineWidth(3);
  leg->AddEntry(&f,"Full profile fit","L");

  TF1 f2("f2","(1.5*(x-125))*(1.5*(x-125))",120,130);
  f2.SetLineColor(kBlue);
  f2.SetLineWidth(3);
  leg->AddEntry(&f2,"Fit fixing nuisance parameter to best fit","L");

  
  TGraph *envelope = new TGraph();
  envelope->SetLineWidth(3);
  envelope->SetLineStyle(4);
  envelope->SetLineColor(kGreen);
  int npoints = 500;
  double x,y;
  for (int p=0; p<npoints; p++){
    double mh = 120.+p*((130.-120.)/double(npoints));
    envelope->SetPoint(p,mh,100.);
  }


  for (int m=122; m<=128; m+=1){
    double yval = f.Eval(m);
    yval *= 1.+0.05*TMath::Abs(m-125.); 
    TF1 *fS = new TF1(Form("fS_%d",m),Form("(1.5*(x-%d))*(1.5*(x-%d))+%6.3f",m,m,yval),120,130);
    fS->SetLineColor(kRed);
    fS->SetLineStyle(kDashed);
    fS->Draw("Lsame");
    if (m==122) {
      leg->AddEntry(fS,"Fits fixing nuisance parameter to arbitary values","L");
    }
    for (int p=0; p<npoints; p++){
      envelope->GetPoint(p,x,y);
      if (fS->Eval(x)<y) envelope->SetPoint(p,x,fS->Eval(x));
    }
  
  }
  leg->AddEntry(envelope,"Envelope","L");
  envelope->Draw("Csame");
  f2.Draw("Lsame");
  f.Draw("Lsame");
  leg->Draw("same");
  c.Update();
  c.Modified();
  c.Print("../concept/envelope_cartoon.pdf");
}
