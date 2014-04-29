void deltaChiSq() {
  gROOT->ProcessLine(".x paperStyle.C");
  //gROOT->SetStyle("Plain");
  TCanvas *canv(new TCanvas());
  gStyle->SetOptStat(0);
  //gStyle->SetPalette(1);

  //TFile *file=TFile::Open("RooDumb1249.root");
  TH1F *h1;
  TH2F *h2;
  TH1D *d1;

  TLine a;
  a.SetLineStyle(2);

  TLatex t;

  // Chi-sq and probability
  TH1F *fDeltaChi2[10];
  TH1F *hDeltaChi2[10];

  for(unsigned j(1);j<10;j++) {
    std::ostringstream sout;
    sout << 320-j;
    fDeltaChi2[j]=new TH1F((std::string("FDeltaChi2_320_")+sout.str()).c_str(),
                           (std::string(";p-value;Change of #chi^{2} for change of DoF to 320")).c_str(),
                           10000,0.0,1.0);
    
    hDeltaChi2[j]=new TH1F((std::string("HDeltaChi2_320_")+sout.str()).c_str(),
                           (std::string(";Change of #chi^{2} for change of DoF to 320;Arbitrary units")).c_str(),
                           1000,0.0,10.0);
    
    for(unsigned i(0);i<10000;i++) {
      double p=0.0001*(i+0.5);
      double c0=TMath::ChisquareQuantile(1-p,320);
      double c1=TMath::ChisquareQuantile(1-p,320-j);
      fDeltaChi2[j]->SetBinContent(i+1,(c0-c1));
      hDeltaChi2[j]->Fill(c0-c1);
    }
  }

  fDeltaChi2[5]->SetLineColor(1);
  fDeltaChi2[5]->GetYaxis()->SetRangeUser(0.0,6.0);
  fDeltaChi2[5]->Draw();
  fDeltaChi2[4]->SetLineColor(2);
  fDeltaChi2[4]->Draw("same");
  fDeltaChi2[3]->SetLineColor(3);
  fDeltaChi2[3]->Draw("same");
  fDeltaChi2[2]->SetLineColor(4);
  fDeltaChi2[2]->Draw("same");
  fDeltaChi2[1]->SetLineColor(6);
  fDeltaChi2[1]->Draw("same");

  t.SetTextColor(6);
  t.DrawLatex(0.1,0.6,"Dof = 319");
  t.SetTextColor(4);
  t.DrawLatex(0.2,1.6,"Dof = 318");
  t.SetTextColor(3);
  t.DrawLatex(0.3,2.6,"Dof = 317");
  t.SetTextColor(2);
  t.DrawLatex(0.4,3.6,"Dof = 316");
  t.SetTextColor(1);
  t.DrawLatex(0.5,4.6,"Dof = 315");

  canv->Print("../correction/DeltaChiSq1.pdf");

  hDeltaChi2[1]->SetLineColor(6);
  hDeltaChi2[1]->GetXaxis()->SetRangeUser(0.0,6.0);
  hDeltaChi2[1]->Draw();
  hDeltaChi2[5]->SetLineColor(1);
  hDeltaChi2[5]->Draw("same");
  hDeltaChi2[4]->SetLineColor(2);
  hDeltaChi2[4]->Draw("same");
  hDeltaChi2[3]->SetLineColor(3);
  hDeltaChi2[3]->Draw("same");
  hDeltaChi2[2]->SetLineColor(4);
  hDeltaChi2[2]->Draw("same");  

  t.SetTextColor(6);
  t.DrawLatex(4.0,900.0,"Dof = 319");
  t.SetTextColor(4);
  t.DrawLatex(4.0,820.0,"Dof = 318");
  t.SetTextColor(3);
  t.DrawLatex(4.0,740.0,"Dof = 317");
  t.SetTextColor(2);
  t.DrawLatex(4.0,660.0,"Dof = 316");
  t.SetTextColor(1);
  t.DrawLatex(4.0,580.0,"Dof = 315");

  canv->Print("../correction/DeltaChiSq2.pdf");
}
