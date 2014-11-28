{

  gROOT->SetBatch(1);
  gROOT->ProcessLine(".x paperStyle.C");
  TFile *fi = TFile::Open("symmetrized_errors_powerlawtoys_compare_corrections.root");


  TCanvas *can = new TCanvas();  
 
  err_hist_pow_c0->SetLineColor(kGreen+2);
  err_hist_pow_c0->SetLineWidth(2);
  err_hist_pow_c0->SetMarkerColor(kGreen+2);
  err_hist_pow_c0->SetMarkerStyle(20);
  err_hist_pow_c0->SetMarkerSize(0.8);

  err_hist_pow_c0->GetXaxis()->SetRangeUser(0.4,0.7);

  // Envelope
  err_hist_env_c1->SetLineWidth(2);
  err_hist_env_c1->SetLineColor(1);
  err_hist_env_c1->SetMarkerColor(1);
  err_hist_env_c1->SetMarkerStyle(23);
  err_hist_env_c1->SetMarkerSize(1.1);

  err_hist_env_cP->SetLineWidth(2);
  err_hist_env_cP->SetLineColor(kBlue);
  err_hist_env_cP->SetMarkerColor(kBlue);
  err_hist_env_cP->SetMarkerStyle(25);
  err_hist_env_cP->SetMarkerSize(0.7);
  
  err_hist_env_c2->SetLineWidth(2);
  err_hist_env_c2->SetLineColor(kMagenta);
  err_hist_env_c2->SetMarkerColor(kMagenta);
  err_hist_env_c2->SetMarkerStyle(28);
  err_hist_env_c2->SetMarkerSize(1.1);
 
  err_hist_pow_c0->Draw("pE");
  err_hist_env_c1->Draw("pEsame");
  err_hist_env_cP->Draw("pEsame");
  err_hist_env_c2->Draw("pEsame");

  err_hist_pow_c0->Draw("histsame");
  err_hist_env_c1->Draw("histsame");
  err_hist_env_cP->Draw("histsame");
  err_hist_env_c2->Draw("histsame");

  TLegend *leg = new TLegend(0.52,0.52,0.87,0.89);
  leg->SetFillColor(0);
  leg->AddEntry(err_hist_pow_c0,"fit power law","PLE");
  leg->AddEntry(err_hist_env_c1,"approx. p-value ","PLE");
  leg->AddEntry(err_hist_env_cP,"p-value ","PLE");
  leg->AddEntry(err_hist_env_c2,"Akaike ","PLE");
  leg->Draw();
  can->RedrawAxis();

  can->SaveAs("../correction/compare_error_magnitude.pdf");
}
