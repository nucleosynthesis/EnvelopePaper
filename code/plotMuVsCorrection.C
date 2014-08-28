void plotMuVsCorrection() {

  gROOT->ProcessLine(".x paperStyle.C");

	const int nPoints = 8;
	double corr[nPoints] = {-1., 0., 0.5, 1., 1.5, 2., 2.5, 3.};

	double mu[nPoints] 			= {0.9 	 			, 1.1 	 			, 0.95 				, 0.95 			, 0.95 			, 0.95 			, 0.95 			, 0.95 				};
	double muLow1[nPoints] 	= {-0.558575 	, -0.600908 	, -0.527051 	, -0.560094 , -0.560094 , -0.560094 , -0.560094	, -0.560094		};
	double muHigh1[nPoints] = {0.373674		, 0.597126 	, 0.67814 		, 0.473881 	, 0.473881 	, 0.473881 	, 0.473881	, 0.473881		};
	double muLow2[nPoints] 	= {-1.10778		, -1.18148 	, -1.05678 		, -1.11012 	, -1.11012	, -1.11012	, -1.11012	, -1.11012		};
	double muHigh2[nPoints] = {0.841168		, 1.21183 		, 1.32824 		, 1.18089 	, 0.98518		, 0.98518		, 0.98518		, 0.98518		};

	TH1F *bestFitH = new TH1F("bestFitH","",nPoints,0,nPoints);
	bestFitH->SetLineColor(kRed);
	bestFitH->SetLineWidth(3);
	for (int p=0; p<nPoints; p++){
		TString name = "";
		if (p==0) name = "p-value";
		else name = Form("%3.1f / dof",corr[p]);
		bestFitH->GetXaxis()->SetBinLabel(p+1,name.Data());
	}

	TGraphAsymmErrors *bestFit = new TGraphAsymmErrors();
	bestFit->SetLineColor(kRed);
	bestFit->SetLineWidth(4);
	bestFit->SetMarkerColor(kRed);
	bestFit->SetMarkerSize(0);

	TGraphAsymmErrors *err1sigma = new TGraphAsymmErrors();
	err1sigma->SetName("err1sigma");
	err1sigma->SetLineColor(kGreen+3);
	err1sigma->SetFillColor(kGreen+3);

	TGraphAsymmErrors *err2sigma = new TGraphAsymmErrors();
	err2sigma->SetName("err2sigma");
	err2sigma->SetLineColor(kYellow);
	err2sigma->SetFillColor(kYellow);

	for (int p=0; p<nPoints; p++){
		bestFitH->SetBinContent(p+1,mu[p]);
		bestFitH->SetBinError(p+1,0.);
		bestFit->SetPoint(p,bestFitH->GetBinCenter(p+1),mu[p]);
		err1sigma->SetPoint(p,bestFitH->GetBinCenter(p+1),mu[p]);
		err2sigma->SetPoint(p,bestFitH->GetBinCenter(p+1),mu[p]);
		bestFit->SetPointError(p,bestFitH->GetBinWidth(p+1)/2.,bestFitH->GetBinWidth(p+1)/2.,0.,0.);
		err1sigma->SetPointError(p,bestFitH->GetBinWidth(p+1)/2.,bestFitH->GetBinWidth(p+1)/2.,TMath::Abs(muLow1[p]),muHigh1[p]);
		err2sigma->SetPointError(p,bestFitH->GetBinWidth(p+1)/2.,bestFitH->GetBinWidth(p+1)/2.,TMath::Abs(muLow2[p]),muHigh2[p]);
	}

	TCanvas *canv = new TCanvas();

	TLegend *leg = new TLegend(0.6,0.7,0.89,0.89);
	leg->SetFillColor(0);
	leg->SetLineColor(0);
	leg->AddEntry(bestFitH,"Fit value","L");
	leg->AddEntry(err1sigma,"#pm 1#sigma error","LF");
	leg->AddEntry(err2sigma,"#pm 2#sigma error","LF");

	bestFitH->GetYaxis()->SetRangeUser(-1.,4.);
	bestFitH->GetXaxis()->LabelsOption("d");
	bestFitH->GetXaxis()->SetLabelSize(0.05);
	bestFitH->GetYaxis()->SetTitle("#mu");
	bestFitH->GetXaxis()->SetTitle("#Lambda correction");
	bestFitH->GetXaxis()->SetTitleOffset(1.2);
	bestFitH->SetMarkerSize(0);

	bestFitH->Draw("AXIS");
	err2sigma->Draw("2same");
	err1sigma->Draw("2same");
	bestFit->Draw("EPsame");
	//bestFitH->Draw("HIST ][ same");
	leg->Draw();
	canv->Print("../correction/correction.pdf");




}
