void plotMuVsCorrection() {

   const double epsilon = 0.05;
  //gROOT->ProcessLine(".x paperStyle.C");

	const int nPoints = 8;
	double corr[nPoints] = {0., 0.5, 1., 1.5, 2., 2.5, 3., -1.};

	double mu[nPoints] 			= {1.1 	 			, 0.95 				, 0.95 			, 0.95 			, 0.95 			, 0.95 			, 0.95 			,	0.9 	 			, 	};
	double muLow1[nPoints] 	= {-0.600908 	, -0.527051 	, -0.560094 , -0.560094 , -0.560094 , -0.560094	, -0.560094	,	-0.558575 	, 	};
	double muHigh1[nPoints] = {0.597126 	, 0.67814 		, 0.473881 	, 0.473881 	, 0.473881 	, 0.473881	, 0.473881	,	0.373674		, 	};
	double muLow2[nPoints] 	= {-1.18148 	, -1.05678 		, -1.11012 	, -1.11012	, -1.11012	, -1.11012	, -1.11012	,	-1.10778		, 	};
	double muHigh2[nPoints] = {1.21183 		, 1.32824 		, 1.18089 	, 0.98518		, 0.98518		, 0.98518		, 0.98518		,	0.841168		, };

	TH1F *bestFitH = new TH1F("bestFitH","",nPoints,0,nPoints);
	bestFitH->SetStats(0);
	bestFitH->SetLineColor(kRed);
	bestFitH->SetLineWidth(3);
	for (int p=0; p<nPoints; p++){
		TString name = Form("%3.1f / dof",corr[p]);;
		if (p==nPoints-1) name = "p-value";
		if (p==2) name += " = app. p-value";
		if (p==4) name += " = Akaike";
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
		bestFit->SetPointError(p,bestFitH->GetBinWidth(p+1)/2.-epsilon,bestFitH->GetBinWidth(p+1)/2.-epsilon,0.,0.);
		err1sigma->SetPointError(p,bestFitH->GetBinWidth(p+1)/2.-epsilon,bestFitH->GetBinWidth(p+1)/2.-epsilon,TMath::Abs(muLow1[p]),muHigh1[p]);
		err2sigma->SetPointError(p,bestFitH->GetBinWidth(p+1)/2.-epsilon,bestFitH->GetBinWidth(p+1)/2.-epsilon,TMath::Abs(muLow2[p]),muHigh2[p]);
	}

	TCanvas *canv = new TCanvas("c","c",800,700);
	canv->SetBottomMargin(0.2);

	TLegend *leg = new TLegend(0.6,0.7,0.89,0.89);
	leg->SetFillColor(0);
	leg->SetLineColor(0);
	leg->AddEntry(bestFitH,"Fit value","L");
	leg->AddEntry(err1sigma,"68.3% Interval","LF");
	leg->AddEntry(err2sigma,"95.4% Interval","LF");

	bestFitH->GetYaxis()->SetRangeUser(-1.,4.);
	bestFitH->GetXaxis()->LabelsOption("d");
	bestFitH->GetXaxis()->SetLabelOffset(0.01);
	bestFitH->GetXaxis()->SetLabelSize(0.05);
	bestFitH->GetXaxis()->SetTitle("#Lambda correction");
	bestFitH->GetXaxis()->SetTitleSize(0.055);
	bestFitH->GetXaxis()->SetTitleOffset(1.7);
	bestFitH->SetMarkerSize(0);
	bestFitH->GetYaxis()->SetTitle("#mu");
	bestFitH->GetYaxis()->SetTitleSize(0.045);
	bestFitH->GetYaxis()->SetLabelSize(0.04);

	TLine *sepLine = new TLine();
	sepLine->SetLineWidth(3);
	sepLine->SetLineStyle(7);
	sepLine->SetLineColor(kBlack);

	bestFitH->Draw("AXIS");
	err2sigma->Draw("2same");
	err1sigma->Draw("2same");
	bestFit->Draw("Esame");
	//bestFitH->Draw("HIST ][ same");
	sepLine->DrawLine(nPoints-1,bestFitH->GetMinimum(),nPoints-1, bestFitH->GetMaximum());
	leg->Draw();
	canv->RedrawAxis();
	canv->Print("../correction/correction.pdf");




}
