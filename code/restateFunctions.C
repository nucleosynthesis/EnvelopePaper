{

	TFile *fi = TFile::Open("envelopews_wsignal_toy1.root");
	RooWorkspace *w = (RooWorkspace*)fi->Get("multipdf");
	RooRealVar *x = w->var("CMS_hgg_mass");
	TFile *fo = new TFile("testws.root","RECREATE");
	RooWorkspace *wnew =(RooWorkspace*) w->Clone();
	 
	// Power laws with coeffs
	RooRealVar p1p1("env_pdf_1_8TeV_powc1_p1","env_pdf_1_8TeV_powc1_p1",-0.5,-10,0);
	RooGenericPdf pow1("env_pdf_1_8TeV_powc1","env_pdf_1_8TeV_powc1","TMath::Power(@0,@1)",RooArgList(*x,p1p1));

	// pow1, pow 3, pow 5
	RooRealVar p3v1("env_pdf_1_8TeV_powc3_f1","env_pdf_1_8TeV_powc3_f1",0.5,-100,100);
	RooRealVar p3p1("env_pdf_1_8TeV_powc3_p1","env_pdf_1_8TeV_powc3_p1",-0.5,-10,0);
	RooRealVar p3p2("env_pdf_1_8TeV_powc3_p2","env_pdf_1_8TeV_powc3_p2",-0.1,-10,0);
	RooGenericPdf pow3("env_pdf_1_8TeV_powc3","env_pdf_1_8TeV_powc3","TMath::Power(@0,@1)+@3*TMath::Power(@0,@2)",RooArgList(*x,p3p1,p3p2,p3v1));

	RooRealVar p5v1("env_pdf_1_8TeV_powc5_f1","env_pdf_1_8TeV_powc5_f1",0.5,-100,100);
	RooRealVar p5v2("env_pdf_1_8TeV_powc5_f2","env_pdf_1_8TeV_powc5_f2",0.2,-100,100);
	RooRealVar p5p1("env_pdf_1_8TeV_powc5_p1","env_pdf_1_8TeV_powc5_p1",-0.5,-10,0);
	RooRealVar p5p2("env_pdf_1_8TeV_powc5_p2","env_pdf_1_8TeV_powc5_p2",-0.1,-10,0);
	RooRealVar p5p3("env_pdf_1_8TeV_powc5_p3","env_pdf_1_8TeV_powc5_p3",-0.2,-10,0);
	RooGenericPdf pow5("env_pdf_1_8TeV_powc5","env_pdf_1_8TeV_powc5","TMath::Power(@0,@1)+@4*TMath::Power(@0,@2)+@5*TMath::Power(@0,@3)",RooArgList(*x,p5p1,p5p2,p5p3,p5v1,p5v2));

	// Power laws with coeffs
	RooRealVar e1p1("env_pdf_1_8TeV_expc1_p1","env_pdf_1_8TeV_expc1_p1",-0.2,-2,0);
	RooGenericPdf RooExponential("env_pdf_1_8TeV_expc1","env_pdf_1_8TeV_expc1","TMath::Exp(@0*@1)",RooArgList(*x,e1p1));

	// exp1, exp 3, exp 5
	RooRealVar e3v1("env_pdf_1_8TeV_expc3_f1","env_pdf_1_8TeV_expc3_f1",0.5,-100,100);
	RooRealVar e3p1("env_pdf_1_8TeV_expc3_p1","env_pdf_1_8TeV_expc3_p1",-0.2,-2,0);
	RooRealVar e3p2("env_pdf_1_8TeV_expc3_p2","env_pdf_1_8TeV_expc3_p2",-0.1,-2,0);
	RooGenericPdf exp3("env_pdf_1_8TeV_expc3","env_pdf_1_8TeV_expc3","TMath::Exp(@0*@1)+@3*TMath::Exp(@0*@2)",RooArgList(*x,e3p1,e3p2,e3v1));

	RooRealVar e5v1("env_pdf_1_8TeV_expc5_f1","env_pdf_1_8TeV_expc5_f1",0.5,-100,100);
	RooRealVar e5v2("env_pdf_1_8TeV_expc5_f2","env_pdf_1_8TeV_expc5_f2",0.2,-100,100);
	RooRealVar e5p1("env_pdf_1_8TeV_expc5_p1","env_pdf_1_8TeV_expc5_p1",-0.2,-2,0);
	RooRealVar e5p2("env_pdf_1_8TeV_expc5_p2","env_pdf_1_8TeV_expc5_p2",-0.1,-2,0);
	RooRealVar e5p3("env_pdf_1_8TeV_expc5_p3","env_pdf_1_8TeV_expc5_p3",-0.05,-2,0);
	RooGenericPdf exp5("env_pdf_1_8TeV_expc5","env_pdf_1_8TeV_expc5","TMath::Exp(@0*@1)+@4*TMath::Exp(@0*@2)+@5*TMath::Exp(@0*@3)",RooArgList(*x,e5p1,e5p2,e5p3,e5v1,e5v2));

	// Laurents with coeffs, pain in the ass to write order is -4,-5,-3,-6,-2,-7
	RooRealVar l1p1("env_pdf_1_8TeV_lauc1_p1","env_pdf_1_8TeV_lauc1_p1",0.4,-10,10);
	RooGenericPdf lau1("env_pdf_1_8TeV_lauc1","env_pdf_1_8TeV_lauc1","TMath::Power(@0,-4) + @1*TMath::Power(@0,-5) ",RooArgList(*x,l1p1));
	
	RooRealVar l2p1("env_pdf_1_8TeV_lauc2_p1","env_pdf_1_8TeV_lauc2_p1",0.4,-10,10);
	RooRealVar l2p2("env_pdf_1_8TeV_lauc2_p2","env_pdf_1_8TeV_lauc2_p2",0.2,-10,10);
	RooGenericPdf lau2("env_pdf_1_8TeV_lauc2","env_pdf_1_8TeV_lauc2","TMath::Power(@0,-4) + @1*TMath::Power(@0,-5) + @2*TMath::Power(@0,-3)",RooArgList(*x,l2p1,l2p2));

	RooRealVar l3p1("env_pdf_1_8TeV_lauc3_p1","env_pdf_1_8TeV_lauc3_p1",0.4,-10,10);
	RooRealVar l3p2("env_pdf_1_8TeV_lauc3_p2","env_pdf_1_8TeV_lauc3_p2",0.2,-10,10);
	RooRealVar l3p3("env_pdf_1_8TeV_lauc3_p3","env_pdf_1_8TeV_lauc3_p3",0.1,-10,10);
	RooGenericPdf lau3("env_pdf_1_8TeV_lauc3","env_pdf_1_8TeV_lauc3","TMath::Power(@0,-4) + @1*TMath::Power(@0,-5) + @2*TMath::Power(@0,-3)+ @3*TMath::Power(@0,-6)",RooArgList(*x,l3p1,l3p2,l3p3));

	RooRealVar l4p1("env_pdf_1_8TeV_lauc4_p1","env_pdf_1_8TeV_lauc4_p1",0.4,-10,10);
	RooRealVar l4p2("env_pdf_1_8TeV_lauc4_p2","env_pdf_1_8TeV_lauc4_p2",0.2,-10,10);
	RooRealVar l4p3("env_pdf_1_8TeV_lauc4_p3","env_pdf_1_8TeV_lauc4_p3",0.1,-10,10);
	RooRealVar l4p4("env_pdf_1_8TeV_lauc4_p4","env_pdf_1_8TeV_lauc4_p4",0.05,-10,10);
	RooGenericPdf lau4("env_pdf_1_8TeV_lauc4","env_pdf_1_8TeV_lauc4","TMath::Power(@0,-4) + @1*TMath::Power(@0,-5) + @2*TMath::Power(@0,-3)+ @3*TMath::Power(@0,-6)+ @4*TMath::Power(@0,-2)",RooArgList(*x,l4p1,l4p2,l4p3,l4p4));

	RooRealVar l5p1("env_pdf_1_8TeV_lauc5_p1","env_pdf_1_8TeV_lauc5_p1",0.4,-10,10);
	RooRealVar l5p2("env_pdf_1_8TeV_lauc5_p2","env_pdf_1_8TeV_lauc5_p2",0.2,-10,10);
	RooRealVar l5p3("env_pdf_1_8TeV_lauc5_p3","env_pdf_1_8TeV_lauc5_p3",0.1,-10,10);
	RooRealVar l5p4("env_pdf_1_8TeV_lauc5_p4","env_pdf_1_8TeV_lauc5_p4",0.05,-10,10);
	RooRealVar l5p5("env_pdf_1_8TeV_lauc5_p5","env_pdf_1_8TeV_lauc5_p5",0.02,-10,10);
	RooGenericPdf lau5("env_pdf_1_8TeV_lauc5","env_pdf_1_8TeV_lauc5","TMath::Power(@0,-4) + @1*TMath::Power(@0,-5) + @2*TMath::Power(@0,-3)+ @3*TMath::Power(@0,-6)+ @4*TMath::Power(@0,-2)+ @5*TMath::Power(@0,-7)",RooArgList(*x,l5p1,l5p2,l5p3,l5p4,l5p5));

	wnew->import(pow1);
	wnew->import(pow3);
	wnew->import(pow5);
	wnew->import(exp1);
	wnew->import(exp3);
	wnew->import(exp5);
	wnew->import(lau1);
	wnew->import(lau2);
	wnew->import(lau3);
	wnew->import(lau4);
	wnew->import(lau5);

	fo->cd();
	wnew->Write();
	fo->Close();
	
}
