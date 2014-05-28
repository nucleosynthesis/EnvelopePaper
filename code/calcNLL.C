{
   TFile *fi = TFile::Open("envelopews.root");
   RooWorkspace *multipdf = fi->Get("multipdf");
   RooRealVar *x    = multipdf->var("CMS_hgg_mass");

   // The data we will use (a Toy dataset);
   RooDataHist *datatoy = multipdf->data("roohist_data_mass_cat1_toy__CMS_hgg_mass");

   // Build signal model
   RooRealVar mean("mean","mean",125); mean.setConstant();
   RooRealVar sigma("sigma","sigma",1.19); sigma.setConstant();
   RooRealVar nsignal_const("nsignal_const","nsignal_const",50.8); nsignal_const.setConstant();

   RooGaussian sig_pdf("gaus","gaus",*x,mean,sigma);
   
   // Should add this to the RooWorkspace
   RooRealVar nbkg("nbkg","nbkg",datatoy->sumEntries(),0,10E6);

   RooRealVar mu("r","r",1,-10,10);
   RooFormulaVar nsig("nsig","nsig","@0*@1",RooArgList(mu,nsignal_const));
   mu.setConstant();
   
   x->setVal(120);  
   RooAbsPdf *bkg_pdf = multipdf->pdf(Form("env_pdf_1_8TeV_pow%d",1));
   RooAddPdf spdf("splusb","splusb",RooArgList(sig_pdf,*bkg_pdf),RooArgList(nsig,nbkg));


   RooArgSet *params = spdf.getParameters(RooArgSet(0));
   RooRealVar *var = (RooRealVar*)params->find("CMS_hgg_mass");
   var->setVal(180);
   std::cout << " New RooRealVar " << std::endl;
     var->Print("");
   std::cout << " THE RooRealVar " << std::endl;
   x->Print("");

   RooAbsReal *nll = spdf.createNLL(*datatoy,RooFit::Optimize(false),RooFit::Verbose(true),RooFit::Extended(false));
//   nll->Print("V");
    
   std::cout << nll->getVal() << std::endl;   // No difference 
}
