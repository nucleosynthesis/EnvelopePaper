// Simple Macro to produce 2NLL ratio value which is similar to chi-squared

double getChisq(RooAbsData &dat, RooAbsPdf &pdf, RooRealVar &var, bool prt=false) {

  // Find total number of events
  double nEvt;
  double nTot=0.0;

  for(int j=0;j<dat.numEntries();j++) {
    dat.get(j);
    nEvt=dat.weight();
    nTot+=nEvt;    
  }

  // Find chi-squared equivalent 2NLL
  //RooRealVar *var=(RooRealVar*)(pdf.getParameters(*dat)->find("CMS_hgg_mass"));
  double totNLL=0.0;
  double prbSum=0.0;

  for(int j=0;j<dat.numEntries();j++) {
    double m=dat.get(j)->getRealValue(var.GetName());
    if ( m < var.getMin() || m > var.getMax())  continue;
    // Find probability density and hence probability
    var.setVal(m);
    double prb=0.25*pdf.getVal(var);
    prbSum+=prb;

    dat.get(j);
    nEvt=dat.weight();
	  
    double mubin=nTot*prb;
    double contrib(0.);
    if (nEvt < 1) contrib = mubin;
    else contrib=mubin-nEvt+nEvt*log(nEvt/mubin);
    totNLL+=contrib;
    
    if(prt) cout << "Bin " << j << " prob = " << prb << " nEvt = " << nEvt << ", mu = " << mubin << " contribution " << contrib << endl;    
  }
  
  totNLL*=2.0;
  if(prt) cout << pdf.GetName() << " nTot = " << nTot << " 2NLL constant = " << totNLL << endl;
  
  return totNLL;
}
