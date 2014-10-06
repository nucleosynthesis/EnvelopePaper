#ifndef PdfModelBuilder_h
#define PdfModelBuilder_h

#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooProduct.h"
#include "RooWorkspace.h"
#include "TFile.h"

using namespace std;
using namespace RooFit;

class PdfModelBuilder {

  public:
    PdfModelBuilder();
    ~PdfModelBuilder();

    void setObsVar(RooRealVar *var);
    void setSignalModifier(RooRealVar *var);
    void setSignalModifierVal(float val);
    void setSignalModifierConstant(bool val);

    void addBkgPdf(string name, bool cache=true);
    void addBkgPdf(RooAbsPdf *pdf, bool cache=true);
    RooAbsPdf* getPdfFromFile(string &prefix);

		void setSignalPdf(RooAbsPdf *pdf, RooRealVar *norm=NULL);
    void setSignalPdfFromMC(RooDataSet *data);
    void makeSBPdfs(bool cache=true);

    map<string,RooAbsPdf*> getBkgPdfs();
    map<string,RooAbsPdf*> getSBPdfs();
    RooAbsPdf *getSigPdf();

    void fitToData(RooAbsData *data, bool bkgOnly=true, bool cache=true, bool print=false);
    void plotPdfsToData(RooAbsData *data, int binning, string name, bool bkgOnly=true, string specificPdfName="");
    void plotToysWithPdfs(string prefix, int binning, bool bkgOnly=true);

    void setSeed(int seed);
    void throwToy(string name, int nEvents, bool bkgOnly=true, bool binned=true, bool poisson=true, bool cache=true);
    map<string,RooAbsData*> getToyData();
		RooAbsData *getToyDataSingle();

    void saveWorkspace(TFile* file);
    void saveWorkspace(string filename);

    RooAbsPdf* returnProfiledBackground(string name);
    double getFractionValue(string name);

    void setPValueCorrection(bool arg=true);
    void setCorrectionValue(double value);

  private:

    map<string,RooAbsPdf*> bkgPdfs;
    map<string,RooAbsPdf*> sbPdfs;
    RooAbsPdf* sigPdf;
    RooAbsReal* sigNorm;
    RooRealVar *bkgYield;
    RooAbsReal *sigYield;

    // this is for the frequentist toy function
    double getChisq(RooAbsData *dat, RooAbsPdf *pdf, bool prt=false);
    map<string,double> cachedNllVals;
    double correctionValue;
    bool isPValueCorrection;

    map<string,RooAbsData*> toyData;
    map<string,RooDataSet*> toyDataSet;
    map<string,RooDataHist*> toyDataHist;

    map<string,RooRealVar*> params;
    map<string,RooFormulaVar*> prods;
    map<string,RooAbsPdf*> utilities;

    RooRealVar *obs_var;
    bool obs_var_set;
    RooRealVar *signalModifier;
    bool signal_modifier_set;
    bool signal_set;
    bool bkgHasFit;
    bool sbHasFit;
    vector<string> cut_strings;

    RooWorkspace *wsCache;

    int verbosity;

};
#endif
