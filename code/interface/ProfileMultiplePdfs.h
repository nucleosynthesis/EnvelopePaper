#ifndef ProfileMultiplePdfs_h
#define ProfileMultiplePdfs_h

#include <iostream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"

#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooAbsData.h"
#include "RooWorkspace.h"

using namespace std;
using namespace RooFit;

class ProfileMultiplePdfs {

  public:
    ProfileMultiplePdfs(bool verbose=true);
    ~ProfileMultiplePdfs();

    void addPdf(RooAbsPdf *pdf, float penaltyTerm=0.);
    void addPdfs(map<string,RooAbsPdf*> pdfs);
    void clearPdfs();
    void printPdfs();
    void setObsVar(RooRealVar *var){
     obs_var=var;
     wspace->import(*obs_var);
     obs_var_set=true;
   }

		void setSavePVal(bool val);

		// using these for envelope paper:
		void makeProfiles(RooAbsData *data, RooRealVar *poi, float poiLow, float poiHigh, int npoints, string name_ext="", bool printProgress=false);
		map<string,TGraph*> returnProfiles();
		void printProfiles();
		void saveProfiles(TFile *outfile);
		void addProfile(TGraph* graph, float corr_val);
		void addProfile(TGraph* graph, TGraph *correction);
		double getCorrection(RooAbsData *data, RooAbsPdf *pdf, int add_pars=0);
		double getChisq(RooAbsData *dat, RooAbsPdf *pdf, bool prt=false);

		//void printProfilesWithNPars();
		void printProfilesWithCorrectionGraph();
		void constructEnvelope(string ext="");
		void checkPointsInGraphs();
		void drawEnvelope(string name, string xaxis, bool shiftToZero=false);
		TGraph* shiftGraph(TGraph *graph, double shiftBy);
		TGraph *getEnvelopeGraph();
		double getEnvelopeBestFitValue();
		string getEnvelopeBestFitName();
		double getEnvelopeBestFitNll();
		double getEnvelopeErrorUp(double sigma);
		double getEnvelopeErrorDn(double sigma);

  private:

    map<string,pair<RooAbsPdf*,float> > listOfPdfs;
		map<string,TGraph*> profileGraphs;
		map<string,TGraph*> correctionGraphs;
		//map<string,pair<TGraph*,float> > profileGraphsWithNPars;
		map<string,pair<TGraph*,TGraph*> > profileGraphsWithCorrectionGraph;

		void printProgressBar(float val, float low, float high, bool asBar=true);
		bool verbose_;
		bool alreadyChecked_;
		bool savePVal_;
		TGraph *envelope_;
		double envelopeMinNll_;
		string envelopeBestFitName_;
		double envelopeBestFitVal_;
    RooArgList *pdfsArgSet;
    		RooRealVar *obs_var;
		bool obs_var_set;
		double globMu;
		RooWorkspace *wspace;
		vector<TGraph*> cleanUp_;
};

#endif
