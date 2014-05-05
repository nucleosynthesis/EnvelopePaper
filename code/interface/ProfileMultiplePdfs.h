#ifndef ProfileMultiplePdfs_h 
#define ProfileMultiplePdfs_h

#include <iostream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TGraph.h"

#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"

using namespace std;
using namespace RooFit;

class ProfileMultiplePdfs {

  public:
    ProfileMultiplePdfs(bool verbose=false);
    ~ProfileMultiplePdfs();

    void addPdf(RooAbsPdf *pdf, float penaltyTerm=0.);
    void addPdfs(map<string,RooAbsPdf*> pdfs);
    void clearPdfs();
    void printPdfs();

		// using these for envelope paper:
		void makeProfiles(RooAbsData *data, RooRealVar *poi, float poiLow, float poiHigh, int npoints, string name_ext="", bool printProgress=false);
		map<string,TGraph*> returnProfiles();
		void printProfiles();
		void saveProfiles(TFile *outfile);
		void addProfile(TGraph* graph, int npars);
		void printProfilesWithNPars();
		void setCorrPdof(double corr);
		void constructEnvelope(string ext="");
		void drawEnvelope(string name, bool shiftToZero=false);
		TGraph* shiftGraph(TGraph *graph, double shiftBy);
		TGraph *getEnvelopeGraph();
		double getEnvelopeBestFitValue();
		string getEnvelopeBestFitName();
		double getEnvelopeErrorUp(double sigma);
		double getEnvelopeErrorDn(double sigma);

  private:
    
    map<string,pair<RooAbsPdf*,float> > listOfPdfs;
		map<string,TGraph*> profileGraphs;
		map<string,pair<TGraph*,int> > profileGraphsWithNPars;

		void printProgressBar(float val, float low, float high, bool asBar=true);
		bool verbose_;
		double corr_pdof_;
		TGraph *envelope_;
		double envelopeMinNll_;
		string envelopeBestFitName_;
		double envelopeBestFitVal_;
    RooArgList *pdfsArgSet;
};
  
#endif
