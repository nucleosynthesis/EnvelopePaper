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
		void addProfile(TGraph* graph, float corr_val);
		void addProfile(TGraph* graph, TGraph *correction);
		//void printProfilesWithNPars();
		void printProfilesWithCorrectionGraph();
		void constructEnvelope(string ext="");
		void checkPointsInGraphs();
		void drawEnvelope(string name, string xaxis, bool shiftToZero=false);
		TGraph* shiftGraph(TGraph *graph, double shiftBy);
		TGraph *getEnvelopeGraph();
		double getEnvelopeBestFitValue();
		string getEnvelopeBestFitName();
		double getEnvelopeErrorUp(double sigma);
		double getEnvelopeErrorDn(double sigma);

  private:
    
    map<string,pair<RooAbsPdf*,float> > listOfPdfs;
		map<string,TGraph*> profileGraphs;
		//map<string,pair<TGraph*,float> > profileGraphsWithNPars;
		map<string,pair<TGraph*,TGraph*> > profileGraphsWithCorrectionGraph;

		void printProgressBar(float val, float low, float high, bool asBar=true);
		bool verbose_;
		bool alreadyChecked_;
		TGraph *envelope_;
		double envelopeMinNll_;
		string envelopeBestFitName_;
		double envelopeBestFitVal_;
    RooArgList *pdfsArgSet;

		vector<TGraph*> cleanUp_;
};
  
#endif
