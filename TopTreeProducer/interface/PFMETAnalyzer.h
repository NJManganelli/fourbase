#ifndef PFMETAnalyzer_h
#define PFMETAnalyzer_h

// system include files
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "../interface/TRootEvent.h"
#include "../interface/TRootPFMET.h"

#include "../interface/METAnalyzer.h"

#include "TClonesArray.h"


class PFMETAnalyzer{
	
public:
	PFMETAnalyzer(const edm::ParameterSet& myConfig, int verbosity);
	~PFMETAnalyzer();
	void Process(const edm::Event& iEvent, TClonesArray* rootMET, edm::EDGetTokenT<pat::METCollection> metToken);

private:
	int verbosity_;
	bool useMC_;
	METAnalyzer* myMETAnalyzer;
};

#endif

