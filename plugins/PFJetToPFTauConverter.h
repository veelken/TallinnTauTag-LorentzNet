#ifndef TallinnTauTag_LorentzNet_PFJetToPFTauConverter_h
#define TallinnTauTag_LorentzNet_PFJetToPFTauConverter_h

#include "FWCore/Framework/interface/stream/EDProducer.h"            // edm::stream::EDProducer<>
#include "FWCore/Framework/interface/Event.h"                        // edm::Event
#include "FWCore/Framework/interface/EventSetup.h"                   // edm::EventSetup
#include "FWCore/ParameterSet/interface/ParameterSet.h"              // edm::ParameterSet
#include "FWCore/Utilities/interface/InputTag.h"                     // edm::InputTag
#include "FWCore/Utilities/interface/EDGetToken.h"                   // edm::EDGetTokenT<>
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h" // edm::ConfigurationDescriptions

#include "RecoTauTag/RecoTau/interface/RecoTauQualityCuts.h"         // reco::tau::RecoTauQualityCuts
#include "RecoTauTag/RecoTau/interface/RecoTauVertexAssociator.h"    // reco::tau::RecoTauVertexAssociator 

#include "DataFormats/JetReco/interface/PFJet.h"                     // reco::PFJet, reco::PFJetRef
#include "DataFormats/JetReco/interface/PFJetCollection.h"           // reco::PFJetCollection

#include <string>                                                    // std::string

namespace reco 
{
  namespace tau 
  {
    class PFJetToPFTauConverter : public edm::stream::EDProducer<>
    {
     public:
      explicit PFJetToPFTauConverter(const edm::ParameterSet& cfg);
      ~PFJetToPFTauConverter() override;

      void 
      produce(edm::Event& evt, const edm::EventSetup& es) override;

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions); 

     private:
      std::string moduleLabel_;

      edm::InputTag srcPFJets_;
      edm::EDGetTokenT<reco::PFJetCollection> tokenPFJets_;

      RecoTauQualityCuts signalQualityCuts_;
      RecoTauVertexAssociator vertexAssociator_;

      int verbosity_;
    };
  }
}

#endif
