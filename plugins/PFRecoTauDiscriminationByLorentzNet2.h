#ifndef TallinnTauTag_LorentzNet_PFRecoTauDiscriminationByLorentzNet2_h
#define TallinnTauTag_LorentzNet_PFRecoTauDiscriminationByLorentzNet2_h

#include "FWCore/Framework/interface/stream/EDProducer.h"            // edm::stream::EDProducer<>
#include "FWCore/Framework/interface/Event.h"                        // edm::Event
#include "FWCore/Framework/interface/EventSetup.h"                   // edm::EventSetup
#include "FWCore/ParameterSet/interface/ParameterSet.h"              // edm::ParameterSet
#include "FWCore/Utilities/interface/InputTag.h"                     // edm::InputTag
#include "FWCore/Utilities/interface/EDGetToken.h"                   // edm::EDGetTokenT<>
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h" // edm::ConfigurationDescriptions

#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"          // ONNXRuntime

#include "RecoTauTag/RecoTau/interface/RecoTauQualityCuts.h"         // reco::tau::RecoTauQualityCuts
#include "RecoTauTag/RecoTau/interface/RecoTauVertexAssociator.h"    // reco::tau::RecoTauVertexAssociator 

#include "DataFormats/TauReco/interface/PFTau.h"                     // reco::PFTau
#include "DataFormats/TauReco/interface/PFTauFwd.h"                  // reco::PFTauCollection

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"     // HepPDT::ParticleDataTable, PDTRecord

#include <string>                                                    // std::string

namespace reco 
{
  namespace tau 
  {
    class PFRecoTauDiscriminationByLorentzNet2 : public edm::stream::EDProducer<edm::GlobalCache<cms::Ort::ONNXRuntime>>
    {
     public:
      explicit PFRecoTauDiscriminationByLorentzNet2(const edm::ParameterSet& cfg, const cms::Ort::ONNXRuntime* cache);
      ~PFRecoTauDiscriminationByLorentzNet2() override;

      static
      std::unique_ptr<cms::Ort::ONNXRuntime>
      initializeGlobalCache(const edm::ParameterSet& cfg);

      static
      void 
      globalEndJob(const cms::Ort::ONNXRuntime* cache);

      void 
      produce(edm::Event& evt, const edm::EventSetup& es) override;

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions); 

     private:
      std::string moduleLabel_;

      edm::ESGetToken<HepPDT::ParticleDataTable, PDTRecord> pdt_;

      edm::InputTag srcPFTaus_;
      edm::EDGetTokenT<reco::PFTauCollection> tokenPFTaus_;

      edm::InputTag srcPFCandidates_;
      edm::EDGetTokenT<reco::PFCandidateCollection> tokenPFCandidates_;
      bool addPFCandidatesOutsideJets_;
      double isolationConeSize_;

      RecoTauQualityCuts signalQualityCuts_;
      RecoTauVertexAssociator vertexAssociator_;

      bool applySignalQualityCuts_;
      double maxChargedPFCand_dZ_;
      unsigned maxPFCandsPerJet_;

      bool usePdgMass_;

      int verbosity_;
    };
  }
}

#endif
