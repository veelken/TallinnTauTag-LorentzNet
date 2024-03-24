#include "TallinnTauTag/LorentzNet/plugins/PFJetToPFTauConverter.h"

#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"      // edm::ParameterSetDescription
#include "FWCore/MessageLogger/interface/MessageLogger.h"               // edm::LogWarning

#include "DataFormats/Common/interface/RefToBase.h"                     // edm::RefToBase
#include "DataFormats/JetReco/interface/Jet.h"                          // reco::Jet
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"    // reco::PFCandidate, reco::PFCandidatePtr
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h" // reco::PFCandidateCollection
#include "DataFormats/TauReco/interface/PFTau.h"                        // reco::PFTau
#include "DataFormats/TauReco/interface/PFTauFwd.h"                     // reco::PFTauCollection
#include "DataFormats/TrackReco/interface/Track.h"                      // reco::Track
#include "DataFormats/VertexReco/interface/Vertex.h"                    // reco::Vertex::Point

#include <algorithm>                                                    // std::sort
#include <iostream>                                                     // std::cout, std::endl
#include <memory>                                                       // std::unique_ptr
#include <vector>                                                       // std::vector

using namespace reco::tau;

typedef std::vector<std::string> vstring;

PFJetToPFTauConverter::PFJetToPFTauConverter(const edm::ParameterSet& cfg) 
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , srcPFJets_(cfg.getParameter<edm::InputTag>("srcPFJets"))
  , tokenPFJets_(consumes<reco::PFJetCollection>(srcPFJets_))
  , signalQualityCuts_(cfg.getParameterSet("qualityCuts").getParameterSet("signalQualityCuts"))
  , vertexAssociator_(cfg.getParameterSet("qualityCuts"), consumesCollector())
  , verbosity_(cfg.getParameter<int>("verbosity"))
{
  std::cout << "<PFJetToPFTauConverter::PFJetToPFTauConverter (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;

  produces<reco::PFTauCollection>();
}

PFJetToPFTauConverter::~PFJetToPFTauConverter()
{}

//-------------------------------------------------------------------------------
namespace 
{
  bool
  isHigherPt(const reco::PFCandidatePtr& pfCandidate1, const reco::PFCandidatePtr& pfCandidate2)
  {
    return pfCandidate1->pt() > pfCandidate2->pt();
  }

  reco::Candidate::LorentzVector
  getSumP4(const std::vector<reco::PFCandidatePtr>& pfCands)
  {
    reco::Candidate::LorentzVector sumP4;
    for ( const reco::PFCandidatePtr& pfCand : pfCands )
    {
      sumP4 += pfCand->p4();
    }
    return sumP4;
  }

  reco::PFCandidatePtr
  getLeadingPFCandPtr(const std::vector<reco::PFCandidatePtr>& pfCands)
  {
    reco::PFCandidatePtr leadingPFCand;
    double leadingPFCandPt = 0.;
    for ( const reco::PFCandidatePtr& pfCand : pfCands )
    {
      if ( pfCand->pt() > leadingPFCandPt )
      {
        leadingPFCand = pfCand;
        leadingPFCandPt = pfCand->pt();
      }
    }
    return leadingPFCand;
  }

  reco::CandidatePtr
  convertToCandidatePtr(const reco::PFCandidatePtr& pfCand)
  {
    return reco::CandidatePtr(pfCand);
  }

  std::vector<reco::CandidatePtr>
  convertToCandidatePtrs(const std::vector<reco::PFCandidatePtr>& pfCands)
  {
    std::vector<reco::CandidatePtr> pfCandPtrs;
    for ( const reco::PFCandidatePtr& pfCand : pfCands )
    {
      pfCandPtrs.push_back(reco::CandidatePtr(pfCand));
    }
    return pfCandPtrs;
  }

  std::vector<reco::PFCandidatePtr>
  getPFCands_of_type(const std::vector<reco::PFCandidatePtr>& pfCands, const std::vector<int>& particleIds)
  {
    std::vector<reco::PFCandidatePtr> pfCands_of_type;
    for ( const reco::PFCandidatePtr& pfCand : pfCands )
    {
      bool isPFCand_of_type = false;
      for ( auto const& particleId : particleIds )
      {
        if ( pfCand->particleId() == particleId )
        {
          isPFCand_of_type = true;
          break;
        }
      }
      if ( isPFCand_of_type ) pfCands_of_type.push_back(pfCand);
    }
    return pfCands_of_type;
  }
} // namespace
//-------------------------------------------------------------------------------

void 
PFJetToPFTauConverter::produce(edm::Event& evt, const edm::EventSetup& es)
{
  std::unique_ptr<reco::PFTauCollection> pfTaus = std::make_unique<reco::PFTauCollection>();

  edm::Handle<reco::PFJetCollection> pfJets;
  evt.getByToken(tokenPFJets_, pfJets);

  vertexAssociator_.setEvent(evt);

  size_t numPFJets = pfJets->size();
  for ( size_t idxPFJet = 0; idxPFJet < numPFJets; ++idxPFJet )
  {
    reco::PFJetRef pfJet(pfJets, idxPFJet);
    if ( verbosity_ >= 1 )
    {
      std::cout << "pfJet #" << pfJet.key() << ":" 
                << " pT = " << pfJet->pt() << ","
                << " eta = " << pfJet->eta() << ","
                << " phi = " << pfJet->phi() << ","
                << " mass = " << pfJet->mass() << std::endl;
    }

    reco::VertexRef vertex = vertexAssociator_.associatedVertex(*pfJet);
    if ( !vertex.get() )
    {
      if ( verbosity_ >= 1 )
      {
        edm::LogWarning("PFJetToPFTauConverter::produce") 
          << "Failed to find vertex --> skipping !!";
      }
      continue;
    }

    const reco::Track* leadTrack = vertexAssociator_.getLeadTrack(*pfJet);
    if ( !leadTrack ) 
    {
      if ( verbosity_ >= 1 )
      {
        edm::LogWarning("PFJetToPFTauConverter::produce") 
          << "Failed to find leading track --> skipping !!";
      }
      continue;
    }

    signalQualityCuts_.setPV(vertex);
    signalQualityCuts_.setLeadTrack(*leadTrack);

    std::vector<reco::PFCandidatePtr> pfJetConstituents = pfJet->getPFConstituents();
    std::vector<reco::PFCandidatePtr> signalPFCands;
    for ( const reco::PFCandidatePtr& pfJetConstituent : pfJetConstituents )
    {
      bool passesSignalQualityCuts = signalQualityCuts_.filterCand(*pfJetConstituent);
      if ( passesSignalQualityCuts )
      {
        signalPFCands.push_back(pfJetConstituent);
      }
    }
     
    // sort PFCandidates by decreasing pT
    std::sort(signalPFCands.begin(), signalPFCands.end(), isHigherPt);

    reco::PFTau pfTau;
    pfTau.setjetRef(edm::RefToBase<reco::Jet>(pfJet));
    reco::Candidate::LorentzVector signalPFCandP4 = getSumP4(signalPFCands);
    pfTau.setP4(signalPFCandP4);
    pfTau.setDecayMode(reco::PFTau::kRareDecayMode);
    pfTau.setCharge(leadTrack->charge());
    int pdgId = ( leadTrack->charge() >= 0. ) ? -15 : +15; 
    pfTau.setPdgId(pdgId);
    pfTau.setVertex(vertex->position());
    pfTau.setleadCand(convertToCandidatePtr(getLeadingPFCandPtr(signalPFCands)));
    std::vector<reco::PFCandidatePtr> signalPFChargedHadrCands = getPFCands_of_type(signalPFCands, { reco::PFCandidate::h, reco::PFCandidate::e, reco::PFCandidate::mu });
    pfTau.setleadChargedHadrCand(convertToCandidatePtr(getLeadingPFCandPtr(signalPFChargedHadrCands)));
    std::vector<reco::PFCandidatePtr> signalPFGammaCands = getPFCands_of_type(signalPFCands, { reco::PFCandidate::gamma });
    pfTau.setleadNeutralCand(convertToCandidatePtr(getLeadingPFCandPtr(signalPFGammaCands)));
    pfTau.setsignalCands(convertToCandidatePtrs(signalPFCands));
    if ( verbosity_ >= 1 )
    {
      std::cout << "pfTau:" 
                << " pT = " << pfTau.pt() << ","
                << " eta = " << pfTau.eta() << ","
                << " phi = " << pfTau.phi() << ","
                << " mass = " << pfTau.mass() << ","
                << " charge = " << pfTau.charge() << ","
                << " decayMode = " << pfTau.decayMode() << std::endl;
    }
    pfTaus->push_back(pfTau);
  }

  evt.put(std::move(pfTaus));
}

void 
PFJetToPFTauConverter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcPFJets", edm::InputTag("ak4PFJets"));
  edm::ParameterSetDescription desc_qualityCuts;
  RecoTauQualityCuts::fillDescriptions(desc_qualityCuts);
  desc.add<edm::ParameterSetDescription>("qualityCuts", desc_qualityCuts);
  desc.add<int>("verbosity", 0);
  descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"                // DEFINE_FWK_MODULE()
DEFINE_FWK_MODULE(PFJetToPFTauConverter);
