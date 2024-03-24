#include "TallinnTauTag/LorentzNet/plugins/PFRecoTauDiscriminationByLorentzNet.h"

#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"          // edm::ParameterSetDescription
#include "FWCore/MessageLogger/interface/MessageLogger.h"                   // edm::LogWarning

#include "DataFormats/Common/interface/RefToPtr.h"                          // edm::refToPtr
#include "DataFormats/JetReco/interface/PFJet.h"                            // reco::PFJet, reco::PFJetRef
#include "DataFormats/JetReco/interface/PFJetCollection.h"                  // reco::PFJetCollection
#include "DataFormats/Math/interface/deltaR.h"                              // reco::deltaR()
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"        // reco::PFCandidate, reco::PFCandidatePtr
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"     // reco::PFCandidateCollection
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"               // reco::PFTauDiscriminator
#include "DataFormats/TrackReco/interface/Track.h"                          // reco::Track
#include "DataFormats/VertexReco/interface/Vertex.h"                        // reco::Vertex::Point

#include "TallinnTauTag/Tools/interface/format_vT.h"                        // format_vfloat

#include <algorithm>                                                        // std::sort
#include <cmath>                                                            // std::fabs, std::log, std::sqrt
#include <iostream>                                                         // std::cout, std::endl
#include <memory>                                                           // std::unique_ptr
#include <vector>                                                           // std::vector

using namespace reco::tau;

typedef std::vector<std::string> vstring;

PFRecoTauDiscriminationByLorentzNet::PFRecoTauDiscriminationByLorentzNet(const edm::ParameterSet& cfg, const cms::Ort::ONNXRuntime* cache) 
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , srcPFTaus_(cfg.getParameter<edm::InputTag>("srcPFTaus"))
  , tokenPFTaus_(consumes<reco::PFTauCollection>(srcPFTaus_))
  , srcPFCandidates_(cfg.getParameter<edm::InputTag>("srcPFCandidates"))
  , tokenPFCandidates_(consumes<reco::PFCandidateCollection>(srcPFCandidates_))
  , addPFCandidatesOutsideJets_(cfg.getParameter<bool>("addPFCandidatesOutsideJets"))
  , isolationConeSize_(cfg.getParameter<double>("isolationConeSize"))
  , signalQualityCuts_(cfg.getParameterSet("qualityCuts").getParameterSet("signalQualityCuts"))
  , vertexAssociator_(cfg.getParameterSet("qualityCuts"), consumesCollector())
  , applySignalQualityCuts_(cfg.getParameter<bool>("applySignalQualityCuts"))
  , maxChargedPFCand_dZ_(cfg.getParameter<double>("maxChargedPFCand_dZ"))
  , maxPFCandsPerJet_(cfg.getParameter<unsigned>("maxPFCandsPerJet"))
  , verbosity_(cfg.getParameter<int>("verbosity"))
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<PFRecoTauDiscriminationByLorentzNet::PFRecoTauDiscriminationByLorentzNet (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  }

  produces<reco::PFTauDiscriminator>();
}

PFRecoTauDiscriminationByLorentzNet::~PFRecoTauDiscriminationByLorentzNet()
{}

std::unique_ptr<cms::Ort::ONNXRuntime>
PFRecoTauDiscriminationByLorentzNet::initializeGlobalCache(const edm::ParameterSet& cfg) 
{
  return std::make_unique<cms::Ort::ONNXRuntime>(cfg.getParameter<edm::FileInPath>("model").fullPath());
}

void 
PFRecoTauDiscriminationByLorentzNet::globalEndJob(const cms::Ort::ONNXRuntime* cache) 
{}

//-------------------------------------------------------------------------------
namespace 
{
  // CV: function helper::init_result_object copied from
  //       RecoTauTag/RecoTau/src/TauDiscriminationProducerBase.cc
  struct helper 
  {
    static std::unique_ptr<reco::PFTauDiscriminator> 
    init_result_object(const edm::Handle<reco::PFTauCollection>& taus) 
    {
      return std::make_unique<reco::PFTauDiscriminator>(edm::RefProd<reco::PFTauCollection>(taus));
    }
  };

  bool
  isHigherPt(const reco::PFCandidatePtr& pfCandidate1, const reco::PFCandidatePtr& pfCandidate2)
  {
    return pfCandidate1->pt() > pfCandidate2->pt();
  }

  bool 
  isSameParticle(const reco::PFCandidate& pfCandidate1, const reco::PFCandidate& pfCandidate2)
  {
    if ( pfCandidate1.particleId()  != pfCandidate2.particleId() ) return false;
    const double epsilon = 1.e-3;
    if ( std::fabs(pfCandidate1.energy() - pfCandidate2.energy()) > epsilon*0.5*(pfCandidate1.energy() + pfCandidate2.energy()) ) return false;
    double dR = deltaR(pfCandidate1.p4(), pfCandidate2.p4());
    if ( dR < epsilon ) return true;
    else                return false;
  }

  int
  sgn(float val) 
  {
    return (0 < val) - (val < 0);
  }

  float
  psi(float val)
  {
    return sgn(val)*std::log(std::fabs(val) + 1.);
  }

  void
  extend(std::vector<float>& values, const std::vector<float>& to_add)
  {
    for ( float value : to_add )
    {
      values.push_back(value);
    }
  }
} // namespace
//-------------------------------------------------------------------------------

void 
PFRecoTauDiscriminationByLorentzNet::produce(edm::Event& evt, const edm::EventSetup& es)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<PFRecoTauDiscriminationByLorentzNet::produce (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
  }

  edm::Handle<reco::PFTauCollection> pfTaus;
  evt.getByToken(tokenPFTaus_, pfTaus);

  std::unique_ptr<reco::PFTauDiscriminator> pfTauDiscriminator = helper::init_result_object(pfTaus);

  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  evt.getByToken(tokenPFCandidates_, pfCandidates);

  vertexAssociator_.setEvent(evt);

  size_t numPFTaus = pfTaus->size();
  for ( size_t idxPFTau = 0; idxPFTau < numPFTaus; ++idxPFTau )
  {
    reco::PFTauRef pfTau(pfTaus, idxPFTau);

    const reco::PFJet* pfJet = dynamic_cast<const reco::PFJet*>(pfTau->jetRef().get());
    if ( verbosity_ >= 1 )
    {
      std::cout << "pfJet #" << pfTau->jetRef().key() << ":" 
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
        edm::LogWarning("PFRecoTauDiscriminationByLorentzNet::analyze") 
          << "Failed to find vertex --> skipping !!";
      }
      continue;
    }

    const reco::Track* leadTrack = vertexAssociator_.getLeadTrack(*pfJet);
    if ( !leadTrack ) 
    {
      if ( verbosity_ >= 1 )
      {
        edm::LogWarning("PFRecoTauDiscriminationByLorentzNet::analyze") 
          << "Failed to find leading track --> skipping !!";
      }
      continue;
    }

    signalQualityCuts_.setPV(vertex);
    signalQualityCuts_.setLeadTrack(*leadTrack);

    std::vector<reco::PFCandidatePtr> pfJetConstituents = pfJet->getPFConstituents();
    std::vector<reco::PFCandidatePtr> associatedPFCandidates;
    for ( const reco::PFCandidatePtr& pfJetConstituent : pfJetConstituents )
    {
      associatedPFCandidates.push_back(pfJetConstituent);
    }

    if ( addPFCandidatesOutsideJets_ )
    {
      size_t numPFCandidates = pfCandidates->size();
      for ( size_t idxPFCandidate = 0; idxPFCandidate < numPFCandidates; ++idxPFCandidate )
      {
        reco::PFCandidateRef pfCandidate(pfCandidates, idxPFCandidate);
        double dR = reco::deltaR(pfCandidate->p4(), pfJet->p4());
        if ( dR < isolationConeSize_ ) 
        {
          bool isPFJetConstituent = false;
          for ( const reco::PFCandidatePtr& pfJetConstituent : pfJetConstituents )
          {
            if ( isSameParticle(*pfCandidate, *pfJetConstituent) )
            {
              isPFJetConstituent = true;
              break;
            }
          }
          if ( !isPFJetConstituent )
          {
            associatedPFCandidates.push_back(edm::refToPtr(pfCandidate));
          }
        } 
      }
    }
     
    // sort PFCandidates by decreasing pT
    std::sort(associatedPFCandidates.begin(), associatedPFCandidates.end(), isHigherPt);

    std::vector<reco::PFCandidatePtr> selectedPFCandidates;
    unsigned numSelectedPFCandidates = 0;
    for ( const reco::PFCandidatePtr& pfCandidate : associatedPFCandidates )
    {
      if ( applySignalQualityCuts_ )
      {
        bool passesSignalQualityCuts = signalQualityCuts_.filterCand(*pfCandidate);
        if ( passesSignalQualityCuts ) continue;
      }
      else
      {
        if ( std::fabs(pfCandidate->charge()) > 0.5 )
        {
          const reco::Track* track = pfCandidate->bestTrack();
          double dZ = track->dz(vertex->position());
          if ( std::fabs(dZ) > maxChargedPFCand_dZ_ ) continue;
        }
      }
      if ( numSelectedPFCandidates < maxPFCandsPerJet_ )
      {
        selectedPFCandidates.push_back(pfCandidate);
        ++numSelectedPFCandidates;
      }
    }
    if ( verbosity_ >= 1 )
    {
      std::cout << "numSelectedPFCandidates = " << numSelectedPFCandidates << std::endl;
    }

    std::vector<float> outputs = { -1., -1. };
    if ( numSelectedPFCandidates >= 1 )
    {
      vstring input_names = { "scalars", "x", "edges0", "edges1", "node_mask" };
      std::vector<float> input_scalars;
      std::vector<float> input_x;
      std::vector<float> input_node_mask;
      for ( int idxPFCandidate = 0; idxPFCandidate < int(numSelectedPFCandidates + 2); ++idxPFCandidate )
      {
        if ( idxPFCandidate == 0 )
        {
          extend(input_scalars, { 0., psi(1.) });
          extend(input_x, { std::sqrt(2.), 0., 0., +1. });
          input_node_mask.push_back(true);
        }
        else if ( idxPFCandidate == 1 )
        {
          extend(input_scalars, { 0., psi(1.) });
          extend(input_x, { std::sqrt(2.), 0., 0., -1. });
          input_node_mask.push_back(true);
        }
        else
        {
          const reco::PFCandidatePtr& pfCandidate = selectedPFCandidates[idxPFCandidate - 2];
          extend(input_scalars, { psi(pfCandidate->p4().mass()), 0. });
          extend(input_x, { float(pfCandidate->energy()), float(pfCandidate->px()), float(pfCandidate->py()), float(pfCandidate->pz()) });
          input_node_mask.push_back(true);
        }
      }
      if ( verbosity_ >= 1 )
      {
        std::cout << "shape(scalars) = " << input_scalars.size() << std::endl;
        std::cout << "scalars = " << format_vfloat(input_scalars) << std::endl;
        std::cout << "shape(x) = " << input_x.size() << std::endl;
        std::cout << "x = " << format_vfloat(input_x) << std::endl;
        std::cout << "shape(node_mask) = " << input_node_mask.size() << std::endl;
        std::cout << "node_mask = " << format_vfloat(input_node_mask) << std::endl;
      }
      std::vector<float> input_edges0;
      std::vector<float> input_edges1;
      for ( int idxPFCandidate1 = 0; idxPFCandidate1 < int(numSelectedPFCandidates + 2); ++idxPFCandidate1 )
      {
        for ( int idxPFCandidate2 = 0; idxPFCandidate2 < int(numSelectedPFCandidates + 2); ++idxPFCandidate2 )
        {
          if ( idxPFCandidate2 == idxPFCandidate1 ) continue;
          input_edges0.push_back(idxPFCandidate1);
          input_edges1.push_back(idxPFCandidate2);
        }
      }
      assert(input_edges0.size() == input_edges1.size());
      if ( verbosity_ >= 1 )
      {
        std::cout << "shape(edges0) = " << input_edges0.size() << std::endl;
        std::cout << "edges0 = " << format_vfloat(input_edges0) << std::endl;
        std::cout << "shape(edges1) = " << input_edges1.size() << std::endl;
        std::cout << "edges1 = " << format_vfloat(input_edges1) << std::endl;
      }
      int64_t num_particles = input_node_mask.size();
      int64_t num_edges = input_edges0.size();
      cms::Ort::FloatArrays input_values = { input_scalars, input_x, input_edges0, input_edges1, input_node_mask };
      const int64_t batch_size = 1;
      std::vector<std::vector<int64_t>> input_shapes = {{ num_particles, 2 }, { num_particles, 4 }, { num_edges }, { num_edges }, { num_particles }};
      std::vector<std::string> output_names = { "output" };
      outputs = globalCache()->run(input_names, input_values, input_shapes, output_names, batch_size)[0];
    }
    if ( verbosity_ >= 1 )
    {
      std::cout << "pfTau:" 
                << " pT = " << pfTau->pt() << ","
                << " eta = " << pfTau->eta() << ","
                << " phi = " << pfTau->phi() << ","
                << " mass = " << pfTau->mass() << ","
                << " charge = " << pfTau->charge() << ","
                << " decayMode = " << pfTau->decayMode() << ","
                << " outputs = {" << outputs[0] << ", " << outputs[1] << "}" << std::endl;
    }

    (*pfTauDiscriminator)[pfTau] = outputs[0];
  }

  evt.put(std::move(pfTauDiscriminator));
}

void 
PFRecoTauDiscriminationByLorentzNet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcPFTaus", edm::InputTag("hpsPFTauProducer"));
  desc.add<edm::InputTag>("srcPFCandidates", edm::InputTag("particleFlow"));
  desc.add<bool>("addPFCandidatesOutsideJets", true);
  desc.add<double>("isolationConeSize", 0.5);
  edm::ParameterSetDescription desc_qualityCuts;
  RecoTauQualityCuts::fillDescriptions(desc_qualityCuts);
  desc.add<edm::ParameterSetDescription>("qualityCuts", desc_qualityCuts);
  desc.add<bool>("applySignalQualityCuts", false);
  desc.add<double>("maxChargedPFCand_dZ", 0.2);
  desc.add<unsigned>("maxPFCandsPerJet", 100);
  desc.add<edm::FileInPath>("model", edm::FileInPath("TallinnTauTag/LorentzNet/data/LorentzNet_tau_2022Dec11.onnx"));
  desc.add<int>("verbosity", 0);
  descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"                // DEFINE_FWK_MODULE()
DEFINE_FWK_MODULE(PFRecoTauDiscriminationByLorentzNet);
