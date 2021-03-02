#ifndef DELPHES2LCIO_DELPHESLCIOCONVERTER_H
#define DELPHES2LCIO_DELPHESLCIOCONVERTER_H

// LCIO
#include "EVENT/MCParticle.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "UTIL/LCRelationNavigator.h"
#include "lcio.h"

// ROOT
#include "TTree.h"
#include "TClonesArray.h"

//Delphes
#include "modules/Delphes.h"

#include <array>
#include <string_view>
#include <string>
#include <vector>

// Delphes output classes
class Muon;
class Electron;
class Photon;

namespace delphes_lcio {

/**
 * Classes that will be stored as reconstructed particle with an attached track
 */
constexpr std::array<std::string_view, 1> RECO_TRACK_OUTPUT = {"Track"};

/**
 * Classes that will be stored as reconstructed particle with an attached cluster
 */
constexpr std::array<std::string_view, 1> RECO_CLUSTER_OUTPUT = {"Tower"};

struct BranchSettings {
  std::string input;
  std::string name;
  std::string className;
};

std::vector<BranchSettings> getBranchSettings(ExRootConfParam /*const&*/treeConf) {
  std::vector<BranchSettings> branches;
  for (int b = 0; b < treeConf.GetSize(); b += 3) {
    BranchSettings branch{treeConf[b].GetString(),
                                        treeConf[b + 1].GetString(),
                                        treeConf[b + 2].GetString()};
    branches.push_back(branch);
  }
  return branches;
}

class OutputSettings;

class DelphesLCIOConverter {
public:
  DelphesLCIOConverter(const std::vector<BranchSettings>& branches,
                       OutputSettings const& outputSettings, double magFieldBz);

  void process(TTree* delphesTree, lcio::LCEventImpl* evt);

  void processEvtHeader(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt, std::string const&);

  void processParticles(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt, std::string const& branch);
  void processTracks(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt, std::string const& branch);
  void processClusters(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt,  std::string const& branch);
  void processJets(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt, std::string const& branch);
  void processPhotons(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt, std::string const& branch) {
    fillReferenceCollection<Photon>(delphesCollection, evt, branch, "photon");
  }

  void processMissingET(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt, std::string const& branch);
  void processScalarHT(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt, std::string const& branch);

  void processMuons(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt, std::string const& branch) {
    fillReferenceCollection<Muon>(delphesCollection, evt, branch, "muon");
  }
  void processElectrons(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt, std::string const& branch) {
    fillReferenceCollection<Electron>(delphesCollection, evt, branch, "electron");
  }

private:

  template<typename DelphesT>
  void fillReferenceCollection(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt, std::string const& branch,
                               const std::string_view type);

  // cannot mark DelphesT as const, because for Candidate* the GetCandidates()
  // method is not marked as const.  template<typename DelphesT>
  template<typename DelphesT>
  lcio::ReconstructedParticleImpl* getMatchingReco(/*const*/ DelphesT* delphesCand) const;

  using ProcessFunction = void (DelphesLCIOConverter::*)(const TClonesArray*, lcio::LCEventImpl*, const std::string&);

  std::vector<BranchSettings> m_branches;
  // std::unordered_map<std::string_view, podio::CollectionBase*> m_collections;
  std::unordered_map<std::string_view, ProcessFunction> m_processFunctions;

  double m_magneticFieldBz; // necessary for determining track parameters

  std::string m_recoCollName;
  std::string m_particleIDName;
  std::string m_mcRecoAssocCollName;

  // reco particle collection, populated by more than one processing function
  lcio::LCCollectionVec* m_recoParticles{nullptr};

  // delphes2lcio populates both, but in this case we only do one to be more
  // comparable to what k4SimDelphes does
  lcio::LCRelationNavigator* m_recoToMCNav{nullptr};
  // lcio::LCRelationNavigator* m_mcToRecoNav{nullptr};

  // map from UniqueIDs (delphes generated particles) to MCParticles
  std::unordered_map<UInt_t, lcio::MCParticle*> m_genParticleIds;
  // map from UniqueIDs (delphes generated particles) to (possibly multiple)
  // ReconstructedParticles
  std::unordered_multimap<UInt_t, lcio::ReconstructedParticleImpl*> m_recoParticleGenIds;
};

} // namespace

#endif
