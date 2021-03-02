#include "DelphesLCIOConverter_new.h"
#include "DelphesLCIOOutputConfiguration.h"
#include "delphesHelpers.h"

#include "IMPL/LCCollectionVec.h"
#include "IMPL/MCParticleImpl.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/TrackImpl.h"
#include "IMPL/TrackStateImpl.h"
#include "IMPL/ClusterImpl.h"
#include "IMPL/ParticleIDImpl.h"
#include "UTIL/PIDHandler.h"

#include "classes/DelphesClasses.h"

namespace delphes_lcio {

// TODO: Take these from HepPDT / HepMC?
constexpr double M_PIPLUS = 0.13957039; // GeV (PDG 2020)
constexpr double M_MU = 0.1056583745; // GeV (PDG 2020)
constexpr double M_ELECTRON = 0.5109989461e-3; // GeV (PDG 2020)

// TODO: Make configurable?
constexpr double trackMass = M_PIPLUS;

/**
 * Order in which the different delphes output classes will be processed.
 * Everything not defined here will not be processed.
 *
 * NOTE: not a configuration parameter. this has to be done in this order to
 * ensure that products required by later stages are producd early enough
 */
constexpr std::array<std::string_view, 10> PROCESSING_ORDER = {
  "Event",
  "GenParticle",
  "Track",
  "Tower",
  "Muon",
  "Electron",
  "Photon",
  "Jet",
  "MissingET",
  "SclalarHT"
};

template<size_t N>
void sortBranchesProcessingOrder(std::vector<BranchSettings>& branches,
                                 std::array<std::string_view, N> const& processingOrder);

void setMotherDaughterRelations(GenParticle const* delphesCand,
                                lcio::MCParticleImpl* particle,
                                lcio::LCCollectionVec* mcParticles);

lcio::TrackImpl* convertTrack(Track const* cand, const double magFieldBz);

/**
 * Simple helper function to make it easier to refactor later
 */
template<typename Container>
inline bool contains(Container const& container, typename Container::value_type const& value)
{
  return std::find(container.cbegin(), container.cend(), value) != container.cend();
}



DelphesLCIOConverter::DelphesLCIOConverter(const std::vector<BranchSettings>& branches,
                                           OutputSettings const& outputSettings, double magFieldBz) :
  m_magneticFieldBz(magFieldBz),
  m_recoCollName(outputSettings.RecoParticleCollectionName),
  m_particleIDName(outputSettings.ParticleIDCollectionName),
  m_mcRecoAssocCollName(outputSettings.MCRecoAssociationCollectionName)
{
  for (const auto& branch : branches) {
    if (contains(PROCESSING_ORDER, branch.className)) {
      m_branches.push_back(branch);
    }
  }
  m_branches.emplace_back(BranchSettings{"Event", "Event", "LHEFEvent"});
  sortBranchesProcessingOrder(m_branches, PROCESSING_ORDER);

  const std::unordered_map<std::string, ProcessFunction> refProcessFunctions = {
    {"Photon", &DelphesLCIOConverter::processPhotons},
    {"Muon", &DelphesLCIOConverter::processMuons},
    {"Electron", &DelphesLCIOConverter::processElectrons}};

  m_processFunctions.emplace("Event", &DelphesLCIOConverter::processEvtHeader);

  for (const auto& branch : m_branches) {
    if (contains(outputSettings.GenParticleCollections, branch.name.c_str())) {
      m_processFunctions.emplace(branch.name, &DelphesLCIOConverter::processParticles);
    }

    if (contains(outputSettings.ReconstructedParticleCollections, branch.name.c_str()) &&
        contains(RECO_TRACK_OUTPUT, branch.className.c_str())) {
      m_processFunctions.emplace(branch.name, &DelphesLCIOConverter::processTracks);
    }

    if (contains(outputSettings.ReconstructedParticleCollections, branch.name.c_str()) &&
        contains(RECO_CLUSTER_OUTPUT, branch.className.c_str())) {
      m_processFunctions.emplace(branch.name, &DelphesLCIOConverter::processClusters);
    }

    if (contains(outputSettings.JetCollections, branch.name.c_str())) {
      m_processFunctions.emplace(branch.name, &DelphesLCIOConverter::processJets);
    }

    if (contains(outputSettings.MuonCollections, branch.name.c_str()) ||
        contains(outputSettings.ElectronCollections, branch.name.c_str()) ||
        contains(outputSettings.PhotonCollections, branch.name.c_str())) {
      m_processFunctions.emplace(branch.name, refProcessFunctions.at(branch.className));
    }

    if (contains(outputSettings.MissingETCollections, branch.name.c_str())) {
      m_processFunctions.emplace(branch.name, &DelphesLCIOConverter::processMissingET);
    }

    if (contains(outputSettings.ScalarHTCollections, branch.name.c_str())) {
      m_processFunctions.emplace(branch.name, &DelphesLCIOConverter::processScalarHT);
    }
  }
}


void DelphesLCIOConverter::process(TTree* delphesTree, lcio::LCEventImpl* evt) {
  m_recoParticles = new lcio::LCCollectionVec(lcio::LCIO::RECONSTRUCTEDPARTICLE);
  evt->addCollection(m_recoParticles, m_recoCollName);

  m_recoToMCNav = new lcio::LCRelationNavigator(lcio::LCIO::RECONSTRUCTEDPARTICLE, lcio::LCIO::MCPARTICLE);
  // m_mcToRecoNav = new lcio::LCRelationNavigator(lcio::LCIO::MCPARTICLE, lcio::LCIO::RECONSTRUCTEDPARTICLE );

  for (const auto& branch : m_branches) {
    // at this point it is not guaranteed that all entries in branch (which follow
    // the input from the delphes card) are also present in the processing
    // functions. Checking this here, basically allows us to skip these checks
    // in all called processing functions, since whatever is accessed there will
    // also be in the collection map, since that is filled with the same keys as
    // the processing function map
    const auto processFuncIt = m_processFunctions.find(branch.name);
    // Also check whether the desired input is actually in the branch
    auto* rootBranch = delphesTree->GetBranch(branch.name.c_str());
    if (processFuncIt != m_processFunctions.end() && rootBranch) {
      auto* delphesCollection = *(TClonesArray**) rootBranch->GetAddress();
      (this->*processFuncIt->second)(delphesCollection, evt, branch.name);
    }
  }

  // Clear the internal maps that hold references to entites that have been put
  // into maps here for internal use only (see #89)
  m_genParticleIds.clear();
  m_recoParticleGenIds.clear();

  evt->addCollection(m_recoToMCNav->createLCCollection(), "RecoMCTruthLink");
  // evt->addCollection(m_mcToRecoNav->createLCCollection(), "MCTruthRecoLink");

  m_recoParticles = nullptr;

  delete m_recoToMCNav;
  // delete m_mcToRecoNav;
}

void DelphesLCIOConverter::processEvtHeader(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt, std::string const&) {
  auto* delphesEvt = static_cast<LHEFEvent*>(delphesCollection->At(0));
  evt->setRunNumber(0);
  evt->setEventNumber(delphesEvt->Number);
  evt->setWeight(delphesEvt->Weight);
}


void DelphesLCIOConverter::processParticles(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt, std::string const& branch)
{
  auto* collection = new lcio::LCCollectionVec(lcio::LCIO::MCPARTICLE);
  evt->addCollection(collection, branch);

  for (int iCand = 0; iCand < delphesCollection->GetEntries(); ++iCand) {
    auto* delphesCand = static_cast<GenParticle*>(delphesCollection->At(iCand));
    auto* cand = new lcio::MCParticleImpl;
    collection->addElement(cand);

    cand->setCharge(delphesCand->Charge);
    cand->setPDG(delphesCand->PID);
    cand->setMass(delphesCand->Mass);
    cand->setGeneratorStatus(delphesCand->Status);
    const double mom[3] = {delphesCand->Px, delphesCand->Py, delphesCand->Pz};
    cand->setMomentum(mom);
    const double vertex[3] = {delphesCand->X, delphesCand->Y, delphesCand->Z};
    cand->setVertex(vertex);

    if (const auto [it, inserted] = m_genParticleIds.emplace(delphesCand->GetUniqueID(), cand); !inserted) {
      std::cerr << "**** WARNING: UniqueID " << delphesCand->GetUniqueID() << " is already used by MCParticle with id: " << it->second << std::endl;
    }
  }

  // mother-daughter relations
  for (int iCand = 0; iCand < delphesCollection->GetEntries(); ++iCand) {
    const auto* delphesCand = static_cast<GenParticle*>(delphesCollection->At(iCand));
    auto cand = static_cast<lcio::MCParticleImpl*>(collection->at(iCand));

    setMotherDaughterRelations(delphesCand, cand, collection);
  }
}

void DelphesLCIOConverter::processTracks(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt, std::string const& branch) {
  auto* trackCollection = new lcio::LCCollectionVec(lcio::LCIO::TRACK);
  evt->addCollection(trackCollection, branch);
 
  for (auto iCand = 0; iCand < delphesCollection->GetEntries(); ++iCand) {
    auto* delphesCand = static_cast<Track*>(delphesCollection->At(iCand));
    auto* track = convertTrack(delphesCand, m_magneticFieldBz);
    trackCollection->addElement(track);

    auto* cand = new lcio::ReconstructedParticleImpl;
    cand->setCharge(delphesCand->Charge);
    const auto momentum = delphesCand->P4();
    cand->setEnergy(momentum.E());
    float mom[3] = {(float) momentum.Px(), (float) momentum.Py(), (float) momentum.Pz()};
    cand->setMomentum(mom);
    cand->setMass(trackMass);

    cand->addTrack(track);

    const auto genId = delphesCand->Particle.GetUniqueID();
    if (const auto genIt = m_genParticleIds.find(genId); genIt != m_genParticleIds.end()) {
      // m_mcToRecoNav->addRelation(genIt->second, cand, 1.0);
      m_recoToMCNav->addRelation(cand, genIt->second, 1.0);
    }

    m_recoParticleGenIds.emplace(genId, cand);
  }

}

void DelphesLCIOConverter::processClusters(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt, std::string const& branch) {
  auto* clusterCollection = new lcio::LCCollectionVec(lcio::LCIO::CLUSTER);
  evt->addCollection(clusterCollection, branch);

  for (auto iCand = 0; iCand < delphesCollection->GetEntries(); ++iCand) {
    auto* delphesCand = static_cast<Tower*>(delphesCollection->At(iCand));
    auto* cluster = new lcio::ClusterImpl;
    clusterCollection->addElement(cluster);

    cluster->setEnergy(delphesCand->E);
    // TODO: more info. But not set in k4SimDelhpes either yet

    auto* cand = new lcio::ReconstructedParticleImpl;
    cand->setEnergy(delphesCand->E);
    const auto momentum = delphesCand->P4(); // NOTE: assuming massless here!
    const float mom[3] = {(float) momentum.Px(), (float) momentum.Py(), (float) momentum.Pz()};
    cand->setMomentum(mom);

    cand->addCluster(cluster);

    for (const auto genId : getAllParticleIDs(delphesCand)) {
      if (const auto genIt = m_genParticleIds.find(genId); genIt != m_genParticleIds.end()) {
        m_recoToMCNav->addRelation(cand, genIt->second, 1.0);
      }

      m_recoParticleGenIds.emplace(genId, cand);
    }

  }
}

void DelphesLCIOConverter::processJets(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt, std::string const& branch) {
  auto* jetCollection = new lcio::LCCollectionVec(lcio::LCIO::RECONSTRUCTEDPARTICLE);
  evt->addCollection(jetCollection, branch);

  lcio::PIDHandler pidHandler(jetCollection);
  auto jetParamId = pidHandler.addAlgorithm("JetParameters", {"BTag", "TauTag"});

  for (auto iCand = 0; iCand < delphesCollection->GetEntries(); ++iCand) {
    auto* delphesCand = static_cast<Jet*>(delphesCollection->At(iCand));
    auto* jet = new lcio::ReconstructedParticleImpl;
    jetCollection->addElement(jet);

    // NOTE: Filling the jet with the information delievered by Delphes, which
    // is not necessarily the same as the sum of its constituents (filled below)
    jet->setCharge(delphesCand->Charge);
    jet->setMass(delphesCand->Mass);
    const auto momentum = delphesCand->P4();
    jet->setEnergy(momentum.E());
    const float mom[3] = {(float) momentum.Px(), (float) momentum.Py(), (float) momentum.Pz()};
    jet->setMomentum(mom);

    pidHandler.setParticleID(jet, 0, 0, 1., jetParamId, {
      (float) delphesCand->BTag,
      (float) delphesCand->TauTag
      });

    const auto& constituents = delphesCand->Constituents;
    for (auto iConst = 0; iConst < constituents.GetEntries(); ++iConst) {
      // TODO: Can we do better than Candidate here?
      auto* constituent = static_cast<Candidate*>(constituents.At(iConst));
      if (auto* matchedReco = getMatchingReco(constituent)) {
        jet->addParticle(matchedReco);
      } else {
        std::cerr << "**** WARNING: No matching ReconstructedParticle was found for a Jet constituent" << std::endl;
      }
    }
  }

}

void DelphesLCIOConverter::processMissingET(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt, std::string const& branch) {
  auto* collection = new lcio::LCCollectionVec(lcio::LCIO::RECONSTRUCTEDPARTICLE);
  evt->addCollection(collection, branch);

  auto* delphesCand = static_cast<MissingET*>(delphesCollection->At(0));
  auto* cand = new lcio::ReconstructedParticleImpl;
  collection->addElement(cand);
  const auto momentum = delphesCand->P4(); // NOTE: assuming massless here!
  const float mom[3] = {(float) momentum.Px(), (float) momentum.Py(), (float) momentum.Pz()};
  cand->setMomentum(mom);
  cand->setEnergy(momentum.E());
}

void DelphesLCIOConverter::processScalarHT(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt, std::string const& branch) {
  auto* collection = new lcio::LCCollectionVec(lcio::LCIO::PARTICLEID);
  evt->addCollection(collection, branch);

  auto* delphesCand = static_cast<ScalarHT*>(delphesCollection->At(0));

  auto* cand = new lcio::ParticleIDImpl;
  collection->addElement(cand);
  cand->addParameter(delphesCand->HT);
}

template<typename DelphesT>
void DelphesLCIOConverter::fillReferenceCollection(const TClonesArray* delphesCollection, lcio::LCEventImpl* evt, std::string const& branch, const std::string_view type) {
  auto* collection = new lcio::LCCollectionVec(lcio::LCIO::RECONSTRUCTEDPARTICLE);
  collection->setSubset(true);
  evt->addCollection(collection, branch);

  for (auto iCand = 0; iCand < delphesCollection->GetEntries(); ++iCand) {
    auto* delphesCand = static_cast<DelphesT*>(delphesCollection->At(iCand));

    if (auto* matchedReco = getMatchingReco(delphesCand)) {
      collection->addElement(matchedReco);

      // if we have an electron or muon we update the mass as well here
      if constexpr (std::is_same_v<DelphesT, Muon>) {
        matchedReco->setMass(M_MU);
      } else if constexpr (std::is_same_v<DelphesT, Electron>) {
        matchedReco->setMass(M_ELECTRON);
      }

      // If we have a charge available, also set it
      if constexpr (!std::is_same_v<DelphesT, Photon>) {
        matchedReco->setCharge(delphesCand->Charge);
      }

    } else {
      std::cerr << "**** WARNING: No matching ReconstructedParticle was found for a Delphes " << type << std::endl;
    }
  }
}



template<size_t N>
void sortBranchesProcessingOrder(std::vector<BranchSettings>& branches,
                                 std::array<std::string_view, N> const& processingOrder)
{
  std::unordered_map<std::string_view, size_t> classSortIndices;
  for (size_t i = 0; i < processingOrder.size(); ++i) {
    classSortIndices.emplace(processingOrder[i], i);
  }

  const auto endIt = classSortIndices.end();
  std::sort(branches.begin(), branches.end(),
            [&classSortIndices, endIt] (const auto& left, const auto& right) {
              const auto leftIt = classSortIndices.find(left.className);
              const auto rightIt = classSortIndices.find(right.className);

              // if we have the class defined in the processing order use the
              // corresponding index, otherwise use one definitely not inside
              // the processing order
              return (leftIt != endIt ? leftIt->second : N + 1) < (rightIt != endIt ? rightIt->second : N + 1);
            });
}


void setMotherDaughterRelations(GenParticle const* delphesCand,
                                lcio::MCParticleImpl* particle,
                                lcio::LCCollectionVec* mcParticles)
{
  // NOTE: it is in general probably not possible to handle all the different
  // possibilities that are present in the different readers. So, for now we are
  // going to follow the pythia documentation (with some additional sanity
  // checks, to avoid some crashs, in case the input is buggy) in the hope that
  // it will work for most inputs
  // Pythia documentation: http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html

  // Mothers
  // Only set parents if not accessing out of bounds
  const auto safeSetParent = [&mcParticles, &particle] (int index) {
    if (index < mcParticles->size()) {
      particle->addParent(static_cast<lcio::MCParticle*>(mcParticles->getElementAt(index)));
    }
  };

  // If M1 == -1, then this particle has no mother, so we only handle cases
  // where there is at least one
  if (delphesCand->M1 > -1) {
    // case 3, only one mother
    if (delphesCand->M2 == -1) {
      safeSetParent(delphesCand->M1);
    }
    if (delphesCand->M2 > -1){
      // case 6, two mothers
      if (delphesCand->M2 < delphesCand->M1) {
        safeSetParent(delphesCand->M1);
        safeSetParent(delphesCand->M2);

      } else {
        //  cases 2, 5 (and 4 without checking the status), mothers in a range
        for (auto iMother = delphesCand->M1; iMother <= delphesCand->M2; ++iMother) {
          safeSetParent(iMother);
        }
      }
    }
  }

  // Daughters are already handled by addParent in LCIO
}

lcio::TrackImpl* convertTrack(Track const* cand, const double magFieldBz) {
  auto* track = new lcio::TrackImpl;
  auto* trackState = new lcio::TrackStateImpl;

  // Delphes does not really provide any information that would go into the
  // track itself. But some information can be used to at least partially
  // populate a TrackState
  trackState->setD0(cand->D0);
  trackState->setZ0(cand->DZ);

  // Delphes calculates this from the momentum 4-vector at the track
  // perigee so this should be what we want. Note that this value
  // only undergoes smearing in the TrackSmearing module but e.g.
  // not in the MomentumSmearing module
  trackState->setPhi(cand->Phi);
  // Same thing under different name in Delphes
  trackState->setTanLambda(cand->CtgTheta);

  // Only do omega when there is actually a magnetic field.
  double varOmega = 0;
  if (magFieldBz) {
    // conversion to have omega in [1/mm]
    constexpr double a = c_light * 1e3 * 1e-15;

    trackState->setOmega(a * magFieldBz / cand->PT * std::copysign(1.0, cand->Charge));
    // calculate variation using simple error propagation, assuming
    // constant B-field -> relative error on pT is relative error on omega
    const double omega = trackState->getOmega();
    varOmega = cand->ErrorPT * cand->ErrorPT / cand->PT / cand->PT * omega * omega;
  }

  lcio::FloatVec covMatrix(15);

  // TODO: fix order of this or use full cov-matrix below
  covMatrix[0] = cand->ErrorD0 * cand->ErrorD0;
  covMatrix[5] = cand->ErrorPhi * cand->ErrorPhi;
  covMatrix[9] = varOmega;
  covMatrix[12] = cand->ErrorDZ * cand->ErrorDZ;
  covMatrix[14] = cand->ErrorCtgTheta * cand->ErrorCtgTheta;

  // TODO: need newer version of Delphes here

  // auto covaFB = cand->CovarianceMatrix();
  // covMatrix[1]  = covaFB(1,0) *scale1 * scale0;
  // covMatrix[2]  = covaFB(1,1) *scale1 * scale1;

  // covMatrix[3]  = covaFB(2,0) *scale2 * scale0;
  // covMatrix[4]  = covaFB(2,1) *scale2 * scale1;
  // covMatrix[5]  = covaFB(2,2) *scale2 * scale2;

  // covMatrix[6]  = covaFB(3,0) *scale3 * scale0;
  // covMatrix[7]  = covaFB(3,1) *scale3 * scale1;
  // covMatrix[8]  = covaFB(3,2) *scale3 * scale2;
  // covMatrix[9]  = covaFB(3,3) *scale3 * scale3;

  // covMatrix[10] = covaFB(4,0) *scale4 * scale0;
  // covMatrix[11] = covaFB(4,1) *scale4 * scale1;
  // covMatrix[12] = covaFB(4,2) *scale4 * scale2;
  // covMatrix[13] = covaFB(4,3) *scale4 * scale3;
  // covMatrix[14] = covaFB(4,4) *scale4 * scale4;

  trackState->setCovMatrix(covMatrix);
  track->addTrackState(trackState);

  return track;
}

template<typename DelphesT>
lcio::ReconstructedParticleImpl* DelphesLCIOConverter::getMatchingReco(DelphesT* delphesCand) const
{
  // Here we have to do some work to actually match the Delphes candidate to
  // the correct edm4hep::ReconstructedParticle because the possibility exists
  // that more than one ReconstructedParticle point to the same UniqueID. Here
  // we are NOT interested in the physics interpretation of such a case, but
  // only want to identify the correct ReconstructedParticle to which we
  // should point. To do the matching we compare the 4-momenta of the stored
  // edm4hep::ReconstructedParticle associated to the GenParticle UniqueID and
  // take the FIRST good match. Since the delphes candidate originates from
  // either a Track or a Tower, there should always be exactly one such good
  // match.
  for (const auto genId : getAllParticleIDs(delphesCand)) {
    const auto [recoBegin, recoEnd] = m_recoParticleGenIds.equal_range(genId);
    for (auto it = recoBegin; it != recoEnd; ++it) {
      // Handling slightly different member names for delphes depending on
      // whether we are still working with Candidates or the actual output
      // classes already
      if constexpr(std::is_same_v<DelphesT, Candidate>) {
        if (equalP4(getP4(it->second), delphesCand->Momentum)) {
          return it->second;
        } else if (equalP4(getP4(it->second), delphesCand->Momentum, 1e-2, false)) {
          // std::cout << "**** DEBUG: Kinematic matching successful after dropping energy matching and dropping momentum matching to percent level" << std::endl;
          return it->second;
        }
      } else {
        if (equalP4(getP4(it->second), delphesCand->P4())) {
          return it->second;
        }
      }
    }
  }

  return nullptr;
}


}
