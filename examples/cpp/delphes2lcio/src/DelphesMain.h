#include "DelphesLCIOConverter_new.h"
#include "DelphesLCIOOutputConfiguration.h"
#include "DelphesInputReader.h"

#include "IOIMPL/LCFactory.h"


#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootConfReader.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

#include <csignal>
#include <iostream>
#include <stdexcept>
#include <memory>

static bool interrupted = false;
void SignalHandler(int /*si*/) {
  interrupted = true;
}


int runConverter(int argc, char* argv[], DelphesInputReader& inputReader) {
  using namespace delphes_lcio;

  // We can't make this a unique_ptr because it interferes with whatever ROOT is
  // doing under the hood to clean up
  auto* modularDelphes = new Delphes("Delphes");
  const auto outputFile = inputReader.init(modularDelphes, argc, argv);
  if (outputFile.empty()) {
    std::cerr << inputReader.getUsage() << std::endl;
    return 1;
  }

  signal(SIGINT, SignalHandler);
  try {
    auto confReader = std::make_unique<ExRootConfReader>();
    confReader->ReadFile(argv[1]);
    modularDelphes->SetConfReader(confReader.get());

    const auto branches = getBranchSettings(confReader->GetParam("TreeWriter::Branch"));
    const auto lcioOutputSettings = getEDM4hepOutputSettings(argv[2]);
    DelphesLCIOConverter converter(branches,
                                   lcioOutputSettings,
                                   confReader->GetDouble("ParticlePropagator::Bz", 0));

    auto lcWriter = std::unique_ptr<lcio::LCWriter>(lcio::LCFactory::getInstance()->createLCWriter());
    lcWriter->open(outputFile, lcio::LCIO::WRITE_NEW);

    // has to happen before InitTask
    TObjArray* allParticleOutputArray = modularDelphes->ExportArray("allParticles");
    TObjArray* stableParticleOutputArray = modularDelphes->ExportArray("stableParticles");
    TObjArray* partonOutputArray = modularDelphes->ExportArray("partons");

    modularDelphes->InitTask();
    modularDelphes->Clear();

    const int maxEvents = confReader->GetInt("::MaxEvents", 0);
    ExRootProgressBar progressBar(-1);
    Int_t eventCounter = 0;
    for (Int_t entry = 0;
         !inputReader.finished() && (maxEvents > 0 ?  entry < maxEvents : true) && !interrupted;
         ++entry) {

      if (!inputReader.readEvent(modularDelphes,
                                 allParticleOutputArray,
                                 stableParticleOutputArray,
                                 partonOutputArray)) {
        break;
      }

      auto* evt = new lcio::LCEventImpl;

      modularDelphes->ProcessTask();

      converter.process(inputReader.converterTree(), evt);
      lcWriter->writeEvent(evt);
      delete evt;

      modularDelphes->Clear();
      progressBar.Update(eventCounter, eventCounter);
      eventCounter++;
    }

    progressBar.Update(eventCounter, eventCounter, true);
    progressBar.Finish();
    modularDelphes->Finish();
    lcWriter->close();
  } catch (std::runtime_error& e) {
    std::cerr << "** ERROR: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
