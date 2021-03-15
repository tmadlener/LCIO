#ifndef BENCHMARK_INSTRUMENTATION
#define BENCHMARK_INSTRUMENTATION 0
#endif

#include "DelphesLCIOConverter_new.h"
#include "DelphesLCIOOutputConfiguration.h"
#include "DelphesInputReader.h"

#include "IOIMPL/LCFactory.h"

#include "podio/BenchmarkRecorder.h"

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
  using namespace podio::benchmark;

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

#if BENCHMARK_INSTRUMENTATION
    BenchmarkRecorder benchmarkRecorder(outputFile + ".bench.root");
    auto& setupTree = benchmarkRecorder.addTree("setup_times", {"constructor", "open_file", "close"});
    auto& eventTree = benchmarkRecorder.addTree("event_times", {"write_event"});
    auto& loopTree = benchmarkRecorder.addTree("non_io_times",
                                               {"loop", "read", "process_delphes", "process_convert", "write"});

    const auto constStart = ClockT::now();
#endif
    auto lcWriter = std::unique_ptr<lcio::LCWriter>(lcio::LCFactory::getInstance()->createLCWriter());
#if BENCHMARK_INSTRUMENTATION
    const auto constEnd = ClockT::now();
    setupTree.recordTime("constructor", constEnd - constStart);
    // Cannot figure out how to use run_void_member_timed with LCWriter here,
    // because open is overloaded
    const auto openStart = ClockT::now();
#endif
    lcWriter->open(outputFile, lcio::LCIO::WRITE_NEW);
#if BENCHMARK_INSTRUMENTATION
    const auto openEnd = ClockT::now();
    setupTree.recordTime("open_file", openEnd - openStart);
#endif
   
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
#if BENCHMARK_INSTRUMENTATION
      const auto loopStartTime = podio::benchmark::ClockT::now();
#endif
      if (!inputReader.readEvent(modularDelphes,
                                 allParticleOutputArray,
                                 stableParticleOutputArray,
                                 partonOutputArray)) {
        break;
      }
#if BENCHMARK_INSTRUMENTATION
      const auto readEndTime = podio::benchmark::ClockT::now();
#endif
      modularDelphes->ProcessTask();
#if BENCHMARK_INSTRUMENTATION
      const auto delphesProcessEndTime = podio::benchmark::ClockT::now();
#endif
      auto evt = std::make_unique<lcio::LCEventImpl>();
      converter.process(inputReader.converterTree(), evt.get());
#if BENCHMARK_INSTRUMENTATION
      const auto convertEndTime = podio::benchmark::ClockT::now();
#endif

#if BENCHMARK_INSTRUMENTATION
      eventTree.recordTime("write_event",
                           run_void_member_timed(*lcWriter, &lcio::LCWriter::writeEvent, evt.get()));
#else
      lcWriter->writeEvent(evt.get());
#endif
      modularDelphes->Clear();
      progressBar.Update(eventCounter, eventCounter);
      eventCounter++;
#if BENCHMARK_INSTRUMENTATION
      const auto loopEndTime = podio::benchmark::ClockT::now();
      loopTree.recordTime("loop", loopEndTime - loopStartTime);
      loopTree.recordTime("read", readEndTime - loopStartTime);
      loopTree.recordTime("process_delphes", delphesProcessEndTime - readEndTime);
      loopTree.recordTime("process_convert", convertEndTime - delphesProcessEndTime);
      loopTree.recordTime("write", loopEndTime - convertEndTime);
      loopTree.Fill();
      eventTree.Fill();
#endif
    }

    progressBar.Update(eventCounter, eventCounter, true);
    progressBar.Finish();
    modularDelphes->Finish();
#if BENCHMARK_INSTRUMENTATION
    setupTree.recordTime("close", run_void_member_timed(*lcWriter, &lcio::LCWriter::close));
    setupTree.Fill();
#else
    lcWriter->close();
#endif
  } catch (std::runtime_error& e) {
    std::cerr << "** ERROR: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
