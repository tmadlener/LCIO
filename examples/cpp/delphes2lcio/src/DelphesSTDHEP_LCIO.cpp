#include "DelphesSTDHEPInputReader.h"
#include "DelphesMain.h"

int main(int argc, char* argv[]) {
  DelphesSTDHEPInputReader inputReader{};
  runConverter(argc, argv, inputReader);
}
