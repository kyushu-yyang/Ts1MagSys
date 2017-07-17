#include <iostream>
#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XCOMETConstruction.h"

using namespace Quench;

int main(int argc, char** argv)
{
 
  XQuenchLogger* log = XQuenchLogger::GetInstance();
  log->Start(XQuenchLogger::DEBUG, "quenchlogger.log");
  //

  XCOMETConstruction* comet = new XCOMETConstruction();

  try {
    comet->SetTime(0.*sec, 100.*sec, 1.e-3*msec);
    //comet->SetDisplayStep(10);
    comet->SetCurrent(2700.*Amp);
    comet->SetHotSpot(40/2, 1, 2);
    comet->SetDumpResistor(0.185*Ohm);
    comet->SetInductance(12.69);
    comet->SetThreshold(0.1);
    comet->SetDetectTime(0.1*sec);

    comet->ConstructTs1a();
    comet->ConstructTs1b();
    comet->Begin();
    comet->Run();
    comet->End();
  }
  catch (XQuenchExcept except) {
    std::cerr << "Error: " << except.what() << std::endl;
    delete comet;
  }

  return 0;
}
