#include "LKLogger.h"

void run_conv()
{
    lk_logger("data/run_conv.log");

    auto run = new LKRun();
    run -> AddPar("config_conv.mac");
    run -> AddDetector(new TexAT2());
    run -> Add(new LKMFMConversionTask());

    run -> Init();
    run -> Run();
}
