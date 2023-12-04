#include <iostream>
using namespace std;

#include "LKMFMConversionTask.h"
#include "MMChannel.h"

#include "GSpectra.h"
#include "GNetServerRoot.h"

LKMFMConversionTask::LKMFMConversionTask()
:LKTask("LKMFMConversionTask","")
{
}

bool LKMFMConversionTask::Init() 
{
    fChannelArray = fRun -> RegisterBranchA("RawData", "MMChannel", 1000);

    fMode             = fPar -> GetParInt("RunMode");
    fReadtype         = fPar -> GetParInt("ReadType");
    fConverterPort    = fPar -> GetParInt("ConverterPort");
    fBucketSize       = fPar -> GetParInt("fBucketSize");
    fScalerMode       = fPar -> GetParInt("fScalerMode");
    fD2pMode          = fPar -> GetParInt("2pMode");
    fUpdatefast       = fPar -> GetParInt("UpdateFast");
    //numfiles         = fPar -> GetParInt("NumberofFiles");
    //watcherPort      = fPar -> GetParInt("watcherPort");
    //CoBoServerPort   = fPar -> GetParInt("CoBoServerPort");
    //MutantServerPort = fPar -> GetParInt("MutantServerPort");
    //energymethod     = fPar -> GetParInt("EnergyFindingMethod");
    //readrw           = fPar -> GetParInt("ReadResponseWaveformFlag");
    //RootConvert      = fPar -> GetParInt("RootConvertEnable");
    //DrawWaveform     = fPar -> GetParInt("DrawWaveformEnable");
    //cleantrack       = fPar -> GetParInt("CleanTrackEnable");
    //DrawTrack        = fPar -> GetParInt("DrawTrackEnable");
    //SkipEvents       = fPar -> GetParInt("SkipEvents");
    //firstEventNo     = fPar -> GetParInt("firstEventNo");
    //IgnoreMM         = fPar -> GetParBool("IgnoreMicromegas");
    //mfmfilename      = fPar -> GetParString("MFMFileName");
    //watcherIP        = fPar -> GetParString("watcherIP");
    //mapChanToMM      = fPar -> GetParString("ChanToMMMapFileName");
    //mapChanToSi      = fPar -> GetParString("ChanToSiMapFileName");
    //rwfilename       = fPar -> GetParString("ResponseWaveformFileName");
    //supdatefast      = fPar -> GetParString("UpdateFast");

    fFrameBuilder = new LKFrameBuilder(fConverterPort);
    fFrameBuilder -> SetBucketSize(fBucketSize);
    fFrameBuilder -> Init(fMode,fD2pMode);
    fFrameBuilder -> SetReadMode(fMode);
    fFrameBuilder -> SetReadType(fReadtype);
    fFrameBuilder -> SetScaler(fScalerMode);
    fFrameBuilder -> Set2pMode(fD2pMode);
    fFrameBuilder -> SetUpdateSpeed(fUpdatefast);
    fFrameBuilder -> SetChannelArray(fChannelArray);

    fBuffer = (char *) malloc (matrixSize);
    lk_info << "Opening file stream." << endl;
    fFileStream.open(infname.c_str(),std::ios::binary | std::ios::in);
    if(!fFileStream) {
        lk_error << "Could not open input file!" << std::endl;
        return false;
    }
    lk_info << "Reading block size is " << matrixSize << endl;

    return true;
}

void LKMFMConversionTask::Exec(Option_t*)
{
    if (fFileStream.eof()) {
        lk_warning << "end of MFM file!" << endl;
        fRun -> SignalEndOfRun();
        return;
    }

    int filebuffer = 0;
    while (true)
    {
        fFileStream.seekg(filebuffer, std::ios_base::cur);
        fFileStream.read(fBuffer,matrixSize);
        filebuffer += matrixSize;

        if(!fFileStream.eof()) {
            try {
                ++fCountAddDataChunk;
                fFrameBuilder -> addDataChunk(fBuffer,fBuffer+matrixSize);
            }catch (const std::exception& e){
                lk_error << "Error occured from " << fCountAddDataChunk << "-th addDataChunk()" << endl;
                e_cout << e.what() << endl;
                fRun -> SignalEndOfRun();
                return;
            }
        }
        else if(fFileStream.gcount()>0) {
            try {
                ++fCountAddDataChunk;
                fFrameBuilder -> addDataChunk(fBuffer,fBuffer+fFileStream.gcount());
            }catch (const std::exception& e){
                lk_error << "Error occured from LAST " << fCountAddDataChunk << "-th addDataChunk()" << endl;
                e_cout << e.what() << endl;
                return;
            }
        }
    }
}

bool LKMFMConversionTask::EndOfRun()
{
    return true;
}
