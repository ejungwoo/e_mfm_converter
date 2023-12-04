#ifndef ATMFMCONVERSIONTASK_HH
#define ATMFMCONVERSIONTASK_HH

#include "LKRun.h"
#include "LKTask.h"
#include "LKFrameBuilder.h"

/*
 * AT-TPC MFM conversion class
 * Imported from MFMHistServer
 */
class LKMFMConversionTask : public LKTask
{
    public:
        LKMFMConversionTask();
        virtual ~LKMFMConversionTask() {};

        bool Init();
        void Exec(Option_t*);
        bool EndOfRun();

        const int kOnline = 0;
        const int kReadMFM = 1;
        const int kReadList = 3;
        const int kReadExpNameNo = 11;
        const int kReadType2 = 2;
        const int kReadType4 = 4;
        const int kReadType10 = 10;
        const int kReadType13 = 13;
        const int kReadType14 = 13;

    private:
        int maxit = 0;
        int currit = 0;
        int percent = 0;
        int size_buffer;
        char *fBuffer;
        int oflcnt=0;

        //size_t const matrixSize = 4*68*512*sizeof(double);
        size_t const matrixSize = 512;

        // init parameters
        int fMode = 1;
        int fReadtype;
        //int numfiles;
        //int watcherPort;
        int fConverterPort;
        //int CoBoServerPort;
        //int MutantServerPort;
        int fBucketSize;
        //int energymethod;
        //int readrw;
        //int RootConvert;
        int fScalerMode;
        int fD2pMode;
        int fUpdatefast;
        //int DrawWaveform;
        //int cleantrack;
        //int DrawTrack;
        //int SkipEvents;
        //int firstEventNo;
        //bool IgnoreMM;
        //TString mfmfilename;
        //TString watcherIP;
        //TString mapChanToMM;
        //TString mapChanToSi;
        //TString rwfilename;
        //TString supdatefast;
        //TString goodEventList;

        LKFrameBuilder* fFrameBuilder; // convServer

        // global parameters in main.cc
        //static HistServer* histServer;
        //static int rsize_buffer=512;
        static string inrfname;
        static int maxinrfidx;
        static const int maxfileno=2000;
        static string listinrfname[maxfileno];
        static string inrtname;
        static string inoutrfname; // output root file name
        static string listinoutrfname[maxfileno];
        static string maxinoutrfidx;
        static string inoutrtname;
        //static char keyinput='r';
        //static thread tgetkey;

        // parameter related to files
        string infname;
        string minfname;
        string outrfname;
        string outrtname;
        string slistinfname;
        string listinfname[maxfileno];
        int maxinfidx=0;
        int infidx=0;

        //TString infname;

        ifstream fFileStream;

        int fCountAddDataChunk = 0;

        TClonesArray *fChannelArray = nullptr;

    ClassDef(LKMFMConversionTask, 1)
};

#endif
