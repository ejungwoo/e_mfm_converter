#ifndef ATMFMCONVERSIONTASK_HH
#define ATMFMCONVERSIONTASK_HH

#include "LKTask.hpp"
#include "ATMFMFrameBuilder.hpp"

/*
 * AT-TPC MFM conversion class
 * Imported from MFMHistServer
 */
class ATMFMConversionTask : public LKTask
{
    public:
        ATMFMConversionTask();
        virtual ~ATMFMConversionTask() {};

        enum ReadType {
            kOnline = 0,
            kReadMFM = 1,
            kReadList = 3,
            kReadExpNameNo = 11,
            kReadType2 = 2,
            kReadType4 = 4,
            kReadType10 = 10,
            kReadType13 = 13,
            kReadType14 = 13,
        };

        bool Init();
        void Exec(Option_t*);
        bool EndOfRun();

    private:
        bool   checkport(int usPort);
        string _executeShellCommand(string command);

        int maxit = 0;
        int currit = 0;
        int percent = 0;
        int size_buffer;
        char *fBuffer;
        int oflcnt=0;

        // init parameters
        int mode = 1;
        ReadType readtype;
        int numfiles;
        int watcherPort;
        int converterPort;
        int CoBoServerPort;
        int MutantServerPort;
        int BucketSize;
        int energymethod;
        int readrw;
        //int RootConvert;
        int ScalerMode;
        int d2pMode;
        int updatefast;
        int DrawWaveform;
        int cleantrack;
        int DrawTrack;
        int SkipEvents;
        int firstEventNo;
        bool IgnoreMM;
        TString mfmfilename;
        TString watcherIP;
        TString mapChanToMM;
        TString mapChanToSi;
        TString rwfilename;
        TString supdatefast;
        TString goodEventList;

        // global parameters in main.cc
        static HistServer* histServer;
        static HistServer* fFrameBuilder; // convServer
        static int rsize_buffer=512;
        static string inrfname;
        static int maxinrfidx;
        static const int maxfileno=2000;
        static string listinrfname[maxfileno];
        static string inrtname;
        static string inoutrfname; // output root file name
        static string listinoutrfname[maxfileno];
        static string maxinoutrfidx;
        static string inoutrtname;
        static char keyinput='r';
        static thread tgetkey;

        // parameter related to files
        int currfiles = 1;
        string infname;
        string minfname;
        string outrfname;
        string outrtname;
        string slistinfname;
        string listinfname[maxfileno];
        int maxinfidx=0;
        int infidx=0;

        TString infname;

        ifstream fFileStreamForEventLoop;

        TClonesArray *fChannelArray = nullptr;

    ClassDef(ATMFMConversionTask, 1)
};

#endif
