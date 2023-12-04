#ifndef LKFRAMEBUILDER_H
#define LKFRAMEBUILDER_H

#include "mfm/FrameBuilder.h"
#include <map>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TTree.h>
#include <TKey.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TGraph2D.h>
#include <TProfile.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TF2.h>
#include <TRandom.h>
#include <TMath.h>
#include <TCutG.h>
#include <string.h>
#include <TError.h>
#include <fstream>

using namespace std;
class GSpectra;
class GNetServerRoot;

class WaveForms {
    public:
        WaveForms();
        ~WaveForms();
        Bool_t isRejected;
        vector<vector<Bool_t>> hasSignal; // To let us know if there was a signal that fired
        vector<Bool_t> hasHit; // To let us know if there was any signal other than FPNs
        vector<Bool_t> hasFPN; // To let us know if there was any FPN signal
        vector<Bool_t> doneFPN; // To let us know if there was any FPN signal averaged
        vector<vector<Bool_t>> isOverflow; // To let us know if the signal was overflow.
        Bool_t hasOverflow; // To let us know if the event has overflow channels.,
        vector<vector<Int_t>> isDecay;
        Bool_t hasImplant;
        Bool_t hasDecay;
        UInt_t frameIdx;
        UInt_t decayIdx;
        UInt_t ICenergy;
        UInt_t Sienergy;
        Double_t SiX;
        Double_t SiY;
        Double_t SiZ;
        UInt_t EstripL;
        UInt_t EstripR;
        UInt_t coboIdx;
        UInt_t asadIdx;
        vector<vector<vector<UInt_t>>> waveform; // To save the waveform for each Aget, Channel and Bucket
        vector<vector<vector<Int_t>>> corrwaveform; // To save the corrected waveform by the averaged FPN waveform
        vector<vector<UInt_t>> fpnwaveform; // To save the waveform for each Aget, FPN Channel and Bucket
        vector<vector<UInt_t>> energy; // Digitized energy value
        vector<vector<UInt_t>> time; // time value
        vector<vector<UInt_t>> baseline; // baseline value
        vector<vector<Double_t>> PSDIntegral; // To save ratio of peak to integral.
        vector<vector<Double_t>> PSDRatio; // To save ratio of peak to integral.
};

class MM_Track {
    public:
        MM_Track();
        ~MM_Track();
        UInt_t hasTrack; // total number of hits in track
        UInt_t hasTrackChain; // total number of chain hits in track
        UInt_t hasTrackStrip; // total number of strip hits in track
        Bool_t hasProjectile; // incoming projectile flag 
        Bool_t hasEjectile; // outgoing ejectile flag 
        Bool_t hasRecoil; // outgoing recoil flag 
        Bool_t hasOverflow; // overflow flag
        Bool_t hasImplant; // implant flag
        Bool_t hasDecay; // decay flag
        Bool_t L1B[150][170];//Is L1B
        UInt_t ICenergy; //IC energy
        UInt_t Sienergy; //Silicon energy
        Double_t SiX;
        Double_t SiY;
        Double_t SiZ;
        UInt_t decay[150][170]; // decay flag by pixel X and Y, 0=implant, 1=decay
        UInt_t pixel[150][170]; // hit pattern by pixel X and Y
        UInt_t energy[150][170]; // Digitized energy value
        UInt_t time[150][170]; // time value
        Double_t avgposx[170]; // Average x position value
        Double_t avgposy[170]; // Average y position value
        Double_t sumenergy[170]; // Average energy value
        Double_t sum2penergy[4]; // Sum of all 2p energies, 0=bc,1=bl,2=br,4=all
        UInt_t coloridx[150][170]; // color value
        Double_t posx[150][170]; // Position X by pixel X
        Double_t posy[150][170]; // Position Y by pixel Y
        Double_t posz[150][170]; // Position Z by time and drift velocity
        Double_t posxerr[150][170]; // Position X Error by pixel X
        Double_t posyerr[150][170]; // Position Y Error by pixel Y
        Double_t poszerr[150][170]; // Position Z Error by time and drift velocity
};

class Si_Track {
    public:
        Si_Track();
        ~Si_Track();
        UInt_t hasTrack; // track flag
        UInt_t hasFWC; // Forward Si Center flag
        UInt_t hasFWL; // Forward Si Left flag
        UInt_t hasFWR; // Forward Si Right flag
        UInt_t hasX6L; // X6 Left flag
        UInt_t hasX6R; // X6 Right flag
        UInt_t hasX6BL; // X6 Bottom Left flag
        UInt_t hasX6BR; // X6 Bottom Right flag
        Bool_t hasProjectile; // incoming projectile flag 
        Bool_t hasEjectile; // outgoing ejectile flag 
        Bool_t hasRecoil; // outgoing recoil flag 
        UInt_t pixel[150][170]; // hit pattern by pixel X and Y
        UInt_t energy[150][170]; // Digitized energy value
        UInt_t time[150][170]; // time value
        UInt_t coloridx[150][170]; // color value
        Double_t posx[150][170]; // Position X by pixel X
        Double_t posy[150][170]; // Position Y by pixel Y
        Double_t posz[150][170]; // Position Z by time and drift velocity
        Double_t posxerr[150][170]; // Position X Error by pixel X
        Double_t posyerr[150][170]; // Position Y Error by pixel Y
        Double_t poszerr[150][170]; // Position Z Error by time and drift velocity
        UInt_t hasICE; // IC energy
        UInt_t hasICT; // IC time
        UInt_t hasCsIE[9]; // CsI energy
        UInt_t hasCsIT[9]; // CsI time
        UInt_t agetid; // agetid
        UInt_t chanid; // chanid
};

class MapChanToMM {
    public:
        MapChanToMM();
        ~MapChanToMM();
        UInt_t pxidx[4][4][64];
        UInt_t pxidy[4][4][64];
        Double_t pxposx[4][4][64];
        Double_t pxposy[4][4][64];
};

class MapChanToSi {
    public:
        MapChanToSi();
        ~MapChanToSi();
        UInt_t pxidx[4][4][64];
        UInt_t pxidy[4][4][64];
        Double_t pxposx[4][4][64];
        Double_t pxposy[4][4][64];
};

class MapChanToX6 {
    public:
        MapChanToX6();
        ~MapChanToX6();
        UInt_t X6flag[4][4][68];
        UInt_t X6det[4][4][68];
        UInt_t X6strip[4][4][68];
        Double_t X6posx[300][8]; //[det][strip]
        Double_t X6posy[300][8]; //[det][strip]
        Double_t X6posz[300][8]; //[det][strip]

        UInt_t CsICTnum[4][4][68];
        UInt_t CsI_X6det[4][4][68];
        UInt_t CsI_X6ud[4][4][68];
};

class LKFrameBuilder : public mfm::FrameBuilder {
    public:
        void SetChannelArray(TClonesArray *channelArray) { fChannelArray = channelArray; }

    private:
        TClonesArray *fChannelArray = nullptr;

    public:
        LKFrameBuilder(int);
        ~LKFrameBuilder();
        void processFrame(mfm::Frame & frame);
        void InitWaveforms();
        void ResetWaveforms();
        void Init(int mode, int d2pmode); //converting mode=1, histogram mode=2
        void RootWOpenFile(string & outputFileName, string & outputTreeName);
        void SetUpdateSpeed(int updatefast);
        void RootWInit();
        void RootWReset();
        void RootWConvert();
        void RootWriteEvent();
        void RootWCloseFile();
        void RootROpenFile(string & inputFileName, string & inputTreeName);
        void RootRInit();
        void RootRInitWaveforms();
        void RootFindEvent();
        void RootReadEvent();
        void RootRResetWaveforms();
        void RootRReset();
        void RootRCloseFile();
        bool RootRRecovered();
        void RootRWOpenFile(string & outputFileName, string & outputTreeName);
        void RootReadWriteEvent();
        void RootRWInit();
        void RootRWriteEvent();
        void RootRWCloseFile();
        void RootRWReset();
        void RootRHInit();
        void RootRHInitWaveforms();
        void RootRHResetWaveforms();
        void RootRHReset();
        void RootReadHistEvent();
        UInt_t GetRootFileSize();
        void SetRootConverter(int flag);
        void SetIgnoreMM(int flag);
        void SetDrawWaveform(int flag);
        void SetCleanTrack(int flag);
        void SetDrawTrack(int flag);
        void SetSkipEvents(int flag);
        void SetfirstEventNo(int flag);
        void SetEnergyMethod(int flag);
        void SetReadRF(int flag);
        int  GetForceReadTree();
        void ValidateEvent(mfm::Frame & frame);
        void ValidateFrame(mfm::Frame& frame);
        void Event(mfm::Frame & frame);
        void UnpackFrame(mfm::Frame& frame);
        void SetReadMode(int flag);
        int GetReadMode();
        void SetReadType(int flag);
        void SetScaler(int flag);
        void Set2pMode(int flag);
        void SetHistMode();
        void decodeMuTanTFrame(mfm::Frame & frame);
        void decodeCoBoTopologyFrame(mfm::Frame& frame);
        void SetBucketSize(int BucketSize);
        void ReadMapChanToMM(string filename);
        void ReadMapChanToSi(string filename);
        void ReadMapChanToX6();
        void ReadGoodEventList(string filename);
        void ReadResponseWaveform(string filename);
        void SetResponseWaveform();
        Bool_t IsResponseSample(Int_t type, Int_t cobo, Int_t asad, Int_t aget, Int_t chan);
        void SetResponseSample(Int_t type, Int_t evtno, Int_t buckcut, Int_t buckwidth, Int_t cobo, Int_t asad, Int_t aget, Int_t chan, Int_t rep, Int_t iter, Int_t boost);
        void GetMaxResponseWaveform();
        void GetSigmaResponseWaveform();
        Double_t ShaperF_GET1(Double_t *x, Double_t *p);
        Double_t ShaperF_GET2(Double_t *x, Double_t *p);
        Double_t ShaperF_MSCF(Double_t *x, Double_t *p);
        void GetAverageFPN(Int_t decayIdx, Int_t cobo, Int_t asad, Int_t aget);
        void GetCorrWaveform(Int_t decayIdx, Int_t cobo, Int_t asad, Int_t aget, Int_t chan);
        Double_t GetBaseline(Int_t decayIdx, Int_t cobo, Int_t asad, Int_t aget, Int_t chan);
        void GetEnergyTime(Int_t decayIdx, Int_t cobo, Int_t asad, Int_t aget, Int_t chan);
        //void GetEnergybyFitWaveform(Int_t type, Int_t cobo, Int_t asad, Int_t aget, Int_t chan, Int_t maxAt, Int_t peakAt, Int_t *maxValuedec, Int_t *maxValueBucketdec); //using ratio method
        void GetEnergybyFitWaveform(Int_t type, Int_t decayIdx, Int_t cobo, Int_t asad, Int_t aget, Int_t chan, Int_t peakAt); //using ShaperF method
        void DrawWaveForm(Int_t decayIdx, Int_t cobo, Int_t asad, Int_t aget, Int_t chan);
        void FillTrack();
        void FindBoxCorner();
        void FindX6Hits();
        void DrawSiDetector();
        void ReplaceEnergy();
        void ReplaceEnergybyRatio();
        Int_t GetSumEnergy(Int_t px, Int_t py);
        void FillDecayFlag();
        void Sum2pEnergy();
        void ResetHitPattern();
        void DrawHitPattern(Int_t decayIdx, UInt_t cobo);
        void WaveformShapeFilter(Int_t decayIdx, UInt_t cobo);
        void WaveletFilter(Int_t decayIdx, UInt_t cobo);
        void DrawPSDFilter(Int_t decayIdx, UInt_t cobo);
        void DrawdEvsE();
        void DrawSiEvsCsIE();
        void CleanTrack();
        void FilldEvsE();
        void DrawTrack();
        void ChangeTrackHistTitle();
        void DrawSumEnergyTrack();
        void DrawTrack2pMode();
        void HoughTransform();
        void GetTrackPosYLimit();
        void GetXYZTrack();
        void InitTrack();
        void InitMutantScaler();
        void ResetTrack();
        void ResetTrackHist();
        void ResetdEvsE();
    protected:
        GSpectra* spectra_;
        GNetServerRoot* serv_;
        unsigned int triggerrate[2];
        int maxcobo;
        int maxasad;
        WaveForms* waveforms;
        WaveForms* rwaveforms[2][3];
        MM_Track* mm_tracks;
        MM_Track* mm_tracks_2p;
        Si_Track* si_tracks;
        MapChanToMM* mapchantomm;
        MapChanToSi* mapchantosi;
        MapChanToX6* mapchantox6;
        TH1D* hWaveForm[64];
        TH1D* hCorrWaveForm[64];
        TH1D* hCorrWaveFormDec[64];
        TH1D* hCorrWaveFormFit[64];
        TH2D* hCorrWaveFormRDF[64];
        TH2D* hFPNWaveFormAll[64];
        TH2D* hWaveFormbyEvent[64];
        TH2D* hCorrWaveFormbyEvent[64];
        TH2D* hWaveFormIC;
        TSpectrum* sCorrWaveForm[64];
        TH2D* hMM_TrackAll;
        TH2D* hMMICenergyvsMME;
        TH2D* hMM_EstripL;
        TH2D* hMM_EstripR;
        TH2D* hMM_EstripAll;
        TH2D* hMM_Track[64];
        TH2D* hMM_TrackDecay[64];
        TH2D* hMM_TrackDecay2[64];
        TGraph2D* gMM_TrackDecay[64];
        TH2D* hMM_TrackvsE[64];
        TH2D* hMM_TrackPosAll;
        TH2D* hMM_TrackPos[64];
        TGraph2D* gMM_TrackDecayPos[64];
        TH2D* hMM_TrackdEvsEALL[3];
        TH2D* hMM_TrackdE1vsSiE;
        TH2D* hMM_TrackdE1vsdE2;
        TH2D* hMM_TrackdEvsE[64][3];
        TH2D* hMM_TrackPosHough[64][3];
        TH2D* hMM_TrackXZ[64];
        TH2D* hMM_TrackYZ[64];
        TH2D* hMM_TrackPosXYAll;
        TH2D* hMM_TrackPosXYTAll;
        TH2D* hMM_TrackPosXZAll;
        TH2D* hMM_TrackPosYZAll;
        TH2D* hMM_TrackPosXY[64];
        TH2D* hMM_TrackPosXZ[64];
        TH2D* hMM_TrackPosYZ[64];
        TH2D* hMM_TimevsPxIDX[64];
        TH2D* hMM_TimevsPxIDXPosAll;
        TH2D* hMM_TimevsPxIDYAll;
        TH2D* hMM_TimevsPxIDXPos[64];
        TH2D* hMM_TimevsPxIDXPosHough[64][3];
        TH2D* hMM_TimevsPxIDY[64];
        TH2D* hMM_TimevsPxIDYPosAll;
        TH2D* hMM_TimevsPxIDYPos[64];
        TH2D* hMM_TimevsPxIDYPosHough[64][3];
        TH2D* hMM_EnergyvsPxIDYALL;
        TH2D* hMM_EnergyvsPxIDY[64];
        TH2D* hMM_SiEvsCsIEAll;
        TH2D* hMM_SiEvsCsIE[10][11];
        TH2D* hMM_X6EvsCsIEAll;
        TH2D* hMM_X6EvsCsIE[30];
        TH1I* hMM_SumEnergyvsPxIDY[64];
        TH1I* hMM_SumEnergyMaxIDY;
        TH1I* hMM_SumEnergyLastIDY;
        TH2D* hMM_PSDIntegralvsMax[64];
        TH2D* hMM_PSDRatiovsMax[64];
        TH2D* hMM_PSDRatio[64];
        TH2D* hMM_Pa_vs_Energy2D;
        TH2D* hMM_Pa_vs_Time2D;
        TH2D* hSi_Pa_vs_Energy2D;
        TH2D* hSi_Pa_vs_Time2D;
        TH1D* hMM_D2PTime;
        TH1D* hMM_Sum2pEnergy[4];
        TH1D* hMM_Time[64];
        TH1D* hMM_Energy[64];
        TH1D* hGET_EALL[64];
        TH1D* hGET_E[64][68];
        TH1I* hGET_HitPattern;
        TH2I* hGET_EHitPattern2D;
        TH2I* hGET_THitPattern2D;
        TH2I* hGET_Si16x16HitPattern2D;
        TH2I* hGET_SiForwardHitPattern2D;
        TH1I* hGET_SiForwardMult;
        TH2I* hGET_SivsCsIHitPattern2D;
        TH2I* hGET_DecEHitPattern2D;
        TH2I* hGET_DecTHitPattern2D;
        TH2I* hGET_FitEHitPattern2D;
        TH2I* hGET_FitTHitPattern2D;
        TH2I* hGET_DecdEMaxvsEMax2D;
        TH2I* hGET_DecdTMaxvsTMax2D;
        TH2I* hGET_FitdEMaxvsEMax2D;
        TH2I* hGET_FitdTMaxvsTMax2D;
        TH2D* hGET_EHitPattern[64];
        TH2D* hGET_THitPattern[64];
        TH2D* hGET_ERHitPattern[64];
        TH1D* hMM_a0[4];
        TH1D* hMM_a1[4];
        TH1D* hMM_b0[4];
        TH1D* hMM_b1[4];
        TH2D* hBM_LR;
        TH2D* hBM_Esum1vsPos;
        TH2D* hBM_Esum2vsEsum1;
        TH2D* hX6_LRAll;
        TH2D* hX6_EsumvsPosAll;
        TH2D* hX6_EsumvsPosAllLeft;
        TH2D* hX6_EsumvsPosAllRight;
        TH1D* hX6_EbackAll;
        TH2D* hX6_LR[48][8];
        TH2D* hX6_EsumvsPos[48][8];
        TH1D* hX6_Eback[48][4];
        TH2D* hX6_LeftHitPattern;
        TH2D* hX6_BottomHitPattern;
        TH2D* hX6_RightHitPattern;
        TH2D* hX6_OhmicHitPattern;
        TRandom* rpos;
        TCutG* cut_pinSiEvsCsIE;
        TCutG* cut_pinX6EvsCsIE;
        TCutG* cut_dtinSiEvsCsIE;
        TCutG* cut_inSiEvsCsIE;
        ofstream X6out;
        const static int maxevtno=100000000;
        bool evtmask[maxevtno];

        int bucketmax;
        int readmode;
        int readtype;
        int enableroot;
        int forcereadtree;
        int enablecleantrack;
        int energymethod;
        int readrw;
        int ignoremm;
        int enabledraww;
        int enableskipevent;
        int firsteventno;
        int enabletrack;
        int enablescaler;
        int enable2pmode;
        int enablehist;
        int enableupdatefast;
        int hasDrawn[16];
        UInt_t bmpos; //BM Position
        UInt_t bmsum1; //BM Energy by position ratio
        UInt_t bmsum2; //BM Energy by sig_sum
        UInt_t eventID; //weventIdx;
        UInt_t weventTime;
        UInt_t goodsievt;
        UInt_t goodmmevt;
        UInt_t prevgoodmmevt;
        UInt_t goodicevt;
        UInt_t prevgoodicevt;
        UInt_t goodevtidx;
        UInt_t goodsicsievt;
        UInt_t goodx6evt;
        UInt_t goodx6csievt;
        UInt_t goodsicsievtidx;
        UInt_t goodsicsipevt;
        UInt_t goodsicsipevtidx;
        UInt_t reventIdx;
        UInt_t reventTime;
        UInt_t rd2ptime;
        UInt_t rtstmp;
        UInt_t FirsteventIdx;
        UInt_t LasteventIdx;
        Bool_t IsFirstevent;
        Bool_t IsDecayEvt;
        UInt_t coboIC;
        UInt_t asadIC;
        UInt_t chanIC;
        UInt_t valueIC_min;
        Int_t L1Aflag;
        UInt_t IsTrig[16];
        int evtcounter;
        int goodevtcounter;
        int gatedevtcounter;
        int badevtcounter;
        int framecounter;
        int printed;
        int hitcounter;
        int prevmidx;
        int mutantcounter;
        vector<UInt_t> mevtidx; // eventidx value
        vector<UInt_t> mtstmp; // timestamp value
        vector<UInt_t> md2ptime; // d2ptime value
        double scaler1;
        double scaler1start;
        double scaler1end;
        double scaler2;
        double scaler2start;
        double scaler2end;
        double scaler3;
        double scaler3start;
        double scaler3end;
        double scaler4;
        double scaler4start;
        double scaler4end;
        double scaler5;
        double scaler5start;
        double scaler5end;
        ofstream scalerout;
        uint32_t scevent;
        uint32_t scevtidx;
        uint64_t sctstmp;
        double driftv;
        double rx;
        double ry;
        double rz;
        double rt;

        Int_t maxValue;
        Int_t maxValueBucket;
        Int_t maxValuedec;
        Int_t maxValueBucketdec;
        Double_t maxValuefit;
        Double_t maxValueBucketfit;
        UInt_t mm_minenergy;
        UInt_t mm_maxenergy;
        UInt_t mm_mintime;
        UInt_t mm_maxtime;
        UInt_t si_minenergy;
        UInt_t si_maxenergy;
        UInt_t si_mintime;
        UInt_t si_maxtime;
        TFile* fInputFile;
        TTree* fInputTree;
        TFile* fOutputFile;
        UInt_t fInputFileSize;
        TTree* fOutputTree;
        Int_t fNumberEvents;
        Float_t fTimePerBin;
        Int_t wGETMul;
        Int_t wGETHit;
        Int_t wGETEventIdx;
        Int_t wGETTimeStamp;
        Int_t wGETD2PTime;
        Int_t wGETFrameNo[4352];
        Int_t wGETDecayNo[4352];
        Int_t wGETL1Aflag[4352];
        Float_t wGETTime[4352];
        Float_t wGETEnergy[4352];
        Int_t wGETCobo[4352];
        Int_t wGETAsad[4352];
        Int_t wGETAget[4352];
        Int_t wGETChan[4352];
        Int_t wGETWaveformX[4352][512];
        Int_t wGETWaveformY[4352][512];
        Int_t wfmaxvalue[4352];
        Int_t wfminvalue[4352];
        Int_t wfbaseline[4352];
        Int_t wfdvalue[4352];
        Int_t rGETMul;
        Int_t rGETHit;
        Int_t rGETEventIdx;
        Int_t rGETTimeStamp;
        Int_t rGETD2PTime;
        Int_t rGETFrameNo[4352];
        Int_t rGETDecayNo[4352];
        Int_t rGETL1Aflag[4352];
        Float_t rGETTime[4352];
        Float_t rGETEnergy[4352];
        Int_t rGETCobo[4352];
        Int_t rGETAsad[4352];
        Int_t rGETAget[4352];
        Int_t rGETChan[4352];
        Int_t rGETWaveformX[4352][512];
        Int_t rGETWaveformY[4352][512];
        Double_t posenergysum[170];
        Double_t energysum[170];
        Double_t stripenergysum[170];
        Double_t chainenergysum[170][3];
        Int_t sumcounter[170];
        Int_t hMM_TrackCounter1[3];
        Int_t hMM_TrackCounter2[3];
        Int_t hMM_TrackdECounter[3];
        Double_t hMM_dE1[3];
        Double_t hMM_dE2[3];
        Double_t hMM_E[3];
        Int_t hMM_cornerfound[3];
        Int_t hMM_MinPxX[3];
        Int_t hMM_MinPxY[3];
        Int_t hMM_MaxPxX[3];
        Int_t hMM_MaxPxY[3];
        Double_t hMM_MinX[3];
        Double_t hMM_MinY[3];
        Double_t hMM_MaxX[3];
        Double_t hMM_MaxY[3];
        Double_t hMM_Slope[3];
        Double_t hMM_Offset[3];
        //int csi_channo[9]={2,7,10,16,19,25,28,33,36};
        int csi_channo[9]={2,7,16,10,19,25,33,28,36};
        int csi_map[64];
        Double_t trackposymin[3];
        Double_t trackposymax[3];
        int maxrespsample = 0;
        Int_t rftype;
        TF1* fResponse[10];
        TH1D* hResponse[10];
        Double_t response[10][512];
        Int_t responsesample[10][7];
        UInt_t maxresponse[10];
        UInt_t responsesigma[10];
        UInt_t deconvrep[10];
        UInt_t deconviter[10];
        UInt_t deconvboost[10];
        UInt_t decoffset;
        Double_t corrwaveformdec[512];
        Double_t corrwaveformfit[512];
};

#endif
