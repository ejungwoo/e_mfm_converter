#include "LKFrameBuilder.h"
#include "GSpectra.h"
#include "GNetServerRoot.h"
#include "WaveletNew.h"
#include "mfm/BitField.h"
#include "mfm/Field.h"
#include "mfm/Frame.h"
#include "mfm/FrameBuilder.h"
#include "mfm/FrameDictionary.h"
#include "mfm/Item.h"
#include <sstream>
#include <cstdio>
#include <boost/utility/binary.hpp>
const int no_cobos=3;
using namespace std;

#include "GETChannel.hpp"

WaveForms::WaveForms() {
}

WaveForms::~WaveForms() {
}

MM_Track::MM_Track() {
}

MM_Track::~MM_Track() {
}

Si_Track::Si_Track() {
}

Si_Track::~Si_Track() {
}

MapChanToMM::MapChanToMM() {
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            for(int k=0;k<64;k++){
                pxidx[i][j][k]=-1;
                pxidy[i][j][k]=-1;
            }
        }
    }
}

MapChanToMM::~MapChanToMM() {
}

MapChanToSi::MapChanToSi() {
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            for(int k=0;k<64;k++){
                pxidx[i][j][k]=-1;
                pxidy[i][j][k]=-1;
            }
        }
    }
}

MapChanToSi::~MapChanToSi() {
}

MapChanToX6::MapChanToX6() {
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            for(int k=0;k<68;k++){
                X6flag[i][j][k]=-1;
                X6det[i][j][k]=-1;
                X6strip[i][j][k]=-1;

                CsICTnum[i][j][k]=-1;
                CsI_X6det[i][j][k]=-1;
                CsI_X6ud[i][j][k]=-1;
            }
        }
    }
    for(int i=0; i<300; i++) {
        for(int j=0; j<8; j++) {
            X6posx[i][j]=-1;
            X6posy[i][j]=-1;
            X6posz[i][j]=-1;
        }
    }
}

MapChanToX6::~MapChanToX6() {
}

LKFrameBuilder::LKFrameBuilder(int port) {
    spectra_ = new GSpectra();
    serv_ = new GNetServerRoot(port,spectra_);

    serv_->StartServer();
}

LKFrameBuilder::~LKFrameBuilder() {
    delete serv_;
    delete spectra_;
}

void LKFrameBuilder::processFrame(mfm::Frame &frame)
{
    lk_debug << "[processFrame]" << endl;
    //XXX
    //cout << "isBlobFrame=" << frame.header().isBlobFrame() << ", frameType=" << frame.header().frameType() << endl;
    if (frame.header().isBlobFrame())
    {
        if (frame.header().frameType() == 0x7)
        {
            decodeCoBoTopologyFrame(frame);
        }
        else if (frame.header().frameType()  == 0x8) {
            decodeMuTanTFrame(frame);
        }
    }
    else
    {
        ValidateEvent(frame);
        goodsievt=1;
        goodicevt=1;

        if(goodsievt==1)
        {
            goodmmevt=1;
            //XXX
            Event(frame);
        }
    }
}

void LKFrameBuilder::Init(int mode, int d2pmode) {
    //cout<<"INIT HIST SERVER"<<endl;
    framecounter = 0;
    enableroot = 0;
    forcereadtree = 0;
    enablecleantrack = 0;
    energymethod = 0;
    readrw = 0;
    enable2pmode = d2pmode;
    enablehist = 0;
    enableupdatefast = 0;
    maxasad = 4;
    FirsteventIdx = 0;
    IsFirstevent = true;
    LasteventIdx = 0;
    weventIdx = 0;
    weventTime = 0;
    coboIC=2;
    asadIC=1;
    chanIC=2;
    valueIC_min=3500;
    decoffset=500; //ofset value for Deconvolution Method
    mm_minenergy = 0;
    //mm_maxenergy = 4096;
    //mm_minenergy = 0; //for AllChannelReading
    //mm_minenergy = 0;
    mm_maxenergy = 4096;
    //mm_mintime = 130;
    //mm_maxtime = 380;
    mm_mintime = 0;
    mm_maxtime = 350;
    si_minenergy = 0;
    si_maxenergy = 5000;
    si_mintime = 0;
    si_maxtime = 512;
    //driftv = 0.0025; // mm/ns, for P5 at 100 Torr and -800V
    driftv = 0.017; // mm/ns, for CO2 at 20 Torr and -800V
    fTimePerBin = 1000.; // 100 ns per bin

    mapchantomm = new MapChanToMM();
    mapchantosi = new MapChanToSi();
    mapchantox6 = new MapChanToX6();
    for(int i=0; i<maxasad ; i++) {
        for(int j=0; j<4 ; j++) {
            for(int k=0; k<64 ; k++) {
                mapchantomm->pxidx[i][j][k]=0;
                mapchantomm->pxidy[i][j][k]=0;
                mapchantosi->pxidx[i][j][k]=0;
                mapchantosi->pxidy[i][j][k]=0;
            }
        }
    }

    cout<<"MODE: "<<mode<<endl;

    if(mode==1){
        mfm::FrameDictionary::instance().addFormats("/usr/local/get/share/get-bench/format/CoboFormats.xcfg");
        waveforms = new WaveForms();
        InitWaveforms();
    }
}

void LKFrameBuilder::InitWaveforms() {
    waveforms->waveform.resize(maxasad*4);
    waveforms->hasSignal.resize(maxasad*4);
    waveforms->hasHit.resize(maxasad*4); //default set to false
    waveforms->hasFPN.resize(maxasad*4); //default set to false
    waveforms->doneFPN.resize(maxasad*4); //default set to false
    for(int i=0;i<maxasad;i++){
        for(int j=0;j<4;j++){
            waveforms->waveform[i*4+j].resize(68);
            waveforms->hasSignal[i*4+j].resize(68); //default set to false
            for(int k=0;k<68;k++){
                waveforms->waveform[i*4+j][k].resize(bucketmax);
            }
        }
    }
}

void LKFrameBuilder::ResetWaveforms() {
    for(int i=0;i<maxasad;i++){
        for(int j=0;j<4;j++){
            if(waveforms->hasHit[i*4+j] || waveforms->hasFPN[i*4+j]){
                for(int k=0;k<68;k++){
                    if(waveforms->hasSignal[i*4+j][k]){
                        waveforms->hasSignal[i*4+j][k] = false;
                        fill(waveforms->waveform[i*4+j][k].begin(),
                                waveforms->waveform[i*4+j][k].end(),0);
                    }
                }
            }
            waveforms->hasHit[i*4+j] = false;
            waveforms->hasFPN[i*4+j] = false;
            waveforms->doneFPN[i*4+j] = false;

        }
    }
    waveforms->ICenergy = 0;
    waveforms->Sienergy = 0;
    waveforms->SiX = 0;
    waveforms->SiY = 0;
    waveforms->SiZ = 0;
    waveforms->isRejected = 0;
}

void LKFrameBuilder::ValidateEvent(mfm::Frame& frame)
{
    lk_debug << "[ValidateEvent]" << endl;
    if(frame.header().isLayeredFrame()) {
        unique_ptr<mfm::Frame> subFrame;
        for(int i = 0;i<frame.itemCount();i++) {
            try{
                try{
                    subFrame = frame.frameAt(i);
                }catch (const std::exception& e){
                    lk_debug << "error" << endl;
                    cout << e.what() << endl;
                    return;
                }
                if((*subFrame.get()).itemCount()>0){ //Make sure we have data
                    ValidateFrame(*subFrame.get());
                }else{
                    //cout << "1no subframe" << endl;
                }
            }catch (const std::exception& e){
                cout << e.what() << endl;
                return;
            }
        }
    } else {
        if(frame.itemCount()>0){ //Make sure we have data
            cout << "Not Layered Frame" << endl;
            ValidateFrame(frame);
        }
    }
}

void LKFrameBuilder::ValidateFrame(mfm::Frame& frame)
{
    UInt_t coboIdx = frame.headerField("coboIdx").value<UInt_t>();
    UInt_t asadIdx = frame.headerField("asadIdx").value<UInt_t>();
    UInt_t ceventIdx = (Int_t)frame.headerField("eventIdx").value<UInt_t>();

    //unsigned short SiMask = BOOST_BINARY( 111011111111111011111111110111111111111111111111101111111111011111111111 );
    ULong_t SiMask = 0xFEFFDFFFFFBFF7FF; //LSB 8byte
    //unsigned short MMMask = BOOST_BINARY( 111011111111111011111111110111111111111111111111101111111111011111111111 );
    ULong_t MMMask = 0xFEFFDFFFFFBFF7FF; //LSB 8byte
    //unsigned short MMMask1 = BOOST_BINARY( 000000000000000000000000000000000000000111111111101111111111011111111111 );
    ULong_t MMMask1 = 0x3FFBFF7FF; //LSB 8byte
    //unsigned short IcMask = BOOST_BINARY( 000000100000000000000000000000000000000000000000000000000000000000000000 );
    ULong_t IcMask = 0x13; // Channel 19

    //cout << coboIdx << " " << asadIdx << "(LSB):" << endl;
    ULong_t hitPat_0 = frame.headerField(32u, 8u).value<ULong_t>();
    ULong_t hitPat_1 = frame.headerField(41u, 8u).value<ULong_t>();
    ULong_t hitPat_2 = frame.headerField(50u, 8u).value<ULong_t>();
    ULong_t hitPat_3 = frame.headerField(59u, 8u).value<ULong_t>();
    //cout << "HitPattern 0: " << hex << hitPat_0 << endl;
    //cout << "HitPattern 1: " << hex << hitPat_1 << endl;
    //cout << "HitPattern 2: " << hex << hitPat_2 << endl;
    //cout << "HitPattern 3: " << hex << hitPat_3 << endl;
    //cout << "===========================" << endl;
    //000100000000000100000000001000000000000000000000010000000000100000000000
    //Si and IC data validation
    if(coboIdx==1 && asadIdx==0 && ((hitPat_0&SiMask)||(hitPat_1&SiMask)||(hitPat_2&SiMask)||(hitPat_3&SiMask))){
        goodsievt=1;
        goodevtidx = ceventIdx;
    }

    //MM and IC data validation
    //if(coboIdx==0 && ((hitPat_0&MMMask)||(hitPat_1&MMMask)||(hitPat_2&MMMask)||(hitPat_3&MMMask)))
    if(coboIdx==0 && ((asadIdx==2||asadIdx==3) && hitPat_2&MMMask1) && (asadIdx==0 && (hitPat_0&MMMask)) && (asadIdx<2 && ((hitPat_1&MMMask)||(hitPat_2&MMMask)||(hitPat_3&MMMask))))
    {
        goodmmevt=1;
        goodevtidx = ceventIdx;
    }

    //cout << coboIdx << " " << asadIdx << "(USB):" << endl;
    hitPat_0 = frame.headerField(31u, 1u).value<ULong_t>();
    hitPat_1 = frame.headerField(40u, 1u).value<ULong_t>();
    hitPat_2 = frame.headerField(49u, 1u).value<ULong_t>();
    hitPat_3 = frame.headerField(58u, 1u).value<ULong_t>();
    //cout << "HitPattern 0: " << hex << hitPat_0 << endl;
    //cout << "HitPattern 1: " << hex << hitPat_1 << endl;
    //cout << "HitPattern 2: " << hex << hitPat_2 << endl;
    //cout << "HitPattern 3: " << hex << hitPat_3 << endl;
    //cout << "===========================" << endl;
    //000100000000000100000000001000000000000000000000010000000000100000000000

    //Si and IC data validation
    //unsigned short SiMask = BOOST_BINARY( 111011111111111011111111110111111111111111111111101111111111011111111111 );
    SiMask = 0xEF; //USB 1byte
    if(coboIdx==1 && asadIdx==0 && ((hitPat_0&SiMask)||(hitPat_1&SiMask)||(hitPat_2&SiMask)||(hitPat_3&SiMask))){
        goodsievt=1;
        goodevtidx = ceventIdx;
    }

    //MM and IC data validation
    MMMask = 0xEF; //USB 1byte
    //if(coboIdx==0 && ((hitPat_0&MMMask)||(hitPat_1&MMMask)||(hitPat_2&MMMask)||(hitPat_3&MMMask)))
    if(coboIdx==0 && (asadIdx==0 && (hitPat_0&MMMask)) && (asadIdx<2 && ((hitPat_1&MMMask)||(hitPat_2&MMMask)||(hitPat_3&MMMask))))
    {
        goodmmevt=1;
        goodevtidx = ceventIdx;
    }
    if(coboIdx==1 && asadIdx==1 && (hitPat_0&IcMask)){
        goodicevt=1;
    }
}

void LKFrameBuilder::Event(mfm::Frame& frame)
{
    //XXX
    waveforms->frameIdx = 0;
    waveforms->decayIdx = 0;
    if(frame.header().isLayeredFrame()) {
        for(int i = 0;i<frame.itemCount();i++) {
            //cout << "1isLayered=" << frame.header().isLayeredFrame() << ", itemCount=" << frame.itemCount() << ", frameIndex=" << i << endl;
            waveforms->frameIdx = 0;
            waveforms->decayIdx = 0;
            try{
                auto_ptr<mfm::Frame> subFrame = frame.frameAt(i);
                if((*subFrame.get()).itemCount()>0){ //Make sure we have data
                    //cout << "2isLayered=" << frame.header().isLayeredFrame() << ", itemCount=" << frame.itemCount() << ", frameIndex=" << i << endl;
                    waveforms->frameIdx = i;
                    UnpackFrame(*subFrame.get());
                    RootWConvert();
                    ResetWaveforms();
                }else{
                    //cout << "2no subframe" << endl;
                }
            }catch (const std::exception& e){
                cout << e.what() << endl;
                return;
            }
        }
    } else {
        try{
            if(frame.itemCount()>0){ //Make sure we have data
                //cout << "Not Layered Frame" << endl;
                UnpackFrame(frame);
                RootWConvert();
                ResetWaveforms();
            }else{
                cout << "3no subframe" << endl;
            }
        }catch (const std::exception& e){
            cout << e.what() << endl;
            return;
        }
    }
}

void LKFrameBuilder::UnpackFrame(mfm::Frame& frame)
{
    //XXX
    UInt_t coboIdx = frame.headerField("coboIdx").value<UInt_t>();
    UInt_t asadIdx = frame.headerField("asadIdx").value<UInt_t>();
    Int_t prevweventIdx = weventIdx;
    UInt_t frameSize = frame.headerField("frameSize").value<UInt_t>();
    UInt_t itemSize = frame.headerField("itemSize").value<UInt_t>();
    weventIdx = (Int_t)frame.headerField("eventIdx").value<UInt_t>();
    UInt_t prevweventTime = weventTime;
    weventTime = (UInt_t)frame.headerField("eventTime").value<UInt_t>();
    //cout<<"weventIdx=" << weventIdx << ", prevweventIdx+1=" << (prevweventIdx+1) << endl;
    //if(weventIdx>prevweventIdx+1) return;
    if(IsFirstevent){
        cout<<"FIRST EVENT!!!!!"<<endl;
        FirsteventIdx = weventIdx;
        IsFirstevent=false;
    }

    //if(enableskipevent==1 && weventIdx<firsteventno) return;
    //if(enableskipevent==2 && weventIdx<maxevtno){
    //if(!evtmask[weventIdx]) return;
    //}else if(enableskipevent==2 && weventIdx>=maxevtno){
    //  //cout << Form("weventIdx >= %d!",maxevtno) << endl;
    //}
    //else cout << "weventIdx " << weventIdx << endl;

    if((prevweventIdx != weventIdx))
    {
        if(readmode==1)
        {
            if(wGETMul>0){
                wGETEventIdx = prevweventIdx;
                if(prevgoodmmevt==1 || 1){
                    cout << "In data: " << wGETEventIdx << "enable2pmode: "<<enable2pmode << endl;
                    if(enable2pmode==1){
                        //cout<<"mevtidx.size() = "<<mevtidx.size()<<endl;
                        for(int i=prevmidx;i<mevtidx.size();i++){
                            //cout<<mevtidx.at(i)<<"\t"<<wGETEventIdx<<endl;
                            if(mevtidx.at(i)==wGETEventIdx){
                                wGETD2PTime = md2ptime.at(i);
                                wGETTimeStamp = mtstmp.at(i);
                                //cout << "In data: " << i << " " << wGETEventIdx << " " << wGETD2PTime << " " << wGETTimeStamp << endl;
                                //cout << "In data: " << i << " " << wGETEventIdx << " " << wGETD2PTime << " " << wGETTimeStamp << " " << prevweventTime << endl;
                                prevmidx = i;
                            }
                        }
                    }
                    //XXX
                    RootWriteEvent();
                    RootWReset();
                    prevgoodmmevt = goodmmevt;
                    prevgoodicevt = goodicevt;
                    goodsievt=0;
                    goodmmevt=0;
                    goodicevt=0;
                    goodevtidx=0;
                }else{
                    wGETEventIdx = prevweventIdx;
                    RootWReset();
                    prevgoodmmevt = goodmmevt;
                    prevgoodicevt = goodicevt;
                    goodsievt=0;
                    goodmmevt=0;
                    goodicevt=0;
                    goodevtidx=0;
                }
                if((wGETEventIdx-FirsteventIdx)%50==0){
                    if(wGETEventIdx==FirsteventIdx){
                        cout << Form("M%d:Starting from Event No. %d.",readmode,wGETEventIdx) << endl;
                    }else{
                        cout << Form("M%d:Upto Event No. %d Converted..(from Event No. %d)",readmode,wGETEventIdx, FirsteventIdx) << endl;
                    }
                }
                wGETEventIdx = weventIdx;
                wGETMul=0;
                wGETHit=0;
                for(int i=0;i<16;i++) IsTrig[i]=0;
            }
        }
    }

    lx_info<<"Frame type: "<<frame.header().frameType()<<"\teventIdx"<<weventIdx<<endl;

    //  if(weventIdx<3350){
    //    return;
    //  }

    mfm::Item item = frame.itemAt(0u);
    mfm::Field field = item.field("");

    const size_t numSamples = frame.itemCount();

    const size_t numChannels = 68u;
    const size_t numChips = 4u;
    vector<uint32_t> chanIdx_(numChips, 0u);
    vector<uint32_t> buckIdx_(numChips, 0u);
    if(frame.header().frameType() == 1u) {
        mfm::Item item = frame.itemAt(0u);
        mfm::Field field = item.field("");
        mfm::BitField agetIdxField = field.bitField("agetIdx");
        mfm::BitField chanIdxField = field.bitField("chanIdx");
        mfm::BitField buckIdxField = field.bitField("buckIdx");
        mfm::BitField sampleValueField = field.bitField("sample");

        if(enable2pmode==1){
            if(IsTrig[coboIdx*4+asadIdx]>0){
                IsTrig[coboIdx*4+asadIdx]=2;
                waveforms->decayIdx = 1;
            }else{
                IsTrig[coboIdx*4+asadIdx]=1;
                waveforms->decayIdx = 0;
            }
        }else{
            IsTrig[coboIdx*4+asadIdx]=0;
            waveforms->decayIdx = 0;
        }
        //cout << "Type:" << frame.header().frameType() << " " <<  weventIdx << "  " << frameSize << " " << coboIdx << " " << asadIdx << " " << waveforms->frameIdx << " " << waveforms->decayIdx << " " << frame.itemCount() << endl;
        int lastaget=-1;
        for(UInt_t i=0; i<frame.itemCount(); i++) {
            item = frame.itemAt(i);
            field = item.field(field);
            agetIdxField = field.bitField(agetIdxField);
            chanIdxField = field.bitField(chanIdxField);
            buckIdxField = field.bitField(buckIdxField);
            sampleValueField =field.bitField(sampleValueField);

            UInt_t agetIdx = agetIdxField.value<UInt_t>();
            //if(agetIdx==lastaget) break;
            lastaget=agetIdx;
            if(agetIdx<0 || agetIdx>3) return;
            UInt_t chanIdx = chanIdxField.value<UInt_t>();
            UInt_t buckIdx = buckIdxField.value<UInt_t>();
            UInt_t sampleValue = sampleValueField.value<UInt_t>();

            waveforms->coboIdx = coboIdx;
            waveforms->asadIdx = asadIdx;

            if((ignoremm==0) || (ignoremm==1 && coboIdx>0)){ // skip MM waveform data
                waveforms->hasSignal[asadIdx*4+agetIdx][chanIdx] = true;
                //      cout<<"BEEP "<<coboIdx<<"\t"<<asadIdx<<"\t"<<agetIdx<<"\t"<<chanIdx<<endl;
                if((chanIdx==11||chanIdx==22||chanIdx==45||chanIdx==56)) waveforms->hasFPN[asadIdx*4+agetIdx] = true;
                else waveforms->hasHit[asadIdx*4+agetIdx] = true;
                waveforms->waveform[asadIdx*4+agetIdx][chanIdx][buckIdx] = sampleValue;
                //cout << "Type:" << frame.header().frameType() << " " <<  weventIdx << "  " << frameSize << " " << coboIdx << " " << asadIdx << " " << agetIdx << " " << chanIdx << " " << buckIdx <<  " " << sampleValue << endl;
            }
        }

    }
    else if (frame.header().frameType() == 2u) {
        if(enable2pmode==1){
            if(IsTrig[coboIdx*4+asadIdx]>0){
                IsTrig[coboIdx*4+asadIdx]=2;
                waveforms->decayIdx = 1;
            }else{
                IsTrig[coboIdx*4+asadIdx]=1;
                waveforms->decayIdx = 0;
            }
        }else{
            IsTrig[coboIdx*4+asadIdx]=0;
            waveforms->decayIdx = 0;
        }
        //cout << "Type:" << frame.header().frameType() << " " <<  weventIdx << " " << frameSize << " " << coboIdx << " " << asadIdx << " " << waveforms->frameIdx << " " << waveforms->decayIdx << " " << endl;
        // Decode first item
        if(numSamples>0) {
            mfm::Item item = frame.itemAt(0u);
            mfm::Field field = item.field("");
            mfm::BitField agetIdxField = field.bitField("agetIdx");
            mfm::BitField sampleValueField = field.bitField("sample");

            const uint32_t agetIdx = agetIdxField.value<uint32_t>();
            const uint32_t sampleValue = sampleValueField.value<uint32_t>();

            waveforms->coboIdx = coboIdx;
            waveforms->asadIdx = asadIdx;

            if((ignoremm==0) || (ignoremm==1 && coboIdx>0)){ // skip MM waveform data
                waveforms->hasSignal[asadIdx*4+agetIdx][chanIdx_[agetIdx]] = true;
                if((chanIdx_[agetIdx]==11||chanIdx_[agetIdx]==22||chanIdx_[agetIdx]==45||chanIdx_[agetIdx]==56)) waveforms->hasFPN[asadIdx*4+agetIdx] = true;
                else waveforms->hasHit[asadIdx*4+agetIdx] = true;
                waveforms->waveform[asadIdx*4+agetIdx][chanIdx_[agetIdx]][buckIdx_[agetIdx]] = sampleValue;
            }
        }
        // Loop through other items
        for(size_t itemId=1; itemId<numSamples; ++itemId) {
            item = frame.itemAt(itemId);
            field = item.field(field);
            mfm::BitField agetIdxField = field.bitField(agetIdxField);
            mfm::BitField sampleValueField = field.bitField(sampleValueField);

            if((ignoremm==0) || (ignoremm==1 && coboIdx>0)){ // skip MM waveform data
                const uint32_t agetIdx = agetIdxField.value<uint32_t>();
                if(chanIdx_[agetIdx]>=numChannels) {
                    chanIdx_[agetIdx] = 0u;
                    buckIdx_[agetIdx]++;
                }
                const uint32_t sampleValue = sampleValueField.value<uint32_t>();

                //if(itemId>57120) cout << itemId << " " << numSamples << " " << coboIdx << " " << asadIdx << " " << agetIdx << " " << chanIdx_[agetIdx] << " " << buckIdx_[agetIdx] << " " << sampleValue << endl;

                waveforms->hasSignal[asadIdx*4+agetIdx][chanIdx_[agetIdx]] = true;
                if((chanIdx_[agetIdx]==11||chanIdx_[agetIdx]==22||chanIdx_[agetIdx]==45||chanIdx_[agetIdx]==56)) waveforms->hasFPN[asadIdx*4+agetIdx] = true;
                else waveforms->hasHit[asadIdx*4+agetIdx] = true;
                waveforms->waveform[asadIdx*4+agetIdx][chanIdx_[agetIdx]][buckIdx_[agetIdx]] = sampleValue;
                chanIdx_[agetIdx]++;
            }
        }
    }
    else {
        cout << "Frame type " << frame.header().frameType() << " not found" << endl;
        return;
    }
}

void LKFrameBuilder::RootFindEvent()
{
    UInt_t reventIdx1 = 0;
    UInt_t reventIdx2 = 0;
    Int_t rgetmul1 = 0;
    Int_t rgetmul2 = 0;
    Int_t coboIdx1 = 0;
    Int_t asadIdx1 = 0;
    Int_t agetIdx1 = 0;
    Int_t chanIdx1 = 0;
    Int_t coboIdx2 = 0;
    Int_t asadIdx2 = 0;
    Int_t agetIdx2 = 0;
    Int_t chanIdx2 = 0;
    UInt_t oldbucketmax = bucketmax;
    RootRInit();
    for(int i=0;i<fNumberEvents;i++){
        if(readmode>0){
            fInputFile->cd();
            fInputTree->GetEntry(i);
        }
        reventIdx1 = rGETEventIdx;
        rgetmul1 = rGETMul;
        if(reventIdx1>0){
            fInputTree->GetEntry(i);
            coboIdx1 = rGETCobo[0];
            asadIdx1 = rGETAsad[0];
            agetIdx1 = rGETAget[0];
            chanIdx1 = rGETChan[0];
            //cout << "Idx1:" << reventIdx1 << " " << rgetmul1 << " " << coboIdx1 << " " << asadIdx1 << " " << agetIdx1 << " " << chanIdx1 << endl;
            for(int k=(i+1);k<fNumberEvents;k++){
                if(readmode>0){
                    fInputFile->cd();
                    fInputTree->GetEntry(k);
                }
                reventIdx2 = rGETEventIdx;
                if(reventIdx1 == reventIdx2){
                    for(Int_t l=0; l<rGETMul; l++){
                        coboIdx2 = rGETCobo[l];
                        asadIdx2 = rGETAsad[l];
                        agetIdx2 = rGETAget[l];
                        chanIdx2 = rGETChan[l];
                        if(chanIdx2==11 || chanIdx2==22 || chanIdx2==45 || chanIdx2==56){
                            continue;
                        }else{
                            //            cout << "Idx2:" << reventIdx2 << " " << rGETMul << " " << coboIdx2 << " " << asadIdx2 << " " << agetIdx2 << " " << chanIdx2 << endl;
                        }
                    }
                }
            }
        }
    }
}

void LKFrameBuilder::RootReadEvent()
{
    Int_t frameIdx = 0;
    Int_t decayIdx = 0;
    L1Aflag = 0;
    Int_t coboIdx = 0;
    Int_t asadIdx = 0;
    Int_t agetIdx = 0;
    Int_t chanIdx = 0;
    UInt_t percent = 0;
    UInt_t Entries = 1;
    UInt_t oldbucketmax = bucketmax;
    if(enable2pmode==1){
        bucketmax = bucketmax/2;
        cout << "bucketmax set to " << bucketmax << endl;
    }
    ostringstream s4out;
    ostringstream s5out;
    ofstream f4out;
    ofstream f5out;
    RootRInit();
    if(readmode>0){
        if(fNumberEvents>0){
            Entries = fNumberEvents;
            percent = Entries/10;
        }
        SetResponseWaveform();
        GetMaxResponseWaveform();
        GetSigmaResponseWaveform();
    }
    for(int i=0;i<fNumberEvents;i++){
        if(readmode>0){
            fInputFile->cd();
            fInputTree->GetEntry(i);
            if(RootRRecovered()){
                if(LasteventIdx>=reventIdx) continue;
                else LasteventIdx = reventIdx;
            }else{
                LasteventIdx = reventIdx;
            }
        }
        reventIdx = rGETEventIdx;
        goodsicsievt=0;
        goodsicsipevt=0;
        if(IsFirstevent){
            //cout<<"SET FIRST EVENT"<<endl;
            FirsteventIdx = reventIdx;
            evtcounter=0;
            goodevtcounter=0;
            gatedevtcounter=0;
            badevtcounter=0;
            IsFirstevent=false;
            IsDecayEvt=false;
            RootRResetWaveforms();
        }else{
            if(enable2pmode==0){
                evtcounter++;
            }else if(enable2pmode==1 && IsDecayEvt==true){
                evtcounter++;
            }
            IsDecayEvt=false;
            RootRResetWaveforms();
        }

        if(enableskipevent==1 && reventIdx<firsteventno) continue;
        if(enableskipevent==2 && reventIdx<maxevtno){
            if(!evtmask[reventIdx]) continue;
        }else if(enableskipevent==2 && reventIdx>=maxevtno){
            cout << Form("reventIdx >= %d!",maxevtno) << endl;
        }

        //if((reventIdx-FirsteventIdx)%1000==0)
        if((reventIdx-FirsteventIdx)%50==0)
        {
            if(reventIdx==FirsteventIdx){
                //cout << Form("M%d:Starting from Event No. %d.",readmode,reventIdx) << endl;
            }else{
                //cout << Form("M%d:Upto Event No. %d Processed (from %d)..",readmode,reventIdx, FirsteventIdx) << endl;
            }
        }

        //cout << rGETEventIdx << " " << rGETD2PTime << " " << rGETMul << endl;
        rd2ptime = rGETD2PTime;
        //cout<<rd2ptime<<endl;
        rtstmp = rGETTimeStamp;
        if(rd2ptime>30000000) continue;
        int goodX6counter=0;
        int goodFWfcounter=0;
        int goodFWbcounter=0;

        int L0time=0;
        if(enable2pmode==1){
            if(rd2ptime>0) hMM_D2PTime->Fill(rd2ptime);
            //else cout << reventIdx << " " << rd2ptime << " " << rGETMul << endl;

            printed=0;

            for(Int_t i=0; i<rGETMul; i++){
                frameIdx = rGETFrameNo[i];
                decayIdx = rGETDecayNo[i];
                coboIdx = rGETCobo[i];
                asadIdx = rGETAsad[i];
                agetIdx = rGETAget[i];
                chanIdx = rGETChan[i];

                if(ignoremm==1 && coboIdx==0) continue;// skip MM waveform data
                if(decayIdx==1) IsDecayEvt=true;
                if(coboIdx==coboIC && asadIdx==asadIC && chanIdx==chanIC){
                    for(int j=0;j<bucketmax;j++){
                        if(rGETWaveformY[i][j] < valueIC_min){
                            //            L1Aflag = frameIdx%2;
                        }
                    }
                    //cout << reventIdx << " " << coboIdx << " " << asadIdx << " " << chanIdx << " " << frameIdx << " " << L1Aflag << " " << decayIdx << endl;
                }else if(coboIdx==0){
                    for(int j=0;j<bucketmax;j++){
                        if(printed==0 && rGETWaveformY[i][j]>4000){
                            //cout << "MM2: " << reventIdx << " " << coboIdx << " " << asadIdx << " " << chanIdx << " " << frameIdx << " " << decayIdx << endl;
                            printed++;
                        }
                    }
                }
            }
        }

        for(Int_t i=0; i<rGETMul; i++){
            decayIdx = rGETDecayNo[i];
            coboIdx = rGETCobo[i];
            asadIdx = rGETAsad[i];
            agetIdx = rGETAget[i];
            chanIdx = rGETChan[i];

            if(coboIdx==2 && agetIdx!=3) goodX6counter++;
            if(coboIdx==1 && asadIdx==0 && agetIdx==0) goodFWfcounter++;
            if(coboIdx==1 && asadIdx==0 && agetIdx==1) goodFWbcounter++;
            if(coboIdx==1 && asadIdx==1 && agetIdx==3 && chanIdx==19) {
                int maxval=0;
                for(int t=0;t<100;t++) {
                    if(rGETWaveformY[i][t]>maxval) {
                        maxval=rGETWaveformY[i][t];
                        L0time=t;
                    }
                }
            }

            if(ignoremm==1 && coboIdx==0) continue;// skip MM waveform data
            wfmaxvalue[i]=-10000;
            wfminvalue[i]=10000;
            wfbaseline[i] = 0;
            Int_t baselinecounter = 0;
            wfdvalue[i] = 0;
            for(int j=0;j<bucketmax;j++){
                if(rGETWaveformY[i][j]>0 && rGETWaveformY[i][j]<4095){
                    if(wfmaxvalue[i]<rGETWaveformY[i][j]) wfmaxvalue[i] = rGETWaveformY[i][j];
                    if(wfminvalue[i]>rGETWaveformY[i][j]) wfminvalue[i] = rGETWaveformY[i][j];
                    if(j>=20+decayIdx*256 && j<30+decayIdx*256){
                        wfbaseline[i] += rGETWaveformY[i][j];
                        baselinecounter++;
                    }
                }
            }
            if(baselinecounter>0) wfbaseline[i] /= baselinecounter;
            wfdvalue[i] = wfmaxvalue[i] - wfminvalue[i];
            //if(coboIdx==0 && (wfdvalue[i]>(wfmaxvalue[i]-wfbaseline[i]+200))) wfdvalue[i] = -1;
            //if(coboIdx==1 && asadIdx==0 && (agetIdx==1||agetIdx==3) && (wfdvalue[i]>(wfmaxvalue[i]-wfbaseline[i]+100))) wfdvalue[i] = -1;
            //cout << reventIdx << ", asadIdx=" << asadIdx << ", wfdvalue[" << i << "]= " << wfdvalue[i] << ",  wfbaseline[" << i << "]= " << wfbaseline[i] << " " << (wfmaxvalue[i]-wfbaseline[i]+300) << endl;
        }
        goodsievt=0;
        //if((goodX6counter!=3 && (goodFWbcounter+goodFWfcounter)!=2)||(goodFWbcounter!=1 && goodFWfcounter!=1)) {evtcounter--; continue;} // for X6 and Fwd Si. condition
        //if(goodX6counter!=3 || (L0time<29 || L0time>33)) {evtcounter--; return;}
        //if(goodX6counter<3) {evtcounter--; return;}
        goodsievt=1;
        //std::cout<<"Good event no: "<<reventIdx<< ", Mult=" << rGETMul << std::endl;
        for(Int_t i=0; i<rGETMul; i++){
            //if(wfdvalue[i]>40)
            if(wfdvalue[i]>0)
            {
                decayIdx = rGETDecayNo[i];
                frameIdx = rGETFrameNo[i];
                coboIdx = rGETCobo[i];
                asadIdx = rGETAsad[i];
                agetIdx = rGETAget[i];
                chanIdx = rGETChan[i];
                if(ignoremm==1 && coboIdx==0) continue;// skip MM waveform data
                if(enable2pmode==1){
                    //          if(coboIdx==coboIC && asadIdx==asadIC && chanIdx==chanIC){
                    //              decayIdx=0;
                    //          }else{
                    //            if(frameIdx%2 == L1Aflag){
                    //              decayIdx=0;
                    //            }else{
                    //              decayIdx=1;
                    //            }
                    //          }
                }else{
                    //          decayIdx=0;
                }
                rwaveforms[decayIdx][coboIdx]->decayIdx = decayIdx;
                rwaveforms[decayIdx][coboIdx]->hasHit[asadIdx*4+agetIdx] = true;
                if((chanIdx==11||chanIdx==22||chanIdx==45||chanIdx==56)){
                    rwaveforms[decayIdx][coboIdx]->hasFPN[asadIdx*4+agetIdx] = true;
                }
                //cout << reventIdx << " " << coboIdx << " " << asadIdx << " " << agetIdx << " " << chanIdx << " " << rwaveforms[decayIdx][coboIdx]->isDecay[asadIdx*4+agetIdx][chanIdx] << endl;
                if(enable2pmode==0){
                    for(int j=0;j<bucketmax;j++){
                        rwaveforms[decayIdx][coboIdx]->waveform[asadIdx*4+agetIdx][chanIdx][j] = rGETWaveformY[i][j];
                        //cout << coboIdx << " " << asadIdx << " " << agetIdx << " " << j << " " << rwaveforms[decayIdx][coboIdx]->waveform[asadIdx*4+agetIdx][chanIdx][j] << " " << rGETWaveformY[i][j] << endl;
                    }
                }else{
                    //if(decayIdx==1) cout << "Yo:" << reventIdx << " " << decayIdx << " " << coboIdx << " " << asadIdx << " " << agetIdx << " " << chanIdx << " " << rwaveforms[decayIdx][coboIdx]->isDecay[asadIdx*4+agetIdx][chanIdx] << endl;
                    for(int j=0;j<bucketmax;j++){
                        rwaveforms[decayIdx][coboIdx]->waveform[asadIdx*4+agetIdx][chanIdx][j] = rGETWaveformY[i][j];
                        //if(coboIdx==1 && rwaveforms[decayIdx][coboIdx]->waveform[asadIdx*4+agetIdx][chanIdx][j]>0)cout << j << " " << rwaveforms[decayIdx][coboIdx]->waveform[asadIdx*4+agetIdx][chanIdx][j] << " " << rGETWaveformY[i][j-decayIdx*256] << endl;
                    }
                }
                /*
                   if(coboIdx==1&&asadIdx==0){
                   s4out << reventIdx << " " << coboIdx << " " << asadIdx << " " << agetIdx << " " << chanIdx << " 4" << endl;
                   }else if(coboIdx==1&&asadIdx==1&&chanIdx==33){
                   s5out << reventIdx << " " << coboIdx << " " << asadIdx << " " << agetIdx << " " << chanIdx << " 5" << endl;
                   }
                 */
            }
        }

        /*
           if(s4out.str().size()!=0){
           f4out.open("silabels.txt", std::ofstream::out|std::ofstream::app);
           f4out << s4out.str() << endl;
           f4out.close();
           s4out.clear();
           s4out.str("");
           }
           if(s5out.str().size()!=0){
           f5out.open("iclabels.txt", std::ofstream::out|std::ofstream::app);
           f5out << s5out.str() << endl;
           f5out.close();
           s5out.clear();
           s5out.str("");
           }
         */
        //cout << "Done with reading events from root tree" << endl;

        Int_t lcwaveforms[6][512];
        Int_t lccounter[6];
        Int_t dchanIdx=0;
        Int_t Beam_med = 180; //low
        Int_t Beam_er = 150; //low
        Int_t Bi = Beam_med-Beam_er;
        Int_t Bf = Beam_med+Beam_er;
        Int_t Beam_window = 40; //event window
        Int_t Particle_window = 150; //event window
        Int_t btime[6];
        Int_t bestbtime;
        Int_t bmax[6];
        for(Int_t i=0;i<6;i++){
            for(Int_t buck=0; buck<512; buck++) lcwaveforms[i][buck]=0;
            lccounter[i]=0;
            bmax[i]=0;
        }
        for(Int_t i=0; i<rGETMul; i++){
            if(wfdvalue[i]>300){
                coboIdx = rGETCobo[i];
                if(ignoremm==1 && coboIdx==0) continue;// skip MM waveform data
                if(coboIdx==0){
                    asadIdx = rGETAsad[i];
                    agetIdx = rGETAget[i];
                    chanIdx = rGETChan[i];
                    if(chanIdx==11 || chanIdx==22 || chanIdx==45 || chanIdx==56) continue; // We want to skip the FPN channels.
                    if(chanIdx<11) {
                        dchanIdx = chanIdx;
                    }else if(chanIdx>11 && chanIdx<22) {
                        dchanIdx = chanIdx - 1;
                    }else if(chanIdx>22 && chanIdx<45) {
                        dchanIdx = chanIdx - 2;
                    }else if(chanIdx>45 && chanIdx<56) {
                        dchanIdx = chanIdx - 3;
                    }else if(chanIdx>56) {
                        dchanIdx = chanIdx - 4;
                    }
                    if(mapchantomm->pxidx[asadIdx][agetIdx][dchanIdx]==64){
                        for(Int_t buck=Bi; buck<Bf; buck++) lcwaveforms[0][buck]+=rwaveforms[decayIdx][coboIdx]->waveform[asadIdx*4+agetIdx][chanIdx][buck];
                        lccounter[0]++;
                    }else if(mapchantomm->pxidx[asadIdx][agetIdx][dchanIdx]==65){
                        for(Int_t buck=Bi; buck<Bf; buck++) lcwaveforms[1][buck]+=rwaveforms[decayIdx][coboIdx]->waveform[asadIdx*4+agetIdx][chanIdx][buck];
                        lccounter[1]++;
                    }else if(mapchantomm->pxidx[asadIdx][agetIdx][dchanIdx]==66){
                        for(Int_t buck=Bi; buck<Bf; buck++) lcwaveforms[2][buck]+=rwaveforms[decayIdx][coboIdx]->waveform[asadIdx*4+agetIdx][chanIdx][buck];
                        lccounter[2]++;
                    }else if(mapchantomm->pxidx[asadIdx][agetIdx][dchanIdx]==67){
                        for(Int_t buck=Bi; buck<Bf; buck++) lcwaveforms[3][buck]+=rwaveforms[decayIdx][coboIdx]->waveform[asadIdx*4+agetIdx][chanIdx][buck];
                        lccounter[3]++;
                    }else if(mapchantomm->pxidx[asadIdx][agetIdx][dchanIdx]==68){
                        for(Int_t buck=Bi; buck<Bf; buck++) lcwaveforms[4][buck]+=rwaveforms[decayIdx][coboIdx]->waveform[asadIdx*4+agetIdx][chanIdx][buck];
                        lccounter[4]++;
                    }else if(mapchantomm->pxidx[asadIdx][agetIdx][dchanIdx]==69){
                        for(Int_t buck=Bi; buck<Bf; buck++) lcwaveforms[5][buck]+=rwaveforms[decayIdx][coboIdx]->waveform[asadIdx*4+agetIdx][chanIdx][buck];
                        lccounter[5]++;
                    }
                }
            }
        }
        if(enable2pmode==0){
            for(Int_t pxidx=0;pxidx<6;pxidx++){
                if(lccounter[pxidx]>0){
                    for(Int_t buck=Bi; buck<Bf; buck++){
                        lcwaveforms[pxidx][buck]/=lccounter[pxidx];
                    }
                }
            }
            for(Int_t pxidx=0;pxidx<6;pxidx++){
                if(lccounter[pxidx]>0){
                    for(Int_t buck=Bi; buck<Bf; buck++){
                        if(bmax[pxidx]<lcwaveforms[pxidx][buck]){
                            bmax[pxidx]=lcwaveforms[pxidx][buck];
                            btime[pxidx]=buck;
                        }
                    }
                }
            }
            bestbtime = 0;
            for(Int_t pxidx=0;pxidx<6;pxidx++){
                if(lccounter[pxidx]>0){
                    if(TMath::Abs(Beam_med-bestbtime)>TMath::Abs(Beam_med-btime[pxidx])) bestbtime=btime[pxidx];
                }
            }

            //for(Int_t pxidx=0;pxidx<6;pxidx++){
            //  if(lccounter[pxidx]>0){
            //    cout << reventIdx << " " << pxidx << " " << btime[pxidx] << " " << bestbtime << " " << bmax[pxidx] << endl;
            //  }
            //}
        }

        for(Int_t i=0; i<rGETMul; i++){
            if(wfdvalue[i]>40){
                frameIdx = rGETFrameNo[i];
                decayIdx = rGETDecayNo[i];
                coboIdx = rGETCobo[i];
                asadIdx = rGETAsad[i];
                agetIdx = rGETAget[i];
                chanIdx = rGETChan[i];
                if(ignoremm==1 && coboIdx==0) continue;// skip MM waveform data
                rwaveforms[decayIdx][coboIdx]->isDecay[asadIdx*4+agetIdx][chanIdx]=frameIdx%2;
                if(enable2pmode==1){
                    if(frameIdx%2 == L1Aflag){
                        decayIdx=0;
                    }else{
                        decayIdx=1;
                    }
                }else{
                    decayIdx=0;
                }
                rwaveforms[decayIdx][coboIdx]->frameIdx = frameIdx;
                rwaveforms[decayIdx][coboIdx]->decayIdx = decayIdx;
                //cout << coboIdx << " " << asadIdx << " " << agetIdx << " " << chanIdx <<endl;
                //      if(decayIdx==1) cout<<"Why do we see them here but not after?"<<endl;
                GetAverageFPN(decayIdx,coboIdx,asadIdx,agetIdx);

                if(chanIdx<11) {
                    dchanIdx = chanIdx;
                }else if(chanIdx>11 && chanIdx<22) {
                    dchanIdx = chanIdx - 1;
                }else if(chanIdx>22 && chanIdx<45) {
                    dchanIdx = chanIdx - 2;
                }else if(chanIdx>45 && chanIdx<56) {
                    dchanIdx = chanIdx - 3;
                }else if(chanIdx>56) {
                    dchanIdx = chanIdx - 4;
                }
                bestbtime=180;
                if(coboIdx==0 && mapchantomm->pxidx[asadIdx][agetIdx][dchanIdx]>=64 && mapchantomm->pxidx[asadIdx][agetIdx][dchanIdx]<=69){
                    mm_mintime = bestbtime-Beam_window;
                    mm_maxtime = bestbtime+Beam_window;
                }else{
                    mm_mintime = bestbtime-Particle_window;
                    mm_maxtime = bestbtime+Particle_window;
                }
                GetEnergyTime(decayIdx,coboIdx,asadIdx,agetIdx,chanIdx);
                if(coboIdx==1&&asadIdx==0&&rwaveforms[decayIdx][coboIdx]->energy[asadIdx*4+agetIdx][chanIdx]>3000){
                    rwaveforms[decayIdx][coboIdx]->hasHit[asadIdx*4+agetIdx] = false;
                    rwaveforms[decayIdx][0]->isRejected = true;
                    rwaveforms[decayIdx][1]->isRejected = true;
                }
                if(enablehist==1 && rwaveforms[decayIdx][0]->isRejected == false){
                }
                if(enablehist==1){
                    //cout<<rwaveforms[decayIdx][coboIdx]->energy[asadIdx*4+agetIdx][chanIdx]<<endl;
                    hGET_EHitPattern2D->Fill(coboIdx*2000+asadIdx*500+agetIdx*100+chanIdx,rwaveforms[decayIdx][coboIdx]->energy[asadIdx*4+agetIdx][chanIdx]);
                    hGET_THitPattern2D->Fill(coboIdx*2000+asadIdx*500+agetIdx*100+chanIdx,rwaveforms[decayIdx][coboIdx]->time[asadIdx*4+agetIdx][chanIdx]);
                }
                //WaveletFilter(decayIdx,coboIdx);
                //cout << "Done with energy/time" << endl;
                if(enabledraww==1){
                    DrawWaveForm(decayIdx,coboIdx,asadIdx,agetIdx,chanIdx);
                    //          cout << "Done with drawing waveform for " << reventIdx<< endl;
                }
                ResetHitPattern();
                DrawHitPattern(decayIdx, coboIdx);
            }
        }
        if(enabledraww==1){
            hWaveFormbyEvent[evtcounter%16]->SetTitle(Form("hWaveFormbyEvent(EvtNo=%d);ADC Channel;Counts [D2PTime=%d usec]",reventIdx,int(rd2ptime/1000)));
            hCorrWaveFormbyEvent[evtcounter%16]->SetTitle(Form("hCorrWaveFormbyEvent(EvtNo=%d);ADC Channel;Counts [ D2PTime=%d usec]",reventIdx,int(rd2ptime/1000)));
        }
        //cout << "Done with energy/time" << endl;
        //WaveformShapeFilter(coboIdx);
        //DrawPSDFilter(coboIdx);
        //cout << "Done with filters" << endl;

        goodx6csievt=0;
        goodx6evt=0;
        FindX6Hits();
        //if(goodx6csievt==0) { evtcounter--; continue; }
        //if(goodx6evt==0) { evtcounter--; continue; }
        //else {cout << "Check out the vigru!!" << endl; }
        //if(si_tracks->hasX6L>0 || si_tracks->hasX6BL>0 || si_tracks->hasX6R>0 || si_tracks->hasX6BR>0)
        if(si_tracks->hasFWC>0 || 1)
        {
            ResetTrackHist();
            FillTrack();
            //cout << "Done with filling tracks" << endl;
            if(enabletrack==1){
                FindBoxCorner();
                //cout << "Done with finding corners using a box method" << endl;
                //ReplaceEnergy();
                //ReplaceEnergybyRatio();
                //cout << "Done with replacing energy of strips/chain from slop information" << endl;

                //if(mm_tracks->hasTrack>0 && si_tracks->hasTrack>0)
                if(mm_tracks->hasTrack>0)
                {
                    //if((!mm_tracks->hasOverflow)||
                    //(mm_tracks->hasOverflow && si_tracks->agetid==0 &&
                    //si_tracks->chanid!=1 && si_tracks->chanid!=20))
                    if(!mm_tracks->hasOverflow)
                    {
                        DrawSiEvsCsIE();
                        //cout << "Done with drawing SiEvsCsIE" << endl;
                        goodsicsipevt=1;
                        if(IsDecayEvt==true || true){
                            //cout<<"Decay event"<<endl;
                            if(goodsicsipevt==1 || 1){
                                goodsicsipevtidx++;
                                if(enablecleantrack==1) CleanTrack();
                                if(enable2pmode==1){
                                    FillDecayFlag();
                                    cout << "Done with filling decay flags" << endl;
                                    Sum2pEnergy();
                                }
                                //cout << "Done with cleaning tracks" << endl;
                                FilldEvsE();
                                //cout << "Done with filling dE vs E" << endl;
                                //cout << "Done with resetting track histograms" << endl;
                                ChangeTrackHistTitle();
                                //cout << "M2:Good Event Found! (Idx=" << goodsicsipevtidx << ", EvtNo=" << reventIdx << ")" << endl;
                                DrawTrack();
                                //cout << "Done with drawing tracks" << endl;
                                DrawSumEnergyTrack();
                                //cout << "Done with drawing sum energy tracks" << endl;
                                DrawdEvsE();
                                //cout << "Done with drawing dEvsE" << endl;
                                if(enable2pmode==1){
                                    cout<<"Drawing decay track"<<endl;
                                    DrawTrack2pMode();
                                }else{
                                    //HoughTransform();
                                    //cout << "Done with Hough Transformation" << endl;
                                    //GetTrackPosYLimit();
                                    //cout << "Done with getting track posy min/max" << endl;
                                    //GetXYZTrack();
                                    //cout << "Done with getting track" << endl;
                                }
                                if(goodx6csievt>0) {
                                    cout<<"Good X6 CsI event"<<endl;
                                    //      goodevtcounter++;
                                }
                                goodevtcounter++;

                            }
                        }
                        ResetTrack();
                        //cout << "Done with resetting tracks" << endl;
                    }else{
                        badevtcounter++;
                        ResetTrack();
                    }
                }else{
                    badevtcounter++;
                    ResetTrack();
                }
                //cout << "Done with finding tracks" << endl;
            }
            ResetdEvsE();
        }
    }
    bucketmax = oldbucketmax;
}

void LKFrameBuilder::RootReadWriteEvent() {
  Int_t frameIdx = 0;
  Int_t decayIdx = 0;
  L1Aflag = 0;
  Int_t coboIdx = 0;
  Int_t asadIdx = 0;
  Int_t agetIdx = 0;
  Int_t chanIdx = 0;
  UInt_t percent = 0;
  UInt_t Entries = 1;
  UInt_t oldbucketmax = bucketmax;
  if(enable2pmode==1){
    bucketmax = bucketmax/2;
    cout << "bucketmax set to " << bucketmax << endl;
  }
  ostringstream s4out;
  ostringstream s5out;
  ofstream f4out;
  ofstream f5out;
  RootRInit();
  RootRWInit();
  if(fNumberEvents>0){
    Entries = fNumberEvents;
    percent = Entries/10;
  }
  SetResponseWaveform();
  GetMaxResponseWaveform();
  GetSigmaResponseWaveform();
  for(int i=1;i<fNumberEvents;i++){
    fInputFile->cd();
    fInputTree->GetEntry(i);
    reventIdx = rGETEventIdx;
    if(RootRRecovered()){
      if(LasteventIdx>=reventIdx) continue;
      else LasteventIdx = reventIdx;
    }else{
      LasteventIdx = reventIdx;
    }
    if(IsFirstevent){
      FirsteventIdx = reventIdx;
      evtcounter=0;
      goodevtcounter=0;
      gatedevtcounter=0;
      badevtcounter=0;
      IsFirstevent=false;
      RootRResetWaveforms();
    }else{
      evtcounter++;
      RootRResetWaveforms();
    }
    if((reventIdx-FirsteventIdx)%100==0){
      if(reventIdx==FirsteventIdx){
        cout << Form("M3:Starting from Event No. %d.",reventIdx) << endl;
      }else{
        cout << Form("M3:Upto Event No. %d Processed (from %d)..",reventIdx, FirsteventIdx) << endl;
      }
    }

    //cout << rGETEventIdx << " " << rGETD2PTime << " " << rGETMul << endl;
    rd2ptime = rGETD2PTime;
    rtstmp = rGETTimeStamp;
    if(enable2pmode==1){
      if(rd2ptime>0) hMM_D2PTime->Fill(rd2ptime);
      else cout << reventIdx << " " << rd2ptime << " " << rGETMul << endl;

      printed=0;
      for(Int_t i=0; i<rGETMul; i++){
        frameIdx = rGETFrameNo[i];
        decayIdx = rGETDecayNo[i];
        coboIdx = rGETCobo[i];
        asadIdx = rGETAsad[i];
        agetIdx = rGETAget[i];
        chanIdx = rGETChan[i];
        if(coboIdx==coboIC && asadIdx==asadIC && chanIdx==chanIC){
          for(int j=0;j<bucketmax;j++){
            if(rGETWaveformY[i][j] < valueIC_min){
//            L1Aflag = frameIdx%2;
            }
          }
//        cout << reventIdx << " " << coboIdx << " " << asadIdx << " " << chanIdx << " " << frameIdx << " " << L1Aflag << " " << decayIdx << endl;
        }else if(coboIdx==0){
          for(int j=0;j<bucketmax;j++){
            if(printed==0 && rGETWaveformY[i][j]>4000){
              cout << "MM2: " << reventIdx << " " << coboIdx << " " << asadIdx << " " << chanIdx << " " << frameIdx << " " << decayIdx << endl;
              printed++;
            }
          }
        }
      }
    }

    for(Int_t i=0; i<rGETMul; i++){
      decayIdx = rGETDecayNo[i];
      coboIdx = rGETCobo[i];
      if(ignoremm==1 && coboIdx==0) continue;// skip MM waveform data
      wfmaxvalue[i]=-10000;
      wfminvalue[i]=10000;
      wfbaseline[i] = 0;
      Int_t baselinecounter = 0;
      wfdvalue[i] = 0;
      for(int j=0;j<bucketmax;j++){
        if(rGETWaveformY[i][j]>0 && rGETWaveformY[i][j]<4095){
          if(wfmaxvalue[i]<rGETWaveformY[i][j]) wfmaxvalue[i] = rGETWaveformY[i][j];
          if(wfminvalue[i]>rGETWaveformY[i][j]) wfminvalue[i] = rGETWaveformY[i][j];
          if(j>=20+decayIdx*256 && j<30+decayIdx*256){
            wfbaseline[i] += rGETWaveformY[i][j];
            baselinecounter++;
          }
        }
      }
      if(baselinecounter>0) wfbaseline[i] /= baselinecounter;
      wfdvalue[i] = wfmaxvalue[i] - wfminvalue[i];
      //if(coboIdx==0 && (wfdvalue[i]>(wfmaxvalue[i]-wfbaseline[i]+200))) wfdvalue[i] = -1;
      //if(coboIdx==1 && asadIdx==0 && (agetIdx==1||agetIdx==3) && (wfdvalue[i]>(wfmaxvalue[i]-wfbaseline[i]+100))) wfdvalue[i] = -1;
      //cout << reventIdx << ", wfdvalue[" << i << "]= " << wfdvalue[i] << ",  wfbaseline[" << i << "]= " << wfbaseline[i] << " " << (wfmaxvalue[i]-wfbaseline[i]+300) << endl;
    }
    for(Int_t i=0; i<rGETMul; i++){
      //cout << "wfdvalue[" << i << "]= " << wfdvalue[i] << endl;
      //if(wfdvalue[i]>40)
      if(wfdvalue[i]>0)
      {
        decayIdx = rGETDecayNo[i];
        frameIdx = rGETFrameNo[i];
        coboIdx = rGETCobo[i];
        asadIdx = rGETAsad[i];
        agetIdx = rGETAget[i];
        chanIdx = rGETChan[i];
        if(ignoremm==1 && coboIdx==0) continue;// skip MM waveform data
        if(enable2pmode==1){
          if(coboIdx==coboIC && asadIdx==asadIC && chanIdx==chanIC){
//              decayIdx=0;
          }else{
            if(frameIdx%2 == L1Aflag){
              decayIdx=0;
            }else{
              decayIdx=1;
            }
          }
        }else{
          decayIdx=0;
        }
        rwaveforms[decayIdx][coboIdx]->decayIdx = decayIdx;
        if((chanIdx==11||chanIdx==22||chanIdx==45||chanIdx==56)){
                rwaveforms[decayIdx][coboIdx]->hasFPN[asadIdx*4+agetIdx] = true;
        }
        //cout << reventIdx << " " << coboIdx << " " << asadIdx << " " << agetIdx << " " << chanIdx << " " << rwaveforms[decayIdx][coboIdx]->isDecay[asadIdx*4+agetIdx][chanIdx] << endl;
        if(enable2pmode==0){
          for(int j=0;j<bucketmax;j++){
            rwaveforms[decayIdx][coboIdx]->waveform[asadIdx*4+agetIdx][chanIdx][j] = rGETWaveformY[i][j];
          }
        }else{
          //if(decayIdx==1) cout << "Yo:" << reventIdx << " " << decayIdx << " " << coboIdx << " " << asadIdx << " " << agetIdx << " " << chanIdx << " " << rwaveforms[decayIdx][coboIdx]->isDecay[asadIdx*4+agetIdx][chanIdx] << endl;
          for(int j=0;j<bucketmax;j++){
            rwaveforms[decayIdx][coboIdx]->waveform[asadIdx*4+agetIdx][chanIdx][j] = rGETWaveformY[i][j];
            //cout << j << " " << rwaveforms[decayIdx][coboIdx]->waveform[asadIdx*4+agetIdx][chanIdx][j] << " " << rGETWaveformY[i][j-decayIdx*256] << endl;
          }
        }

/*
        if(coboIdx==1&&asadIdx==0){
          s4out << reventIdx << " " << coboIdx << " " << asadIdx << " " << agetIdx << " " << chanIdx << " 4" << endl;
        }else if(coboIdx==1&&asadIdx==1&&chanIdx==33){
          s5out << reventIdx << " " << coboIdx << " " << asadIdx << " " << agetIdx << " " << chanIdx << " 5" << endl;
        }
*/
      }
    }

/*
    if(s4out.str().size()!=0){
      f4out.open("silabels.txt", std::ofstream::out|std::ofstream::app);
      f4out << s4out.str() << endl;
      f4out.close();
      s4out.clear();
      s4out.str("");
    }
    if(s5out.str().size()!=0){
      f5out.open("iclabels.txt", std::ofstream::out|std::ofstream::app);
      f5out << s5out.str() << endl;
      f5out.close();
      s5out.clear();
      s5out.str("");
    }
*/
    //cout << "Done with reading events from root tree(" << rGETMul << ")" << endl;

    wGETEventIdx = rGETEventIdx;
    wGETD2PTime = rGETD2PTime;
    wGETTimeStamp = rGETTimeStamp;
    for(Int_t i=0; i<rGETMul; i++){
      wfbaseline[i]=0;
      //cout << "wfdvalue[" << i << "]= " << wfdvalue[i] << endl;
      if(wfdvalue[i]>40){
        frameIdx = rGETFrameNo[i];
        decayIdx = rGETDecayNo[i];
        coboIdx = rGETCobo[i];
        asadIdx = rGETAsad[i];
        agetIdx = rGETAget[i];
        chanIdx = rGETChan[i];
        if(ignoremm==1 && coboIdx==0) continue;// skip MM waveform data
        maxValueBucket = 100+decayIdx*256;
        wfbaseline[i] = GetBaseline(decayIdx,coboIdx,asadIdx,agetIdx,chanIdx);
        //cout << reventIdx << ", wfdvalue[" << i << "]= " << wfdvalue[i] << ",  wfbaseline[" << i << "]= " << wfbaseline[i] << endl;
        if(wfdvalue[i]>(wfmaxvalue[i]-wfbaseline[i]+300)) wfdvalue[i] = -1;
      }
      if(wfdvalue[i]>40){
        if(enable2pmode==1){
          if(frameIdx%2 == L1Aflag){
            decayIdx=0;
          }else{
            decayIdx=1;
          }
        }else{
          decayIdx=0;
        }
        rwaveforms[decayIdx][coboIdx]->frameIdx = frameIdx;
        rwaveforms[decayIdx][coboIdx]->decayIdx = decayIdx;
        GetAverageFPN(decayIdx,coboIdx,asadIdx,agetIdx);
        GetEnergyTime(decayIdx,coboIdx,asadIdx,agetIdx,chanIdx);
        //WaveletFilter(decayIdx,coboIdx);
        wGETFrameNo[i]=frameIdx;
        wGETDecayNo[i]=decayIdx;
        wGETL1Aflag[i]=L1Aflag;
        wGETCobo[i]=coboIdx;
        wGETAsad[i]=asadIdx;
        wGETAget[i]=agetIdx;
        wGETChan[i]=chanIdx;
        wGETEnergy[i]=rwaveforms[decayIdx][coboIdx]->energy[asadIdx*4+agetIdx][chanIdx];
        wGETTime[i]=rwaveforms[decayIdx][coboIdx]->time[asadIdx*4+agetIdx][chanIdx];
        wGETMul++;
        if(chanIdx!=11&&chanIdx!=22&&chanIdx!=45&&chanIdx!=56) wGETHit++;
        //cout << "Done with energy/time" << endl;
      }
    }
    RootRWriteEvent();
    RootRWReset();
  }
  bucketmax = oldbucketmax;
}

void LKFrameBuilder::RootReadHistEvent()
{
    Int_t frameIdx = 0;
    Int_t decayIdx = 0;
    L1Aflag = 0;
    Int_t coboIdx = 0;
    Int_t asadIdx = 0;
    Int_t agetIdx = 0;
    Int_t chanIdx = 0;
    UInt_t percent = 0;
    UInt_t Entries = 1;
    UInt_t oldbucketmax = bucketmax;
    if(enable2pmode==1){
        bucketmax = bucketmax/2;
        cout << "bucketmax set to " << bucketmax << endl;
    }
    ostringstream s4out;
    ostringstream s5out;
    ofstream f4out;
    ofstream f5out;
    RootRHInit();
    if(fNumberEvents>0){
        Entries = fNumberEvents;
        percent = Entries/10;
    }
    for(int i=1;i<fNumberEvents;i++){
        fInputFile->cd();
        fInputTree->GetEntry(i);
        reventIdx = rGETEventIdx;
        goodsicsievt=0;
        goodsicsipevt=0;
        if(RootRRecovered()){
            if(LasteventIdx>=reventIdx) continue;
            else LasteventIdx = reventIdx;
        }else{
            LasteventIdx = reventIdx;
        }
        if(IsFirstevent){
            FirsteventIdx = reventIdx;
            evtcounter=0;
            goodevtcounter=0;
            gatedevtcounter=0;
            badevtcounter=0;

            IsFirstevent=false;
            RootRHResetWaveforms();
        }else{
            evtcounter++;
            RootRHResetWaveforms();
        }
        if((reventIdx-FirsteventIdx)%2000==0){
            if(reventIdx==FirsteventIdx){
                cout << Form("M4: Starting from Event No. %d.",reventIdx) << endl;
            }else{
                cout << Form("M4: Upto Event No. %d Processed (from %d)..",reventIdx, FirsteventIdx) << endl;
            }
        }

        //cout << rGETEventIdx << " " << rGETD2PTime << " " << rGETMul << endl;
        rd2ptime = rGETD2PTime;
        rtstmp = rGETTimeStamp;
        if(enable2pmode==1){
            if(rd2ptime>0) hMM_D2PTime->Fill(rd2ptime);
            else cout << reventIdx << " " << rd2ptime << " " << rGETMul << endl;
            printed=0;
        }

        for(Int_t i=0; i<rGETMul; i++){
            frameIdx = rGETFrameNo[i];
            decayIdx = rGETDecayNo[i];
            coboIdx = rGETCobo[i];
            asadIdx = rGETAsad[i];
            agetIdx = rGETAget[i];
            chanIdx = rGETChan[i];
            if(ignoremm==1 && coboIdx==0) continue;// skip MM waveform data
            L1Aflag = rGETL1Aflag[i];
            rwaveforms[decayIdx][coboIdx]->frameIdx = frameIdx;
            rwaveforms[decayIdx][coboIdx]->decayIdx = decayIdx;
            rwaveforms[decayIdx][coboIdx]->hasHit[asadIdx*4+agetIdx] = true;

            rwaveforms[decayIdx][coboIdx]->hasSignal[asadIdx*4+agetIdx][chanIdx] = true;
            rwaveforms[decayIdx][coboIdx]->energy[asadIdx*4+agetIdx][chanIdx] = rGETEnergy[i];
            rwaveforms[decayIdx][coboIdx]->time[asadIdx*4+agetIdx][chanIdx] = rGETTime[i];
            if(coboIdx==1 && asadIdx==1 && agetIdx==1 && chanIdx==2) {
                rwaveforms[decayIdx][coboIdx]->ICenergy = rwaveforms[decayIdx][coboIdx]->energy[asadIdx*4+agetIdx][chanIdx];
            }
            //cout << "Done with energy/time" << endl;
            DrawHitPattern(decayIdx, coboIdx);
            //cout << "Done with drawhitpattern" << endl;
        }

        FindX6Hits();

        bmpos = 0;
        bmsum1 = 0;
        bmsum2 = 0;
        if(rwaveforms[0][1]->energy[1*4+2][19]>20 && rwaveforms[0][1]->energy[1*4+2][25]>20){
            bmpos = (rwaveforms[0][1]->energy[1*4+2][19]-rwaveforms[0][1]->energy[1*4+2][25])/(rwaveforms[0][1]->energy[1*4+2][19]+rwaveforms[0][1]->energy[1*4+2][25]);
            bmsum1 = (rwaveforms[0][1]->energy[1*4+2][19]+rwaveforms[0][1]->energy[1*4+2][25]);
            if(rwaveforms[0][1]->energy[1*4+2][33]>20) bmsum2 = rwaveforms[0][1]->energy[1*4+2][33];
            hBM_LR->Fill(rwaveforms[0][1]->energy[1*4+2][19],rwaveforms[0][1]->energy[1*4+2][25]);
            hBM_Esum1vsPos->Fill(bmpos, bmsum1);
            hBM_Esum2vsEsum1->Fill(bmsum1, bmsum2);
        }

        FillTrack();
        //cout << "Done with filling tracks" << endl;

        if(enabletrack==1){
            FindBoxCorner();
            //cout << "Done with finding corners using a box method" << endl;
            //ReplaceEnergy();
            //ReplaceEnergybyRatio();
            //cout << "Done with replacing energy of strips/chain from slop information" << endl;

            //if(mm_tracks->hasTrack>0 && si_tracks->hasTrack>0)
            if(mm_tracks->hasTrack>2)
            {
                //if((!mm_tracks->hasOverflow)||
                //(mm_tracks->hasOverflow && si_tracks->agetid==0 &&
                //si_tracks->chanid!=1 && si_tracks->chanid!=20))
                if(!mm_tracks->hasOverflow)
                {
                    DrawSiEvsCsIE();
                    //cout << "Done with drawing SiEvsCsIE" << endl;
                    goodsicsipevt=1;
                    if(goodsicsipevt==1){
                        goodsicsipevtidx++;
                        if(enablecleantrack==1) CleanTrack();
                        if(enable2pmode==1){
                            FillDecayFlag();
                            cout << "Done with filling decay flags" << endl;
                            Sum2pEnergy();
                        }
                        //cout << "Done with cleaning tracks" << endl;
                        FilldEvsE();
                        //cout << "Done with filling dE vs E" << endl;
                        ResetTrackHist();
                        //cout << "Done with resetting track histograms" << endl;
                        ChangeTrackHistTitle();
                        //cout << "M4:Good Proton Found! (Idx=" << goodsicsipevtidx << ", EvtNo=" << reventIdx << ")" << endl;
                        DrawTrack();
                        //cout << "M4:Done with drawing tracks" << endl;
                        DrawSumEnergyTrack();
                        //cout << "Done with drawing sum energy tracks" << endl;
                        DrawdEvsE();
                        //cout << "Done with drawing dEvsE" << endl;
                        if(enable2pmode==1){
                            DrawTrack2pMode();
                        }else{
                            //HoughTransform();
                            //cout << "Done with Hough Transformation" << endl;
                            //GetTrackPosYLimit();
                            //cout << "Done with getting track posy min/max" << endl;
                            //GetXYZTrack();
                            //cout << "Done with getting track" << endl;
                        }
                        goodevtcounter++;
                    }
                    ResetTrack();
                    //cout << "Done with resetting tracks" << endl;
                }else{
                    badevtcounter++;
                    ResetTrack();
                }
            }else{
                badevtcounter++;
                ResetTrack();
            }
            ResetdEvsE();
            //cout << "Done with finding tracks" << endl;
        }
    }
    bucketmax = oldbucketmax;
}

//We will get the averaged FPN waveform for each Aget.
void LKFrameBuilder::GetAverageFPN(Int_t decayIdx, Int_t cobo, Int_t asad, Int_t aget)
{
    if(rwaveforms[decayIdx][cobo]->hasFPN[asad*4+aget] && !rwaveforms[decayIdx][cobo]->doneFPN[asad*4+aget]){
        for(UInt_t buck=0;buck<bucketmax;buck++){
            rwaveforms[decayIdx][cobo]->fpnwaveform[asad*4+aget][buck] = rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][11][buck];
            rwaveforms[decayIdx][cobo]->fpnwaveform[asad*4+aget][buck] += rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][22][buck];
            rwaveforms[decayIdx][cobo]->fpnwaveform[asad*4+aget][buck] += rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][45][buck];
            rwaveforms[decayIdx][cobo]->fpnwaveform[asad*4+aget][buck] += rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][56][buck];
            rwaveforms[decayIdx][cobo]->fpnwaveform[asad*4+aget][buck] = rwaveforms[decayIdx][cobo]->fpnwaveform[asad*4+aget][buck]/4;
        }
        //rwaveforms[decayIdx][cobo]->doneFPN[asad*4+aget] = true;
    }
}

void LKFrameBuilder::DrawWaveForm(Int_t decayIdx, Int_t cobo, Int_t asad, Int_t aget, Int_t chan)
{
    Int_t maxyval=0;
    Int_t baseline = rwaveforms[decayIdx][cobo]->baseline[asad*4+aget][chan];
    if(chan==11 || chan==22 || chan==45 || chan==56) return; // We want to skip the FPN channels.
    if(!rwaveforms[decayIdx][cobo]->hasSignal[asad*4+aget][chan]) return; // skip signals that did not fire
    //shaped waveform type
    if(cobo==0){
        rftype=0;
    }else if(cobo==1){
        if(asad==0){
            rftype=2;
        }else if(asad==1){
            rftype=1;
        }
    }

    Int_t mintime=0;
    Int_t maxtime=0;
    if(cobo==0){
        mintime=mm_mintime;
        maxtime=mm_maxtime;
    }else{
        mintime=si_mintime;
        maxtime=si_maxtime;
    }
    if(cobo==1 && asad>0) {//For the non-Si have different times
        mintime=0;
        maxtime=512;
    }

    //if(cobo>=0 && rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]>mm_minenergy && rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]<mm_maxenergy && rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan]>mintime+decayIdx*20 && rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan]<maxtime+decayIdx*256)
    if((cobo==0 && rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]>mm_minenergy && rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]<mm_maxenergy && rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan]>mintime+decayIdx*20 && rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan]<maxtime+decayIdx*256)||cobo>0)
    {
        hWaveForm[cobo*maxasad*4+asad*4+aget]->Reset();
        hCorrWaveForm[cobo*maxasad*4+asad*4+aget]->Reset();
        hCorrWaveFormDec[cobo*maxasad*4+asad*4+aget]->Reset();
        hCorrWaveFormFit[cobo*maxasad*4+asad*4+aget]->Reset();
        hCorrWaveFormRDF[cobo*maxasad*4+asad*4+aget]->Reset();
        //for(int i=0;i<=512;i++) hWaveForm[cobo*maxasad*4+asad*4+aget]->SetBinContent(i,0);
        //for(int i=0;i<=512;i++){
        //  hCorrWaveForm[cobo*maxasad*4+asad*4+aget]->SetBinContent(i,0);
        //  hCorrWaveFormDec[cobo*maxasad*4+asad*4+aget]->SetBinContent(i,0);
        //  hCorrWaveFormFit[cobo*maxasad*4+asad*4+aget]->SetBinContent(i,0);
        //  for(int j=0;j<=512;j++){
        //    hCorrWaveFormRDF[cobo*maxasad*4+asad*4+aget]->SetBinContent(i,j,0);
        //  }
        //}
        //cout << evtcounter << " " << cobo << " " << asad << " " << aget << " " << chan << " " << rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][0] << endl;

        //for(UInt_t buck=0;buck<bucketmax;buck++)
        for(UInt_t buck=mintime;buck<maxtime;buck++)
        {
            if(cobo>=0){
                hWaveForm[cobo*maxasad*4+asad*4+aget]->Fill(buck,rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck]);
                hCorrWaveForm[cobo*maxasad*4+asad*4+aget]->Fill(buck,rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck]);
                //if(buck<responsesample[0][6]) hCorrWaveForm[cobo*maxasad*4+asad*4+aget]->Fill(buck,response[rftype][buck]-decoffset+baseline);
                hCorrWaveFormDec[cobo*maxasad*4+asad*4+aget]->Fill(buck,corrwaveformdec[buck]-decoffset);
                hCorrWaveFormFit[cobo*maxasad*4+asad*4+aget]->Fill(buck,corrwaveformfit[buck]);
                hCorrWaveFormRDF[cobo*maxasad*4+asad*4+aget]->Fill(buck,baseline);
                hCorrWaveFormRDF[cobo*maxasad*4+asad*4+aget]->Fill(buck,rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck]);
                hCorrWaveFormRDF[cobo*maxasad*4+asad*4+aget]->Fill(buck,response[rftype][buck],10);
                hCorrWaveFormRDF[cobo*maxasad*4+asad*4+aget]->Fill(buck,corrwaveformdec[buck],50);
                hCorrWaveFormRDF[cobo*maxasad*4+asad*4+aget]->Fill(buck,corrwaveformfit[buck],100);
                hWaveFormbyEvent[evtcounter%16]->Fill(buck+decayIdx*256,rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck]);
                //        cout<<"BEEP"<<endl;
                //      if(evtcounter%16>1)cout<<buck<<"\t"<<rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck]<<endl;
                hCorrWaveFormbyEvent[evtcounter%16]->Fill(buck+decayIdx*256,rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck]-baseline);
                //hCorrWaveFormbyEvent[evtcounter%16]->Fill(buck+decayIdx*256,response[rftype][buck]-baseline);
                if(cobo==1&&asad==0&&aget==0&&chan==0) hFPNWaveFormAll[cobo*maxasad*4+asad*4+aget]->Fill(buck,rwaveforms[decayIdx][cobo]->fpnwaveform[asad*4+aget][buck]);
                if(cobo==1&&asad==1&&aget==2&&chan==19) hWaveFormIC->Fill(buck+decayIdx*256,rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck]);
            }
        }
        hWaveForm[cobo*maxasad*4+asad*4+aget]->SetTitle(Form("hWaveForm_%d_%d_%d_%d(EvtNo=%d);ADC Channel;Counts",cobo,asad,aget,chan,reventIdx));
        hCorrWaveForm[cobo*maxasad*4+asad*4+aget]->SetTitle(Form("hCorrWaveForm_%d_%d_%d_%d(EvtNo=%d);ADC Channel;Counts",cobo,asad,aget,chan,reventIdx));
        hCorrWaveFormDec[cobo*maxasad*4+asad*4+aget]->SetTitle(Form("hCorrWaveFormDec_%d_%d_%d_%d(EvtNo=%d);ADC Channel;Counts",cobo,asad,aget,chan,reventIdx));
        hCorrWaveFormFit[cobo*maxasad*4+asad*4+aget]->SetTitle(Form("hCorrWaveFormFit_%d_%d_%d_%d(EvtNo=%d);ADC Channel;Counts",cobo,asad,aget,chan,reventIdx));
        hCorrWaveFormRDF[cobo*maxasad*4+asad*4+aget]->SetTitle(Form("hCorrWaveFormFit_%d_%d_%d_%d(EvtNo=%d);ADC Channel;Counts (bl=%d",cobo,asad,aget,chan,reventIdx,baseline));
        hWaveForm[cobo*maxasad*4+asad*4+aget]->SetOption("hist");
        hCorrWaveForm[cobo*maxasad*4+asad*4+aget]->SetOption("hist");
        hCorrWaveFormDec[cobo*maxasad*4+asad*4+aget]->SetOption("hist");
        hCorrWaveFormFit[cobo*maxasad*4+asad*4+aget]->SetOption("hist");
        UInt_t nfound = 0;
        //UInt_t nfound = sCorrWaveForm[cobo*maxasad*4+asad*4+aget]->Search(hCorrWaveForm[cobo*maxasad*4+asad*4+aget],2,"",0.4);
        if(nfound>100){
            rwaveforms[decayIdx][cobo]->hasSignal[asad*4+aget][chan] = 0;
        }
    }
}

void LKFrameBuilder::GetCorrWaveform(Int_t decayIdx, Int_t cobo, Int_t asad, Int_t aget, Int_t chan)
{
    Int_t minvalue=10000;
    if(enable2pmode==1) decayIdx = rwaveforms[decayIdx][cobo]->decayIdx;
    for(UInt_t buck=0;buck<bucketmax;buck++){
        // Flip
        if(cobo==0) {
            rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck] =  rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck] - rwaveforms[decayIdx][cobo]->fpnwaveform[asad*4+aget][buck]; // Subtract the averaged FPN value (Negative Polarity)

        }
        if(cobo==2) {//JEB
            if(aget!=0){
                rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck] =  -1*(rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck] - rwaveforms[decayIdx][cobo]->fpnwaveform[asad*4+aget][buck]); // Subtract the averaged FPN value (Positive polarity)
            }
            if(aget==0) {
                rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck] =  rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck] - rwaveforms[decayIdx][cobo]->fpnwaveform[asad*4+aget][buck]; // Subtract the averaged FPN value (Negative Polarity)
            }
        }
        if(cobo==1) {
            if(asad==0) {
                if(aget==0 || aget==2) {
                    //rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck] =  -1*(rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck] - rwaveforms[decayIdx][cobo]->fpnwaveform[asad*4+aget][buck]); // Subtract the averaged FPN value (Positive polarity)
                    rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck] =  -1*(rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck] - 4096); // Subtract the averaged FPN value (Positive polarity)
                    //rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck] =  rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck] - rwaveforms[decayIdx][cobo]->fpnwaveform[asad*4+aget][buck]; // Subtract the averaged FPN value (Negative Polarity)
                }else{
                    //rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck] =  -1*(rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck] - rwaveforms[decayIdx][cobo]->fpnwaveform[asad*4+aget][buck]); // Subtract the averaged FPN value (Positive polarity)
                    rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck] =  rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck] - rwaveforms[decayIdx][cobo]->fpnwaveform[asad*4+aget][buck]; // Subtract the averaged FPN value (Negative Polarity)
                }
            }else if(asad==2 || asad==3) {
                if(aget==2 || aget==3) {
                    rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck] =  -1*(rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck] - 4096); // Subtract the averaged FPN value (Positive polarity)
                }
                else {
                    rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck] =  rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck] - rwaveforms[decayIdx][cobo]->fpnwaveform[asad*4+aget][buck]; // Subtract the averaged FPN value (Negative Polarity)
                }
            }else{
                rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck] =  -1*(rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck] - rwaveforms[decayIdx][cobo]->fpnwaveform[asad*4+aget][buck]); // Subtract the averaged FPN value (Positive polarity)
                //rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck] =  rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck] - rwaveforms[decayIdx][cobo]->fpnwaveform[asad*4+aget][buck]; // Subtract the averaged FPN value (Negative Polarity)
            }
        }
        //    if(cobo==0 && (asad==0||asad==1) && aget>1) {
        //      rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck] =  -1*(rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck] - rwaveforms[decayIdx][cobo]->fpnwaveform[asad*4+aget][buck]); // Subtract the averaged FPN value (Positive polarity)
        //      }
        //      cout<<buck<<"\t"<<cobo<<"\t"<<asad<<"\t"<<aget<<"\t"<<rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck]<<endl;

        if(minvalue>rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck]){
            minvalue = rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck];
        }
    }
    if(minvalue<0){
        for(UInt_t buck=0;buck<bucketmax;buck++){
            rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck] -= minvalue;
        }
    }
}

Double_t LKFrameBuilder::GetBaseline(Int_t decayIdx, Int_t cobo, Int_t asad, Int_t aget, Int_t chan)
{
    Double_t baseline;
    Double_t dbaseline;
    TF1* fWaveform;
    TH1D* hWaveform;
    Double_t sigma=10;

    fWaveform = new TF1("fWaveform", "gaus", 0, 4096);
    fWaveform->SetParameter(2,sigma);
    fWaveform->SetParLimits(2, sigma*0.1, sigma*2);
    hWaveform = new TH1D("hWaveform","hWaveform",1024,0,4096);
    for(int buck=0;buck<512;buck++){
        hWaveform->Fill(rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck]);
    }
    hWaveform->Fit("fWaveform","QN"); //Q=Quiet and N=Do not draw
    dbaseline = TMath::Abs(fWaveform->GetMaximumX() - hWaveform->GetMaximumBin());
    if(dbaseline<100){
        baseline = fWaveform->GetMaximumX();
    }else{
        baseline = hWaveform->GetMaximumBin();
    }
    delete hWaveform;
    delete fWaveform;

    UInt_t bgmin=50+decayIdx*256;
    UInt_t bgmax=bucketmax+decayIdx*256-50;
    UInt_t backgroundSum=0;
    UInt_t backgroundSamples=0;

    for(UInt_t buck=0;buck<bucketmax;buck++){
        if(maxValueBucket>bgmin && maxValueBucket<bgmax){
            if(buck>(bgmin-30) && buck<(bgmin-20)){
                backgroundSum += rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck];
                backgroundSamples++;
            }else if((buck>(bgmax+20) && buck<(bgmax+30))){
                backgroundSum += rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck];
                backgroundSamples++;
            }
        }else if(maxValueBucket<bgmin){
            if((buck>bgmax && buck<(bgmax+20))){
                backgroundSum += rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck];
                backgroundSamples++;
            }
        }else if(maxValueBucket>bgmax){
            if(buck>(bgmin-20) && buck<bgmin){
                backgroundSum += rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck];
                backgroundSamples++;
            }
        }
    }
    baseline = (backgroundSamples>0) ? round(backgroundSum/backgroundSamples) : 0;

    return baseline;
}

void LKFrameBuilder::GetEnergyTime(Int_t decayIdx, Int_t cobo, Int_t asad, Int_t aget, Int_t chan)
{
    if(chan==11 || chan==22 || chan==45 || chan==56) return; // We want to skip the FPN channels.
    Int_t maxvalue=10000;
    Int_t minvalue=10000;
    Int_t corrwaveformintegral = 0;
    Double_t psdratio = 0;
    Int_t baseline = 0;

    corrwaveformintegral = 0;
    psdratio = 0;

    maxValue=-1000;
    maxValueBucket=0;
    maxValuedec=-1000;
    maxValueBucketdec=0;
    maxValuefit=-1000;
    maxValueBucketfit=0;
    rftype=0;

    GetCorrWaveform(decayIdx,cobo,asad,aget,chan);

    Int_t mintime=0;
    Int_t maxtime=0;
    if(cobo==0){
        mintime=mm_mintime;
        maxtime=mm_maxtime;
    }else{
        mintime=si_mintime;
        maxtime=si_maxtime;
    }
    if(cobo==1 && asad>0) {//For the non-Si have different times
        mintime=0;
        maxtime=512;
    }

    // Find the maximum
    Int_t *samplecorrwf = rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan].data();
    //for(UInt_t buck=0;buck<bucketmax;buck++)
    for(UInt_t buck=mintime;buck<maxtime;buck++)
    {
        //   if(cobo==1 && asad==1 && aget==3 && chan==19) std::cout<<buck<<"\t"<<samplecorrwf[buck]<<std::endl;
        if(maxValue<samplecorrwf[buck]) {
            maxValue = samplecorrwf[buck];
            maxValueBucket = buck;
        }
        //if(maxValue<rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck]) {
        //  maxValue = rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck];
        //  maxValueBucket = buck;
        //}
    }

    //Find baseline from corrwaveform
    baseline = GetBaseline(decayIdx,cobo,asad,aget,chan);

    //if(maxValue>3500) energymethod=2;
    //cout << cobo << " " << asad << " " << aget << " " << chan << " " << maxValue << endl;

    if(energymethod>0){
        for(int i=0;i<512;i++){
            corrwaveformdec[i]=0;
            corrwaveformfit[i]=0;
        }

        //shaped waveform type
        if(cobo==0){
            rftype=2;
        }else if(cobo==1 || cobo==2){
            if(asad==0){
                rftype=2;
            }else if(asad==1){
                rftype=1;
            }
        }
        int min_val=405;
        for(UInt_t buck=0;buck<bucketmax;buck++){
            corrwaveformdec[buck] = rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck]-baseline+decoffset;
            if(corrwaveformdec[buck]<min_val && corrwaveformdec[buck]>0) min_val=corrwaveformdec[buck];
            //    cout<<"IN1\t"<<buck<<"\t"<<corrwaveformdec[buck]<<"\t"<<min_val<<endl;

        }
        min_val+=50;
        for(UInt_t buck=0;buck<bucketmax;buck++) {
            corrwaveformdec[buck]-=min_val;
            if(buck<100) corrwaveformdec[buck]=0;
            //      cout<<"IN CORR RESPONSE\t"<<buck<<"\t"<<corrwaveformdec[buck]<<"\t"<<response[rftype][buck]<<endl;
        }
        //    cout << "Deconvolution Result: " << sCorrWaveForm[cobo*maxasad*4+asad*4+aget]->Deconvolution(corrwaveformdec,response[rftype],bucketmax,deconvrep[rftype],deconviter[rftype],deconvboost[rftype]) << endl;
        sCorrWaveForm[cobo*maxasad*4+asad*4+aget]->Deconvolution(corrwaveformdec,response[rftype],bucketmax,deconvrep[rftype],deconviter[rftype],deconvboost[rftype]);
        for(UInt_t buck=responsesample[rftype][6]*0.8+decayIdx*256;buck<bucketmax+decayIdx*256;buck++){
            //  cout<<buck<<"\t"<<corrwaveformdec[buck]<<endl;
            // Find the maximum
            if(maxValuedec<corrwaveformdec[buck]) {
                maxValuedec = corrwaveformdec[buck];
                maxValueBucketdec = buck;
            }
        }
        maxValuedec = rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][maxValueBucketdec];

        GetEnergybyFitWaveform(rftype,decayIdx,cobo,asad,aget,chan,maxValueBucket); //using ShaperF_Fusion method
        Int_t dBucket = TMath::Abs(maxValueBucket-maxValueBucketdec);
        if(maxValue>1000 && dBucket>2 && dBucket<100){
            //GetEnergybyFitWaveform(rftype,decayIdx,cobo,asad,aget,chan,maxValueBucket); //using ratio method
            //GetEnergybyFitWaveform(rftype,decayIdx,cobo,asad,aget,chan,maxValueBucket); //using ShaperF method
            //GetEnergybyFitWaveform(2,decayIdx,cobo,asad,aget,chan,maxValueBucket); //using ShaperF_Fusion method
        }else{
            //GetEnergybyFitWaveform(rftype,decayIdx,cobo,asad,aget,chan,maxValueBucket);
            //GetEnergybyFitWaveform(rftype,decayIdx,cobo,asad,aget,chan,maxValueBucket); //using ShaperF method
            //GetEnergybyFitWaveform(2,decayIdx,cobo,asad,aget,chan,maxValueBucket); //using ShaperF_Fusion method
            //maxValuefit = maxValuedec;
        }
        //cout << cobo << " " << asad << " " << aget << " " << chan << " " << maxValue << " " << maxValuedec << " " << maxValuefit << " " << dBucket << endl;

        if(enablehist==1){
            hGET_DecEHitPattern2D->Fill(cobo*2000+asad*500+aget*100+chan,maxValuedec-maxValue);
            hGET_DecTHitPattern2D->Fill(cobo*2000+asad*500+aget*100+chan,maxValueBucketdec-maxValueBucket);
            hGET_FitEHitPattern2D->Fill(cobo*2000+asad*500+aget*100+chan,maxValuefit-maxValuedec);
            hGET_FitTHitPattern2D->Fill(cobo*2000+asad*500+aget*100+chan,maxValueBucketfit-maxValueBucketdec);
            hGET_DecdEMaxvsEMax2D->Fill(maxValue,maxValue-maxValuedec);
            hGET_DecdTMaxvsTMax2D->Fill(maxValueBucket,maxValueBucket-maxValueBucketdec);
            hGET_FitdEMaxvsEMax2D->Fill(maxValue,maxValue-maxValuefit);
            hGET_FitdTMaxvsTMax2D->Fill(maxValueBucket,maxValueBucket-maxValueBucketfit);
        }
        if(energymethod==1){
            if(maxValue<maxValuedec){
                maxValue = maxValuedec;
                maxValueBucket = maxValueBucketdec;
            }
        }else if(energymethod==2){
            if(maxValue<maxValuefit){
                maxValue = maxValuefit;
                maxValueBucket = maxValueBucketfit;
            }
        }
    }

    corrwaveformintegral=0;
    maxValue = maxValue-baseline;
    maxValueBucket = maxValueBucket;
    for(UInt_t buck=0;buck<bucketmax;buck++){
        corrwaveformintegral+=(rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck]-baseline);
    }
    psdratio = ((Double_t)maxValue/(Double_t)corrwaveformintegral);

    //if(maxValue>mm_minenergy && maxValue<mm_maxenergy && maxValueBucket>mintime && maxValueBucket<maxtime)
    if((cobo==0&&maxValue>mm_minenergy && maxValue<mm_maxenergy && maxValueBucket>mintime && maxValueBucket<maxtime)||cobo>0)
    {
        rwaveforms[decayIdx][cobo]->hasHit[asad*4+aget] = true;
        rwaveforms[decayIdx][cobo]->hasSignal[asad*4+aget][chan] = true;
        rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan] = maxValue;
        if(decayIdx>0) cout<<"DECAY MODE TIME"<<maxValueBucket<<"\t"<<decayIdx<<endl;
        rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan] = maxValueBucket;// - decayIdx*256;
        rwaveforms[decayIdx][cobo]->baseline[asad*4+aget][chan] = baseline;
        rwaveforms[decayIdx][cobo]->PSDIntegral[asad*4+aget][chan] = corrwaveformintegral;
        rwaveforms[decayIdx][cobo]->PSDRatio[asad*4+aget][chan] = psdratio;
        if(enable2pmode==1){
            if((rwaveforms[decayIdx][cobo]->decayIdx)==0){
                rwaveforms[decayIdx][cobo]->hasImplant = true;
                //        rwaveforms[decayIdx][cobo]->isDecay[asad*4+aget][chan] = 0;
            }else{
                rwaveforms[decayIdx][cobo]->hasDecay = true;
                //        rwaveforms[decayIdx][cobo]->isDecay[asad*4+aget][chan] = 1;
            }

        }
        if(cobo==1 && asad==1 && aget==1 && chan==2) {
            rwaveforms[decayIdx][cobo]->ICenergy = rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan];
        }
    }
}

/*
Int_t LKFrameBuilder::GetEnergybyFitWaveform(Int_t type, Int_t decayIdx, Int_t cobo, Int_t asad, Int_t aget, Int_t chan, Int_t maxAt, Int_t peakAt, Int_t *maxValuedec, Int_t *maxvalueBucketdec){
  Int_t bucketmin = peakAt - maxresponse[type];
  Int_t bucketgap = peakAt - maxAt;
  Int_t bucketmax = peakAt + responsesigma[type]*1.5;

  Double_t ratio = rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][maxAt-1] / response[type][maxresponse[type]-(peakAt-maxAt)-1];
  Double_t chisq=0;
  Double_t ndof=0;
  Double_t minchisq=1E9;
  Double_t tryratio=0;
  Double_t bestratio=ratio;
  Int_t maxstep=40;
  Double_t stepsize=0.05;
  for(int i=0;i<maxstep;i++){
    chisq = 0;
    for(int buck=0;buck<bucketmax;buck++){
      if((buck+bucketmin)<bucketmax){
        if(buck<maxresponse[type]-bucketgap || buck>maxresponse[type]+bucketgap){
          tryratio = ratio*(0.5+i*stepsize);
          chisq += TMath::Power(rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck+bucketmin] - response[type][buck]*tryratio,2);
          ndof += rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck+bucketmin];
        }
      }
    }

    chisq = chisq/ndof;
    if(minchisq>chisq){
      minchisq = chisq;
      bestratio = tryratio;
    }
    if(reventIdx==7861 && cobo==0 && asad==2 && aget==1 && chan==42){
      cout << i << " " << tryratio << " " << chisq << " " << minchisq << " " << bestratio << endl;
    }
    //delete hScaledResponse;
  }
  for(int buck=0;buck<bucketmax;buck++){
    if((buck+bucketmin)<bucketmax){
      corrwaveformfit[buck+bucketmin] = response[type][buck]*bestratio;
    }
  }
  if(response[type][maxresponse[type]]*bestratio>0 && response[type][maxresponse[type]]*bestratio>100000){
    *maxValuedec = (Int_t)response[type][maxresponse[type]]*bestratio;
  }else{
    *maxValuedec = 0;
  }
  //hWaveform->SetTitle(Form("minchisq=%f, bestratio=%f, maxValuedec=%d",minchisq,bestratio,*maxValuedec));
  //delete hWaveform;
}
*/

void LKFrameBuilder::GetEnergybyFitWaveform(Int_t type, Int_t decayIdx, Int_t cobo, Int_t asad, Int_t aget, Int_t chan, Int_t peakAt)
{
    Double_t tau, power;
    TF1* fWaveform;
    TH1D* hWaveform;
    Int_t bucketmin = peakAt - responsesigma[type]*1.5;
    Int_t bucketmax = peakAt + responsesigma[type]*1.5;
    Int_t baseline = rwaveforms[decayIdx][cobo]->baseline[asad*4+aget][chan];
    Int_t saturated = 0;
    Int_t minbuck = 0;

    if(type<3){
        hWaveform = new TH1D("hWaveform","hWaveform",512,0,512);
        for(int buck=0;buck<512;buck++){
            hWaveform->Fill(buck,rwaveforms[decayIdx][cobo]->corrwaveform[asad*4+aget][chan][buck]);
            if(rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck]==0){
                hWaveform->SetBinError(buck,0.0);
                if(saturated==0) minbuck = buck;
                saturated++;
            }else if(rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][buck]==4095){
                hWaveform->SetBinError(buck,0.0);
                if(saturated==0) minbuck = buck;
                saturated++;
            }
        }

        if(type==0){
            //fWaveform = new TF1("fWaveform", this, &LKFrameBuilder::ShaperF_GET1, bucketmin, bucketmax, 6);
            fWaveform = new TF1("fWaveform", this, &LKFrameBuilder::ShaperF_GET1, 0, 512, 6);
            fWaveform->SetParNames("offset", "amplitude", "peakAt", "sigma", "power", "p2");
            fWaveform->SetParameters(baseline, 500, peakAt, responsesigma[type], 3, 0.2);
            if(saturated>0){
                peakAt = minbuck + saturated/2;
                fWaveform->SetParameter(2, peakAt);
            }
            fWaveform->FixParameter(3, responsesigma[type]);
            fWaveform->SetParLimits(1, 0, 10000);
        }else if(type==1) {
            fWaveform = new TF1("fWaveform", this, &LKFrameBuilder::ShaperF_MSCF, bucketmin, bucketmax, 6);
            fWaveform->SetParNames("offset", "amplitude", "peakAt", "sigma", "power", "p2");
            fWaveform->SetParameters(baseline, 2000, peakAt, responsesigma[type], 3, 0.2);
            if(saturated>0){
                peakAt = minbuck + saturated/2;
                fWaveform->SetParameter(2, peakAt);
            }
            fWaveform->FixParameter(3, responsesigma[type]);
            fWaveform->SetParLimits(1, 0, 10000);
        }else if(type==2){
            //fWaveform = new TF1("fWaveform", this, &LKFrameBuilder::ShaperF_GET2, bucketmin, bucketmax, 5);
            fWaveform = new TF1("fWaveform", this, &LKFrameBuilder::ShaperF_GET2, 0, 512, 5);
            fWaveform->SetParNames("offset", "amplitude", "peakAt", "sigma", "power");
            if(cobo==0){
                tau = 4.5; power = 2.3e-3;
            }else if(cobo==0){
                tau = 15.5; power = 3.2e-4;
            }
            fWaveform->SetParameters(baseline, 500, peakAt, tau, power);
            if(saturated>0){
                peakAt = minbuck + saturated/2;
                fWaveform->SetParameter(2, peakAt);
            }
            fWaveform->SetParLimits(1, 0, 10000);
            fWaveform->SetParLimits(3, tau*0.96, tau*1.04);
            fWaveform->SetParLimits(4, power*0.96, power*1.1);
        }

        hWaveform->Fit("fWaveform","WWQN+"); //Q=Quiet and N=Do not draw
        maxValuefit = fWaveform->GetMaximum();
        maxValueBucketfit = fWaveform->GetMaximumX();
        for(int buck=0;buck<512;buck++){
            corrwaveformfit[buck] = fWaveform->Eval(buck);
            //if(reventIdx==104970&&cobo==0&&asad==0&&aget==0&&chan==64) cout << cobo << " " << asad << " " << aget << " " << chan << " " << corrwaveformfit[buck] << " " << maxValuefit << " " << maxValueBucketfit << " " << maxValueBucket << endl;
        }
        delete hWaveform;
        delete fWaveform;
    }else{
        cout << Form("No wavefunction type defined: %d",type) << endl;
        maxValuefit = 0;
        maxValueBucketfit = 0;
    }
}

void LKFrameBuilder::ResetHitPattern() {
    hGET_THitPattern[goodevtcounter%16]->Reset();
    hGET_EHitPattern[goodevtcounter%16]->Reset();
    hGET_ERHitPattern[goodevtcounter%16]->Reset();
    /*
       for(int i=0;i<2000;i++){
       for(int j=0;j<512;j++){
       hGET_THitPattern[goodevtcounter%16]->SetBinContent(i,j,0);
       hGET_EHitPattern[goodevtcounter%16]->SetBinContent(i,j,0);
       }
       for(int j=0;j<100;j++){
       hGET_ERHitPattern[goodevtcounter%16]->SetBinContent(i,j,0);
       }
       }
     */
}

void LKFrameBuilder::DrawHitPattern(Int_t decayIdx, UInt_t cobo)
{
    Double_t psdratio = 0;
    std::vector<Int_t> si16x16front;
    Int_t frontid=0;
    Int_t dchan=0;
    std::vector<Int_t> si16x16back;
    std::vector<Int_t> sifront;
    std::vector<Int_t> siback;
    std::vector<Int_t> csifront;
    std::vector<Int_t> sibackenergy;
    std::vector<Int_t> csifrontenergy;
    Int_t backid=0;
    Int_t sifwmult=0;
    for(UInt_t asad=0; asad<maxasad; asad++) {
        for(UInt_t aget=0; aget<4; aget++) {

            if(rwaveforms[decayIdx][cobo]->hasHit[asad*4+aget]){
                for(UInt_t chan=0; chan<68; chan++) {

                    if(!rwaveforms[decayIdx][cobo]->hasSignal[asad*4+aget][chan]) continue; // skip signals that did not fire

                    if(rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]>4095) rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan] = 4095;
                    if(readmode==2){
                        hGET_EHitPattern[goodevtcounter%16]->Fill(asad*4*100+aget*100+chan,rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);
                        hGET_THitPattern[goodevtcounter%16]->Fill(asad*4*100+aget*100+chan,rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan]);
                        psdratio = rwaveforms[decayIdx][cobo]->PSDRatio[asad*4+aget][chan];
                        hGET_ERHitPattern[goodevtcounter%16]->Fill(asad*4*100+aget*100+chan,psdratio);
                    }
                    hGET_E[cobo*maxasad*4+asad*4+aget][chan]->Fill(rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);
                    hGET_EALL[cobo*maxasad*4+asad*4+aget]->Fill(rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);
                    if(readmode==4){
                        hGET_EHitPattern2D->Fill(cobo*2000+asad*500+aget*100+chan,rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);
                        hGET_THitPattern2D->Fill(cobo*2000+asad*500+aget*100+chan,rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan]);
                    }
                    if(cobo==1&&asad==0){
                        if(chan<11) {
                            dchan = chan;
                        }else if(chan>11 && chan<22) {
                            dchan = chan - 1;
                        }else if(chan>22 && chan<45) {
                            dchan = chan - 2;
                        }else if(chan>45 && chan<56) {
                            dchan = chan - 3;
                        }else if(chan>56) {
                            dchan = chan - 4;
                        }

                        if(aget==0){
                            //hGET_SiForwardHitPattern2D->Fill(mapchantosi->pxidx[asad][aget][chan],mapchantosi->pxidy[asad][aget][chan]);
                            sifwmult++;
                        }

                        //if((aget==0||aget==2)&&rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]>0)
                        if((aget==0)&&rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]>0)
                        {
                            sifront.push_back(aget*100+chan);
                        }
                        //if((aget==1||aget==3)&&rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]>350)
                        //if(aget==1||aget==3)
                        if(aget==1)
                        {
                            siback.push_back(aget*100+chan);
                            sibackenergy.push_back(rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);
                            //              cout << "si_back:" << asad << " " << aget << " " << rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan] << endl;
                        }
                        if(rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]>0){
                            if(aget==2){
                                frontid = dchan/2;
                                si16x16front.push_back(frontid);
                            }else if(aget==3){
                                backid = dchan/2;
                                if(backid>7) backid = 8 - backid + 15;
                                else backid = backid;
                                //if(backid<8) backid += 8;
                                //else backid = 15 - backid;
                                //if(backid>8) backid = (15-backid)+8;
                                si16x16back.push_back(backid);
                            }
                        }
                    }
                    if(cobo==1&&asad==1){
                        if((aget==0)&&rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]>400){
                            csifront.push_back(aget*100+chan);
                            csifrontenergy.push_back(rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);
                            //              cout << reventIdx<<" csi:" << asad << " " << aget << " " << rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan] << endl;
                        }
                    }
                }
                if(readmode==2){
                    hGET_EHitPattern[goodevtcounter%16]->SetTitle(Form("hGET_EHitPattern_%d(EvtNo=%d);asad*400+aget*100+chan;ADC Channel (ch)",goodevtcounter%16,reventIdx));
                    hGET_THitPattern[goodevtcounter%16]->SetTitle(Form("hGET_THitPattern_%d(EvtNo=%d);asad*400+aget*100+chan;ADC Channel (ch)",goodevtcounter%16,reventIdx));
                    hGET_ERHitPattern[goodevtcounter%16]->SetTitle(Form("hGET_ERHitPattern_%d(EvtNo=%d);asad*400+aget*100+chan;Ratio",goodevtcounter%16,reventIdx));
                }
            }
        }
    }

    if(sifwmult>0){
        hGET_SiForwardMult->Fill(sifwmult);
    }

    for(int i=0;i<si16x16front.size();i++){
        for(int j=0;j<si16x16back.size();j++){
            hGET_Si16x16HitPattern2D->Fill(si16x16front.at(i),si16x16back.at(j));
            hGET_SiForwardHitPattern2D->Fill((si16x16front.at(i)/8)+4,(si16x16back.at(j)/8));
        }
    }

    for(int i=0;i<sifront.size();i++){
        for(int j=0;j<csifront.size();j++){
            hGET_SivsCsIHitPattern2D->Fill(sifront.at(i),csifront.at(j)); //front
        }
    }
    //cout<<"Siback size = "<<siback.size()<<"\tCsI size: "<<csifront.size()<<endl;
    if(siback.size()>0 && csifront.size()>0){
        for(int i=0;i<siback.size();i++){
            for(int j=0;j<csifront.size();j++){
                csi_map[2]=1;
                csi_map[7]=2;
                csi_map[10]=3;
                csi_map[16]=4;
                csi_map[19]=5;
                csi_map[25]=6;
                csi_map[28]=7;
                csi_map[33]=8;
                csi_map[36]=9;
                csi_map[41]=10;
                if(siback.at(i)<105 && csifront.at(j)<264) hMM_SiEvsCsIE[siback.at(i)-100][csi_map[csifront.at(j)-200]]->Fill(csifrontenergy.at(j),sibackenergy.at(i));
                if(siback.at(i)>133 && siback.at(i)<138 && csifront.at(j)<264) hMM_SiEvsCsIE[siback.at(i)-129][csi_map[csifront.at(j)-200]]->Fill(csifrontenergy.at(j),sibackenergy.at(i));
                //if(siback.at(i)==104 && csifront.at(j)==233) hMM_SiEvsCsIE[1]->Fill(csifrontenergy.at(j),sibackenergy.at(i));
                //if(reventIdx==670644) cout << siback.at(i) << " " << sibackenergy.at(i) << " " << csifront.at(j) << " " << csifrontenergy.at(j) << endl;

                //if(sibackenergy.at(i)>100 && sibackenergy.at(i)<2500 && csifrontenergy.at(j)>200 && cut_pinSiEvsCsIE->IsInside(csifrontenergy.at(j),sibackenergy.at(i)))
                if(sibackenergy.at(i)>100 && sibackenergy.at(i)<2500 && csifrontenergy.at(j)>200)
                {
                    //if((siback.at(i)==100 && csifront.at(j)==216) ||
                    //        (siback.at(i)==101 && csifront.at(j)==210) ||
                    //        (siback.at(i)==102 && csifront.at(j)==219) ||
                    //        (siback.at(i)==103 && csifront.at(j)==225) ||
                    //        (siback.at(i)==104 && csifront.at(j)==233) ||
                    //        (siback.at(i)==134 && csifront.at(j)==202) ||
                    //        (siback.at(i)==135 && csifront.at(j)==207) ||
                    //        (siback.at(i)==136 && csifront.at(j)==236) ||
                    //        (siback.at(i)==137 && csifront.at(j)==228))
                    if(sibackenergy.at(i)>100 && sibackenergy.at(i)<2500 && csifrontenergy.at(j)>200){
                        hGET_SivsCsIHitPattern2D->Fill(siback.at(i),csifront.at(j)); //back
                    }
                    //cout << reventIdx << " " << siback.at(i) << " " << sibackenergy.at(i) << " " << csifront.at(j) << " " << csifrontenergy.at(j) << endl;

                    hMM_SiEvsCsIEAll->Fill(csifrontenergy.at(j),sibackenergy.at(i));
                    if(cut_pinSiEvsCsIE->IsInside(csifrontenergy.at(j),sibackenergy.at(i))) si_tracks->hasFWC++;
                }
            }
        }
        siback.clear();
        csifront.clear();
        csifrontenergy.clear();
        sibackenergy.clear();
    }
}

void LKFrameBuilder::WaveletFilter(Int_t decayIdx, UInt_t cobo)
{
    std::vector<Double_t> widths;
    widths.push_back(8);
    //widths.push_back(32); //for 512 timebucket
    widths.push_back(4); //for 256 timebucket
    WaveletNew* wavelet;

    for(UInt_t asad=0; asad<maxasad; asad++) {
        for(UInt_t aget=0; aget<4; aget++) {
            if(!rwaveforms[decayIdx][cobo]->hasHit[asad*4+aget]) continue; // skip aget that did not fire
            for(UInt_t chan=0; chan<68; chan++) {
                if(chan==11 || chan==22 || chan==45 || chan==56) continue; // We want to skip the FPN channels.
                if(!rwaveforms[decayIdx][cobo]->hasSignal[asad*4+aget][chan]) continue; // skip signals that did not fire
                //if(cobo>0)
                if(cobo==0 && rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]>mm_minenergy && rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]<mm_maxenergy)
                {
                    std::vector<double> wfinput;
                    //std::vector<double> wfinput(rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan].begin(),rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan].end());

                    Int_t decayIdx = rwaveforms[decayIdx][cobo]->decayIdx;
                    for(int i=0;i<bucketmax;i++){
                        //for(int i=0;i<bucketmax;i++)
                        wfinput.push_back(static_cast<double>(rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][i])/rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);
                        // wfinput[i] /= rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan];
                    }
                    WaveletNew* wavelet = new WaveletNew(wfinput, widths, true);
                    wavelet->CalcCWTFast();
                    std::vector<Double_t> scale = wavelet->GetScale();
                    std::vector<Double_t> pa = wavelet->GetPa();
                    if(enablehist==1){
                        hMM_Pa_vs_Energy2D->Fill(rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan],pa[1]/pa[0]);
                        hMM_Pa_vs_Time2D->Fill(rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan],pa[1]/pa[0]);
                    }
                    delete wavelet;
                }else if(cobo==1 && rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]>si_minenergy && rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]<si_maxenergy){
                    std::vector<double> wfinput;
                    //std::vector<double> wfinput(rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan].begin(),rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan].end());

                    Int_t decayIdx = rwaveforms[decayIdx][cobo]->decayIdx;
                    for(int i=0;i<bucketmax;i++){
                        //for(int i=0;i<bucketmax;i++)
                        wfinput.push_back(static_cast<double>(rwaveforms[decayIdx][cobo]->waveform[asad*4+aget][chan][i])/rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);
                        // wfinput[i] /= rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan];
                    }
                    WaveletNew* wavelet = new WaveletNew(wfinput, widths, true);
                    wavelet->CalcCWTFast();
                    std::vector<Double_t> scale = wavelet->GetScale();
                    std::vector<Double_t> pa = wavelet->GetPa();
                    if(enablehist==1){
                        hSi_Pa_vs_Energy2D->Fill(rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan],pa[1]/pa[0]);
                        hSi_Pa_vs_Time2D->Fill(rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan],pa[1]/pa[0]);
                    }
                    delete wavelet;
                }
            }
        }
    }
}

void LKFrameBuilder::WaveformShapeFilter(Int_t decayIdx, UInt_t cobo)
{
    Int_t maxValue = -100000;
    UInt_t maxValueBucket = 0;
    Int_t projybegin,projyend;
    Int_t eentriesfh,emeanfh,elimitfh,eintegralfh;
    Int_t tentriesfh,tmeanfh,tlimitfh,tintegralfh;
    Int_t eentriessh,emeansh,elimitsh,eintegralsh;
    Int_t tentriessh,tmeansh,tlimitsh,tintegralsh;
    Int_t bin1,bin2;
    ostringstream s1out;
    ostringstream s2out;
    ostringstream s3out;
    ofstream f1out;
    ofstream f2out;
    ofstream f3out;
    for(UInt_t asad=0; asad<maxasad; asad++) {
        for(UInt_t aget=0; aget<4; aget++) {
            if(!rwaveforms[decayIdx][cobo]->hasHit[asad*4+aget]) continue; // skip aget that did not fire
            if((asad==2||asad==3)&&(aget==0||aget==1)){
                projybegin=asad*4*100+aget*100;
                projyend=asad*4*100+aget*100+32;
                eentriesfh = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetEntries();
                emeanfh = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetMean();
                elimitfh = 50;
                bin1 = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetXaxis()->FindBin(emeanfh-elimitfh);
                bin2 = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetXaxis()->FindBin(emeanfh+elimitfh);
                eintegralfh = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->Integral(bin1,bin2);
                tentriesfh = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetEntries();
                tmeanfh = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetMean();
                tlimitfh = 2;
                bin1 = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetXaxis()->FindBin(tmeanfh-tlimitfh);
                bin2 = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetXaxis()->FindBin(tmeanfh+tlimitfh);
                tintegralfh = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->Integral(bin1,bin2);
                projybegin=asad*4*100+aget*100+33;
                projyend=asad*4*100+aget*100+67;
                eentriessh = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetEntries();
                emeansh = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetMean();
                elimitsh = 50;
                bin1 = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetXaxis()->FindBin(emeansh-elimitsh);
                bin2 = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetXaxis()->FindBin(emeansh+elimitsh);
                eintegralsh = hGET_EHitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->Integral(bin1,bin2);
                tentriessh = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetEntries();
                tmeansh = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetMean();
                tlimitsh = 2;
                bin1 = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetXaxis()->FindBin(tmeansh-tlimitsh);
                bin2 = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->GetXaxis()->FindBin(tmeansh+tlimitsh);
                tintegralsh = hGET_THitPattern[goodevtcounter%16]->ProjectionY("projy",projybegin,projyend)->Integral(bin1,bin2);
            }
            for(UInt_t chan=0; chan<68; chan++) {
                if(chan==11 || chan==22 || chan==45 || chan==56) continue; // We want to skip the FPN channels.
                if(!rwaveforms[decayIdx][cobo]->hasSignal[asad*4+aget][chan]) continue; // skip signals that did not fire

                maxValue = rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan];
                maxValueBucket = rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan];
                if((asad==2||asad==3)&&(aget==0||aget==1) &&
                        (!rwaveforms[decayIdx][cobo]->isOverflow[asad*4+aget][chan])){
                    if(chan<33 && eintegralfh*3>eentriesfh && tintegralfh*3>tentriesfh){
                        if( maxValue>=(emeanfh-elimitfh) &&
                                maxValue<=(emeanfh+elimitfh) &&
                                maxValueBucket>=(tmeanfh-tlimitfh) &&
                                maxValueBucket<=(tmeanfh+tlimitfh)){
                            //                rwaveforms[decayIdx][cobo]->hasSignal[asad*4+aget][chan] = 0;
                            //continue; // skip fake signals by saturated signal
                            //s2out << reventIdx << " " << cobo << " " << asad << " " << aget << " " << chan << " 2" << endl;
                        }
                    }else if(chan>32 && eintegralsh*3>eentriessh && tintegralsh*3>tentriessh){
                        if( maxValue>=(emeansh-elimitsh) &&
                                maxValue<=(emeansh+elimitsh) &&
                                maxValueBucket>=(tmeansh-tlimitsh) &&
                                maxValueBucket<=(tmeansh+tlimitsh)){
                            //                rwaveforms[decayIdx][cobo]->hasSignal[asad*4+aget][chan] = 0;
                            //continue; // skip fake signals by saturated signal
                            //s2out << reventIdx << " " << cobo << " " << asad << " " << aget << " " << chan << " 2" << endl;
                        }
                    }
                }else if((cobo==0)&&(asad==0||asad==1)&&(!rwaveforms[decayIdx][cobo]->isOverflow[asad*4+aget][chan])){
                    //s1out << reventIdx << " " << cobo << " " << asad << " " << aget << " " << chan << " 1" << endl;
                }
            }
        }
    }
    /*
       if(s1out.str().size()!=0){
       f1out.open("goodlabels.txt", std::ofstream::out|std::ofstream::app);
       f1out << s1out.str() << endl;
       f1out.close();
       s1out.clear();
       s1out.str("");
       }
       if(s2out.str().size()!=0){
       f2out.open("badlabels.txt", std::ofstream::out|std::ofstream::app);
       f2out << s2out.str() << endl;
       f2out.close();
       s2out.clear();
       s2out.str("");
       }
     */
}

void LKFrameBuilder::DrawPSDFilter(Int_t decayIdx, UInt_t cobo)
{
    Double_t ratiocut=0.0397+0.00412; // centroid + 1sigma for good data
    Int_t maxValue = -100000;
    Int_t corrwaveformintegral = 0;
    Double_t psdratio = 0;
    Int_t baseline = 0;
    for(UInt_t asad=0; asad<maxasad; asad++) {
        for(UInt_t aget=0; aget<4; aget++) {
            if(!rwaveforms[decayIdx][cobo]->hasHit[asad*4+aget]) continue; // skip aget that did not fire
            for(UInt_t chan=0; chan<68; chan++) {
                if(chan==11 || chan==22 || chan==45 || chan==56) continue; // We want to skip the FPN channels.
                if(!rwaveforms[decayIdx][cobo]->hasSignal[asad*4+aget][chan]) continue; // skip signals that did not fire

                if((asad==2||asad==3)&&(aget==0||aget==1)&&(rwaveforms[decayIdx][cobo]->hasOverflow)){
                    maxValue = rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan];
                    corrwaveformintegral = rwaveforms[decayIdx][cobo]->PSDIntegral[asad*4+aget][chan];
                    psdratio = rwaveforms[decayIdx][cobo]->PSDRatio[asad*4+aget][chan];
                    //((maxValue<300 && psdratio>0.05) || (maxValue>300 && psdratio>ratiocut)))
                    if(rwaveforms[decayIdx][cobo]->isOverflow[asad*4+aget][chan])
                    {
                        hMM_PSDIntegralvsMax[0]->Fill(maxValue,corrwaveformintegral);
                        hMM_PSDRatiovsMax[0]->Fill(maxValue,psdratio);
                        hMM_PSDRatio[0]->Fill(0.0,psdratio);
                    }else{
                        hMM_PSDIntegralvsMax[1]->Fill(maxValue,corrwaveformintegral);
                        hMM_PSDRatiovsMax[1]->Fill(maxValue,psdratio);
                        hMM_PSDRatio[0]->Fill(1.0,psdratio);
                    }
                }else if((asad==2||asad==3)&&(aget==0||aget==1)){
                    maxValue = rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan];
                    baseline = rwaveforms[decayIdx][cobo]->baseline[asad*4+aget][chan];
                    corrwaveformintegral = rwaveforms[decayIdx][cobo]->PSDIntegral[asad*4+aget][chan];
                    psdratio = rwaveforms[decayIdx][cobo]->PSDRatio[asad*4+aget][chan];
                    hMM_PSDIntegralvsMax[2]->Fill(maxValue,corrwaveformintegral);
                    hMM_PSDRatiovsMax[2]->Fill(maxValue,psdratio);
                    hMM_PSDRatio[0]->Fill(2.0,psdratio);
                }
            }
        }
    }
}

void LKFrameBuilder::FillTrack()
{
    UInt_t dchan = 0;
    UInt_t spxidx = 0;
    UInt_t spxidy = 0;
    int rgidx=0;
    Int_t decayIdxMax=1;
    if(enable2pmode==1) decayIdxMax=2;

    for(Int_t decayIdx=0; decayIdx<decayIdxMax; decayIdx++) {
        for(UInt_t cobo=0; cobo<no_cobos; cobo++) {
            if(cobo==1){
                mm_tracks->ICenergy=rwaveforms[decayIdx][cobo]->ICenergy;
                mm_tracks->Sienergy=rwaveforms[decayIdx][cobo]->Sienergy;
                mm_tracks->SiX=rwaveforms[decayIdx][cobo]->SiX;
                mm_tracks->SiY=rwaveforms[decayIdx][cobo]->SiY;
                mm_tracks->SiZ=rwaveforms[decayIdx][cobo]->SiZ;
            }
            for(UInt_t asad=0; asad<maxasad; asad++) {
                for(UInt_t aget=0; aget<4; aget++) {
                    if(rwaveforms[decayIdx][cobo]->isRejected) continue; // skip a bad event from saturated Si det signals
                    if(!rwaveforms[decayIdx][cobo]->hasHit[asad*4+aget]) continue; // skip agets that did not have any fire
                    for(UInt_t chan=0; chan<68; chan++) {
                        if(chan==11 || chan==22 || chan==45 || chan==56) continue; // We want to skip the FPN channels.
                        if(!rwaveforms[decayIdx][cobo]->hasSignal[asad*4+aget][chan]) continue; // skip signals that did not fire

                        //cout << "Good ones: " << reventIdx << " " << cobo << " " << asad << " " << aget << " " << chan << endl;
                        if(cobo==0
                                && rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]>mm_minenergy
                                && rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]<mm_maxenergy
                                && rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan]>mm_mintime
                                && rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan]<mm_maxtime)
                        {
                            if(chan<11) {
                                dchan = chan;
                            }else if(chan>11 && chan<22) {
                                dchan = chan - 1;
                            }else if(chan>22 && chan<45) {
                                dchan = chan - 2;
                            }else if(chan>45 && chan<56) {
                                dchan = chan - 3;
                            }else if(chan>56) {
                                dchan = chan - 4;
                            }

                            //if((asad==2||asad==3)&&(aget==0||aget==1)&&(rwaveforms[decayIdx][cobo]->hasOverflow[asad*4+aget]))
                            if(rwaveforms[decayIdx][cobo]->hasOverflow)
                            {
                                mm_tracks->hasOverflow=true;
                            }
                            //if(cobo==0
                            //        && asad==1
                            //        && aget==0
                            //        && rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan]>220
                            //        && rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan]<280) mm_tracks->hasTrack+=1;
                            mm_tracks->hasTrack+=1;
                            spxidx = mapchantomm->pxidx[asad][aget][dchan];
                            spxidy = mapchantomm->pxidy[asad][aget][dchan];

                            //if(rwaveforms[decayIdx][cobo]->decayIdx==1 || 1)
                            //  cout<<rwaveforms[decayIdx][cobo]->decayIdx<<"\t"<<spxidx<<"\t"<<spxidy<<"\t"
                            //  <<rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan]<<endl;

                            if(cobo==0 && spxidx!=-1 && spxidy!=-1) {//JEB
                                mm_tracks_2p->pixel[spxidx][spxidy]+=1;
                                if(rwaveforms[decayIdx][cobo]->isDecay[asad*4+aget][chan]==0)
                                    hMM_TrackDecay[goodevtcounter%16]->Fill(spxidx,spxidy,rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);
                                if(rwaveforms[decayIdx][cobo]->isDecay[asad*4+aget][chan]==1)
                                    hMM_TrackDecay2[goodevtcounter%16]->Fill(spxidx,spxidy,rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);
                                //if(rwaveforms[decayIdx][cobo]->decayIdx==1)
                                //  hMM_TrackDecay2[goodevtcounter%16]->Fill(spxidx,spxidy,rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);
                            }
                            if(cobo==0 && spxidx!=-1 && spxidy==-1) {//JEB chains
                                for(int i=0;i<128;i++) {
                                    if(rwaveforms[decayIdx][cobo]->isDecay[asad*4+aget][chan]==0)
                                        hMM_TrackDecay[goodevtcounter%16]->Fill(spxidx,2*i+1,rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);
                                    if(rwaveforms[decayIdx][cobo]->isDecay[asad*4+aget][chan]==1)
                                        hMM_TrackDecay2[goodevtcounter%16]->Fill(spxidx,2*i+1,rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);

                                }
                            }

                            if(cobo==0 && spxidx==-1 && spxidy!=-1) {//JEB strips
                                if(rwaveforms[decayIdx][cobo]->isDecay[asad*4+aget][chan]==0)
                                    hMM_TrackDecay[goodevtcounter%16]->Fill(spxidx,spxidy,rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);
                                if(rwaveforms[decayIdx][cobo]->isDecay[asad*4+aget][chan]==1)
                                    hMM_TrackDecay2[goodevtcounter%16]->Fill(spxidx,spxidy,rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);
                                if(asad==2) {
                                    for(int i=0;i<63;i++) {
                                        if(rwaveforms[decayIdx][cobo]->decayIdx==0)
                                            hMM_TrackDecay[goodevtcounter%16]->Fill(i,spxidy,rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);
                                        if(rwaveforms[decayIdx][cobo]->decayIdx==1)
                                            hMM_TrackDecay2[goodevtcounter%16]->Fill(i,spxidy,rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);
                                    }
                                }
                                if(asad==3) {
                                    for(int i=70;i<140;i++) {
                                        if(rwaveforms[decayIdx][cobo]->decayIdx==0)
                                            hMM_TrackDecay[goodevtcounter%16]->Fill(i,spxidy,rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);
                                        if(rwaveforms[decayIdx][cobo]->decayIdx==1)
                                            hMM_TrackDecay2[goodevtcounter%16]->Fill(i,spxidy,rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]);
                                    }
                                }
                            }
                            if(spxidx!=-1 && spxidy!=-1){
                                mm_tracks->pixel[spxidx][spxidy]+=decayIdx*100-1*(decayIdx-1);
                                mm_tracks->time[spxidx][spxidy]=rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan];
                                mm_tracks->energy[spxidx][spxidy]=rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan];
                                mm_tracks->posx[spxidx][spxidy]=spxidx*3.4-67*3.4 + 3.4/2;
                                mm_tracks->posy[spxidx][spxidy]=spxidy*1.7 + 1.7/2;
                                mm_tracks->posxerr[spxidx][spxidy]=3.4/2;
                                mm_tracks->posyerr[spxidx][spxidy]=1.7/2;
                                mm_tracks->poszerr[spxidx][spxidy]=fTimePerBin*driftv/2;
                                rx=rpos->Uniform(-1,1)*mm_tracks->posxerr[spxidx][spxidy];
                                ry=rpos->Uniform(-1,1)*mm_tracks->posyerr[spxidx][spxidy];
                                rz=rpos->Uniform(-1,1)*mm_tracks->poszerr[spxidx][spxidy];
                                //if(spxidx>63 && spxidx<68)
                                //cout << spxidx << " " << spxidy << ", posx=" <<  mm_tracks->posx[spxidx][spxidy]
                                //<< ", posy=" << mm_tracks->posy[spxidx][spxidy] << endl;
                                if(asad<2) mm_tracks->coloridx[spxidx][spxidy]=1;
                                else mm_tracks->coloridx[spxidx][spxidy]=2;

                                if(mm_tracks->energy[spxidx][spxidy]>0){
                                    posenergysum[spxidy]+=(mm_tracks->posx[spxidx][spxidy]+rx)*(mm_tracks->energy[spxidx][spxidy]);
                                    energysum[spxidy]+=mm_tracks->energy[spxidx][spxidy];
                                    sumcounter[spxidy]++;
                                }
                            }
                            if(spxidx==-1){ //strip
                                mm_tracks->hasTrackStrip+=1;
                                if(asad==2){
                                    for(int i=0;i<64;i++){
                                        mm_tracks->pixel[i][spxidy]+=decayIdx*100-1*(decayIdx-1);
                                        mm_tracks->time[i][spxidy]=rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan];
                                        mm_tracks->energy[i][spxidy]=rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan];
                                        mm_tracks->posx[i][spxidy]=i*1.7-64*1.7-3*3.4 + 1.7/2;
                                        mm_tracks->posy[i][spxidy]=spxidy*1.7 + 1.7/2;
                                        mm_tracks->posxerr[i][spxidy]=1.7/2;
                                        mm_tracks->posyerr[i][spxidy]=1.7/2;
                                        mm_tracks->poszerr[i][spxidy]=fTimePerBin*driftv/2;
                                        mm_tracks->coloridx[i][spxidy]=3;
                                        rx=rpos->Uniform(-1,1)*mm_tracks->posxerr[i][spxidy];
                                        ry=rpos->Uniform(-1,1)*mm_tracks->posyerr[i][spxidy];
                                        rz=rpos->Uniform(-1,1)*mm_tracks->poszerr[i][spxidy];
                                    }
                                    rgidx=1;
                                    if(spxidy%2==0) {
                                        if(spxidy>=112 && spxidy<116) {
                                            hMM_TrackCounter1[rgidx]++;
                                            hMM_dE1[rgidx]+=rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan];
                                        }
                                        //hMM_TrackCounter2[rgidx]++;
                                        //hMM_dE2[rgidx]+=rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan];
                                    }
                                }else if(asad==3){
                                    for(int i=0;i<64;i++){
                                        mm_tracks->pixel[i+70][spxidy]+=decayIdx*100-1*(decayIdx-1);
                                        mm_tracks->time[i+70][spxidy]=rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan];
                                        mm_tracks->energy[i+70][spxidy]=rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan];
                                        mm_tracks->posx[i+70][spxidy]=i*1.7 + 3*3.4 + 1.7/2;
                                        mm_tracks->posy[i+70][spxidy]=spxidy*1.7 + 1.7/2;
                                        mm_tracks->posxerr[i+70][spxidy]=1.7/2;
                                        mm_tracks->posyerr[i+70][spxidy]=1.7/2;
                                        mm_tracks->poszerr[i+70][spxidy]=fTimePerBin*driftv/2;
                                        mm_tracks->coloridx[i+70][spxidy]=3;
                                        rx=rpos->Uniform(-1,1)*mm_tracks->posxerr[i+70][spxidy];
                                        ry=rpos->Uniform(-1,1)*mm_tracks->posyerr[i+70][spxidy];
                                        rz=rpos->Uniform(-1,1)*mm_tracks->poszerr[i+70][spxidy];
                                    }
                                    rgidx=2;
                                    if(spxidy%2==0) {
                                        if(spxidy>=112 && spxidy<116) {
                                            hMM_TrackCounter1[rgidx]++;
                                            hMM_dE1[rgidx]+=rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan];
                                        }
                                        //hMM_TrackCounter2[rgidx]++;
                                        //hMM_dE2[rgidx]+=rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan];
                                    }
                                }
                            }
                            if(spxidy==-1){ //chain
                                mm_tracks->hasTrackChain+=1;
                                for(int i=0;i<64;i++){
                                    mm_tracks->pixel[spxidx][i*2+1]+=decayIdx*100-1*(decayIdx-1);
                                    mm_tracks->pixel[spxidx][i*2]+=decayIdx*100-1*(decayIdx-1); //Fill strip pixel as well
                                    mm_tracks->time[spxidx][i*2+1]=rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan];
                                    mm_tracks->energy[spxidx][i*2+1]=rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan];
                                    if(spxidx<64) mm_tracks->posx[spxidx][i*2+1]=(spxidx-64)*1.7 - 3*3.4 + 1.7/2;
                                    else mm_tracks->posx[spxidx][i*2+1]=(spxidx-70)*1.7 + 3*3.4 + 1.7/2;
                                    mm_tracks->posy[spxidx][i*2+1]=(i*2+1)*1.7;
                                    mm_tracks->posxerr[spxidx][i*2+1]=1.7/2;
                                    mm_tracks->posyerr[spxidx][i*2+1]=1.7/2;
                                    mm_tracks->poszerr[spxidx][i*2+1]=fTimePerBin*driftv/2;
                                    mm_tracks->coloridx[spxidx][i*2+1]=4;
                                    rx=rpos->Uniform(-1,1)*mm_tracks->posxerr[spxidx][i*2+1];
                                    ry=rpos->Uniform(-1,1)*mm_tracks->posyerr[spxidx][i*2+1];
                                    rz=rpos->Uniform(-1,1)*mm_tracks->poszerr[spxidx][i*2+1];
                                }
                            }
                        }

                        if(cobo==1
                                && rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]>si_minenergy
                                && rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]<si_maxenergy)
                        {
                            //if(asad==0
                            //&& aget==3
                            //&& ( chan==25 || chan==5 || chan==18 || chan==12 || chan==37 || chan==43 || chan==31 || chan==50))
                            if(asad==0 && (aget==0||aget==2))
                            {
                                si_tracks->hasTrack+=1;
                                si_tracks->agetid=aget;
                                si_tracks->chanid=chan;
                                spxidx = mapchantosi->pxidx[asad][aget][chan];
                                spxidy = mapchantosi->pxidy[asad][aget][chan];
                                //if(spxidx>0 && spxidy>0){
                                si_tracks->pixel[spxidx][spxidy]+=decayIdx*100-1*(decayIdx-1);
                                si_tracks->time[spxidx][spxidy]=rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan];
                                si_tracks->energy[spxidx][spxidy]=rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan];
                                //cout << "Si: " << spxidx << " " << spxidy << " " << si_tracks->energy[spxidx][spxidy] << endl;
                                //if(spxidx>=48&&spxidx<60) si_tracks->coloridx[spxidx][spxidy]=5;
                                //if(spxidx>=60&&spxidx<72) si_tracks->coloridx[spxidx][spxidy]=6;
                                //if(spxidx>=72&&spxidx<84) si_tracks->coloridx[spxidx][spxidy]=7;
                                //}
                            }
                            if(asad==1 && aget==1 && chan==2){
                                si_tracks->hasICE = rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan];
                                si_tracks->hasICT = rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan];
                            }
                            if(asad==1 && aget==2 && rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan]>0){
                                for(int i=0;i<9;i++){
                                    if(chan == csi_channo[i]){
                                        si_tracks->hasCsIE[i] = rwaveforms[decayIdx][cobo]->energy[asad*4+aget][chan];
                                        si_tracks->hasCsIT[i] = rwaveforms[decayIdx][cobo]->time[asad*4+aget][chan];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void LKFrameBuilder::FindBoxCorner()
{
    int rgidx=0;
    for(int i=0;i<3;i++){
        hMM_cornerfound[i]=0;
        hMM_MinPxX[i]=0;
        hMM_MinPxY[i]=0;
        hMM_MaxPxX[i]=0;
        hMM_MaxPxY[i]=0;
        hMM_MinX[i]=1000;
        hMM_MinY[i]=1000;
        hMM_MaxX[i]=-1000;
        hMM_MaxY[i]=-1000;
        hMM_Slope[i]=0;
    }

    //if(mm_tracks->hasTrack>0 && si_tracks->hasTrack>0)
    if(mm_tracks->hasTrack>0)
    {
        for(int i=0;i<150;i++){ // pixel_x
            for(int j=0;j<140;j++){ // pixel_y
                if(i>=64 && i<70 && mm_tracks->pixel[i][j]>0){
                    rgidx=0;
                    hMM_cornerfound[rgidx]=1;
                    if(hMM_MinX[rgidx]>mm_tracks->posx[i][j]){
                        hMM_MinX[rgidx] = mm_tracks->posx[i][j];
                        hMM_MinPxX[rgidx] = i;
                    }
                    if(hMM_MinY[rgidx]>mm_tracks->posy[i][j]){
                        hMM_MinY[rgidx] = mm_tracks->posy[i][j];
                        hMM_MinPxY[rgidx] = j;
                    }
                    if(hMM_MaxX[rgidx]<mm_tracks->posx[i][j]){
                        hMM_MaxX[rgidx] = mm_tracks->posx[i][j];
                        hMM_MaxPxX[rgidx] = i;
                    }
                    if(hMM_MaxY[rgidx]<mm_tracks->posy[i][j]){
                        hMM_MaxY[rgidx] = mm_tracks->posy[i][j];
                        hMM_MaxPxY[rgidx] = j;
                    }
                }else if((i>=70) && mm_tracks->pixel[i][j]>1){
                    rgidx=1;
                    hMM_cornerfound[rgidx]=1;
                    if(hMM_MinX[rgidx]>mm_tracks->posx[i][j]){
                        hMM_MinX[rgidx] = mm_tracks->posx[i][j];
                        hMM_MinPxX[rgidx] = i;
                    }
                    if(hMM_MinY[rgidx]>mm_tracks->posy[i][j]){
                        hMM_MinY[rgidx] = mm_tracks->posy[i][j];
                        hMM_MinPxY[rgidx] = j;
                    }
                    if(hMM_MaxX[rgidx]<mm_tracks->posx[i][j]){
                        hMM_MaxX[rgidx] = mm_tracks->posx[i][j];
                        hMM_MaxPxX[rgidx] = i;
                    }
                    if(hMM_MaxY[rgidx]<mm_tracks->posy[i][j]){
                        hMM_MaxY[rgidx] = mm_tracks->posy[i][j];
                        hMM_MaxPxY[rgidx] = j;
                    }
                }else if((i<64) && mm_tracks->pixel[i][j]>1){
                    rgidx=2;
                    hMM_cornerfound[rgidx]=1;
                    if(hMM_MinX[rgidx]>mm_tracks->posx[i][j]){
                        hMM_MinX[rgidx] = mm_tracks->posx[i][j];
                        hMM_MinPxX[rgidx] = i;
                    }
                    if(hMM_MinY[rgidx]>mm_tracks->posy[i][j]){
                        hMM_MinY[rgidx] = mm_tracks->posy[i][j];
                        hMM_MinPxY[rgidx] = j;
                    }
                    if(hMM_MaxX[rgidx]<mm_tracks->posx[i][j]){
                        hMM_MaxX[rgidx] = mm_tracks->posx[i][j];
                        hMM_MaxPxX[rgidx] = i;
                    }
                    if(hMM_MaxY[rgidx]<mm_tracks->posy[i][j]){
                        hMM_MaxY[rgidx] = mm_tracks->posy[i][j];
                        hMM_MaxPxY[rgidx] = j;
                    }
                }
            }
        }
    }
    for(int i=0;i<3;i++){
        if(hMM_cornerfound[i]>0){
            hMM_Slope[i] = (hMM_MaxY[i]-hMM_MinY[i]) / (hMM_MaxX[i]-hMM_MinX[i]);
            hMM_Offset[i] = hMM_MaxY[i]-hMM_Slope[i]*hMM_MaxX[i];
            //cout << (evtcounter%16) << " " << evtcounter << " " << reventIdx << " " << i << " " << hMM_Slope[i] << endl;
        }
    }
}

void LKFrameBuilder::FindX6Hits()
{
    Int_t frameIdx = 0;
    Int_t decayIdx = 0;
    L1Aflag = 0;
    Int_t coboIdx = 0;
    Int_t asadIdx = 0;
    Int_t agetIdx = 0;
    Int_t chanIdx = 0;
    //X6out.open("X6_good_event.txt", std::ofstream::out|std::ofstream::app);
    X6out.open("X6_proton_event.txt", std::ofstream::out|std::ofstream::app);
    Int_t goodCsIhit = 0;

    Int_t num_X6 = 300;
    int Eleft[num_X6][8];
    int Eright[num_X6][8];
    int Eback[num_X6][4];
    int CsIE[num_X6][2];
    for(int i=0;i<num_X6;i++){
        for(int j=0;j<8;j++){
            Eleft[i][j]=0;
            Eright[i][j]=0;
        }
        for(int j=0;j<4;j++){
            Eback[i][j]=0;
        }
        for(int j=0;j<2;j++){
            CsIE[i][j]=0;
        }
    }

    for(Int_t i=0; i<rGETMul; i++){ // find El, Er, Eb, CsIE
        frameIdx = rGETFrameNo[i];
        decayIdx = rGETDecayNo[i];
        coboIdx = rGETCobo[i];
        asadIdx = rGETAsad[i];
        agetIdx = rGETAget[i];
        chanIdx = rGETChan[i];
        if(coboIdx==0 || coboIdx==1) continue;// skip MM waveform data

        Int_t csiloc = 0;
        Int_t csidet = 0;
        Int_t sidet = 0;
        Int_t silr = 0;
        Int_t sistrip = 0;
        //Int_t sichan = 0;
        //if(chanIdx<11) {
        //    sichan = chanIdx;
        //}else if(chanIdx>11 && chanIdx<22) {
        //    sichan = chanIdx - 1;
        //}else if(chanIdx>22 && chanIdx<45) {
        //    sichan = chanIdx - 2;
        //}else if(chanIdx>45 && chanIdx<56) {
        //    sichan = chanIdx - 3;
        //}else if(chanIdx>56) {
        //    sichan = chanIdx - 4;
        //}
        if(coboIdx==2 && !(chanIdx==11 || chanIdx==22 || chanIdx==45 || chanIdx==56)) {
            if(asadIdx==1 && agetIdx==3) { //CsI
                if(rwaveforms[decayIdx][coboIdx]->energy[asadIdx*4+agetIdx][chanIdx]>100){
                    csidet = mapchantox6->CsI_X6det[asadIdx][agetIdx][chanIdx];
                    csiloc = mapchantox6->CsI_X6ud[asadIdx][agetIdx][chanIdx];
                    CsIE[csidet][csiloc]=rwaveforms[decayIdx][coboIdx]->energy[asadIdx*4+agetIdx][chanIdx];
                    goodCsIhit++;
                }
            }
            else if(agetIdx==1 || agetIdx==2) { //Junction
                if(rwaveforms[decayIdx][coboIdx]->energy[asadIdx*4+agetIdx][chanIdx]>100){
                    silr = (mapchantox6->X6strip[asadIdx][agetIdx][chanIdx]+1)%2;
                    sidet = mapchantox6->X6det[asadIdx][agetIdx][chanIdx];
                    sistrip = (int)(mapchantox6->X6strip[asadIdx][agetIdx][chanIdx]+1)/2-1;

                    if(silr==0) Eleft[sidet][sistrip]=rwaveforms[decayIdx][coboIdx]->energy[asadIdx*4+agetIdx][chanIdx];
                    if(silr==1) Eright[sidet][sistrip]=rwaveforms[decayIdx][coboIdx]->energy[asadIdx*4+agetIdx][chanIdx];
                }
            }
            else if(agetIdx==0) { //Ohmic
                if(rwaveforms[decayIdx][coboIdx]->energy[asadIdx*4+agetIdx][chanIdx]>100){
                    sidet = mapchantox6->X6det[asadIdx][agetIdx][chanIdx];
                    sistrip = mapchantox6->X6strip[asadIdx][agetIdx][chanIdx];
                    //cout << sidet << " " << sistrip <<  endl;
                    Eback[sidet][sistrip]=rwaveforms[decayIdx][coboIdx]->energy[asadIdx*4+agetIdx][chanIdx];
                }
            }
        }
    }

    // X6 junction hit pattern
    Double_t Esum=0;
    Double_t Esub=0;
    Double_t Epos=0;
    Double_t posx_X6=0;
    Double_t posy_X6=0;
    Double_t posz_X6=0;
    Double_t X6posz_offset=191.7; //distance from window to MM active area
    Int_t reduceDet = 0;
    for(int i=0;i<num_X6;i++){
        for(int j=0;j<8;j++){
            if(Eleft[i][j]>0 && Eright[i][j]>0){
                Esum = Eleft[i][j]+Eright[i][j];
                Esub = Eleft[i][j]-Eright[i][j];
                Epos = 75/2*(Esub/Esum+1);
                posx_X6 = mapchantox6 -> X6posx[i][j];
                posy_X6 = mapchantox6 -> X6posy[i][j];
                posz_X6 = mapchantox6 -> X6posz[i][j];
                //cout << "Epos" << Epos << endl;
                //hX6_LR[i][j]->Fill(Eright[i][j],Eleft[i][j]);
                //hX6_EsumvsPos[i][j]->Fill(Epos,Esum);
                hX6_LRAll->Fill(Eright[i][j],Eleft[i][j]);
                hX6_EsumvsPosAll->Fill(Epos,Esum);
                if(Esum>1100){
                    if(i>=0 && i<=4) hX6_EsumvsPosAllLeft->Fill(Epos,Esum);
                    else if(i>=8 && i<=12) hX6_EsumvsPosAllRight->Fill(Epos,Esum);

                    if(posy_X6==0) //side
                    {
                        if(posx_X6<0) hX6_LeftHitPattern -> Fill(Epos-75/2,posz_X6);
                        else hX6_RightHitPattern -> Fill(Epos-75/2,posz_X6);
                        posy_X6 = Epos-75/2;
                        //cout << Epos << " " << X6posz[i][j] << endl;
                    }
                    else //bottom
                    {
                        if(posx_X6<0){
                            hX6_BottomHitPattern -> Fill(posx_X6+Epos-75/2,posz_X6);
                            posx_X6 = posx_X6+Epos-75/2;
                        }else{
                            hX6_BottomHitPattern -> Fill(posx_X6-Epos+75/2,posz_X6);
                            posx_X6 = posx_X6-Epos+75/2;
                        }
                    }

                    if((int)i/10==1) si_tracks->hasX6L++;
                    if((int)i/10==2) si_tracks->hasX6R++;
                    if((int)i/100==1) si_tracks->hasX6BL++;
                    if((int)i/100==2) si_tracks->hasX6BR++;
                    if((i>=0 && i<=4)||(i>=8 && i<=12)){
                        rwaveforms[decayIdx][coboIdx]->Sienergy = Esum;
                        rwaveforms[decayIdx][coboIdx]->SiX = posx_X6;
                        rwaveforms[decayIdx][coboIdx]->SiY = posy_X6;
                        rwaveforms[decayIdx][coboIdx]->SiZ = posz_X6-X6posz_offset;
                        //cout << "X6 particles! " << reventIdx << " " << i << ", Front: " << j << " " << Eleft[i][j] << " " << Eright[i][j] << " " << Epos << " " << Esum << " " << si_tracks->hasX6L << " " << si_tracks->hasX6BL << " " << si_tracks->hasX6R << " " << si_tracks->hasX6BR << " " << X6posx[i][j] << " " << X6posy[i][j] << " " << X6posz[i][j] << " " << X6posz[i][j]-X6posz_offset << endl;
                        //X6out << reventIdx << endl;
                    }
                }
            }
        }
    }

    // X6 Ohmic E && CsI PID
    for(int i=0;i<num_X6;i++){
        for(int j=0;j<4;j++){
            if(Eback[i][j]>0){
                //hX6_Eback[i][j]->Fill(Eback[i][j]);
                hX6_EbackAll->Fill(Eback[i][j]);
                goodx6evt++;
                //if(i==15 || i==16 || i==17 || i==25 || i==26 || i==27) goodx6evt++;
                //cout << "X6 particles! " << reventIdx << " " << i << ", Back: " << j << " " << Eback[i][j] << endl;
                if(goodCsIhit>0) {
                    for(int k=0; k<2; k++){
                        //cout << "proton X6 particles! " << reventIdx << " " << i << ", Back: " << j << " " << Eback[i][j] << " " << CsIE[i][k] << endl;
                        if(CsIE[i][k]>0 && cut_pinX6EvsCsIE->IsInside(Eback[i][j],CsIE[i][k])) {
                            hMM_X6EvsCsIEAll -> Fill(CsIE[i][k],Eback[i][j]);
                            goodx6csievt++;
                        }
                    }
                }
            }
        }
    }

    // Ohmic hit pattern
    for(int i=0;i<num_X6;i++){
        for(int j=1;j<5;j++){
            if(Eback[i][j]>0){
                if((int)i/10==1){ //LS
                    if(i%10<5) hX6_OhmicHitPattern -> Fill(-11+(5-j),2+i%10);
                    else hX6_OhmicHitPattern -> Fill(-11+(5-j),-5+i%10);
                }
                else if((int)i/10==2){ //RS
                    hX6_OhmicHitPattern -> Fill(11-(5-j),7-i%10);
                }
                else if((int)i/100==1){ //LB
                    if(i%100==1) hX6_OhmicHitPattern -> Fill(-1,8-(5-j));
                    else if(i%100<5) hX6_OhmicHitPattern -> Fill(-1-(5-j),2+i%100);
                    else hX6_OhmicHitPattern -> Fill(-1-(5-j),-5+i%100);
                }
                else if((int)i/100==2){ //RB
                    if(i%100==4) hX6_OhmicHitPattern -> Fill(1,8-(5-j));
                    else if(i%100<4) hX6_OhmicHitPattern -> Fill(1+(5-j),7-i%100);
                    else hX6_OhmicHitPattern -> Fill(1+(5-j),8-i%100);
                }
            }
        }
    }

    //if(si_tracks->hasFWC>0){
    //  X6out << reventIdx << endl;
    //}
    //X6out.close();
}

//hMM_SiEvsCsIEAll->Fill(csifrontenergy.at(j),sibackenergy.at(i));

void LKFrameBuilder::DrawSiDetector()
{
    int rgidx=0;

    //if(mm_tracks->hasTrack>0 && si_tracks->hasTrack>0)
    if(mm_tracks->hasTrack>3)
    {
        for(int i=0;i<150;i++){ // pixel_x
            if(i>=70){
                rgidx=1;
            }else if(i<64){
                rgidx=2;
            }
            if(hMM_cornerfound[rgidx]>0){
                hMM_Slope[rgidx];
            }
        }
    }
}

void LKFrameBuilder::ReplaceEnergy()
{
    int rgidx=0;
    Double_t tan60deg = TMath::Tan(60*TMath::DegToRad());
    Double_t tan30deg = TMath::Tan(30*TMath::DegToRad());
    Double_t tan45deg = TMath::Tan(45*TMath::DegToRad());
    Double_t gainratio = 900/450;
    Double_t doubleenergy = 0;

    //if(mm_tracks->hasTrack>0 && si_tracks->hasTrack>0)
    if(mm_tracks->hasTrack>3)
    {
        for(int i=0;i<150;i++){ // pixel_x
            if(i>=70){
                rgidx=1;
            }else if(i<64){
                rgidx=2;
            }
            if(hMM_cornerfound[rgidx]>0){
                //if(TMath::Abs(hMM_Slope[rgidx])>tan45deg){
                //  for(int j=0;j<70;j++){ // pixel_y
                //    doubleenergy = mm_tracks->energy[i][j*2];
                //    mm_tracks->energy[i][j*2+1] = doubleenergy/gainratio;
                //  }
                //}else if(TMath::Abs(hMM_Slope[rgidx])<tan45deg){
                //  for(int j=0;j<70;j++){ // pixel_y
                //    doubleenergy = mm_tracks->energy[i][j*2+1];
                //    mm_tracks->energy[i][j*2] = doubleenergy*gainratio;
                //  }
                //}
                if(TMath::Abs(hMM_Slope[rgidx])>tan45deg){
                    for(int j=0;j<70;j++){ // pixel_y
                        mm_tracks->energy[i][j*2+1] = mm_tracks->energy[i][j*2];
                    }
                }else if(TMath::Abs(hMM_Slope[rgidx])<tan45deg){
                    for(int j=0;j<70;j++){ // pixel_y
                        mm_tracks->energy[i][j*2] = mm_tracks->energy[i][j*2+1];
                    }
                }
                //if(TMath::Abs(hMM_Slope[rgidx])>tan60deg){
                //  for(int j=0;j<70;j++){ // pixel_y
                //    mm_tracks->energy[i][j*2+1] = mm_tracks->energy[i][j*2];
                //  }
                //}else if(TMath::Abs(hMM_Slope[rgidx])<tan30deg){
                //  for(int j=0;j<70;j++){ // pixel_y
                //    mm_tracks->energy[i][j*2] = mm_tracks->energy[i][j*2+1];
                //  }
                //}else{
                //  for(int j=0;j<70;j++){ // pixel_y
                //    if(mm_tracks->pixel[i][j*2]>0 && mm_tracks->pixel[i][j*2+1]>0){
                //      mm_tracks->energy[i][j*2] = (mm_tracks->energy[i][j*2]+mm_tracks->energy[i][j*2+1])/2;
                //    }else if(mm_tracks->pixel[i][j*2]>0 && mm_tracks->pixel[i][j*2+1]==0){
                //      mm_tracks->energy[i][j*2] = mm_tracks->energy[i][j*2];
                //    }else if(mm_tracks->pixel[i][j*2]==0 && mm_tracks->pixel[i][j*2+1]>0){
                //      mm_tracks->energy[i][j*2] = mm_tracks->energy[i][j*2+1];
                //    }
                //  }
                //}
            }
        }
    }
}

void LKFrameBuilder::ReplaceEnergybyRatio()
{
    Int_t rgidx=0;
    Int_t maxenergy=0;
    Double_t tan60deg = TMath::Tan(60*TMath::DegToRad());
    Double_t tan30deg = TMath::Tan(30*TMath::DegToRad());
    Double_t tan45deg = TMath::Tan(45*TMath::DegToRad());

    for(int i=0;i<170;i++){
        stripenergysum[i]=0;
        for(int j=0;j<3;j++){
            chainenergysum[i][3]=0;
        }
    }

    //if(mm_tracks->hasTrack>0 && si_tracks->hasTrack>0)
    if(mm_tracks->hasTrack>3)
    {
        for(int i=0;i<150;i++){ // pixel_x
            if(i>=70){
                rgidx=1;
            }else if(i<64){
                rgidx=2;
            }
            for(int j=0;j<140;j++){ // pixel_y
                if(j%2==0){
                    if(i>0) stripenergysum[i] += mm_tracks->energy[i-1][j];
                    stripenergysum[i] += mm_tracks->energy[i][j];
                    stripenergysum[i] += mm_tracks->energy[i+1][j];
                }else{
                    if(j>0) chainenergysum[j][rgidx] += mm_tracks->energy[i][j-1];
                    chainenergysum[j][rgidx] += mm_tracks->energy[i][j];
                    chainenergysum[j][rgidx] += mm_tracks->energy[i][j+1];
                }
            }
        }
    }

    //if(mm_tracks->hasTrack>0 && si_tracks->hasTrack>0)
    if(mm_tracks->hasTrack>0)
    {
        for(int i=0;i<150;i++){ // pixel_x
            if(i>=70){
                rgidx=1;
            }else if(i<64){
                rgidx=2;
            }
            if(hMM_cornerfound[rgidx]>0){
                if(TMath::Abs(hMM_Slope[rgidx])>tan45deg){
                    for(int j=0;j<70;j++){ // pixel_y
                        if(stripenergysum[i]>0){
                            //mm_tracks->energy[i][j*2+1] = mm_tracks->energy[i][j*2+1]*mm_tracks->energy[i][j*2]/stripenergysum[i];
                            maxenergy = GetSumEnergy(i,j*2);
                            mm_tracks->energy[i][j*2+1] = mm_tracks->energy[i][j*2+1]*maxenergy/stripenergysum[i];
                        }else mm_tracks->energy[i][j*2+1] = 0;
                    }
                }else if(TMath::Abs(hMM_Slope[rgidx])<tan45deg){
                    for(int j=0;j<70;j++){ // pixel_y
                        if(chainenergysum[i]>0){
                            //mm_tracks->energy[i][j*2] = mm_tracks->energy[i][j*2]*mm_tracks->energy[i][j*2]/chainenergysum[j*2+1][rgidx];
                            maxenergy = GetSumEnergy(i,j*2+1);
                            mm_tracks->energy[i][j*2] = mm_tracks->energy[i][j*2]*maxenergy/chainenergysum[j*2+1][rgidx];
                        }else mm_tracks->energy[i][j*2] = 0;
                    }
                }
                //if(TMath::Abs(hMM_Slope[rgidx])>tan60deg){
                //  for(int j=0;j<70;j++){ // pixel_y
                //    if(stripenergysum[i]>0){
                //      //mm_tracks->energy[i][j*2+1] = mm_tracks->energy[i][j*2+1]*mm_tracks->energy[i][j*2]/stripenergysum[i];
                //      maxenergy = GetSumEnergy(i,j*2);
                //      mm_tracks->energy[i][j*2+1] = mm_tracks->energy[i][j*2+1]*maxenergy/stripenergysum[i];
                //    }else mm_tracks->energy[i][j*2+1] = 0;
                //  }
                //}else if(TMath::Abs(hMM_Slope[rgidx])<tan30deg){
                //  for(int j=0;j<70;j++){ // pixel_y
                //    if(chainenergysum[i]>0){
                //      //mm_tracks->energy[i][j*2] = mm_tracks->energy[i][j*2]*mm_tracks->energy[i][j*2]/chainenergysum[j*2+1][rgidx];
                //      maxenergy = GetSumEnergy(i,j*2+1);
                //      mm_tracks->energy[i][j*2] = mm_tracks->energy[i][j*2]*maxenergy/chainenergysum[j*2+1][rgidx];
                //    }else mm_tracks->energy[i][j*2] = 0;
                //  }
                //}else{
                //  for(int j=0;j<70;j++){ // pixel_y
                //    if(mm_tracks->pixel[i][j*2]>0 && mm_tracks->pixel[i][j*2+1]>0){
                //      mm_tracks->energy[i][j*2] = (mm_tracks->energy[i][j*2]+mm_tracks->energy[i][j*2+1])/2;
                //    }else if(mm_tracks->pixel[i][j*2]>0 && mm_tracks->pixel[i][j*2+1]==0){
                //      mm_tracks->energy[i][j*2] = mm_tracks->energy[i][j*2];
                //    }else if(mm_tracks->pixel[i][j*2]==0 && mm_tracks->pixel[i][j*2+1]>0){
                //      mm_tracks->energy[i][j*2] = mm_tracks->energy[i][j*2+1];
                //    }
                //  }
                //}
            }
        }
    }
}

Int_t LKFrameBuilder::GetSumEnergy(Int_t px, Int_t py){
    Int_t maxenergy=mm_tracks->energy[px][py];
    if(px>0) maxenergy+=mm_tracks->energy[px-1][py];
    maxenergy+=mm_tracks->energy[px+1][py];

    return maxenergy;
}

void LKFrameBuilder::CleanTrack()
{
    //if(mm_tracks->hasTrack>0 && si_tracks->hasTrack>0)
    if(mm_tracks->hasTrack>0)
    {
        for(int i=0;i<150;i++){ // pixel_x
            for(int j=0;j<140;j++){ // pixel_y
                if((i<64 || i>=70) && mm_tracks->pixel[i][j]>0){
                    if(mm_tracks->pixel[i][j+1]>0){
                        //if(mm_tracks->time[i][j]==mm_tracks->time[i][j+1]) //for fast drift velocity
                        if(mm_tracks->time[i][j]<=((mm_tracks->time[i][j+1])+5) && mm_tracks->time[i][j]>=((mm_tracks->time[i][j+1])-5))//for slow drift velocity
                        {
                            //mm_tracks->pixel[i][j] = mm_tracks->pixel[i][j]; //for fast drift velocity
                        }else{
                            mm_tracks->pixel[i][j] = 0;
                        }
                    }else{
                        mm_tracks->pixel[i][j] = 0;
                    }
                }
            }
        }
    }
}

void LKFrameBuilder::FillDecayFlag()
{
    UInt_t dchan = 0;
    UInt_t spxidx = 0;
    UInt_t spxidy = 0;
    Int_t decayIdxMax=1;
    if(enable2pmode==1) decayIdxMax=2;

    for(int i=0;i<150;i++) {for(int j=0;j<170;j++) { mm_tracks->L1B[i][j]=false;}}

    for(Int_t decayIdx=0; decayIdx<decayIdxMax; decayIdx++) {
        for(UInt_t cobo=0; cobo<no_cobos; cobo++) {
            for(UInt_t asad=0; asad<maxasad; asad++) {
                for(UInt_t aget=0; aget<4; aget++) {
                    for(UInt_t chan=0; chan<68; chan++) {
                        //cout << "FDF:" << cobo << " " << asad << " " << aget << " " << chan << " " << rwaveforms[decayIdx][cobo]->isDecay[asad*4+aget][chan] << endl;
                        if(decayIdx==1) mm_tracks->hasDecay=true;
                        if(chan<11) {
                            dchan = chan;
                        }else if(chan>11 && chan<22) {
                            dchan = chan - 1;
                        }else if(chan>22 && chan<45) {
                            dchan = chan - 2;
                        }else if(chan>45 && chan<56) {
                            dchan = chan - 3;
                        }else if(chan>56) {
                            dchan = chan - 4;
                        }

                        spxidx = mapchantomm->pxidx[asad][aget][dchan];
                        spxidy = mapchantomm->pxidy[asad][aget][dchan];
                        //  mm_tracks->L1B[spxidx][spxidy]=rwaveforms[decayIdx][cobo]->isDecay[asad*4+aget][chan];
                        if(spxidx!=-1 && spxidy!=-1 && mm_tracks->pixel[spxidx][spxidy]>(0+decayIdx*99)){
                            mm_tracks->decay[spxidx][spxidy]=decayIdx;
                        }
                        if(spxidx==-1){ //strip
                            if(asad==2){
                                for(int i=0;i<64;i++){
                                    if(mm_tracks->pixel[i][spxidy]>(0+decayIdx*99)){
                                        mm_tracks->decay[i][spxidy]=decayIdx;
                                    }
                                }
                            }else if(asad==3){
                                for(int i=0;i<64;i++){
                                    if(mm_tracks->pixel[i+70][spxidy]>(0+decayIdx*99)){
                                        mm_tracks->decay[i+70][spxidy]=decayIdx;
                                    }
                                }
                            }
                        }
                        if(spxidy==-1){ //chain
                            for(int i=0;i<64;i++){
                                if(mm_tracks->pixel[spxidx][i*2+1]>(0+decayIdx*99)){
                                    mm_tracks->decay[spxidx][i*2+1]=decayIdx;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void LKFrameBuilder::GetTrackPosYLimit()
{
    int rgidx=0;
    //if(mm_tracks->hasTrack>0 && si_tracks->hasTrack>0)
    if(mm_tracks->hasTrack>3)
    {
        for(rgidx=0;rgidx<3;rgidx++){
            trackposymin[rgidx] = 10000;
            trackposymax[rgidx] = -10000;
        }
        for(int i=0;i<150;i++){
            for(int j=0;j<170;j++){
                if(mm_tracks->pixel[i][j]>0){
                    if(i>63 && i<71) rgidx=0;
                    if(i<64) rgidx=1;
                    if(i>70) rgidx=2;
                    if(trackposymin[rgidx]>mm_tracks->posy[i][j]) trackposymin[rgidx] = mm_tracks->posy[i][j];
                    if(trackposymax[rgidx]<mm_tracks->posy[i][j]) trackposymax[rgidx] = mm_tracks->posy[i][j];
                }
            }
        }
    }
}

void LKFrameBuilder::HoughTransform()
{
    double theta;
    double costheta;
    double sintheta;
    double radiusxy;
    double radiusxt;
    double radiusyt;
    int rgidx=0;
    int avgflag=0;
    //if(mm_tracks->hasTrack>3 && si_tracks->hasTrack>0)
    if(mm_tracks->hasTrack>3)
    {
        //if(mm_tracks->hasTrack>3 && mm_tracks->hasTrackChain>10 && mm_tracks->hasTrackStrip>10)
        for(int i=0;i<150;i++)
        {
            for(int j=0;j<170;j++){
                if(mm_tracks->pixel[i][j]>0){
                    if(i>63 && i<71) rgidx=0;
                    if(i<64) rgidx=1;
                    if(i>70) rgidx=2;
                    for(theta=0;theta<180;theta+=0.1){
                        costheta = TMath::Cos(theta*TMath::DegToRad());
                        sintheta = TMath::Sin(theta*TMath::DegToRad());
                        if(rgidx==0){
                            radiusxy = (mm_tracks->avgposx[j])*costheta+(mm_tracks->avgposy[j])*sintheta;
                        }else if(rgidx>0){
                            radiusxy = (mm_tracks->posx[i][j])*costheta+(mm_tracks->posy[i][j])*sintheta;
                        }
                        if(radiusxy>-300 && radiusxy<300) hMM_TrackPosHough[goodevtcounter%16][rgidx]->Fill(theta,radiusxy);
                        radiusxt = (mm_tracks->posx[i][j])*costheta+(mm_tracks->time[i][j])*sintheta;
                        if(radiusxt>-400 && radiusxt<600) hMM_TimevsPxIDXPosHough[goodevtcounter%16][rgidx]->Fill(theta,radiusxt);
                        radiusyt = (mm_tracks->posy[i][j])*costheta+(mm_tracks->time[i][j])*sintheta;
                        if(radiusyt>-400 && radiusyt<600) hMM_TimevsPxIDYPosHough[goodevtcounter%16][rgidx]->Fill(theta,radiusyt);
                    }
                }
            }
        }
    }
}

void LKFrameBuilder::GetXYZTrack()
{
    int maxbinxy=0;
    int maxbinxt=0;
    int maxbinyt=0;
    int xbin=0;
    int ybin=0;
    int zbin=0;
    double thetaxy=0;
    double thetaxt=0;
    double thetayt=0;
    double costhetaxy=0;
    double sinthetaxy=0;
    double costhetaxt=0;
    double sinthetaxt=0;
    double costhetayt=0;
    double sinthetayt=0;
    double radiusxy=0;
    double radiusxt=0;
    double radiusyt=0;
    double stddevxy=0;
    double minstddevxy=100000;
    double minthetayt=0;
    double stddevyt=0;
    double minstddevyt=100000;
    double minthetaxt=0;
    double stddevxt=0;
    double minstddevxt=100000;
    int entries=0;
    double maxr=0;
    double maxtheta=0;
    int binmax=0;
    int maximum=0;
    //int minentries[3]={20,10,10};
    int minentries[3]={0,0,0};
    //if(mm_tracks->hasTrack>3 && si_tracks->hasTrack>0)
    if(mm_tracks->hasTrack>3)
    {
        //if(mm_tracks->hasTrack>0 && mm_tracks->hasTrackChain>10 && mm_tracks->hasTrackStrip>10)
        for(int rgidx=0;rgidx<3;rgidx++)
        {
            minstddevxy=100000;
            minstddevyt=100000;
            minstddevxt=100000;
            entries = hMM_TrackPosHough[goodevtcounter%16][rgidx]->ProjectionY("projy",450,451)->GetEntries();
            if(entries>minentries[rgidx]){
                for(int i=0;i<1800;i++){
                    if(hMM_TrackPosHough[goodevtcounter%16][rgidx]->ProjectionY("projy",i,i+1)->GetEntries()==entries){
                        stddevxy = hMM_TrackPosHough[goodevtcounter%16][rgidx]->ProjectionY("projy",i,i+1)->GetStdDev();
                        if(minstddevxy>=stddevxy){
                            minstddevxy = stddevxy;
                            thetaxy = hMM_TrackPosHough[goodevtcounter%16][rgidx]->GetXaxis()->GetBinCenter(i);
                            ybin = hMM_TrackPosHough[goodevtcounter%16][rgidx]->ProjectionY("projy",i,i+1)->GetMaximumBin();
                            radiusxy = hMM_TrackPosHough[goodevtcounter%16][rgidx]->GetYaxis()->GetBinCenter(ybin);
                        }

                        stddevxt = hMM_TimevsPxIDXPosHough[goodevtcounter%16][rgidx]->ProjectionY("projy",i,i+1)->GetStdDev();
                        if(minstddevxt>stddevxt){
                            minstddevxt = stddevxt;
                            thetaxt = hMM_TimevsPxIDXPosHough[goodevtcounter%16][rgidx]->GetXaxis()->GetBinCenter(i);
                            xbin = hMM_TimevsPxIDXPosHough[goodevtcounter%16][rgidx]->ProjectionY("projy",i,i+1)->GetMaximumBin();
                            radiusxt = hMM_TimevsPxIDXPosHough[goodevtcounter%16][rgidx]->GetYaxis()->GetBinCenter(xbin);
                        }
                        stddevyt = hMM_TimevsPxIDYPosHough[goodevtcounter%16][rgidx]->ProjectionY("projy",i,i+1)->GetStdDev();
                        if(minstddevyt>stddevyt){
                            minstddevyt = stddevyt;
                            thetayt = hMM_TimevsPxIDYPosHough[goodevtcounter%16][rgidx]->GetXaxis()->GetBinCenter(i);
                            ybin = hMM_TimevsPxIDYPosHough[goodevtcounter%16][rgidx]->ProjectionY("projy",i,i+1)->GetMaximumBin();
                            radiusyt = hMM_TimevsPxIDYPosHough[goodevtcounter%16][rgidx]->GetYaxis()->GetBinCenter(ybin);
                        }
                    }
                }

                costhetaxy = TMath::Cos(thetaxy*TMath::DegToRad());
                sinthetaxy = TMath::Sin(thetaxy*TMath::DegToRad());
                double ifzeroxy = sinthetaxy*costhetaxy;
                if(ifzeroxy==0){
                    thetaxy+=0.1;
                    costhetaxy = TMath::Cos(thetaxy*TMath::DegToRad());
                    sinthetaxy = TMath::Sin(thetaxy*TMath::DegToRad());
                }
                maxbinxt = hMM_TimevsPxIDXPosHough[goodevtcounter%16][rgidx]->GetMaximumBin();
                costhetaxt = TMath::Cos(thetaxt*TMath::DegToRad());
                sinthetaxt = TMath::Sin(thetaxt*TMath::DegToRad());
                double ifzeroxt = sinthetaxt*costhetaxt;
                if(ifzeroxt==0){
                    thetaxt+=0.1;
                    costhetaxt = TMath::Cos(thetaxt*TMath::DegToRad());
                    sinthetaxt = TMath::Sin(thetaxt*TMath::DegToRad());
                }
                maxbinyt = hMM_TimevsPxIDYPosHough[goodevtcounter%16][rgidx]->GetMaximumBin();
                costhetayt = TMath::Cos(thetayt*TMath::DegToRad());
                sinthetayt = TMath::Sin(thetayt*TMath::DegToRad());
                double ifzeroyt = sinthetayt*costhetayt;
                if(ifzeroyt==0){
                    thetayt+=0.1;
                    costhetayt = TMath::Cos(thetayt*TMath::DegToRad());
                    sinthetayt = TMath::Sin(thetayt*TMath::DegToRad());
                }
                int i = 0; // pixel_x
                double posxt = 0; //pos_xt
                double posx = 0; //pos_x
                double posy = 0; //pos_y
                double posz = 0; //pos_z
                double timemm = 0; // time_mm
                double timez = 0; // time_z
                double timesi = 0; // time_si
                int mmcounts=0;
                int maxmmcounts=0;
                for(int i=0;i<150;i++){ // pixel_x
                    for(int j=140;j<170;j++){ // pixel_y
                        if(si_tracks->time[i][j]>0) timesi = si_tracks->time[i][j]*fTimePerBin;
                    }
                }
                posx = (radiusxy+300*sinthetaxy)/costhetaxy;
                if(rgidx==0 && TMath::Abs(posx)>5){
                    for(int i=65;i<69;i++){ // pixel_x
                        mmcounts=0;
                        for(int j=0;j<128;j++){ // pixel_y
                            if(mm_tracks->pixel[i][j]>0) mmcounts++;
                        }
                        if(mmcounts>maxmmcounts){
                            radiusxy = (i-67)*3.4+3.4/2;
                            maxmmcounts = mmcounts;
                        }
                    }
                    thetaxy=0;
                    costhetaxy = 1;
                    sinthetaxy = 0;
                }
                for(int j=-300;j<400;j++){ // pos_y (mm)
                    //posxt = (radiusxt - sinthetaxt*(radiusyt-j*costhetayt)/sinthetayt)/costhetaxt;
                    posx = (radiusxy-j*sinthetaxy)/costhetaxy;
                    timemm = (radiusyt-j*costhetayt)/sinthetayt*fTimePerBin;
                    //timez = timemm - timesi;
                    //posz = 70 - (timez * driftv);
                    posz = 0;
                    timez = timemm;
                    //if(j==0) cout << a0 << " " << a1 << " " << posx << " " << b0 << " " << b1 << " " << timemm << " " << timesi << " " << timez << " " << posz << endl;
                    if(j>=-164 && j<-100) hMM_TrackPosXZ[(164+j)]->Fill(posx,timez);
                    if(j>=-205 && j<-195) hMM_TrackPosXZAll->Fill(posx,timez);
                    for(int k=0;k<10;k++){
                        //hMM_TrackPos[goodevtcounter%16]->Fill(posxt,j);
                        //hMM_TrackPosXY[goodevtcounter%16]->Fill(posxt,j);
                    }
                    if(j<=trackposymax[rgidx]){
                        hMM_TrackPos[goodevtcounter%16]->Fill(posx,j);
                        hMM_TrackPosXY[goodevtcounter%16]->Fill(posx,j);
                        hMM_TrackPosYZ[goodevtcounter%16]->Fill(j,posz);
                        /*
                           if(j<0 && ((rgidx==0) || (posx<0 && rgidx==1) || (posx>0 && rgidx==2))){
                           hMM_TrackPos[goodevtcounter%16]->Fill(posx,j);
                           hMM_TrackPosXY[goodevtcounter%16]->Fill(posx,j);
                           hMM_TrackPosYZ[goodevtcounter%16]->Fill(j,posz);

                        //if(posxt<150){
                        //  hMM_TrackPosXYTAll->Fill(posxt,j);
                        //  timemm = (radiusxt-posxt*costhetaxt)/sinthetaxt;
                        //  hMM_TimevsPxIDXPos[goodevtcounter%16]->Fill(posxt,timemm);
                        //}
                        }else if(j>=0 && ((posx>=-10.2 && posx<=10.2 && rgidx==0) || (posx<-10.2 && rgidx==1) || (posx>10.2 && rgidx==2))){
                        hMM_TrackPos[goodevtcounter%16]->Fill(posx,j);
                        hMM_TrackPosXY[goodevtcounter%16]->Fill(posx,j);
                        hMM_TrackPosYZ[goodevtcounter%16]->Fill(j,posz);
                        //if(posxt<150){
                        //  hMM_TrackPosXYTAll->Fill(posxt,j);
                        //  timemm = (radiusxt-posxt*costhetaxt)/sinthetaxt;
                        //  hMM_TimevsPxIDXPos[goodevtcounter%16]->Fill(posxt,timemm);
                        //}
                        }
                         */
                    }
                    hMM_TrackPosXYAll->Fill(posx,j);
                    hMM_TrackPosYZAll->Fill(j,posz);
                    timemm = (radiusyt-j*costhetayt)/sinthetayt;
                    hMM_TimevsPxIDYPos[goodevtcounter%16]->Fill(j,timemm);
                    //hMM_TrackXZ[goodevtcounter%16]->Fill(i,timez);
                    //hMM_TrackYZ[goodevtcounter%16]->Fill(j,timez);
                    //hMM_TrackXZ[goodevtcounter%16]->Fill(i,posz);
                    //hMM_TrackYZ[goodevtcounter%16]->Fill(j,posz);
                }
                hMM_a0[rgidx]->Fill(costhetaxy/sinthetaxy*-1);
                hMM_a1[rgidx]->Fill(radiusxy/sinthetaxy);
                hMM_b0[rgidx]->Fill((costhetayt/sinthetayt*-1)*fTimePerBin);
                hMM_b1[rgidx]->Fill(((radiusyt-260*sinthetayt)/costhetayt)*fTimePerBin);
                if(rgidx>0){
                    hMM_a0[3]->Fill(costhetaxy/sinthetaxy*-1);
                    hMM_a1[3]->Fill(radiusxy/sinthetaxy);
                    hMM_b0[3]->Fill((costhetayt/sinthetayt*-1)*fTimePerBin);
                    hMM_b1[3]->Fill(((radiusyt-260*sinthetayt)/costhetayt)*fTimePerBin);
                }
            }
        }
    }
}

void LKFrameBuilder::FilldEvsE()
{
    //if(mm_tracks->hasTrack>0 && si_tracks->hasTrack>0 && si_tracks->hasICE>0)
    //if(mm_tracks->hasTrack>0 && si_tracks->hasTrack>0)
    int detno=0;
    int rgidx=0;
    if(mm_tracks->hasTrack>3){
        for(int i=0;i<10;i++){
            for(int j=0;j<4;j++){
                if(si_tracks->pixel[i][j]>0){
                    //if(si_tracks->energy[i][j]>0 &&(si_tracks->energy[i][j]>hMM_E[rgidx])) hMM_E[rgidx] = si_tracks->energy[i][j];
                    if(si_tracks->energy[i][j]>0){
                        //cout << "Si: " << i << " " << j << " " << si_tracks->energy[i][j] << endl;
                        if(i<2 && j<2) detno=9;
                        else if(i<2 && j>=2) detno=8;
                        else if((i>=2&&i<4) && j<2) detno=7;
                        else if((i>=2&&i<4) && j>=2) detno=6;
                        else if((i>=4&&i<6) && j<2) detno=5;
                        else if((i>=4&&i<6) && j>=2) detno=4;
                        else if((i>=6&&i<8) && j<2) detno=3;
                        else if((i>=6&&i<8) && j>=2) detno=2;
                        else if((i>=8&&i<10) && j<2) detno=1;
                        else if((i>=8&&i<10) && j>=2) detno=0;

                        if(detno==0||detno==1||detno==8||detno==9){
                            hMM_E[0] += si_tracks->energy[i][j];
                        }else if(detno==2||detno==3||detno==4||detno==6||detno==7){
                            hMM_E[1] += si_tracks->energy[i][j];
                        }else{
                            hMM_E[2] += si_tracks->energy[i][j];
                        }
                    }
                }
            }
        }
        for(int j=0;j<170;j++){
            if(sumcounter[j]>0){
                mm_tracks->avgposx[j] = posenergysum[j]/energysum[j];
                mm_tracks->sumenergy[j] = energysum[j];
            }
        }
        for(int i=64;i<71;i++){
            for(int j=0;j<170;j++){
                if(mm_tracks->pixel[i][j]>0){
                    mm_tracks->avgposy[j] = mm_tracks->posy[i][j];
                }
            }
        }
    }
}

void LKFrameBuilder::Sum2pEnergy()
{
    int rgidx=0;
    if(mm_tracks->hasTrack>3 && (mm_tracks->hasDecay)){
        for(int i=0;i<150;i++){
            for(int j=0;j<170;j++){
                //Micromega
                if(mm_tracks->decay[i][j]>0){
                    if(i>63 && i<71) rgidx=0;
                    if(i<64) rgidx=1;
                    if(i>70) rgidx=2;
                    mm_tracks->sum2penergy[rgidx] += mm_tracks->energy[i][j];
                    mm_tracks->sum2penergy[3] += mm_tracks->energy[i][j];
                }
            }
        }
        //cout << "sum = " << mm_tracks->sum2penergy << endl;;
        for(int i=0;i<4;i++){
            if(mm_tracks->sum2penergy[i]>0){
                hMM_Sum2pEnergy[i]->Fill(mm_tracks->sum2penergy[i]);
            }
        }
        if(mm_tracks->sum2penergy[3]>0){
            //hMM_TrackvsE[goodevtcounter%16]->SetTitle(Form("Single Event Track vs E;Cell ID X [SumE(bc,bl,br,all)=%d,%d,%d,%d];Cell ID Y [EventNo=%d, goodevtcounter=%d]",mm_tracks->sum2penergy[0],mm_tracks->sum2penergy[1],mm_tracks->sum2penergy[2],mm_tracks->sum2penergy[3],reventIdx,goodevtcounter));
            //hMM_TrackvsE[gatedevtcounter%16]->SetTitle(Form("Single Event Track vs E;Cell ID X;Cell ID Y [EventNo=%d, goodevtcounter=%d]",reventIdx,goodevtcounter));
        }
        hMM_TrackvsE[goodevtcounter%16]->SetTitle(Form("Single Event Track vs E;Cell ID X;Cell ID Y [EventNo=%d, goodevtcounter=%d]",reventIdx,goodevtcounter));
    }
}

void LKFrameBuilder::DrawSumEnergyTrack()
{
    int rgidx=0;
    int sumstrip[170];
    int laststrip=0;
    int stripcounter[170];
    int fusionflag=0;
    for(int i=0;i<170;i++){
        sumstrip[i] = 0;
        stripcounter[i] = 0;
    }
    if(mm_tracks->hasTrack>0){
        for(int i=0;i<150;i++){
            for(int j=0;j<170;j++){
                //Micromega
                if(mm_tracks->pixel[i][j]>0){
                    if(i<64 || i>69){
                    }else{
                        if(mm_tracks->energy[i][j]>0){
                            //sumstrip[j]+=mm_tracks->energy[i][j];
                            //stripcounter[j]++;
                        }
                    }
                    if(i>=64 && i<70){
                        if(mm_tracks->energy[i][j]>0){
                            sumstrip[j]+=mm_tracks->energy[i][j];
                            stripcounter[j]++;
                        }
                    }
                }
                if(j==0 && i==70) {
                    hMMICenergyvsMME->Fill(mm_tracks->ICenergy,sumstrip[0]);
                }
            }
        }

        for(int j=0;j<170;j++){
            if(sumstrip[j]>0 && stripcounter[j]>0){
                hMM_SumEnergyvsPxIDY[goodevtcounter%16]->Fill(j,(1.0*((double)sumstrip[j])/stripcounter[j]));
                //cout << j << " " << sumstrip[j] << " " << stripcounter[j] << endl;
                if(sumstrip[j]>100){
                    laststrip = j;
                }
            }
        }
        hMM_SumEnergyMaxIDY->Fill(hMM_SumEnergyvsPxIDY[goodevtcounter%16]->GetMaximumBin());
        hMM_SumEnergyLastIDY->Fill(laststrip);
        //hMM_SumEnergyvsPxIDY[goodevtcounter%16]->Fill(130,mm_tracks->ICenergy);
        hMM_SumEnergyvsPxIDY[goodevtcounter%16]->Fill(131,mm_tracks->Sienergy);

        for(int j=0;j<170;j++){
            if(sumstrip[j]>0 && stripcounter[j]>0){
                if(1.0*((double)sumstrip[j])/stripcounter[j]>1300 && (j>1 && j<127)) fusionflag=1;
                //cout << j << " " << sumstrip[j] << " " << stripcounter[j] << endl;
            }
        }
        if(fusionflag==1){
            for(int i=0;i<150;i++){
                for(int j=0;j<170;j++){
                    if(mm_tracks->pixel[i][j]>0){
                        //hMM_TrackAll->Fill(i,j);
                        //hMM_TrackvsE[gatedevtcounter%16]->Fill(i,j,mm_tracks->energy[i][j]);
                    }
                }
            }
            for(int j=0;j<170;j++){
                if(sumstrip[j]>0 && stripcounter[j]>0){
                    //hMM_SumEnergyvsPxIDY[gatedevtcounter%16]->Fill(j,(1.0*((double)sumstrip[j])/stripcounter[j]));
                }
            }
            gatedevtcounter++;
        }
    }

    int detno=0;
    for(int i=0;i<10;i++){
        for(int j=0;j<4;j++){
            if(si_tracks->pixel[i][j]>0){
                if(i<2 && j<2) detno=9;
                else if(i<2 && j>=2) detno=8;
                else if((i>=2&&i<4) && j<2) detno=7;
                else if((i>=2&&i<4) && j>=2) detno=6;
                else if((i>=4&&i<6) && j<2) detno=5;
                else if((i>=4&&i<6) && j>=2) detno=4;
                else if((i>=6&&i<8) && j<2) detno=3;
                else if((i>=6&&i<8) && j>=2) detno=2;
                else if((i>=8&&i<10) && j<2) detno=1;
                else if((i>=8&&i<10) && j>=2) detno=0;

                if(si_tracks->hasCsIE[detno]>100){
                    hMM_SumEnergyvsPxIDY[goodevtcounter%16]->Fill(132,si_tracks->energy[i][j]);
                    hMM_SumEnergyvsPxIDY[goodevtcounter%16]->Fill(134,si_tracks->hasCsIE[detno]);
                    hGET_SiForwardHitPattern2D->Fill(i,j);
                }
            }
        }
    }
}

void LKFrameBuilder::ChangeTrackHistTitle()
{
    TString title0 = Form("[D2PTime=%d usec, FrameNo~%d, EventNo=%d, goodevtcounter=%d]",int(rd2ptime/1000),(reventIdx-FirsteventIdx+2),reventIdx,goodevtcounter);
    hMM_Track[goodevtcounter%16]->SetTitle(Form("Single Event Track by channels %s;Cell ID X;Cell ID Y",title0.Data()));
    hMM_TrackvsE[goodevtcounter%16]->SetTitle(Form("Single Event of (X6:%d %d %d %d) (%d) Track vs E;Cell ID X;Cell ID Y [goodevtcounter=%d]",si_tracks->hasX6L,si_tracks->hasX6R,si_tracks->hasX6BL,si_tracks->hasX6BR,reventIdx,goodevtcounter));
    hMM_TrackDecay[goodevtcounter%16]->SetTitle(Form("Single Event Decay Track by channels %s;Cell ID X;Cell ID Y",title0.Data()));
    hMM_TrackDecay2[goodevtcounter%16]->SetTitle(Form("Single Event Decay Track by channels %s;Cell ID X;Cell ID Y",title0.Data()));
    hMM_TrackPos[goodevtcounter%16]->SetTitle(Form("Single Event Track by Position %s;Pos X (mm,title0.Data()));Pos Y (mm)",title0.Data()));
    hMM_TrackPosXY[goodevtcounter%16]->SetTitle(Form("Single Event Track Pos Y vs Pos X %s;Pos X (mm,title0.Data()));Pos Y (mm)",title0.Data()));
    hMM_TimevsPxIDX[goodevtcounter%16]->SetTitle(Form("Single Event Time vs Pixel X %s;Cell ID X;Time Bucket",title0.Data()));
    hMM_TimevsPxIDXPos[goodevtcounter%16]->SetTitle(Form("Single Event Time vs Pos X %s;Pos X (mm,title0.Data()));Time Bucket",title0.Data()));
    hMM_TimevsPxIDY[goodevtcounter%16]->SetTitle(Form("Single Event Time vs Pixel Y %s;Cell ID Y;Time Bucket",title0.Data()));
    hMM_TimevsPxIDYPos[goodevtcounter%16]->SetTitle(Form("Single Event Time vs Pos Y %s;Pos Y (mm,title0.Data()));Time Bucket",title0.Data()));
    hMM_EnergyvsPxIDY[goodevtcounter%16]->SetTitle(Form("Single Event Energy vs Pixel Y %s;Cell ID Y;Energy",title0.Data()));
    hMM_SumEnergyvsPxIDY[goodevtcounter%16]->SetTitle(Form("Single Event Sum Energy vs Pixel Y %s;Cell ID Y",reventIdx,title0.Data()));
}

void LKFrameBuilder::DrawTrack()
{
    int rgidx=0;
    int fusionflag=0;
    int EstripL=0,EstripR=0;
    //if(mm_tracks->hasTrack>0 && si_tracks->hasTrack>0 && si_tracks->hasICE>0)
    //if(mm_tracks->hasTrack>0 && si_tracks->hasTrack>0)
    if(mm_tracks->hasTrack>0)
    {
        //cout<<"BOOP"<<endl;
        if(mm_tracks->SiX<0) hMM_TrackvsE[goodevtcounter%16]->Fill(1,mm_tracks->SiZ/1.75,mm_tracks->Sienergy);
        if(mm_tracks->SiX>0) hMM_TrackvsE[goodevtcounter%16]->Fill(140,mm_tracks->SiZ/1.75,mm_tracks->Sienergy);
        //for(int i=0;i<150;i++)
        //for(int j=0;j<170;j++)
        for(int i=10;i<116;i++) // ignore the left 10 and right 10 pixels
        {
            for(int j=10;j<115;j++){ // ignore the first 10 pixels
                //Micromega
                if(mm_tracks->pixel[i][j]>0){
                    if(i>63 && i<71) rgidx=0;
                    if(i<64) rgidx=1;
                    if(i>70) rgidx=2;
                    hMM_Time[goodevtcounter%16]->Fill(mm_tracks->time[i][j]);
                    hMM_Energy[goodevtcounter%16]->Fill(mm_tracks->energy[i][j]);
                    hMM_TrackAll->Fill(i,j);
                    hMM_TrackPosAll->Fill(mm_tracks->posx[i][j]+rx,mm_tracks->posy[i][j]+ry);
                    hMM_TimevsPxIDXPosAll->Fill(mm_tracks->posx[i][j]+rx,mm_tracks->time[i][j]);
                    hMM_TimevsPxIDYPosAll->Fill(mm_tracks->posy[i][j]+ry,mm_tracks->time[i][j]);
                    hMM_TrackvsE[goodevtcounter%16]->Fill(i,j,mm_tracks->energy[i][j]);
                    hMM_Track[goodevtcounter%16]->Fill(i,j,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
                    hMM_TrackPos[goodevtcounter%16]->Fill(mm_tracks->posx[i][j]+rx,mm_tracks->posy[i][j]+ry,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
                    if(i<64 || i>69){
                        if(j%2!=0){ //chain
                            hMM_TimevsPxIDX[goodevtcounter%16]->Fill(i,mm_tracks->time[i][j],mm_tracks->coloridx[i][j]);
                            hMM_TimevsPxIDXPos[goodevtcounter%16]->Fill(mm_tracks->posx[i][j]+rx,mm_tracks->time[i][j],mm_tracks->coloridx[i][j]);
                            hMM_EnergyvsPxIDY[goodevtcounter%16]->Fill(j,mm_tracks->energy[i][j],mm_tracks->coloridx[i][j]);
                        }else{ //strip
                            hMM_TimevsPxIDYAll->Fill(j,mm_tracks->time[i][j]);
                            hMM_TimevsPxIDY[goodevtcounter%16]->Fill(j,mm_tracks->time[i][j],mm_tracks->coloridx[i][j]);
                            hMM_TimevsPxIDYPos[goodevtcounter%16]->Fill(mm_tracks->posy[i][j]+ry,mm_tracks->time[i][j],mm_tracks->coloridx[i][j]);
                            hMM_EnergyvsPxIDY[goodevtcounter%16]->Fill(j,mm_tracks->energy[i][j],mm_tracks->coloridx[i][j]);
                        }
                    }else{
                        hMM_TimevsPxIDX[goodevtcounter%16]->Fill(i,mm_tracks->time[i][j],mm_tracks->coloridx[i][j]);
                        hMM_TimevsPxIDXPos[goodevtcounter%16]->Fill(mm_tracks->posx[i][j]+rx,mm_tracks->time[i][j],mm_tracks->coloridx[i][j]);
                        //hMM_TimevsPxIDYAll->Fill(j,mm_tracks->time[i][j]);
                        hMM_TimevsPxIDY[goodevtcounter%16]->Fill(j,mm_tracks->time[i][j],mm_tracks->coloridx[i][j]);
                        hMM_TimevsPxIDYPos[goodevtcounter%16]->Fill(mm_tracks->posy[i][j]+ry,mm_tracks->time[i][j],mm_tracks->coloridx[i][j]);
                        hMM_EnergyvsPxIDY[goodevtcounter%16]->Fill(j,mm_tracks->energy[i][j],mm_tracks->coloridx[i][j]);
                        hMM_SumEnergyvsPxIDY[goodevtcounter%16]->Fill(j,mm_tracks->energy[i][j]);
                    }
                    if(i>=64 && i<70){
                        hMM_EnergyvsPxIDYALL->Fill(j,mm_tracks->energy[i][j]);
                    }
                    if(i==0 && j>=120 && j%2==0) {
                        EstripL+=mm_tracks->energy[i][j];
                    }
                    if(i==72 && j>=120 && j%2==0) {
                        EstripR+=mm_tracks->energy[i][j];
                    }

                }
                //Si
                if(si_tracks->pixel[i][j]>0){
                    hMM_TrackAll->Fill(si_tracks->chanid,141+si_tracks->agetid);
                    //hMM_Time[goodevtcounter%16]->Fill(si_tracks->time[i][j]);
                    //hMM_Energy[goodevtcounter%16]->Fill(si_tracks->energy[i][j]);
                    for(int l=0;l<6;l++){
                        for(int m=0;m<6;m++){
                            //hMM_TrackAll->Fill(i+l,j+m);
                        }
                    }
                    //for(int k=0;k<si_tracks->coloridx[i][j];k++)
                    for(int k=0;k<1;k++)
                    {
                        for(int l=0;l<6;l++){
                            for(int m=0;m<6;m++){
                                //hMM_Track[goodevtcounter%16]->Fill(i+l,j+m);
                                //hMM_Track[goodevtcounter%16]->Fill(si_tracks->chanid,141+si_tracks->agetid);
                                hMM_TimevsPxIDX[goodevtcounter%16]->Fill(i+l,si_tracks->time[i][j]);
                                hMM_TimevsPxIDY[goodevtcounter%16]->Fill(j+m,si_tracks->time[i][j]);
                                hMM_EnergyvsPxIDY[goodevtcounter%16]->Fill(j+m,si_tracks->energy[i][j]);
                            }
                        }
                    }
                }
            }
        }

        if(EstripL>0){
            hMM_EstripL->Fill(mm_tracks->ICenergy,EstripL);
            hMM_EstripAll->Fill(mm_tracks->ICenergy,EstripL);
        }
        if(EstripR>0){
            hMM_EstripR->Fill(mm_tracks->ICenergy,EstripR);
            hMM_EstripAll->Fill(mm_tracks->ICenergy,EstripR);
        }
        //hMM_SumEnergyvsPxIDY[goodevtcounter%16]->Fill(130,mm_tracks->ICenergy);

        //if((EstripL>0 && EstripL<8000 && sumstrip[0]>700) || (EstripR>0 && EstripR<8000 && sumstrip[0]>700))
        /*
        for(int i=64;i<71;i++){
            for(int j=0;j<170;j++){
                if(mm_tracks->pixel[i][j]>0){
                    hMM_TrackPos[goodevtcounter%16]->Fill(mm_tracks->avgposx[j],mm_tracks->posy[i][j],20);
                }
            }
        }
        */
    }
}

void LKFrameBuilder::DrawTrack2pMode()
{
    int rgidx=0;
    int gtrackidx=0;
    double posz=0;
    double timez=0;
    double poszoffset=178; //tb=178 is the beam position=0
    cout<<mm_tracks->hasTrack<<"\ttrackhits"<<endl;
    //if(mm_tracks->hasTrack>0 && si_tracks->hasTrack>0 && si_tracks->hasICE>0)
    //if(mm_tracks->hasTrack>0 && si_tracks->hasTrack>0)
    if(mm_tracks->hasTrack>3)
    {
        gMM_TrackDecayPos[goodevtcounter%16]->SetPoint(gtrackidx++,-150,-300,-150);
        gMM_TrackDecayPos[goodevtcounter%16]->SetPoint(gtrackidx++,150,400,150);
        for(int i=0;i<150;i++){
            for(int j=0;j<170;j++){
                //Micromega
                if(mm_tracks->L1B[i][j]==1){
                    cout<<"I am drawing: "<<i<<"\t"<<j<<endl;
                }
                if(mm_tracks->pixel[i][j]>0){
                    if(i>63 && i<71) rgidx=0;
                    if(i<64) rgidx=1;
                    if(i>70) rgidx=2;
                    if(mm_tracks->posx[i][j]>-150 && mm_tracks->posx[i][j]<150){
                        posz = (poszoffset-(mm_tracks->time[i][j]-mm_tracks->decay[i][j]*256))*fTimePerBin*driftv;
                        timez = bucketmax - (mm_tracks->time[i][j]-mm_tracks->decay[i][j]*256);
                        hMM_TrackPosXYAll->Fill(mm_tracks->posx[i][j]+rx,mm_tracks->posy[i][j]+ry);
                        hMM_TrackPosYZAll->Fill(mm_tracks->posy[i][j]+ry,posz+rz);
                        hMM_TrackPosXZAll->Fill(mm_tracks->posx[i][j]+rx,posz+rz);
                        if(i<64 || i>69){
                            if(j%2!=0){ //chain
                                hMM_TrackXZ[goodevtcounter%16]->Fill(i,timez,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
                                hMM_TrackPosXZ[goodevtcounter%16]->Fill(mm_tracks->posx[i][j]+rx,posz+rz,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
                            }else{ //strip
                                hMM_TrackYZ[goodevtcounter%16]->Fill(j,timez,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
                                hMM_TrackPosYZ[goodevtcounter%16]->Fill(mm_tracks->posy[i][j]+ry,posz+rz,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
                            }
                        }else{
                            hMM_TrackXZ[goodevtcounter%16]->Fill(i,timez,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
                            hMM_TrackYZ[goodevtcounter%16]->Fill(j,timez,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
                            hMM_TrackPosXZ[goodevtcounter%16]->Fill(mm_tracks->posx[i][j]+rx,posz+rz,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
                            hMM_TrackPosYZ[goodevtcounter%16]->Fill(mm_tracks->posy[i][j]+ry,posz+rz,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
                        }
                        if(mm_tracks->L1B[i][j]==1){
                            cout<<"I am drawing: "<<i<<"\t"<<j<<endl;
                            //hMM_TrackDecay[goodevtcounter%16]->Fill(i,j,mm_tracks->decay[i][j]*100+mm_tracks->coloridx[i][j]);
                            gMM_TrackDecay[goodevtcounter%16]->SetPoint(gtrackidx,i,j,(int)(posz/1.7));
                            gMM_TrackDecayPos[goodevtcounter%16]->SetPoint(gtrackidx,mm_tracks->posx[i][j]+rx,mm_tracks->posy[i][j]+ry,posz+rz);
                            gtrackidx++;
                        }
                    }
                }
            }
        }
    }
}

void LKFrameBuilder::InitTrack()
{
    rpos = new TRandom();
    ResetTrack();
}

void LKFrameBuilder::DrawSiEvsCsIE()
{
    int rgidx=0;
    int detno=0;
    if(mm_tracks->hasTrack>3){
        for(int i=0;i<10;i++){
            for(int j=0;j<4;j++){
                if(si_tracks->pixel[i][j]>0){
                    if(i<2 && j<2) detno=9;
                    else if(i<2 && j>=2) detno=8;
                    else if((i>=2&&i<4) && j<2) detno=7;
                    else if((i>=2&&i<4) && j>=2) detno=6;
                    else if((i>=4&&i<6) && j<2) detno=5;
                    else if((i>=4&&i<6) && j>=2) detno=4;
                    else if((i>=6&&i<8) && j<2) detno=3;
                    else if((i>=6&&i<8) && j>=2) detno=2;
                    else if((i>=8&&i<10) && j<2) detno=1;
                    else if((i>=8&&i<10) && j>=2) detno=0;

                    if(si_tracks->hasCsIE[detno]>100){
                        if(detno==0||detno==1||detno==8||detno==9){
                            //hMM_SiEvsCsIE[0]->Fill(si_tracks->hasCsIE[detno],si_tracks->energy[i][j]);
                        }else if(detno==2||detno==3||detno==4||detno==6||detno==7){
                            //hMM_SiEvsCsIE[1]->Fill(si_tracks->hasCsIE[detno],si_tracks->energy[i][j]);
                            //if(si_tracks->energy[i][j]>400 && cut_dtinSiEvsCsIE->IsInside(si_tracks->hasCsIE[detno],si_tracks->energy[i][j])) goodsicsipevt=1;
                        }else{
                            //hMM_SiEvsCsIE[2]->Fill(si_tracks->hasCsIE[detno],si_tracks->energy[i][j]);
                            //if(si_tracks->energy[i][j]>400 && cut_pinSiEvsCsIE->IsInside(si_tracks->hasCsIE[detno],si_tracks->energy[i][j])) goodsicsipevt=1;
                            //if(si_tracks->energy[i][j]>400 && cut_dtinSiEvsCsIE->IsInside(si_tracks->hasCsIE[detno],si_tracks->energy[i][j])) goodsicsipevt=1;
                        }
                        //if(cut_inSiEvsCsIE->IsInside(si_tracks->hasCsIE[detno],si_tracks->energy[i][j])) goodsicsievt=1;
                        //if(si_tracks->energy[i][j]>400) goodsicsipevt=1;
                    }
                }
            }
        }
    }
}

void LKFrameBuilder::ResetdEvsE(){
    for(int i=0;i<3;i++){
        hMM_TrackCounter1[i]=0;
        hMM_TrackCounter2[i]=0;
        hMM_dE1[i]=0;
        hMM_dE2[i]=0;
        hMM_E[i]=0;
    }
}

void LKFrameBuilder::DrawdEvsE(){
    int rgidx=0;
    for(int j=112;j<128;j++){
        if(mm_tracks->sumenergy[j]>0){
            hMM_TrackCounter2[rgidx]++;
            if(j<116){
                hMM_TrackCounter1[rgidx]++;
                hMM_dE1[rgidx]+=mm_tracks->sumenergy[j];
            }
            hMM_dE2[rgidx]+=mm_tracks->sumenergy[j];
        }
    }
    for(int rgidx=1;rgidx<3;rgidx++){
        if(hMM_cornerfound[rgidx]>0){
            for(int j=hMM_MinPxY[rgidx];j<hMM_MaxPxY[rgidx];j++){
                if(j%2==0){
                    hMM_dE2[rgidx]+=mm_tracks->energy[hMM_MinPxX[rgidx]][j];
                    //cout << reventIdx << " " << hMM_MinPxX[rgidx] << " " << j << " " << mm_tracks->energy[hMM_MinPxX[rgidx]][j] << endl;
                }
            }
            for(int i=hMM_MinPxX[rgidx];i<hMM_MaxPxX[rgidx];i++){
                if(hMM_MinPxY[rgidx]%2==0){
                    hMM_dE2[rgidx]+=mm_tracks->energy[i][hMM_MinPxY[rgidx]+1];
                    //cout << reventIdx << " " << i << " " << hMM_MinPxY[rgidx] << " "  << mm_tracks->energy[i][hMM_MinPxY[rgidx]+1] << endl;
                }else{
                    hMM_dE2[rgidx]+=mm_tracks->energy[i][hMM_MinPxY[rgidx]];
                    //cout << reventIdx << " " << i << " " << hMM_MinPxY[rgidx] << " "  << mm_tracks->energy[i][hMM_MinPxY[rgidx]] << endl;
                }
            }
            //hMM_TrackCounter2[rgidx] = hMM_MaxPxY[rgidx]-hMM_MinPxY[rgidx];
        }
    }

    for(rgidx=0;rgidx<3;rgidx++){
        if(hMM_E[rgidx]>0 && hMM_dE2[rgidx]>0 && rgidx>0 && hMM_cornerfound[rgidx]>0){
            double angle = atan(hMM_Slope[rgidx]);
            hMM_dE1[rgidx]=hMM_dE1[rgidx]*TMath::Cos(angle)/hMM_TrackCounter1[rgidx];
            //hMM_dE2[rgidx]=hMM_dE2[rgidx]/hMM_TrackCounter2[rgidx];
            //hMM_dE2[rgidx]=hMM_dE2[rgidx]*TMath::Cos(angle)/hMM_TrackCounter2[rgidx];
            hMM_dE2[rgidx]=hMM_dE2[rgidx]/TMath::Sqrt(TMath::Power(hMM_MaxY[rgidx]-hMM_MinY[rgidx],2)+TMath::Power(hMM_MaxX[rgidx]-hMM_MinX[rgidx],2));
        }else if(hMM_E[rgidx]>0 && hMM_dE2[rgidx]>0){
            hMM_dE1[rgidx]=hMM_dE1[rgidx]/hMM_TrackCounter1[rgidx];
            hMM_dE2[rgidx]=hMM_dE2[rgidx]/hMM_TrackCounter2[rgidx];
        }
    }

    for(rgidx=0;rgidx<3;rgidx++){
        if(hMM_E[rgidx]>0 && hMM_dE2[rgidx]>0){
            hMM_TrackdEvsEALL[rgidx]->Fill(hMM_E[rgidx],hMM_dE2[rgidx]);
            hMM_TrackdE1vsSiE->Fill(hMM_E[rgidx],hMM_dE1[rgidx]);
            hMM_TrackdE1vsdE2->Fill(hMM_dE2[rgidx],hMM_dE1[rgidx]);
            //cout << rgidx << " " << hMM_E[rgidx] << " " << hMM_dE2[rgidx] << endl;;
        }
    }
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            if(hMM_E[i]>0 && hMM_dE2[j]>0){
                hMM_TrackdEvsEALL[i]->Fill(hMM_E[i],hMM_dE2[j]);
            }
        }
    }
}

void LKFrameBuilder::ResetTrack(){
    rx = 0;
    ry = 0;
    rz = 0;
    rt = 0;
    mm_tracks->ICenergy=0;
    mm_tracks->Sienergy=0;
    mm_tracks->SiX=0;
    mm_tracks->SiY=0;
    mm_tracks->SiZ=0;
    mm_tracks->hasTrack=0;
    mm_tracks->hasTrackChain=0;
    mm_tracks->hasTrackStrip=0;
    mm_tracks->hasProjectile=0;
    mm_tracks->hasEjectile=0;
    mm_tracks->hasRecoil=0;
    mm_tracks->hasDecay=false;
    for(int i=0;i<4;i++) mm_tracks->sum2penergy[i] = 0;
    si_tracks->hasTrack=0;
    si_tracks->hasFWC=0;
    si_tracks->hasFWL=0;
    si_tracks->hasFWR=0;
    si_tracks->hasX6L=0;
    si_tracks->hasX6R=0;
    si_tracks->hasX6BL=0;
    si_tracks->hasX6BR=0;
    si_tracks->hasProjectile=0;
    si_tracks->hasEjectile=0;
    si_tracks->hasRecoil=0;
    si_tracks->hasICE=0;
    si_tracks->hasICT=0;
    for(int i=0;i<9;i++) si_tracks->hasCsIE[i]=0;
    for(int i=0;i<9;i++) si_tracks->hasCsIT[i]=0;
    si_tracks->agetid=0;
    si_tracks->chanid=0;
    for(int i=0;i<150;i++){
        for(int j=0;j<170;j++){
            mm_tracks->L1B[i][j]=0;
            mm_tracks->decay[i][j]=0;
            mm_tracks->pixel[i][j]=0;
            mm_tracks->energy[i][j]=0;
            mm_tracks->time[i][j]=0;
            mm_tracks->coloridx[i][j]=0;
            si_tracks->pixel[i][j]=0;
            si_tracks->energy[i][j]=0;
            si_tracks->time[i][j]=0;
            si_tracks->coloridx[i][j]=0;
        }
    }
    for(int j=0;j<170;j++){
        mm_tracks->avgposx[j] = 0;
        mm_tracks->sumenergy[j] = 0;
        posenergysum[j]=0;
        energysum[j]=0;
        sumcounter[j]=0;
    }
}

void LKFrameBuilder::ResetTrackHist()
{
    gMM_TrackDecay[goodevtcounter%16]->Clear();
    gMM_TrackDecayPos[goodevtcounter%16]->Clear();
    hMM_Track[goodevtcounter%16]->Reset();
    if(enable2pmode==1) hMM_TrackDecay[goodevtcounter%16]->Reset();
    if(enable2pmode==1) hMM_TrackDecay2[goodevtcounter%16]->Reset();
    hMM_TrackvsE[goodevtcounter%16]->Reset();
    hMM_TimevsPxIDX[goodevtcounter%16]->Reset();
    hMM_TimevsPxIDY[goodevtcounter%16]->Reset();
    hMM_TrackXZ[goodevtcounter%16]->Reset();
    hMM_EnergyvsPxIDY[goodevtcounter%16]->Reset();
    hMM_SumEnergyvsPxIDY[goodevtcounter%16]->Reset();
    hMM_TimevsPxIDYPos[goodevtcounter%16]->Reset();
    hMM_TimevsPxIDXPos[goodevtcounter%16]->Reset();
    hMM_TrackYZ[goodevtcounter%16]->Reset();
    hMM_TrackPos[goodevtcounter%16]->Reset();
    hMM_TrackPosXY[goodevtcounter%16]->Reset();
    hMM_TrackPosYZ[goodevtcounter%16]->Reset();
    //hMM_TrackdEvsE[goodevtcounter%16][rgidx]->Reset();
    for(int rgidx=0;rgidx<3;rgidx++){
        hMM_TrackPosHough[goodevtcounter%16][rgidx]->Reset();
        hMM_TimevsPxIDXPosHough[goodevtcounter%16][rgidx]->Reset();
        hMM_TimevsPxIDYPosHough[goodevtcounter%16][rgidx]->Reset();
    }

    /*
    for(int i=0;i<150;i++){
        for(int j=0;j<170;j++){
            hMM_Track[goodevtcounter%16]->SetBinContent(i,j,0);
            if(enable2pmode==1) hMM_TrackDecay[goodevtcounter%16]->SetBinContent(i,j,0);
            hMM_TrackvsE[goodevtcounter%16]->SetBinContent(i,j,0);
            //hMM_TrackvsE[gatedevtcounter%16]->SetBinContent(i,j,0);
        }
        for(int j=0;j<512;j++){
            hMM_TimevsPxIDX[goodevtcounter%16]->SetBinContent(i,j,0);
            hMM_TimevsPxIDY[goodevtcounter%16]->SetBinContent(i,j,0);
            hMM_TrackXZ[goodevtcounter%16]->SetBinContent(i,j,0);
        }
        for(int j=0;j<4000;j++){
            hMM_EnergyvsPxIDY[goodevtcounter%16]->SetBinContent(i,j,0);
        }
        hMM_SumEnergyvsPxIDY[goodevtcounter%16]->SetBinContent(i,0);
        //hMM_SumEnergyvsPxIDY[gatedevtcounter%16]->SetBinContent(i,0);
    }
    for(int i=0;i<700;i++){
        for(int j=0;j<512;j++){
            hMM_TimevsPxIDYPos[goodevtcounter%16]->SetBinContent(i,j,0);
        }
    }
    for(int i=0;i<300;i++){
        for(int j=0;j<512;j++){
            hMM_TimevsPxIDXPos[goodevtcounter%16]->SetBinContent(i,j,0);
        }
    }
    for(int i=0;i<170;i++){
        for(int j=0;j<512;j++){
            hMM_TrackYZ[goodevtcounter%16]->SetBinContent(i,j,0);
        }
    }
    for(int i=0;i<300;i++){
        for(int j=0;j<700;j++){
            hMM_TrackPos[goodevtcounter%16]->SetBinContent(i,j,0);
            hMM_TrackPosXY[goodevtcounter%16]->SetBinContent(i,j,0);
        }
    }
    for(int i=0;i<700;i++){
        for(int j=0;j<400;j++){
            hMM_TrackPosYZ[goodevtcounter%16]->SetBinContent(i,j,0);
        }
    }
    for(int rgidx=0;rgidx<3;rgidx++){
        for(int i=0;i<512;i++){
            for(int j=0;j<512;j++){
                //hMM_TrackdEvsE[goodevtcounter%16][rgidx]->SetBinContent(i,j,0);
            }
        }
        for(int i=0;i<1800;i++){
            for(int j=0;j<600;j++){
                hMM_TrackPosHough[goodevtcounter%16][rgidx]->SetBinContent(i,j,0);
            }
            for(int j=0;j<1000;j++){
                hMM_TimevsPxIDXPosHough[goodevtcounter%16][rgidx]->SetBinContent(i,j,0);
                hMM_TimevsPxIDYPosHough[goodevtcounter%16][rgidx]->SetBinContent(i,j,0);
            }
        }
    }
    */

    TString title0 = Form("[D2PTime=%d usec, FrameNo~%d, EventNo=%d, goodevtcounter=%d]",int(rd2ptime/1000),(reventIdx-FirsteventIdx+2),reventIdx,goodevtcounter);
    hMM_Track[goodevtcounter%16]->SetTitle(Form("Single Event Track by channels;Cell ID X %s;Cell ID Y",title0.Data()));
    hMM_TrackvsE[goodevtcounter%16]->SetTitle(Form("Single Event Track vs E;Cell ID X;Cell ID Y [EventNo=%d, goodevtcounter=%d]",reventIdx,goodevtcounter));
    hMM_TrackDecay[goodevtcounter%16]->SetTitle(Form("Single Event Decay Track by channels;Cell ID X %s;Cell ID Y",title0.Data()));
    hMM_TrackDecay2[goodevtcounter%16]->SetTitle(Form("Single Event Decay Track by channels;Cell ID X %s;Cell ID Y",title0.Data()));
    hMM_TrackPos[goodevtcounter%16]->SetTitle(Form("Single Event Track by Position;Pos X (mm) %s;Pos Y (mm)",title0.Data()));
    hMM_TrackPosXY[goodevtcounter%16]->SetTitle(Form("Single Event Track Pos Y vs Pos X;Pos X (mm) %s;Pos Y (mm)",title0.Data()));
    hMM_TimevsPxIDX[goodevtcounter%16]->SetTitle(Form("Single Event Time vs Pixel X;Cell ID X %s;Time Bucket",title0.Data()));
    hMM_TimevsPxIDXPos[goodevtcounter%16]->SetTitle(Form("Single Event Time vs Pos X;Pos X (mm) %s;Time Bucket",title0.Data()));
    hMM_TimevsPxIDY[goodevtcounter%16]->SetTitle(Form("Single Event Time vs Pixel Y;Cell ID Y %s;Time Bucket",title0.Data()));
    hMM_TimevsPxIDYPos[goodevtcounter%16]->SetTitle(Form("Single Event Time vs Pos Y;Pos Y (mm) %s;Time Bucket",title0.Data()));
    hMM_EnergyvsPxIDY[goodevtcounter%16]->SetTitle(Form("Single Event Energy vs Pixel Y;Cell ID Y %s;Energy",title0.Data()));
    hMM_SumEnergyvsPxIDY[goodevtcounter%16]->SetTitle(Form("Single Event Sum Energy vs Pixel Y %s;Cell ID Y;Energy",title0.Data()));
    hMM_SumEnergyvsPxIDY[goodevtcounter%16]->SetFillStyle(3001);
    hMM_SumEnergyvsPxIDY[goodevtcounter%16]->SetFillColor(kMagenta+3);
    hMM_SumEnergyvsPxIDY[goodevtcounter%16]->SetLineColor(1);
    hMM_SumEnergyvsPxIDY[goodevtcounter%16]->SetOption("hist");
    //hMM_SumEnergyvsPxIDY[gatedevtcounter%16]->SetTitle(Form("Single Event Sum Energy vs Pixel Y;Cell ID Y;Energy [EventNo=%d]",reventIdx));
    //hMM_SumEnergyvsPxIDY[gatedevtcounter%16]->SetFillStyle(3001);
    //hMM_SumEnergyvsPxIDY[gatedevtcounter%16]->SetFillColor(kMagenta+3);
    //hMM_SumEnergyvsPxIDY[gatedevtcounter%16]->SetLineColor(1);
    //hMM_SumEnergyvsPxIDY[gatedevtcounter%16]->SetOption("hist");

    hMM_Time[goodevtcounter%16]->Reset();
    hMM_Energy[goodevtcounter%16]->Reset();

    //for(int j=0;j<512;j++){
    //    hMM_Time[goodevtcounter%16]->SetBinContent(j,0);
    //}
    //for(int j=0;j<4000;j++){
    //    hMM_Energy[goodevtcounter%16]->SetBinContent(j,0);
    //}
    hMM_Time[goodevtcounter%16]->SetTitle(Form("Single Event Time;Time Bucket [FrameNo~%d, EventNo=%d, goodevtcounter=%d, D2PTime=%d usec]",(reventIdx-FirsteventIdx+2),reventIdx,goodevtcounter,int(rd2ptime/1000)));
    hMM_Energy[goodevtcounter%16]->SetTitle(Form("Single Event Energy;Energy [FrameNo~%d, EventNo=%d, goodevtcounter=%d, D2PTime=%d usec]",(reventIdx-FirsteventIdx+2),reventIdx,goodevtcounter,int(rd2ptime/1000)));
    //for(int j=0;j<170;j++){
    //    pMM_Track[goodevtcounter%16]->SetBinContent(j,0);
    //}
}

void LKFrameBuilder::decodeCoBoTopologyFrame(mfm::Frame& frame) {
    if (frame.header().frameType() == 0x7){
        // Print meta-data
        //cout << "Topology:" << endl;
        // Skip 8 Byte header
        istream & stream = frame.serializer().inputStream(8u);
        uint8_t value;
        stream.read(reinterpret_cast<char *>(&value), 1lu);
        //cout << " coboIdx=" << (uint16_t) value << ", ";
        stream.read(reinterpret_cast<char *>(&value), 1lu);
        //cout << " asadMask=" << hex << showbase << (uint16_t) value << dec << ", ";
        stream.read(reinterpret_cast<char *>(&value), 1lu);
        //cout << " 2pMode=" << boolalpha << (bool) value << endl;
    }
}

void LKFrameBuilder::InitMutantScaler()
{
    printed=0;
    prevmidx=0;
    mevtidx.resize(0);
    md2ptime.resize(0);
    mtstmp.resize(0);
    scaler1=0;
    scaler1start=0;
    scaler1end=0;
    scaler2=0;
    scaler2start=0;
    scaler2end=0;
    scaler3=0;
    scaler3start=0;
    scaler3end=0;
    scaler4=0;
    scaler4start=0;
    scaler4end=0;
    scaler5=0;
    scaler5start=0;
    scaler5end=0;
    remove("scalers.txt");
}

void LKFrameBuilder::SetReadMode(int flag){
    readmode = flag;
}

int LKFrameBuilder::GetReadMode(){
    return readmode;
}

void LKFrameBuilder::SetReadType(int flag){
    readtype = flag;
}

void LKFrameBuilder::SetScaler(int flag){
    enablescaler = flag;
}

void LKFrameBuilder::Set2pMode(int flag){
    enable2pmode = flag;
}

void LKFrameBuilder::SetHistMode(){
    enablehist = 1;
}

void LKFrameBuilder::SetUpdateSpeed(int flag){
    enableupdatefast = flag;
}


void LKFrameBuilder::decodeMuTanTFrame(mfm::Frame & frame)
{
    mutantcounter++;
    if (frame.header().frameType() == 0x8)
    {
        try
        {
            if(enable2pmode==1){
                mevtidx.push_back(frame.headerField(14u, 4u).value< uint32_t >());
                md2ptime.push_back(frame.headerField(60u, 4u).value< uint32_t >()*10);
                mtstmp.push_back(frame.headerField(8u, 6u).value< uint64_t >());
                cout << "decoded: " << mevtidx.at(mevtidx.size()-1) << " " << md2ptime.at(md2ptime.size()-1) << " " << mtstmp.at(mtstmp.size()-1) << endl;
            }else{
                scaler1 = frame.headerField(48u, 4u).value< uint32_t >();
                scaler2 = frame.headerField(52u, 4u).value< uint32_t >();
                scaler3 = frame.headerField(56u, 4u).value< uint32_t >();
                scaler4 = frame.headerField(40u, 4u).value< uint32_t >();
                scaler5 = frame.headerField(44u, 4u).value< uint32_t >();
                scevent = frame.headerField(24u, 4u).value< uint32_t >();
                scevtidx = frame.headerField(14u, 4u).value< uint32_t >();
                sctstmp = frame.headerField(8u, 6u).value< uint64_t >();
                //if(scaler3>0 && (int)scaler3%1000000==0)
                if(mutantcounter==100)
                {
                    scaler1end = scaler1;
                    scaler2end = scaler2;
                    scaler3end = scaler3;
                    cout << " Scaler 1 Rate = " << 100*(scaler1end-scaler1start)/(scaler3end-scaler3start) << ", " << (scaler1end-scaler1start) << " " << (scaler3end-scaler3start) << " " << scaler1end << " " << scaler1start << " " << scaler3end << " " << scaler3start << " or " << (100*scaler1end/scaler3end) << endl;
                    cout << " Scaler 2 Rate(Live) = " << 100*(scaler2end-scaler2start)/(scaler3end-scaler3start) << ", " << scaler2end << " " << scaler2start << " " << scaler3end << " " << scaler3start << " or " << (100*scaler2end/scaler3end) << endl;
                    if(enablescaler==1){
                        scalerout.open("scalers.txt", std::ofstream::out|std::ofstream::app);
                        scalerout << sctstmp << " " << scevent << " " << scevtidx << " " << scaler1 << " " << scaler2 << " " << scaler3 << endl;
                        scalerout.close();
                    }
                    mutantcounter=0;
                }
                //else if(scaler3>0 && (int)scaler3%1000000==1)
                else if(mutantcounter==1)
                {
                    scaler1start = scaler1;
                    scaler2start = scaler2;
                    scaler3start = scaler3;
                }
                if(enablescaler==0){
                    printed++;
                    if(printed>15){
                        printed=0;
                        cout << "Mutant: ";
                        cout << " TSTMP=" << frame.headerField(8u, 6u).value< uint64_t >() << ", ";
                        cout << " EVT_NO=" << frame.headerField(14u, 4u).value< uint32_t >() << ", ";
                        cout << " TRIG=" << hex << showbase << frame.headerField(18u, 2u).value< uint16_t >() << dec << ", ";
                        cout << " MULT_A="  << frame.headerField(20u, 2u).value< uint16_t >() << ", ";
                        cout << " MULT_B="  << frame.headerField(22u, 2u).value< uint16_t >() << ", ";
                        cout << " L0_EVT="  << frame.headerField(24u, 4u).value< uint32_t >() << ", ";
                        cout << " L1A_EVT="  << frame.headerField(28u, 4u).value< uint32_t >() << ", ";
                        cout << " L1B_EVT="  << frame.headerField(32u, 4u).value< uint32_t >() << ", ";
                        cout << " L2_EVT="  << frame.headerField(36u, 4u).value< uint32_t >() << ", ";
                        cout << " SCA1="  << frame.headerField(48u, 4u).value< uint32_t >() << ", ";
                        cout << " SCA2="  << frame.headerField(52u, 4u).value< uint32_t >() << ", ";
                        cout << " SCA3="  << frame.headerField(56u, 4u).value< uint32_t >() << ", ";
                        cout << " SCA4="  << frame.headerField(40u, 4u).value< uint32_t >() << ", ";
                        cout << " SCA5="  << frame.headerField(44u, 4u).value< uint32_t >() << ", ";
                        cout << " D2P="  << frame.headerField(60u, 4u).value< uint32_t >() << endl;
                        //cout << frame.headerField(8u, 6u).value< uint64_t >() << " ";
                        //cout << frame.headerField(14u, 4u).value< uint32_t >() << " ";
                        //cout << frame.headerField(24u, 4u).value< uint32_t >() << " ";
                        //cout << frame.headerField(28u, 4u).value< uint32_t >() << " ";
                        //cout << frame.headerField(32u, 4u).value< uint32_t >() << " ";
                        //cout << frame.headerField(36u, 4u).value< uint32_t >() << " ";
                        cout << scaler1 << " " << scaler2 << " " << scaler3 << ", IC_Rate=" << Form("%5.2f",(100*scaler1/scaler3)) << "pps, Live=" << Form("%5.2f",((100*scaler2)/scaler3)) << "\%" << endl;
                        //double beamrate = scaler3/2000;
                        //cout << "beamrate: " << scevent << " " << beamrate << endl;
                        /*
                           cout << "Mutant:" << endl;
                           cout << " TIMESTAMP=" << frame.headerField(8u, 6u).value< uint32_t >() << endl;
                           cout << " EVENT_NUMBER=" << frame.headerField(14u, 4u).value< uint32_t >() << endl;
                           cout << " TRIGGER_INFO=" << hex << showbase << frame.headerField(18u, 2u).value< uint16_t >() << dec << endl;
                           cout << " MULTIPLICITY_A_MEM="  << frame.headerField(20u, 2u).value< uint16_t >() << endl;
                           cout << " MULTIPLICITY_B_MEM="  << frame.headerField(22u, 2u).value< uint16_t >() << endl;
                           cout << " L0_EVT_COUNTER="  << frame.headerField(24u, 4u).value< uint32_t >() << endl;
                           cout << " L1A_EVT_COUNTER="  << frame.headerField(28u, 4u).value< uint32_t >() << endl;
                           cout << " L1B_EVT_COUNTER="  << frame.headerField(32u, 4u).value< uint32_t >() << endl;
                           cout << " L2_EVT_COUNTER="  << frame.headerField(36u, 4u).value< uint32_t >() << endl;
                           cout << " SCALER1_REG="  << frame.headerField(48u, 4u).value< uint32_t >() << endl;
                           cout << " SCALER2_REG="  << frame.headerField(52u, 4u).value< uint32_t >() << endl;
                           cout << " SCALER3_REG="  << frame.headerField(56u, 4u).value< uint32_t >() << endl;
                           cout << " SCALER4_REG="  << frame.headerField(40u, 4u).value< uint32_t >() << endl;
                           cout << " SCALER5_REG="  << frame.headerField(44u, 4u).value< uint32_t >() << endl;
                           cout << " D2P_TIME="  << frame.headerField(60u, 4u).value< uint32_t >() << endl;
                         */
                    }
                }
            }
        } catch (const mfm::Exception &)
        {}
    }
}

void LKFrameBuilder::RootWOpenFile(string & outputFileName, string & outputTreeName)
{
    return;
    cout<<"FIRST EVENT?:\t"<<IsFirstevent<<endl;

    if(IsFirstevent){
        fOutputFile = new TFile(outputFileName.c_str(),"recreate");
        fOutputTree = new TTree(outputTreeName.c_str(),"Experimental Event Data");
    }
    else{
        cout<<"Update rootwopenfile "<<outputFileName.c_str()<<endl;

        fOutputFile = new TFile(outputFileName.c_str(),"update");
        fOutputTree = (TTree *)fOutputFile->Get(outputTreeName.c_str());
    }
    if(enableupdatefast==0){
        // void TTree::SetAutoSave	(	Long64_t 	autos = -300000000	)	
        //
        // In case of a program crash, it will be possible to recover the data in the tree up to the last AutoSave point.
        // This function may be called before filling a TTree to specify when the branch buffers and TTree header are flushed to disk as part of TTree::Fill().
        // The default is -300000000, ie the TTree will write data to disk once it exceeds 300 MBytes.
        // CASE 1: If fAutoSave is positive the watermark is reached when a multiple of fAutoSave entries have been filled.
        // CASE 2: If fAutoSave is negative the watermark is reached when -fAutoSave bytes can be written to the file.
        // CASE 3: If fAutoSave is 0, AutoSave() will never be called automatically as part of TTree::Fill().
        fOutputTree->SetAutoSave(-1000000);
    }
    //cout << Form("Created new TTree named %s in %s (new=%d)",outputTreeName.c_str(),outputFileName.c_str(),IsFirstevent) << endl;
}

void LKFrameBuilder::RootWInit()
{
    //XXX
    return;
    //cout<<"Root W Init"<<endl;
    wGETMul = 0;
    wGETHit = 0;
    wGETEventIdx = 0;
    wGETD2PTime = 0;
    wGETTimeStamp = 0;
    for(int i=0;i<4352;i++){
        wGETFrameNo[i] = 0;
        wGETDecayNo[i] = 0;
        wGETEnergy[i] = 0;
        wGETCobo[i] = 0;
        wGETAsad[i] = 0;
        wGETAget[i] = 0;
        wGETChan[i] = 0;
        for(int j=0;j<512;j++){
            wGETWaveformX[i][j] = 0;
            wGETWaveformY[i][j] = 0;
        }
    }
    if(readmode==1){
        fOutputTree->Branch("mmMul",&wGETMul,"mmMul/I");
        fOutputTree->Branch("mmHit",&wGETHit,"mmHit/I");
        fOutputTree->Branch("mmEventIdx",&wGETEventIdx,"mmEventIdx/I");
        fOutputTree->Branch("mmD2PTime",&wGETD2PTime,"mmD2PTime/I");
        fOutputTree->Branch("mmTimeStamp",&wGETTimeStamp,"mmTimeStamp/I");
        fOutputTree->Branch("mmFrameNo",wGETFrameNo,"mmFrameNo[mmMul]/I");
        fOutputTree->Branch("mmDecayNo",wGETDecayNo,"mmDecayNo[mmMul]/I");
        fOutputTree->Branch("mmCobo",wGETCobo,"mmCobo[mmMul]/I");
        fOutputTree->Branch("mmAsad",wGETAsad,"mmAsad[mmMul]/I");
        fOutputTree->Branch("mmAget",wGETAget,"mmAget[mmMul]/I");
        fOutputTree->Branch("mmChan",wGETChan,"mmChan[mmMul]/I");
        fOutputTree->Branch("mmTime",wGETTime,"mmTime[mmMul]/F");
        fOutputTree->Branch("mmEnergy",wGETEnergy,"mmEnergy[mmMul]/F");
        fOutputTree->Branch("mmWaveformX",wGETWaveformX,Form("mmWaveformX[mmMul][%d]/I",bucketmax)); // To record the waveform timebucket
        fOutputTree->Branch("mmWaveformY",wGETWaveformY,Form("mmWaveformY[mmMul][%d]/I",bucketmax)); // To record the waveform samplevalue
    }
}

void LKFrameBuilder::RootWReset()
{
    for(int i=0;i<=wGETMul;i++){
        wGETFrameNo[i] = 0;
        wGETDecayNo[i] = 0;
        wGETTime[i] = 0;
        wGETEnergy[i] = 0;
        wGETCobo[i] = 0;
        wGETAsad[i] = 0;
        wGETAget[i] = 0;
        wGETChan[i] = 0;
        for(int j=0;j<bucketmax;j++){
            wGETWaveformX[i][j] = 0;
            wGETWaveformY[i][j] = 0;
        }
    }
}

void LKFrameBuilder::RootWConvert()
{
    //XXX
    UInt_t frameIdx = waveforms->frameIdx;
    UInt_t decayIdx = waveforms->decayIdx;
    //  if(decayIdx==1) cout<<"RootWConvert L1B"<<endl;
    UInt_t coboIdx = waveforms->coboIdx;
    UInt_t asadIdx = waveforms->asadIdx;

    //auto header = (GETFrameHeader *) fFrameHeader -> ConstructedAt(0);
    //header -> SetFrameIdx(frameIdx);
    //header -> SetDecayIdx(decayIdx);


    int countPad = 0;
    for(UInt_t asad=0; asad<maxasad; asad++) {
        for(UInt_t aget=0; aget<4; aget++) {
            if(!waveforms->hasHit[asad*4+aget]) continue; // Skip no fired agets.
            for(UInt_t chan=0; chan<68; chan++) {

                if(!waveforms->hasSignal[asad*4+aget][chan]) continue; // Skip no fired channels.
                if(coboIdx>=0){
                    if(readmode==1){

                        auto channel = (GETChannel *) fChannelArray -> ConstructedAt(countPad++);
                        channel -> SetCobo(coboIdx);
                        channel -> SetAsad(asad);
                        channel -> SetAget(aget);
                        channel -> SetChan(chan);
                        channel -> SetTime(0);
                        channel -> SetEnergy(0);
                        channel -> SetWaveform(waveforms->waveform[asad*4+aget][chan]);

                        /*
                        wGETFrameNo[wGETMul] = frameIdx;
                        wGETDecayNo[wGETMul] = decayIdx;
                        wGETTime[wGETMul] = 0;
                        wGETEnergy[wGETMul] = 0;
                        wGETCobo[wGETMul] = coboIdx;
                        wGETAsad[wGETMul] = asad;
                        wGETAget[wGETMul] = aget;
                        wGETChan[wGETMul] = chan;
                        for(int i=0;i<bucketmax;i++){
                            wGETWaveformX[wGETMul][i] = i;
                            wGETWaveformY[wGETMul][i] = waveforms->waveform[asad*4+aget][chan][i];
                        }
                        //cout << weventIdx << " " << wGETMul << " " << coboIdx << " " << asad << " " << aget << " " << chan << endl;
                        //if(weventIdx==1412) cout << wGETMul << " " << wGETFrameNo[wGETMul] << " " << wGETCobo[wGETMul] << " " << wGETAsad[wGETMul] << " " << wGETAget[wGETMul] << " " << wGETChan[wGETMul] << endl;
                        */
                        wGETMul++;
                        if(chan!=11&&chan!=22&&chan!=45&&chan!=56) wGETHit++;
                    }
                }
            }
        }
    }
}

void LKFrameBuilder::RootWriteEvent(){
    if(wGETMul>0){
        //cout << "Writing data: " << wGETEventIdx << " " << wGETMul << endl;
        fOutputFile->cd();
        fOutputTree->Fill();
        if(enableupdatefast==1){
            if(framecounter%16==0){
                fOutputTree->Write();
//  XXX
//  void TDirectoryFile::SaveSelf(Bool_t force=kFALSE)
//
//  Save Directory keys and header.
//  If the directory has been modified (fModified set), write the keys and the directory header.
//  This function assumes the cd is correctly set.
//
//  It is recommended to use this function in the following situation:
//      Assume a process1 using a directory in Update mode
//  - New objects or modified objects have been written to the directory.
//  - You do not want to close the file.
//  - You want your changes be visible from another process2 already connected to this directory in read mode.
//  - Call this function.
//  - In process2, use TDirectoryFile::ReadKeys to refresh the directory.
                fOutputFile->SaveSelf();
                if((fOutputFile->GetSize())>5000000){
                    enableupdatefast=0;
                }
                //cout << (int)(framecounter/16) << " times written, Size = " << fOutputFile->GetSize() << endl;
            }
        }else{
            if(framecounter%128==0){
                fOutputTree->Write();
                fOutputFile->SaveSelf();
            }
        }
        framecounter++;
        wGETMul = 0;
        wGETHit = 0;
    }
}

void LKFrameBuilder::RootWCloseFile(){
    //fOutputFile->cd();
    fOutputFile->Close();
    IsFirstevent = true;
}

bool LKFrameBuilder::RootRRecovered(){
    return fInputFile->TestBit(TFile::kRecovered);
}

UInt_t LKFrameBuilder::GetRootFileSize(){
    fInputFileSize = fInputFile->GetSize();
    return fInputFileSize;
}

void LKFrameBuilder::RootROpenFile(string & inputFileName, string & inputTreeName){
    size_t f = inputFileName.find("online");
    if(f!=std::string::npos) gErrorIgnoreLevel = kFatal;
    fInputFile = new TFile(inputFileName.c_str(),"read");
    fInputTree = (TTree*)fInputFile->Get(inputTreeName.c_str());
    if(fInputTree == nullptr){
        TKey *keyobj = (TKey *)fInputFile->GetKey(inputTreeName.c_str());
        short cycleno = keyobj->GetCycle();
        inputTreeName += ";"+ to_string(cycleno-1);
        fInputTree = (TTree*)fInputFile->Get(inputTreeName.c_str());
        forcereadtree =1;
    }
    fInputFileSize = fInputFile->GetSize();
    if(f!=std::string::npos) gErrorIgnoreLevel = kError;
}

void LKFrameBuilder::RootRInit(){
    if(readmode==0){
        fNumberEvents = 1;
    }else{
        fInputTree->SetBranchAddress("mmMul",&rGETMul);
        fInputTree->SetBranchAddress("mmHit",&rGETHit);
        fInputTree->SetBranchAddress("mmEventIdx",&rGETEventIdx);
        fInputTree->SetBranchAddress("mmD2PTime",&rGETD2PTime);
        fInputTree->SetBranchAddress("mmTimeStamp",&rGETTimeStamp);
        fInputTree->SetBranchAddress("mmFrameNo",rGETFrameNo);
        fInputTree->SetBranchAddress("mmDecayNo",rGETDecayNo);
        fInputTree->SetBranchAddress("mmCobo",rGETCobo);
        fInputTree->SetBranchAddress("mmAsad",rGETAsad);
        fInputTree->SetBranchAddress("mmAget",rGETAget);
        fInputTree->SetBranchAddress("mmChan",rGETChan);
        fInputTree->SetBranchAddress("mmTime",rGETTime);
        fInputTree->SetBranchAddress("mmEnergy",rGETEnergy);
        fInputTree->SetBranchAddress("mmWaveformX",rGETWaveformX);
        fInputTree->SetBranchAddress("mmWaveformY",rGETWaveformY);
        /*
           TBranch *branch;
           branch = fInputTree->GetBranch("mmMul");
           branch->SetAutoDelete(kTRUE);
           branch = fInputTree->GetBranch("mmHit");
           branch->SetAutoDelete(kTRUE);
           branch = fInputTree->GetBranch("mmEventIdx");
           branch->SetAutoDelete(kTRUE);
           branch = fInputTree->GetBranch("mmD2PTime");
           branch->SetAutoDelete(kTRUE);
           branch = fInputTree->GetBranch("mmFrameNo");
           branch->SetAutoDelete(kTRUE);
           branch = fInputTree->GetBranch("mmDecayNo");
           branch->SetAutoDelete(kTRUE);
           branch = fInputTree->GetBranch("mmCobo");
           branch->SetAutoDelete(kTRUE);
           branch = fInputTree->GetBranch("mmAsad");
           branch->SetAutoDelete(kTRUE);
           branch = fInputTree->GetBranch("mmAget");
           branch->SetAutoDelete(kTRUE);
           branch = fInputTree->GetBranch("mmChan");
           branch->SetAutoDelete(kTRUE);
           branch = fInputTree->GetBranch("mmTime");
           branch->SetAutoDelete(kTRUE);
           branch = fInputTree->GetBranch("mmEnergy");
           branch->SetAutoDelete(kTRUE);
           branch = fInputTree->GetBranch("mmWaveformX");
           branch->SetAutoDelete(kTRUE);
           branch = fInputTree->GetBranch("mmWaveformY");
           branch->SetAutoDelete(kTRUE);
         */
        rGETMul = 0;
        rGETHit = 0;
        rGETEventIdx = 0;
        rGETD2PTime = 0;
        rGETTimeStamp = 0;
        goodx6csievt=0;
        for(int i=0;i<4352;i++){
            rGETFrameNo[i] = 0;
            rGETDecayNo[i] = 0;
            rGETL1Aflag[i] = 0;
            rGETEnergy[i] = 0;
            rGETTime[i] = 0;
            rGETCobo[i] = 0;
            rGETAsad[i] = 0;
            rGETAget[i] = 0;
            rGETChan[i] = 0;
            for(int j=0;j<512;j++){
                rGETWaveformX[i][j] = 0;
                rGETWaveformY[i][j] = 0;
            }
        }
        fNumberEvents = fInputTree->GetEntries();
    }
    RootRInitWaveforms();
}

void LKFrameBuilder::RootRInitWaveforms() {
    Int_t decayIdxMax=1;
    if(enable2pmode==1) decayIdxMax=2;

    for(Int_t decayIdx=0; decayIdx<decayIdxMax; decayIdx++) {
        for(int cobo=0;cobo<no_cobos;cobo++){
            rwaveforms[decayIdx][cobo]->hasSignal.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->hasHit.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->hasFPN.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->doneFPN.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->isOverflow.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->hasOverflow = false;
            rwaveforms[decayIdx][cobo]->isDecay.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->hasImplant = false;
            rwaveforms[decayIdx][cobo]->hasDecay = false;
            rwaveforms[decayIdx][cobo]->coboIdx = 0;
            rwaveforms[decayIdx][cobo]->asadIdx = 0;
            rwaveforms[decayIdx][cobo]->waveform.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->corrwaveform.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->fpnwaveform.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->energy.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->time.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->baseline.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->PSDIntegral.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->PSDRatio.resize(maxasad*4);
            for(int i=0;i<maxasad;i++){
                for(int j=0;j<4;j++){
                    rwaveforms[decayIdx][cobo]->hasSignal[i*4+j].resize(68);
                    rwaveforms[decayIdx][cobo]->isOverflow[i*4+j].resize(68);
                    rwaveforms[decayIdx][cobo]->isDecay[i*4+j].resize(68);
                    rwaveforms[decayIdx][cobo]->waveform[i*4+j].resize(68);
                    rwaveforms[decayIdx][cobo]->corrwaveform[i*4+j].resize(68);
                    rwaveforms[decayIdx][cobo]->energy[i*4+j].resize(68);
                    rwaveforms[decayIdx][cobo]->time[i*4+j].resize(68);
                    rwaveforms[decayIdx][cobo]->baseline[i*4+j].resize(68);
                    rwaveforms[decayIdx][cobo]->PSDIntegral[i*4+j].resize(68);
                    rwaveforms[decayIdx][cobo]->PSDRatio[i*4+j].resize(68);
                    for(int k=0;k<68;k++){
                        rwaveforms[decayIdx][cobo]->waveform[i*4+j][k].resize(bucketmax);
                        rwaveforms[decayIdx][cobo]->corrwaveform[i*4+j][k].resize(bucketmax);
                    }
                    rwaveforms[decayIdx][cobo]->fpnwaveform[i*4+j].resize(bucketmax);
                }
            }
        }
    }
}

void LKFrameBuilder::RootRResetWaveforms() {
    Int_t decayIdxMax=1;
    if(enable2pmode==1) decayIdxMax=2;

    for(Int_t decayIdx=0; decayIdx<decayIdxMax; decayIdx++) {
        for(int cobo=0;cobo<no_cobos;cobo++){
            rwaveforms[decayIdx][cobo]->hasOverflow = false;
            rwaveforms[decayIdx][cobo]->hasImplant = false;
            rwaveforms[decayIdx][cobo]->hasDecay = false;
            rwaveforms[decayIdx][cobo]->isRejected = false;
            for(int i=0;i<maxasad;i++){
                for(int j=0;j<4;j++){
                    if(rwaveforms[decayIdx][cobo]->hasHit[i*4+j] || rwaveforms[decayIdx][cobo]->hasFPN[i*4+j]){
                        for(int k=0;k<68;k++){
                            if(!rwaveforms[decayIdx][cobo]->hasSignal[i*4+j][k]) continue;
                            rwaveforms[decayIdx][cobo]->hasSignal[i*4+j][k] = false;
                            rwaveforms[decayIdx][cobo]->isOverflow[i*4+j][k] = false;
                            rwaveforms[decayIdx][cobo]->isDecay[i*4+j][k] = 0;
                            rwaveforms[decayIdx][cobo]->energy[i*4+j][k] = 0.0;
                            rwaveforms[decayIdx][cobo]->time[i*4+j][k] = 0.0;
                            rwaveforms[decayIdx][cobo]->baseline[i*4+j][k] = 0.0;
                            rwaveforms[decayIdx][cobo]->PSDIntegral[i*4+j][k] = 0.0;
                            rwaveforms[decayIdx][cobo]->PSDRatio[i*4+j][k] = 0.0;
                            fill(rwaveforms[decayIdx][cobo]->waveform[i*4+j][k].begin(),
                                    rwaveforms[decayIdx][cobo]->waveform[i*4+j][k].end(),0);
                            fill(rwaveforms[decayIdx][cobo]->corrwaveform[i*4+j][k].begin(),
                                    rwaveforms[decayIdx][cobo]->corrwaveform[i*4+j][k].end(),0);
                        }
                    }
                    if(rwaveforms[decayIdx][cobo]->hasFPN[i*4+j]){
                        fill(rwaveforms[decayIdx][cobo]->fpnwaveform[i*4+j].begin(),
                                rwaveforms[decayIdx][cobo]->fpnwaveform[i*4+j].end(),0);
                    }
                    rwaveforms[decayIdx][cobo]->hasHit[i*4+j] = false;
                    rwaveforms[decayIdx][cobo]->hasFPN[i*4+j] = false;
                    rwaveforms[decayIdx][cobo]->doneFPN[i*4+j] = false;
                }
            }
        }
    }

    if(enabledraww==1){
        hWaveFormbyEvent[evtcounter%16]->Reset();
        hCorrWaveFormbyEvent[evtcounter%16]->Reset();
        //for(int i=0;i<=512;i++){
        //  for(int j=0;j<=512;j++){
        //    hWaveFormbyEvent[evtcounter%16]->SetBinContent(i,j,0);
        //    hCorrWaveFormbyEvent[evtcounter%16]->SetBinContent(i,j,0);
        //  }
        //}
    }
}

void LKFrameBuilder::RootRReset(){
    rGETHit = 0;
    rGETEventIdx = 0;
    rGETD2PTime = 0;
    rGETTimeStamp = 0;
    for(int i=0;i<=rGETMul;i++){
        rGETFrameNo[i] = 0;
        rGETDecayNo[i] = 0;
        rGETL1Aflag[i] = 0;
        rGETEnergy[i] = 0;
        rGETTime[i] = 0;
        rGETCobo[i] = 0;
        rGETAsad[i] = 0;
        rGETAget[i] = 0;
        rGETChan[i] = 0;
        for(int j=0;j<512;j++){
            rGETWaveformX[i][j] = 0;
            rGETWaveformY[i][j] = 0;
        }
    }
    rGETMul = 0;
}

void LKFrameBuilder::RootRCloseFile(){
    fInputFile->cd();
    fInputFile->Close();
}

void LKFrameBuilder::RootRWOpenFile(string & outputFileName, string & outputTreeName){
    if(IsFirstevent){
        fOutputFile = new TFile(outputFileName.c_str(),"recreate");
        fOutputTree = new TTree(outputTreeName.c_str(),"Experimental Event Data");
    }
    else{
        fOutputFile = new TFile(outputFileName.c_str(),"update");
        fOutputTree = (TTree *)fOutputFile->Get(outputTreeName.c_str());
    }
    if(enableupdatefast==0){
        fOutputTree->SetAutoSave(-1000);
    }
    cout << Form("Created new TTree named %s in %s (new=%d)",outputTreeName.c_str(),outputFileName.c_str(),IsFirstevent) << endl;
}

void LKFrameBuilder::RootRWInit(){
    wGETMul = 0;
    wGETHit = 0;
    wGETEventIdx = 0;
    wGETD2PTime = 0;
    wGETTimeStamp = 0;
    for(int i=0;i<4352;i++){
        wGETFrameNo[i] = 0;
        wGETDecayNo[i] = 0;
        wGETL1Aflag[i] = 0;
        wGETEnergy[i] = 0;
        wGETCobo[i] = 0;
        wGETAsad[i] = 0;
        wGETAget[i] = 0;
        wGETChan[i] = 0;
    }
    fOutputTree->Branch("mmMul",&wGETMul,"mmMul/I");
    fOutputTree->Branch("mmHit",&wGETHit,"mmHit/I");
    fOutputTree->Branch("mmEventIdx",&wGETEventIdx,"mmEventIdx/I");
    fOutputTree->Branch("mmD2PTime",&wGETD2PTime,"mmD2PTime/I");
    fOutputTree->Branch("mmTimeStamp",&wGETTimeStamp,"mmTimeStamp/I");
    fOutputTree->Branch("mmL1Aflag",wGETL1Aflag,"mmL1Aflag[mmMul]/I");
    fOutputTree->Branch("mmFrameNo",wGETFrameNo,"mmFrameNo[mmMul]/I");
    fOutputTree->Branch("mmDecayNo",wGETDecayNo,"mmDecayNo[mmMul]/I");
    fOutputTree->Branch("mmCobo",wGETCobo,"mmCobo[mmMul]/I");
    fOutputTree->Branch("mmAsad",wGETAsad,"mmAsad[mmMul]/I");
    fOutputTree->Branch("mmAget",wGETAget,"mmAget[mmMul]/I");
    fOutputTree->Branch("mmChan",wGETChan,"mmChan[mmMul]/I");
    fOutputTree->Branch("mmTime",wGETTime,"mmTime[mmMul]/F");
    fOutputTree->Branch("mmEnergy",wGETEnergy,"mmEnergy[mmMul]/F");
}

void LKFrameBuilder::RootRWriteEvent(){
    if(wGETMul>0){
        //cout << "Writing data: " << wGETEventIdx << " " << wGETMul << endl;
        fOutputFile->cd();
        fOutputTree->Fill();
        if(enableupdatefast==1){
            if(framecounter%16==0){
                fOutputTree->Write();
                fOutputFile->SaveSelf();
                if((fOutputFile->GetSize())>5000000){
                    enableupdatefast=0;
                }
                //cout << (int)(framecounter/16) << " times written, Size = " << fOutputFile->GetSize() << endl;
            }
        }else{
            if(framecounter%1600==0){
                fOutputTree->Write();
                fOutputFile->SaveSelf();
            }
        }
        framecounter++;
        wGETMul = 0;
        wGETHit = 0;
    }
}

void LKFrameBuilder::RootRWCloseFile(){
    fOutputFile->cd();
    fOutputFile->Close();
}

void LKFrameBuilder::RootRWReset(){
    for(int i=0;i<=wGETMul;i++){
        wGETFrameNo[i] = 0;
        wGETDecayNo[i] = 0;
        wGETL1Aflag[i] = 0;
        wGETTime[i] = 0;
        wGETEnergy[i] = 0;
        wGETCobo[i] = 0;
        wGETAsad[i] = 0;
        wGETAget[i] = 0;
        wGETChan[i] = 0;
    }
    wGETMul=0;
    wGETHit=0;
}

void LKFrameBuilder::RootRHInit(){
    fInputTree->SetBranchAddress("mmMul",&rGETMul);
    fInputTree->SetBranchAddress("mmHit",&rGETHit);
    fInputTree->SetBranchAddress("mmEventIdx",&rGETEventIdx);
    fInputTree->SetBranchAddress("mmD2PTime",&rGETD2PTime);
    fInputTree->SetBranchAddress("mmTimeStamp",&rGETTimeStamp);
    fInputTree->SetBranchAddress("mmFrameNo",rGETFrameNo);
    fInputTree->SetBranchAddress("mmDecayNo",rGETDecayNo);
    fInputTree->SetBranchAddress("mmCobo",rGETCobo);
    fInputTree->SetBranchAddress("mmAsad",rGETAsad);
    fInputTree->SetBranchAddress("mmAget",rGETAget);
    fInputTree->SetBranchAddress("mmChan",rGETChan);
    fInputTree->SetBranchAddress("mmTime",rGETTime);
    fInputTree->SetBranchAddress("mmEnergy",rGETEnergy);

    rGETMul = 0;
    rGETHit = 0;
    rGETEventIdx = 0;
    rGETD2PTime = 0;
    rGETTimeStamp = 0;
    for(int i=0;i<4352;i++){
        rGETFrameNo[i] = 0;
        rGETDecayNo[i] = 0;
        rGETL1Aflag[i] = 0;
        rGETEnergy[i] = 0;
        rGETTime[i] = 0;
        rGETCobo[i] = 0;
        rGETAsad[i] = 0;
        rGETAget[i] = 0;
        rGETChan[i] = 0;
    }
    fNumberEvents = fInputTree->GetEntries();
    RootRHInitWaveforms();
}

void LKFrameBuilder::RootRHInitWaveforms() {
    Int_t decayIdxMax=1;
    if(enable2pmode==1) decayIdxMax=2;

    for(Int_t decayIdx=0; decayIdx<decayIdxMax; decayIdx++) {
        for(int cobo=0;cobo<no_cobos;cobo++){
            rwaveforms[decayIdx][cobo]->hasSignal.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->hasHit.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->hasFPN.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->doneFPN.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->isOverflow.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->hasOverflow = false;
            rwaveforms[decayIdx][cobo]->isDecay.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->hasImplant = false;
            rwaveforms[decayIdx][cobo]->hasDecay = false;
            rwaveforms[decayIdx][cobo]->coboIdx = 0;
            rwaveforms[decayIdx][cobo]->asadIdx = 0;
            rwaveforms[decayIdx][cobo]->energy.resize(maxasad*4);
            rwaveforms[decayIdx][cobo]->time.resize(maxasad*4);
            for(int i=0;i<maxasad;i++){
                for(int j=0;j<4;j++){
                    rwaveforms[decayIdx][cobo]->hasSignal[i*4+j].resize(68);
                    rwaveforms[decayIdx][cobo]->isOverflow[i*4+j].resize(68);
                    rwaveforms[decayIdx][cobo]->isDecay[i*4+j].resize(68);
                    rwaveforms[decayIdx][cobo]->energy[i*4+j].resize(68);
                    rwaveforms[decayIdx][cobo]->time[i*4+j].resize(68);
                }
            }
        }
    }
}

void LKFrameBuilder::RootRHResetWaveforms()
{
    Int_t decayIdxMax=1;
    if(enable2pmode==1) decayIdxMax=2;

    for(Int_t decayIdx=0; decayIdx<decayIdxMax; decayIdx++) {
        for(int cobo=0;cobo<no_cobos;cobo++){
            rwaveforms[decayIdx][cobo]->hasOverflow = false;
            rwaveforms[decayIdx][cobo]->hasImplant = false;
            rwaveforms[decayIdx][cobo]->hasDecay = false;
            for(int i=0;i<maxasad;i++){
                for(int j=0;j<4;j++){
                    if(rwaveforms[decayIdx][cobo]->hasHit[i*4+j] || rwaveforms[decayIdx][cobo]->hasFPN[i*4+j]){
                        for(int k=0;k<68;k++){
                            if(!rwaveforms[decayIdx][cobo]->hasSignal[i*4+j][k]) continue;
                            rwaveforms[decayIdx][cobo]->hasSignal[i*4+j][k] = false;
                            rwaveforms[decayIdx][cobo]->isOverflow[i*4+j][k] = false;
                            rwaveforms[decayIdx][cobo]->isDecay[i*4+j][k] = 0;
                            rwaveforms[decayIdx][cobo]->energy[i*4+j][k] = 0.0;
                            rwaveforms[decayIdx][cobo]->time[i*4+j][k] = 0.0;
                        }
                    }
                    rwaveforms[decayIdx][cobo]->hasHit[i*4+j] = false;
                    rwaveforms[decayIdx][cobo]->hasFPN[i*4+j] = false;
                    rwaveforms[decayIdx][cobo]->doneFPN[i*4+j] = false;
                }
            }
        }
    }
}

void LKFrameBuilder::RootRHReset()
{
    rGETHit = 0;
    rGETEventIdx = 0;
    rGETD2PTime = 0;
    rGETTimeStamp = 0;
    for(int i=0;i<=rGETMul;i++){
        rGETFrameNo[i] = 0;
        rGETDecayNo[i] = 0;
        rGETL1Aflag[i] = 0;
        rGETEnergy[i] = 0;
        rGETTime[i] = 0;
        rGETCobo[i] = 0;
        rGETAsad[i] = 0;
        rGETAget[i] = 0;
        rGETChan[i] = 0;
    }
    rGETMul = 0;
}

void LKFrameBuilder::SetBucketSize(int BucketSize){
    bucketmax = BucketSize;
}

void LKFrameBuilder::SetRootConverter(int flag){
    enableroot = flag; // this is not used anywere...
}

int LKFrameBuilder::GetForceReadTree(){
    return forcereadtree;
}

void LKFrameBuilder::SetIgnoreMM(int flag){
    ignoremm = flag;
}

void LKFrameBuilder::SetDrawWaveform(int flag){
    enabledraww = flag;
}

void LKFrameBuilder::SetCleanTrack(int flag){
    enablecleantrack = flag;
}

void LKFrameBuilder::SetDrawTrack(int flag){
    enabletrack = flag;
}

void LKFrameBuilder::SetSkipEvents(int flag){
    enableskipevent = flag;
}

void LKFrameBuilder::SetfirstEventNo(int flag){
    firsteventno = flag;
}

void LKFrameBuilder::ReadGoodEventList(string filename){
    cout << "goodevtlistfilename is " << filename << endl;
    int goodevtno=0;
    ifstream goodevtlist;
    goodevtlist.open(filename.data());
    if(goodevtlist.fail()==true){
        cerr<<"The goodevtlist file wasn't opened!"<<endl;
    }
    while(goodevtlist.good()){
        goodevtlist >> goodevtno;
        if(goodevtno<maxevtno){
            evtmask[goodevtno]=true;
        }else{
            cout << "goodevtno >= 10^8!" << endl;
        }
    }
    goodevtlist.close();
}

void LKFrameBuilder::ReadMapChanToMM(string filename){
    ifstream MapEData;
    UInt_t asadid;
    UInt_t agetid;
    UInt_t dchanid;
    UInt_t spxidx;
    UInt_t spxidy;

    MapEData.open(filename.data());
    if(MapEData.fail()==true){
        cerr<<"The ChanToMM_Map file wasn't opened!"<<endl;
    }
    while(MapEData.good()){
        MapEData >> asadid >> agetid >> dchanid >> spxidx >> spxidy;
        mapchantomm->pxidx[asadid][agetid][dchanid]=spxidx;
        mapchantomm->pxidy[asadid][agetid][dchanid]=spxidy;
    }
    MapEData.close();
}

void LKFrameBuilder::ReadMapChanToSi(string filename){
    ifstream MapEData;
    UInt_t asadid;
    UInt_t agetid;
    UInt_t chanid;
    UInt_t spxidx;
    UInt_t spxidy;

    MapEData.open(filename.data());
    if(MapEData.fail()==true){
        cerr<<"The ChanToSi_Map file wasn't opened!"<<endl;
    }
    while(MapEData.good()){
        MapEData >> asadid >> agetid >> chanid >> spxidx >> spxidy;
        mapchantosi->pxidx[asadid][agetid][chanid]=spxidx;
        mapchantosi->pxidy[asadid][agetid][chanid]=spxidy;
    }
    MapEData.close();
}

void LKFrameBuilder::ReadMapChanToX6(){
    UInt_t tX6asad;
    UInt_t tX6aget;
    UInt_t tX6chan;
    UInt_t tX6flag;
    UInt_t tX6det;
    UInt_t tX6strip;
    ifstream fmapchantoX6;
    fmapchantoX6.open("mapchantoX6_CRIB.txt");
    if(fmapchantoX6.fail()==true){
        cerr<<"The mapchantoX6 file wasn't opened!"<<endl;
    }else{
        //cout << "mapchantoX6 open" <<endl;
    }
    while(fmapchantoX6.good()){
        fmapchantoX6 >> tX6asad >> tX6aget >> tX6chan >> tX6flag >> tX6det >> tX6strip;
        mapchantox6->X6flag[tX6asad][tX6aget][tX6chan]=tX6flag;
        mapchantox6->X6det[tX6asad][tX6aget][tX6chan]=tX6det;
        mapchantox6->X6strip[tX6asad][tX6aget][tX6chan]=tX6strip;
    }
    fmapchantoX6.close();

    Double_t tX6x, tX6y, tX6z;
    ifstream mapdimtoX6;
    mapdimtoX6.open("mapdimtoX6_CRIB.txt");
    if(mapdimtoX6.fail()==true){
        cerr<<"The mapdimtoX6 file wasn't opened!"<<endl;
    }else{
        //cout << "mapdimtoX6 open" <<endl;
    }
    while(mapdimtoX6.good()){
        mapdimtoX6 >> tX6det >> tX6strip >> tX6x >> tX6y >> tX6z;
        mapchantox6 -> X6posx[tX6det][tX6strip] = tX6x;
        mapchantox6 -> X6posy[tX6det][tX6strip] = tX6y;
        mapchantox6 -> X6posz[tX6det][tX6strip] = tX6z;
    }
    mapdimtoX6.close();

    Int_t tX6CsIasad, tX6CsIaget, tX6CsIchan, tX6CsIRCflag, tX6CsIloc, tX6CsIref;
    ifstream mapchantoX6CsI;
    mapchantoX6CsI.open("mapchantoX6CsI_CRIB.txt");
    if(mapchantoX6CsI.fail()==true){
        cerr<<"The mapchantoX6CsI file wasn't opened!"<<endl;
    }else{
        //cout << "mapchantoX6CsI open" <<endl;
    }
    while(mapchantoX6CsI.good()){
        mapchantoX6CsI >> tX6CsIasad >> tX6CsIaget >> tX6CsIchan >> tX6CsIRCflag >> tX6CsIloc >> tX6CsIref;
        mapchantox6 -> CsICTnum[tX6CsIasad][tX6CsIaget][tX6CsIchan] = tX6CsIRCflag;
        mapchantox6 -> CsI_X6det[tX6CsIasad][tX6CsIaget][tX6CsIchan] = tX6CsIloc;
        mapchantox6 -> CsI_X6ud[tX6CsIasad][tX6CsIaget][tX6CsIchan] = tX6CsIref;
    }
    mapchantoX6CsI.close();
}

void LKFrameBuilder::SetEnergyMethod(int flag){
    energymethod = flag;
}

void LKFrameBuilder::SetReadRF(int flag){
    readrw = flag;
}

void LKFrameBuilder::ReadResponseWaveform(string filename){
    ifstream ResponseData;
    UInt_t type;
    UInt_t timebucket;
    UInt_t amplitude;

    ResponseData.open(filename.data());
    if(ResponseData.fail()==true){
        cerr<<"The ResponseData file wasn't opened!"<<endl;
    }else{
        cout << "The ResponseData file: " << filename.data() <<endl;
    }
    while(ResponseData.good()){
        ResponseData >> type >> timebucket >> amplitude;
        response[type][timebucket] = amplitude;
        deconvrep[type] = 15;
        deconviter[type] = 15;
        deconvboost[type] = 1;
    }

    Int_t min_val=4095;
    for(int i=0;i<timebucket;i++) {
        if(response[0][i]<min_val && response[0][i]>0) min_val=response[0][i]; //find the min val
    }
    for(int i=0;i<timebucket;i++) {
        response[0][i]-=min_val;
        if(response[0][i]<0) response[0][i]=0;
    }
    cout<<"Min value: "<<min_val<<endl;
    ResponseData.close();
}

void LKFrameBuilder::SetResponseSample(Int_t type, Int_t evtno, Int_t buckcut, Int_t buckwidth, Int_t cobo, Int_t asad, Int_t aget, Int_t chan, Int_t rep, Int_t iter, Int_t boost)
{
    responsesample[type][0] = evtno;
    responsesample[type][1] = cobo;
    responsesample[type][2] = asad;
    responsesample[type][3] = aget;
    responsesample[type][4] = chan;
    responsesample[type][5] = buckcut;
    responsesample[type][6] = buckwidth;
    deconvrep[type] = rep;
    deconviter[type] = iter;
    deconvboost[type] = boost;
    if(maxrespsample<=type) maxrespsample = type;
}

Bool_t LKFrameBuilder::IsResponseSample(Int_t type, Int_t cobo, Int_t asad, Int_t aget, Int_t chan){
    if(responsesample[type][1] == cobo && responsesample[type][2] == asad && responsesample[type][3] == aget && responsesample[type][4] == chan) return true;
    return false;
}

void LKFrameBuilder::SetResponseWaveform(){
    Int_t decayIdx = 0;
    Int_t coboIdx = 0;
    Int_t asadIdx = 0;
    Int_t agetIdx = 0;
    Int_t chanIdx = 0;
    Int_t buckbegin;
    Int_t buckend;
    Int_t maxreventIdx=0;
    Int_t baseline=0;
    for(int j=0;j<=maxrespsample;j++){
        if(maxreventIdx <= responsesample[j][0]){
            maxreventIdx = responsesample[j][0];
        }
    }
    for(int i=0;i<fNumberEvents;i++){
        fInputFile->cd();
        RootRInit();
        fInputTree->GetEntry(i);
        reventIdx = rGETEventIdx;
        if(reventIdx > maxreventIdx) return;
        for(int j=0;j<=maxrespsample;j++){
            if(responsesample[j][0] == reventIdx){
                for(Int_t i=0; i<rGETMul; i++){
                    decayIdx = rGETDecayNo[i];
                    coboIdx = rGETCobo[i];
                    asadIdx = rGETAsad[i];
                    agetIdx = rGETAget[i];
                    chanIdx = rGETChan[i];
                    if((ignoremm==0) || (ignoremm==1 && coboIdx>0)) continue;// skip MM waveform data
                    if((chanIdx==11||chanIdx==22||chanIdx==45||chanIdx==56)){
                        rwaveforms[decayIdx][coboIdx]->hasFPN[asadIdx*4+agetIdx] = true;
                    }
                    for(int j=0;j<bucketmax;j++){
                        rwaveforms[decayIdx][coboIdx]->waveform[asadIdx*4+agetIdx][chanIdx][j] = rGETWaveformY[i][j];
                    }
                }

                for(Int_t i=0; i<rGETMul; i++){
                    decayIdx = rGETDecayNo[i];
                    coboIdx = rGETCobo[i];
                    asadIdx = rGETAsad[i];
                    agetIdx = rGETAget[i];
                    chanIdx = rGETChan[i];
                    if((ignoremm==0) || (ignoremm==1 && coboIdx>0)) continue;// skip MM waveform data
                    if(chanIdx!=11 && chanIdx!=22 && chanIdx!=45 && chanIdx!=56){ // We want to skip the FPN channels.
                        //cout << j << " " << reventIdx << " " << coboIdx << " " << asadIdx << " " << agetIdx << " " << chanIdx << endl;
                        if(IsResponseSample(j,coboIdx,asadIdx,agetIdx,chanIdx)){
                            buckbegin = responsesample[j][5];
                            buckend = responsesample[j][5] + responsesample[j][6];
                            GetAverageFPN(decayIdx,coboIdx,asadIdx,agetIdx);
                            GetCorrWaveform(decayIdx,coboIdx,asadIdx,agetIdx,chanIdx);
                            baseline = GetBaseline(decayIdx,coboIdx,asadIdx,agetIdx,chanIdx);
                            cout << "rwtype eventNo CoBo AsAd Aget Chan Bucket Amplitude Baseline" << endl;
                            for(int buck=buckbegin;buck<buckend;buck++){
                                response[j][buck-buckbegin] = rwaveforms[decayIdx][coboIdx]->corrwaveform[asadIdx*4+agetIdx][chanIdx][buck]-baseline+decoffset;
                                cout << j << " " << reventIdx << " " << coboIdx << " " << asadIdx << " " << agetIdx << " " << chanIdx << " " << buck-buckbegin << " " << response[j][buck-buckbegin] << " " << baseline << endl;
                            }
                        }
                    }
                }
            }
        }
    }
}

void LKFrameBuilder::GetMaxResponseWaveform(){
    Int_t maxamplitude;
    for(int i=0;i<=maxrespsample;i++){
        maxamplitude=-10000;
        for(int buck=0;buck<512;buck++){
            if(maxamplitude<response[i][buck]){
                maxamplitude = response[i][buck];
                maxresponse[i] = buck;
            }
        }
    }
}

void LKFrameBuilder::GetSigmaResponseWaveform(){
    for(int i=0;i<=maxrespsample;i++){
        hResponse[i] = new TH1D(Form("hResponse_%d",i),Form("hResponse_%d",i),512,0,512);
        if(i==0){
            fResponse[i] = new TF1(Form("fResponse_%d",i), this, &LKFrameBuilder::ShaperF_GET1, 0, 512, 6);
            fResponse[i]->SetParNames("offset", "amplitude", "peakAt", "sigma", "power", "p2");
            fResponse[i]->SetParameters(10, 500, 80, 15, 3, 0.2);
            fResponse[i]->SetParLimits(1, 0, 10000);
            fResponse[i]->SetParLimits(3, 0, 50);
        }else if(i>0) {
            fResponse[i] = new TF1(Form("fResponse_%d",i), this, &LKFrameBuilder::ShaperF_MSCF, 0, 512, 6);
            fResponse[i]->SetParNames("offset", "amplitude", "peakAt", "sigma", "power", "p2");
            fResponse[i]->SetParameters(10, 2000, 80, 25, 3, 0.2);
            fResponse[i]->SetParLimits(1, 0, 10000);
            fResponse[i]->SetParLimits(3, 0, 50);
        }
        /*
           fResponse[i] = new TF1(Form("fResponse_%d",i),conv, 0,512);
           for(int ii=0; ii<p2[0]; ii++) {
           f_conv->FixParameter(ii, p2[ii]);
           }
         */

        for(int buck=0;buck<512;buck++){
            hResponse[i]->Fill(buck,response[i][buck]);
        }
        hResponse[i]->Fit(Form("fResponse_%d",i),"QN"); //Q=Quiet and N=Do not draw
        responsesigma[i] = fResponse[i]->GetParameter(3);
        cout << "Sigma[" << i << "] = " << responsesigma[i] << endl;
    }
}

Double_t LKFrameBuilder::ShaperF_GET1(Double_t *x, Double_t *p) {
    Double_t semiGaus = p[1]*22.68113723*TMath::Exp(-3.0*(x[0] - p[2])/p[3])*sin((x[0] - p[2])/p[3])*pow((x[0] - p[2])/p[3], p[4]);
    return (x[0] >= p[2]) ? p[0] + x[0]*p[5] + semiGaus : p[0] + x[0]*p[5];
    /*
       Double_t val=0.;
       for(int i=0; i<(p[0]-2)/2; i++) {
       if(!(x[0]<p[2*i+2+1] || x[0]>1000000.))
       val += p[2*i+2] * 22.68113723 * exp(-3.*(x[0]-p[2*i+2+1])/p[1]) * sin((x[0]-p[2*i+2+1])/p[1]) * pow((x[0]-p[2*i+2+1])/p[1], 3);
       }
       return(val);
     */
}

Double_t LKFrameBuilder::ShaperF_MSCF(Double_t *x, Double_t *p) {
    Double_t semiGaus = p[1]*21.928*TMath::Exp(-3.0*(x[0] - p[2])/p[3])*sin((x[0] - p[2])/p[3])*pow((x[0] - p[2])/p[3], p[4]);
    return (x[0] >= p[2]) ? p[0] + x[0]*p[5] + semiGaus : p[0] + x[0]*p[5];
}

Double_t LKFrameBuilder::ShaperF_GET2(Double_t *x, Double_t *p) {
    Double_t semiGaus = p[0] + 0.5*p[1]*(1.0 + TMath::TanH( (x[0]-p[2])/p[3])) * TMath::Exp(-p[4]*((x[0]-p[2])*(x[0]-p[2])));
    return semiGaus;
}
