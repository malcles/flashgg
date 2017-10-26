// -*- C++ -*-
//
// Package:    Commissioning/outputAnalyzerDiMuHighMass
// Class:      outputAnalyzerDiMuHighMass
//
/**\class outputAnalyzerDiMuHighMass outputAnalyzerDiMuHighMass.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Julie Malcles
//         Created:  Tue, 02 Dec 2014 10:57:22 GMT
//
//


#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//for miniAOD:
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "flashgg/MicroAOD/interface/VertexSelectorBase.h"
#include "flashgg/DataFormats/interface/VertexCandidateMap.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "TTree.h"

// **********************************************************************

// define the structures used to create tree branches and fill the trees


// **********************************************************************

using namespace flashgg;

using namespace std;
using namespace edm;
using namespace reco;
using namespace math;

struct BSInfo {
    float x;
    float y;
    float z;
    float sigmaz;

};

struct MyInfo {
    int nVtx;
    int nVtxNoMu;
    int nPU;
    float mcWeight;
    float zTrue;
    float zRecoTrue;
    float zRecoTrueNoMu;
    float zNextRecoTrue;
    float zNextRecoTrueNoMu;
    float zRecoFromMu;
    float zRecoFromMuNoMu;        
    int nTracks;
    int nTracksNoMu;        
    int muonsConsistentNoMu;
    int muonsConsistentWithMu;
    float diMuonMass;
    float diMuonPt;
    float zChosenWithMu;
    float zChosenNoMu;
    float dzMatchingWithMu;
    float dzMatchingNoMu;
    float mvaProbWithMu;
    float mvaProbNoMu;
    float probWithMu;
    float probNoMu;
    float logSumPt2WithMu;
    float logSumPt2NoMu;
    float ptBalWithMu;
    float ptBalNoMu;
    float ptAsymWithMu;
    float ptAsymNoMu;
    float dZ1WithMu;
    float dZ1NoMu;
    float dZ2WithMu;
    float dZ2NoMu;

    int nTracksRecoTrue;
    int nTracksRecoTrueNoMu;
    int nTracksChosenWithMu;
    int nTracksChosenNoMu;
    int nTracksZerothWithMu;
    int nTracksZerothNoMu;
    float ptRecoTrue;
    float ptRecoTrueNoMu;
    float ptChosenWithMu;
    float ptChosenNoMu;
    float ptZerothWithMu;
    float ptZerothNoMu;
    float zZerothWithMu;
    float zZerothNoMu;

};
struct MuonInfo {
    float eta;
    float pt;
    float vx;
    float vy;
    float vz;
    float dB;
    float edB;
    float sip;
    float dzBS;
    float dxyBS;
};

struct SignalInfo {

    int nvertex;
    int ndipho;
    int dipho_index;

    float LogSumPt2;
    float PtBal;
    float PtAsym;
    float NConv;
    float PullConv;
    float DZtrue;

};

struct BackgroundInfo {

    int nvertex;
    int ndipho;
    int dipho_index;

    float LogSumPt2;
    float PtBal;
    float PtAsym;
    float NConv;
    float PullConv;
    float DZtrue;

};

// **********************************************************************

class outputAnalyzerDiMuHighMass : public edm::EDAnalyzer
{
public:
    explicit outputAnalyzerDiMuHighMass( const edm::ParameterSet & );
    ~outputAnalyzerDiMuHighMass();
    
    static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );
    
    
private:
    
    edm::Service<TFileService> fs_;
    
    virtual void beginJob() override;
    virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
    virtual void endJob() override;

    void initEventStructure();
   

    void printMuon( MuonInfo inf ,string name);
    void getMCTruthVertexCoord(  const std::vector<edm::Ptr<reco::GenParticle>>& gens , float& vx, float&vy, float&vz );
    int getRecoClosestToTrueVertexIndex(  const std::vector<edm::Ptr<reco::GenParticle>>& gens, const vector<edm::Ptr<reco::Vertex> > & vertices, double dzMatch = 50. );
    int getRecoNextToClosestToTrueVertexIndex(  const std::vector<edm::Ptr<reco::GenParticle>>& gens , const vector<edm::Ptr<reco::Vertex> > & vertices,  Ptr<reco::Vertex>  closest , double dzMatch =50. );
    int getRecoWithMuonsVertexIndex( const vector<edm::Ptr<reco::Vertex> > & vertices, Ptr<pat::Muon> mu1, Ptr<pat::Muon> mu2, int& muonsConsistent, double &dzFound , double dzMatch = 0.2 );
    int getMatchingVertexIndexNoMu( const vector<edm::Ptr<reco::Vertex> > & verticesNoMu, Ptr<reco::Vertex> vertexTrueWithMu, double dzMatch = 0.5);
    int getClosestVertexIndex( const vector<edm::Ptr<reco::Vertex> > & vertices, Ptr<reco::Vertex> vertexTrue, double dzMatch = 0.5);


    edm::EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
    edm::EDGetToken genEventInfoToken_;
    edm::EDGetTokenT<View<pat::Muon> > patMuonToken_;
    edm::EDGetTokenT<edm::View<reco::Vertex> >               vertexToken_;
    edm::EDGetTokenT<edm::View<reco::Vertex> >               vertexTokenNoMu_;
    edm::EDGetTokenT<edm::View<reco::Vertex> >               vertexTokenWithMu_;
    EDGetTokenT< VertexCandidateMap > vertexCandidateMapToken_;
    EDGetTokenT< VertexCandidateMap > vertexCandidateMapNoMuToken_;
    unique_ptr<VertexSelectorBase> vertexSelector_;
    //unique_ptr<VertexSelectorBase> vertexSelectorZero_;
    edm::EDGetTokenT<reco::BeamSpot > beamSpotToken_;
  edm::EDGetTokenT<edm::View<PileupSummaryInfo> >  PileUpToken_;

    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    std::map<std::string, unsigned int> triggerIndices_;

    bool isData = false;
    bool relaxDZ=true;
    bool verbose = true;
    bool verbose2 = false;
    bool muonEventMC = false;
    int  nMuonInAccMC=0;
    int  nZInAccMC =0;
    int noMatch=0;
    int nPass=0;
    
    int nPatMuonReco=0;
    int nZReco=0;

    bool passTrigger = false;
    std::vector<std::string> pathName_;
    //bool doReweightBS_;
    // double BSData_;
    //double BSMC_;

    MyInfo myInfo;
    MuonInfo mu1Info;
    MuonInfo mu2Info;
    BSInfo bsInfo;
    TTree *vtxTree;

    TTree *signalTree;
    TTree *backgroundTree;

    BackgroundInfo bkgInfo;
    SignalInfo sigInfo;
};

// ******************************************************************************************


//
// constructors and destructor
//
outputAnalyzerDiMuHighMass::outputAnalyzerDiMuHighMass( const edm::ParameterSet &iConfig ):
    genParticleToken_(consumes<View<reco::GenParticle> >(iConfig.getUntrackedParameter<InputTag> ("GenParticleTag", InputTag("prunedGenParticles")))),
    genEventInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter <edm::InputTag> ("genEventInfoProduct"))),
    patMuonToken_( consumes<View<pat::Muon> >( iConfig.getParameter<InputTag>( "patMuonTag" ) ) ),
    vertexToken_( consumes<View<reco::Vertex> >( iConfig.getUntrackedParameter<InputTag> ( "VertexTag", InputTag( "offlineSlimmedPrimaryVertices" ) ) ) ),
    vertexTokenNoMu_( consumes<View<reco::Vertex> >( iConfig.getUntrackedParameter<InputTag> ( "VertexTagNoMu", InputTag( "offlinePrimaryVerticesNoMu" ) ) ) ),
    vertexTokenWithMu_( consumes<View<reco::Vertex> >( iConfig.getUntrackedParameter<InputTag> ( "VertexTagWithMu", InputTag( "offlinePrimaryVerticesWithMu" ) ) ) ),
    vertexCandidateMapToken_( consumes<VertexCandidateMap>( iConfig.getParameter<InputTag>( "VertexCandidateMapTag" ) ) ),
    vertexCandidateMapNoMuToken_( consumes<VertexCandidateMap>( iConfig.getParameter<InputTag>( "VertexCandidateMapTagNoMu" ) ) ),
    beamSpotToken_( consumes<reco::BeamSpot >( iConfig.getUntrackedParameter<InputTag>( "BeamSpotTag", InputTag( "offlineBeamSpot" ) ) ) ),
    PileUpToken_( consumes<View<PileupSummaryInfo> >( iConfig.getUntrackedParameter<InputTag> ( "PileUpTag", InputTag( "addPileupInfo" ) ) ) ),
    triggerBits_( consumes<edm::TriggerResults>( iConfig.getParameter<edm::InputTag>( "bits" ) ) )
{
        const std::string &VertexSelectorName = iConfig.getParameter<std::string>( "VertexSelectorName" );
        vertexSelector_.reset( FlashggVertexSelectorFactory::get()->create( VertexSelectorName, iConfig ) );
        //        const std::string &VertexSelectorNameZero = iConfig.getParameter<std::string>( "VertexSelectorNameZero" );
        //        vertexSelectorZero_.reset( FlashggVertexSelectorFactory::get()->create( VertexSelectorNameZero, iConfig ) );
        pathName_ = iConfig.getParameter<std::vector<std::string>>( "pathName" ) ; 
        //        doReweightBS_        = iConfig.getParameter<bool>( "doReweightBS" );
        //BSData_           = iConfig.getParameter<double>( "BSData" ); 
        //BSMC_           = iConfig.getParameter<double>( "BSMC" );  
        

}


outputAnalyzerDiMuHighMass::~outputAnalyzerDiMuHighMass()
{


}



void
outputAnalyzerDiMuHighMass::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
{

    passTrigger=false;

    Handle<View<pat::Muon> >  patMuons;
    iEvent.getByToken( patMuonToken_, patMuons );

    Handle<View<reco::Vertex> > primaryVertices;
    iEvent.getByToken( vertexToken_, primaryVertices );

    Handle<View<reco::Vertex> > primaryVerticesNoMu;
    iEvent.getByToken( vertexTokenNoMu_, primaryVerticesNoMu );

    Handle<View<reco::Vertex> > primaryVerticesWithMu;
    iEvent.getByToken( vertexTokenWithMu_, primaryVerticesWithMu );

    Handle<VertexCandidateMap> vertexCandidateMap;
    iEvent.getByToken( vertexCandidateMapToken_, vertexCandidateMap );

    Handle<VertexCandidateMap> vertexCandidateMapNoMu;
    iEvent.getByToken( vertexCandidateMapNoMuToken_, vertexCandidateMapNoMu );

    Handle<reco::BeamSpot> recoBeamSpotHandle;
    iEvent.getByToken( beamSpotToken_, recoBeamSpotHandle );
    math::XYZPoint beamSpot;
    double sigmaZ=0;
    if( recoBeamSpotHandle.isValid() ) {
        beamSpot = recoBeamSpotHandle->position();
        sigmaZ = recoBeamSpotHandle->sigmaZ();

    } else {
        cout << " WARNING NO VALID BEAM SPOT: this should not happen!" << endl;
    }

    edm::Handle<edm::TriggerResults> triggerBits;
    iEvent.getByToken( triggerBits_, triggerBits );

    initEventStructure();

    isData = iEvent.isRealData();

    // Pass HLT

    if( !triggerBits.isValid() ) {
        LogDebug( "" ) << "TriggerResults product not found - returning result=false!";
        return;
    }

    const edm::TriggerNames &triggerNames = iEvent.triggerNames( *triggerBits );

    triggerIndices_.clear();
    for( unsigned int i = 0; i < triggerNames.triggerNames().size(); i++ ) {
        std::string trimmedName = HLTConfigProvider::removeVersion( triggerNames.triggerName( i ) );
        triggerIndices_[trimmedName] = triggerNames.triggerIndex( triggerNames.triggerName( i ) );
    }

    passTrigger = false;

    for (std::map<std::string, unsigned int>::const_iterator cit = triggerIndices_.begin(); cit != triggerIndices_.end(); cit++) {
        if (triggerBits->accept(cit->second)) {

            cout << " accept " << " " << cit->first << " second : " << cit->second << endl;
            std::vector<std::string>::const_iterator it = find(pathName_.begin(), pathName_.end(), cit->first);
            
            if (it != pathName_.end()){
                passTrigger = true;
                cout << "passTrigger = true" << endl;
            }
        }
    }

    if(!passTrigger) return;

    
    int trueRecoVtxIdx=-1;
    int trueRecoVtxIdxNoMu=-1;
    int nextTrueRecoVtxIdx=-1;
    int nextTrueRecoVtxIdxNoMu=-1;
    float zTrueVtxIdx=-999.;
    float yTrueVtxIdx=-999.;
    float xTrueVtxIdx=-999.;

    Ptr<reco::Vertex> vertexTrueReco;
    Ptr<reco::Vertex> vertexTrueRecoNoMu;

    Ptr<reco::Vertex> nextVertexTrueReco;
    Ptr<reco::Vertex> nextVertexTrueRecoNoMu;

    if( ! isData ) {

        double maxDZ=50;

        Handle<View<reco::GenParticle> > genParticles;
        iEvent.getByToken(genParticleToken_,genParticles);
        const std::vector<edm::Ptr<reco::GenParticle>>& gens = genParticles->ptrs();
      
        Handle<GenEventInfoProduct> genEventInfo;
        iEvent.getByToken(genEventInfoToken_,genEventInfo);
        myInfo.mcWeight = genEventInfo->weight();

        const vector<edm::Ptr<reco::Vertex> > & vtxsStdProd = primaryVertices->ptrs();
        const vector<edm::Ptr<reco::Vertex> > & vtxs = primaryVerticesWithMu->ptrs();
        const vector<edm::Ptr<reco::Vertex> > & vtxsnomu = primaryVerticesNoMu->ptrs();

        if(vtxsStdProd.size()!=vtxs.size()) cout<< "ERROR DIFFERENT SIZES "<< vtxsStdProd.size()<<" "<<vtxs.size()<< endl;

        trueRecoVtxIdx= getRecoClosestToTrueVertexIndex( gens, vtxs , maxDZ );
        trueRecoVtxIdxNoMu= getRecoClosestToTrueVertexIndex( gens, vtxsnomu , maxDZ );

        if(trueRecoVtxIdx!=-1){
            vertexTrueReco=vtxs[trueRecoVtxIdx];        
            nextTrueRecoVtxIdx= getRecoNextToClosestToTrueVertexIndex( gens, vtxs ,vertexTrueReco, maxDZ );
        } else cout<<" NO TRUE RECO"<<endl;
        if(trueRecoVtxIdxNoMu!=-1){
            vertexTrueRecoNoMu=vtxsnomu[trueRecoVtxIdxNoMu];
            nextTrueRecoVtxIdxNoMu= getRecoNextToClosestToTrueVertexIndex( gens, vtxsnomu , vertexTrueRecoNoMu, maxDZ );
        } else cout<<" NO TRUE RECO NoMu "<<vtxsnomu.size()<<" "<<vtxs.size()<<endl;
        
        if(nextTrueRecoVtxIdx!=-1) nextVertexTrueReco=vtxs[nextTrueRecoVtxIdx]; 
        if(nextTrueRecoVtxIdxNoMu!=-1) nextVertexTrueRecoNoMu=vtxsnomu[nextTrueRecoVtxIdxNoMu]; 

        getMCTruthVertexCoord(  gens , xTrueVtxIdx, yTrueVtxIdx,zTrueVtxIdx );

        cout<<" ZTRUE="<<zTrueVtxIdx<< endl;
        if(trueRecoVtxIdx!=-1) cout<<" ZTRUERECO="<<vertexTrueReco->z()<< endl;

        Handle<View< PileupSummaryInfo> > PileupInfos;
        iEvent.getByToken( PileUpToken_, PileupInfos ); 
        
        float  pun = 0;
        
        // pileup info
        for( unsigned int PVI = 0; PVI < PileupInfos->size(); ++PVI ) {
            Int_t pu_bunchcrossing = PileupInfos->ptrAt( PVI )->getBunchCrossing();
            if( pu_bunchcrossing == 0 ) {
                pun = PileupInfos->ptrAt( PVI )->getPU_NumInteractions();
            }
        }
        
        myInfo.nPU=pun;
    }
    
    
    nPatMuonReco=0;
  
    for(unsigned int j = 0 ; j < patMuons->size(); j++){

        Ptr<pat::Muon> pat_muon = patMuons->ptrAt(j);
        
        if(!pat_muon->isLooseMuon() || pat_muon->pt()<10. ) continue;
        if(!pat_muon->innerTrack().isNonnull()) continue;
        if(!pat_muon->globalTrack().isNonnull()) continue;
        if(pat_muon->globalTrack()->normalizedChi2() > 10. ) continue;
        if(pat_muon->globalTrack()->hitPattern().numberOfValidMuonHits() <= 0) continue;
        if(pat_muon->numberOfMatchedStations() <= 1) continue;
        if(pat_muon->innerTrack()->hitPattern().numberOfValidPixelHits() <= 0) continue;
        if(pat_muon->innerTrack()->hitPattern().trackerLayersWithMeasurement() <= 5) continue;
        if(pat_muon->isolationR03().sumPt/pat_muon->pt()>0.05) continue;

        // cout<<" PAT muon #"<<j<<" pt="<<pat_muon->pt()<<" eta="<<pat_muon->eta()<<" "<< pat_muon->isLooseMuon()<<endl;
        if(verbose) cout<<"JM PATSEL muon #"<<j<<" pt="<<pat_muon->pt()<<" eta="<<pat_muon->eta()<<" "<< pat_muon->isLooseMuon()<<endl;
        nPatMuonReco++;
    }
  
    Ptr<pat::Muon> pat_muon1;
    Ptr<pat::Muon> pat_muon2;
    float mass=0;
    float pt=0;

    if( nPatMuonReco>=2 ){
        
        float smallerDM=999.;
        float Zmass=91.19;
        Ptr<pat::Muon> pat_muontmp1;
        Ptr<pat::Muon> pat_muontmp2;
  
        for(unsigned int j = 0 ; j < patMuons->size(); j++){
            
            pat_muontmp1= patMuons->ptrAt(j);      

            if(!pat_muontmp1->isLooseMuon() || pat_muontmp1->pt()<10. ) continue;
            if(!pat_muontmp1->innerTrack().isNonnull()) continue;
            if(!pat_muontmp1->globalTrack().isNonnull()) continue;
            if(pat_muontmp1->globalTrack()->normalizedChi2() > 10. ) continue;
            if(pat_muontmp1->globalTrack()->hitPattern().numberOfValidMuonHits() <= 0) continue;
            if(pat_muontmp1->numberOfMatchedStations() <= 1) continue;
            if(pat_muontmp1->innerTrack()->hitPattern().numberOfValidPixelHits() <= 0) continue;
            if(pat_muontmp1->innerTrack()->hitPattern().trackerLayersWithMeasurement() <= 5) continue;
            if(pat_muontmp1->isolationR03().sumPt/pat_muontmp1->pt()>0.05) continue;
            
            
            for(unsigned int i = j+1 ; i < patMuons->size(); i++){
                
                pat_muontmp2= patMuons->ptrAt(i);
                if(!pat_muontmp2->isLooseMuon() || pat_muontmp2->pt()<10. ) continue;
                if(!pat_muontmp2->innerTrack().isNonnull()) continue;
                if(!pat_muontmp2->globalTrack().isNonnull()) continue;
                if(pat_muontmp2->globalTrack()->normalizedChi2() > 10. ) continue;
                if(pat_muontmp2->globalTrack()->hitPattern().numberOfValidMuonHits() <= 0) continue;
                if(pat_muontmp2->numberOfMatchedStations() <= 1) continue;
                if(pat_muontmp2->innerTrack()->hitPattern().numberOfValidPixelHits() <= 0) continue;
                if(pat_muontmp2->innerTrack()->hitPattern().trackerLayersWithMeasurement() <= 5) continue;
                if(pat_muontmp2->isolationR03().sumPt/pat_muontmp2->pt()>0.05) continue;
                
                XYZTLorentzVector p4Sum;
                p4Sum += pat_muontmp1->p4();
                p4Sum += pat_muontmp2->p4();
                
                double masstmp = p4Sum.M();
                double pttmp = p4Sum.pt();
                double DM=fabs(masstmp-Zmass);
                if(DM<smallerDM){
                    pat_muon1=pat_muontmp1;
                    pat_muon2=pat_muontmp2;
                    mass=masstmp;
                    pt=pttmp;
                    smallerDM=DM;
                }
                //cout<<" JM In Loop dimuon invariant mass:"<<masstmp<<" "<<DM<<" "<< endl;
            }
        }    
        
        if(verbose) cout<<"JM dimuon invariant mass:"<<mass<< endl;
        myInfo.diMuonMass=mass;
        myInfo.diMuonPt=pt;
    }
    bool pass=false;
    if( mass>50 && mass<1000 && nPatMuonReco>=2 ){
        if(verbose) cout<< "JM DIMUON CANDIDATE" << endl;
        pass=true;
    }
    
    if(pass){ 
        cout <<"pass"<< endl;
        cout<<" PV size:"<< primaryVertices->size()<<" NoMu size:"<< primaryVerticesNoMu->size()<<" WithMu size:"<< primaryVerticesWithMu->size()<< endl;
        
        myInfo.nTracks=0;
        myInfo.nTracksNoMu=0;
        
        for( unsigned int i = 0 ; i < primaryVerticesWithMu->size() ; i++ ) {
            Ptr<reco::Vertex> v=primaryVerticesWithMu->ptrAt(i);         
            myInfo.nTracks+=v->nTracks();
            if(verbose2) std::cout << "With Mu vtx "<< i 
                                   << " n " << std::setw(3) << v->nTracks()
                                   << " z "  << std::setw(6) << v->position().z() 
                                   << " dz " << std::setw(6) << v->zError()
                                   << std::endl;
        }   
        
        for( unsigned int i = 0 ; i < primaryVerticesNoMu->size() ; i++ ) {
            Ptr<reco::Vertex> v=primaryVerticesNoMu->ptrAt(i);         
            
            myInfo.nTracksNoMu+=v->nTracks();
            if(verbose2) std::cout << "NoMu vtx "<< i 
                                   << " n " << std::setw(3) << v->nTracks()
                                   << " z "  << std::setw(6) << v->position().z() 
                                   << " dz " << std::setw(6) << v->zError()
                                   << std::endl;
        }   
        
        
        
        const vector<edm::Ptr<reco::Vertex> > & vtxNoMu = primaryVerticesNoMu->ptrs();
        const vector<edm::Ptr<reco::Vertex> > & vtxWithMu = primaryVerticesWithMu->ptrs();
       

        myInfo.nVtx=vtxWithMu.size();
        myInfo.nVtxNoMu=vtxNoMu.size();
        
        bsInfo.x=beamSpot.x();
        bsInfo.y=beamSpot.y();
        bsInfo.z=beamSpot.z();
        bsInfo.sigmaz=sigmaZ;
        
        

        // test with muons:

        // double T1=sqrt(mu1->vx()*mu1->vx()+mu1->vy()*mu1->vy());
        // double T2=sqrt(mu2->vx()*mu2->vx()+mu2->vy()*mu2->vy()); 

        mu1Info.eta=pat_muon1->eta();
        mu1Info.pt=pat_muon1->pt();
        mu1Info.vx=pat_muon1->vx();
        mu1Info.vy=pat_muon1->vy();
        mu1Info.vz=pat_muon1->vz();
        mu1Info.dB=pat_muon1->dB(pat::Muon::BS3D);
        mu1Info.edB=pat_muon1->edB(pat::Muon::BS3D);
        mu1Info.sip=mu1Info.dB/mu1Info.edB;

        if(pat_muon1->innerTrack().isNonnull()){
            mu1Info.dzBS=pat_muon1->innerTrack()->dz(beamSpot);
            mu1Info.dxyBS=pat_muon1->innerTrack()->dxy(beamSpot);
        }

        mu2Info.eta=pat_muon2->eta();
        mu2Info.pt=pat_muon2->pt();
        mu2Info.vx=pat_muon2->vx();
        mu2Info.vy=pat_muon2->vy();
        mu2Info.vz=pat_muon2->vz();
        mu2Info.dB=pat_muon2->dB(pat::Muon::BS3D);
        mu2Info.edB=pat_muon2->edB(pat::Muon::BS3D);
        mu2Info.sip=mu2Info.dB/mu2Info.edB;

        if(pat_muon2->innerTrack().isNonnull()){
            mu2Info.dzBS=pat_muon2->innerTrack()->dz(beamSpot);
            mu2Info.dxyBS=pat_muon2->innerTrack()->dxy(beamSpot);
        }

        printMuon(mu1Info,"mu1");
        printMuon(mu2Info,"mu2");
        
        Ptr<reco::Vertex> pvxWithMu = vertexSelector_->select( pat_muon1,pat_muon2, vtxWithMu, *vertexCandidateMap, beamSpot); 
        float valWithMu;
        vertexSelector_->getInfoFromLastSelection(valWithMu);
        vector<float> allWithMu;
        vertexSelector_->getAllInfoFromLastSelection(allWithMu);

        Ptr<reco::Vertex> pvxNoMu = vertexSelector_->select( pat_muon1,pat_muon2, vtxNoMu, *vertexCandidateMapNoMu, beamSpot); 
        float valNoMu;
        vertexSelector_->getInfoFromLastSelection(valNoMu);
        vector<float> allNoMu;
        vertexSelector_->getAllInfoFromLastSelection(allNoMu);
        
        myInfo.zChosenNoMu=pvxNoMu->z();
        myInfo.zChosenWithMu=pvxWithMu->z();

        myInfo.ptChosenWithMu=pvxWithMu->p4().pt();
        myInfo.nTracksChosenWithMu=pvxWithMu->nTracks();

        myInfo.ptChosenNoMu=pvxNoMu->p4().pt();
        myInfo.nTracksChosenNoMu=pvxNoMu->nTracks();
       
        myInfo.mvaProbWithMu=valWithMu;
        myInfo.mvaProbNoMu=valNoMu;

        //        Ptr<reco::Vertex> pvxZeroWithMu = vertexSelectorZero_->select( pat_muon1,pat_muon2, vtxWithMu, *vertexCandidateMap, beamSpot); 
        //Ptr<reco::Vertex> pvxZeroNoMu = vertexSelectorZero_->select( pat_muon1,pat_muon2, vtxNoMu, *vertexCandidateMapNoMu, beamSpot); 
        Ptr<reco::Vertex> pvxZeroWithMu = vtxWithMu.at(0);
        Ptr<reco::Vertex> pvxZeroNoMu = vtxNoMu.at(0);

        myInfo.zZerothNoMu=pvxZeroNoMu->z();
        myInfo.ptZerothNoMu=pvxZeroNoMu->p4().pt();
        myInfo.nTracksZerothNoMu=pvxZeroNoMu->nTracks();

        myInfo.zZerothWithMu=pvxZeroWithMu->z();
        myInfo.ptZerothWithMu=pvxZeroWithMu->p4().pt();
        myInfo.nTracksZerothWithMu=pvxZeroWithMu->nTracks();

        if(allWithMu.size()==6 && allNoMu.size()==6){
            myInfo.logSumPt2WithMu=allWithMu.at(1);
            myInfo.logSumPt2NoMu=allNoMu.at(1);
            myInfo.ptBalWithMu=allWithMu.at(2);
            myInfo.ptBalNoMu=allNoMu.at(2);
            myInfo.ptAsymWithMu=allWithMu.at(3);
            myInfo.ptAsymNoMu=allNoMu.at(3);
            myInfo.dZ1WithMu=allWithMu.at(4);
            myInfo.dZ1NoMu=allNoMu.at(4);
            myInfo.dZ2WithMu=allWithMu.at(5);
            myInfo.dZ2NoMu=allNoMu.at(5);
        }

        myInfo.probWithMu=(1.+(-0.344)-(-0.091)+(-0.234)-(-0.186))+(-0.344)*valWithMu+(-0.091)*valWithMu*valWithMu+(-0.234)*valWithMu*valWithMu*valWithMu+(-0.186)*valWithMu*valWithMu*valWithMu*valWithMu;
        myInfo.probNoMu=(1.+(-0.344)-(-0.091)+(-0.234)-(-0.186))+(-0.344)*valNoMu+(-0.091)*valNoMu*valNoMu+(-0.234)*valNoMu*valNoMu*valNoMu+(-0.186)*valNoMu*valNoMu*valNoMu*valNoMu;

        int muonsConsistentWithMu;
        int muonsConsistentNoMu;
        float dz=0.2;
        if(relaxDZ) dz=50;
        double dZFoundNoMu=-999.;
        double dZFoundWithMu=-999.;
        
        int iFromMuonsWithMu=getRecoWithMuonsVertexIndex( vtxWithMu, pat_muon1, pat_muon2,muonsConsistentWithMu,dZFoundWithMu,dz);       
        
        Ptr<reco::Vertex> vFromMuonsWithMuon;
        if(iFromMuonsWithMu!=-1) vFromMuonsWithMuon=vtxWithMu[iFromMuonsWithMu];
        
        int iFromMuonsNoMu=getRecoWithMuonsVertexIndex( vtxNoMu, pat_muon1, pat_muon2,muonsConsistentNoMu ,dZFoundNoMu,dz);       
      
        myInfo.dzMatchingWithMu=dZFoundWithMu;
        myInfo.dzMatchingNoMu=dZFoundNoMu;


        Ptr<reco::Vertex> vFromMuonsNoMuon;
        if(iFromMuonsNoMu!=-1)vFromMuonsNoMuon=vtxNoMu[iFromMuonsNoMu];
        
        myInfo.zTrue = zTrueVtxIdx;
        
        if(trueRecoVtxIdx!=-1){
            myInfo.zRecoTrue=vertexTrueReco->z();
            myInfo.ptRecoTrue=vertexTrueReco->p4().pt();
            myInfo.nTracksRecoTrue=vertexTrueReco->nTracks();

        }else{
            myInfo.zRecoTrue=-999;
            myInfo.ptRecoTrue=-999;
        }
        if(nextTrueRecoVtxIdx!=-1) myInfo.zNextRecoTrue=nextVertexTrueReco->z();
        else  myInfo.zNextRecoTrue=-999;

        if(trueRecoVtxIdxNoMu !=-1){
            myInfo.zRecoTrueNoMu=vertexTrueRecoNoMu->z(); 
            myInfo.ptRecoTrueNoMu=vertexTrueRecoNoMu->p4().pt();
            myInfo.nTracksRecoTrueNoMu=vertexTrueRecoNoMu->nTracks();
        }else{
            myInfo.zRecoTrueNoMu=-999;
            myInfo.ptRecoTrueNoMu=-999;
        }
        if(nextTrueRecoVtxIdxNoMu!=-1) myInfo.zNextRecoTrueNoMu=nextVertexTrueRecoNoMu->z();
        else myInfo.zNextRecoTrueNoMu=-999;

        if(iFromMuonsWithMu!=-1) myInfo.zRecoFromMu=vFromMuonsWithMuon->z();
        if(iFromMuonsNoMu!=-1) myInfo.zRecoFromMuNoMu=vFromMuonsNoMuon->z();    

        myInfo.muonsConsistentNoMu=muonsConsistentNoMu;
        myInfo.muonsConsistentWithMu=muonsConsistentWithMu;
        
        vtxTree->Fill();

        // Right vertex:
        unsigned int trueVtxIndex = iFromMuonsNoMu;
        int trueVtxSortedIndexI = vertexSelector_->getSortedIndexFromLastSelection( trueVtxIndex );
       
        
        // Fill Signal Info
        if( trueVtxSortedIndexI!=-1) {
            unsigned int trueVtxSortedIndex= (unsigned int) trueVtxSortedIndexI;
            vector<float> allTrue;
            vertexSelector_->getAllInfoFromLastSelectionForVtxIdx( allTrue , trueVtxSortedIndex );
            if(allTrue.size()==4){
                sigInfo.LogSumPt2 = allTrue[1];
                sigInfo.PtBal   = allTrue[2];
                sigInfo.PtAsym  = allTrue[3];
                sigInfo.NConv  =  0;
                sigInfo.PullConv  =  0;
                sigInfo.DZtrue=primaryVerticesNoMu->ptrAt(trueVtxIndex)->position().z()-primaryVerticesWithMu->ptrAt(iFromMuonsWithMu)->position().z();
            }
        }

        signalTree->Fill();

        vector<int>	pvVecNoTrue;
        for( unsigned int i = 0 ; i < primaryVerticesNoMu->size() ; i++ ) {
            if( i != trueVtxIndex ) { pvVecNoTrue.push_back( i ); }
        }

        int irand = -999;
        if( pvVecNoTrue.size() > 1 ) { irand = rand() % pvVecNoTrue.size(); }

        int randVtxIndex = -999;
        if( irand != -999 ) { randVtxIndex = pvVecNoTrue[irand]; }

        int randVtxSortedIndexI = vertexSelector_->getSortedIndexFromLastSelection( randVtxIndex );

        if( randVtxSortedIndexI!=-1) {
            unsigned int randVtxSortedIndex=(unsigned int) randVtxSortedIndexI;

            vector<float> allRand;
            vertexSelector_->getAllInfoFromLastSelectionForVtxIdx( allRand , randVtxSortedIndex );
            if(allRand.size()==4){
                bkgInfo.LogSumPt2 = allRand[1];
                bkgInfo.PtBal   = allRand[2];
                bkgInfo.PtAsym  = allRand[3];
                bkgInfo.NConv  =  0;
                bkgInfo.PullConv  =  0;
                bkgInfo.DZtrue=primaryVerticesNoMu->ptrAt(randVtxIndex)->position().z()-primaryVerticesWithMu->ptrAt(iFromMuonsWithMu)->position().z();
            }
        }
        backgroundTree->Fill();
        
    }


}


void
outputAnalyzerDiMuHighMass::beginJob()
{

    vtxTree = fs_->make<TTree>( "vtxTree", "vtxTree" );
       
        
    vtxTree->Branch( "BSx",&bsInfo.x,"BSx/F");
    vtxTree->Branch( "BSy",&bsInfo.y,"BSy/F");
    vtxTree->Branch( "BSz",&bsInfo.z,"BSz/F");
    vtxTree->Branch( "BSsigmaz",&bsInfo.sigmaz,"BSsigmaz/F");
    vtxTree->Branch( "mcWeight", &myInfo.mcWeight, "mcWeight/I" );
    vtxTree->Branch( "nVtx", &myInfo.nVtx, "nVtx/I" );
    vtxTree->Branch( "nVtxNoMu", &myInfo.nVtxNoMu, "nVtxNoMu/I" );
    vtxTree->Branch( "nPU", &myInfo.nPU, "nPU/I" );
    vtxTree->Branch( "nTracks", &myInfo.nTracks, "nTracks/I" );
    vtxTree->Branch( "nTracksNoMu", &myInfo.nTracksNoMu, "nTracksNoMu/I" );
    vtxTree->Branch( "zTrue", &myInfo.zTrue, "zTrue/F" );
    vtxTree->Branch( "zRecoTrue", &myInfo.zRecoTrue, "zRecoTrue/F" );
    vtxTree->Branch( "zRecoTrueNoMu", &myInfo.zRecoTrueNoMu, "zRecoTrueNoMu/F" );
    vtxTree->Branch( "zNextRecoTrue", &myInfo.zNextRecoTrue, "zNextRecoTrue/F" );
    vtxTree->Branch( "zNextRecoTrueNoMu", &myInfo.zNextRecoTrueNoMu, "zNextRecoTrueNoMu/F" );
    vtxTree->Branch( "zRecoFromMu", &myInfo.zRecoFromMu, "zRecoFromMu/F" );
    vtxTree->Branch( "zRecoFromMuNoMu", &myInfo.zRecoFromMuNoMu, "zRecoFromMuNoMu/F" );
    vtxTree->Branch( "muonsConsistentNoMu", &myInfo.muonsConsistentNoMu, "muonsConsistentNoMu/I" );
    vtxTree->Branch( "muonsConsistentWithMu", &myInfo.muonsConsistentWithMu, "muonsConsistentWithMu/I" );
    vtxTree->Branch( "diMuonMass", &myInfo.diMuonMass,"diMuonMass/F");
    vtxTree->Branch( "diMuonPt", &myInfo.diMuonPt,"diMuonPt/F");
    vtxTree->Branch( "zChosenWithMu",&myInfo.zChosenWithMu,"zChosenWithMu/F");
    vtxTree->Branch( "zChosenNoMu",&myInfo.zChosenNoMu,"zChosenNoMu/F");
    vtxTree->Branch( "zZerothWithMu",&myInfo.zZerothWithMu,"zZerothWithMu/F");
    vtxTree->Branch( "zZerothNoMu",&myInfo.zZerothNoMu,"zZerothNoMu/F");
    vtxTree->Branch( "mvaProbWithMu",&myInfo.mvaProbWithMu,"mvaProbWithMu/F");
    vtxTree->Branch( "mvaProbNoMu",&myInfo.mvaProbNoMu,"mvaProbNoMu/F");
    vtxTree->Branch( "probWithMu",&myInfo.probWithMu,"probWithMu/F");
    vtxTree->Branch( "probNoMu",&myInfo.probNoMu,"probNoMu/F");
    vtxTree->Branch( "logSumPt2WithMu",&myInfo.logSumPt2WithMu,"logSumPt2WithMu/F");
    vtxTree->Branch( "logSumPt2NoMu",&myInfo.logSumPt2NoMu,"logSumPt2NoMu/F");
    vtxTree->Branch( "ptBalWithMu",&myInfo.ptBalWithMu,"ptBalWithMu/F");
    vtxTree->Branch( "ptBalNoMu",&myInfo.ptBalNoMu,"ptBalNoMu/F");
    vtxTree->Branch( "ptAsymWithMu",&myInfo.ptAsymWithMu,"ptAsymWithMu/F");
    vtxTree->Branch( "ptAsymNoMu",&myInfo.ptAsymNoMu,"ptAsymNoMu/F");
    vtxTree->Branch( "dZ1WithMu",&myInfo.dZ1WithMu,"dZ1WithMu/F");
    vtxTree->Branch( "dZ1NoMu",&myInfo.dZ1NoMu,"dZ1NoMu/F");
    vtxTree->Branch( "dZ2WithMu",&myInfo.dZ2WithMu,"dZ2WithMu/F");
    vtxTree->Branch( "dZ2NoMu",&myInfo.dZ2NoMu,"dZ2NoMu/F");


    vtxTree->Branch( "nTracksRecoTrue",&myInfo.nTracksRecoTrue,"nTracksRecoTrue/I");
    vtxTree->Branch( "nTracksRecoTrueNoMu",&myInfo.nTracksRecoTrueNoMu,"nTracksRecoTrueNoMu/I");
    vtxTree->Branch( "nTracksChosenWithMu",&myInfo.nTracksChosenWithMu,"nTracksChosenWithMu/I");
    vtxTree->Branch( "nTracksChosenNoMu",&myInfo.nTracksChosenNoMu,"nTracksChosenNoMu/I");
    vtxTree->Branch( "nTracksZerothWithMu",&myInfo.nTracksZerothWithMu,"nTracksZerothWithMu/I");
    vtxTree->Branch( "nTracksZerothNoMu",&myInfo.nTracksZerothNoMu,"nTracksZerothNoMu/I");
    vtxTree->Branch( "ptRecoTrue",&myInfo.ptRecoTrue,"ptRecoTrue/F");
    vtxTree->Branch( "ptRecoTrueNoMu",&myInfo.ptRecoTrueNoMu,"ptRecoTrueNoMu/F");
    vtxTree->Branch( "ptChosenWithMu",&myInfo.ptChosenWithMu,"ptChosenWithMu/F");
    vtxTree->Branch( "ptChosenNoMu",&myInfo.ptChosenNoMu,"ptChosenNoMu/F");
    vtxTree->Branch( "ptZerothWithMu",&myInfo.ptZerothWithMu,"ptZerothWithMu/F");
    vtxTree->Branch( "ptZerothNoMu",&myInfo.ptZerothNoMu,"ptZerothNoMu/F");

    vtxTree->Branch( "dzMatchingWithMu",&myInfo.dzMatchingWithMu,"dzMatchingWithMu/F");
    vtxTree->Branch( "dzMatchingNoMu",&myInfo.dzMatchingNoMu,"dzMatchingNoMu/F");

    vtxTree->Branch( "mu1eta",&mu1Info.eta,"mu1eta/F");
    vtxTree->Branch( "mu1pt",&mu1Info.pt,"mu1pt/F");
    vtxTree->Branch( "mu1vx",&mu1Info.vx,"mu1vx/F");
    vtxTree->Branch( "mu1vy",&mu1Info.vy,"mu1vy/F");
    vtxTree->Branch( "mu1vz",&mu1Info.vz,"mu1vz/F");
    vtxTree->Branch( "mu1dB",&mu1Info.dB,"mu1dB/F");
    vtxTree->Branch( "mu1edB",&mu1Info.edB,"mu1edB/F");
    vtxTree->Branch( "mu1sip",&mu1Info.sip,"mu1sip/F");
    vtxTree->Branch( "mu1dzBS",&mu1Info.dzBS,"mu1dzBS/F");
    vtxTree->Branch( "mu1dxyBS",&mu1Info.dxyBS,"mu1dxyBS/F");
    
    vtxTree->Branch( "mu2eta",&mu2Info.eta,"mu2eta/F");
    vtxTree->Branch( "mu2pt",&mu2Info.pt,"mu2pt/F");
    vtxTree->Branch( "mu2vx",&mu2Info.vx,"mu2vx/F");
    vtxTree->Branch( "mu2vy",&mu2Info.vy,"mu2vy/F");
    vtxTree->Branch( "mu2vz",&mu2Info.vz,"mu2vz/F");
    vtxTree->Branch( "mu2dB",&mu2Info.dB,"mu2dB/F");
    vtxTree->Branch( "mu2edB",&mu2Info.edB,"mu2edB/F");
    vtxTree->Branch( "mu2sip",&mu2Info.sip,"mu2sip/F");
    vtxTree->Branch( "mu2dzBS",&mu2Info.dzBS,"mu2dzBS/F");
    vtxTree->Branch( "mu2dxyBS",&mu2Info.dxyBS,"mu2dxyBS/F");
    
    signalTree = fs_->make<TTree>( "signalTree", "per-diphoton tree" );
        
    signalTree->Branch( "BSx",&bsInfo.x,"BSx/F");
    signalTree->Branch( "BSy",&bsInfo.y,"BSy/F");
    signalTree->Branch( "BSz",&bsInfo.z,"BSz/F");
    signalTree->Branch( "BSsigmaz",&bsInfo.sigmaz,"BSsigmaz/F");
    signalTree->Branch( "mcWeight", &myInfo.mcWeight, "mcWeight/I" );
    signalTree->Branch( "nPU", &myInfo.nPU, "nPU/I" );
    signalTree->Branch( "mu1eta",&mu1Info.eta,"mu1eta/F");
    signalTree->Branch( "mu1pt",&mu1Info.pt,"mu1pt/F");
    signalTree->Branch( "mu2eta",&mu2Info.eta,"mu2eta/F");
    signalTree->Branch( "mu2pt",&mu2Info.pt,"mu2pt/F");
    signalTree->Branch( "diMuonMass", &myInfo.diMuonMass,"diMuonMass/F");
    signalTree->Branch( "nvertex", &myInfo.nVtxNoMu, "nvertex/I" );
    signalTree->Branch( "ndipho", &sigInfo.ndipho, "ndipho/I" );
    signalTree->Branch( "dipho_index", &sigInfo.dipho_index, "dipho_index/I" );
    signalTree->Branch( "LogSumPt2", &sigInfo.LogSumPt2, "LogSumPt2/F" );
    signalTree->Branch( "PtBal", &sigInfo.PtBal, "PtBal/F" );
    signalTree->Branch( "PtAsym", &sigInfo.PtAsym, "PtAsym/F" );
    signalTree->Branch( "DZtrue", &sigInfo.DZtrue, "DZtrue/F" );
    signalTree->Branch( "NConv", &sigInfo.NConv, "NConv/F" );
    signalTree->Branch( "PullConv", &sigInfo.PullConv, "PullConv/F" );

    backgroundTree = fs_->make<TTree>( "backgroundTree", "per-diphoton tree" );
 
    backgroundTree->Branch( "BSx",&bsInfo.x,"BSx/F");
    backgroundTree->Branch( "BSy",&bsInfo.y,"BSy/F");
    backgroundTree->Branch( "BSz",&bsInfo.z,"BSz/F");
    backgroundTree->Branch( "BSsigmaz",&bsInfo.sigmaz,"BSsigmaz/F");
    backgroundTree->Branch( "mcWeight", &myInfo.mcWeight, "mcWeight/I" );
    backgroundTree->Branch( "nPU", &myInfo.nPU, "nPU/I" );
    backgroundTree->Branch( "mu1eta",&mu1Info.eta,"mu1eta/F");
    backgroundTree->Branch( "mu1pt",&mu1Info.pt,"mu1pt/F");
    backgroundTree->Branch( "mu2eta",&mu2Info.eta,"mu2eta/F");
    backgroundTree->Branch( "mu2pt",&mu2Info.pt,"mu2pt/F");
    backgroundTree->Branch( "diMuonMass", &myInfo.diMuonMass,"diMuonMass/F");
    backgroundTree->Branch( "nvertex",  &myInfo.nVtxNoMu, "nvertex/I" );
    backgroundTree->Branch( "ndipho", &bkgInfo.ndipho, "ndipho/I" );
    backgroundTree->Branch( "dipho_index", &bkgInfo.dipho_index, "dipho_index/I" );
    backgroundTree->Branch( "LogSumPt2", &bkgInfo.LogSumPt2, "LogSumPt2/F" );
    backgroundTree->Branch( "PtBal", &bkgInfo.PtBal, "PtBal/F" );
    backgroundTree->Branch( "PtAsym", &bkgInfo.PtAsym, "PtAsym/F" );
    backgroundTree->Branch( "DZtrue", &bkgInfo.DZtrue, "DZtrue/F" );
    backgroundTree->Branch( "NConv", &bkgInfo.NConv, "NConv/F" );
    backgroundTree->Branch( "PullConv", &bkgInfo.PullConv, "PullConv/F" );

}

void
outputAnalyzerDiMuHighMass::endJob()
{
    cout<<" No Match:"<<noMatch<<" Pass:"<<nPass<< endl;
}

void
outputAnalyzerDiMuHighMass::initEventStructure()
{

    bsInfo.x=-999;
    bsInfo.y=-999;
    bsInfo.z=-999;
    bsInfo.sigmaz=-999;
        
    myInfo.mcWeight=1;
    myInfo.diMuonMass=-999;
    myInfo.diMuonPt=-999;

    myInfo.nVtx=0;
    myInfo.nVtxNoMu=0;
    myInfo.nPU=0;
    
    myInfo.zTrue =-999;
        
    myInfo.zRecoTrue=-999;
    myInfo.zRecoTrueNoMu=-999;
    myInfo.zNextRecoTrue=-999;
    myInfo.zNextRecoTrueNoMu=-999;
    myInfo.zRecoFromMu=-999;
    myInfo.zRecoFromMuNoMu=-999;
    myInfo.zChosenNoMu=-999;
    myInfo.zChosenWithMu=-999;
    myInfo.zZerothNoMu=-999;
    myInfo.zZerothWithMu=-999;
    myInfo.mvaProbWithMu=-999;
    myInfo.mvaProbNoMu=-999;
    myInfo.probWithMu=-999;
    myInfo.probNoMu=-999;
    myInfo.dzMatchingWithMu=-999;
    myInfo.dzMatchingNoMu=-999;

    myInfo.logSumPt2WithMu=-999;
    myInfo.logSumPt2NoMu=-999;
    myInfo.ptBalWithMu=-999;
    myInfo.ptBalNoMu=-999;
    myInfo.ptAsymWithMu=-999;
    myInfo.ptAsymNoMu=-999;
    myInfo.dZ1WithMu=-999;
    myInfo.dZ1NoMu=-999;
    myInfo.dZ2WithMu=-999;
    myInfo.dZ2NoMu=-999;
    

    myInfo.nTracksRecoTrue=-999;
    myInfo.nTracksRecoTrueNoMu=-999;
    myInfo.nTracksChosenWithMu=-999;
    myInfo.nTracksChosenNoMu=-999;
    myInfo.nTracksZerothWithMu=-999;
    myInfo.nTracksZerothNoMu=-999;
    myInfo.ptRecoTrue=-999;
    myInfo.ptRecoTrueNoMu=-999;
    myInfo.ptChosenWithMu=-999;
    myInfo.ptChosenNoMu=-999;
    myInfo.ptZerothWithMu=-999;
    myInfo.ptZerothNoMu=-999;
    

    mu1Info.eta=-999;
    mu1Info.pt=-999;
    mu1Info.vx=-999;
    mu1Info.vy=-999;
    mu1Info.vz=-999;
    mu1Info.dB=-999;
    mu1Info.edB=-999;
    mu1Info.sip=-999;
    mu1Info.dzBS=-999;
    mu1Info.dxyBS=-999;

    mu2Info.eta=-999;
    mu2Info.pt=-999;
    mu2Info.vx=-999;
    mu2Info.vy=-999;
    mu2Info.vz=-999;
    mu2Info.dB=-999;
    mu2Info.edB=-999;
    mu2Info.sip=-999;
    mu2Info.dzBS=-999;
    mu2Info.dxyBS=-999;
    
    myInfo.muonsConsistentNoMu=0;
    myInfo.muonsConsistentWithMu=0;

    sigInfo.nvertex = -999;
    sigInfo.ndipho = -999;
    sigInfo.dipho_index = -999;

    sigInfo.LogSumPt2  = -999;
    sigInfo.PtBal  = -999;
    sigInfo.PtAsym  = -999;
    sigInfo.DZtrue  = -999;
    sigInfo.NConv  = -999;
    sigInfo.PullConv  = -999;


    bkgInfo.nvertex = -999;
    bkgInfo.ndipho = -999;
    bkgInfo.dipho_index = -999;

    bkgInfo.LogSumPt2  = -999;
    bkgInfo.PtBal  = -999;
    bkgInfo.PtAsym  = -999;
    bkgInfo.DZtrue  = -999;
    bkgInfo.NConv  = -999;
    bkgInfo.PullConv  = -999;

}

int outputAnalyzerDiMuHighMass::getMatchingVertexIndexNoMu( const vector<edm::Ptr<reco::Vertex> > & verticesNoMu, Ptr<reco::Vertex> vertexTrueWithMu, double dzMatch){

    int nTracksTrue=vertexTrueWithMu->nTracks();
    int nTracksToFind=nTracksTrue-2;
    if(nTracksToFind<=0) return -1;

    float dZMin=999;
    int iNoMu=-1;
    for( unsigned int iv = 0; iv < verticesNoMu.size(); iv++ ) {
        if(verticesNoMu[iv]->nTracks()!=(unsigned) nTracksToFind){
            //cout<<" Not passing nTracks cut "<<verticesNoMu[iv]->nTracks()<<" "<<nTracksToFind<< endl;
            continue;
        }
        float dz=fabs(verticesNoMu[iv]->z()-vertexTrueWithMu->z());
        if(dz<dZMin){
            dZMin=dz;
            iNoMu=iv;
        }        
        cout<<" zNoMu="<<verticesNoMu[iv]->z()<<" zMu="<<vertexTrueWithMu->z()<<" dZMin="<<dZMin<<endl;
    }
  
    if(dZMin<dzMatch) return iNoMu;

    cout<<" Not passing MATCH cut "<<dZMin<< endl;
    return -1;
}
    


int outputAnalyzerDiMuHighMass::getRecoWithMuonsVertexIndex( const vector<edm::Ptr<reco::Vertex> > & vertices, Ptr<pat::Muon> mu1, Ptr<pat::Muon> mu2, int& muonsConsistent, double &dzFound , double dzMatch)
{
    muonsConsistent=1;

    double IP1=mu1->vz();
    double IP2=mu2->vz(); 

    double average=0.5*(IP1+IP2);
        
    double mu1sip=fabs(mu1->dB(pat::Muon::BS3D)/mu1->edB(pat::Muon::BS3D));
    double mu2sip=fabs(mu2->dB(pat::Muon::BS3D)/mu2->edB(pat::Muon::BS3D));

    double w1=mu2sip/(mu1sip+mu2sip);
    double w2=mu1sip/(mu1sip+mu2sip);


    double weightedSum=(w1*IP1+w2*IP2)/(w1+w2);

    cout<<" IP1="<<IP1<<" IP2="<<IP2<<" weightedSum="<<weightedSum<<" average="<< average<<" w1="<<w1<<"w2="<<w2<<" "<<w1/(w1+w2)<<" "<<w2/(w1+w2)<<endl;
   
    
    int  ivMatch = 0;
    int  ivMatch1 = 0;
    int  ivMatch2 = 0;
    double dzMin = 999;
    double dzMin1 = 999;
    double dzMin2 = 999;
    
    for( unsigned int iv = 0; iv < vertices.size(); iv++ ) {
        double dz = fabs( vertices[iv]->z() - average );
        double dz1 = fabs( vertices[iv]->z() - IP1 );
        double dz2 = fabs( vertices[iv]->z() - IP2 );
        if( dz < dzMin ) {
            ivMatch = iv;
            dzMin   = dz;
        } 
        if( dz1 < dzMin1 ) {
            ivMatch1 = iv;
            dzMin1   = dz1;
        }
        if( dz2 < dzMin2 ) {
            ivMatch2 = iv;
            dzMin2   = dz2;
        }
    }
    
    
    if(ivMatch!=ivMatch1 || ivMatch!=ivMatch2){
        muonsConsistent=0;
        cout<< "MUONS INCONSISTENT -- dzMin="<<dzMin<<" dzMin1="<<dzMin1<<" dzMin2="<<dzMin2<< endl;
        
        if(dzMin1<dzMin2) {
            dzMin=dzMin1;
            ivMatch=ivMatch1;
        }else{
            dzMin=dzMin2;
            ivMatch=ivMatch2;
        }
    }

    dzFound=dzMin;
    if( dzMin < dzMatch ) { return ivMatch; }
    return -1;
}


int outputAnalyzerDiMuHighMass::getClosestVertexIndex( const vector<edm::Ptr<reco::Vertex> > & vertices, Ptr<reco::Vertex> vertexTrue, double dzMatch){


    float dZMin=999;
    int iVtx=-1;
    for( unsigned int iv = 0; iv < vertices.size(); iv++ ) {

        float dz=fabs(vertices[iv]->z()-vertexTrue->z());
        if(dz<dZMin){
            dZMin=dz;
            iVtx=iv;
        }        
    }
    
    if(dZMin<dzMatch) return iVtx;
    
    cout<<" Not passing MATCH cut "<<dZMin<< endl;
    return -1;
}
    

void outputAnalyzerDiMuHighMass::getMCTruthVertexCoord(  const std::vector<edm::Ptr<reco::GenParticle>>& gens , float& vx, float&vy, float&vz )
{

    for( unsigned int genLoop = 0 ; genLoop < gens.size(); genLoop++ ) {
        
        if( fabs( gens[genLoop]->pdgId() ) < 10 || fabs( gens[genLoop]->pdgId() ) == 23 ) {
            //hardVertex.SetCoordinates( gens[genLoop]->vx(), gens[genLoop]->vy(), gens[genLoop]->vz() );
            vx=gens[genLoop]->vx(); vy=gens[genLoop]->vy(); vz=gens[genLoop]->vz();
            
            break;
        }
    }
}
int outputAnalyzerDiMuHighMass::getRecoClosestToTrueVertexIndex(  const std::vector<edm::Ptr<reco::GenParticle>>& gens , const vector<edm::Ptr<reco::Vertex> > & vertices, double dzMatch )
{

    reco::Vertex::Point hardVertex( 0, 0, 0 );

    for( unsigned int genLoop = 0 ; genLoop < gens.size(); genLoop++ ) {

        if( fabs( gens[genLoop]->pdgId() ) < 10 || fabs( gens[genLoop]->pdgId() ) == 23 ) {
            hardVertex.SetCoordinates( gens[genLoop]->vx(), gens[genLoop]->vy(), gens[genLoop]->vz() );
            break;
        }
    }

    int  ivMatch = 0;
    double dzMin = 999;
    for( unsigned int iv = 0; iv < vertices.size(); iv++ ) {
        double dz = fabs( vertices[iv]->z() - hardVertex.z() );
        if( dz < dzMin ) {
            ivMatch = iv;
            dzMin   = dz;
        }
    }


    if( dzMin < dzMatch ) { return ivMatch; }

    return -1;
}

int outputAnalyzerDiMuHighMass::getRecoNextToClosestToTrueVertexIndex(  const std::vector<edm::Ptr<reco::GenParticle>>& gens , const vector<edm::Ptr<reco::Vertex> > & vertices,  Ptr<reco::Vertex>  closest , double dzMatch )
{

    reco::Vertex::Point hardVertex( 0, 0, 0 );

    for( unsigned int genLoop = 0 ; genLoop < gens.size(); genLoop++ ) {

        if( fabs( gens[genLoop]->pdgId() ) < 10 || fabs( gens[genLoop]->pdgId() ) == 23 ) {
            hardVertex.SetCoordinates( gens[genLoop]->vx(), gens[genLoop]->vy(), gens[genLoop]->vz() );
            break;
        }
    }
    double dzClosest= fabs( closest->z() - hardVertex.z() );

    int  ivMatch = 0;
    double dzMin = 999;
    for( unsigned int iv = 0; iv < vertices.size(); iv++ ) {
        double dz = fabs( vertices[iv]->z() - hardVertex.z() );
        if(dz<=dzClosest) continue;
        if( dz < dzMin ) {
            ivMatch = iv;
            dzMin   = dz;
        }
    }


    if( dzMin < dzMatch ) { return ivMatch; }

    return -1;
}


void
outputAnalyzerDiMuHighMass::printMuon( MuonInfo inf, string name ){

    cout<<name.c_str() <<
        ": eta="<< inf.eta<<
        " pt="<< inf.pt<<
        " vx="<< inf.vx<<
        " vy="<< inf.vy<<
        " vz="<< inf.vz<<
        " dB="<< inf.dB<<
        " edB="<< inf.edB<<
        " sip="<< inf.sip<<
        " dzBS="<< inf.dzBS<<
        " dxyBS="<< inf.dxyBS<<endl;

}
void
outputAnalyzerDiMuHighMass::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault( desc );
}


DEFINE_FWK_MODULE( outputAnalyzerDiMuHighMass );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

