// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TFile.h"
#include "Math/VectorUtil.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "TTree.h"

using namespace std;
using namespace edm;

class ElectronValidator : public EDAnalyzer
{
    public:
        explicit ElectronValidator( const edm::ParameterSet & );
        ~ElectronValidator();

    private:
        virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
        virtual void beginJob() override;
        //virtual void endJob() override;

        bool verbose_;
        bool applyCuts_;
        edm::EDGetTokenT<View<pat::Electron> > electronToken_;
        edm::EDGetTokenT<View<reco::Vertex> > vertexToken_;
        edm::EDGetTokenT<reco::ConversionCollection> convToken_;
        edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
        edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
        edm::EDGetTokenT<View<pat::PackedGenParticle> > genToken_;
        edm::EDGetTokenT<edm::ValueMap<bool> > eleMVAMediumIdMapToken_;
        edm::EDGetTokenT<edm::ValueMap<bool> > eleMVATightIdMapToken_;
        edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
        edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
        edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
        edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken_;

        EffectiveAreas _effectiveAreas;

        TTree *electronTree_;

        float _Rho;
        edm::InputTag rhoFixedGrid_;

        Int_t nElectrons_;
        Int_t evtnum_;

        std::vector<Float_t> pt_;
        std::vector<Float_t> eta_;
        std::vector<Float_t> phi_;

        std::vector<Float_t> mvaValue_;

        std::vector<Int_t> passLooseId_;
        std::vector<Int_t> passMediumId_;
        std::vector<Int_t> passTightId_;

        std::vector<Int_t> passMVAMediumId_;
        std::vector<Int_t> passMVATightId_;

        std::vector<Float_t> isolation_;
        std::vector<Float_t> charged_hadron_isolation_;
        std::vector<Float_t> neutral_hadron_isolation_;
        std::vector<Float_t> photon_isolation_;
        std::vector<Float_t> pile_up_;

        std::vector<Int_t> missing_inner_hit_;
        std::vector<Int_t> has_matched_conversion_; 

        std::vector<Int_t> genmatch_;

};

void ElectronValidator::beginJob(){
    //electron_validation_ = file->make<TH1F>("electron.root", "RECREATE");

}

ElectronValidator::ElectronValidator( const ParameterSet &iConfig ):
    electronToken_( consumes<View<pat::Electron> >( iConfig.getParameter<InputTag>( "electronTag" ) ) ),
    vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag>( "vertexTag" ) ) ),
    convToken_( consumes<reco::ConversionCollection>( iConfig.getParameter<InputTag>( "convTag" ) ) ),
    beamSpotToken_( consumes<reco::BeamSpot >( iConfig.getParameter<InputTag>( "beamSpotTag" ) ) ),
    mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap" ) ) ),
    genToken_(consumes<View<pat::PackedGenParticle > >( iConfig.getParameter<InputTag>( "genTag" ) ) ), 
    eleMVAMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVAMediumIdMap"))),
    eleMVATightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMVATightIdMap"))),
    eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
    eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
    eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
    eleHEEPIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap"))),
    _effectiveAreas( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath())
{
    applyCuts_ = iConfig.getUntrackedParameter<bool>( "ApplyCuts", true );
    verbose_ = iConfig.getUntrackedParameter<bool>( "verbose", false );
    //	eventrhoToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("Rho"));
    rhoFixedGrid_  = iConfig.getParameter<edm::InputTag>( "rhoFixedGridCollection" );

    edm::Service<TFileService> fs;
    electronTree_ = fs->make<TTree> ("ElectronTree", "Electron data");

    //  electronTree_->Branch("run"        ,  &run_     , "run/I");
    //  electronTree_->Branch("lumi"       ,  &lumi_    , "lumi/I");
    electronTree_->Branch("evtnum"     ,  &evtnum_  , "evtnum/I");

    electronTree_->Branch("nEle",  &nElectrons_ , "nEle/I");
    electronTree_->Branch("pt"  ,  &pt_    );
    electronTree_->Branch("eta" ,  &eta_ );
    electronTree_->Branch("phi" ,  &phi_ );

    electronTree_->Branch("mvaValue" ,  &mvaValue_ );

    electronTree_->Branch("passLooseId" ,  &passLooseId_ );
    electronTree_->Branch("passMediumId" ,  &passMediumId_ );
    electronTree_->Branch("passTightId"  ,  &passTightId_ );

    electronTree_->Branch("passMVAMediumId" ,  &passMVAMediumId_ );
    electronTree_->Branch("passMVATightId"  ,  &passMVATightId_ );

    electronTree_->Branch("isolation" , &isolation_);
    electronTree_->Branch("charged_hadron_isolation" , &charged_hadron_isolation_);
    electronTree_->Branch("neutral_hadron_isolation" , &neutral_hadron_isolation_);
    electronTree_->Branch("photon_isolation" , &photon_isolation_);
    electronTree_->Branch("pile_up" , &pile_up_);

    electronTree_->Branch("missing_inner_hit", &missing_inner_hit_);
    electronTree_->Branch("has_matched_conversion", &has_matched_conversion_);

    electronTree_->Branch("genmatch",&genmatch_);
}


ElectronValidator::~ElectronValidator(){}

void ElectronValidator::analyze( const Event &evt, const EventSetup & ){

    //using namespace edm;

    Handle<View<pat::Electron> >  pelectrons;
    evt.getByToken( electronToken_, pelectrons );
    //	const PtrVector<pat::Electron> pelectronPointers = pelectrons->ptrVector();

    _Rho = 0;
    Handle<double> rhoHandle;
    evt.getByLabel( rhoFixedGrid_, rhoHandle );
    _Rho = *rhoHandle;


    Handle<View<reco::Vertex> >  vtxs;
    evt.getByToken( vertexToken_, vtxs );
    //	const PtrVector<reco::Vertex> vertexPointers = vtxs->ptrVector();

    Handle<reco::ConversionCollection> convs;
    evt.getByToken( convToken_, convs );
    //		const PtrVector<reco::Conversion> convPointers = convs->ptrVector();

    Handle<reco::BeamSpot> recoBeamSpotHandle;
    evt.getByToken( beamSpotToken_, recoBeamSpotHandle );
    math::XYZPoint vertexPoint;
    if( recoBeamSpotHandle.isValid() )
    { vertexPoint = recoBeamSpotHandle->position(); }

    edm::Handle<edm::ValueMap<float> > mvaValues;
    evt.getByToken(mvaValuesMapToken_,mvaValues);

    Handle<View<pat::PackedGenParticle> > genParticle;
    evt.getByToken(genToken_,genParticle); 

    edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
    edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
    edm::Handle<edm::ValueMap<bool> > tight_id_decisions; 
    edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
    evt.getByToken(eleLooseIdMapToken_ ,loose_id_decisions);
    evt.getByToken(eleMediumIdMapToken_,medium_id_decisions);
    evt.getByToken(eleTightIdMapToken_,tight_id_decisions);
    evt.getByToken(eleHEEPIdMapToken_ ,heep_id_decisions);

    edm::Handle<edm::ValueMap<bool> > mediumMVA_wp;
    edm::Handle<edm::ValueMap<bool> > tightMVA_wp; 
    evt.getByToken(eleMVAMediumIdMapToken_,mediumMVA_wp);
    evt.getByToken(eleMVATightIdMapToken_,tightMVA_wp);

    vector<pat::PackedGenParticle> matchedGen;

       pt_.clear();
       eta_.clear();
       phi_.clear();
       mvaValue_.clear();
       passLooseId_.clear();
       passMediumId_.clear();
       passTightId_ .clear();
       passMVAMediumId_.clear();
       passMVATightId_.clear();
       isolation_.clear();
       charged_hadron_isolation_.clear();
       neutral_hadron_isolation_.clear();
       photon_isolation_.clear();
       pile_up_.clear();
       missing_inner_hit_.clear();
       has_matched_conversion_.clear(); 
       int nelectrons=0;  

       evtnum_ = evt.id().event();

    for( unsigned int elecIndex = 0; elecIndex < pelectrons->size(); elecIndex++ ) {
        Ptr<pat::Electron> pelec = pelectrons->ptrAt( elecIndex );
        flashgg::Electron felec = flashgg::Electron( *pelec );
        // double nontrigmva_ = -999999;
        double dR = 999;
        nelectrons++;
        for( unsigned int genIndex = 0; genIndex < genParticle->size(); genIndex++){
            Ptr<pat::PackedGenParticle> gen = genParticle->ptrAt(genIndex);
            //cout << ElectronValidator::GenElectronMatchedToBoson(*gen) << endl;
        if ( (*gen).energy() == 0) continue; 

        if ( fabs((*gen).pdgId()) != 11 || !((*gen).isPromptFinalState())) continue;

                double dRtmp = ROOT::Math::VectorUtil::DeltaR( pelec->p4(), gen->p4() );

                if( dRtmp < dR ){
                    dR = dRtmp;

                    Ptr<pat::PackedGenParticle> particle = gen;
                }


        } 

        int genmatch = 0;

        if( dR < 0.1 ) { genmatch = 1; }

        if (genmatch != 1) continue;

        genmatch_.push_back(genmatch);

        float Aeff = 0;

        float pelec_eta = fabs( pelec->superCluster()->eta() );

        double nontrigmva = (*mvaValues)[pelec];

        if(pelec->pt() < 5.0 )continue;

//        if (genmatch==1){
//        cout << "dR " << dR << endl;
//        cout << "match" << endl;
//        cout << "mvaValue " <<nontrigmva << endl;
//        cout << "pt " << pelec->pt() << endl;
       
        mvaValue_ .push_back( nontrigmva );
        pt_  .push_back( pelec->pt() );
        eta_ .push_back( pelec_eta );
        phi_ .push_back( pelec->phi() );

        bool passLooseId =  (*loose_id_decisions)[pelec];
        bool passMediumId = (*medium_id_decisions)[pelec];
        bool passTightId = (*tight_id_decisions)[pelec];

        bool passMVAMediumId = (*mediumMVA_wp)[pelec];
        bool passMVATightId = (*tightMVA_wp)[pelec];

        passLooseId_  .push_back(passLooseId);
        passMediumId_ .push_back(passMediumId);
        passTightId_  .push_back(passTightId);

        passMVAMediumId_ .push_back(passMVAMediumId);
        passMVATightId_  .push_back(passMVATightId);

        //        if( applyCuts_ && ( pelec_eta > 2.5 || ( pelec_eta > 1.442 && pelec_eta < 1.566 ) ) ) { continue; }

        const reco::GsfElectron::PflowIsolationVariables& pfIso = pelec->pfIsolationVariables();
        float chad = pfIso.sumChargedHadronPt;
        float nhad = pfIso.sumNeutralHadronEt;
        float pho = pfIso.sumPhotonEt;
        Aeff = _effectiveAreas.getEffectiveArea( pelec_eta );   
        //cout << "Aff " << Aeff << endl;
        float pile_up = _Rho*Aeff;
        //cout << "Rho " << _Rho << endl;
        //cout << "pile_up " << pile_up << endl; 
        float iso = chad + std::max(0.0f, nhad + pho - _Rho*Aeff);

        isolation_               .push_back(iso);
        charged_hadron_isolation_.push_back(chad);
        neutral_hadron_isolation_.push_back(nhad);
        photon_isolation_        .push_back(pho);    
        pile_up_                 .push_back(pile_up); 

        bool has_matched_conversion = ConversionTools::hasMatchedConversion(*pelec,convs,vertexPoint);
        int missing_inner_hit = pelec->gsfTrack()->hitPattern().numberOfHits( reco::HitPattern::MISSING_INNER_HITS ); 

        has_matched_conversion_.push_back(has_matched_conversion);
        missing_inner_hit_.push_back(missing_inner_hit);
//        }

        nElectrons_ = nelectrons;
    }

  electronTree_->Fill();

}

DEFINE_FWK_MODULE( ElectronValidator );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
