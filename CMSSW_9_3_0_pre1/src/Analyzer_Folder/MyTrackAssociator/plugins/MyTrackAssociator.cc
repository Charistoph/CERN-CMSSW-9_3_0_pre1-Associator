  // -*- C++ -*-
  //
  // Package:    SimTracker/TrackAssociatorProducers
  // Class:      QuickTrackAssociatorByHitsProducer
  //
  /**\class QuickTrackAssociatorByHitsProducer QuickTrackAssociatorByHitsProducer.cc SimTracker/TrackAssociatorProducers/plugins/QuickTrackAssociatorByHitsProducer.cc

  Description: [one line class summary]

  Implementation:
      [Notes on implementation]
  */
  //
  // Original Author:  Christoph Bernkopf
  //         Created:  Mon, 31 Jul 2017 10:50:34 GMT
  //
  //


  // system include files
  #include <memory>
  #include <iostream>

  // user include files
  #include "FWCore/Framework/interface/Frameworkfwd.h"
  #include "FWCore/Framework/interface/one/EDAnalyzer.h"
  #include "FWCore/Framework/interface/ESHandle.h"

  #include "FWCore/Framework/interface/Event.h"
  #include "FWCore/Framework/interface/MakerMacros.h"

  #include "FWCore/ParameterSet/interface/ParameterSet.h"
  #include "FWCore/Utilities/interface/EDGetToken.h"

  //--------------------------------------------------------
  #include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
  #include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
  // wenn lokal gespeichert
  #include "Analyzer_Folder/MyTrackAssociator/interface/QuickTrackAssociatorByHitsImpl.h"
  // wenn von CMSSW genommen
  //#include "SimTracker/TrackAssociatorProducers/plugins/QuickTrackAssociatorByHitsImpl.h"

  #include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
  #include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

  // Import missing track data
  #include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
  #include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
  #include "DataFormats/EgammaReco/interface/ElectronSeed.h"
  #include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"

  // TTree include
  #include "TTree.h"
  #include "FWCore/ServiceRegistry/interface/Service.h"
  #include "CommonTools/UtilAlgos/interface/TFileService.h"

// GsfTrackToVtx imports
  #include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
  #include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

  #include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
  #include "MagneticField/Engine/interface/MagneticField.h"

  #include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
  #include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"

// Brauch ich die?
  #include "FWCore/MessageLogger/interface/MessageLogger.h"

  //
  // class declaration
  //
  namespace {
  }

  class MyTrackAssociator : public edm::one::EDAnalyzer<edm::one::SharedResources> {
    public:
      explicit MyTrackAssociator(const edm::ParameterSet&);
      ~MyTrackAssociator();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
      virtual void beginJob() override;
  //    virtual void analyze(edm::StreamID, edm::Event&, const edm::EventSetup&) override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      edm::ParameterSet makeHitAssociatorParameters(const edm::ParameterSet&);

  // ----------associator member data ---------------------------
      TrackerHitAssociator::Config trackerHitAssociatorConfig_;
      edm::EDGetTokenT<ClusterTPAssociation> cluster2TPToken_;
      double qualitySimToReco_;
      double puritySimToReco_;
      double pixelHitWeight_;
      double cutRecoToSim_;
      QuickTrackAssociatorByHitsImpl::SimToRecoDenomType simToRecoDenominator_;
      bool threeHitTracksAreSpecial_;
      bool useClusterTPAssociation_;
      bool absoluteNumberOfHits_;

  // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<TrajectorySeed> > TrajectorySeedToken_;
      edm::EDGetTokenT<TrackingParticleCollection> tpToken_;
      edm::EDGetTokenT<edm::View<reco::GsfTrack> > GsfTrackCollectionToken_;
      edm::EDGetTokenT<edm::View<reco::Track> > TrackCollectionToken_;

  // GsfTrackToVtx imports
      edm::ESHandle<GlobalTrackingGeometry> trackingGeometryHandle_;
      edm::ESHandle<MagneticField> magneticFieldHandle_;

  // ----------counting variables ---------------------------
      int indexEvent;
      int assocseedfound;
      int assoctrackfound;
      double seedsuccessrate;
      double tracksuccessrate;

  // ----------TTree Varibs ---------------------------
      TTree * track_tree;
      int track_varib_nr;
      int ml_varib_nr;
      float gsf_track[9];
      float seed_assoc_track[9];
      float track_assoc_track[9];
      float ml_track[28];
  };

  //
  // constants, enums and typedefs
  //

  //
  // static data member definitions
  //

  //
  // constructors and destructor
  //

  MyTrackAssociator::MyTrackAssociator(const edm::ParameterSet& iConfig):
    trackerHitAssociatorConfig_(makeHitAssociatorParameters(iConfig), consumesCollector()),
    qualitySimToReco_( iConfig.getParameter<double>( "Quality_SimToReco" ) ),
    puritySimToReco_( iConfig.getParameter<double>( "Purity_SimToReco" ) ),
    pixelHitWeight_( iConfig.getParameter<double>( "PixelHitWeight" ) ),
    cutRecoToSim_( iConfig.getParameter<double>( "Cut_RecoToSim" ) ),
    threeHitTracksAreSpecial_( iConfig.getParameter<bool>( "ThreeHitTracksAreSpecial" ) ),
    useClusterTPAssociation_( iConfig.getParameter<bool>( "useClusterTPAssociation" ) ),
    absoluteNumberOfHits_( iConfig.getParameter<bool>( "AbsoluteNumberOfHits" ) ){

      indexEvent = 0;
      assocseedfound = 0;
      assoctrackfound = 0;
      seedsuccessrate = 0;
      tracksuccessrate = 0;
      track_varib_nr = 9;
      ml_varib_nr = 28;

      TrajectorySeedToken_ = consumes<edm::View<TrajectorySeed> >(edm::InputTag("electronMergedSeeds"));
      tpToken_ = consumes<TrackingParticleCollection>(edm::InputTag("tpSelection"));
      GsfTrackCollectionToken_ = consumes<edm::View<reco::GsfTrack> >(edm::InputTag("electronGsfTracks"));
      TrackCollectionToken_ = consumes<edm::View<reco::Track> >(edm::InputTag("electronGsfTracks"));

      usesResource("TFileService");

  }

  MyTrackAssociator::~MyTrackAssociator(){

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

  }

  //
  // member functions
  //

  // Set up the parameter set for the hit associator
  edm::ParameterSet
  MyTrackAssociator::makeHitAssociatorParameters(const edm::ParameterSet& iConfig) {
   edm::ParameterSet hitAssociatorParameters;
   hitAssociatorParameters.addParameter<bool>( "associatePixel", iConfig.getParameter<bool>("associatePixel") );
   hitAssociatorParameters.addParameter<bool>( "associateStrip", iConfig.getParameter<bool>("associateStrip") );
   // This is the important one, it stops the hit associator searching through the list of sim hits.
   // I only want to use the hit associator methods that work on the hit IDs (i.e. the uint32_t trackId
   // and the EncodedEventId eventId) so I'm not interested in matching that to the PSimHit objects.
   hitAssociatorParameters.addParameter<bool>("associateRecoTracks",true);
   // add these new ones to allow redirection of inputs:
   hitAssociatorParameters.addParameter<edm::InputTag>( "pixelSimLinkSrc", iConfig.getParameter<edm::InputTag>("pixelSimLinkSrc") );
   hitAssociatorParameters.addParameter<edm::InputTag>( "stripSimLinkSrc", iConfig.getParameter<edm::InputTag>("stripSimLinkSrc") );

   return hitAssociatorParameters;
  }

  // ------------ method called for each event  ------------
  void
  MyTrackAssociator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
  //MyTrackAssociator::analyze(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const
  {
    using namespace edm;

// GsfTrackToVtx handle code
    iSetup.get<GlobalTrackingGeometryRecord>().get(trackingGeometryHandle_);
    iSetup.get<IdealMagneticFieldRecord>().get(magneticFieldHandle_);
    MultiTrajectoryStateTransform mtst(&*trackingGeometryHandle_,&*magneticFieldHandle_);

    edm::Handle<edm::View<reco::GsfTrack> > gsfTrackHandle;
    iEvent.getByToken(GsfTrackCollectionToken_, gsfTrackHandle);
// GsfTrackToVtx handle code end

// handles for impl = QuickTrackAssociatorByHitsImpl
    const ClusterTPAssociation *clusterAssoc = nullptr;
    std::unique_ptr<TrackerHitAssociator> trackAssoc;
    trackAssoc = std::make_unique<TrackerHitAssociator>(iEvent, trackerHitAssociatorConfig_);

    edm::Handle<edm::View<TrajectorySeed> > TrajectorySeedHandle;
    iEvent.getByToken(TrajectorySeedToken_, TrajectorySeedHandle);

    edm::Handle<TrackingParticleCollection> tpHandle;
    iEvent.getByToken(tpToken_,tpHandle);

    edm::Handle<edm::View<reco::GsfTrack> > GsfTrackCollectionHandle;
    iEvent.getByToken(GsfTrackCollectionToken_, GsfTrackCollectionHandle);

    edm::Handle<edm::View<reco::Track> > TrackCollectionHandle;
    iEvent.getByToken(TrackCollectionToken_, TrackCollectionHandle);

    std::cout << "\n" << "------------------------------------------" << "\n" << "\n"
              << "--- Output Prints of MyTrackAssociator ---" << "\n" << "\n"
              << "#TrajectorySeeds = " << TrajectorySeedHandle->size() << "\n"
              << "#TrackingParticles = " << tpHandle->size() << "\n"
              << "#RecoTracks = " << TrackCollectionHandle->size() << "\n" << std::endl;

  // Associator Funktion
    auto impl = std::make_unique<QuickTrackAssociatorByHitsImpl>(iEvent.productGetter(),
        std::move(trackAssoc),
        clusterAssoc,
        absoluteNumberOfHits_,
        qualitySimToReco_,
        puritySimToReco_,
        pixelHitWeight_,
        cutRecoToSim_,
        threeHitTracksAreSpecial_,
        simToRecoDenominator_);

    reco::RecoToSimCollectionSeed mySeedToSim = impl->associateRecoToSim(TrajectorySeedHandle,tpHandle);

    reco::RecoToSimCollection myTrackToSim = impl->associateRecoToSim(TrackCollectionHandle,tpHandle);

// Associator loop
    for ( size_t j=0; j< GsfTrackCollectionHandle->size() ; ++j ) {
        const reco::GsfTrack& gsfTrack = GsfTrackCollectionHandle->at(j);

        ++indexEvent;

        const edm::RefToBase<TrajectorySeed>& mySeedRef = gsfTrack.seedRef();
        reco::RecoToSimCollectionSeed::const_iterator iassocseed = mySeedToSim.find(mySeedRef);

        edm::RefToBase<reco::Track> seedRef(TrackCollectionHandle,j);
        reco::RecoToSimCollection::const_iterator iassoctrack = myTrackToSim.find(seedRef);

        std::cout << "GsfTrackCollectionHandle->size() = " << GsfTrackCollectionHandle->size() << std::endl;
        std::cout << "TrajectorySeedHandle->size() = " << TrajectorySeedHandle->size() << std::endl;
        std::cout << "TrackCollectionHandle->size() = " << TrackCollectionHandle->size() << std::endl;

        for (int k = 0; k < track_varib_nr; k++){
            gsf_track[k] = 0;
            seed_assoc_track[k] = 0;
            track_assoc_track[k] = 0;
        }

        for (int k = 0; k < ml_varib_nr; k++){
            ml_track[k] = 0;
        }

        std::cout << "all track set to 0 worked! Loop Nr = " << j << std::endl;

        gsf_track[0] = gsfTrack.pt();
        gsf_track[1] = gsfTrack.phi();
        gsf_track[2] = gsfTrack.eta();
        gsf_track[3] = gsfTrack.charge();
        gsf_track[4] = gsfTrack.dxy();
        gsf_track[5] = gsfTrack.dz();
        gsf_track[6] = gsfTrack.numberOfValidHits();

        std::cout << "gsf_track fill worked!" << std::endl;

        if (iassocseed != mySeedToSim.end()){
//          std::cout << "\n" << "if (iassocseed != mySeedToSim.end()){" << std::endl;

            std::cout << "Sim to reco seed found!" << std::endl;
            std::cout << "SIZE (*iassocseed).val.size() = " << (*iassocseed).val.size() << std::endl;

            size_t kmax = 0;
            double qmax = -1.;

            for (size_t i = 0; i < (*iassocseed).val.size(); i++) {

                std::cout << "loop entered" << std::endl;
                std::cout << "i = " << i << ", (*iassocseed).val[i].second = " << (*iassocseed).val[i].second << std::endl;

                if ((*iassocseed).val[i].second > qmax){
//                    std::cout << "qmax if entered" << std::endl;
                    kmax = i;
                    qmax = (*iassocseed).val[i].second;
                    std::cout << "qmax = " << qmax << std::endl;
                }
                else {
                    std::cout << "assocseed not filled, bad quality!" << std::endl;
                }
            }
//            std::cout << "qmax #2 = " << qmax << std::endl;

            if ( qmax>0.) {
                std::cout << "qmax < 0 found!" << std::endl;

                gsf_track[7] = float(qmax);

                const edm::Ref<TrackingParticleCollection> tref_seed = (*iassocseed).val[kmax].first;
                seed_assoc_track[0] = tref_seed->pt();
                seed_assoc_track[1] = tref_seed->phi();
                seed_assoc_track[2] = tref_seed->eta();
                seed_assoc_track[3] = tref_seed->charge();
                seed_assoc_track[6] = tref_seed->numberOfTrackerLayers();
                std::cout << "seed_assoc_track writen!" << std::endl;

//                std::cout << "(*iassocseed).val[j].second = " << (*iassocseed).val[j].second << "\n"
//                << "tref_seed->pt() = " << tref_seed->pt() << "\n"
//                << "tref_seed->phi() = " << tref_seed->phi() << "\n"
//                << "tref_seed->eta() = " << tref_seed->eta() << "\n"
//                << "tref_seed->charge() = " << tref_seed->charge() << "\n"
//                << "tref_seed->x^2+y^2 = " << tref_seed->vx()*tref_seed->vx()+tref_seed->vy()*tref_seed->vy() << "\n"
//                << "tref_seed->numberOfTrackerLayers() = " << tref_seed->numberOfTrackerLayers()
//                << std::endl;

                ++assocseedfound;
                std::cout << "assocseedfound # increased!" << "\n" << std::endl;

            }
            else {
                std::cout << "No sim to reco seed!" << "\n" << std::endl;
                gsf_track[7] = -1;
            }
        }

        if (iassoctrack != myTrackToSim.end()){

            std::cout << "Sim to reco track found!" << std::endl;
            std::cout << "#myTrackToSim Size = " << myTrackToSim.size() << std::endl;
//            std::cout << "std::typeid((*iassoctrack).first).name() = " << typeid(*iassoctrack).name() << "\n" << std::endl;
  //          std::cout << "#iassoctrack qual = " << (*iassoctrack).val[i].second << std::endl;
            std::cout << "(*iassoctrack).val.size() = " << (*iassoctrack).val.size() << std::endl;

            size_t kmax = 0;
            double qmax = -1.;

            for (size_t i = 0; i < (*iassoctrack).val.size(); i++) {

                if ((*iassoctrack).val[i].second > qmax){
                    kmax = i;
                    qmax = (*iassoctrack).val[i].second;
                }
                else {
                    std::cout << "iassoctrack not filled, bad quality!" << std::endl;
                }
            }

            if ( qmax>0.) {
                std::cout << "qmax < 0 found!" << std::endl;

                gsf_track[8] = float(qmax);

                const edm::Ref<TrackingParticleCollection> tref_track = (*iassoctrack).val[kmax].first;
                track_assoc_track[0] = tref_track->pt();
                track_assoc_track[1] = tref_track->phi();
                track_assoc_track[2] = tref_track->eta();
                track_assoc_track[3] = tref_track->charge();
                track_assoc_track[6] = tref_track->numberOfTrackerLayers();
                std::cout << "track_assoc_track writen!" << std::endl;

//                std::cout
//                << "tref_track->pt() = " << tref_track->pt() << "\n"
//                << "tref_track->px() = " << tref_track->px() << "\n"
//                << "tref_track->py() = " << tref_track->py() << "\n"
//                << "tref_track->pz() = " << tref_track->pz() << "\n"
//                << "tref_track->vertex() = " << tref_track->vertex() << "\n"
//                << "tref_track->vx() = " << tref_track->vx() << "\n"
//                << "tref_track->vy() = " << tref_track->vy() << "\n"
//                << "tref_track->vz() = " << tref_track->vz() << "\n"
//                << "tref_track->phi() = " << tref_track->phi() << "\n"
//                << "tref_track->eta() = " << tref_track->eta() << "\n"
//                << "tref_track->charge() = " << tref_track->charge() << "\n"
//                << "tref_track->numberOfTrackerLayers() = " << tref_track->numberOfTrackerLayers() << "\n"
//                << std::endl;

//                ml_track[0] = tref_track->pt();
//                ml_track[0] = tref_track->px();
//                ml_track[0] = tref_track->py();
//                ml_track[0] = tref_track->pz();
//                ml_track[0] = tref_track->px();
//                ml_track[0] = tref_track->py();
//                ml_track[0] = tref_track->pz();
//                ml_track[1] = tref_track->phi();
//                ml_track[2] = tref_track->eta();
//                ml_track[3] = tref_track->charge();
//                ml_track[6] = tref_track->numberOfTrackerLayers();

                // GsfTrackToVtx loop code
                for ( edm::View<reco::GsfTrack>::const_iterator igsf=gsfTrackHandle->begin();
                igsf!=gsfTrackHandle->end(); ++igsf ) {

                    if ((tref_track->vx()*tref_track->vx()+tref_track->vy()*tref_track->vy())<0.0625) {

                        TrajectoryStateOnSurface innTSOS = mtst.innerStateOnSurface(*igsf);
                        // Changed from GlobalPoint(0.,0.,0.)
                        TrajectoryStateOnSurface vtxTSOS = mtst.extrapolatedState(innTSOS,GlobalPoint(tref_track->vx(),tref_track->vy(),tref_track->vz()));
                        std::cout << "GsfTrackToVtx code" << "\n"
                        << "innTSOS comp.size = " << innTSOS.components().size() << "\n"
                        << "vtxTSOS comp.size = " << vtxTSOS.components().size() << "\n"
//                        << "vtxTSOS localParameters().vector() = " << vtxTSOS.localParameters().vector() << "\n"
//                        << "vtxTSOS localError().matrix() = " << vtxTSOS.localError().matrix() << "\n"
                        << std::endl;

//                        std::cout
//                        << "Gsf pt = " << igsf->pt() << "\n"
//                        << "Gsf px = " << igsf->px() << "\n"
//                        << "Gsf py = " << igsf->py() << "\n"
//                        << "Gsf pz = " << igsf->pz() << "\n"
//                        << std::endl;
//
//                        std::cout
//                        << "tref_track->pt() = " << tref_track->pt() << "\n"
//                        << "tref_track->px() = " << tref_track->px() << "\n"
//                        << "tref_track->py() = " << tref_track->py() << "\n"
//                        << "tref_track->pz() = " << tref_track->pz() << "\n"
//                        << std::endl;

                        std::cout << "Gsf component loop starts:" << "\n" << std::endl;

                        for (size_t ic=0; ic<vtxTSOS.components().size(); ++ic ) {
//                              std::cout << vtxTSOS.components()[ic].localParameters().vector() << std::endl;

                              LocalPoint assocp(vtxTSOS.surface().toLocal(GlobalPoint(tref_track->px(),tref_track->py(),tref_track->pz())));
//                              std::cout << "toLocal seed g p -> l p = " << assocp.x() << " " << assocp.y() << " " << assocp.z() << std::endl;

                              LocalPoint assocv(vtxTSOS.surface().toLocal(GlobalPoint(tref_track->vx(),tref_track->vy(),tref_track->vz())));
//                              std::cout << "toLocal seed g v -> l v = " << assocv.x() << " " << assocv.y() << " " << assocv.z() << std::endl;

                              std::cout << "\n"
                              << "nc = vtxTSOS.components().size() = " << vtxTSOS.components().size() << "\n"
                              << "ic = " << ic << "\n"
                              << "vtxTSOS.components()[ic].weight() = " << vtxTSOS.components()[ic].weight() << "\n"

                              << "q_sim/p_sim, px_local_sim/pz_local_sim, py_local_sim/pz_local_sim, vx_local_sim, vy_local_sim" << "\n"
                              << "vtxTSOS localParameters().vector() = " << vtxTSOS.components()[ic].localParameters().vector() << "\n"
                              << "vtxTSOS localParameters().matrix() = " << vtxTSOS.components()[ic].localError().matrix() << "\n"
                              << "sim                                  "
                              << tref_track->charge()/tref_track->p() << ", " << assocp.x()/assocp.z() << ", " << -assocp.y()/assocp.z() << ", "
                              << assocv.x() << ", " << assocv.y() << "\n"
                              << "tref_track->charge()/tref_track->p() assocp.x()/assocp.z() -assocp.y()/assocp.z() assocv.x() assocv.y()" << "\n"
                              << std::endl;
                        }

//                        LocalPoint lp(vtxTSOS.surface().toLocal(GlobalPoint(0.,0.,0.)));
//                        std::cout << "   " << lp.x() << " " << lp.y() << " " << lp.z() << std::endl;
//                        GlobalVector gy(vtxTSOS.surface().toGlobal(LocalVector(0.,1.,0.)));
//                        std::cout << "   " << gy.x() << " " << gy.y() << " " << gy.z() << std::endl;
                    }
                    else {
                        std::cout << "Gsf Vertex too far from origin, no Gaussian Mix produced" << "\n" << std::endl;
                    }

                    ++assoctrackfound;
                    std::cout << "assoctrackfound # increased!" << "\n" << std::endl;
                }
            }
        }
        else {
            std::cout << "No sim to reco track!" << "\n" << std::endl;
            gsf_track[8] = -1;
        }

  // to have quality information on both branches
        seed_assoc_track[7] = gsf_track[7];
        seed_assoc_track[8] = gsf_track[8];
        track_assoc_track[7] = gsf_track[7];
        track_assoc_track[8] = gsf_track[8];

        track_tree->Fill();
        std::cout << "all fills worked!" << "\n" << std::endl;
    }

    seedsuccessrate = float(assocseedfound) / float(indexEvent);
    tracksuccessrate = float(assoctrackfound) / float(indexEvent);

    std::cout << "indexEvent = " << indexEvent <<"\n"
    << "assocseedfound = " <<assocseedfound << "\n"
    << "p found = " << seedsuccessrate << "\n"
    << "assoctrackfound = " <<assoctrackfound << "\n"
    << "p found = " << tracksuccessrate << "\n"
    <<"\n" << "------------------------------------------" << "\n" << std::endl;

  }

  //------------------------------------------------------------------------------

  // ------------ method called once each job just before starting event loop  ------------
  void
  MyTrackAssociator::beginJob()
  {

    using namespace edm;

  // initialize tree
    edm::Service<TFileService> fs;
    track_tree = fs->make<TTree>("track_associator_tree","Associator tree with branches" );
    track_tree->Branch("gsf_track", &gsf_track, "gsf_track[9]/F");
    track_tree->Branch("seed_assoc_track", &seed_assoc_track, "seed_assoc_track[9]/F");
    track_tree->Branch("track_assoc_track", &track_assoc_track, "track_assoc_track[9]/F");
    track_tree->Branch("ml_track", &ml_track, "ml_track[28]/F");
  }

  // ------------ method called once each job just after ending the event loop  ------------
  void
  MyTrackAssociator::endJob() {
  }

  // ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
  void
  MyTrackAssociator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
   //The following says we do not know what parameters are allowed so do no validation
   // Please change this to state exactly what you do use, even if it is no parameters
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
  }

  //define this as a plug-in
  DEFINE_FWK_MODULE(MyTrackAssociator);

