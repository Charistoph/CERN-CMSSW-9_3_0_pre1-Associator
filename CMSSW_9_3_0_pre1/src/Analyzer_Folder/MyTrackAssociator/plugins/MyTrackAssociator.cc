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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

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

// ----------counting variables ---------------------------
    int indexEvent;
    int assocseedfound;
    double successrate;

// ----------TTree Varibs ---------------------------
    TTree * track_tree;
    int track_varib_nr;
    float gsf_track[8];
// sts = seed to sim
    float assoc_track[8];

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
    successrate = 0;
    track_varib_nr = 7;

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

// handles for impl = QuickTrackAssociatorByHitsImpl
  const ClusterTPAssociation *clusterAssoc = nullptr;
  std::unique_ptr<TrackerHitAssociator> trackAssoc;
//  if(useClusterTPAssociation_)  {
//    edm::Handle<ClusterTPAssociation> clusterAssocHandle;
//    iEvent.getByToken(cluster2TPToken_,clusterAssocHandle);
//    clusterAssoc = clusterAssocHandle.product();
//  }
//  else {
    // If control got this far then either useClusterTPAssociation_ was false or getting the cluster
    // to TrackingParticle association from the event failed. Either way I need to create a hit associator.
    trackAssoc = std::make_unique<TrackerHitAssociator>(iEvent, trackerHitAssociatorConfig_);
//  }

  edm::Handle<edm::View<TrajectorySeed> > TrajectorySeedHandle;
  iEvent.getByToken(TrajectorySeedToken_, TrajectorySeedHandle);

  edm::Handle<TrackingParticleCollection> tpHandle;
  iEvent.getByToken(tpToken_,tpHandle);

  edm::Handle<edm::View<reco::GsfTrack> > GsfTrackCollectionHandle;
  iEvent.getByToken(GsfTrackCollectionToken_, GsfTrackCollectionHandle);

  edm::Handle<edm::View<reco::Track> > TrackCollectionHandle;
  iEvent.getByToken(TrackCollectionToken_, TrackCollectionHandle);

  std::cout << " " << "\n"
            << "--- Output Prints of MyTrackAssociator ---" << "\n" << "\n"
            << "#TrajectorySeeds = " << TrajectorySeedHandle->size() << "\n"
            << "#TrackingParticles = " << tpHandle->size() << "\n" << "\n" << std::endl;

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

// Testing Seed
// REF - Smart Pointer
//  edm::RefToBase<TrajectorySeed> seedRef(TrajectorySeedHandle,0);
//  reco::RecoToSimCollectionSeed::const_iterator iassocseed = mySeedToSim.find(seedRef);

// Test Prints
//  std::cout << "------------------------------" << "\n" << "\n" << "#mySeedToSim Size = " << mySeedToSim.size() << "\n"
//            << "std::typeid((*iassocseed).first).name() = " << typeid(*iassocseed).name() << "\n"
//            << "#iassocseed size = " << (*iassocseed).val.size() << "\n"
//            << "\n" << "------------------------------" << "\n" << "\n"
//            << "iassocseed Loop " << std::endl;

//  std::cout << "\n" << "pt, phi, eta, charge, vertex, pdgId, #TRLayers, qual" << " " << std::endl;

  for ( size_t j=0; j< GsfTrackCollectionHandle->size() ; ++j ) {
      const reco::GsfTrack& gsfTrack = GsfTrackCollectionHandle->at(j);

      ++indexEvent;

      const edm::RefToBase<TrajectorySeed>& mySeedRef = gsfTrack.seedRef();
      reco::RecoToSimCollectionSeed::const_iterator iassocseed = mySeedToSim.find(mySeedRef);

      edm::RefToBase<reco::Track> seedRef(TrackCollectionHandle,j);
      reco::RecoToSimCollection::const_iterator iassoctrack = myTrackToSim.find(seedRef);

      // Wert der Varib in Tree datenstruktur kopieren
      //          for (int k = 0; track_varib_nr; ++k){
      //            std::cout << "Set to 0 loop! " << k << std::endl;
      //            assoc_track[k] = 0;
      //          }
      assoc_track[0] = 0;
      assoc_track[1] = 0;
      assoc_track[2] = 0;
      assoc_track[3] = 0;
      assoc_track[4] = 0;
      assoc_track[5] = 0;
      assoc_track[6] = 0;
      assoc_track[7] = 0;

      std::cout << "assoc_track set to 0 worked! " << std::endl;

      assoc_track[0] = gsfTrack.pt();
      assoc_track[1] = gsfTrack.phi();
      assoc_track[2] = gsfTrack.eta();
      assoc_track[3] = gsfTrack.charge();
      assoc_track[4] = gsfTrack.dxy();
      assoc_track[5] = gsfTrack.dz();
      assoc_track[6] = gsfTrack.numberOfValidHits();

      if (iassocseed != mySeedToSim.end()){
//        for ( size_t i=0; i< (*iassocseed).val.size(); ++i ) {

          std::cout << "Sim to reco found" << "\n" << std::endl;

          if ((*iassocseed).val[j].second == 1){
            assoc_track[7] = float((*iassocseed).val[j].second);
            track_tree->Fill();
            std::cout << "assoc filled! " << std::endl;
          }

          else {
            std::cout << "not filled, bad quality! " << std::endl;
          }

          ++assocseedfound;

//        }
      }

      if (iassoctrack != myTrackToSim.end()){

        std::cout << "#myTrackToSim Size = " << myTrackToSim.size() << "\n"
        << "std::typeid((*iassoctrack).first).name() = " << typeid(*iassoctrack).name() << "\n" << std::endl;
        std::cout << "#iassoctrack qual = " << (*iassoctrack).val[j].second << std::endl;
        std::cout << "iassoctrack qual works." << "\n" << std::endl;

      }

      else {

        std::cout << "No sim to reco! pt = " << gsfTrack.pt() << "\n" << std::endl;
        assoc_track[7] = -1;
        track_tree->Fill();
        std::cout << "gsf Fill worked! " << std::endl;

      }
  }

  successrate = float(assocseedfound) / float(indexEvent);

  std::cout << "indexEvent = " << indexEvent <<"\n"
  << "assocseedfound = " <<assocseedfound << "\n"
  << "p found = " << successrate << "\n"
  <<"\n" << "------------------------------" << "\n" << std::endl;

}

//------------------------------------------------------------------------------

// ------------ method called once each job just before starting event loop  ------------
void
MyTrackAssociator::beginJob()
{

  using namespace edm;

// initialize tree
  edm::Service<TFileService> fs;
  track_tree = fs->make<TTree>("track_associator_tree","Associator tree with one branch" );
  track_tree->Branch("assoc_track", &assoc_track, "assoc_track[8]/F");

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

