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
//#include "DataFormats/TrackReco/interface/TrackBase.h"

// REF
//#include "DataFormats/Common/interface/Ref.h"
//#include "DataFormats/Common/interface/RefVector.h"
//#include "DataFormats/Common/interface/RefToBase.h"
//#include "DataFormats/Common/interface/RefToBaseVector.h"

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
    int indexEvent;

    edm::EDGetTokenT<edm::View<TrajectorySeed> > TrajectorySeedToken_;

    edm::EDGetTokenT<TrackingParticleCollection> tpToken_;

    edm::EDGetTokenT<edm::View<reco::GsfTrack> > GsfTrackCollectionToken_;

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

    TrajectorySeedToken_ = consumes<edm::View<TrajectorySeed> >(edm::InputTag("electronMergedSeeds"));

    tpToken_ = consumes<TrackingParticleCollection>(edm::InputTag("tpSelection"));

    GsfTrackCollectionToken_ = consumes<edm::View<reco::GsfTrack> >(edm::InputTag("electronGsfTracks"));

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

  std::cout << " " << "\n"
            << "--- Output Prints of MyTrackAssociator ---" << "\n" << "\n"
            << "#TrackingParticles = " << tpHandle->size() << "\n"
            << "#TrajectorySeeds = " << TrajectorySeedHandle->size() << "\n" << "\n" << std::endl;


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

// Testing Seed
// REF - Smart Pointer
//  edm::RefToBase<TrajectorySeed> seedRef(TrajectorySeedHandle,0);
//  reco::RecoToSimCollectionSeed::const_iterator iassoc = mySeedToSim.find(seedRef);

// Test Prints
//  std::cout << "------------------------------" << "\n" << "\n" << "#mySeedToSim Size = " << mySeedToSim.size() << "\n"
//            << "std::typeid((*iassoc).first).name() = " << typeid(*iassoc).name() << "\n"
//            << "#iassoc size = " << (*iassoc).val.size() << "\n"
//            << "\n" << "------------------------------" << "\n" << "\n"
//            << "iassoc Loop " << std::endl;

std::cout << "\n" << "pt, phi, eta, charge, vertex, pdgId, #TRLayers, qual" << " " << std::endl;

for ( size_t j=0; j< GsfTrackCollectionHandle->size() ; ++j ) {
    const reco::GsfTrack& gsfTrack = GsfTrackCollectionHandle->at(j);

    const edm::RefToBase<TrajectorySeed>& mySeedRef = gsfTrack.seedRef();
    reco::RecoToSimCollectionSeed::const_iterator iassoc = mySeedToSim.find(mySeedRef);

    std::cout << "GSF Track " << j << "\n"<< std::endl;

    if (iassoc != mySeedToSim.end()){
      for ( size_t i=0; i< (*iassoc).val.size(); ++i ) {
        
        std::cout << "iassoc test" << i << "\n"<< std::endl;

        const edm::Ref<TrackingParticleCollection> tref = (*iassoc).val[i].first;
        double qual = (*iassoc).val[i].second;

        std::cout << "\n" << "Event: " << i << " "
        << tref->pt() << " "
        << tref->phi() << " "
        << tref->eta() << " "
        << tref->charge() << " "
        << tref->vertex() << " "
        << tref->pdgId() << " "
        // The method matchedHit() has been deprecated. Use numberOfTrackerLayers() instead.
        << tref->numberOfTrackerLayers() << " "
        << qual << "\n" << std::endl;
      }
    }

}

  std::cout << "\n" << "------------------------------" << "\n" << std::endl;

}

//------------------------------------------------------------------------------

// ------------ method called once each job just before starting event loop  ------------
void
MyTrackAssociator::beginJob()
{
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

