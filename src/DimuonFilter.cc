//
// author: thomas.mccauley@cern.ch
//

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include <iostream>
#include <string>
#include <fstream>

class DimuonFilter : public edm::EDFilter 
{
public: 
  explicit DimuonFilter(const edm::ParameterSet&);
  ~DimuonFilter();

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
 
  edm::InputTag muonInputTag_;

  double minMuonPt_;
  double maxMuonEta_;
  double invariantMassMin_;
  double invariantMassMax_;

  std::ofstream csvOut_;
  std::string csvFileName_;
};

DimuonFilter::DimuonFilter(const edm::ParameterSet& iConfig)
  : muonInputTag_(iConfig.getParameter<edm::InputTag>("muonInputTag")),
    minMuonPt_(iConfig.getParameter<double>("minMuonPt")),
    maxMuonEta_(iConfig.getParameter<double>("maxMuonEta")),
    invariantMassMin_(iConfig.getParameter<double>("invariantMassMin")),
    invariantMassMax_(iConfig.getParameter<double>("invariantMassMax")),
    csvFileName_(iConfig.getParameter<std::string>("csvFileName"))
{
  csvOut_.open(csvFileName_.c_str());
}

DimuonFilter::~DimuonFilter()
{}

bool
DimuonFilter::filter(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<reco::MuonCollection> muons;
  event.getByLabel(muonInputTag_, muons);

  if ( ! muons.isValid() )
  {
    std::cerr<<"DimuonFilter: invalid collection"<<std::endl;
    return false;
  }
  
  int charge;
  double pt = 0.0;
  double eta = 0.0;
  double phi;
  double energy, px, py, pz;

  double last_energy = 0.0;
  double last_px = 0.0;
  double last_py = 0.0;
  double last_pz = 0.0;

  double last_pt = 0.0;
  double last_eta = 0.0;
  double last_phi = 0.0;

  // Only examine if there are precisely 2 muons in the event (for simplicity)
  if ( muons->size() != 2 )
    return false;

  int last_charge = 0;
  int combined_charge = 0;

  bool last_was_tracker = false;
  bool last_was_global = false;
  bool this_is_tracker = false;
 
  for ( reco::MuonCollection::const_iterator it = muons->begin(), end = muons->end(); 
        it != end; ++it) 
  {
    // We are only looking for tracker or global muons, not stand-alone
    if ( ! ((*it).isTrackerMuon() || (*it).isGlobalMuon()) )
      return false;

    if  ( (*it).isGlobalMuon() ) // Global muon                                                             
    {
      pt = (*it).globalTrack()->pt();
      phi = (*it).globalTrack()->phi();
      eta = (*it).globalTrack()->eta();
    }
    
    if ( (*it).isTrackerMuon() ) // Tracker muon                                                                      
    {
      pt = (*it).track()->pt();
      phi = (*it).track()->phi();
      eta = (*it).track()->eta();
    }
  
    if ( pt < minMuonPt_ )
      return false;
    if ( fabs(eta) > maxMuonEta_ )
      return false;

    charge = (*it).charge();
    
    if ( last_charge == 0 ) // i.e. this is the first of the pair of muons
    {
      last_charge = charge;
      last_energy = (*it).energy();
      last_px = (*it).px();
      last_py = (*it).py();
      last_pz = (*it).pz();

      last_pt = pt;
      last_eta = eta;
      last_phi = phi;

      if ( (*it).isTrackerMuon() )
      {
        last_was_global = false;
        last_was_tracker = true;
      }
      
      if ( (*it).isGlobalMuon() )
      {
        last_was_tracker = false;
        last_was_global = true;
      }
      
    }
    
    else // we are on the second muon of the pair and can compare charge and calculate invariant mass
    {
      combined_charge = last_charge*charge;
  
      if ( combined_charge < 0 )
        return false;
    
      energy = (*it).energy();
      px = (*it).px();
      py = (*it).py();
      pz = (*it).pz();

      double E  = last_energy + energy;
      double PX = last_px + px;
      double PY = last_py + py;
      double PZ = last_pz + pz;

      double m = E*E;
      m -= (PX*PX + PY*PY + PZ*PZ);
      m = sqrt(m);

      if ( m < invariantMassMin_ || m > invariantMassMax_ )
        return false;

      std::string this_type;

      if ( (*it).isTrackerMuon() )
      {
        this_is_tracker = true;
        this_type = "T";
      }
      
      if ( (*it).isGlobalMuon() )
      {
        this_is_tracker = false; // a global actually contains both, so make it false
        this_type = "G";
      }
      
      std::string last_type;

      if ( last_was_tracker )
        last_type = "T";
      if ( last_was_global )
        last_type = "G";

      // If both are tracker then reject
      if ( last_was_tracker && this_is_tracker )
        return false; // TT
      
      csvOut_<< event.id().run() <<","<< event.id().event() <<","
	     << last_type <<","
             << last_energy <<","<< last_px <<","<< last_py <<","<< last_pz <<","
             << last_pt <<","<< last_eta <<","<< last_phi <<","<< last_charge <<","
             << this_type <<","
	     << energy <<","<< px <<","<< py <<","<< pz <<","
             << pt <<","<< eta <<","<< phi <<","<< charge <<","
             << m <<std::endl;
    }
  }

  return true;
}

void 
DimuonFilter::beginJob()
{
  csvOut_<<"Run,Event,Type1,E1,px1,py1,pz1,pt1,eta1,phi1,Q1,Type2,E2,px2,py2,pz2,pt2,eta2,phi2,Q2,M"<<std::endl;
}

void 
DimuonFilter::endJob() 
{
  csvOut_.close();
}

DEFINE_FWK_MODULE(DimuonFilter);
