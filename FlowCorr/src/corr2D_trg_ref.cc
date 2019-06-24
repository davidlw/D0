// ntracks done
// d0 selections done
// d0 daughter rejection done
// correlation same event done
// correlation mixed events

#include <iostream>
#include <vector>
#include <map>
#include <list>

#include "TFile.h"
#include "TChain.h"
#include "Event.h"
#include "TTree.h"
#include "TFileCollection.h"
#include "TCollection.h"
#include "THashList.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVector3.h"
#include "TString.h"

#include "myAnaConsts.h"

using namespace std;

void setBranchStatus(Event*);
bool checkBranchStatus(Event*);

bool passGoodTrack(Event*, const unsigned int&);
inline bool passGoodVtx(Event* event);
inline bool passD0Selections(Event*, const int&, const int&);
bool passD0PreSelections(Event*, const int&);
bool passD0KinematicCuts(Event*, const int&);
inline bool passD0MVA(Event*, const int&, const int&);

int main(int argc, char** argv)
{

   TH1::SetDefaultSumw2(true);

   if(argc!=2) {
      std::cerr << "The number of arguments is wrong" << std::endl;
      return -1;
   }
   
   string datalist(argv[1]);
   std::cout << datalist << std::endl;

   TFile fout(TString::Format("fout_ref_%s.root", datalist.c_str()), "recreate");

   TChain *chain_d0 = new TChain("npd0ana1/VertexCompositeNtuple");
   TChain *chain_tracks = new TChain("track_ana/trackTree");

   TFileCollection* fcData = new TFileCollection(datalist.c_str(), "", datalist.c_str());

   chain_d0->AddFileInfoList(fcData->GetList());
   std::cout << "d0 ready" << std::endl;

   chain_tracks->AddFileInfoList(fcData->GetList());
   std::cout << "tracks ready" << std::endl;

   Event* evt = new Event(chain_d0, chain_tracks);
   setBranchStatus(evt);
   if(!checkBranchStatus(evt)) return 0;

   // declare hists
   TH1D* hMult;


   TH2D* hSignal_Ref;
   TH2D* hBackground_Ref;

   hMult = new TH1D("hMult", "", 600, 0, 600);

   hSignal_Ref = new TH2D("hSignal_Ref",
            "", ana::nEtaBin, ana::etaBegin, ana::etaEnd,
            ana::nPhiBin, ana::phiBegin, ana::phiEnd);
   hBackground_Ref = new TH2D("hBackground_Ref",
            "", ana::nEtaBin, ana::etaBegin, ana::etaEnd,
            ana::nPhiBin, ana::phiBegin, ana::phiEnd);

   // declare vectors
   vector<TVector3> pVect_d0;

   vector<TVector3> pVect_dau1_d0;

   vector<TVector3> pVect_dau2_d0;

   vector<TVector3> pVect_ref;
   list<vector<TVector3>> pVectList_ref[ana::nZ_Vtx_];
   vector<float> effVect_ref;
   list<vector<float>> effVectList_ref[ana::nZ_Vtx_];

   // start timing
   TStopwatch ts;
   ts.Start();

   // loop
   //
   // temporary vectors
   TVector3 p_dau1(0, 0, 0), p_dau2(0, 0, 0), p_d0(0, 0, 0), p_ref(0, 0, 0), p_temp(0, 0, 0);

   std::cout << evt->GetEntries() << std::endl;
   long int nentries = evt->GetEntries();
   int percent = 0;
   long int skip = 0;
   for(long int ientry=0; ientry<nentries; ientry++){
      if(pVect_ref.size()!=0){
         cout << "did not clear the vector pVect_ref.size" << endl;
         continue;
      }
      if(effVect_ref.size()!=0){
         cout << "did not clear the vector effVect_ref.size" << endl;
         continue;
      }
      if(pVect_d0.size()!=0){
         cout << "did not clear the vector pVect_d0.size" << endl;
         continue;
      }
      if(pVect_dau1_d0.size()!=0){
         cout << "did not clear the vector pVect_dau1_d0.size" << endl;
         continue;
      }
      if(pVect_dau2_d0.size()!=0){
         cout << "did not clear the vector pVect_dau2_d0.size" << endl;
         continue;
      }

      int current = (int)(ientry+1)*100/nentries;
      if( current == 5*percent){
         std::cout << 5*percent++ << "\% percents completed" << std::endl;
      }
      auto bytes = evt->GetEntry(ientry);
      if(bytes == 0 || bytes == -1) {
      // std::cout << "vertex unmatched " << ientry << std::endl;
         skip++;
         continue;
      }

      // good vertex check
      if(!passGoodVtx(evt)) continue;
      auto iz = ana::findZVtxBin(evt->BestVtxZ());
      if(iz == -1) continue;

      // count number of good tracks per event
      unsigned int nMult_trk_good = 0;
      for(unsigned int itrack=0; itrack<evt->CandSizeTrk(); itrack++){
         // assume all tracks in TTree are good, since error of dz and dxy are not available
         if(passGoodTrack(evt, itrack)) 
            nMult_trk_good++;
      }
      hMult->Fill(nMult_trk_good);

      if(nMult_trk_good<ana::multMax_ && nMult_trk_good>=ana::multMin_){
         for(int id0=0; id0<evt->CandSize(); id0++){
            int ipt = ana::findPtBin(evt->Pt(id0));
            int iy = ana::findYBin(evt->Y(id0));

            if(ipt == -1) continue;
            if(iy == -1) continue;
            if(!passD0Selections(evt, id0, ipt)) continue;

            p_dau1.SetPtEtaPhi(evt->PtD1(id0), evt->etaD1(id0), evt->phiD1(id0));
            p_dau2.SetPtEtaPhi(evt->PtD2(id0), evt->etaD2(id0), evt->phiD2(id0));
            p_d0 = p_dau1 + p_dau2;

            pVect_d0.push_back(p_d0);
            pVect_dau1_d0.push_back(p_dau1);
            pVect_dau2_d0.push_back(p_dau2);
         }
      }else{
         //std::cout << "multiplicity wrong" << std::endl;
         continue;
      }

      for(unsigned int itrack=0; itrack<evt->CandSizeTrk(); itrack++){
         // some cuts are not available
         bool passDzErr = true;
         bool passDxyErr = true;
         bool passTrkPurity = true;
         bool passPt = evt->PtTrk(itrack) > ana::ptMin_ass_ && evt->PtTrk(itrack) < ana::ptMax_ass_;
         bool passPtError = true;
         bool passEta = evt->EtaTrk(itrack) > ana::etaMin_ass_ && evt->EtaTrk(itrack) < ana::etaMax_ass_;
         bool pass_rejecDaughter = true;
         for(unsigned int id0=0; id0<pVect_d0.size(); id0++){
            if(ana::rejectDaughter_){
               p_temp.SetPtEtaPhi(evt->PtTrk(itrack), evt->EtaTrk(itrack), evt->PhiTrk(itrack));
               if(fabs(p_temp.Eta() - pVect_dau1_d0.at(id0).Eta()<0.03)
                     && fabs(p_temp.DeltaPhi(pVect_dau1_d0.at(id0)))<0.03)
                  pass_rejecDaughter = false;
               if(fabs(p_temp.Eta() - pVect_dau2_d0.at(id0).Eta()<0.03)
                     && fabs(p_temp.DeltaPhi(pVect_dau2_d0.at(id0)))<0.03)
                  pass_rejecDaughter = false;
            }
         }

         bool pass_ass_ = passDzErr && passDxyErr && passTrkPurity && passPt &&
                           passPtError && passEta && pass_rejecDaughter;
         if(pass_ass_) {
            p_ref.SetPtEtaPhi(evt->PtTrk(itrack), evt->EtaTrk(itrack), evt->PhiTrk(itrack));
            pVect_ref.push_back(p_ref);
            effVect_ref.push_back(evt->WeightTrk(itrack));
         }
      }

      // calculate signal
      unsigned int nMult_ref = (unsigned int) pVect_ref.size();
      double nMult_eff_ref = 0;
      for(unsigned int iref=0; iref<nMult_ref; iref++){
         nMult_eff_ref += 1./effVect_ref.at(iref);
      }

      for(unsigned int itrg=0; itrg<nMult_ref; itrg++){
         //for(unsigned int iass=itrg+1; iass<nMult_ref; iass++){
         for(unsigned int iass=0; iass<nMult_ref; iass++){
            if(iass == itrg) continue;
            TVector3 pvector_trg = pVect_ref.at(itrg);
            double effweight_trg = effVect_ref.at(itrg);

            TVector3 pvector_ass = pVect_ref.at(iass);
            double effweight_ass = effVect_ref.at(iass);

            double deltaEta = pvector_ass.Eta() - pvector_trg.Eta();
            double deltaPhi = pvector_ass.DeltaPhi(pvector_trg);
            if(deltaPhi>-ana::PI && deltaPhi<-ana::PI/2.) deltaPhi += 2*ana::PI;

            hSignal_Ref->Fill(deltaEta, deltaPhi, 
                  1./nMult_eff_ref/effweight_trg/effweight_ass);
         }
      }
      //

      // mixed events
      // begin filling background
      // if having enough, erase the first one, fill the current, fill the background
      if(pVectList_ref[iz].size() == ana::nMixedEvts){
         for(unsigned int itrg=0; itrg<nMult_ref; itrg++){
            // simultaneously read momentum and efficiency
            auto ievt_eff_ass = effVectList_ref[iz].begin();
            for(auto ievt_p_ass= pVectList_ref[iz].begin(); ievt_p_ass != pVectList_ref[iz].end(); ievt_p_ass++){

               TVector3 pvector_trg = pVect_ref.at(itrg);
               double effweight_trg = effVect_ref.at(itrg);

               unsigned int n_ass = ievt_p_ass->size();
               if(n_ass!=ievt_eff_ass->size()){
                  cout << "wrong" << endl;
                  break;
               }
               for(unsigned int iass=0; iass<n_ass; iass++){
                  TVector3 pvector_ass = ievt_p_ass->at(iass);
                  double effweight_ass = ievt_eff_ass->at(iass);
                  double deltaEta = pvector_ass.Eta() - pvector_trg.Eta();
                  double deltaPhi = pvector_ass.DeltaPhi(pvector_trg);
                  if(deltaPhi>-ana::PI && deltaPhi<-ana::PI/2.) deltaPhi += 2*ana::PI;

                  hBackground_Ref->Fill(deltaEta, deltaPhi, 
                     1./nMult_eff_ref/effweight_trg/effweight_ass);
               }
               ievt_eff_ass++;
            }
         }
      }

      // clear all
      if(pVectList_ref[iz].size() == ana::nMixedEvts){
         // update the list/buffer
         pVectList_ref[iz].erase(pVectList_ref[iz].begin());
         pVectList_ref[iz].push_back(pVect_ref);
         effVectList_ref[iz].erase(effVectList_ref[iz].begin());
         effVectList_ref[iz].push_back(effVect_ref);
      } else{
         // if not having enough events, keep filling the list of momentum vector of associated particles
         // fill the list/buffer
         pVectList_ref[iz].push_back(pVect_ref);
         effVectList_ref[iz].push_back(effVect_ref);
      }
      pVect_d0.clear();
      pVect_dau1_d0.clear();
      pVect_dau2_d0.clear();
      pVect_ref.clear();
      effVect_ref.clear();
   }
   std::cout << "completed loop" << std::endl;
   std::cout << skip << " events are skipped" << std::endl;

   ts.Stop();
   ts.Print();

   fout.cd();

   // start writing output
   TH2D* hCorrected_Ref;

   hCorrected_Ref = (TH2D*) hSignal_Ref->Clone();
   hCorrected_Ref->SetName("hCorrected_Ref");


   int etaBin = hBackground_Ref->GetXaxis()->FindBin(0.);
   int phiBin = hBackground_Ref->GetYaxis()->FindBin(0.);
   hCorrected_Ref->Divide(hBackground_Ref);
   hCorrected_Ref->Scale(hBackground_Ref->GetBinContent(etaBin, phiBin));
   
   hMult->Write();
   hSignal_Ref->Write();
   hBackground_Ref->Write();
   hCorrected_Ref->Write();

   delete hMult;
   delete hSignal_Ref;
   delete hBackground_Ref;
   delete hCorrected_Ref;

   delete evt;


   return 0;
}

void setBranchStatus(Event* evt)
{
   evt->SetBranchStatus("candSize", 1);
   evt->SetBranchStatus("pT", 1);
   evt->SetBranchStatus("mass", 1);
   evt->SetBranchStatus("mva", 1);
   evt->SetBranchStatus("y", 1);
   evt->SetBranchStatus("3DPointingAngle", 1);
   evt->SetBranchStatus("3DDecayLength", 1);
   evt->SetBranchStatus("3DDecayLengthSignificance", 1);
   evt->SetBranchStatus("*D1*", 1);
   evt->SetBranchStatus("*D2*", 1);
   evt->SetBranchStatus("*Daugther1", 1); // mistype daughter... the writer of TTree
   evt->SetBranchStatus("*Daugther2", 1);
   evt->SetBranchStatus("dedx*", 0);

   evt->SetBranchStatus("tracks.candSizeTRK", 1);
   evt->SetBranchStatus("tracks.pTTRK", 1);
   evt->SetBranchStatus("tracks.etaTRK", 1);
   evt->SetBranchStatus("tracks.phiTRK", 1);
   evt->SetBranchStatus("tracks.weightTRK", 1);
}

bool checkBranchStatus(Event* event)
{
   bool check =
      event->GetBranchStatus("candSize")&&

      event->GetBranchStatus("pT");
      event->GetBranchStatus("mass")&&
      event->GetBranchStatus("mva")&&
      event->GetBranchStatus("y")&&
      event->GetBranchStatus("3DPointingAngle")&&
      event->GetBranchStatus("3DDecayLength")&&
      event->GetBranchStatus("3DDecayLengthSignificance")&&

      event->GetBranchStatus("pTD1")&&
      event->GetBranchStatus("pTerrD1")&&
      event->GetBranchStatus("EtaD1")&&
      event->GetBranchStatus("PhiD1")&&
      event->GetBranchStatus("zDCASignificanceDaugther1")&&
      event->GetBranchStatus("xyDCASignificanceDaugther1")&&
      event->GetBranchStatus("NHitD1")&&
      event->GetBranchStatus("HighPurityDaugther1")&&

      event->GetBranchStatus("pTD2")&&
      event->GetBranchStatus("pTerrD2")&&
      event->GetBranchStatus("EtaD2")&&
      event->GetBranchStatus("PhiD2")&&
      event->GetBranchStatus("zDCASignificanceDaugther2")&&
      event->GetBranchStatus("xyDCASignificanceDaugther2")&&
      event->GetBranchStatus("NHitD2")&&
      event->GetBranchStatus("HighPurityDaugther2")&&

      event->GetBranchStatus("tracks.bestvtxX")&&
      event->GetBranchStatus("tracks.bestvtxY")&&
      event->GetBranchStatus("tracks.bestvtxZ")&&
      event->GetBranchStatus("tracks.candSizeTRK")&&
      event->GetBranchStatus("tracks.pTTRK")&&
      event->GetBranchStatus("tracks.etaTRK")&&
      event->GetBranchStatus("tracks.phiTRK")&&
      event->GetBranchStatus("tracks.weightTRK");

   return check;
}

inline bool passD0Selections(Event* event, const int& icand, const int& ipt)
{
   if(!passD0PreSelections(event, icand)) return false;
   if(!passD0KinematicCuts(event, icand)) return false;
   if(!passD0MVA(event, icand, ipt)) return false;
   return true;
}

bool passD0PreSelections(Event* event, const int& icand)
{
   bool passPointingAngle = std::fabs(event->PointingAngle3D(icand)) < 1;
   bool passTrkEta = std::fabs(event->etaD1(icand)) < 2.4 && std::fabs(event->etaD2(icand)) < 2.4;
   bool passTrkPt = event->PtD1(icand) > 0.7 && event->PtD2(icand) > 0.7;
   bool passTrkPtErr = event->PtErrD1(icand)/event->PtD1(icand) < 0.1 && event->PtErrD2(icand)/event->PtD2(icand) < 0.1;
   bool passTrkPurity = event->highPurityD1(icand) && event->highPurityD2(icand);
   bool passTrkNhits = event->nHitD1(icand) >=11 && event->nHitD2(icand) >=11; 
   bool passDeltaEta = std::fabs(event->etaD1(icand) - event->etaD2(icand)) < 1;

   if(passPointingAngle && passTrkEta && passDeltaEta
         && passTrkPt && passTrkPtErr && passTrkPurity && passTrkNhits
         ) return true;
   return false;
}

bool passD0KinematicCuts(Event* event, const int& icand)
{
   bool passEta = fabs(event->Eta(icand)<1000.);
   return passEta;
}

inline bool passD0MVA(Event* event, const int& icand, const int& ipt)
{
   return event->Mva(icand) > ana::mvaCut[ipt];
}

inline bool passGoodVtx(Event* event)
{
   if(event->BestVtxZ() < -15. || event->BestVtxZ() > 15.) return false;
   return true;
}

bool passGoodTrack(Event* event, const unsigned int& icand)
{
   bool passHighPurity = true;
   bool passDzErr = true;
   bool passDxyErr = true;
   bool passPt = event->PtTrk(icand) > 0.4;
   bool passPtError = true;
   bool passEta = fabs(event->EtaTrk(icand)) < 2.4;
   return passHighPurity && passDzErr && passDxyErr &&
      passPt && passPtError && passEta;
}
