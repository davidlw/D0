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
//bool passD0KinematicCuts(Event*, const int&);
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

   TFile fout(TString::Format("fout%s.root", datalist.c_str()), "recreate");

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
   TH1D* hMult_ass;

   TH1D* hKET_D0[ana::nPt];
   TH1D* hPt_D0[ana::nPt];
   TH1D* hEta_D0[ana::nPt];
   TH1D* hRapidity_D0[ana::nPt];

   TH1D* hMass_D0[ana::nMass][ana::nPt];

   TH1D* hMult_raw_D0[ana::nMass][ana::nPt];
   TH1D* hMult_eff_D0[ana::nMass][ana::nPt];
   map<string, TH2*> hSignal_D0[ana::nMass][ana::nPt];
   map<string, TH2*> hBackground_D0[ana::nMass][ana::nPt];

   hMult = new TH1D("hMult", "", 600, 0, 600);
   hMult_ass = new TH1D("hMult_ass", "", 600, 0, 600);
   for(int ipt=0; ipt<ana::nPt; ipt++){
      hKET_D0[ipt] = new TH1D(Form("hKET_pt%d", ipt), "", 3000, 0, 30);
      hPt_D0[ipt] = new TH1D(Form("hPt_pt%d", ipt), "", 3000, 0, 30);
      hEta_D0[ipt] = new TH1D(Form("hEta_pt%d", ipt), "", 24, -2.4, 2.4);
      hRapidity_D0[ipt] = new TH1D(Form("hRapidity_pt%d", ipt), "", 24, -2.4, 2.4);
      for(int imass=0; imass<ana::nMass; imass++){
         hMass_D0[imass][ipt] = new TH1D(Form("hMassD0_mass%d_pt%d", imass, ipt),
               "", 200, 1.5, 2.5);
         hMult_raw_D0[imass][ipt] = new TH1D(Form("hMult_raw_D0_mass%d_pt%d", imass, ipt),
               "", 50, 0, 50);
         hMult_eff_D0[imass][ipt] = new TH1D(Form("hMult_eff_D0_mass%d_pt%d", imass, ipt),
               "", 50, 0, 50);
         (hSignal_D0[imass][ipt])["largeDCA"] = new TH2D(Form("hSignal_mass%d_pt%d_largeDCA", imass, ipt),
               "", ana::nEtaBin, ana::etaBegin, ana::etaEnd,
               ana::nPhiBin, ana::phiBegin, ana::phiEnd);
         (hBackground_D0[imass][ipt])["largeDCA"] = new TH2D(Form("hBackground_mass%d_pt%d_largeDCA", imass, ipt),
               "", ana::nEtaBin, ana::etaBegin, ana::etaEnd,
               ana::nPhiBin, ana::phiBegin, ana::phiEnd);
         (hSignal_D0[imass][ipt])["smallDCA"] = new TH2D(Form("hSignal_mass%d_pt%d_smallDCA", imass, ipt),
               "", ana::nEtaBin, ana::etaBegin, ana::etaEnd,
               ana::nPhiBin, ana::phiBegin, ana::phiEnd);
         (hBackground_D0[imass][ipt])["smallDCA"] = new TH2D(Form("hBackground_mass%d_pt%d_smallDCA", imass, ipt),
               "", ana::nEtaBin, ana::etaBegin, ana::etaEnd,
               ana::nPhiBin, ana::phiBegin, ana::phiEnd);
      }
   }

   // declare vectors
   vector<TVector3> pVect_trg_d0[ana::nMass][ana::nPt];

   vector<TVector3> pVect_dau1_d0[ana::nMass][ana::nPt];

   vector<TVector3> pVect_dau2_d0[ana::nMass][ana::nPt];

   vector<TVector3> pVect_ass;
   list<vector<TVector3>> pVectList_ass[ana::nZ_Vtx_];
   vector<float> effVect_ass;
   list<vector<float>> effVectList_ass[ana::nZ_Vtx_];

   // start timing
   TStopwatch ts;
   ts.Start();

   // loop
   //
   // temporary vectors
   TVector3 p_dau1(0, 0, 0), p_dau2(0, 0, 0), p_d0(0, 0, 0), p_ass(0, 0, 0);

   std::cout << evt->GetEntries() << std::endl;
   long int nentries = evt->GetEntries();
   int percent = 0;
   long int skip = 0;
   for(long int ientry=0; ientry<nentries; ientry++){

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
      if(!passGoodVtx(evt)) continue;
      auto iz = ana::findZVtxBin(evt->BestVtxZ());
      if(iz == -1) continue;

      // count number of good tracks per event
      unsigned int nMult_ass_good = 0;
      /*
      for(unsigned int itrack=0; itrack<evt->CandSizeTrk(); itrack++){
         // assume all tracks in TTree are good, since error of dz and dxy are not available
         if(passGoodTrack(evt, itrack)) 
            nMult_ass_good++;
      }
      hMult->Fill(nMult_ass_good);
      */

      nMult_ass_good = evt->nTrkOffline();

      if(nMult_ass_good<ana::multMax_ && nMult_ass_good>=ana::multMin_){
         for(int id0=0; id0<evt->CandSize(); id0++){
            int imass = ana::findMassBin(evt->Mass(id0));
            int ipt = ana::findPtBin(evt->Pt(id0));
            int iy = ana::findYBin(evt->Y(id0));

            if(imass == -1) continue;
            if(ipt == -1) continue;
            if(iy == -1) continue;
            if(!passD0Selections(evt, id0, ipt)) continue;

            p_dau1.SetPtEtaPhi(evt->PtD1(id0), evt->etaD1(id0), evt->phiD1(id0));
            p_dau2.SetPtEtaPhi(evt->PtD2(id0), evt->etaD2(id0), evt->phiD2(id0));
            p_d0 = p_dau1 + p_dau2;

            double effks = 1.0;

            hMass_D0[imass][ipt]->Fill(evt->Mass(id0));
            hPt_D0[ipt]->Fill(evt->Pt(id0), 1./effks);
            hEta_D0[ipt]->Fill(p_d0.Eta(), 1./effks);
            hRapidity_D0[ipt]->Fill(evt->Y(id0), 1./effks);
            double KET = sqrt(pow(evt->Mass(id0), 2) + pow(evt->Pt(id0), 2)
                  - evt->Mass(id0));
            hKET_D0[ipt]->Fill(KET, 1./effks);

            pVect_trg_d0[imass][ipt].push_back(p_d0);
            pVect_dau1_d0[imass][ipt].push_back(p_dau1);
            pVect_dau2_d0[imass][ipt].push_back(p_dau2);
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
         bool pass_ass_ = passDzErr && passDxyErr && passTrkPurity && passPt &&
                           passPtError && passEta;

         if(pass_ass_) {
            p_ass.SetPtEtaPhi(evt->PtTrk(itrack), evt->EtaTrk(itrack), evt->PhiTrk(itrack));
            pVect_ass.push_back(p_ass);
            effVect_ass.push_back(evt->WeightTrk(itrack));
         }
      }

      // calculate signal
      unsigned int nMult_ass = (unsigned int) pVect_ass.size();
      hMult_ass->Fill(nMult_ass);
      if((int)nMult_ass_good!=evt->nTrkOffline()) {
         std::cout << "unmatched nTrackOffline" << std::endl;
         std::cout << nMult_ass_good << std::endl;
         std::cout << evt->nTrkOffline() << std::endl;
      }
      
      unsigned int nMult_trg_raw_d0[ana::nMass][ana::nPt] = {0}; // Ntrig for mass & pt bins
      double nMult_trg_eff_d0[ana::nMass][ana::nPt] = {0.}; // eff corrected Ntrig for mass & pt bins

      for(int imass=0; imass<ana::nMass; imass++){
         for(int ipt=0; ipt<ana::nPt; ipt++){
            unsigned int nMult_trg_d0 = (unsigned int) pVect_trg_d0[imass][ipt].size();
            for(unsigned int id0=0; id0<nMult_trg_d0; id0++){
               if(ipt!=ana::findPtBin(pVect_trg_d0[imass][ipt].at(id0).Pt())) std::cout << "pT bin error" << std::endl;
               double effks = 1.0;
               nMult_trg_raw_d0[imass][ipt] += 1;
               nMult_trg_eff_d0[imass][ipt] += 1./effks;
            }
            hMult_raw_D0[imass][ipt]->Fill(nMult_trg_raw_d0[imass][ipt]);
            hMult_eff_D0[imass][ipt]->Fill(nMult_trg_eff_d0[imass][ipt]);

            for(unsigned int id0=0; id0<nMult_trg_d0; id0++){
               if(ipt!=ana::findPtBin(pVect_trg_d0[imass][ipt].at(id0).Pt())) std::cout << "pT bin error" << std::endl;
               double effks = 1.0;
               for(unsigned int iass=0; iass<nMult_ass; iass++){
                  TVector3 pvector_ass = pVect_ass.at(iass);
                  double effweight_ass = effVect_ass.at(iass);
                  if(ana::rejectDaughter_){
                     if(fabs(pvector_ass.Eta() - pVect_dau1_d0[imass][ipt].at(id0).Eta())<0.03
                           && fabs(pvector_ass.DeltaPhi(pVect_dau1_d0[imass][ipt].at(id0)))<0.03)
                        continue;
                     if(fabs(pvector_ass.Eta() - pVect_dau2_d0[imass][ipt].at(id0).Eta())<0.03
                           && fabs(pvector_ass.DeltaPhi(pVect_dau2_d0[imass][ipt].at(id0)))<0.03)
                        continue;
                  }

                  double deltaEta = pvector_ass.Eta() - pVect_trg_d0[imass][ipt].at(id0).Eta();
                  double deltaPhi = pvector_ass.DeltaPhi(pVect_trg_d0[imass][ipt].at(id0));
                  if(deltaPhi>-ana::PI && deltaPhi<-ana::PI/2.) deltaPhi += 2*ana::PI;

                  const double DCA = evt->DecayL3D(id0) * sin(evt->PointingAngle3D(id0));
                  string strDCA = ana::findDCA(DCA);
                  (hSignal_D0[imass][ipt])[strDCA]->Fill(deltaEta, deltaPhi, 
                        1./nMult_trg_eff_d0[imass][ipt]/effks/effweight_ass);
               }
            }
         }
      }
      //


      // mixed events
      // begin filling background
      // if having enough, erase the first one, fill the current, fill the background
      if(pVectList_ass[iz].size() == ana::nMixedEvts){
         for(int imass=0; imass<ana::nMass; imass++){
            for(int ipt=0; ipt<ana::nPt; ipt++){
               for(unsigned int id0=0; id0<nMult_trg_raw_d0[imass][ipt]; id0++){
                  // simultaneously read momentum and efficiency
                  auto ievt_eff_ass = effVectList_ass[iz].begin();
                  for(auto& ievt_p_ass : pVectList_ass[iz]){

                     if(ipt!=ana::findPtBin(pVect_trg_d0[imass][ipt].at(id0).Pt())) std::cout << "pT bin error" << std::endl;
                     double effks = 1.0;
                     unsigned int n_ass = ievt_p_ass.size();
                     if(n_ass!=ievt_eff_ass->size()){
                        cout << "wrong" << endl;
                        break;
                     }
                     for(unsigned int iass=0; iass<n_ass; iass++){
                        TVector3 pvector_ass = ievt_p_ass.at(iass);
                        double effweight_ass = ievt_eff_ass->at(iass);
                        if(ana::rejectDaughter_){
                           if(fabs(pvector_ass.Eta() - pVect_dau1_d0[imass][ipt].at(id0).Eta())<0.03
                                 && fabs(pvector_ass.Phi() - pVect_dau1_d0[imass][ipt].at(id0).Phi())<0.03)
                              continue;
                           if(fabs(pvector_ass.Eta() - pVect_dau2_d0[imass][ipt].at(id0).Eta())<0.03
                                 && fabs(pvector_ass.Phi() - pVect_dau2_d0[imass][ipt].at(id0).Phi())<0.03)
                              continue;
                        }
                        double deltaEta = pvector_ass.Eta() - pVect_trg_d0[imass][ipt].at(id0).Eta();
                        double deltaPhi = pvector_ass.DeltaPhi(pVect_trg_d0[imass][ipt].at(id0));
                        if(deltaPhi>-ana::PI && deltaPhi<-ana::PI/2.) deltaPhi += 2*ana::PI;

                        const double DCA = evt->DecayL3D(id0) * sin(evt->PointingAngle3D(id0));
                        string strDCA = ana::findDCA(DCA);
                        (hBackground_D0[imass][ipt])[strDCA]->Fill(deltaEta, deltaPhi, 
                           1./nMult_trg_eff_d0[imass][ipt]/effks/effweight_ass);
                     }
                     ievt_eff_ass++;
                  }
               }
            }
         }
      }

      // clear all and refresh list/buffer
      if(pVectList_ass[iz].size() == ana::nMixedEvts){
         // update the list/buffer
         pVectList_ass[iz].erase(pVectList_ass[iz].begin());
         pVectList_ass[iz].push_back(pVect_ass);
         effVectList_ass[iz].erase(effVectList_ass[iz].begin());
         effVectList_ass[iz].push_back(effVect_ass);
      }else{
      // if not having enough events, keep filling the list of momentum vector of associated particles
         // fill the list/buffer
         pVectList_ass[iz].push_back(pVect_ass);
         effVectList_ass[iz].push_back(effVect_ass);
      }

      for(int imass=0; imass<ana::nMass; imass++){
         for(int ipt=0; ipt<ana::nPt; ipt++){
            pVect_trg_d0[imass][ipt].clear();
            pVect_dau1_d0[imass][ipt].clear();
            pVect_dau2_d0[imass][ipt].clear();
         }
      }
      pVect_ass.clear();
      effVect_ass.clear();
   }
   std::cout << "completed loop" << std::endl;
   std::cout << skip << " events are skipped" << std::endl;

   ts.Stop();
   ts.Print();

   fout.cd();

   // start writing output
   TH2D* hCorrected_D0_largeDCA[ana::nMass][ana::nPt];
   TH2D* hCorrected_D0_smallDCA[ana::nMass][ana::nPt];
   for(int imass=0; imass<ana::nMass; imass++){
      for(int ipt=0; ipt<ana::nPt; ipt++){
         hCorrected_D0_largeDCA[imass][ipt] = (TH2D*) hSignal_D0[imass][ipt].at("largeDCA")->Clone();
         hCorrected_D0_largeDCA[imass][ipt]->SetName(Form("hCorrected_D0_mass%d_pt%d_largeDCA", imass, ipt));

         hCorrected_D0_smallDCA[imass][ipt] = (TH2D*) hSignal_D0[imass][ipt].at("smallDCA")->Clone();
         hCorrected_D0_smallDCA[imass][ipt]->SetName(Form("hCorrected_D0_mass%d_pt%d_smallDCA", imass, ipt));

         int etaBin = hBackground_D0[imass][ipt].at("largeDCA")->GetXaxis()->FindBin(0.);
         int phiBin = hBackground_D0[imass][ipt].at("largeDCA")->GetYaxis()->FindBin(0.);
         hCorrected_D0_largeDCA[imass][ipt]->Divide(hBackground_D0[imass][ipt].at("largeDCA"));
         hCorrected_D0_largeDCA[imass][ipt]->Scale(hBackground_D0[imass][ipt].at("largeDCA")->GetBinContent(etaBin, phiBin));

         hCorrected_D0_smallDCA[imass][ipt]->Divide(hBackground_D0[imass][ipt].at("smallDCA"));
         hCorrected_D0_smallDCA[imass][ipt]->Scale(hBackground_D0[imass][ipt].at("smallDCA")->GetBinContent(etaBin, phiBin));
      }
   }
   
   hMult->Write();
   hMult_ass->Write();
   for(int ipt=0; ipt<ana::nPt; ipt++){
      hKET_D0[ipt]->Write(); 
      hPt_D0[ipt]->Write();
      hEta_D0[ipt]->Write();
      hRapidity_D0[ipt]->Write();
      for(int imass=0; imass<ana::nMass; imass++){
         hMass_D0[imass][ipt]->Write();
         hMult_raw_D0[imass][ipt]->Write();
         hMult_eff_D0[imass][ipt]->Write();
         for(auto& h : hSignal_D0[imass][ipt]) h.second->Write();
         for(auto& h : hBackground_D0[imass][ipt]) h.second->Write();
         hCorrected_D0_largeDCA[imass][ipt]->Write();
         hCorrected_D0_smallDCA[imass][ipt]->Write();
      }
   }

   delete hMult;
   delete hMult_ass;
   for(int ipt=0; ipt<ana::nPt; ipt++){
      delete hKET_D0[ipt]; 
      delete hPt_D0[ipt];
      delete hEta_D0[ipt];
      delete hRapidity_D0[ipt];
      for(int imass=0; imass<ana::nMass; imass++){
         delete hMass_D0[imass][ipt];
         delete hMult_raw_D0[imass][ipt];
         delete hMult_eff_D0[imass][ipt];
         for(auto& h : hSignal_D0[imass][ipt]) delete h.second;
         for(auto& h : hBackground_D0[imass][ipt]) delete h.second;
         delete hCorrected_D0_largeDCA[imass][ipt];
         delete hCorrected_D0_smallDCA[imass][ipt];
      }
   }

   delete evt;


   return 0;
}

void setBranchStatus(Event* evt)
{
   evt->SetBranchStatus("Ntrkoffline", 1);
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
   //if(!passD0KinematicCuts(event, icand)) return false;
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

//bool passD0KinematicCuts(Event* event, const int& icand)
//{
   // eta is not available in the TTree
   //bool passEta = fabs(event->Eta(icand)<1000.);
//   bool passEta = fabs(event->Eta(icand)<1.5);
//   return passEta;
//}

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
