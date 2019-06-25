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
#include "TComplex.h"

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

   TFile fout(TString::Format("fout_%s.root", datalist.c_str()), "recreate");

   TChain *chain_d0 = new TChain("d0ana/VertexCompositeNtuple");
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
   TH1D* hQ4Cos_D0[ana::nMass][ana::nPt];
   TH1D* hQ4Sin_D0[ana::nMass][ana::nPt];
   TH1D* hQ2Cos_D0[ana::nMass][ana::nPt];
   TH1D* hQ2Sin_D0[ana::nMass][ana::nPt];

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
         (hQ4Cos_D0[imass][ipt]) = new TH1D(Form("hQ4Cos_mass%d_pt%d", imass, ipt),
               "", 4000, -1., 1.);
         (hQ4Sin_D0[imass][ipt]) = new TH1D(Form("hQ4Sin_mass%d_pt%d", imass, ipt),
               "", 4000, -1., 1.);
         (hQ2Cos_D0[imass][ipt]) = new TH1D(Form("hQ2Cos_mass%d_pt%d", imass, ipt),
               "", 4000, -1., 1.);
         (hQ2Sin_D0[imass][ipt]) = new TH1D(Form("hQ2Sin_mass%d_pt%d", imass, ipt),
               "", 4000, -1., 1.);
      }
   }

   // declare vectors
   vector<TVector3> pVect_trg_d0[ana::nMass][ana::nPt];
   vector<TVector3> pVect_ass;
   vector<float> effVect_ass;

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

            double sumcosn_d0=0;
            double sumsinn_d0=0;
//            double sumweight_d0=0;
            double sumcosn_ass=0;
            double sumsinn_ass=0;
//            double sumweight_ass=0;
            double sumcosn_d0_p=0;
            double sumsinn_d0_p=0;
            double sumweight_d0_p=0;
            double sumcosn_ass_p=0;
            double sumsinn_ass_p=0;
            double sumweight_ass_p=0;
            double sumcosn_d0_m=0;
            double sumsinn_d0_m=0;
            double sumweight_d0_m=0;
            double sumcosn_ass_m=0;
            double sumsinn_ass_m=0;
            double sumweight_ass_m=0;

            for(unsigned int id0=0; id0<nMult_trg_d0; id0++){
               if(ipt!=ana::findPtBin(pVect_trg_d0[imass][ipt].at(id0).Pt())) std::cout << "pT bin error" << std::endl;
               double effks = 1.0;

               if(pVect_trg_d0[imass][ipt].at(id0).Eta()>0)
               {
                 sumcosn_d0_p += cos(2*pVect_trg_d0[imass][ipt].at(id0).Phi())/effks;
                 sumsinn_d0_p += sin(2*pVect_trg_d0[imass][ipt].at(id0).Phi())/effks;
                 sumweight_d0_p += 1.0;      
               }
               if(pVect_trg_d0[imass][ipt].at(id0).Eta()<0)
               {
                 sumcosn_d0_m += cos(2*pVect_trg_d0[imass][ipt].at(id0).Phi())/effks;
                 sumsinn_d0_m += sin(2*pVect_trg_d0[imass][ipt].at(id0).Phi())/effks;
                 sumweight_d0_m += 1.0;
               }              

               sumcosn_d0 = sumcosn_d0_p + sumcosn_d0_m;
               sumsinn_d0 = sumsinn_d0_p + sumsinn_d0_m;
//               sumweight_d0 = sumweight_d0_p + sumweight_d0_m;
            }

            for(unsigned int iass=0; iass<nMult_ass; iass++){
               TVector3 pvector_ass = pVect_ass.at(iass);
//               double effweight_ass = effVect_ass.at(iass);
               double effweight_ass = 1.0;

               if(pvector_ass.Eta()>0)
               {
                 sumcosn_ass_p += cos(2*pvector_ass.Phi())/effweight_ass;
                 sumsinn_ass_p += sin(2*pvector_ass.Phi())/effweight_ass;
                 sumweight_ass_p += 1.0/effweight_ass;
               }
               if(pvector_ass.Eta()<0)
               {
                 sumcosn_ass_m += cos(2*pvector_ass.Phi())/effweight_ass;
                 sumsinn_ass_m += sin(2*pvector_ass.Phi())/effweight_ass;
                 sumweight_ass_m += 1.0/effweight_ass;
               }
 
               sumcosn_ass = sumcosn_ass_p + sumcosn_ass_m;
               sumsinn_ass = sumsinn_ass_p + sumsinn_ass_m;
//               sumweight_ass = sumweight_ass_p + sumweight_ass_m;
            }

            TComplex q_d0_p(sumcosn_d0_p,sumsinn_d0_p); 
            TComplex q_ass_p(sumcosn_ass_p,sumsinn_ass_p);
            TComplex q_d0_m(sumcosn_d0_m,sumsinn_d0_m);
            TComplex q_ass_m(sumcosn_ass_m,sumsinn_ass_m);
            TComplex q_d0(sumcosn_d0,sumsinn_d0);
            TComplex q_ass(sumcosn_ass,sumsinn_ass);

            TComplex Q4_p = q_d0_p*q_ass_p*TComplex::Conjugate(q_ass_m)*TComplex::Conjugate(q_ass_m);
            TComplex Q4_m = q_d0_m*q_ass_m*TComplex::Conjugate(q_ass_p)*TComplex::Conjugate(q_ass_p);
            TComplex Q2_p = q_d0_p*TComplex::Conjugate(q_ass_m);
            TComplex Q2_m = q_d0_m*TComplex::Conjugate(q_ass_p);                                   
            TComplex Q4 = Q4_p + Q4_m;
            TComplex Q2 = Q2_p + Q2_m;

            double weight4_tot = sumweight_d0_p*sumweight_ass_p*sumweight_ass_m*sumweight_ass_m + sumweight_d0_m*sumweight_ass_m*sumweight_ass_p*sumweight_ass_p;
            double weight2_tot = sumweight_d0_p*sumweight_ass_m + sumweight_d0_m*sumweight_ass_p;

            hQ4Cos_D0[imass][ipt]->Fill(Q4.Re()/weight4_tot,weight4_tot);
            hQ4Sin_D0[imass][ipt]->Fill(Q4.Im()/weight4_tot,weight4_tot);
            hQ2Cos_D0[imass][ipt]->Fill(Q2.Re()/weight2_tot,weight2_tot);
            hQ2Sin_D0[imass][ipt]->Fill(Q2.Im()/weight2_tot,weight2_tot);
         }
      }

      for(int imass=0; imass<ana::nMass; imass++){
         for(int ipt=0; ipt<ana::nPt; ipt++){
            pVect_trg_d0[imass][ipt].clear();
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
         hQ4Cos_D0[imass][ipt]->Write();
         hQ4Sin_D0[imass][ipt]->Write();
         hQ2Cos_D0[imass][ipt]->Write();
         hQ2Sin_D0[imass][ipt]->Write();
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
         delete hQ4Cos_D0[imass][ipt];
         delete hQ4Sin_D0[imass][ipt];
         delete hQ2Cos_D0[imass][ipt];
         delete hQ2Sin_D0[imass][ipt];
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
