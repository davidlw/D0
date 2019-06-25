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

   TH1D* hQ4Cos_Ref;
   TH1D* hQ4Sin_Ref;
   TH1D* hQ2Cos_Ref;
   TH1D* hQ2Sin_Ref;

   hMult = new TH1D("hMult", "", 600, 0, 600);
   hMult_ass = new TH1D("hMult_ass", "", 600, 0, 600);
   hQ4Cos_Ref = new TH1D("hQ4Cos_Ref","", 4000, -1., 1.);
   hQ4Sin_Ref = new TH1D("hQ4Sin_Ref","", 4000, -1., 1.);
   hQ2Cos_Ref = new TH1D("hQ2Cos_Ref","", 4000, -1., 1.);
   hQ2Sin_Ref = new TH1D("hQ2Sin_Ref","", 4000, -1., 1.);

   // declare vectors
   vector<TVector3> pVect_ass;
   vector<float> effVect_ass;

   // start timing
   TStopwatch ts;
   ts.Start();

   // loop
   //
   // temporary vectors
   TVector3 p_ass(0, 0, 0);

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
      
      double sumcosn_ass=0;
      double sumsinn_ass=0;
      double sumcosn_ass_p=0;
      double sumsinn_ass_p=0;
      double sumweight_ass_p=0;
      double sumcosn_ass_m=0;
      double sumsinn_ass_m=0;
      double sumweight_ass_m=0;

      for(unsigned int iass=0; iass<nMult_ass; iass++){
         TVector3 pvector_ass = pVect_ass.at(iass);
//       double effweight_ass = effVect_ass.at(iass);
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
//         sumweight_ass = sumweight_ass_p + sumweight_ass_m;
      }

      TComplex q_ass_p(sumcosn_ass_p,sumsinn_ass_p);
      TComplex q_ass_m(sumcosn_ass_m,sumsinn_ass_m);
      TComplex q_ass(sumcosn_ass,sumsinn_ass);

      TComplex Q4 = q_ass_p*q_ass_p*TComplex::Conjugate(q_ass_m)*TComplex::Conjugate(q_ass_m);
      TComplex Q2 = q_ass_p*TComplex::Conjugate(q_ass_m);

      double weight4_tot = sumweight_ass_p*sumweight_ass_p*sumweight_ass_m*sumweight_ass_m;
      double weight2_tot = sumweight_ass_p*sumweight_ass_m;

      hQ4Cos_Ref->Fill(Q4.Re()/weight4_tot,weight4_tot);
      hQ4Sin_Ref->Fill(Q4.Im()/weight4_tot,weight4_tot);
      hQ2Cos_Ref->Fill(Q2.Re()/weight2_tot,weight2_tot);
      hQ2Sin_Ref->Fill(Q2.Im()/weight2_tot,weight2_tot);

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
   hQ4Cos_Ref->Write();
   hQ4Sin_Ref->Write();
   hQ2Cos_Ref->Write();
   hQ2Sin_Ref->Write();

   delete hMult;
   delete hMult_ass;
   delete hQ4Cos_Ref;
   delete hQ4Sin_Ref;
   delete hQ2Cos_Ref;
   delete hQ2Sin_Ref;

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
