#ifndef Event_cxx
#define Event_cxx

#include <iostream>

#include "TMath.h"

#include "Event.h"

Event::Event(TTree *d0Collection, TTree *trackCollection) : fChain_D0(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (d0Collection == 0 ) {
      std::cout << "empty tree of D0 candidates" << std::endl;
   }else if(trackCollection == 0){
      std::cout << "empty tree of tracks" << std::endl;
   }else {
      Init(d0Collection, trackCollection);
   }
}

Event::~Event()
{
   if (!fChain_D0) return;
   delete fChain_D0->GetCurrentFile();
}

Int_t Event::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain_D0) return 0;
   Int_t bytes = fChain_D0->GetEntry(entry);
   if(bytes==0) return 0;
   if(bytes==-1) return -1;
   /*
   if( TMath::Abs(bestvtxX_trk - bestvtxX_D0) > 1e-5
      || TMath::Abs(bestvtxY_trk - bestvtxY_D0) > 1e-5
      || TMath::Abs(bestvtxZ_trk - bestvtxZ_D0) > 1e-5
         ){
      std::cout << bestvtxZ_D0 << " " << bestvtxZ_trk << std::endl;
      std::cout << "unmatched" << std::endl;
      return -1;
   }
   */
   return bytes;
}

Long64_t Event::GetEntries()
{
   if(fChain_D0 == 0 || fChain_track == 0) return 0; 
   if(fChain_D0->GetEntries() != fChain_track->GetEntries()) return 0;
   return fChain_D0->GetEntries();
}

void Event::Init(TTree *d0Collection, TTree *trackCollection)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!d0Collection) return;
   if (!trackCollection) return;

   fChain_D0 = d0Collection;
   fChain_D0->SetMakeClass(1);

   fChain_track = trackCollection;
   fChain_track->SetMakeClass(1);

   fChain_D0->AddFriend(fChain_track, "tracks");

   fChain_D0->SetBranchAddress("Ntrkoffline", &Ntrkoffline, &b_Ntrkoffline);
   fChain_D0->SetBranchAddress("Npixel", &Npixel, &b_Npixel);
   fChain_D0->SetBranchAddress("HFsumET", &HFsumET, &b_HFsumET);
   fChain_D0->SetBranchAddress("bestvtxX", &bestvtxX_D0, &b_bestvtxX_D0);
   fChain_D0->SetBranchAddress("bestvtxY", &bestvtxY_D0, &b_bestvtxY_D0);
   fChain_D0->SetBranchAddress("bestvtxZ", &bestvtxZ_D0, &b_bestvtxZ_D0);
   fChain_D0->SetBranchAddress("candSize", &candSize, &b_candSize);
   fChain_D0->SetBranchAddress("pT", pT, &b_pT);
   fChain_D0->SetBranchAddress("y", y, &b_y);
   fChain_D0->SetBranchAddress("phi", phi, &b_phi);
   fChain_D0->SetBranchAddress("mass", mass, &b_mass);
   fChain_D0->SetBranchAddress("mva", mva, &b_mva);
   fChain_D0->SetBranchAddress("flavor", flavor, &b_flavor);
   fChain_D0->SetBranchAddress("VtxProb", VtxProb, &b_VtxProb);
   fChain_D0->SetBranchAddress("3DCosPointingAngle", m3DCosPointingAngle, &b_3DCosPointingAngle);
   fChain_D0->SetBranchAddress("3DPointingAngle", m3DPointingAngle, &b_3DPointingAngle);
   fChain_D0->SetBranchAddress("2DCosPointingAngle", m2DCosPointingAngle, &b_2DCosPointingAngle);
   fChain_D0->SetBranchAddress("2DPointingAngle", m2DPointingAngle, &b_2DPointingAngle);
   fChain_D0->SetBranchAddress("3DDecayLengthSignificance", m3DDecayLengthSignificance, &b_3DDecayLengthSignificance);
   fChain_D0->SetBranchAddress("3DDecayLength", m3DDecayLength, &b_3DDecayLength);
   fChain_D0->SetBranchAddress("2DDecayLengthSignificance", m2DDecayLengthSignificance, &b_2DDecayLengthSignificance);
   fChain_D0->SetBranchAddress("2DDecayLength", m2DDecayLength, &b_2DDecayLength);
   fChain_D0->SetBranchAddress("zDCASignificanceDaugther1", zDCASignificanceDaugther1, &b_zDCASignificanceDaugther1);
   fChain_D0->SetBranchAddress("xyDCASignificanceDaugther1", xyDCASignificanceDaugther1, &b_xyDCASignificanceDaugther1);
   fChain_D0->SetBranchAddress("NHitD1", NHitD1, &b_NHitD1);
   fChain_D0->SetBranchAddress("HighPurityDaugther1", HighPurityDaugther1, &b_HighPurityDaugther1);
   fChain_D0->SetBranchAddress("pTD1", pTD1, &b_pTD1);
   fChain_D0->SetBranchAddress("pTerrD1", pTerrD1, &b_pTerrD1);
   fChain_D0->SetBranchAddress("EtaD1", EtaD1, &b_EtaD1);
   fChain_D0->SetBranchAddress("PhiD1", PhiD1, &b_PhiD1);
   fChain_D0->SetBranchAddress("dedxHarmonic2D1", dedxHarmonic2D1, &b_dedxHarmonic2D1);
   fChain_D0->SetBranchAddress("zDCASignificanceDaugther2", zDCASignificanceDaugther2, &b_zDCASignificanceDaugther2);
   fChain_D0->SetBranchAddress("xyDCASignificanceDaugther2", xyDCASignificanceDaugther2, &b_xyDCASignificanceDaugther2);
   fChain_D0->SetBranchAddress("NHitD2", NHitD2, &b_NHitD2);
   fChain_D0->SetBranchAddress("HighPurityDaugther2", HighPurityDaugther2, &b_HighPurityDaugther2);
   fChain_D0->SetBranchAddress("pTD2", pTD2, &b_pTD2);
   fChain_D0->SetBranchAddress("pTerrD2", pTerrD2, &b_pTerrD2);
   fChain_D0->SetBranchAddress("EtaD2", EtaD2, &b_EtaD2);
   fChain_D0->SetBranchAddress("PhiD2", PhiD2, &b_PhiD2);
   fChain_D0->SetBranchAddress("dedxHarmonic2D2", dedxHarmonic2D2, &b_dedxHarmonic2D2);

   fChain_D0->SetBranchAddress("tracks.RunNb", &RunNb, &b_RunNb);
   fChain_D0->SetBranchAddress("tracks.LSNb", &LSNb, &b_LSNb);
   fChain_D0->SetBranchAddress("tracks.EventNb", &EventNb, &b_EventNb);
   fChain_D0->SetBranchAddress("tracks.bestvtxX", &bestvtxX_trk, &b_bestvtxX_trk);
   fChain_D0->SetBranchAddress("tracks.bestvtxY", &bestvtxY_trk, &b_bestvtxY_trk);
   fChain_D0->SetBranchAddress("tracks.bestvtxZ", &bestvtxZ_trk, &b_bestvtxZ_trk);
   fChain_D0->SetBranchAddress("tracks.candSizeTRK", &candSizeTRK, &b_candSizeTRK);
   fChain_D0->SetBranchAddress("tracks.pTTRK", pTTRK, &b_pTTRK);
   fChain_D0->SetBranchAddress("tracks.etaTRK", etaTRK, &b_etaTRK);
   fChain_D0->SetBranchAddress("tracks.phiTRK", phiTRK, &b_phiTRK);
   fChain_D0->SetBranchAddress("tracks.weightTRK", weightTRK, &b_weightTRK);

   fChain_D0->SetBranchStatus("*", 0);
   fChain_D0->SetBranchStatus("*bestvtx*", 1);
   fChain_D0->SetBranchStatus("tracks.bestvtxX", 1);
   fChain_D0->SetBranchStatus("tracks.bestvtxY", 1);
   fChain_D0->SetBranchStatus("tracks.bestvtxZ", 1);

   Notify();
}

void Event::SetBranchStatus(const char* bname, Bool_t status)
{
   fChain_D0->SetBranchStatus(bname, status);
}

Bool_t Event::GetBranchStatus(const char* bname)
{
   return fChain_D0->GetBranchStatus(bname);
}

Bool_t Event::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}
#endif
