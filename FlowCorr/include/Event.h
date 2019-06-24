//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun 17 16:23:04 2019 by ROOT version 6.06/01
// from TTree Event/Event
// found on file: /storage1/users/wl33/D0Trees/Data/pPb_Tree_D0_default_BDTCutNewAndTrack_b1_v2.root
//////////////////////////////////////////////////////////

#ifndef Event_h
#define Event_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TObject.h>

// Header file for the classes stored in the TTree if any.

class Event : public TObject{

public:
   Event(TTree *d0Collection=0, TTree *trackCollection=0);
   virtual ~Event();
   virtual Long64_t GetEntries();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual void     Init(TTree *d0Collection, TTree *trackCollection);
   virtual void     SetBranchStatus(const char* bname, Bool_t status);
   virtual Bool_t   GetBranchStatus(const char* bname);
   virtual Bool_t   Notify();

   virtual Int_t     nTrkOffline()  const {return Ntrkoffline;}
   virtual Int_t     nPixel()       const {return Npixel;}
   virtual Float_t   HFSumET()      const {return HFsumET;}
   virtual Float_t   BestVtxX()     const {return bestvtxX_D0;}
   virtual Float_t   BestVtxY()     const {return bestvtxY_D0;}
   virtual Float_t   BestVtxZ()     const {return bestvtxZ_D0;}
   virtual Int_t     CandSize()     const {return candSize;}
   virtual Float_t   Pt(const unsigned int icand)           const {return pT[icand];}
   virtual Float_t   Y(const unsigned int icand)            const {return y[icand];}
   virtual Float_t   Mass(const unsigned int icand)         const {return mass[icand];}
   virtual Float_t   Mva(const unsigned int icand)          const {return mva[icand];}
   virtual Float_t   Flavor(const unsigned int icand)       const {return flavor[icand];}
   virtual Float_t   Eta(const unsigned int icand)          const {return eta[icand];}
   virtual Float_t   vtxProb(const unsigned int icand)      const {return VtxProb[icand];}
   virtual Float_t   CosPointingAngle3D(const unsigned int icand) 
                                    const {return m3DCosPointingAngle[icand];}
   virtual Float_t  PointingAngle3D(const unsigned int icand) 
                                    const {return m3DPointingAngle[icand];}
   virtual Float_t  CosPointingAngle2D(const unsigned int icand) 
                                    const {return m2DCosPointingAngle[icand];}
   virtual Float_t  PointingAngle2D(const unsigned int icand) 
                                    const {return m2DPointingAngle[icand];}
   virtual Float_t  DecayLSig3D(const unsigned int icand)  const {return m3DDecayLengthSignificance[icand];}
   virtual Float_t  DecayL3D(const unsigned int icand)     const {return m3DDecayLength[icand];}
   virtual Float_t  DecayLSig2D(const unsigned int icand)  const {return m2DDecayLengthSignificance[icand];}
   virtual Float_t  DecayL2D(const unsigned int icand)     const {return m2DDecayLength[icand];}
   virtual Float_t  zDCASigD1(const unsigned int icand)    const {return zDCASignificanceDaugther1[icand];}
   virtual Float_t  xyDCASigD1(const unsigned int icand)   const {return xyDCASignificanceDaugther1[icand];}
   virtual Float_t  nHitD1(const unsigned int icand)       const {return NHitD1[icand];}
   virtual Bool_t   highPurityD1(const unsigned int icand) const {return HighPurityDaugther1[icand];}
   virtual Float_t  PtD1(const unsigned int icand)         const {return pTD1[icand];}
   virtual Float_t  PtErrD1(const unsigned int icand)      const {return pTerrD1[icand];}
   virtual Float_t  etaD1(const unsigned int icand)        const {return EtaD1[icand];}
   virtual Float_t  phiD1(const unsigned int icand)        const {return PhiD1[icand];}
   virtual Float_t  DedxD1(const unsigned int icand)       const {return dedxHarmonic2D1[icand];}
   virtual Float_t  zDCASigD2(const unsigned int icand)    const {return zDCASignificanceDaugther2[icand];}
   virtual Float_t  xyDCASigD2(const unsigned int icand)   const {return xyDCASignificanceDaugther2[icand];}
   virtual Float_t  nHitD2(const unsigned int icand)       const {return NHitD2[icand];}
   virtual Bool_t   highPurityD2(const unsigned int icand) const {return HighPurityDaugther2[icand];}
   virtual Float_t  PtD2(const unsigned int icand)         const {return pTD2[icand];}
   virtual Float_t  PtErrD2(const unsigned int icand)      const {return pTerrD2[icand];}
   virtual Float_t  phiD2(const unsigned int icand)        const {return PhiD2[icand];}
   virtual Float_t  etaD2(const unsigned int icand)        const {return EtaD2[icand];}
   virtual Float_t  DedxD2(const unsigned int icand)       const {return dedxHarmonic2D2[icand];}

   virtual UInt_t    runNb()        const {return RunNb;}
   virtual UInt_t    lsNb()         const {return LSNb;}
   virtual UInt_t    eventNb()      const {return EventNb;}
   virtual UInt_t    CandSizeTrk()  const {return candSizeTRK;}
   virtual Float_t   PtTrk(const unsigned int icand)        const {return (float)pTTRK[icand]/100.;}
   virtual Float_t   EtaTrk(const unsigned int icand)       const {return (float)etaTRK[icand]/100.;}
   virtual Float_t   PhiTrk(const unsigned int icand)       const {return (float)phiTRK[icand]/100.;}
   virtual Float_t   WeightTrk(const unsigned int icand)    const {return (float)weightTRK[icand]/100.;}

protected :
   TTree          *fChain_D0;   //!pointer to the analyzed TTree or TChain of D0 candidates
   TTree          *fChain_track;   //!pointer to the analyzed TTree or TChain of tracks

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Ntrkoffline;
   Int_t           Npixel;
   Float_t         HFsumET;
   Float_t         bestvtxX_D0;
   Float_t         bestvtxY_D0;
   Float_t         bestvtxZ_D0;
   Int_t           candSize;
   Float_t         pT[50];   //[candSize]
   Float_t         y[50];   //[candSize]
   Float_t         eta[50]; //[cadSize]
   Float_t         phi[50];   //[candSize]
   Float_t         mass[50];   //[candSize]
   Float_t         mva[50];   //[candSize]
   Float_t         flavor[50];   //[candSize]
   Float_t         VtxProb[50];   //[candSize]
   Float_t         m3DCosPointingAngle[50];   //[candSize]
   Float_t         m3DPointingAngle[50];   //[candSize]
   Float_t         m2DCosPointingAngle[50];   //[candSize]
   Float_t         m2DPointingAngle[50];   //[candSize]
   Float_t         m3DDecayLengthSignificance[50];   //[candSize]
   Float_t         m3DDecayLength[50];   //[candSize]
   Float_t         m2DDecayLengthSignificance[50];   //[candSize]
   Float_t         m2DDecayLength[50];   //[candSize]
   Float_t         zDCASignificanceDaugther1[50];   //[candSize]
   Float_t         xyDCASignificanceDaugther1[50];   //[candSize]
   Float_t         NHitD1[50];   //[candSize]
   Bool_t          HighPurityDaugther1[50];   //[candSize]
   Float_t         pTD1[50];   //[candSize]
   Float_t         pTerrD1[50];   //[candSize]
   Float_t         EtaD1[50];   //[candSize]
   Float_t         PhiD1[50];   //[candSize]
   Float_t         dedxHarmonic2D1[50];   //[candSize]
   Float_t         zDCASignificanceDaugther2[50];   //[candSize]
   Float_t         xyDCASignificanceDaugther2[50];   //[candSize]
   Float_t         NHitD2[50];   //[candSize]
   Bool_t          HighPurityDaugther2[50];   //[candSize]
   Float_t         pTD2[50];   //[candSize]
   Float_t         pTerrD2[50];   //[candSize]
   Float_t         EtaD2[50];   //[candSize]
   Float_t         PhiD2[50];   //[candSize]
   Float_t         dedxHarmonic2D2[50];   //[candSize]


   UInt_t          RunNb;
   UInt_t          LSNb;
   UInt_t          EventNb;
   Float_t         bestvtxX_trk;
   Float_t         bestvtxY_trk;
   Float_t         bestvtxZ_trk;
   UInt_t          candSizeTRK;
   UInt_t          pTTRK[553];   //[candSizeTRK]
   Short_t         etaTRK[553];   //[candSizeTRK]
   Short_t         phiTRK[553];   //[candSizeTRK]
   UInt_t          weightTRK[553];   //[candSizeTRK]

   // List of branches
   TBranch        *b_Ntrkoffline;   //!
   TBranch        *b_Npixel;   //!
   TBranch        *b_HFsumET;   //!
   TBranch        *b_bestvtxX_D0;   //!
   TBranch        *b_bestvtxY_D0;   //!
   TBranch        *b_bestvtxZ_D0;   //!
   TBranch        *b_candSize;   //!
   TBranch        *b_pT;   //!
   TBranch        *b_y;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_mass;   //!
   TBranch        *b_mva;   //!
   TBranch        *b_flavor;   //!
   TBranch        *b_VtxProb;   //!
   TBranch        *b_3DCosPointingAngle;   //!
   TBranch        *b_3DPointingAngle;   //!
   TBranch        *b_2DCosPointingAngle;   //!
   TBranch        *b_2DPointingAngle;   //!
   TBranch        *b_3DDecayLengthSignificance;   //!
   TBranch        *b_3DDecayLength;   //!
   TBranch        *b_2DDecayLengthSignificance;   //!
   TBranch        *b_2DDecayLength;   //!
   TBranch        *b_zDCASignificanceDaugther1;   //!
   TBranch        *b_xyDCASignificanceDaugther1;   //!
   TBranch        *b_NHitD1;   //!
   TBranch        *b_HighPurityDaugther1;   //!
   TBranch        *b_pTD1;   //!
   TBranch        *b_pTerrD1;   //!
   TBranch        *b_EtaD1;   //!
   TBranch        *b_PhiD1;   //!
   TBranch        *b_dedxHarmonic2D1;   //!
   TBranch        *b_zDCASignificanceDaugther2;   //!
   TBranch        *b_xyDCASignificanceDaugther2;   //!
   TBranch        *b_NHitD2;   //!
   TBranch        *b_HighPurityDaugther2;   //!
   TBranch        *b_pTD2;   //!
   TBranch        *b_pTerrD2;   //!
   TBranch        *b_EtaD2;   //!
   TBranch        *b_PhiD2;   //!
   TBranch        *b_dedxHarmonic2D2;   //!

   TBranch        *b_RunNb;   //!
   TBranch        *b_LSNb;   //!
   TBranch        *b_EventNb;   //!
   TBranch        *b_bestvtxX_trk;   //!
   TBranch        *b_bestvtxY_trk;   //!
   TBranch        *b_bestvtxZ_trk;   //!
   TBranch        *b_candSizeTRK;   //!
   TBranch        *b_pTTRK;   //!
   TBranch        *b_etaTRK;   //!
   TBranch        *b_phiTRK;   //!
   TBranch        *b_weightTRK;   //!
};
#endif
