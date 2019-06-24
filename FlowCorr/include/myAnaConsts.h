#ifndef __MYANACONST_h__
#define __MYANACONST_h__

#include <iostream>
#include <string>
#include "TMath.h"

namespace ana{

   const double PI = TMath::Pi();

   const int nEtaBin = 33;
   const double etaBegin = -4.95;
   const double etaEnd = 4.95;

   const int nPhiBin = 32;
//   const double phiBegin = -(0.5-1.0/32)*PI;
//   const double phiEnd = (1.5-1.0/32)*PI;
   const double phiBegin = -0.5*PI;
   const double phiEnd = 1.5*PI;

   const double etaMin_ass_ = -2.4;
   const double etaMax_ass_ = 2.4;
   const double ptMin_ass_ = 0.3;
   const double ptMax_ass_ = 3.0;
   const double multMin_ = 185;
   const double multMax_ = 250;
   const bool   rejectDaughter_ = true;
   const double d0_eta_ = 9999.0;
   const double d0_dau_eta_ = 1.5;
   const int    d0_dau_nhit_ = 11;
   const double d0_dau_pterr_ = 0.1;
   const double d0_dau_pt_ = 0.7;
   const double d0_rapidity_min_ = -1.;
   const double d0_rapidity_max_ = 1;

   const int nMass = 14;
   const int nPt = 3;
   const int nY = 1;
   const double massbin[nMass+1] = {1.70, 1.74, 1.78, 1.8, 1.82, 1.84, 1.85, 1.86, 1.865, 
      1.87, 1.88, 1.9, 1.92, 1.96, 2.0};
   //const double massbin[nMass+1] = {1.86, 1.865};
   const double ptbin[nPt+1] = {2., 4., 6., 8.};
   const double ybin[nY+1] = {0., 1.};
   const double mvaCut[nPt] = {
      0.58, // 2 - 4 GeV
      0.46, // 4 - 6 GeV
      0.28  // 6 - 8 GeV
   };

   const double dcaSep = 0.008;

   const unsigned int nMixedEvts = 10;

   const unsigned int nZ_Vtx_ = 15;
   const double z_vtx_bin_[nZ_Vtx_+1] = {
      -15, -13, -11, -9, -7, -5, -3, -1,
      1, 3, 5, 7, 9, 11, 13, 15
   };

   std::string findDCA(const double&);
   int findMassBin(const double&);
   int findPtBin(const double&);
   int findYBin(const double&);
   int findZVtxBin(const float&);
};

#endif
