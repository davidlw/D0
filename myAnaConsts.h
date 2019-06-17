#ifndef _MYANACONSTS_h
#define _MYANACONSTS_h
#include <string>
namespace ana{
   const int nuofY = 2;
   const float ybin[nuofY+1] = {0, 1.0, 2.0};
   const std::string whichtree[] = {"d0ana", "npd0ana", "npd0ana1", "npd0ana2"};

   const double pTMin = 6.0;
   const double pTMax = 7.0;
   const double yMin  = 0.0;
   const double yMax  = 1.0;

   const double mvaMin = 0.26;
   const double mvaMax = 0.40;
   const double mvaStep = 0.02;

   const int nuofpt = 9;
//   const double ptbin[nuofpt+1]={1.5, 2.4, 3.0, 3.5, 4.2, 5.0, 6.0, 7.0, 8.0}; // GeV
   const double ptbin[nuofpt+1]={1.5, 2.2, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10., 20.}; // GeV

   const int nuofDca = 12;
   const double dcaBin[nuofDca+1] = {0, 0.002, 0.004, 0.006, 0.008, 0.012, 0.016, 0.022, 0.028, 0.036, 0.044, 0.06, 0.08}; // cm
   const int nDca = 80;
   const double dcaMin = 0;
   const double dcaMax = 0.08;

   const double fit_range_low = 1.7;
   const double fit_range_high = 2.0;
   const double D0_mass = 1.8648;

   // variables related to side band subtraction 
   // make sure the variables located at the bin edge, and within mass range [1.7, 2.0], 
   // since we use both TH1 Integral and TF1 Integral
   // after convincing they satisfy the requirement, I add the the offset to get the correct bin of histogram
   const double sigma = 0.015;
   //const double left_side_min = 1.7; 
   const double left_side_min = 1.76; 
   const double left_side_max = 1.82;
   const double peak_min = 1.865 - 3*0.015; // 1.82 
   const double peak_max = 1.865 + 3*0.015; // 1.91
   const double right_side_min = 1.91;
   const double right_side_max = 2.0;
}


int whichPt(const double pT){
   for(int ipt=0; ipt<ana::nuofpt; ipt++){
      if(ana::ptbin[ipt]<pT && pT<ana::ptbin[ipt+1]) return ipt;
   }
   return -1;
} 

int whichY(const double y){
   for(int iy=0; iy<ana::nuofY; iy++){
      if(ana::ybin[iy]<fabs(y) && fabs(y)<ana::ybin[iy+1]) return iy;
   }
   return -1;
} 
#endif
