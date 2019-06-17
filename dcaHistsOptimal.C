#include "myAnaConsts.h"
#include "massfitting.C"
#include <iostream>
#include <memory>

TH1D* hDcaMCAll(TH3*, TH1&, const double&, const std::string&);
TH1D* hDcaData(TH3*, TH1&, const double&, const std::string&, const double&, const double&);
void hDcaDataSignal(TH3*, std::map<std::string, TH3*>, std::map<std::string, TH3*>, TH1&, const double&, const std::string&);

void dcaHistsOptimal(int mode = 2)
//void dcaHistsOptimal(int mode = 1)
{
   TH1::SetDefaultSumw2();

   //TGaxis::SetMaxDigits(3);
   gStyle->SetOptStat(0);
   //TFile* f1 = new TFile(Form("%s_hists.root", ana::whichtree[mode].c_str()));
   std::unique_ptr<TFile> f1 = std::unique_ptr<TFile>(new TFile(Form("%s_hists.root", ana::whichtree[mode].c_str())));
   std::map<std::string, TH3*> hDcaVsMassAndMvaPD0;
   std::map<std::string, TH3*> hDcaVsMassAndMvaNPD0;
   hDcaVsMassAndMvaPD0["h_match_unswap"] = (TH3D*) f1->Get("hDcaVsMassAndMvaPD0");
   hDcaVsMassAndMvaPD0["h_match_all"] = (TH3D*) f1->Get("hDcaVsMassAndMvaPD0_All");
   hDcaVsMassAndMvaNPD0["h_match_unswap"] = (TH3D*) f1->Get("hDcaVsMassAndMvaNPD0");
   hDcaVsMassAndMvaNPD0["h_match_all"] = (TH3D*) f1->Get("hDcaVsMassAndMvaNPD0_All");
   TH3D* hDcaVsMassAndMvaDataD0 = (TH3D*) f1->Get("hDcaVsMassAndMvaDataD0");

   TFile f2(Form("%s_dca_hists.root", ana::whichtree[mode].c_str()), "recreate");

   // create Pic to store pictures
   const std::string dirPic = Form("if [ ! -d \"Pic\" ]; then\n"
                              "    mkdir %sPic \n"
                              "fi", ana::whichtree[mode].c_str());
   gSystem->Exec(dirPic.c_str());

   for(int iMva=0; iMva<7; iMva++){
   //for(int iMva=2; iMva<4; iMva++){

      //create Pic/mva to store pictures
      std::string dirMvaPic(TString::Format( "if [ ! -d \"%sPic/mva%d\" ]; then\n"
                                 "    mkdir %sPic/mva%d \n"
                                 "fi", ana::whichtree[mode].c_str(), iMva, 
                                 ana::whichtree[mode].c_str(), iMva));
      gSystem->Exec(dirMvaPic.c_str());

      double mvaCut = (double)iMva*ana::mvaStep + ana::mvaMin;

      // DCA distribution of Prompt D0 
      TCanvas cPrompt("cPrompt", "", 550, 450);
      cPrompt.SetLogy();
      TH1D hDcaMCPD0("hDcaMCPD0", "hDCaMCPD0", ana::nuofDca, ana::dcaBin);
      TH1D* hDcaMCPD0_Proj = hDcaMCAll(hDcaVsMassAndMvaPD0["h_match_all"], hDcaMCPD0, mvaCut, "hprompt");
      hDcaMCPD0.Draw();
      //hDcaMCPD0_Proj->Scale(hDcaMCPD0.GetMaximum()/hDcaMCPD0_Proj->GetMaximum());
      hDcaMCPD0_Proj->SetLineColor(kRed);
      hDcaMCPD0_Proj->Draw("same");
      cPrompt.Print(Form("%sPic/mva%d/cPrompt.png", ana::whichtree[mode].c_str(), iMva));
      delete hDcaMCPD0_Proj;

      // DCA distribution of Non-Prompt D0
      TCanvas cNonPrompt("cNonPrompt", "", 550, 450);
      cNonPrompt.SetLogy();
      TH1D hDcaMCNPD0("hDcaMCNPD0", "hDCaMCNPD0", ana::nuofDca, ana::dcaBin);
      TH1D* hDcaMCNPD0_Proj = hDcaMCAll(hDcaVsMassAndMvaNPD0["h_match_all"], hDcaMCNPD0, mvaCut, "hnonprompt");
      hDcaMCNPD0.Draw();
      //hDcaMCNPD0_Proj->Scale(hDcaMCNPD0.GetMaximum()/hDcaMCNPD0_Proj->GetMaximum());
      hDcaMCNPD0_Proj->SetLineColor(kRed);
      hDcaMCNPD0_Proj->Draw("same");
      cNonPrompt.Print(Form("%sPic/mva%d/cNonPrompt.png", ana::whichtree[mode].c_str(), iMva));
      delete hDcaMCNPD0_Proj;

      // DCA distribution of signal + swap, by fitting mass;
      TCanvas cData("cData", "", 550, 450);
      TH1D hDcaDataD0("hDcaDataD0", "hDcaDataD0", ana::nuofDca, ana::dcaBin);
      std::string tmpName(Form("%sPic/mva%d/massPerDcaBin", ana::whichtree[mode].c_str(), iMva));
      hDcaDataSignal(hDcaVsMassAndMvaDataD0, hDcaVsMassAndMvaNPD0, hDcaVsMassAndMvaPD0, hDcaDataD0, mvaCut, tmpName);
      cData.cd();
      cData.SetLogy();
      //hDcaDataD0.Scale(1./hDcaDataD0.Integral("width"));
      hDcaDataD0.Draw();
      cData.Print(Form("%sPic/mva%d/cData.png", ana::whichtree[mode].c_str(), iMva));
      
      // DCA distribution of N1*signal + N2*swap, by side band subtraction
      // extract invariant mass
      int mvaBinMin = hDcaVsMassAndMvaNPD0["h_match_all"]->GetYaxis()->FindBin(mvaCut+0.1*ana::mvaStep); // 0.1 offset to make sure the correct bin is returned
      int mvaBinMax = hDcaVsMassAndMvaNPD0["h_match_all"]->GetYaxis()->GetNbins()+1;
      int label = mvaCut*100;
      TH1D* hMassNPD0 = hDcaVsMassAndMvaNPD0["h_match_unswap"]->ProjectionX(Form("hMassMCNPD0mva%d", label), mvaBinMin, mvaBinMax);
      TH1D* hMassNPD0All = hDcaVsMassAndMvaNPD0["h_match_all"]->ProjectionX(Form("hMassMCNPD0Allmva%d", label), mvaBinMin, mvaBinMax);
      TH1D* hMassPD0 = hDcaVsMassAndMvaPD0["h_match_unswap"]->ProjectionX(Form("hMassMCPD0mva%d", label), mvaBinMin, mvaBinMax);
      TH1D* hMassPD0All = hDcaVsMassAndMvaPD0["h_match_all"]->ProjectionX(Form("hMassMCPD0Allmva%d", label), mvaBinMin, mvaBinMax);
      TH1D* hMassData = hDcaVsMassAndMvaDataD0->ProjectionX(Form("hMassDataD0mva%d", label), mvaBinMin, mvaBinMax);


      // TFitResultPtr constructed by std::shared_ptr<TFitResult>, be free of it even if you do not delete it
      // fitting the mass
      // rescale the mass
      
      TFitResultPtr fitResultPtr;
      TF1 fMass;
      if(mode == 0) fMass = massfitting(hMassData, hMassPD0, hMassPD0All, TString::Format("fMassMva%d", iMva), fitResultPtr);
      if(mode == 1) fMass = massfitting(hMassData, hMassNPD0, hMassNPD0All, TString::Format("fMassMva%d", iMva), fitResultPtr);
      if(mode == 2) fMass = massfitting(hMassData, hMassNPD0, hMassNPD0All, TString::Format("fMassMva%d", iMva), fitResultPtr);

      // draw the fitting
      drawMassFitting(hMassData, fMass, TString::Format("%sPic/mva%d/MassFittingUnNormalized.png", ana::whichtree[mode].c_str(), iMva), 
            "4 < pT < 5GeV  |y|<1", TString::Format("MVA > %.2f", 0.4+0.02*iMva));

      // https://root.cern/doc/v616/TF1Helper_8cxx_source.html and https://root.cern.ch/doc/master/classTF1.html
      // help you understand how the integral error of TF1 is calculated
      // 1-dimensional function have five parameters, x_low, x_min, array of pars, array of covariance matrix, precision
      // array of pars would replace the original set of pars by last fitting unless no set of pars is passed
      // it is very important to pass in the correct set of pars and covariance matrix
      // explicitly pass the parameters to make sure anything won't go wrong
      // be aware that the fixed parameter would bring effect on the evaluation
      
      // get the covariance matrix
      auto covMat = fitResultPtr->GetCovarianceMatrix();

      // define the background function and then fix or set the parameters
      TF1 bkg("bkg", "[9] + [10]*x + [11]*x*x + [12]*x*x*x"
            "+ 0 *[0]*[1]*[2]*[3]*[4]*[5]*[6]*[7]*[8]", 1.7, 2.0);
      for(int ipar=0; ipar<8+1; ipar++){
         bkg.FixParameter(ipar, fMass.GetParameter(ipar));
      }
      for(int ipar=9; ipar<12+1; ipar++){
         bkg.SetParameter(ipar, fMass.GetParameter(ipar));
      }

      // calculate the yield of background in three regions
      std::map<std::string, double> yield;
      std::map<std::string, double> yieldError;
      yield["left_side"] = bkg.Integral(ana::left_side_min, ana::left_side_max);
      yield["peak"] = bkg.Integral(ana::peak_min, ana::peak_max);
      yield["right_side"] = bkg.Integral(ana::right_side_min, ana::right_side_max);
      yieldError["left_side"] = bkg.IntegralError(ana::left_side_min, ana::left_side_max, 0, covMat.GetMatrixArray());
      yieldError["peak"] = bkg.IntegralError(ana::peak_min, ana::peak_max, 0, covMat.GetMatrixArray());
      yieldError["right_side"] = bkg.IntegralError(ana::right_side_min, ana::right_side_max, 0, covMat.GetMatrixArray());
      for(auto& element : yield) std::cout << element.second << std::endl;
      for(auto& element : yieldError) std::cout << element.second << std::endl;
      
      // calculate the signal-to-sideband ratio
      // propagate the error
      double ratioLeft = yield["peak"] / yield["left_side"];
      double relRatioErrorLeft = std::sqrt(
            std::pow(yieldError["peak"]/yield["peak"], 2) + std::pow(yieldError["left_side"]/yield["left_side"], 2));
      double absRatioErrorLeft = ratioLeft* relRatioErrorLeft;

      double ratioRight = yield["peak"] / yield["right_side"];
      double relRatioErrorRight = std::sqrt(
            std::pow(yieldError["peak"]/yield["peak"], 2) + std::pow(yieldError["right_side"]/yield["right_side"], 2));
      double absRatioErrorRight = ratioRight* relRatioErrorRight;

      // extract dca of left side band 
      TCanvas cLeftSide("cLeftSide", "", 550, 450);
      cLeftSide.SetLogy();
      TH1D hDcaLeft("hDcaLeftSide", "hDcaLeftSide", ana::nuofDca, ana::dcaBin);
      TH1D* hDcaLeft_Proj = hDcaData(hDcaVsMassAndMvaDataD0, hDcaLeft, mvaCut, "hDcaLeftSide_Proj", ana::left_side_min, ana::left_side_max);
      hDcaLeft.Draw();
      hDcaLeft_Proj->Scale(hDcaLeft.GetMaximum()/hDcaLeft_Proj->GetMaximum());
      hDcaLeft_Proj->SetLineColor(kRed);
      hDcaLeft_Proj->Draw("SAME");
      cLeftSide.Print(Form("%sPic/mva%d/cLeftSide.png", ana::whichtree[mode].c_str(), iMva));
      hDcaLeft_Proj->Delete();

      // extract dca of right side band 
      TCanvas cRightSide("cRightSide", "", 550, 450);
      cRightSide.SetLogy();
      TH1D hDcaRight("hDcaRightSide", "hDcaRightSide", ana::nuofDca, ana::dcaBin);
      TH1D* hDcaRight_Proj = hDcaData(hDcaVsMassAndMvaDataD0, hDcaRight, mvaCut, "hDcaRightSide_Proj", ana::right_side_min, ana::right_side_max);
      hDcaRight.Draw();
      hDcaRight_Proj->Scale(hDcaRight.GetMaximum()/hDcaRight_Proj->GetMaximum());
      hDcaRight_Proj->SetLineColor(kRed);
      hDcaRight_Proj->Draw("SAME");
      cRightSide.Print(Form("%sPic/mva%d/cRightSide.png", ana::whichtree[mode].c_str(), iMva));
      hDcaRight_Proj->Delete();

      // extract dca of peak
      TCanvas cPeak("cPeak", "", 550, 450);
      cPeak.SetLogy();
      TH1D hDcaPeak("hDcaPeak", "hDcaPeak", ana::nuofDca, ana::dcaBin);
      TH1D* hDcaPeak_Proj = hDcaData(hDcaVsMassAndMvaDataD0, hDcaPeak, mvaCut, "hDcaPeak_Proj", ana::peak_min, ana::peak_max);
      hDcaPeak.Draw();
      hDcaPeak_Proj->Scale(hDcaPeak.GetMaximum()/hDcaPeak_Proj->GetMaximum());
      hDcaPeak_Proj->SetLineColor(kRed);
      hDcaPeak_Proj->Draw("SAME");
      cPeak.Print(Form("%sPic/mva%d/cPeak.png", ana::whichtree[mode].c_str(), iMva));
      hDcaPeak_Proj->Delete();

      // scale the left-side by signal-to-sideband ratio
      for(int iDca=0; iDca<ana::nuofDca; iDca++){
         double yield = ratioLeft * hDcaLeft.GetBinContent(iDca+1);
         double yieldErr = std::sqrt(std::pow(absRatioErrorLeft* hDcaLeft.GetBinContent(iDca+1), 2) +
               std::pow(ratioLeft* hDcaLeft.GetBinError(iDca+1), 2));
         hDcaLeft.SetBinContent(iDca+1, yield);
         hDcaLeft.SetBinError(iDca+1, yieldErr);
      }

      // scale the right-side by signal-to-sideband ratio
      for(int iDca=0; iDca<ana::nuofDca; iDca++){
         double yield = ratioRight * hDcaRight.GetBinContent(iDca+1);
         double yieldErr = std::sqrt(std::pow(absRatioErrorRight* hDcaRight.GetBinContent(iDca+1), 2) +
               std::pow(ratioRight* hDcaRight.GetBinError(iDca+1), 2));
         hDcaRight.SetBinContent(iDca+1, yield);
         hDcaRight.SetBinError(iDca+1, yieldErr);
      }

      // average two sides
      TH1D hDcaSides(hDcaLeft);
      hDcaSides.SetName("hDcaSides");
      hDcaSides.Add(&hDcaRight);
      hDcaSides.Scale(0.5);

      // compare dca of peak and side band
      TCanvas cPeakAndSides("cPeakAndSides", "", 550, 450);
      cPeakAndSides.SetLogy();
      hDcaPeak.SetMarkerStyle(20);
      hDcaPeak.Draw("E P");
      hDcaSides.Draw("E SAME");
      cPeakAndSides.Print(Form("%sPic/mva%d/cDcaPeakAndSides.png", ana::whichtree[mode].c_str(), iMva));

      TCanvas cSignal("cSignal", "", 550, 450);
      cSignal.SetLogy();
      TH1D* hDcaPeakSignal = (TH1D*) hDcaPeak.Clone();
      hDcaPeakSignal->SetName("hDcaPeakSignal");
      hDcaPeakSignal->Add(&hDcaSides, -1);
      hDcaPeakSignal->Draw();
      hDcaDataD0.Draw("same");
      cSignal.Print(Form("%sPic/mva%d/cSignalSideBand.png", ana::whichtree[mode].c_str(), iMva));

      hDcaPeakSignal->Write(Form("hDcaDataPeakD0mva%d", label));

      f2.cd();
      fMass.Write();

      hDcaDataD0.Write(Form("hDcaDataD0mva%d", label));
      hDcaMCNPD0.Write(Form("hDcaMCNPD0mva%d", label));
      hDcaMCPD0.Write(Form("hDcaMCPD0mva%d", label));

      hMassNPD0->Write();
      hMassNPD0->Delete();
      hMassNPD0All->Write();
      hMassNPD0All->Delete();
      hMassPD0->Write();
      hMassPD0->Delete();
      hMassPD0All->Write();
      hMassPD0All->Delete();
      hMassData->Write();
      hMassData->Delete();
   }
}

TH1D* hDcaMCAll(TH3* hDcaVsMassAndMva, TH1& hDcaMC, const double& mvaCut, const  std::string& proj_name)
{
   int mvaBinMin = hDcaVsMassAndMva->GetYaxis()->FindBin(mvaCut+0.1*ana::mvaStep); // 0.1 offset to make sure the correct bin is returned
   int mvaBinMax = hDcaVsMassAndMva->GetYaxis()->GetNbins()+1;
   TH1D* hDcaMCProj = hDcaVsMassAndMva->ProjectionZ(proj_name.c_str(), 1, 60, mvaBinMin, mvaBinMax); // x: mass, y: dca 
   const double norm = hDcaMCProj->Integral(0, hDcaMCProj->GetXaxis()->GetNbins()+1); // normalized constant
   hDcaMCProj->Scale(1./ (norm * hDcaMCProj->GetBinWidth(1))); // convert histogram to the P.D.F.
   for(int iDca=0; iDca<ana::nuofDca; iDca++){
      int binlw = hDcaMCProj->FindBin(ana::dcaBin[iDca] + hDcaMCProj->GetBinWidth(1)/2.); // binwidth/2 offset to make sure the correct bin is returned
      int binup = hDcaMCProj->FindBin(ana::dcaBin[iDca+1] - hDcaMCProj->GetBinWidth(1)/2.);  // - binwidth/2 offset due to the integral is taken within [binlw, binup]
      // calculate the contents and errors per dca bin
      double binwidth = ana::dcaBin[iDca+1] - ana::dcaBin[iDca];
      double err = 0.;
      double yield = hDcaMCProj->IntegralAndError(binlw, binup, err, "width"); // note the option "width"
      hDcaMC.SetBinContent(iDca+1, yield/binwidth);
      hDcaMC.SetBinError(iDca+1, err/binwidth);
   }
   return hDcaMCProj;
}

TH1D* hDcaData(TH3* hDcaVsMassAndMva, TH1& hDca, const double& mvaCut, const  std::string& proj_name, 
      const double& massMin, const double& massMax)
{
   int mvaBinMin = hDcaVsMassAndMva->GetYaxis()->FindBin(mvaCut+0.1*ana::mvaStep); // 0.1 offset to make sure the correct bin is returned
   int mvaBinMax = hDcaVsMassAndMva->GetYaxis()->GetNbins()+1;
   int massBinMin = hDcaVsMassAndMva->GetXaxis()->FindBin(massMin + 0.1*hDcaVsMassAndMva->GetXaxis()->GetBinWidth(1)); // 0.1 offset to get the correct bin, make sure massMin located at the bin edge
   int massBinMax = hDcaVsMassAndMva->GetXaxis()->FindBin(massMax - 0.1*hDcaVsMassAndMva->GetXaxis()->GetBinWidth(1)); // -0.1 offset to get the correct bin, make sure massMin located at the bin edge
   TH1D* hDcaDataProj = hDcaVsMassAndMva->ProjectionZ(proj_name.c_str(), massBinMin, massBinMax, mvaBinMin, mvaBinMax); // x: mass, y: dca 
   for(int iDca=0; iDca<ana::nuofDca; iDca++){
      int binlw = hDcaDataProj->FindBin(ana::dcaBin[iDca] + hDcaDataProj->GetBinWidth(1) * 0.1); // binwidth*0.1 offset to make sure the correct bin is returned
      int binup = hDcaDataProj->FindBin(ana::dcaBin[iDca+1] - hDcaDataProj->GetBinWidth(1) * 0.1);  // - binwidth*0.1 offset due to the integral is taken within [binlw, binup]
      // calculate the contents and errors per dca bin
      double binwidth = ana::dcaBin[iDca+1] - ana::dcaBin[iDca];
      double err = 0.;
      double yield = hDcaDataProj->IntegralAndError(binlw, binup, err);
      hDca.SetBinContent(iDca+1, yield/binwidth);
      hDca.SetBinError(iDca+1, err/binwidth);
   }
   
   return hDcaDataProj;
}

void hDcaDataSignal(TH3* hData, std::map<std::string, TH3*> hNP, std::map<std::string, TH3*> hP, TH1& hDataD0, const double& mvaCut, const std::string& tmpName)
{
   int mvaBinMin = hP["h_match_unswap"]->GetYaxis()->FindBin(mvaCut+0.1*ana::mvaStep); // 0.1 offset to make sure the correct bin is returned
   int mvaBinMax = hP["h_match_unswap"]->GetYaxis()->GetNbins()+1; // overflow of y axis, mva
   int mcDcaBinMin = 0; // underflow of the z axis, the DCA
   int mcDcaBinMax = hP["h_match_unswap"]->GetZaxis()->GetNbins()+1; // overflow of z axis, the DCA

   // x: mass, y: mva, z: dca
   TH1D* hMCPromptMass = hP["h_match_unswap"]->ProjectionX("hPromptMass", mvaBinMin, mvaBinMax, mcDcaBinMin, mcDcaBinMax); // projection, y-mva:[cut, overflow], z-DCA:[underflow, overflow]
   TH1D* hMCPromptMassAll = hP["h_match_all"]->ProjectionX("hPromptMassAll", mvaBinMin, mvaBinMax, mcDcaBinMin, mcDcaBinMax);
   TH1D* hMCNonPromptMass = hNP["h_match_unswap"]->ProjectionX("hNonPromptMass", mvaBinMin, mvaBinMax, mcDcaBinMin, mcDcaBinMax);
   TH1D* hMCNonPromptMassAll = hNP["h_match_all"]->ProjectionX("hNonPromptMassAll", mvaBinMin, mvaBinMax, mcDcaBinMin, mcDcaBinMax);

   TH1D* hTemp;
   TCanvas cTemp("cTemp", "", 550, 450);
   for(int iDca=0; iDca<ana::nuofDca; iDca++){
      int binlw = hData->GetZaxis()->FindBin(ana::dcaBin[iDca] + hData->GetZaxis()->GetBinWidth(1) * 0.2); // 0.2 offset to make sure the correct bin is returned
      int binup = hData->GetZaxis()->FindBin(ana::dcaBin[iDca+1] - hData->GetZaxis()->GetBinWidth(1) * 0.2); // -0.2 offset due to the integral is taken within [binlw, binup]
      hTemp = hData->ProjectionX("tmp", mvaBinMin, mvaBinMax, binlw, binup); // projection, y-mva:[cut, overflow], z-DCA:[dca[iDca], dca[iDca+1]]
      hTemp->Scale(1./hTemp->GetBinWidth(1));
      hTemp->GetYaxis()->SetRangeUser(0, hTemp->GetMaximum() * 1.3);
      hTemp->SetTitle(";Mass (GeV);Entries /(5 MeV)");
      TFitResultPtr fitResultPtr;
      TF1 f = massfitting(hTemp, hMCNonPromptMass, hMCNonPromptMassAll, TString::Format("iDca%d", iDca), fitResultPtr);
      drawMassFitting(hTemp, f, TString::Format("%s%d.png", tmpName.c_str(), iDca), "4 < pT < 5GeV  |y|<1",
            TString::Format("%2.f < DCA <%2.fmm", ana::dcaBin[iDca]*1e3, ana::dcaBin[iDca+1]*1e3));

      double binWidth = ana::dcaBin[iDca+1] - ana::dcaBin[iDca];
      hDataD0.SetBinContent(iDca+1, f.GetParameter(0)/binWidth);
      hDataD0.SetBinError(iDca+1, f.GetParError(0)/binWidth);

      // clear pointers
      hTemp->Delete();
   }
   hMCPromptMass->Delete(); 
   hMCPromptMass = nullptr;
   hMCPromptMassAll->Delete();
   hMCPromptMassAll = nullptr;
   hMCNonPromptMass->Delete();
   hMCNonPromptMass = nullptr;
   hMCNonPromptMassAll->Delete();
   hMCNonPromptMassAll = nullptr;
}
