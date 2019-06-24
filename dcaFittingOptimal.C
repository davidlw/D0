#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TNtuple.h"
#include "TFitter.h"
#include "TFitResult.h"

#include "RooPlot.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooDataHist.h"

#include "myAnaConsts.h"
#include "massfitting.C"

#include <vector>

using namespace RooFit;

TH1D* hD0DcaMCPSignal;
TH1D* hD0DcaMCNPSignal;
TH1D* hD0DcaData;
TH1D* hD0DcaDataPeak;

void dcaFractionFitting(TH1* hData, TObjArray& mc, 
      std::map<std::string, double>& frac, std::map<std::string, std::string>& name, 
      std::string& cut, std::string& label)
{
   TFractionFitter fit(hD0DcaData, &mc);
   fit.Constrain(0, 0.0, 1.0);
   fit.Constrain(1, 0.0, 1.0);
   fit.SetRangeX(1, 13);
   auto fitResultPtr = fit.Fit();
   auto covMat = fitResultPtr->GetCovarianceMatrix();
   covMat.Print();

   TCanvas cDca("cDca", "", 550, 450);
   cDca.Divide(1, 2, 0, 0);
   auto pad1 = cDca.cd(1);
   pad1->SetBorderSize(0);
   pad1->SetLeftMargin(0.16);
   pad1->SetBottomMargin(0.);

   gStyle->SetErrorX(0);

   TH1D* result = (TH1D*) fit.GetPlot();
   //hD0DcaData->GetXaxis()->SetRangeUser(0, 0.036);
   hD0DcaData->GetXaxis()->SetRangeUser(0, 0.044);
   hD0DcaData->SetMarkerStyle(20);
   hD0DcaData->SetTitle(";DCA (cm);Normalized dN / d(D0 Dca) (cm^{-1})");
   hD0DcaData->Draw("E P");
   result->SetLineColor(kBlack);
   result->Draw("SAME");

   double fracNPD0, fracErrNPD0;
   double fracPD0, fracErrPD0;
   fit.GetResult(0, fracNPD0, fracErrNPD0);
   fit.GetResult(1, fracPD0, fracErrPD0);

   TH1D* hMCPredctionNPD0 = (TH1D*) fit.GetMCPrediction(0);
   hMCPredctionNPD0->Scale(fracNPD0);
   hMCPredctionNPD0->SetFillColor(kBlue);
   hMCPredctionNPD0->SetFillStyle(3001);
   TH1D* hMCPredctionPD0 = (TH1D*) fit.GetMCPrediction(1);
   hMCPredctionPD0->Scale(fracPD0);
   hMCPredctionPD0->SetFillColor(kRed);
   hMCPredctionPD0->SetFillStyle(1001);

   THStack hs("hs", "");
   hs.Add(hMCPredctionPD0);
   hs.Add(hMCPredctionNPD0);
   hs.Draw("SAME HIST nostack");
   result->Draw("SAME");
   hD0DcaData->Draw("E P SAME ][");

   TLegend lgd(0.55, 0.7, 0.9, 0.9);
   lgd.AddEntry(hD0DcaData, "Data", "p");
   lgd.AddEntry(hMCPredctionNPD0, "Predicted MC B #rightarrow D", "f");
   lgd.AddEntry(hMCPredctionPD0, "Predicted MC Prompt D", "f");
   lgd.AddEntry(result, "Fitting", "l");
   lgd.Draw();

   TLatex ltx;
   ltx.SetTextSize(0.04);

   ltx.DrawLatexNDC(0.35, 0.8, cut.c_str());
   ltx.DrawLatexNDC(0.6, 0.6, label.c_str());

   double realFracNPD0 = fracNPD0/(fracNPD0+fracPD0);

   double realFracNPD0Err = std::sqrt(std::pow(fracErrNPD0/(fracNPD0+fracPD0)*(1-realFracNPD0), 2) 
         + std::pow(fracErrPD0/(fracNPD0+fracPD0)*realFracNPD0, 2)
         + 2*covMat[0][1]*fracPD0/(fracNPD0+fracPD0)*(-1)*fracNPD0/(fracNPD0+fracPD0));
   std::cout << "error of non-prompt: " << realFracNPD0Err << std::endl;

   frac["fracNPD0"] = realFracNPD0;
   frac["fracNPD0Err"] = realFracNPD0Err;

   ltx.SetTextSize(0.035);
   ltx.DrawLatexNDC(0.55, 0.5, Form("B2D Frac. = %.2f+/-%.2f", realFracNPD0, realFracNPD0Err));

   auto pad2 = cDca.cd(2);
   pad2->SetBorderSize(0);
   pad2->SetLeftMargin(0.16);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.16);

   TH1D* ratio = (TH1D*) hData->Clone();
   ratio->SetName("name");
   ratio->GetYaxis()->SetTitle("dN_{obs}/dN_{th}");
   ratio->Divide(result);
   ratio->Draw();

   cDca.Print(name["noLog"].c_str());
   pad1->SetLogy();
   cDca.Print(name["log"].c_str());

   delete ratio;

   std::cout << "chi2/ndf = " << fit.GetChisquare() << "/ " << fit.GetNDF() << std::endl;
   std::cout << "fit probability: " << fit.GetProb() << std::endl;
}

void dcaFittingOptimal(int mode=2, int method = 0)
//void dcaFittingOptimal(int mode=1, int method = 0)
{
   std::unique_ptr<TFile> f1 = std::unique_ptr<TFile>(new TFile(Form("%s_dca_hists.root", ana::whichtree[mode].c_str())));

   const int nuOfMVA = 7;

   double sig[nuOfMVA] = {0};
   double sigErr[nuOfMVA] = {0};
   double sigPeak[nuOfMVA] = {0};
   double sigPeakErr[nuOfMVA] = {0};
   for(int iMva=0; iMva<nuOfMVA; iMva++)
   {
      int label = 40 + 2*iMva;
      hD0DcaData = (TH1D*) f1->Get(Form("hDcaDataD0mva%d", label));
      hD0DcaDataPeak = (TH1D*) f1->Get(Form("hDcaDataPeakD0mva%d", label));
      hD0DcaMCPSignal = (TH1D*) f1->Get(Form("hDcaMCPD0mva%d", label));
      hD0DcaMCNPSignal = (TH1D*) f1->Get(Form("hDcaMCNPD0mva%d", label));

      for(int iDca=0; iDca<ana::nuofDca; iDca++){
         //hD0DcaMCPSignal->SetBinError(iDca, 0);
         //hD0DcaMCNPSignal->SetBinError(iDca, 0);
      }

      double yield_Data = hD0DcaData->Integral("width");
      hD0DcaData->Scale(1./yield_Data);
      double yield_DataPeak = hD0DcaDataPeak->Integral("width");
      hD0DcaDataPeak->Scale(1./yield_DataPeak);

      TH1D* hMassData = (TH1D*) f1->Get(Form("hMassDataD0mva%d", label));
      TH1D* hMassMCNPD0 = (TH1D*) f1->Get(Form("hMassMCNPD0mva%d", label));
      TH1D* hMassMCNPD0All = (TH1D*) f1->Get(Form("hMassMCNPD0Allmva%d", label));
      TH1D* hMassMCPD0 = (TH1D*) f1->Get(Form("hMassMCPD0mva%d", label));
      TH1D* hMassMCPD0All = (TH1D*) f1->Get(Form("hMassMCPD0Allmva%d", label));

      // create Pic to store pictures
      const std::string dirPic = Form("if [ ! -d \"%sPic\" ]; then\n"
                                 "    mkdir %sPic \n"
                                 "fi", 
                                 ana::whichtree[mode].c_str(), ana::whichtree[mode].c_str());
      gSystem->Exec(dirPic.c_str());

      if(method == 0){
         //create Pic/mva to store pictures
         std::string dirMvaPic = Form( "if [ ! -d \"%sPic/mva%d\" ]; then\n"
                                       "    mkdir %sPic/mva%d \n"
                                       "fi", ana::whichtree[mode].c_str(), iMva, 
                                       ana::whichtree[mode].c_str(), iMva);
         gSystem->Exec(dirMvaPic.c_str());

         TObjArray mc(2);
         mc.Add(hD0DcaMCNPSignal);
         mc.Add(hD0DcaMCPSignal);

         std::map<std::string, double> frac;
         frac["fracNPD0"] = 0;
         frac["fracNPD0Err"] = 0;
         std::map<std::string, std::string> name;
         name["noLog"] = Form("%sPic/mva%d/DcaFitting.png", ana::whichtree[mode].c_str(), iMva);
         name["log"] = Form("%sPic/mva%d/DcaFittingLogScale.png", ana::whichtree[mode].c_str(), iMva);

         std::string cut(TString::Format("MVA > %.2f", 0.4+0.02*iMva));
         std::string label = "4 < pT < 5GeV  |y|<1";

         dcaFractionFitting(hD0DcaData, mc, frac, name, cut, label);

         std::map<std::string, double> fracPeak;
         fracPeak["fracNPD0"] = 0;
         fracPeak["fracNPD0Err"] = 0;
         std::map<std::string, std::string> namePeak;
         namePeak["noLog"] = Form("%sPic/mva%d/PeakDcaFitting.png", ana::whichtree[mode].c_str(), iMva);
         namePeak["log"] = Form("%sPic/mva%d/PeakDcaFittingLogScale.png", ana::whichtree[mode].c_str(), iMva);

         dcaFractionFitting(hD0DcaDataPeak, mc, fracPeak, namePeak, cut, label);

         // fit the invariant mass
         TFitResultPtr fitResultPtr;
         TF1 f = massfitting(hMassData, hMassMCNPD0, hMassMCNPD0All, TString::Format("mva%d", iMva), fitResultPtr);

         // draw the invariant mass
         hMassData->GetYaxis()->SetTitle("Entries /(5 MeV)");
         drawMassFitting(hMassData, f, TString::Format("%sPic/mva%d/MassFitting.png", ana::whichtree[mode].c_str(), iMva), 
               "4 < pT < 5GeV  |y|<1", TString::Format("MVA > %.2f", 0.4+0.02*iMva));

         // calculate the significance
         TF1 signal("signal", "[0]* [5] * (" "[4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))"
               "+ (1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))" ")", 1.7, 2.0);
         for(int ipar=0; ipar<6+1; ipar++){
            signal.FixParameter(ipar, f.GetParameter(ipar));
         }
         signal.ReleaseParameter(0);
         signal.ReleaseParameter(1);
         signal.ReleaseParameter(6);

         TF1 swapAndBkg("swap", "[0]*((1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])))"
               " + " "[9] + [10]*x + [11]*x*x + [12]*x*x*x"
               "+0 *[1]*[2]*[3]*[4]", 1.7, 2.0);
         for(int ipar=0; ipar<12+1; ipar++){
            swapAndBkg.SetParameter(ipar, f.GetParameter(ipar));
         }
         swapAndBkg.FixParameter(2, f.GetParameter(2));
         swapAndBkg.FixParameter(3, f.GetParameter(3));
         swapAndBkg.FixParameter(4, f.GetParameter(4));
         swapAndBkg.FixParameter(5, f.GetParameter(5));
         swapAndBkg.FixParameter(7, f.GetParameter(7));
         swapAndBkg.FixParameter(8, f.GetParameter(8));

         double max = signal.GetMaximum(1.7, 2.0);
         double xmin = signal.GetX(max/2., 1.84, 1.865);
         double xmax = signal.GetX(max/2., 1.865, 1.9);

         double mass_bin_width = hD0DcaData->GetBinWidth(1);
         double yieldSignal = signal.Integral(xmin, xmax) / mass_bin_width;
         double errSignal = signal.IntegralError(xmin, xmax, 0, fitResultPtr->GetCovarianceMatrix().GetMatrixArray()) / mass_bin_width;
         double yieldSwapAndBkg = swapAndBkg.Integral(xmin, xmax) / mass_bin_width;
         double errSwapAndBkg = swapAndBkg.IntegralError(xmin, xmax, 0, fitResultPtr->GetCovarianceMatrix().GetMatrixArray()) / mass_bin_width;

         //std::cout << yieldSignal << "+/-" << errSignal << std::endl;
         //std::cout << yieldSwapAndBkg << "+/-" << errSwapAndBkg << std::endl;

         // sig = frac*s/sqrt(s+b)
         // err_sig = err_frac*s/sqrt(s+b) cross err_s*(frac/sqrt(s+b) - 0.5*frac*s/(s+b)^1.5) cross err_b * 0.5 * frac*s/(s+b)^1.5
         //std::cout << frac["fracNPD0"] << " +/- " << frac["fracNPD0Err"] << std::endl;
         sig[iMva] = signal.Integral(xmin, xmax) * frac["fracNPD0"] / sqrt(f.Integral(xmin, xmax)) / std::sqrt(mass_bin_width);
         sigErr[iMva] = std::sqrt(
                  std::pow(frac["fracNPD0Err"]/frac["fracNPD0"]*sig[iMva] ,2) + 
                  std::pow(sig[iMva]*(errSignal/yieldSignal - 0.5*errSignal/(yieldSignal + yieldSwapAndBkg)), 2) +
                  std::pow(errSwapAndBkg/(yieldSignal + yieldSwapAndBkg) * 0.5 *sig[iMva], 2)
               );
         sigPeak[iMva] = signal.Integral(xmin, xmax) * fracPeak["fracNPD0"] / sqrt(f.Integral(xmin, xmax)) / std::sqrt(mass_bin_width);
         sigPeakErr[iMva] = std::sqrt(
                  std::pow(fracPeak["fracNPD0Err"]/fracPeak["fracNPD0"]*sigPeak[iMva] ,2) + 
                  std::pow(sigPeak[iMva]*(errSignal/yieldSignal - 0.5*errSignal/(yieldSignal + yieldSwapAndBkg)), 2) +
                  std::pow(errSwapAndBkg/(yieldSignal + yieldSwapAndBkg) * 0.5 *sigPeak[iMva], 2)
               );


         TCanvas cSignal("cSignal", "", 550, 450);
         cSignal.SetLogy();
         cSignal.SetLeftMargin(0.16);
         cSignal.SetBottomMargin(0.16);
         hD0DcaDataPeak->Scale(hD0DcaData->GetMaximum()/hD0DcaDataPeak->GetMaximum());
         hD0DcaData->Draw();
         hD0DcaDataPeak->Draw("same");
         hD0DcaData->SetLineColor(kRed);
         hD0DcaDataPeak->SetLineColor(kBlue);
         hD0DcaDataPeak->SetMarkerStyle(kOpenSquare);
         TLegend lgdDca(0.6, 0.75, 0.9, 0.9);
         lgdDca.AddEntry(hD0DcaData, "mass fitting", "lp");
         lgdDca.AddEntry(hD0DcaDataPeak, "sideband", "lp");
         lgdDca.Draw();
         cSignal.Print(Form("%sPic/mva%d/DcaComparison.png", ana::whichtree[mode].c_str(), iMva));
      }
      if(method == 1){
         // the following cannot give correct results
         TCanvas cTemp("cTemp", "", 550, 450);
         RooRealVar x("x", "x", 0, 0.036);
         RooRealVar fracNPD0("fracNPD0", "fracNPD0", 0.5, 0., 1.);
         RooDataHist hNP("hNP", "hNP", RooArgList(x), hD0DcaMCNPSignal);
         RooDataHist hP("hP", "hP", RooArgList(x), hD0DcaMCPSignal);
         RooDataHist hData("hData", "hData", RooArgList(x), hD0DcaData);
         RooHistPdf fNP("fNP", "fNP", x, hNP, 0);
         RooHistPdf fP("fP", "fP", x, hP, 0);
         RooAddPdf sum("sum", "NPD0+PD0", RooArgList(fNP, fP), RooArgList(fracNPD0));
         sum.fitTo(hData, SumW2Error(kTRUE));
         RooArgSet *parList = sum.getParameters(hData);
         parList->Print("v");
         cTemp.cd();
         RooPlot* xframe = x.frame();
         hData.plotOn(xframe);
         sum.plotOn(xframe);
         xframe->Draw();
         cTemp.Print("test.png");
      }
   }

   FILE *fp;
   fp = fopen(Form("%s_sig.txt", ana::whichtree[mode].c_str()), "w");
   for(int iMva=0; iMva<nuOfMVA; iMva++){
      fprintf(fp, "sig_mass_fitting: %lf +/- %lf\n", sig[iMva], sigErr[iMva]);
      fprintf(fp, "sig_side_band:    %lf +/- %lf\n", sigPeak[iMva], sigPeakErr[iMva]);
   }

   // draw the significance plot
   double mva[nuOfMVA];
   for(int iMva=0; iMva<nuOfMVA; iMva++){
      mva[iMva] = 0.4 + ana::mvaStep * iMva;
   }

   TCanvas cSig("cSig", "", 550, 450);

   double mva_x[nuOfMVA];
   double mvaPeak_x[nuOfMVA];
   double mva_ex1[nuOfMVA];
   double mva_ex2[nuOfMVA];
   double mvaPeak_ex1[nuOfMVA];
   double mvaPeak_ex2[nuOfMVA];
   for(int iMva=0; iMva<nuOfMVA; iMva++){
      double mva_x1 = 0.4 + iMva*ana::mvaStep;
      double mva_x2 = 0.4 + (iMva+1)*ana::mvaStep;
      mva_x[iMva] = mva_x1 + 0.5*ana::mvaStep;
      mvaPeak_x[iMva] = mva_x1 + 0.7*ana::mvaStep;
      mva_ex1[iMva] = mva_x[iMva] - mva_x1;
      mva_ex2[iMva] = mva_x2 - mva_x[iMva];
      mvaPeak_ex1[iMva] = mvaPeak_x[iMva] - mva_x1;
      mvaPeak_ex2[iMva] = mva_x2 - mvaPeak_x[iMva];
   }

   TGraphAsymmErrors gSig(nuOfMVA, mva_x, sig, mva_ex1, mva_ex2, sigErr, sigErr);
   TGraphAsymmErrors gSigPeak(nuOfMVA, mvaPeak_x, sigPeak, mvaPeak_ex1, mvaPeak_ex2, sigPeakErr, sigPeakErr);
   gSig.SetMarkerStyle(kFullCircle);
   gSigPeak.SetMarkerStyle(kOpenSquare);
   gSig.SetMarkerColor(kBlue);
   gSigPeak.SetMarkerColor(kRed);
   TMultiGraph multiGraph;
   multiGraph.Add(&gSig);
   multiGraph.Add(&gSigPeak);
   multiGraph.Draw("P SAME");
   auto hist = multiGraph.GetHistogram();
   hist->SetTitle(";MVA;S_{B2D}/#sqrt{S_{D}+B}");
   const double ymin = 0;
   const double ymax = 2.0* (*std::max_element(sig, sig+nuOfMVA));
   hist->GetXaxis()->SetRangeUser(0.38, 0.56);
   hist->GetYaxis()->SetRangeUser(ymin, ymax);
   TLatex ltx;
   ltx.DrawLatexNDC(0.5, 0.2, "4 < p_{T} < 5GeV  |y|<1");
   TLegend lgdSig(0.6, 0.75, 0.9, 0.9);
   lgdSig.AddEntry(&gSig, "mass fitting", "p");
   lgdSig.AddEntry(&gSigPeak, "sideband", "p");
   lgdSig.Draw();
   cSig.Print(TString::Format("%sPic/Sig.png", ana::whichtree[mode].c_str()));

}
