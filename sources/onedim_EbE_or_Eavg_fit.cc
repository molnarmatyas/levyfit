// Why can't you see I am modifying this file???
// Compile: make exe/onedim_EbE_or_Eavg_fit.exe
// Run also from path: exe/onedim_Ebe_or_Eavg_fit.exe **args, e.g. exe/onedim_EbE_or_Eavg_fit.exe 11 "27" 1 10000 1 1
// "path": parent directory of levyfit, analysed, figs directories
// make sure there is a results directory within levyfit dir!
// make sure to move .root input files to "analysed" directory (if not already there)!

#include <Levy_proj_reader.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <complex>
#include <vector>
#include <list>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <sstream>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TNamed.h>
#include <TLegend.h>
#define SQR(x)  ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#include <TLatex.h>

using namespace std;

// Global variables
const double Mass2_pi = 0.019479835;
const int NKT = 10;
const int NFRAME = 3;
const int NCH = 2;
const int NCENT = 10; // number of centrality classes
const char* centleg[NCENT+2] = {"0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-100","all","0-10"};
const int colors[NKT] = {632,416,600,432,616,400,8,9,42,46};
const double ktbins[NKT + 1] = {0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675};
//const double ktbins[NKT + 1] = {0.175,0.200,0.225,0.250,0.275,0.300,0.325,0.350,0.375,0.400,0.425,0.450,0.475,0.500,0.525,0.550,0.575};
const int linestyles[NCENT] = {1,9,7,1};
const int markerstyles[NCENT] = {20,21,34,24};

const float rightmargins[9] = {0.,0.,0.05,0.,0.,0.05,0.,0.,0.05};
const float leftmargins[9]  = {0.2,0.,0.,0.2,0.,0.,0.2,0.,0.};
const float topmargins[9]  = {0.15,0.15,0.15,0.,0.,0.,0.,0.,0.};
const float bottommargins[9]  = {0.15,0.15,0.15,0.,0.,0.,0.2,0.2,0.2};
//const float bottommargins[9]  = {0.,0.,0.,0.,0.,0.,0.2,0.2,0.2};
const char* frames[3] = {"lcms","pcms","lab"};
const char* osl_labels[3] = {"out","side","long"};

// !!! change to actual path
//const char* path = "/home/starelte/OneDrive/KutatÁsás/EPOS4/Dr_fromEPOS";
const char* path = "..";
//const char* path = "/mnt/c/Users/MolnarMatyas/OneDrive-elte.hu/KutatÁsás/EPOS4/Dr_fromEPOS";

TH1* histograms[NKT][NFRAME];
Levy_reader* myLevy_reader;

const double rfitmax = 100.; // TODO make adjustable for kT/cent/... & incr/decr if fit does not converge
const double rfitmin = 5.;

int thiskt = 0;
int thisframe = 0;
int NDF = 0;

// Define the fit function
double fitFunction(const double *x, const double *par)
{
  double alpha = par[0];
  double R = par[1];
  double N = par[2];
  double Rcc = (R*pow(2.,1./alpha));
  return (N/Rcc/Rcc/Rcc)*(myLevy_reader->getValue_3d(alpha, x[0]/Rcc));
  //return (N/(Rcc*Rcc*Rcc))*(myLevy_reader->getValue_3d(alpha, x[0]/Rcc));
}

// The chi-square function to minimize
double chiSquare(const double *params)
{
  NDF = 0;
  double chi2 = 0.0;

  TH1* hist = (TH1F*)histograms[thiskt][thisframe]->Clone();
  for (int bin = 1; bin <= hist->GetNbinsX(); ++bin)
  {
    double x = hist->GetBinCenter(bin);
    double ex = hist->GetBinError(bin);
    if(x > rfitmax) continue;
    if(x < rfitmin) continue;
    double observed = hist->GetBinContent(bin);
    double expected = fitFunction(&x, params);
    if (ex > 0)
    {
      chi2 += pow((observed - expected)/ex, 2.);
    }
    NDF++;
  }
  NDF -= 3;
  return chi2;
}
// The log-likelihood function to minimize
double logLikelihood(const double *params)
{
  NDF = 0;
  double logL = 0.0;
  double integral = histograms[thiskt][thisframe]->Integral();
  for (int ibin = 1; ibin <= histograms[thiskt][thisframe]->GetNbinsX(); ++ibin)
  {
    double x = histograms[thiskt][thisframe]->GetXaxis()->GetBinCenter(ibin);
    if(x > rfitmax) continue;
    if(x < rfitmin) continue;
    double binVolume = ( histograms[thiskt][thisframe]->GetXaxis()->GetBinWidth(ibin) );
      
    double observed = histograms[thiskt][thisframe]->GetBinContent(ibin);
    double expected = fitFunction(&x, params)*binVolume*integral*(x*x*4.*M_PI);
    if(expected <= 0) continue; // Avoid log(0) or negative expected values
    if(observed != 0)
      logL += expected + observed*log(observed/expected) - observed;
    NDF++;
  }
  NDF -= 3;
  return 2.0 * logL;
}

int main(int argc, char *argv[])
{
  // Default args
  int ICENT=-1;
  const char* energy = "200GeV";
  int NFILEMAX = 10;//2;
  int NEVT = 25;//1000;
  int NEVT_AVG = 100; // number of events to be averaged over for each fit
  bool IsUrQMD = false;
  // Args overridden
  if(argc < 2 || argc > 7)
  {
    cerr << "Correct usage: " << endl;
    cerr << "exe/programname.exe ICENT energy NFILEMAX NEVT NEVT_AVG IsUrQMD[bool:true/false]" << endl;
    cerr << "(NEVT is number of events per 'file', NEVT_AVG is number of events over to be averaged)" << endl;
    return 1;
  }
  if(argc > 1)
  {
    ICENT = (int)atoi(argv[1]);
  }
  if(argc > 2)
  {
    energy = argv[2];
  }
  if(argc > 3)
  {
    NFILEMAX = (int)atoi(argv[3]);
  }
  if(argc > 4)
  {
    NEVT = (int)atoi(argv[4]);
  }
  if(argc > 5)
  {
    NEVT_AVG = (int)atoi(argv[5]);
    if(NEVT_AVG < 1)
    {
      cerr << "NEVT_AVG must be at least 1!" << endl;
      return 1;
    }
  }
  if(argc == 7)
  {
    IsUrQMD = static_cast<bool>(atoi(argv[6]));
  }

  // Correct indexing of centrality classes' array - implemented from lcmsonly_epos_pair_source.C macro 
  if(ICENT==NCENT || ICENT < 0)
  {
    ICENT = NCENT; // centleg element: "all"
    cout << "WARNING: all centrality classes together." << endl;
  }else if(ICENT > NCENT) // 0-10% together
  {
    ICENT = NCENT + 1;
    cout << "WARNING: 0-10 percent centrality together." << endl;
  }else{
    cout << centleg[ICENT] << " percent centrality class to be fitted." << endl;
  }

  cout << "about to create levy reader." << endl;
  myLevy_reader = new Levy_reader(Form("%s/levyfit/tables/levy_values_moreprecise.dat",path));
  cout << "levy reader created." << endl;

  TH1* alphahist[NKT];
  TH1* Rhist[NKT];
  TH2* alpha_vs_R[NKT];
  TH2F* alpha_vs_R_all = new TH2F("alpha_vs_R_all","",100,0.6,1.8,100,2.5,11.);
  // Temporary storages
  // For writing the vectors into TGraphAsymmErrors
  Double_t ktbin_centers[NKT], xLow[NKT], xHigh[NKT];//, binWidths[NKT];
  for(int ikt = 0; ikt < NKT; ikt++)
  {
    ktbin_centers[ikt] = (ktbins[ikt] + ktbins[ikt + 1]) / 2.;
    //binWidths[ikt] = ktbins[ikt + 1] - ktbins[ikt];
    double binWidth = ktbins[ikt + 1] - ktbins[ikt];
    xLow[ikt] = binWidth / 2.;
    xHigh[ikt] = binWidth / 2.;
  }
  Double_t alpha_vec[NKT]={0};
  Double_t alpha_errup_vec[NKT]={0};
  Double_t alpha_errdn_vec[NKT]={0};
  Double_t R_vec[NKT]={0};
  Double_t R_errup_vec[NKT]={0};
  Double_t R_errdn_vec[NKT]={0};
  Double_t N_vec[NKT]={0};
  Double_t N_errup_vec[NKT]={0};
  Double_t N_errdn_vec[NKT]={0};
  // For collecting from temporary storages
  int Ngoodfits = 0;
  std::vector<TGraphAsymmErrors*> alpha_vs_KT_all;
  std::vector<TGraphAsymmErrors*> N_vs_KT_all;
  std::vector<TGraphAsymmErrors*> R_vs_KT_all;
  /*
  vector<double> alpha_vec;
  vector<double> alpha_err_vec;
  vector<double> R_vec;
  vector<double> R_err_vec;
  vector<double> N_vec;
  vector<double> N_err_vec;
  */
  //vector<double> chi2_vec;
  for(int ikt = 0; ikt < NKT; ikt++)
  {
    alphahist[ikt] = new TH1F(Form("alphahist_ikt%i",ikt),"",100,0.,2.);
    Rhist[ikt] = new TH1F(Form("Rhist_ikt%i",ikt),"",100,2.5,15.);
    alpha_vs_R[ikt] = new TH2F(Form("alpha_vs_R_ikt%i",ikt),"",100,0.,2.,100,2.5,11.);
  }
  // Open the file and get the histograms
  const char* isPathUrqmd = IsUrQMD ? "UrQMD" : "EPOS";
  TFile *file = TFile::Open(Form("%s/analysed/%s_3d_source_%scent_all_%s.root", path,isPathUrqmd,centleg[ICENT],energy));
  if (!file)
  {
    std::cerr << "Error opening file" << std::endl;
    return 1;
  }

  //TH3F* sourcehists[NKT][NFRAME];
  TCanvas* canvas = new TCanvas("c1", "", 600, 600);//1200, 800); // Compare with Yan

  bool firstplot = true;
  for(int ifile = 1; ifile < NFILEMAX + 1; ifile++) // sorry, the files coming out from the simulation were indexed from 1, I feel too lazy to change everything
  {
    //int firstfile = 0;
    //if(!(ifile%20 == 1 || ifile%20 == 2)) continue;
    for(int ievt = 0; ievt < NEVT; ievt++)
    {
      //int firstevt = 0;
      for (int ikt = 0; ikt < NKT; ++ikt)
      {
        cout << "ifile,ievt,ikt: " << ifile << "," << ievt << "," << ikt << "," << endl;
        TLatex Tl; Tl.SetTextSize(0.04);//Tl.SetTextFont(43); Tl.SetTextSize(15);
        Tl.SetNDC(kTRUE);
    
        for(int iframe = 0; iframe < 3; iframe++)
        {
          thiskt = ikt;
          if(iframe != thisframe) continue;

          TH1* temp_rhohist = nullptr;

          if(IsUrQMD)
          {
            // Add both charges together...
            // pi- pi-
            // -------
            // Form the histogram name
            TString histName = Form("D_%s_ev%d_ch%d_KT%d", frames[iframe], ievt, 0, ikt);
            // Read the histogram
            temp_rhohist = dynamic_cast<TH1D*>(file->Get(histName));
            // Check if histogram was successfully retrieved
            if (!temp_rhohist)
            {
              std::cerr << "Error: Histogram " << histName << " not found in the file." << std::endl;
              continue;
            }

            // pi+ pi+
            // -------
            // Form the histogram name
            histName = Form("D_%s_ev%d_ch%d_KT%d", frames[iframe], ievt, 1, ikt);
            // Read the histogram
            temp_rhohist->Add(dynamic_cast<TH1D*>(file->Get(histName)));
            if (file->Get(histName) == nullptr)
            {
              std::cerr << "Error: Histogram " << histName << " not found in the file." << std::endl;
              continue;
            }
          }
          else
          {
            // Both charges already added...
            // Form the histogram name
            TString histName = Form("pion_pair_source_avg_%s_ifile%i_ievt%i_ikt%i_sqrtrho2", frames[iframe], ifile, ievt, ikt);
            // Read the histogram
            temp_rhohist = dynamic_cast<TH1F*>(file->Get(histName));
          
            // Check if histogram was successfully retrieved
            if (!temp_rhohist)
            {
              std::cerr << "Error: Histogram " << histName << " not found in the file." << std::endl;
              continue;
            }
          }

          //  --- AVERAGING logic ---
          // Only do resetting if first one of (NEVT_AVG) averaged
          if((ievt % NEVT_AVG == 0) || NEVT_AVG == 1)
          {
            if(histograms[ikt][iframe]) histograms[ikt][iframe]->Reset();
            histograms[ikt][iframe] = (TH1F*)temp_rhohist->Clone();
            if(NEVT_AVG != 1) continue; // if no averaging, proceed to fitting right away
          }
          else
          {
            histograms[ikt][iframe]->Add(temp_rhohist); // add for averaging
            if(ievt % NEVT_AVG != NEVT_AVG - 1)
            {
              continue; // do until last one of averaging
            }
          }
          
          // Print some information about the histogram
          int entries = histograms[ikt][iframe]->GetEntries();
          std::cout << "Processed histogram with " << entries << " entries (= pairs, I guess)." << std::endl;
          // TODO based on that, with a global variable, the averaging could stop at a given number of pairs as well (though there will be some fluctuation, as we already can only loop over evts)
          //sourcehist->Scale(1.0 / sourcehist->Integral(1,sourcehist->GetNbinsX()));
          cout << "Integral w overflow after normalization: " << histograms[ikt][iframe]->Integral(0,histograms[ikt][iframe]->GetNbinsX()+1) << endl;
    
          // Create the minimizer
          ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        
          // Set the minimizer properties
          minimizer->SetMaxFunctionCalls(10000);
          minimizer->SetMaxIterations(10000);
          minimizer->SetTolerance(0.001);
        
          // Create the function to be minimized
          ROOT::Math::Functor functor(&logLikelihood, 3);
          //ROOT::Math::Functor functor(&chiSquare, 3);
          minimizer->SetFunction(functor);
        
          // Set the initial parameters and their steps
          minimizer->SetLimitedVariable(0, "alpha", 1.6, 0.01, 0.5, 2.0);
          //minimizer->SetFixedVariable(0, "alpha", 1.3);
          minimizer->SetLimitedVariable(1, "R", 8.5, 0.01, 0., 15.);
          //minimizer->SetFixedVariable(3, "R_long", 20.0);
          //minimizer->SetFixedVariable(4, "N_out", 1.);
          //minimizer->SetFixedVariable(5, "N_side", 1.);
          //minimizer->SetFixedVariable(6, "N_long", 1.);
          minimizer->SetVariable(2, "N", 1., 0.01);
        
          // Minimize the function
          minimizer->Minimize();
          //minimizer->ProvidesError(); // !!!
          //minimizer->Hesse(); // !!!
                  
          // Print the results
          const double *results = minimizer->X();
          const double *errors = minimizer->Errors();
          const double chi2val = minimizer->MinValue();
          std::cout << "alpha: "  << results[0] << std::endl;
          std::cout << "R: "  << results[1] << std::endl;
          std::cout << "N "       << results[2] << std::endl;
          std::cout << "chi2: "   << chi2val << std::endl;

          // Find out hits in a given range vs. number of bins
          int nBins = histograms[ikt][iframe]->GetNbinsX();
          Int_t binsInRange = 0;
          Double_t hitsInRange = 0;
          Double_t xMin=1.; Double_t xMax=100.;
          for(Int_t i=1; i<=nBins; i++)
          {
            Double_t binCenter = histograms[ikt][iframe]->GetBinCenter(i);
            if(binCenter >= xMin && binCenter <= xMax)
            {
              binsInRange++;
              hitsInRange += histograms[ikt][iframe]->GetBinContent(i);
            }
          }
          // Take care of "bad" fits
          cout << "CovMatrixStatus: " << minimizer->CovMatrixStatus() << endl;
          cout << "Status: " << minimizer->Status() << endl;
          double confidence = chi2val / NDF;          
          if(minimizer->CovMatrixStatus() != 3 || minimizer->Status() > 1 || 
          results[0] < 0.55 || results[0] > 1.95 || results[1] < 0.05 || results[1] > 14.95 ||
          (NEVT_AVG==1 && confidence < 0.01) || static_cast<Double_t>(binsInRange) * 0.5 > hitsInRange) // FIXME * 3 maybe too strict?; confidence < 0.01 should not be checked when averaged!!!
          {
            NDF=0; // prbably not needed, but just in case
            cout << "Bad fit, skipping..." << endl;
            delete minimizer;
            continue;
          }
          
          /*
          gStyle->SetLabelSize(0.06,"Y");
          canvas->SetLogx(1);
          canvas->SetLogy(1);
          canvas->SetRightMargin(0.01);
          canvas->SetLeftMargin(0.15);
          canvas->SetTopMargin(0.01);
          canvas->SetBottomMargin(0.15);
          */
          gStyle->SetLabelSize(0.04, "XY"); // Set label size for both axes
          gStyle->SetTitleSize(0.05, "XY"); // Set title size for both axes
          gStyle->SetTitleOffset(1.2, "X"); // Adjust X-axis title offset
          gStyle->SetTitleOffset(1.5, "Y"); // Adjust Y-axis title offset
          
          canvas->SetLogx(1);
          canvas->SetLogy(1);
          canvas->SetRightMargin(0.05);
          canvas->SetLeftMargin(0.15);
          canvas->SetTopMargin(0.10);
          canvas->SetBottomMargin(0.15);

          // Here, among others, creating D(rho) from rho hists (normalising to 1)
          double integral = histograms[ikt][iframe]->Integral();//0,histograms[ikt][iframe][iosl]->GetNbinsX()+1);
          for (int x = 1; x <= histograms[ikt][iframe]->GetNbinsX(); ++x)
          {
            double content = histograms[ikt][iframe]->GetBinContent(x);
            double error = histograms[ikt][iframe]->GetBinError(x);
            double binVolume = histograms[ikt][iframe]->GetXaxis()->GetBinWidth(x);
            histograms[ikt][iframe]->SetBinContent(x, content / binVolume);
            histograms[ikt][iframe]->SetBinError(x, error / binVolume);
          }
          TF1* f_r2 = new TF1("f_r2","1./x/x/4./pi",0.,1.e8);
          histograms[ikt][iframe]->Multiply(f_r2,1.);
          histograms[ikt][iframe]->Scale(1.0 / integral);
          delete f_r2;
          histograms[ikt][iframe]->SetStats(0);
          histograms[ikt][iframe]->SetTitle("");
          //histograms[ikt][iframe]->SetTitle(Form("EPOS4 pion-pion pair source, k_{T} #in [%.2f, %.2f] GeV/c, #sqrt{s_{NN}} = %s, %s%%", ktbins[ikt], ktbins[ikt + 1], energy, centleg[ICENT]));
          histograms[ikt][iframe]->GetXaxis()->SetRangeUser(0.2,5000.);//0.1,1000.); // Compare with Yan
          //histograms[ikt][iframe]->GetXaxis()->SetTitleSize(0.1);
          //histograms[ikt][iframe]->GetYaxis()->SetTitleSize(0.1);
          //histograms[ikt][iframe]->GetYaxis()->SetLabelSize(0.08);
          histograms[ikt][iframe]->GetYaxis()->SetTitle("D(#rho)");
          //histograms[ikt][iframe]->GetYaxis()->SetTitleOffset(0.9);
          histograms[ikt][iframe]->GetXaxis()->SetTitle("#rho [fm]");
          //histograms[ikt][iframe]->GetXaxis()->SetTitleOffset(0.7);
          histograms[ikt][iframe]->SetMinimum(5.e-13); // -12 FIXED
          histograms[ikt][iframe]->SetMaximum(1.e-2);
          
          histograms[ikt][iframe]->GetXaxis()->SetTitleSize(0.05);
          histograms[ikt][iframe]->GetYaxis()->SetTitleSize(0.05);
          histograms[ikt][iframe]->GetXaxis()->SetLabelSize(0.04);
          histograms[ikt][iframe]->GetYaxis()->SetLabelSize(0.04);
          histograms[ikt][iframe]->GetYaxis()->SetTitleOffset(1.5);
          histograms[ikt][iframe]->GetXaxis()->SetTitleOffset(1.2);
          
          TF1* f_levyfunc = new TF1(Form("levyfunc%i",iframe), fitFunction, 0.1, 5000., 3);
          f_levyfunc->SetParNames("alpha","R","norm");
          f_levyfunc->SetLineStyle(1);
          double alpha = results[0];
          double dalpha = errors[0];
          double R = results[1];
          double dR = errors[1];
          double N = results[2];
          double dN = errors[2];

          alphahist[ikt]->Fill(alpha);
          Rhist[ikt]->Fill(R);
          alpha_vs_R[ikt]->Fill(alpha,R);
          alpha_vs_R_all->Fill(alpha,R);

          // Saving to arrays (then vectors):
          alpha_vec[ikt] = alpha;
          alpha_errdn_vec[ikt] = dalpha;
          alpha_errup_vec[ikt] = dalpha;
          //minimizer->GetMinosError(0, alpha_errdn_vec[ikt], alpha_errup_vec[ikt]);
          R_vec[ikt] = R;
          R_errdn_vec[ikt] = dR;
          R_errup_vec[ikt] = dR;
          //minimizer->GetMinosError(1, R_errdn_vec[ikt], R_errup_vec[ikt]);
          N_vec[ikt] = N;
          N_errup_vec[ikt] = dN;
          N_errdn_vec[ikt] = dN;
          //minimizer->GetMinosError(2, N_errdn_vec[ikt], N_errup_vec[ikt]);

          f_levyfunc->SetParameters(alpha,R,N);
          //histograms[ikt][iframe]->GetXaxis()->SetLabelSize(0.08);
          //histograms[ikt][iframe]->GetXaxis()->SetLabelOffset(-0.03);
          histograms[ikt][iframe]->Draw("pe");
          f_levyfunc->SetLineStyle(2);
          f_levyfunc->DrawCopy("same");
          f_levyfunc->SetLineStyle(1);
          f_levyfunc->SetRange(rfitmin,rfitmax);
          //f_levyfunc->SetTitle(Form("EPOS4 pion-pion pair source, k_{T}#in [%.2f, %.2f] GeV/c, #sqrt{s_{NN}}=%s, %s%%", ktbins[ikt],ktbins[ikt+1],energy,centleg[ICENT]));
          f_levyfunc->DrawCopy("same");
          double conflev = TMath::Prob(chi2val, NDF);
          /*
          Tl.SetTextSize(10);
          Tl.DrawLatex(0.58-0.05, 0.78, Form("#chi^{2}/NDF = %1.0f/%i", chi2val, NDF));
          Tl.DrawLatex(0.58-0.05, 0.73, Form("conf.lev. = %1.5f", conflev));
          Tl.DrawLatex(0.58-0.05, 0.68, Form("R = (%1.2f #pm %1.2f) fm",R, dR));
          Tl.DrawLatex(0.58-0.05, 0.63, Form("#alpha = %1.2f #pm %1.2f",alpha, dalpha));
          Tl.DrawLatex(0.58-0.05, 0.58, Form("N = %1.2f #pm %1.2f",N, dN));
          */
          Tl.SetTextSize(0.04);
          Tl.DrawLatexNDC(0.58, 0.83, Form("R = (%.2f #pm %.2f) fm", R, dR));
          Tl.DrawLatexNDC(0.58, 0.78, Form("#alpha = %.2f #pm %.2f", alpha, dalpha));
          Tl.DrawLatexNDC(0.58, 0.73, Form("N = %.2f #pm %.2f", N, dN));
          Tl.DrawLatexNDC(0.58, 0.68, Form("#chi^{2}/NDF = %.1f / %d", chi2val, NDF));
          Tl.DrawLatexNDC(0.58, 0.63, Form("C.L. = %.2f%%", 100 * conflev));
          
          delete minimizer;
          delete f_levyfunc;

          // Add a title using TLatex for better customization
          const char* isPathUrqmdTitle = IsUrQMD ? "UrQMD" : "EPOS4"; // FIXME prbly isPathUrQMD would do as well ("EPOS" instead of "EPOS4", who cares)
          TLatex title;
          title.SetTextAlign(12);  //centered
          title.SetTextSize(0.03);
          title.SetTextFont(42);
          title.SetNDC(true);
          title.DrawLatex(0.04, 0.95, Form("%s pion-pion pair source, k_{T}#in[%.3f, %.3f] GeV/c,#sqrt{s_{NN}} = %s, %s%%",
                                 isPathUrqmdTitle, ktbins[ikt], ktbins[ikt + 1], energy, centleg[ICENT]));


          //Tl.SetTextSize(20);
          //Tl.DrawLatex(0.05, 0.94, Form("EPOS4 200 GeV %s%% AuAu PION PAIR SOURCE, #LTm_{T}#GT = %1.3f, LCMS", centleg[ICENT],sqrt(Mass2_pi+SQR(0.5*(ktbins[ikt]+ktbins[ikt+1]))) ) );
          //if(ievt == 10 && ifile == 1) canvas->SaveAs(Form("%s/figs/fitting/%s/onedsource_cent%s_%s_ifile%i_ievt%i_ikt%i_ich0.png", path, frames[thisframe], centleg[ICENT], energy, ifile, ievt, ikt));
          if(firstplot==true) // && ievt%10 == 0 && ikt == 2 && entries > 100 // !!! FIXME from if(true) 
          {
            canvas->SaveAs(Form("%s/figs/fitting/%s/%s_onedsource_cent%s_%s_ifile%i_ievt%i_ikt%i_ich0.png", path, isPathUrqmd, frames[thisframe], centleg[ICENT], energy, ifile, ievt, ikt));
            firstplot = false;
          }
        } // end of iframe loop
        canvas->Clear();
      } // end of ikt loop
      
      // new TGraphAsymmErrors + put in the vector
      TGraphAsymmErrors* alpha_vs_kt = new TGraphAsymmErrors(NKT, ktbin_centers, alpha_vec, xLow, xHigh, alpha_errdn_vec, alpha_errup_vec);
      alpha_vs_kt->SetTitle(Form("#alpha(K_{T}), #sqrt{s_{NN}}=%s;K_{T} (GeV/c);#alpha",energy));
      alpha_vs_kt->SetName(Form("alpha_vs_kt_%d", Ngoodfits));
      TGraphAsymmErrors* R_vs_kt = new TGraphAsymmErrors(NKT, ktbin_centers, R_vec, xLow, xHigh, R_errdn_vec, R_errup_vec);
      R_vs_kt->SetTitle(Form("#R(K_{T}), #sqrt{s_{NN}}=%s;K_{T} (GeV/c);#R",energy));
      R_vs_kt->SetName(Form("R_vs_kt_%d", Ngoodfits));
      TGraphAsymmErrors* N_vs_kt = new TGraphAsymmErrors(NKT, ktbin_centers, N_vec, xLow, xHigh, N_errdn_vec, N_errup_vec);
      N_vs_kt->SetTitle(Form("#N(K_{T}), #sqrt{s_{NN}}=%s;K_{T} (GeV/c);#N",energy));
      N_vs_kt->SetName(Form("N_vs_kt_%d", Ngoodfits));

      alpha_vs_KT_all.push_back(alpha_vs_kt);
      R_vs_KT_all.push_back(R_vs_kt);
      N_vs_KT_all.push_back(N_vs_kt);
      
      Ngoodfits++;
    } // end of ievt loop
  } // end of ifile loop
  delete canvas;
  delete myLevy_reader;
  cout << "about to close input file." << endl;
  file->Close();
  delete file;
  cout << "input file closed, delete done." << endl;

  TFile* file_output = new TFile(Form("./results/%s_onedfitresults_%s_cent%s_%s.root", isPathUrqmd, frames[thisframe],centleg[ICENT],energy), "RECREATE"); // Compare with Yan, add:  _Yan
  cout << "output file created." << endl;
  file_output->cd();
  cout << "output file cd() done." << endl;

  // Writing out to files
  for(int ikt = 0; ikt < NKT; ikt++)
  {
    cout << "ikt: " << ikt << endl;
    //alphahist[ikt]->Write(); // these two will not be used is alpha_vs_R_all and the TGraphAsymmErrors vectors are written out
    //Rhist[ikt]->Write();
    alpha_vs_R[ikt]->Write();
  }
  alpha_vs_R_all->Write();
  cout << "histograms written." << endl;

  for(size_t i=0; i<alpha_vs_KT_all.size(); i++)
  {
    //cout << "i: " << i << endl;
    alpha_vs_KT_all[i]->Write();
    R_vs_KT_all[i]->Write();
    N_vs_KT_all[i]->Write();
  }
  cout << "vectors in form of TGraphAsymmErrors written." << endl;
  file_output->Close();
  return 0;
}

