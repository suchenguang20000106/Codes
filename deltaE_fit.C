#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include <iostream>

using namespace RooFit;
using namespace std;


void deltaE()
{
 TTree* tree = (TTree*) gDirectory->Get("selected");
 RooRealVar deltaE("deltaE", "deltaE", -0.1, 0.1);
 RooDataSet data("data", "data of deltaE", tree, deltaE);
 
 RooPlot* Eframe = deltaE.frame();
 data.plotOn(Eframe,Name("data"));

// mframe->Draw(); 

 RooRealVar mean1("mean1", "Mean of deltaE",0., -0.1, 0.1);
 RooRealVar sigma1("sigma1", "sigma of deltaE",0, 0.1);
 RooRealVar mean2("mean2", "Mean of deltaE",-0.03, -0.1, 0.1);
 RooRealVar sigma2("sigma2", "sigma of deltaE",0, 0.1);
 RooRealVar a("a", "a",0.2, 0, 1);
 RooRealVar b("b", "b",0.5, 0, 1);
 //RooRealVar c("c", "c", 0, 1);


 RooRealVar coefficient1("coefficient1", "coe for two guasses", 0, 1);
 RooRealVar coefficient2("coefficient2", "coe for guasses and polynominal", 0, 1);

 RooGaussian gauss1("gauss1", "gauss1(deltaE, mean, sigma)", deltaE, mean1, sigma1);
 RooGaussian gauss2("gauss2", "gauss2(deltaE, mean, sigma)", deltaE, mean2, sigma2);
 RooPolynomial bkg("bkg","back ground",deltaE, RooArgSet(a,b));
// RooGenericPdf sec_poly("sec_poly",  " deltaE * deltaE + a * deltaE + b", RooArgList(c, deltaE, a, b));
 RooAddPdf model("model", "model for deltaE", RooArgList(bkg, gauss1, gauss2), RooArgList(coefficient2, coefficient1));
// RooAddPdf model("model", "model for deltaE", RooArgList( gauss1, gauss2), coefficient1);

 model.fitTo(data);
 model.plotOn(Eframe,Name("model"));
 model.plotOn(Eframe, Components(gauss1),LineStyle(kDashed),LineColor(kRed));
 model.plotOn(Eframe, Components(gauss2),LineStyle(kDashed),LineColor(kGreen));
 model.plotOn(Eframe, Components(bkg),LineStyle(kDashed), LineColor(kOrange)); 

 mean1.Print();
 sigma1.Print();
 mean2.Print();
 sigma2.Print();


 TCanvas* c1 = new TCanvas();

 cout << "chi^2 = " <<  Eframe->chiSquare("model", "data",8) << endl;
 Eframe->GetYaxis()->SetTitle("Events"); 
 Eframe->GetXaxis()->SetTitle("deltaE(GeV)");
 Eframe->SetTitle("Distribution of deltaE");
 Eframe->Draw();

}
