/*****************************************************************************************
 *****************************************************************************************
 **                                                                                     **
 ** Questo programma esegue il fit della funzione che descrive la figura di diffrazione **
 ** prodotta da una fenditura lineare. Due funzioni accessorie permettono la visualiz-  **
 ** zazione dei dati sperimentali (utile per stimare i valori iniziali delle grandezze  **
 ** coinvolte) e di produrre il grafico della funzione da fittare.                      **
 **                                                                                     **
 *****************************************************************************************
  C.M. 8 Giugno 2015
  C.M. 4 Febbraio 2016, aggiunta del fondo par[5]
  C.M. 14 Marzo 2024, aggiornamento linea 48, char --> Tstring
 *****************************************************************************************
Si possono utilizzare questi comandi da terminale

> root -l --> Avvia ROOT
root[0] .L Fit_Diffrazione.cpp  --> Carica il file Fit_Diffrazione.cpp
root[1] mydata("nomefile.txt") --> costruisce il grafico dei dati sperimentali
root[2] myfunc() --> visualizza la funzione con valori di default
root[3] myfunc(background, normalizzazione, lambda, larghezza fenditura, shift lungo x, distanza fenditura - schermo) --> visualizza la funzione con valori forniti dall'utente
root[4] myfit("nomefile.txt", background, normalizzazione, lambda, larghezza fenditura, shift lungo x, distanza fenditura - schermo) --> esegue il fit

Tutte le grandezze sono in metri.
 ****************************************************************************************/

/*
PARAMETRI FUNZIONE E FIT
*/
#define BKG 0.
#define I_0 0.6
#define LAMBDA 650.E-9
#define D_PICCOLO 0.15e-3
#define X_0 0.0585
#define LUNGHEZZA 0.954
#define D_GRANDE 1.E-3

Double_t SoloInterferenza(double *x, double *par)
{
   double dx = x[0] - par[1];
   double sin_theta = dx / TMath::Sqrt(dx * dx + par[2] * par[2]);
   double cos_beta = cos(TMath::Pi() * par[6] * sin_theta / par[3]);
   return par[5] + par[4] * cos_beta * cos_beta;
}
Double_t SoloDiffrazione(double *x, double *par)
{
   double dx = x[0] - par[1];
   double sin_theta = dx / TMath::Sqrt(dx * dx + par[2] * par[2]);
   double arg = TMath::Pi() / par[3] * par[0] * sin_theta;
   double sinc_arg = sin(arg) / arg;
   return par[5] + par[4] * sinc_arg * sinc_arg;
}

Double_t Young(double *x, double *par)
{
   double dx = x[0] - par[1];
   double sin_theta = dx / TMath::Sqrt(dx * dx + par[2] * par[2]);
   double arg = TMath::Pi() / par[3] * par[0] * sin_theta;
   double sinc_arg = sin(arg) / arg;
   double cos_beta = cos(TMath::Pi() * par[6] * sin_theta / par[3]);

   return par[5] + par[4] * cos_beta * cos_beta * sinc_arg * sinc_arg;
}

void myfunc(double bkg = BKG, double I0 = I_0, double lambda = LAMBDA, double d = D_PICCOLO, double x0 = X_0, double L = LUNGHEZZA, double D = D_GRANDE)
{
   TF1 *f1 = new TF1("myfunc", Young, x0 - 0.03, x0 + 0.03, 7);
   f1->SetParameter(0, d);
   f1->SetParameter(1, x0);
   f1->SetParameter(2, L);
   f1->SetParameter(3, lambda);
   f1->SetParameter(4, I0);
   f1->SetParameter(5, bkg);
   f1->SetParameter(6, D);
   f1->SetParName(0, "Larghezza fenditura");
   f1->SetParName(1, "shift lungo x");
   f1->SetParName(2, "Distanza fenditura - schermo");
   f1->SetParName(3, "Lambda");
   f1->SetParName(4, "Normalizzazione");
   f1->SetParName(5, "fondo");
   f1->SetParName(6, "distanza tra fenditure");
   f1->Draw("same");
   //f1->Draw();
}

void mydata(TString fname = " ")
{
   TGraphErrors *data = new TGraphErrors(fname, "%lg %lg %lg");
   data->Draw("AP");
   data->SetLineColor(4);
   data->SetMarkerColor(4);
   data->SetTitle("Figura di diffrazione");
   data->GetXaxis()->SetTitle("Posizione [m]");
   data->GetYaxis()->SetTitle("Intensit#grave{a} [u.a.]");
   data->GetXaxis()->CenterTitle(true);
   data->GetXaxis()->CenterTitle(true);
}

void myfit(TString fname = " ", double bkg = BKG, double I0 = I_0, double lambda = LAMBDA, double d = D_PICCOLO, double x0 = X_0, double L = LUNGHEZZA, double D = D_GRANDE)
{
   TGraphErrors *data = new TGraphErrors(fname, "%lg %lg %lg");
   TF1 *f1 = (TF1 *)gROOT->GetFunction("myfunc");
   f1->SetParameter(0, d);
   f1->SetParameter(1, x0);
   f1->SetParameter(2, L);
   f1->SetParameter(3, lambda);
   f1->SetParameter(4, I0);
   f1->SetParameter(5, bkg);
   f1->SetParameter(6, D);

   f1->FixParameter(0, d);
   f1->FixParameter(2, L);

   f1->SetParLimits(3, 600e-9, 700e-9);
   // f1->SetParLimits(1,x0-0.001,x0+0.001);
   // f1->SetParLimits(4,I0-10., I0+10.);
   // f1->SetParLimits(5,bkg-5., bkg+5.);

   data->Fit("myfunc", "R");
   data->Draw("AP");
   data->SetLineColor(9);
   data->SetMarkerStyle(8);
   data->SetMarkerColor(9);
   data->SetMarkerSize(1.5);
   f1->SetLineWidth(1);
   f1->Draw("same");
   data->SetTitle("Figura di diffrazione");
   data->GetXaxis()->SetTitle("Posizione [m]");
   data->GetYaxis()->SetTitle("Intensit#grave{a} [u.a.]");
   data->GetXaxis()->CenterTitle(true);
   data->GetXaxis()->CenterTitle(true);
   TLegend *leg = new TLegend(.6, .7, .9, .9);
   leg->SetTextSize(0.04);
   leg->SetFillColor(0); // fill color is white
   leg->AddEntry(data, "L=0.954m, d=99#mum", "p");
   leg->AddEntry(f1, "fit", "l");
   leg->Draw();
}
