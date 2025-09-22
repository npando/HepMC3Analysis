#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TAxis.h"
#include "TLine.h"
#include <cstdio>
#include <cstring>



// If the directory exits overwrite it 
static bool LooksSafeDir(const char* d) {
  if (!d || !*d) return false;
  if (!std::strcmp(d,"/"))  return false;
  if (!std::strcmp(d,"."))  return false;
  if (!std::strcmp(d,"..")) return false;
  return true;
}

// Try a few flexible naming patterns for CC histograms (no hadron tag)
static TH1D* GetCC(TFile* f, const char* base, const char* cat) {
  if (!f) return 0;
  char n1[256], n2[256], n3[256];
  // preferred: base_cat
  std::snprintf(n1, sizeof(n1), "%s_%s", base, cat);
  // spaces: base cat
  std::snprintf(n2, sizeof(n2), "%s %s", base, cat);
  // mixed legacy: base_ cat  (underscore then space)
  std::snprintf(n3, sizeof(n3), "%s_ %s", base, cat);

  TH1D* h = (TH1D*)f->Get(n1);
  if (!h) h = (TH1D*)f->Get(n2);
  if (!h) h = (TH1D*)f->Get(n3);
  return h;
}

// Try flexible naming patterns for NC-proton histograms (hadron tag 'p')
static TH1D* GetNCp(TFile* f, const char* base, const char* cat) {
  if (!f) return 0;
  char n1[256], n2[256], n3[256], n4[256], n5[256];
  // preferred: base_p_cat
  std::snprintf(n1, sizeof(n1), "%s_p_%s", base, cat);
  // spaces: base p cat
  std::snprintf(n2, sizeof(n2), "%s %s %s", base, "p", cat);
  // mixed: base_p cat
  std::snprintf(n3, sizeof(n3), "%s_p %s", base, cat);
  // mixed: base p_cat
  std::snprintf(n4, sizeof(n4), "%s %s_%s", base, "p", cat);
  // rare legacy: base _p_ cat (extra underscores)
  std::snprintf(n5, sizeof(n5), "%s _p_ %s", base, cat);

  TH1D* h = (TH1D*)f->Get(n1);
  if (!h) h = (TH1D*)f->Get(n2);
  if (!h) h = (TH1D*)f->Get(n3);
  if (!h) h = (TH1D*)f->Get(n4);
  if (!h) h = (TH1D*)f->Get(n5);
  return h;
}

static void StyleHist(TH1* h, Color_t c, int m=20) {
  if (!h) return;
  h->SetLineColor(c);
  h->SetMarkerColor(c);
  h->SetLineWidth(2);
  h->SetMarkerStyle(m);
  h->SetMarkerSize(0.8);
}

void ratio_cc_nc_minimal(const char* ccfile="cc.root",
                         const char* ncfile="nc.root",
                         const char* outdir="ratios",
                         const char* category="total")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  if (LooksSafeDir(outdir)) {
    // wipe and recreate
    TString cmd = TString::Format("rm -rf -- '%s'", outdir);
    gSystem->Exec(cmd);
  }
  gSystem->mkdir(outdir, kTRUE);

  TFile* fcc = TFile::Open(ccfile);
  if (!fcc || fcc->IsZombie()) { printf("ERROR: cannot open CC file %s\n", ccfile); return; }
  TFile* fnc = TFile::Open(ncfile);
  if (!fnc || fnc->IsZombie()) { printf("ERROR: cannot open NC file %s\n", ncfile); fcc->Close(); delete fcc; return; }

  // Variable bases and x-axis titles (adjust to your exact set)
  struct Fam { const char* base; const char* xtitle; };
  const Fam fams[] = {
    {"hIn_nu_E",        "E_{#nu}^{in} [MeV]"},
    {"hOut_nu_pz",      "p_{z}^{#nu,out} [MeV]"},
    {"hOut_h_cosTheta", "cos#theta(hadron)"},
    {"homega",          "#omega = E_{#nu}^{in} - E_{#nu}^{out} [MeV]"},
    {"hnu_pmag",        "|#vec{p}_{#nu}^{out}| [MeV]"},
    {"hhad_pmag",       "|#vec{p}_{h}| [MeV]"},
    {"hQ2",             "Q^{2} [MeV^{2}]"},
    {"htheta_nu_h",     "#theta(#nu^{out}, h) [rad]"},
    {"hdeltapt",        "#delta p_{T} [MeV]"},
    {"halphat",         "#delta #alpha_{T} [rad]"},
    {"hphit",           "#delta #phi_{T} [rad]"},
    {"hEnu_out",        "E_{#nu}^{out} [MeV]"}
  };
  const int NFAM = sizeof(fams)/sizeof(fams[0]);

  for (int k=0; k<NFAM; ++k) {
    const char* base = fams[k].base;

    // Fetch numerator (CC) and denominator (NC proton)
    TH1D* hCC  = GetCC (fcc, base, category);
    TH1D* hNCp = GetNCp(fnc, base, category);


    if (!hCC || !hNCp) {
      // print a concise diagnostic
      printf("SKIP %-16s : CC=%s  NCp=%s\n",
             base,
             hCC  ? "OK" : "missing",
             hNCp ? "OK" : "missing");
      continue;
    }

    // Check bin compatibility
    if (hCC->GetNbinsX() != hNCp->GetNbinsX()
        || hCC->GetXaxis()->GetXmin() != hNCp->GetXaxis()->GetXmin()
        || hCC->GetXaxis()->GetXmax() != hNCp->GetXaxis()->GetXmax()) {
      printf("MISMATCH bins for %s: cannot form ratio safely.\n", base);
      continue;
    }

    // Prepare ratio: clone CC, divide by NC proton
    if(hCC->Integral()  != 0) hCC->Scale(1.0 / hCC->Integral());
    if(hNCp->Integral() != 0) hNCp->Scale(1.0 / hNCp->Integral());

    TH1D* hRatio = (TH1D*)hCC->Clone(TString::Format("%s_ratio_%s", base, category));
    hRatio->SetDirectory(nullptr);
    hRatio->Sumw2(); // ensure proper errors (in case source lacked Sumw2)


    // Do the division. If you want binomial/Clopper-Pearson style, use "B". For weighted counts, plain divide:
    hRatio->Divide(hNCp); // or hRatio->Divide(hCC, hNCp, 1.0, 1.0, "B");

    // Style and draw
    TCanvas* c = new TCanvas(TString::Format("c_ratio_%s_%s", base, category),
                             TString::Format("ratio %s %s", base, category), 900, 650);
    StyleHist(hRatio, kBlue+2);
    hRatio->GetXaxis()->SetTitle(fams[k].xtitle);
    hRatio->GetYaxis()->SetTitle("CC / NC_{p}");
    hRatio->SetMinimum(0); // tweak if you prefer
    hRatio->Draw("E1");

    // Optional: draw a reference line at 1
    auto line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0,
                          hRatio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kGray+2);
    line->Draw("same");

    // Save
    TString outBase = TString::Format("%s/ratio_%s_%s", outdir, base, category);
    c->SaveAs(outBase + ".png");
    c->SaveAs(outBase + ".pdf");

    delete line;
    delete c;
    delete hRatio;
  }

  fcc->Close(); delete fcc;
  fnc->Close(); delete fnc;
}