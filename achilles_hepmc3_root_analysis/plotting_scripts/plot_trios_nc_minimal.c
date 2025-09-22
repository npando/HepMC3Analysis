// plot_trios_nc_minimal.C
// Draws {total, qe, intf} on one canvas for each (variable base, hadron in {p,n})
// for histograms named base_had_cat. Optionally wipes the output directory.
//
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TAxis.h"
#include <cstdio>
#include <cstring>

// --- safety guard for wiping directories ---
static bool LooksSafeDir(const char* d) {
  if (!d || !*d) return false;
  if (!std::strcmp(d,"/"))  return false;
  if (!std::strcmp(d,"."))  return false;
  if (!std::strcmp(d,"..")) return false;
  return true;
}

// Get hist with several tolerant name patterns
static TH1D* GetHistFlex(TFile* f, const char* base, const char* had, const char* cat) {
  if (!f) return 0;
  char nm1[256], nm2[256], nm3[256];
  std::snprintf(nm1, sizeof(nm1), "%s_%s_%s", base, had, cat); // preferred: base_had_cat
  std::snprintf(nm2, sizeof(nm2), "%s %s %s", base, had, cat); // spaces:   base had cat
  std::snprintf(nm3, sizeof(nm3), "%s_%s %s", base, had, cat); // mixed:    base_had cat
  TH1D* h = (TH1D*)f->Get(nm1);
  if (!h) h = (TH1D*)f->Get(nm2);
  if (!h) h = (TH1D*)f->Get(nm3);
  return h;
}

static void StyleCurve(TH1* h, Color_t c) {
  if (!h) return;
  h->SetLineColor(c);
  h->SetMarkerColor(c);
  h->SetLineWidth(3);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.7);
}

void plot_trios_nc_minimal(const char* infile="my_output.root",
                           const char* outdir="plots",
                           bool wipeOutDir=true)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0); // hide big title box

  if (wipeOutDir && LooksSafeDir(outdir)) {
    TString cmd = TString::Format("rm -rf -- '%s'", outdir);
    gSystem->Exec(cmd);
  }
  gSystem->mkdir(outdir, kTRUE);

  TFile* f = TFile::Open(infile);
  if (!f || f->IsZombie()) { printf("ERROR: cannot open %s\n", infile); return; }

  // Bases and x-axis titles
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

  const char* hads[2] = {"p","n"};                 // proton / neutron
  const char* cats[3] = {"total","qe","intf"};     // three categories

  for (int ih=0; ih<2; ++ih) { // loop over hadron type
    const char* had = hads[ih];

    for (int k=0; k<NFAM; ++k) {
      // fetch trio
      TH1D* h[3] = {0,0,0};
      for (int ic=0; ic<3; ++ic) h[ic] = GetHistFlex(f, fams[k].base, had, cats[ic]);

      if (!h[0] && !h[1] && !h[2]) {
        // not found → skip quietly
        continue;
      }

      // pick frame (first non-null)
      TH1D* frame = h[0] ? h[0] : (h[1] ? h[1] : h[2]);

      // style
      StyleCurve(h[0], kBlack);   // total
      StyleCurve(h[1], kRed+1);   // qe
      StyleCurve(h[2], kGreen+2); // intf

      // canvas name: base_had
      char cname[256]; std::snprintf(cname, sizeof(cname), "c_%s_%s", fams[k].base, had);
      TCanvas* c = new TCanvas(cname, cname, 900, 650);

      // y-range
      double ymax = 0.0;
      for (int ic=0; ic<3; ++ic) if (h[ic] && h[ic]->GetMaximum() > ymax) ymax = h[ic]->GetMaximum();
      if (ymax <= 0) ymax = 1.0;
      ymax *= 1.2;

      frame->GetXaxis()->SetTitle(fams[k].xtitle);
      frame->GetYaxis()->SetTitle("Events (AU)"); // adjust if you flux-average / differentialize
      frame->SetMinimum(0.0);
      frame->SetMaximum(ymax);

      frame->Draw("E1");
      if (h[1]) h[1]->Draw("E1 same");
      if (h[2]) h[2]->Draw("E1 same");

      // small transparent legend
      TLegend* leg = new TLegend(0.65,0.80,0.85,0.90);
      leg->SetTextSize(0.03);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      if (h[0]) leg->AddEntry(h[0], "total", "lpe");
      if (h[1]) leg->AddEntry(h[1], "QE",    "lpe");
      if (h[2]) leg->AddEntry(h[2], "Interference", "lpe");
      leg->Draw();

      // save: outdir/base_had.{png,pdf}
      TString outBase = TString::Format("%s/%s_%s", outdir, fams[k].base, had);
      c->SaveAs(outBase + ".png");

      delete leg;
      delete c;
    }
  }

  f->Close();
  delete f;
}