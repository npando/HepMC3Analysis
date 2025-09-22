// plot_trios.C
// Draws each trio (total, qe, intf) on a single canvas for every histogram family in your file.

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TAxis.h"
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

namespace {

TH1D* GetHistFlex(TFile* f, const std::string& base, const std::string& cat) {
  if (!f) return nullptr;
  const std::string name_us = base + "_" + cat; // e.g. hQ2_total
  const std::string name_sp = base + " " + cat; // e.g. hQ2 total
  if (auto h = dynamic_cast<TH1D*>(f->Get(name_us.c_str()))) return h;
  if (auto h = dynamic_cast<TH1D*>(f->Get(name_sp.c_str()))) return h;
  return nullptr;
}

void StyleCurve(TH1* h, Color_t c) {
  if (!h) return;
  h->SetLineColor(c);
  h->SetMarkerColor(c);
  h->SetLineWidth(3);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.7);
}

double MaxY(const std::vector<TH1D*>& hs) {
  double m = 0.0;
  for (auto* h : hs) if (h) m = std::max(m, h->GetMaximum());
  return m;
}

} // anon

void plot_trios(const char* infile="my_output.root", const char* outdir="plots")
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(outdir, kTRUE);

  TFile* f = TFile::Open(infile);
  if (!f || f->IsZombie()) { std::cerr << "ERROR: cannot open " << infile << "\n"; return; }

  // Your histogram families:
  struct Fam { std::string base; std::string xtitle; std::string ytitle; bool useErrors; };
  std::vector<Fam> fams = {
    {"hIn_nu_E",        "E_{#nu}^{in} [MeV]",              "Events", false},
    {"hOut_mu_pz",      "p_{z}^{#mu} [MeV]",               "Events", false},
    {"hOut_p_cosTheta", "cos#theta(hadron)",               "Events", false},
    {"homega",          "#omega = E_{#nu}-E_{#mu} [MeV]",  "Events", false},
    {"hmupmag",         "|#vec{p}_{#mu}| [MeV]",           "Events", false},
    {"hppmag",          "|#vec{p}_{h}| [MeV]",             "Events", false},
    {"hQ2",             "Q^{2} [MeV^{2}]",                 "Events", false},
    {"hthetamup",       "#theta(#mu,h) [rad]",             "Events", false},
    {"hdeltapt",        "#delta p_{T} [MeV]",              "Events", false},
    {"halphat",         "#delta #alpha_{T} [rad]",         "Events", false},
    {"hphit",           "#delta #phi_{T} [rad]",           "Events", false},
    {"hEmu",            "E_{#mu}^{out} [MeV]",             "Events", false}
  };

  const std::vector<std::string> cats = {"total","qe","intf"};

  for (const auto& fam : fams) {
    // fetch the three histograms
    std::vector<TH1D*> h(3, nullptr);
    for (int i=0;i<3;i++) {
      h[i] = GetHistFlex(f, fam.base, cats[i]);
    }
    // if none found, skip
    if (!h[0] && !h[1] && !h[2]) {
      std::cerr << "WARNING: none of " << fam.base << " {total,qe,intf} found; skipping.\n";
      continue;
    }

    // choose a frame (first non-null)
    TH1D* frame = h[0] ? h[0] : (h[1] ? h[1] : h[2]);

    // style lines/markers
    StyleCurve(h[0], kBlack);    // total
    StyleCurve(h[1], kRed+1);    // qe
    StyleCurve(h[2], kGreen+2);  // intf

    // canvas
    TString cname = TString::Format("c_%s", fam.base.c_str());
    TCanvas* c = new TCanvas(cname, cname, 900, 650);

    // prepare axes
    double ymax = 1.2 * MaxY(h);
    if (ymax <= 0) ymax = 1.0;

    
    frame->GetXaxis()->SetTitle(fam.xtitle.c_str());
    frame->GetYaxis()->SetTitle(fam.ytitle.c_str());
    frame->SetMinimum(0.0);
    frame->SetMaximum(ymax);

    // draw first
    frame->Draw(fam.useErrors ? "E1" : "hist");
    // overlay the others
    if (h[1]) h[1]->Draw((fam.useErrors ? "E1 same" : "hist same"));
    if (h[2]) h[2]->Draw((fam.useErrors ? "E1 same" : "hist same"));

    // legend
    TLegend* leg = new TLegend(0.65,0.75,0.85,0.88);
    if (h[0]) leg->AddEntry(h[0], "total", "lpe");
    if (h[1]) leg->AddEntry(h[1], "QE",    "lpe");
    if (h[2]) leg->AddEntry(h[2], "Interference", "lpe");
    leg->SetBorderSize(0);
    leg->Draw();

    // save
    TString outBase = TString::Format("%s/%s", outdir, fam.base.c_str());
    c->SaveAs(outBase + ".png");

    delete leg;
    delete c;
  }

  f->Close();
  delete f;
}
