#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/Print.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TVector3.h"
#include "TVector3.h"
#include <iostream>

using namespace HepMC3;

// -------------------- Functions --------------------
double compute_Q2(const HepMC3::ConstGenParticlePtr& nu,
                  const HepMC3::ConstGenParticlePtr& mu,
                  const HepMC3::ConstGenParticlePtr& p) {
  // q = k_in - k_out (here: use outgoing muon only, as in your code)
  auto nu_f_vec = nu->momentum();
  auto mu_f_vec = mu->momentum();

  auto q = nu_f_vec - mu_f_vec;
  double Q2 = -q.m2(); // beware units: MeV^2 if inputs are MeV
  return Q2;
}

double compute_opening_angle(const HepMC3::ConstGenParticlePtr& mu,
                             const HepMC3::ConstGenParticlePtr& p) {
  TVector3 vmu(mu->momentum().px(), mu->momentum().py(), mu->momentum().pz());
  TVector3 vp (p  ->momentum().px(), p  ->momentum().py(), p  ->momentum().pz());
  return vmu.Angle(vp); // radians
}

double compute_deltapt(const HepMC3::ConstGenParticlePtr& mu,
                       const HepMC3::ConstGenParticlePtr& p) {
  TVector3 vmuT(mu->momentum().px(), mu->momentum().py(), 0);
  TVector3 vpT (p  ->momentum().px(), p  ->momentum().py(), 0);
  return (vmuT + vpT).Mag();
}

double compute_delta_alphat(const HepMC3::ConstGenParticlePtr& mu,
                            const HepMC3::ConstGenParticlePtr& p) {
  TVector3 vmuT(mu->momentum().px(), mu->momentum().py(), 0);
  TVector3 vpT (p  ->momentum().px(), p  ->momentum().py(), 0);
  TVector3 delta_pT = vmuT + vpT;

  if (vmuT.Mag() == 0 || delta_pT.Mag() == 0) return 0.0;
  double cosval = (-vmuT.Dot(delta_pT)) / (vmuT.Mag() * delta_pT.Mag());
  // numerical safety
  cosval = std::max(-1.0, std::min(1.0, cosval));
  return std::acos(cosval);
}

double compute_delta_phit(const HepMC3::ConstGenParticlePtr& mu,
                          const HepMC3::ConstGenParticlePtr& p) {
  TVector3 vmuT(mu->momentum().px(), mu->momentum().py(), 0);
  TVector3 vpT (p  ->momentum().px(), p  ->momentum().py(), 0);
  if (vmuT.Mag() == 0 || vpT.Mag() == 0) return 0.0;
  double cosval = (-vmuT.Dot(vpT)) / (vmuT.Mag() * vpT.Mag());
  cosval = std::max(-1.0, std::min(1.0, cosval));
  return std::acos(cosval);
}

// -------------------- main --------------------
int main(int argc, char** argv) {
  if (argc != 3) {
    std::cout << "Usage: " << argv[0]
              << " <HepMC3_input_file> <Output_file_name>\n";
    return 1;
  }

  // Flux for weighting
  TFile* fflux = TFile::Open("fluxes/sbnd_flux.root");
  TH1D* numuFlux_Gev = (TH1D*)fflux->Get("flux_sbnd_numu");

  // Counters
  int    events_parsed = 0;
  double xsec = 0.0;
  double sum_weights = 0.0;

  // Reader
  ReaderAscii input_file(argv[1]);

  // Hist arrays
  TH1D* hIn_nu_E[3];
  TH1D* hOut_mu_pz[3];
  TH1D* hOut_p_cosTheta[3];
  TH1D* homega[3];
  TH1D* hmupmag[3];
  TH1D* hppmag[3];
  TH1D* hQ2[3];
  TH1D* hthetamup[3];
  TH1D* hdeltapt[3];
  TH1D* halphat[3];
  TH1D* hphit[3];
  TH1D* hEmu[3];

  std::string inter[3] = {"total","qe","intf"};

  for (int i = 0; i < 3; ++i) {
    // book
    hIn_nu_E[i]       = new TH1D(("hIn_nu_E_"     + inter[i]).c_str(), ("hIn_nu_E_"     + inter[i]).c_str(), 100,    0,  2500);
    hOut_mu_pz[i]     = new TH1D(("hOut_mu_pz_"   + inter[i]).c_str(), ("hOut_mu_pz_"   + inter[i]).c_str(), 100,    0,  4000);
    hOut_p_cosTheta[i]= new TH1D(("hOut_p_cosTheta_"+ inter[i]).c_str(),("hOut_p_cosTheta_"+ inter[i]).c_str(),100,  -1,     1);
    homega[i]         = new TH1D(("homega "       + inter[i]).c_str(), ("homega "       + inter[i]).c_str(), 100,    0,  1000);
    hmupmag[i]        = new TH1D(("hmupmag "      + inter[i]).c_str(), ("hmupmag "      + inter[i]).c_str(), 100,    0,  2500);
    hppmag[i]         = new TH1D(("hppmag "       + inter[i]).c_str(), ("hppmag "       + inter[i]).c_str(), 100,    0,  2500);
    hQ2[i]            = new TH1D(("hQ2 "          + inter[i]).c_str(), ("hQ2 "          + inter[i]).c_str(), 100,    0, 2000e3);
    hthetamup[i]      = new TH1D(("hthetamup "    + inter[i]).c_str(), ("hthetamup "    + inter[i]).c_str(), 100,    0,   3.5);
    hdeltapt[i]       = new TH1D(("hdeltapt "     + inter[i]).c_str(), ("hdeltapt "     + inter[i]).c_str(), 100,    0,   500);
    halphat[i]        = new TH1D(("halphat "      + inter[i]).c_str(), ("halphat "      + inter[i]).c_str(), 100,    0,   3.5);
    hphit[i]          = new TH1D(("hphit "        + inter[i]).c_str(), ("hphit "        + inter[i]).c_str(), 100,    0,   3.5);
    hEmu[i]           = new TH1D(("hEmu_"         + inter[i]).c_str(), "E_{#mu}^{out};E_{#mu}^{out} [MeV];d#sigma/dE_{#mu}", 80, 0.0, 4000);

    // Sumw2 for proper errors
    hIn_nu_E[i]->Sumw2();
    hOut_mu_pz[i]->Sumw2();
    hOut_p_cosTheta[i]->Sumw2();
    homega[i]->Sumw2();
    hmupmag[i]->Sumw2();
    hppmag[i]->Sumw2();
    hQ2[i]->Sumw2();
    hthetamup[i]->Sumw2();
    hdeltapt[i]->Sumw2();
    halphat[i]->Sumw2();
    hphit[i]->Sumw2();
    hEmu[i]->Sumw2();

    hIn_nu_E[i]->SetCanExtend(TH1::kAllAxes);
    hOut_mu_pz[i]->SetCanExtend(TH1::kAllAxes);
    hOut_p_cosTheta[i]->SetCanExtend(TH1::kAllAxes);
    homega[i]->SetCanExtend(TH1::kAllAxes);
    hmupmag[i]->SetCanExtend(TH1::kAllAxes);
    hppmag[i]->SetCanExtend(TH1::kAllAxes);
    hQ2[i]->SetCanExtend(TH1::kAllAxes);
    hthetamup[i]->SetCanExtend(TH1::kAllAxes);
    hdeltapt[i]->SetCanExtend(TH1::kAllAxes);
    halphat[i]->SetCanExtend(TH1::kAllAxes);
    hphit[i]->SetCanExtend(TH1::kAllAxes);
    hEmu[i]->SetCanExtend(TH1::kAllAxes);
  
  }

  // -------------------- event loop --------------------
  while (!input_file.failed()) {
    GenEvent evt(Units::MEV, Units::MM);
    input_file.read_event(evt);
    if (input_file.failed()) break;

    // xsec and event weight
    auto cs = evt.attribute<GenCrossSection>("GenCrossSection");
    if (cs) xsec = cs->xsec();

    if (evt.weights().empty()) {
      // If no generator weight, treat as 1.0
      // (you can choose to skip instead)
      // continue;
    } else {
      sum_weights += evt.weights()[0];
    }

    if (events_parsed++ == 0) {
      std::cout << "First event:\n";
      Print::listing(evt);
      Print::content(evt);
      if (cs) {
        std::cout << "GenCrossSection: ";
        Print::line(cs);
      }
    }
    if (events_parsed % 10000 == 0) {
      std::cout << "Events parsed: " << events_parsed << "\n";
    }

    // collect final-state lepton & protons; incoming neutrino
    ConstGenParticlePtr lepton_out = nullptr;
    ConstGenParticlePtr neutrino_in = nullptr;

    std::vector<ConstGenParticlePtr> protons_out;
    protons_out.reserve(4);

    for (const auto& part : evt.particles()) {
      // beam parent → incoming neutrino is a child (guard empties)
      if (part->status() == 4) {
        const auto& ch = part->children();
        if (!ch.empty()) neutrino_in = ch[0];
      }

      // only final-state
      if (part->status() != 1) continue;

      if (std::abs(part->pid()) == 13) {
        lepton_out = part;
      } else if (part->pid() == 2212) {
        protons_out.push_back(part);
      }
    }

    // we need: neutrino_in, lepton_out, and at least one proton
    if (!neutrino_in || !lepton_out || protons_out.empty()) {
      // Skip proton-based fills if any of these are missing
      continue;
    }

    // pick highest-energy proton (total energy e())
    ConstGenParticlePtr highE_proton = protons_out[0];
    for (size_t j = 1; j < protons_out.size(); ++j) {
      if (protons_out[j]->momentum().e() > highE_proton->momentum().e()) {
        highE_proton = protons_out[j];
      }
    }

    // kinematics
    const double in_nu_E     = neutrino_in->momentum().e();
    const double out_mu_pz   = lepton_out->momentum().z();
    const double out_mu_E    = lepton_out->momentum().e();
    const double out_p_cosTh = std::cos(highE_proton->momentum().theta());
    const double omega       = in_nu_E - out_mu_E;
    const double muon_pmag   = lepton_out->momentum().rho();
    const double proton_pmag = highE_proton->momentum().rho();
    const double Q2_val      = compute_Q2(neutrino_in, lepton_out, highE_proton);
    const double opening_ang = compute_opening_angle(lepton_out, highE_proton);
    const double deltapt_val = compute_deltapt(lepton_out, highE_proton);
    const double deltaphi_val= compute_delta_phit(lepton_out, highE_proton);
    const double dalphaT_val = compute_delta_alphat(lepton_out, highE_proton);
    const double Emu_MeV     = out_mu_E;

    // event weight = flux(Eν in GeV) * generatorWeight
    double gen_w = evt.weights().empty() ? 1.0 : evt.weights()[0];
    // guard flux histogram binning
    double flux_w = 0.0;
    {
      int bin = numuFlux_Gev->FindBin(in_nu_E / 1000.0); // MeV → GeV
      // clamp bin range to avoid 0 or nbins+1 under/overflow if desired
      if (bin < 1) bin = 1;
      if (bin > numuFlux_Gev->GetNbinsX()) bin = numuFlux_Gev->GetNbinsX();
      flux_w = numuFlux_Gev->GetBinContent(bin);
    }
    const double weight = gen_w;

    

    // fill totals
    hIn_nu_E[0]      ->Fill(in_nu_E,     weight);
    hOut_mu_pz[0]    ->Fill(out_mu_pz,   weight);
    hOut_p_cosTheta[0]->Fill(out_p_cosTh,weight);
    homega[0]        ->Fill(omega,       weight);
    hmupmag[0]       ->Fill(muon_pmag,   weight);
    hppmag[0]        ->Fill(proton_pmag, weight);
    hQ2[0]           ->Fill(Q2_val,      weight);
    hthetamup[0]     ->Fill(opening_ang, weight);
    hdeltapt[0]      ->Fill(deltapt_val, weight);
    halphat[0]       ->Fill(dalphaT_val, weight);
    hphit[0]         ->Fill(deltaphi_val,weight);
    hEmu[0]          ->Fill(Emu_MeV,     weight);

    // by interaction type (proc: 1=qe, 2=intf; else ignore)
    int proc = 10;
    if (auto attr = evt.attribute<IntAttribute>("signal_process_id")) {
      int procID = attr->value();
      if      (procID >= 200 && procID <= 300) proc = 1; // qe
      else if (procID >= 700 && procID <= 800) proc = 2; // intf
    }

    if (proc == 1 || proc == 2) {
      hIn_nu_E[proc]      ->Fill(in_nu_E,     weight);
      hOut_mu_pz[proc]    ->Fill(out_mu_pz,   weight);
      hOut_p_cosTheta[proc]->Fill(out_p_cosTh,weight);
      homega[proc]        ->Fill(omega,       weight);
      hmupmag[proc]       ->Fill(muon_pmag,   weight);
      hppmag[proc]        ->Fill(proton_pmag, weight);
      hQ2[proc]           ->Fill(Q2_val,      weight);
      hthetamup[proc]     ->Fill(opening_ang, weight);
      hdeltapt[proc]      ->Fill(deltapt_val, weight);
      halphat[proc]       ->Fill(dalphaT_val, weight);
      hphit[proc]         ->Fill(deltaphi_val,weight);
      hEmu[proc]          ->Fill(Emu_MeV,     weight);
    }
  } // end event loop

  // normalization & output
  if (sum_weights == 0) {
    std::cerr << "WARNING: sum_weights == 0; skipping scaling.\n";
  } else {
    std::cout << "Normalization factor (xsec/sum_weights) = "
              << xsec / sum_weights << "\n";
    const double sf = xsec / sum_weights;
    for (int j = 0; j < 3; ++j) {
      hIn_nu_E[j]      ->Scale(sf);
      hOut_mu_pz[j]    ->Scale(sf);
      hOut_p_cosTheta[j]->Scale(sf);
      homega[j]        ->Scale(sf);
      hmupmag[j]       ->Scale(sf);
      hppmag[j]        ->Scale(sf);
      hQ2[j]           ->Scale(sf);
      hthetamup[j]     ->Scale(sf);
      hdeltapt[j]      ->Scale(sf);
      halphat[j]       ->Scale(sf);
      hphit[j]         ->Scale(sf);
      hEmu[j]          ->Scale(sf);
    }
  }

  // simple integrals 
  {
    double QEerror = 0, Intferror = 0;
    double QEint   = hOut_p_cosTheta[1]->IntegralAndError(1, hOut_p_cosTheta[1]->GetNbinsX(), QEerror);
    double Intfint = hOut_p_cosTheta[2]->IntegralAndError(1, hOut_p_cosTheta[2]->GetNbinsX(), Intferror);
    std::cout << "Integral of QE   = " << QEint   << " +/- " << QEerror   << "\n";
    std::cout << "Integral of Intf = " << Intfint << " +/- " << Intferror << "\n";
  }

  
  {
    std::unique_ptr<TFile> outfile(new TFile(argv[2], "RECREATE"));
    if (!outfile || outfile->IsZombie()) {
      std::cerr << "ERROR: cannot create output file " << argv[2] << "\n";
      return 1;
    }
    for (int n = 0; n < 3; ++n) {
      hIn_nu_E[n]->Write();
      hOut_mu_pz[n]->Write();
      hOut_p_cosTheta[n]->Write();
      homega[n]->Write();
      hmupmag[n]->Write();
      hppmag[n]->Write();
      hQ2[n]->Write();
      hthetamup[n]->Write();
      hdeltapt[n]->Write();
      halphat[n]->Write();
      hphit[n]->Write();
      hEmu[n] ->Write();
    }
    // outfile auto-closed by unique_ptr dtor
  }

  return 0;
}
