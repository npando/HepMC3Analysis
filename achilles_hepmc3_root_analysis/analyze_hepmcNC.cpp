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

// create a marker for what type of hadron emerges as an outgoing particle

enum HadronType { HPROTON = 0, HNEUTRON = 1 };
const char* HADNAME[2] = {"p", "n"};

//create array for each type of interaction mode

const char* CATNAME[3] = {"total","qe","intf"};




int main(int argc, char **argv) {
  if(argc != 3) {
    std::cout << "Usage: " << argv[0] << " <HepMC3_input_file> <Output_file_name>" << std::endl;
    exit(-1);
  }
  
  //Get flux plot for weighting
  TFile* fflux = TFile::Open("fluxes/sbnd_flux.root");
  TH1D* numuFlux_Gev = (TH1D*)fflux->Get("flux_sbnd_numu");
  
  //Set up variables 
  int events_parsed = 0;
  double xsec = 0;
  double sum_weights = 0;
  
  //Read in file
  ReaderAscii input_file(argv[1]);
  
  //Set up arrays of histograms for events with both protons and neutrons
  //Each array corresponds to a variable we're plotting.
  //E.g.:
  //                   hIn_nu_E - Incoming neutrino energy
  //                   hOut_nu_pz - Outgoing neutrino z-momentum
  //etc.
 
  TH1D* hIn_nu_E[2][3];
  TH1D* hOut_nu_pz[2][3];
  TH1D* hOut_h_cosT[2][3];   // hadron cos(theta)
  TH1D* hOmega[2][3];
  TH1D* hNuPmag[2][3];
  TH1D* hHadPmag[2][3];
  TH1D* hQ2[2][3];
  TH1D* hTheta_nu_h[2][3];
  TH1D* hDpt[2][3];
  TH1D* hAlphaT[2][3];
  TH1D* hPhiT[2][3];
  TH1D* hEnuOut[2][3];   
  
  
  //Book 3 histograms for each variable for each hadron: total events, QE (quasi-elastic) events, and inft (interference) events
  
  for (int had = 0; had < 2; ++had)
   {
    for (int ic = 0; ic < 3; ++ic) {
      // helper to build names that are easy to find later, e.g. hQ2_p_total
      auto nm = [&](const char* base){
        return std::string(base) + "_" + HADNAME[had] + "_" + CATNAME[ic];
      };

      hIn_nu_E[had][ic] = new TH1D(nm("hIn_nu_E").c_str(),    "", 100, 0, 2500);
      hOut_nu_pz[had][ic] = new TH1D(nm("hOut_nu_pz").c_str(),  "", 100, 0, 4000);
      hOut_h_cosT[had][ic] = new TH1D(nm("hOut_h_cosTheta").c_str(),"",100,-1,1);
      hOmega[had][ic] = new TH1D(nm("homega").c_str(),      "", 100, 0, 1000);
      hNuPmag[had][ic] = new TH1D(nm("hnu_pmag").c_str(),    "", 100, 0, 2500);
      hHadPmag[had][ic] = new TH1D(nm("hhad_pmag").c_str(),   "", 100, 0, 2500);
      hQ2[had][ic] = new TH1D(nm("hQ2").c_str(),         "", 100, 0, 2.0e6);
      hTheta_nu_h[had][ic] = new TH1D(nm("htheta_nu_h").c_str(), "", 100, 0, 3.5);
      hDpt[had][ic] = new TH1D(nm("hdeltapt").c_str(),    "", 100, 0, 500);
      hAlphaT[had][ic] = new TH1D(nm("halphat").c_str(),     "", 100, 0, 3.5);
      hPhiT[had][ic] = new TH1D(nm("hphit").c_str(),       "", 100, 0, 3.5);
      hEnuOut[had][ic] = new TH1D(nm("hEnu_out").c_str(),    "", 80,  0, 4000);

      // errors on all
      hIn_nu_E[had][ic]->Sumw2();
      hOut_nu_pz[had][ic]->Sumw2();
      hOut_h_cosT[had][ic]->Sumw2();
      hOmega[had][ic]->Sumw2();
      hNuPmag[had][ic]->Sumw2();
      hHadPmag[had][ic]->Sumw2();
      hQ2[had][ic]->Sumw2();
      hTheta_nu_h[had][ic]->Sumw2();
      hDpt[had][ic]->Sumw2();
      hAlphaT[had][ic]->Sumw2();
      hPhiT[had][ic]->Sumw2();
      hEnuOut[had][ic]->Sumw2();

   
  }
   };
  
  
  //Enter loop through events, will stop running when the file "fails" (when you run out of events)
  while(!input_file.failed()) {
    
    // Read event from input file
    GenEvent evt(Units::MEV, Units::MM);
    input_file.read_event(evt);
    
    // Exit if failed to read next event
    if(input_file.failed()) break;
    
    // Get the cross section estimate and the event weight
    std::shared_ptr<GenCrossSection> cs = evt.attribute<GenCrossSection>("GenCrossSection");
    xsec = cs->xsec();
    sum_weights += evt.weights()[0];

    //Read out the details of the first event (good for debugging or understanding the event structure)
    if(events_parsed++==0) {
      std::cout << " First event: " << std::endl;
      Print::listing(evt);
      Print::content(evt);
      std::cout << " GenCrossSection:  ";
      Print::line(cs);
    }

    //Read out the number of events analyzed (good for seeing your code running)
    if(events_parsed%10000 == 0) {
      std::cout << "Events parsed: " << events_parsed << std::endl;
    }

    //Define the particles we'll get variables from
    ConstGenParticlePtr nu_in=nullptr, nu_out=nullptr, p_out=nullptr, n_out=nullptr;

    // Loop through the particles in the event. First look for the incoming particle
    for (const auto& part : evt.particles()) {
      if (part->status() == 4) {
        const auto& ch = part->children();
        if (!ch.empty() && std::abs(ch[0]->pid()) == 14) 
        nu_in = ch[0];
      }
      // Now look for the outgoing particles
      if (part->status() != 1) continue;
      if (std::abs(part->pid()) == 14) nu_out = part;
      else if (part->pid() == 2212)   p_out  = part;
      else if (part->pid() == 2112)   n_out  = part;
    }

    if (!nu_in || !nu_out) continue;

    //Handle the different cases if there is a proton/ neutron

    int had = -1;
    ConstGenParticlePtr h_out = nullptr;
    if (p_out) { had = HPROTON; h_out = p_out; }
    else if (n_out) { had = HNEUTRON; h_out = n_out; }
    else { continue; }


    //Handle the different cases for each interaction mode
    int cat = 0;
    auto procID = evt.attribute<IntAttribute>("signal_process_id");
    if (procID) {
      int pid = procID->value();
      if      (pid >= 200 && pid <= 300) cat = 1; // qe
      else if (pid >= 700 && pid <= 800) cat = 2; // intf
      else                               cat = 0; // falls back to total only
    }


  // extract kinematics (MeV units in file)
    const double in_nu_E     = nu_in->momentum().e();
    const double out_nu_E    = nu_out->momentum().e();
    const double out_nu_pz   = nu_out->momentum().z();
    const double out_nu_pmag = nu_out->momentum().rho();
    const double omega       = in_nu_E - out_nu_E;

    const auto& hP4         = h_out->momentum();
    const double had_cosT    = std::cos(hP4.theta());
    const double had_pmag    = hP4.rho();

    // angles/ TKI variables with TVector3
    TVector3 vnu(nu_out->momentum().px(), nu_out->momentum().py(), nu_out->momentum().pz());
    TVector3 vh (hP4.px(), hP4.py(), hP4.pz());
    const double theta_nu_h  = vnu.Angle(vh);
    TVector3 vnuT(vnu.X(), vnu.Y(), 0), vhT(vh.X(), vh.Y(), 0);

    const double dpt         = (vnuT + vhT).Mag();
    const double phiT        = std::acos( (-vnuT.Dot(vhT)) / (vnuT.Mag()*vhT.Mag()) );
    const TVector3 dpt_vec   = vnuT + vhT;
    const double alphaT      = std::acos( (-vnuT.Dot(dpt_vec)) / (vnuT.Mag()*dpt_vec.Mag()) );

    // NC lepton-side Q2 (uses outgoing ν)
    const auto q = nu_in->momentum() - nu_out->momentum();
    const double Q2_val = -q.m2(); // still in MeV^2

    // fill TOTAL (cat=0) and CATEGORY (cat=1 or 2 if matched)
    const int cats_to_fill[2] = {0, cat};
    const int nfill = (cat == 0 ? 1 : 2);



    
    //Get event weight = fluxWeight*generatorWeight (remembering to convert the neutrino energy to GeV)
    float weight =  numuFlux_Gev->GetBinContent(numuFlux_Gev->FindBin(in_nu_E/1000.))*evt.weights()[0] ;

    const double Phi_tot = numuFlux_Gev->Integral("width");
    //Fill Histograms and weight


    for (int k = 0; k < nfill; ++k) {
      int ic = cats_to_fill[k];
      hIn_nu_E[had][ic]->Fill(in_nu_E,weight);
      hOut_nu_pz[had][ic]->Fill(out_nu_pz,weight);
      hOut_h_cosT[had][ic]->Fill(had_cosT,weight);
      hOmega[had][ic]->Fill(omega,weight);
      hNuPmag[had][ic]->Fill(out_nu_pmag,weight);
      hHadPmag[had][ic]->Fill(had_pmag,weight);
      hQ2[had][ic]->Fill(Q2_val,weight);
      hTheta_nu_h[had][ic]->Fill(theta_nu_h,weight);
      hDpt[had][ic]->Fill(dpt,weight);
      hAlphaT[had][ic]->Fill(alphaT,weight);
      hPhiT[had][ic]->Fill(phiT,weight);
      hEnuOut[had][ic]->Fill(out_nu_E,weight);
    }
  }

    


  
  std::cout << "Normalization factor (xsec/sum_weights) = " << xsec/sum_weights << std::endl;

  //Scale factors applied to the histograms in order to create cross sections
  const double gen_norm = (sum_weights > 0.0) ? (xsec / sum_weights) : 0.0;
  for (int had = 0; had < 2; ++had) {
    for (int ic = 0; ic < 3; ++ic) {
      hIn_nu_E[had][ic]->Scale(gen_norm);
      hOut_nu_pz[had][ic]->Scale(gen_norm);
      hOut_h_cosT[had][ic]->Scale(gen_norm);
      hOmega[had][ic]->Scale(gen_norm);
      hNuPmag[had][ic]->Scale(gen_norm);
      hHadPmag[had][ic]->Scale(gen_norm);
      hQ2[had][ic]->Scale(gen_norm);
      hTheta_nu_h[had][ic]->Scale(gen_norm);
      hDpt[had][ic]->Scale(gen_norm);
      hAlphaT[had][ic]->Scale(gen_norm);
      hPhiT[had][ic]->Scale(gen_norm);
      hEnuOut[had][ic]->Scale(gen_norm);}
    }

  //Close the input file so it doesn't crash
  input_file.close();


  
  //Create an output file to save the histograms
  TFile* outfile = new TFile(argv[2], "RECREATE");
  for (int had = 0; had < 2; ++had) {
    for (int ic = 0; ic < 3; ++ic) {
      hIn_nu_E[had][ic]->Write();
      hOut_nu_pz[had][ic]->Write();
      hOut_h_cosT[had][ic]->Write();
      hOmega[had][ic]->Write();
      hNuPmag[had][ic]->Write();
      hHadPmag[had][ic]->Write();
      hQ2[had][ic]->Write();
      hTheta_nu_h[had][ic]->Write();
      hDpt[had][ic]->Write();
      hAlphaT[had][ic]->Write();
      hPhiT[had][ic]->Write();
      hEnuOut[had][ic]->Write();
    }
  }
  delete outfile;
  
  return 0;
}