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


int main(int argc, char **argv) {
  
    
    ReaderAscii input_file(argv[1]);

    std::string hist1_name = "Proton Multiplicity";
    std::string hist2_name = "Pion Plus Multiplicity";
    std::string hist3_name = "Pion Neg Multiplicity";
    std::string hist4_name = "Pion Neu Multiplicity";
    std::string hist5_name = "PID Counts";
    
    TH1D* proton_hist = new TH1D( hist1_name.c_str(),  hist1_name.c_str() , 100, -1, 10);
    TH1D* pion_plus_hist = new TH1D( hist2_name.c_str(),  hist2_name.c_str() , 100, -1, 10);
    TH1D* pion_neg_hist = new TH1D( hist3_name.c_str(),  hist3_name.c_str() , 100, -1, 10);
    TH1D* pion_neu_hist = new TH1D( hist4_name.c_str(),  hist4_name.c_str() , 100, -1, 10);
    TH1I* pid_hist = new TH1I(hist5_name.c_str(),hist5_name.c_str(),1,0,1);
    pid_hist->SetCanExtend(TH1::kAllAxes);

    while(!input_file.failed()) {
        GenEvent evt(Units::MEV, Units::MM);

        // Read event from input file
        input_file.read_event(evt);

        // Exit if failed to read next event
        if(input_file.failed()) break;

        // Loop over particles in event, keep the outgoing lepton and all protons
        ConstGenParticlePtr lepton_out = nullptr;
        std::vector<ConstGenParticlePtr> proton_out;
        std::vector<ConstGenParticlePtr> pion_plus_out;
        std::vector<ConstGenParticlePtr> pion_neg_out;
        std::vector<ConstGenParticlePtr> pion_neu_out;

        for(const auto &part : evt.particles()) {
            if(part->status() != 1) continue;
            pid_hist->Fill(Form("%d", part->pid()), 1.0); // Fill by label
            if(std::abs(part->pid()) == 13) lepton_out = part;
            else if(part->pid() == 2212) proton_out.push_back(part);
            else if(part->pid() == 211) pion_plus_out.push_back(part);
            else if(part->pid() == -211) pion_neg_out.push_back(part);
            else if(part->pid() == 111) pion_neu_out.push_back(part);

        }
        proton_hist->Fill(proton_out.size());
        pion_plus_hist->Fill(pion_plus_out.size());
        pion_neg_hist->Fill(pion_neg_out.size());
        pion_neu_hist->Fill(pion_neu_out.size());
    
    }

    
    
    input_file.close();

    TFile* outfile = new TFile(argv[2], "RECREATE");
    proton_hist->Write();
    pion_plus_hist->Write();
    pion_neg_hist->Write();
    pion_neu_hist->Write();
    pid_hist->Write();

    delete outfile;

        
    return 0;
}