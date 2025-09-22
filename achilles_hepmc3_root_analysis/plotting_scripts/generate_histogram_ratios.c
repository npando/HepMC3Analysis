void generate_histogram_ratios(const char* argon_file, const char* carbon_file) {

  TFile *fAr = TFile::Open(argon_file);
  TFile *fC  = TFile::Open(carbon_file);

  auto get = [](TFile* f, const char* name){ return (TH1D*)f->Get(name); };

  TH1D *Ar_tot = get(fAr,"hdeltapt total");
  TH1D *Ar_qe  = get(fAr,"hdeltapt qe");
  TH1D *Ar_int = get(fAr,"hdeltapt intf");

  TH1D *C_tot  = get(fC,"hdeltapt total");
  TH1D *C_qe   = get(fC,"hdeltapt qe");
  TH1D *C_int  = get(fC,"hdeltapt intf");

  // (Optional) assert same binning
  if (!Ar_tot || !C_tot) { std::cout<<"missing hists\n"; return; }

  // Build ratios: Ar/C  (errors propagated via Sumw2 content)
  TH1D *R_tot = (TH1D*)Ar_tot->Clone("R_tot"); R_tot->Divide(C_tot);
  TH1D *R_qe  = (TH1D*)Ar_qe ->Clone("R_qe");  R_qe ->Divide(C_qe );
  TH1D *R_int = (TH1D*)Ar_int->Clone("R_int"); R_int->Divide(C_int);

  // Style
  R_tot->SetLineColor(kBlack);
  R_qe ->SetLineColor(kRed+1);
  R_int->SetLineColor(kGreen+2);
  R_tot->SetLineWidth(2); R_qe->SetLineWidth(2); R_int->SetLineWidth(2);

  // Axes and draw
  TCanvas *c = new TCanvas("c","Ar/C",1000,800);
  R_tot->GetYaxis()->SetTitle("Ratio of Delta PT (Ar/C)");
  R_tot->GetXaxis()->SetTitle("#delta p_{T} [MeV]");
  R_tot->SetMinimum(0.0);
  R_tot->SetMaximum(12.0);  // adjust to taste
  R_tot->Draw("hist");
  R_qe ->Draw("hist same");
  R_int->Draw("hist same");

  auto leg = new TLegend(0.58,0.20,0.88,0.38);
  leg->AddEntry(R_tot,"Total","l");
  leg->AddEntry(R_qe ,"QE","l");
  leg->AddEntry(R_int,"Interference","l");
  leg->Draw();

  c->SaveAs("ratio_Ar_over_C_deltapt.pdf");
}
