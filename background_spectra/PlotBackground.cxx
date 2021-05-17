void Plot(TString BG);
void PlotBackground();

void PlotBackground() {
  
  vector<TString> theBG = {"_all", "_P", "_TC", "_UC", "_C"};
  for (int i=0; i<theBG.size(); i++) {
    Plot( theBG[i] );
  }

}
void Plot(TString BG){


  TFile *f = new TFile(TString::Format("HadronicBackground.spectrum%s.root", BG.Data()).Data());
  TCanvas * C = (TCanvas*)f->Get("c1");
  TH1D * h_had = (TH1D*)C->GetPrimitive("EnergySpectrum");

  f = new TFile(TString::Format("LeptonicBackground.spectrum%s.root",  BG.Data()).Data());
  C = (TCanvas*)f->Get("c1");
  TH1D * h_lep = (TH1D*)C->GetPrimitive("EnergySpectrum");

  f = new TFile(TString::Format("PhotonicBackground.spectrum%s.root",  BG.Data()).Data());
  C = (TCanvas*)f->Get("c1");
  TH1D * h_pho = (TH1D*)C->GetPrimitive("EnergySpectrum");

  f = new TFile(TString::Format("TrappedHadronicBackground.spectrum%s.root",  BG.Data()).Data());
  C = (TCanvas*)f->Get("c1");
  TH1D * h_thad = (TH1D*)C->GetPrimitive("EnergySpectrum");

  f = new TFile(TString::Format("TrappedLeptonicBackground.spectrum%s.root",  BG.Data()).Data());
  C = (TCanvas*)f->Get("c1");
  TH1D * h_tlep = (TH1D*)C->GetPrimitive("EnergySpectrum");
  
  f = new TFile(TString::Format("HadronicBackgroundDecay.spectrum%s.root",  BG.Data()).Data());
  C = (TCanvas*)f->Get("c1");
  TH1D * h_had_decay = (TH1D*)C->GetPrimitive("EnergySpectrum");
       
  f = new TFile(TString::Format("TrappedHadronicBackgroundDecay.spectrum%s.root",  BG.Data()).Data());
  C = (TCanvas*)f->Get("c1");
  TH1D * h_thad_decay = (TH1D*)C->GetPrimitive("EnergySpectrum");



	TCanvas * C_tot = new TCanvas("C_tot","C_tot", 1181, 802);
	C_tot->SetLeftMargin(0.15);
	C_tot->SetBottomMargin(0.15);

	h_had->SetFillStyle(0);
	h_had->SetLineColor(kRed);
	h_had->SetLineWidth(2);
	h_had->SetTitle("");
	h_had->GetYaxis()->SetTitle("Flux (#gamma/keV/s)");
	h_had->GetXaxis()->SetRangeUser(10, 1e6);
	h_had->GetYaxis()->SetTitleSize(0.06);
	h_had->GetYaxis()->SetTitleOffset(1.0);
	h_had->GetYaxis()->SetLabelSize(0.05);
	h_had->GetXaxis()->SetTitleSize(0.06);
	h_had->GetXaxis()->SetLabelSize(0.05);
	h_had->GetXaxis()->SetTitleOffset(1.1);

	//Convert the x-axis to MeV instead of keV
	//for (int i = 0; i < h_had->GetNBins(); i++) {
	//	h_had->Get


	double T = 3600; //time in simulated seconds
	double scale = 1.0/3600.; //conversion to event rate (per second)

	//Simulations were for 3600 seconds, so divide by time
	h_had->Scale(scale);
	//Multiply by 1000 to conver to /MeV
	//h_had->Scale(1000.);
	h_had->GetYaxis()->SetRangeUser(2e-6, 2);
	h_had->Draw("hist");

	h_lep->SetFillStyle(0);
	h_lep->SetLineColor(kGreen+1);
	h_lep->SetLineWidth(2);
	h_lep->Scale(scale);
	//h_lep->Scale(1000.);
	h_lep->Draw("hist same");

	h_pho->SetFillStyle(0);
	h_pho->SetLineColor(kMagenta);
	h_pho->SetLineWidth(2);
	h_pho->Scale(scale);
	//h_pho->Scale(1000.);
	h_pho->Draw("hist same");

	h_thad->SetFillStyle(0);
	h_thad->SetLineColor(kBlue);
	h_thad->SetLineWidth(2);
	h_thad->Scale(scale);
	//h_thad->Scale(1000.);
	h_thad->Draw("hist same");

	h_tlep->SetFillStyle(0);
	h_tlep->SetLineColor(kOrange);
	h_tlep->SetLineWidth(2);
	h_tlep->Scale(scale);
	//h_tlep->Scale(1000.);
	h_tlep->Draw("hist same");

        h_had_decay->SetFillStyle(0);
        h_had_decay->SetLineColor(kRed+1);
	h_had_decay->SetLineStyle(2);
        h_had_decay->SetLineWidth(2);
        h_had_decay->Scale(scale);
        //h_pho->Scale(1000.);
        h_had_decay->Draw("hist same");

	h_thad_decay->SetFillStyle(0);
	h_thad_decay->SetLineColor(kBlue+1);
	h_thad_decay->SetLineStyle(2);
	h_thad_decay->SetLineWidth(2);
	h_thad_decay->Scale(scale);
	//h_thad_decay->Scale(1000.);
	h_thad_decay->Draw("hist same");


	TH1D * h_tot = (TH1D*)h_had->Clone("h_tot");
	h_tot->Add(h_lep);
	h_tot->Add(h_pho);
	h_tot->Add(h_thad);
	h_tot->Add(h_tlep);
	h_tot->Add(h_had_decay);
	h_tot->Add(h_thad_decay);
	h_tot->SetLineColor(kBlack);
	h_tot->SetLineWidth(2);
	h_tot->Draw("hist same");

	C_tot->SetLogx();
	C_tot->SetLogy();
	gPad->RedrawAxis();

	TLegend * leg = new TLegend(0.65, 0.59, 0.95, 0.95);
	leg->AddEntry(h_had, "Hadrons", "l");
	leg->AddEntry(h_had_decay, "Hadron Activation", "l");
	leg->AddEntry(h_thad, "Trapped Hadrons", "l");
	leg->AddEntry(h_thad_decay, "Trapped Hadron Activation", "l");
	leg->AddEntry(h_lep, "Leptons", "l");
	leg->AddEntry(h_tlep, "Trapped Leptons", "l");
	leg->AddEntry(h_pho, "Photons", "l");
	leg->AddEntry(h_tot, "Total", "l");
	leg->Draw();

	C_tot->SaveAs(TString::Format("bg_pixel%s.png", BG.Data()).Data() );
	C_tot->SaveAs(TString::Format("bg_pixel%s.pdf", BG.Data()).Data() );
	C_tot->SaveAs(TString::Format("bg_pixel%s.root", BG.Data()).Data() );

	h_tot->SetName(TString::Format("totalBG%s", BG.Data() ).Data() );
	h_tot->SaveAs(TString::Format("totalBG_pixel%s.root", BG.Data()).Data() );

};

