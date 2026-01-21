{
    TCanvas * c = new TCanvas("c", "c", 1000, 600);

    TFile * theFile = new TFile( "pixel_rightOrbit/totalBG_pixel.root");
    //TFile * theFile = new TFile( "base_wrongOrbit/totalBG_DSSD.root");
    theFile->cd();

    TH1D* totalBG_all = (TH1D*)theFile->Get("totalBG_all");
    TH1D* totalBG_P = (TH1D*)theFile->Get("totalBG_P");
    TH1D* totalBG_UC = (TH1D*)theFile->Get("totalBG_UC");
    TH1D* totalBG_TC = (TH1D*)theFile->Get("totalBG_TC");

    totalBG_all->Rebin(10);
    totalBG_all->Scale(0.1);
    totalBG_P->Rebin(10);
    totalBG_P->Scale(0.1);
    totalBG_TC->Rebin(10);
    totalBG_TC->Scale(0.1);
    totalBG_UC->Rebin(10);
    totalBG_UC->Scale(0.1);
    
    totalBG_all->GetXaxis()->SetRangeUser(50, 1e6);

    totalBG_all->SetLineColor(kGray);
    totalBG_all->SetLineWidth(3);
    totalBG_all->DrawCopy("hist l");
    totalBG_P->SetLineColor(kRed);
    totalBG_P->SetLineWidth(3);
    totalBG_P->DrawCopy("hist l same ");
    totalBG_TC->SetLineWidth(3);
    totalBG_UC->SetLineWidth(3);
    totalBG_UC->SetLineColor(kBlue);
    totalBG_TC->SetLineColor(kGreen);
    totalBG_TC->DrawCopy("hist l same ");
    totalBG_UC->DrawCopy("hist l same ");
    c->SetLogy();
    c->SetLogx();

    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.15);
    c->SetTopMargin(0.05);

    TLegend * legend = new TLegend(0.65,0.7,0.95,0.95);
    //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    legend->AddEntry(totalBG_all ,"Total BG","l");
    legend->AddEntry(totalBG_P ,"Fake pair events","l");
    legend->AddEntry(totalBG_TC ,"Fake tracked Compton","l");
    legend->AddEntry(totalBG_UC ,"Fake untracked Compton","l");
    legend->Draw();


    //TFile * theFile = new TFile( "pixel_rightOrbit/totalBG_pixel.root");
    TFile * theFile2 = new TFile( "pixel_EHCut/totalBG_pixel_EHCut.root");
    theFile2->cd();

    theFile2->ls();
    
    TH1D* totalBG_all_CUT = (TH1D*)theFile2->Get("totalBG_all");
    TH1D* totalBG_P_CUT = (TH1D*)theFile2->Get("totalBG_P");
    TH1D* totalBG_UC_CUT = (TH1D*)theFile2->Get("totalBG_UC");
    TH1D* totalBG_TC_CUT = (TH1D*)theFile2->Get("totalBG_TC");

    totalBG_all_CUT->Rebin(10);
    totalBG_all_CUT->Scale(0.1);
    totalBG_P_CUT->Rebin(10);
    totalBG_P_CUT->Scale(0.1);
    totalBG_TC_CUT->Rebin(10);
    totalBG_TC_CUT->Scale(0.1);
    totalBG_UC_CUT->Rebin(10);
    totalBG_UC_CUT->Scale(0.1);
    
    
    totalBG_all_CUT->SetLineStyle(7);
    totalBG_P_CUT->SetLineStyle(7);
    totalBG_UC_CUT->SetLineStyle(7);
    totalBG_TC_CUT->SetLineStyle(7);

    totalBG_all_CUT->SetLineColor(kGray);
    totalBG_all_CUT->SetLineWidth(3);
    totalBG_all_CUT->DrawCopy("hist l same");
    totalBG_P_CUT->SetLineColor(kRed);
    totalBG_P_CUT->SetLineWidth(3);
    totalBG_P_CUT->DrawCopy("hist l same ");
    totalBG_TC_CUT->SetLineWidth(3);
    totalBG_UC_CUT->SetLineWidth(3);
    totalBG_UC_CUT->SetLineColor(kBlue);
    totalBG_TC_CUT->SetLineColor(kGreen);
    totalBG_TC_CUT->DrawCopy("hist l same ");
    totalBG_UC_CUT->DrawCopy("hist l same ");

    
    legend->AddEntry(totalBG_all_CUT ,"after horizon cut","l");
//    legend->AddEntry(totalBG_P_CUT ,"Fake pair events","l");
//    legend->AddEntry(totalBG_TC_CUT ,"Fake tracked Compton","l");
//    legend->AddEntry(totalBG_UC_CUT ,"Fake untracked Compton","l");
//    legend->Draw();

    c->SaveAs("AMEGO-X-Pixel-BG_EHCut.png");

}
