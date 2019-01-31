void test_recoE_for_consistency(){

    TFile *f_cv = new TFile("/uboone/data/users/alister1/numuSelection/files/selected_events/selectedEvents_cv.root");
    TTree* t_cv = (TTree*)f_cv->Get("simulation");

    TFile *f_dirt = new TFile("/uboone/data/users/alister1/numuSelection/files/selected_events/selectedEvents_dirt.root");
    TTree* t_dirt = (TTree*)f_dirt->Get("simulation");

    TFile *f_offbeam = new TFile("/uboone/data/users/alister1/numuSelection/files/selected_events/selectedEvents_offbeam.root");
    TTree* t_offbeam = (TTree*)f_offbeam->Get("offbeam");

    TFile *f_onbeam = new TFile("/uboone/data/users/alister1/numuSelection/files/selected_events/selectedEvents_onbeam.root");
    TTree* t_onbeam = (TTree*)f_onbeam->Get("onbeam");


    TH1D* h_mc = new TH1D("h_mc", ";;", 25, 0, 3);
    TH1D* h_dirt = new TH1D("h_dirt", ";;", 25, 0, 3);
    TH1D* h_onbeam = new TH1D("h_onbeam", ";;", 25, 0, 3);
    TH1D* h_offbeam = new TH1D("h_offbeam", ";;", 25, 0, 3);

    t_cv->Draw("reconstructed_neutrino_energy_calib >> h_mc");
    t_dirt->Draw("reconstructed_neutrino_energy_calib >> h_dirt");
    t_onbeam->Draw("reconstructed_neutrino_energy_calib >> h_onbeam");
    t_offbeam->Draw("reconstructed_neutrino_energy_calib >> h_offbeam");

    h_mc->Sumw2();
    h_mc->Scale(0.489/1.96824);

    h_dirt->Sumw2();
    h_dirt->Scale(0.489/0.450955);

    h_offbeam->Sumw2();
    h_offbeam->Scale(10905211./77219137.);

    h_onbeam->GetYaxis()->SetRangeUser(0, h_onbeam->GetMaximum()*1.3);
    h_onbeam->DrawClone("p");
    h_offbeam->DrawClone("histsame");
    h_dirt->DrawClone("histsame");
    h_offbeam->Add(h_mc, 1);
    h_offbeam->Add(h_dirt, 1);
    h_offbeam->Draw("histsame");
    h_onbeam->SetMarkerStyle(20);
    h_onbeam->SetMarkerSize(0.6);
    h_onbeam->Draw("samep");

}
