void make_cv_detvar_syst_plots(){

    std::string typeup = "dic";
    std::string typedown = "dic";
    bool isDoubleVar = false;
    float data_pot = 4.89e+19;
    float cv_pot = 1.93572e+20;
    float dirt_pot = 1.68796e+21;
    float dv_up_pot = 1.95567e+20;
    float dv_down_pot = 1.96824e+20;

    float onbeam_triggers = 10905211.; 
    float offbeam_triggers = 77219137.;

    int nbins = 25;
    int bin_low = 0;
    int bin_high = 3;

    TFile *f_cv = new TFile("/uboone/data/users/alister1/numuSelection/files/selected_events/selectedEvents_cv.root");
    TTree* t_cv = (TTree*)f_cv->Get("simulation");

    TFile *f_offbeam = new TFile("/uboone/data/users/alister1/numuSelection/files/selected_events/selectedEvents_offbeam.root");
    TTree* t_offbeam = (TTree*)f_offbeam->Get("offbeam");

    TFile *f_dirt = new TFile("/uboone/data/users/alister1/numuSelection/files/selected_events/selectedEvents_dirt.root");
    TTree* t_dirt = (TTree*)f_dirt->Get("simulation");

    TFile *f_dv_up = new TFile(std::string("/uboone/data/users/alister1/numuSelection/files/selected_events/selectedEvents_"+typeup+".root").c_str());
    TTree* t_dv_up = (TTree*)f_dv_up->Get("simulation");

    TFile *f_dv_down;
    TTree* t_dv_down;

    if (isDoubleVar == true){
        f_dv_down= new TFile(std::string("/uboone/data/users/alister1/numuSelection/files/selected_events/selectedEvents_"+typedown+".root").c_str());
        t_dv_down = (TTree*)f_dv_down->Get("simulation");
    }

    TH1D* h_cv = new TH1D("h_cv", ";Reconstructed Neutrino Energy (GeV);Selected Events", nbins, bin_low, bin_high);
    TH1D* h_dirt = new TH1D("h_dirt", ";Reconstructed Neutrino Energy (GeV);", nbins, bin_low, bin_high);
    TH1D* h_offbeam = new TH1D("h_offbeam", ";Reconstructed Neutrino Energy (GeV);", nbins, bin_low, bin_high);
    TH1D* h_dv_up = new TH1D("h_dv_up", ";Reconstructed Neutrino Energy (GeV;", nbins, bin_low, bin_high);
    TH1D* h_dv_down;
    
    if (isDoubleVar == true){
        h_dv_down = new TH1D("h_dv_down", ";Reconstructed Neutrino Energy (GeV;", nbins, bin_low, bin_high);
    }
 
    t_cv->Draw("reconstructed_neutrino_energy_calib >> h_cv");
    t_offbeam->Draw("reconstructed_neutrino_energy_calib >> h_offbeam");
    t_dirt->Draw("reconstructed_neutrino_energy_calib >> h_dirt");
    t_dv_up->Draw("reconstructed_neutrino_energy_calib >> h_dv_up");

    h_dirt->Sumw2();
    h_dirt->Scale(data_pot/dirt_pot);
    h_offbeam->Sumw2();
    h_offbeam->Scale(onbeam_triggers/offbeam_triggers);
    if (isDoubleVar == true){
        t_dv_down->Draw("reconstructed_neutrino_energy_calib >> h_dv_down");
        h_dv_down->Sumw2();
        h_dv_down->Scale(data_pot/dv_down_pot);
        h_dv_down->Add(h_dirt);
        h_dv_down->Add(h_offbeam);
        h_dv_down->SetLineColor(kRed+1);
        h_dv_down->SetLineWidth(2);
    }

    TH1D* h_cv_clone = (TH1D*)h_cv->Clone("h_cv_clone");
 
    h_cv->Sumw2();
    h_cv_clone->Sumw2();
    h_dv_up->Sumw2();

    h_cv->Scale(data_pot/cv_pot);
    h_cv_clone->Scale(data_pot/cv_pot);
    h_dv_up->Scale(data_pot/dv_up_pot);

    h_cv->Add(h_dirt);
    h_cv->Add(h_offbeam);
    h_cv_clone->Add(h_dirt);
    h_cv_clone->Add(h_offbeam);
    h_dv_up->Add(h_dirt);
    h_dv_up->Add(h_offbeam);

    h_cv->SetLineColor(kBlack);
    h_cv->SetLineWidth(2);
    h_cv_clone->SetFillColor(kGray);
    h_cv_clone->SetLineWidth(0);
    h_dv_up->SetLineColor(kGreen+1);
    h_dv_up->SetLineWidth(2);

    h_cv_clone->GetYaxis()->SetRangeUser(0, h_cv->GetMaximum()*1.3);

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);
    c1->cd();

    h_cv_clone->Draw("E2");
    h_dv_up->Draw("samehist");
    
    if (isDoubleVar == true){
        h_dv_down->Draw("samehist");
    }
    h_cv->DrawClone("samehist");

    TLegend *leg = new TLegend(0.6, 0.70, 0.85, 0.85);
    leg->AddEntry(h_cv, "CV MC");
    if (isDoubleVar == true){
        leg->AddEntry(h_dv_up, "+1 #sigma");
        leg->AddEntry(h_dv_down, "-1 #sigma");
    }
    else
        leg->AddEntry(h_dv_up, "Variation");
    leg->Draw("same");

    std::cout << "CV, number of selected events: " << h_cv->Integral(-1, h_cv->GetNbinsX()+1) << std::endl;
    if (isDoubleVar == true){
        std::cout << "detvar down, number of selected events: " << h_dv_down->Integral(-1, h_dv_down->GetNbinsX()+1) << std::endl;
        std::cout << "detvar up, number of selected events: " << h_dv_up->Integral(-1, h_dv_up->GetNbinsX()+1) << std::endl;

        c1->SaveAs(std::string("CV_detvar_"+typedown+"_"+typeup+".pdf").c_str());
    }
    else {
        std::cout << "detvar, number of selected events: " << h_dv_up->Integral(-1, h_dv_up->GetNbinsX()+1) << std::endl;

        c1->SaveAs(std::string("CV_detvar_"+typeup+".pdf").c_str());
    }
}
