Double_t f_gamma(float x, float M){
    return (std::sqrt(1.+std::pow(x/M,2)));
}

Double_t f_beta(float x, float M){
    return std::sqrt(1.-1./(std::pow(f_gamma(x,M),2)));
}

Double_t f_beta_gamma(float x, float M){
    return x/M;
}

Double_t f_tmax(float x, float M, float m_ec2){

    return ( 2 * m_ec2 * std::pow(f_beta_gamma(x,M),2) / (1+ 2 * f_gamma(x,M) * m_ec2/M + std::pow(m_ec2/M,2)));

}

void plot_dedx_func_ke(){

    float M = 105;
    float K = 0.307075;
    float z = 1.;
    float Z = 18.;
    float A = 39.948;
    float m_ec2 = 0.511;
    float I = 10.*Z/1000;

    float Kz2ZoverA = K*z*z*Z/A;

    TF1* mean_dedx = new TF1("mean_dedx", 
            "- [1] * (1./std::pow(f_beta(x, [0]),2)) * ( 0.5* std::log( (2*[2]* f_tmax(x,[0],[2])*std::pow(f_beta_gamma(x,[0]),2))/([3]*[3])) - std::pow(f_beta(x, [0]),2))", 0.01, 1000000);
    mean_dedx->SetParameters(M, Kz2ZoverA, m_ec2, I);

    TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
    c1->SetLogx();
    mean_dedx->Draw();
    //mean_dedx->GetYaxis()->SetRangeUser(0,100);

}
