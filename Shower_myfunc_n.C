struct myfunc {
    std::vector<double> par1 = {8.62369, 5.32678, 12.0, 0.696003, 22.8543};
    std::vector<double> par2 = {83.4757, -0.171518, 1.88588, 0.199225, 23.0441};
    std::vector<double> par3 = {0.0227033, 4.97094, 0.772748, 4.0, 0.0018315};
    std::vector<double> par4 = {0.00790974, 1.56139, 1.39, 0.428328, 0.00115601};
    std::vector<double> par5 = {-0.0849149, 2.58396, 0.419646};
    std::vector<double> par6 = {10, 0.0124098, 0.0780432, 0.182612};
    std::vector<double> par7 = {0.144939, -0.435278, 0.0642399};
    Double_t func_xx(Double_t d){
        if (d < 1.7) return par1[0]*TMath::Vavilov(d - par1[1], par1[2], par1[3]) + par1[4];
        else return par2[0]*TMath::Landau(d - par2[1], par2[2], par2[3]) + par2[4];
    }
    Double_t func_a(Double_t d){
        if (d < 1.39) return par3[0]*TMath::Poisson(par3[1]*(d-par3[2]), par3[3])+par3[4];
        else return par4[0]*TMath::Poisson(par4[1]*(d-par4[2]), par4[3])+par4[4];
    }
    Double_t func_h(Double_t d){
        if (d < 1.4 ) return par5[0]*pow(d,par5[1])+par5[2];
        else if ( d < 3.5) return par6[0]*TMath::Landau(d-par6[1],par6[2],par6[3]);
        else return par7[0]*exp(par7[1]*d+par7[2]);
    }
    Double_t m(Double_t d, Double_t angle){
        Double_t h = func_h(d);
        if ((d < 0.8) || (d > 2.8)) return h;
        else{
            Double_t a = func_a(d);
            Double_t xx = func_xx(d);
            a *= a;
            angle = abs(fmod(angle,45.0) - 45*(((int)(angle/45.0))%2));
            if (angle<xx) return a*angle*angle+h;
            else return a*xx*xx+(a/(45/xx-1))*(xx-45)*(xx-45)-(a/(45/xx-1))*(angle-45)*(angle-45)+h;
        }
    }
}func;

int Shower_myfunc_n(){
    Double_t rangex_min(0),rangex_max(8), rangey_min(0), rangey_max(2);
    int nstepx(50), nstepy(200);
    const double L0 = 2.5667;
    
    TCanvas* c1=new TCanvas("PANDA1","c1",800,600);
    gStyle->SetOptTitle(0);
    gStyle->SetStatX(0.36);
    gStyle->SetStatY(0.88);
    gStyle->SetOptStat(0);
    gStyle->SetLabelFont(42,"xyz");
    gStyle->SetLabelSize(0.05,"xyz");
    gStyle->SetLabelOffset(0.01,"xyz");
    gStyle->SetNdivisions(510,"xyz");
    gStyle->SetTitleFont(42,"xyz");
    gStyle->SetTitleColor(1,"xyz");
    gStyle->SetTitleSize(0.05,"xyz");
    gStyle->SetTitleOffset(1.0,"xyz");
    
    TH3D* h = new TH3D("h3D","3D",200,rangex_min,rangex_max,200,rangey_min,rangey_max,200,0,1);
    //TH2D* h = new TH2D("h3D","3D",50,rangex_min,rangex_max,50,rangey_min,rangey_max);
    h->SetMarkerStyle(7);
    h->SetMarkerColorAlpha(kAzure+3, 0.7);
    h->GetXaxis()->SetTitle("distance");
    h->GetYaxis()->SetTitle("d0");
    h->GetZaxis()->SetTitle("E/E_0");
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->GetZaxis()->CenterTitle();
    double Range = 3;
    for (int i =0;i<2000;i++){
        double xx = L0 * ((rand()/double(RAND_MAX)) - 0.5);
        double yy = L0 * ((rand()/double(RAND_MAX)) - 0.5);
        double angle_0 = abs(57.29578*TMath::ATan(yy/xx));
        double d_0 = sqrt(xx*xx+yy*yy);
        double E_0 = func.m(d_0,angle_0);
        
        for (int j =0;j<2000;j++){
            
            double x = Range* ((rand()/double(RAND_MAX)) - 0.5) + xx + L0 + Range/2;
            double y = Range* (rand()/double(RAND_MAX));
            double angle = abs(57.29578*TMath::ATan(y/x));
            double d = sqrt(x*x+y*y);
            double E = func.m(d,angle);
             h->Fill(d,d_0,E/E_0);
            //cout << d << "," << d_0 << "," << E << "," << E_0 << endl;
        }
    }
    
    h->Draw("SCAT");
    return 0;
}
