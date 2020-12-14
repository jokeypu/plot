Double_t par[8] = {1.22,  3.15,  169.0/1405.0, 0.887, 45.0/1405.0, 0.354 ,0.8, 0.8};
void SetPar(Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5, Double_t p6, Double_t p7){
    par[0] = p0;
    par[1] = p1;
    par[2] = p2;
    par[3] = p3;
    par[4] = p4;
    par[5] = p5;
    par[6] = p6;
    par[7] = p7;
}
Double_t FFF(const Double_t *x){
    Double_t r = sqrt(x[0]*x[0]+x[1]*x[1]);
    return (exp(-par[1]*r)+par[2]*exp(-par[3]*r)+par[4]*exp(-par[5]*r))/pow(r,par[7]);
    //return exp(-1/cos(x[0]))+0*x[1];
}
Double_t MMM(Double_t distance, Double_t angle, Double_t p0 = 1.22, Double_t p1 = 3.15, Double_t p2 = 0.12, Double_t p3 = 0.89, Double_t p4 = 0.032, Double_t p5 = 0.354, Double_t p6 = 0.4, Double_t p7 = 0.8){
    //time_t begin,end;
    //begin = clock();
    SetPar(p0, p1, p2, p3, p4, p5 , p6, p7);
    if ( angle > 90 && angle <= 180 ) angle = 180 - angle;
    else if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
    angle *= TMath::DegToRad();
    Double_t L0 = par[0];
    Double_t x0 = distance*cos(angle), y0 = distance*sin(angle);
    Double_t xp = x0+L0 , xm = x0-L0, yp = y0+L0, ym = y0-L0;
    double a[2] = {xm,ym};
    double b[2] = {xp,yp};
    //const double ERRORLIMIT = 1E+2;
    ROOT::Math::Functor wf(&FFF,2);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE,0.0001,0.0001,1000);
    ig.SetFunction(wf);
    Double_t value = ig.Integral(a,b);
    //end = clock();
    //cout << "TIME:" << end - begin << endl;
    return par[6]*value;
}
int Fit_DigiEnergy_test1(){
    TCanvas* c1=new TCanvas("PANDA1","test1",800,600);
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
    
    TGraph2D *g = new TGraph2D("doc/DigiEnergy_R.txt","%lg %lg %lg");
    g->SetMarkerStyle(7);
    g->SetMarkerColorAlpha(kAzure+3, 0.5);
    g->GetZaxis()->SetTitle("E_{digi}");
    g->GetXaxis()->SetTitle("distance");
    g->GetYaxis()->SetTitle("angle");
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    g->GetZaxis()->CenterTitle();
    
    Double_t distance_cut = 10;
    Double_t *x = g->GetX();
    for (int i = 0; i < g->GetN(); i++) if (x[i] > distance_cut) g->RemovePoint(i);
    
    //TF2* f=new TF2("f2","[4]*Shower_Function.shower_Digi(x,y,1.22,3.2,0,[3],0,0.35)",0,distance_cut,0,45);
    TF2* f=new TF2("f2","MMM(x,y,[0], [1], [2], [3], [4], [5], [6], [7])",0,distance_cut,0,45);
    f->SetParameters( 1.22, 3.15, 0.12, 0.89, 0.032, 0.354, 0.4 , 0.8);
    //f->SetParLimits(0, 1.0, 2.0);
    //f->SetParLimits(1, 0.01,0.5);
    //f->SetParLimits(2, 0.001, 0.05);
    //f->SetParLimits(3, 0, 10);
    
    g->GetXaxis()->SetRangeUser(0,distance_cut);
    g->Draw("p.");
    g->Fit(f,"R");
    cout << f->GetParameter(0) << ", " << f->GetParameter(1) << ", " << f->GetParameter(2) << ", " << f->GetParameter(3) << ", " << f->GetParameter(4) << ", " << f->GetParameter(5) << ", " << f->GetParameter(6) << ", "  << f->GetParameter(7) << endl;
    return 0;
}
