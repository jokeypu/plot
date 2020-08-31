#include "Math/IntegratorMultiDim.h"
#include "Math/Functor.h"
const double L0 = 2.5667;
Double_t myfunc(const Double_t * xx);
Double_t func(TVector3 distance);
int Shower_Integral(){
    TCanvas* c1=new TCanvas("PANDA1","c1",800,600);
    //TCanvas* c2=new TCanvas("PANDA2","c2",800,600);
    gStyle->SetOptTitle(0);
    gStyle->SetStatX(0.36);
    gStyle->SetStatY(0.88);
    gStyle->SetOptStat(1);
    gStyle->SetLabelFont(42,"xyz");
    gStyle->SetLabelSize(0.05,"xyz");
    gStyle->SetLabelOffset(0.01,"xyz");
    gStyle->SetNdivisions(510,"xyz");
    gStyle->SetTitleFont(42,"xyz");
    gStyle->SetTitleColor(1,"xyz");
    gStyle->SetTitleSize(0.05,"xyz");
    gStyle->SetTitleOffset(1.0,"xyz");
    
    TH2D *h = new TH2D("h2d","hist2d",200,0,15,200,0,1.1);
    h->SetMarkerStyle(7);
    h->SetMarkerColorAlpha(kAzure+3, 0.5);
    h->GetXaxis()->SetTitle("distance");
    h->GetYaxis()->SetTitle("E/E0");
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    
    TH1D* h1 = new TH1D("Hist","h1",100,0,1);
    h1->SetLineColor(kBlue);
    h1->SetLineWidth(2);
    h1->GetXaxis()->SetTitle("distance");
    h1->GetYaxis()->SetTitle("Energy");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->CenterTitle();
    
    
    TVector3 test(1,0.0,1);
    cout << func(test) << endl;
    
    vector<double> d;
    for (int i = 0; i < 21; i++){
        d.push_back(-10*L0+i*L0);
    }
    
    for (int k = 0; k < 500; k++){
        double x0 = L0 * ((rand()/double(RAND_MAX)) - 0.5);
        double y0 = L0 * ((rand()/double(RAND_MAX)) - 0.5);
        //double x0 = gRandom->Gaus(0,1);
        //double y0 = gRandom->Gaus(0,1);
        if ((abs(x0) > L0/2) || (abs(y0) > L0/2)) continue;
        TVector3 distance0(x0,0,y0);
        double E0 = -1.0;
        for (int i = 0; i < d.size(); i++){
            for (int j = 0; j < d.size(); j++){
                TVector3 distance(x0+d[i],0,y0+d[j]);
                double E = func(distance);
                if (E>E0) E0 = E;
            }
        }
        for (int i = 0; i < d.size(); i++){
            for (int j = 0; j < d.size(); j++){
                TVector3 distance(x0+d[i],0,y0+d[j]);
                double E = func(distance);
                //if (E<0.002) continue;
                h->Fill(distance.Mag(),E/E0);
            }
        }
        h1->Fill(E0);
    }
    c1->cd();
    h->Draw("SCAT");
    //c2->cd();
    //h1->Draw();
    return 0;
}

/*Double_t myfunc(const Double_t * xx) {
    double x = sqrt(xx[0]*xx[0]+xx[1]*xx[1]);
    Double_t p1(-0.2222), p2(-0.1224), x0(4.439);
    Double_t p0 = 0.03153/0.507112;
    Double_t value;
    if (x <= x0) value = p0/((x-p1) * (x-p1) - p2);
    else{
        Double_t fx0 = p0/((x0-p1) * (x0-p1) - p2);
        Double_t c = fx0 * (x0-p1) * (x0-x) / (p0*1.1512925);
        value = fx0 * pow(10, c);
    }
    return value;
}*/

Double_t myfunc(const Double_t * xx) {
    double x = sqrt(xx[0]*xx[0]+xx[1]*xx[1]);
    Double_t p1(-0.18), p2(-0.12), x0(4.667);
    Double_t p0 = 0.03453/0.507112;
    Double_t value;
    if (x <= x0) value = p0/((x-p1) * (x-p1) - p2);
    else{
        Double_t fx0 = p0/((x0-p1) * (x0-p1) - p2);
        Double_t c = fx0 * (x0-p1) * (x0-x) / (p0*1.1512925);
        value = fx0 * pow(10, c);
    }
    return value;
}

Double_t func(TVector3 distance){
    double L = L0/2;
    TVector3 vz(0.0,0.0,1.0);
    double r = distance.Mag(), alpha = distance.Angle(vz);
    alpha = abs(fmod(alpha,45.0) - ((int)(alpha/45.0))%2);
    double x0 = r*cos(alpha), y0 = r*sin(alpha);
    //if (((x0-L)*(x0+L)<0) && ((y0-L)*(y0+L)<0)) return 1;
    double a[2] = {x0-L,y0-L};
    double b[2] = {x0+L,y0+L};
    const double ERRORLIMIT = 1E-2;
    ROOT::Math::Functor wf(&myfunc,2);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    ig.SetFunction(wf);
    return ig.Integral(a,b);
}
