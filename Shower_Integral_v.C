#include "Math/IntegratorMultiDim.h"
#include "Math/Functor.h"
const double L0 = 2.5667;
Double_t myfunc(const Double_t * xx);
Double_t func(TVector3 distance);
Double_t ff(Double_t x,Double_t a,Double_t x0,Double_t h);
int Shower_Integral_v(){
    
    Double_t rangex_min(0),rangex_max(3), rangey_min(0), rangey_max(90);
    int nstepx(50), nstepy(200);
        
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
    
    TH3D* h = new TH3D("h3D","3D",200,rangex_min,rangex_max,200,rangey_min,rangey_max,200,0,0.5);
    h->SetMarkerStyle(7);
    h->SetMarkerColorAlpha(kAzure+3, 0.7);
    h->GetXaxis()->SetTitle("distance");
    h->GetYaxis()->SetTitle("angle");
    h->GetZaxis()->SetTitle("E");
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->GetZaxis()->CenterTitle();
    
    Double_t stepx = (rangex_max - rangex_min)/nstepx;
    Double_t stepy = (rangey_max - rangey_min)/nstepy;
    Double_t d = rangex_min;
    for (int i = 0;i<nstepx+1;i++){
        cout << "complete: " << 100*i/nstepx << "%" << endl;
        Double_t angle = rangey_min;
        for (int j = 0;j<nstepy+1;j++){
            TVector2 distance;
            distance.SetMagPhi(d,angle*3.14159/180.0);
            TVector3 distance3(0,distance.X(),distance.Y());
            Double_t E = func(distance3);
            h->Fill(d,angle,E);
            angle+=stepy;
        }
        d+=stepx;
    }
    h->Draw("SCAT");
    return 0;
}

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
    alpha = abs(fmod(alpha,45.0) - 45*(((int)(alpha/45.0))%2));
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

Double_t ff(Double_t x,Double_t a,Double_t x0,Double_t h){
    a *= a;
    if (x<x0) return a*x*x+h;
    else return a*x0*x0+(a/(45/x0-1))*(x0-45)*(x0-45)-(a/(45/x0-1))*(x-45)*(x-45)+h;
}
