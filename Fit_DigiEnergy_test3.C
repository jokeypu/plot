struct INTEGRAL {
    Double_t y, ci, L0;
    void SetPar(Int_t index, Double_t par){
        if (index == 1) L0 = par;
        else if (index == 2) ci = par;
        else if (index == 3) y = par;
        else std::cout << "Parameter error!!" << std::endl;
    }
    Double_t func_Int(Double_t a, Double_t b){
        float value(0);
        if (y < 0.1){
            Int_t N = 3+100*(b-a);
            float step = (b-a)/N, k = a+step, m = b-step/2;
            for (float i = k; i < m; i+=step) value += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
            return step*((y*(1-exp(-ci*sqrt(a*a+y*y)))/(a*a+y*y)+y*(1-exp(-ci*sqrt(b*b+y*y)))/(b*b+y*y))/2+value);
        }else{
            Int_t n = 4;
            float step = (b-a)/(2*n), Twostep = 2*step, StepOver2 = step/2;
            float sum1(0), sum2(0), k1 = a + step, k2 = k1 + step, m1 = b - StepOver2, m2 = m1 - step;
            for (float i = k1; i < m1; i += Twostep) sum1 += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
            for (float i = k2; i < m2; i += Twostep) sum2 += y*(1-exp(-ci*sqrt(i*i+y*y)))/(i*i+y*y);
            return step/3*(y*(1-exp(-ci*sqrt(a*a+y*y)))/(a*a+y*y)+y*(1-exp(-ci*sqrt(b*b+y*y)))/(b*b+y*y)+4*sum1+2*sum2);
        }
    }
    Double_t line_Int(Double_t x_n, Double_t y_n){
        SetPar(3,y_n);
        if ( x_n-L0 > 0) return func_Int(x_n-L0,x_n+L0);
        else return 2*func_Int(0,L0-x_n)+func_Int(L0-x_n,x_n+L0);
    }
    Double_t block_Int(Double_t distance, Double_t angle, Double_t ci){
        Double_t x0 = distance*cos(angle), y0 = distance*sin(angle);
        Double_t xp = x0+L0 , xm = x0-L0, yp = y0+L0, ym = y0-L0;
        SetPar(2,ci);
        Int_t w1(1),w3(1);
        if (xm>0) { w1 = -1; if (ym>0) w3 = -1; }
        return w1*line_Int(y0,fabs(xm))+line_Int(y0,xp)+w3*line_Int(x0,fabs(ym))+line_Int(x0,yp);
    }
    Double_t shower_Digi(Double_t distance,Double_t angle, Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5){
        if ( angle > 90 && angle <= 180 ) angle = 180 - angle;
        else if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
        angle *= TMath::DegToRad();
        Double_t p[6] = {p0, p1, p2, p3, p4, p5};
        SetPar(1,p[0]);
        return (block_Int(distance,angle,p[1])/p[1]+p[2]*block_Int(distance,angle,p[3])/p[3]+p[4]*block_Int(distance,angle,p[5])/p[5])/3;
    }
    /*Double_t shower_Digi(Double_t distance,Double_t angle, Double_t L_1,Double_t L_2, Double_t L_3, Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5){
        if ( angle > 90 && angle <= 180 ) angle = 180 - angle;
        else if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
        angle *= TMath::DegToRad();
        Double_t p[6] = {p0, p1, p2, p3, p4, p5};
        SetPar(1,L_1);
        Double_t T1 = block_Int(distance,angle,p[1])/p[1];
        SetPar(1,L_2);
        Double_t T2 = p[2]*block_Int(distance,angle,p[3])/p[3];
        SetPar(1,L_3);
        Double_t T3 = p[4]*block_Int(distance,angle,p[5])/p[5];
        return p[0]*(T1+T2+T3)/3;
    }*/
}Shower_Function;

Double_t par[3] = {1.18772, 2.58453, 0.688391};

void SetPar(Double_t p0, Double_t p1, Double_t p2){
    par[0] = p0;
    par[1] = p1;
    par[2] = p2;
}

Double_t FFF(const Double_t *x){
    Double_t r = sqrt(x[0]*x[0]+x[1]*x[1]);
    //return (exp(-par[1]*r)+par[2]*exp(-par[3]*r)+par[4]*exp(-par[5]*r))/pow(r,par[7]);
    return exp(-par[1]*r)/pow(r,par[2]);
    //return exp(-1/cos(x[0]))+0*x[1];
}
Double_t MMM(Double_t distance, Double_t angle, Double_t p0, Double_t p1, Double_t p2){
    //time_t begin,end;
    //begin = clock();
    SetPar(p0,p1,p2);
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
    return value;
}

Double_t fff(Double_t x,Double_t A,Double_t b, Double_t xp, Double_t yp, Double_t k1, Double_t k2){
    //k1 *= (1-exp(-A*pow(x,b)));
    //k2 *= (1-exp(-A*pow(x,b)));
     
    return (TMath::Log10(yp)+k1*xp*TMath::LogE())*exp(-k1*x)+(TMath::Log10(yp)+k2*xp*TMath::LogE())*exp(-k2*x);
}

Double_t FABC(Double_t x,Double_t A, Double_t p1, Double_t p2, Double_t c1, Double_t c2){
    p2 *= (1-exp(-A*pow(x,3)));
    c2 *= (1-exp(-A*pow(x,3)));
    return p1*exp(-p2*x)+c1*exp(-c2*x);
}

int Fit_DigiEnergy_test3(){
    std::string dir_name = "WorkData_1Gamma_A7_E1.0_OR";
    std::string file_name= "doc/WorkData_1Gamma_A7_E1.0_OR_R.txt";
    Int_t NO_Angle = 7;
    Double_t Energy = 1.0;
    //title[20] = "doc/"+ dir_name +"_R.txt"
    ostringstream out1,out2;
    out1 << NO_Angle;
    out2 << fixed << setprecision(1) << Energy;
    string str_NO_Angle = out1.str(), str_Energy = out2.str();
    
    std::ifstream in_file;
    in_file.open(file_name,std::ios::in);
    
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
    
    TGraph *g = new TGraph();
    g->SetMarkerStyle(7);
    g->SetMarkerColorAlpha(kAzure+3, 0.5);
    g->GetYaxis()->SetTitle("E_{digi}");
    g->GetXaxis()->SetTitle("distance");
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    
    string str;
    Int_t N = 0;
    Double_t distance_cut = 10;
    while (std::getline(in_file, str)) {
        std::stringstream strStream(str);
        float distance, angle, energy;
        strStream >> distance >> angle >> energy;
        //if (angle>10 || angle<0) continue;
        if (distance > distance_cut) continue;
        g->SetPoint(N,distance,TMath::Log(energy));
        //g->SetPoint(N,distance,energy);
        N++;
    }
    g->Draw("AP.");
    
    //TF1* f=new TF1("f1","[0]*exp(-[1]*x)+[2]*exp(-[3]*x)",1.2,distance_cut);
    TF1* f=new TF1("f1","TMath::Log(FABC(x,[0],[1],[2],[3],[4]))",0,distance_cut);
    //TF1* f=new TF1("f1","[2]*FABC(x,[0],[1],1,[3],[4],[5])",0,distance_cut);
    f->SetParameters(0.119652, 0.78, 1.8, 0.0198265, 0.088);
    //f->SetParameters( 4, 2, 2, 2,1, 0.2674);
    //f->SetParLimits(4, 2, 3);
    //f->SetParLimits(5, 0.1, 0.5);
    //f->SetParameters( 9.88,2.63, 0.031,0.2674);
    //f->SetParLimits(0, 2, 8);
    //f->SetParLimits(1, 0.1, 3);
    //f->SetParLimits(3, 0.001, 0.1);
    //f->SetParLimits(4, 0.1, 1);
    //f->SetParLimits(5, 0.00001, 0.01);
    //f->SetParLimits(6, 0.0001, 0.1);
    //f->SetParameters();
    //f->SetParLimits(0, 0.8, 1.5);
    g->Fit(f,"R");
    //f->SetParLimits(0, 5, 50);
    //f->SetParLimits(1, 1, 5);
    //f->SetParLimits(2, 0, 1);
    //f->SetParLimits(3, 0.001, 1);
    f->Draw("SAME");
    if (f->GetParameter(2)>f->GetParameter(4))
    cout << f->GetParameter(0) << ", " << f->GetParameter(2) << ", " << (f->GetParameter(3))/(f->GetParameter(1)) << ", " << f->GetParameter(4) << "         ," << f->GetParameter(1) << endl;
    else cout << f->GetParameter(0) << ", " << f->GetParameter(4) << ", " << (f->GetParameter(1))/(f->GetParameter(3)) << ", " << f->GetParameter(2) << "         ," << f->GetParameter(1) << endl;
    
    TF1 *f1=new TF1("f1","[0]*exp(-[1]*x)",0,8);
    f1->SetLineWidth(2);
    f1->SetLineColor(kGreen);
    f1->SetLineStyle(2);
    f1->SetParameters(f->GetParameter(1),f->GetParameter(2));
    f1->SetParLimits(0, 5, 50);
    f1->SetParLimits(1, 1, 5);
    
    TF1 *f2=new TF1("f2","[0]*exp(-[1]*x)",0,8);
    f2->SetLineWidth(2);
    f2->SetLineStyle(2);
    f2->SetLineColor(kCyan);
    f2->SetParameters(f->GetParameter(3),f->GetParameter(4));
    f2->SetParLimits(0, 0, 1);
    f2->SetParLimits(1, 0.001, 1);
    
    //c1->SetLogy();
    f1->Draw("SAME");
    f2->Draw("SAME");

    in_file.close();
    return 0;
}
