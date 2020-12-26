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

Double_t FABC(Double_t x,Double_t A, Double_t p1, Double_t p2, Double_t c1, Double_t c2){
    p2 *= (1-exp(-A*pow(x,3)));
    c2 *= (1-exp(-A*pow(x,3)));
    return p1*exp(-p2*x)+c1*exp(-c2*x);
}

int Fit_DigiEnergy(std::string dir_name, const char title[20], Int_t NO_Angle, Double_t Energy){
    //title[20] = "doc/"+ dir_name +"_R.txt"
    ostringstream out1,out2;
    out1 << NO_Angle;
    out2 << fixed << setprecision(1) << Energy;
    string str_NO_Angle = out1.str(), str_Energy = out2.str();
    std::string out_name = "doc/A"+str_NO_Angle+"_par.txt";
    std::ofstream par_file;
    par_file.open(out_name,std::ios::app);
    
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
    
    std::string file_name(title);
    //std::string file_name = "doc/WorkData_1Gamma_A7_E1.0_OR_R.txt";
    std::ifstream in_file;
    in_file.open(file_name,std::ios::in);
    
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
        g->SetPoint(N,distance,pow(energy,1.0/3.0));
        //g->SetPoint(N,distance,energy);
        N++;
    }
    
    TF1* f=new TF1("f1","pow(FABC(x,[0],[1],[2],[3],[4]),1.0/3.0)",0,distance_cut);
    //TF1* f=new TF1("f1","[2]*FABC(x,[0],[1],1,[3],[4],[5])",0,distance_cut);
    f->SetParameters(0.18, 0.8, 1.6, 0.0087, 0.055);
    g->Draw("AP.");
    g->Fit(f,"R");
    
    if (f->GetParameter(2)>f->GetParameter(4))
    par_file << str_Energy << " " << f->GetParameter(0) << " " << f->GetParameter(2) << " " << (f->GetParameter(3))/(f->GetParameter(1)) << " " << f->GetParameter(4) << endl;
    else par_file << str_Energy << " " <<f->GetParameter(0) << " " << f->GetParameter(4) << ", " << (f->GetParameter(1))/(f->GetParameter(3)) << " " << f->GetParameter(2) << endl;

    //TF1* ff=new TF1("ff","FABC(x,[0],[1],[2],[3],[4])",0,distance_cut);
    //ff->SetParameters(f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3),f->GetParameter(4));
    //g->Draw("AP.");
    //ff->Draw("SAME");
    c1->SetLogy();
    TString picture_name= "doc/A"+str_NO_Angle+"_FitPicture/A"+str_NO_Angle+"_E"+str_Energy+"_FitPar.png";
    c1->Print(picture_name);
    /*TGraph2D *g = new TGraph2D(title,"%lg %lg %lg");
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
    for (int i = 0; i < g->GetN(); i++) if (x[i] > distance_cut) g->RemovePoint(i);*/
    
    //TF2* f=new TF2("f2","[4]*Shower_Function.shower_Digi(x,y,1.22,3.2,0,[3],0,0.35)",0,distance_cut,0,45);
    /*TF2* f=new TF2("f2","[6]*Shower_Function.shower_Digi(x,0,[0], [1],  [2],[3],  [4],[5])",0,distance_cut,0,45);
    f->SetParameters( 1.22,  3.15,  0.12,0.887, 0.032,0.035, 100);
    f->SetParLimits(0, 0.8, 1.5);
    f->SetParLimits(1, 1, 10);
    f->SetParLimits(3, 0.1, 1);
    f->SetParLimits(5, 0.01, 0.1);
    //TF2* f=new TF2("f2","[1]*exp(-[0]*x)+0*y",0,distance_cut,0,45);
    //f->SetParameters( 1.25,  100);
    //TF2* f=new TF2("f2","Shower_Function.shower_Digi(x,y,[0],[1],[2],  [3],  [4], [5],[6], [7],[8])",0,distance_cut,0,90);
    //f->SetParameters( 1.22, 1.22, 1.22,   0.3,   3.15,  0.12,0.887,   0.032,0.354);
    //f->SetParLimits(4, 2.5, 3.5);
    //f->SetParLimits(6, 0.5,1);
    //f->SetParLimits(8, 0.1, 0.5);
    
    //f->SetParameters( 1.22, 3.2,   0.1,0.887,    0.025,0.354,  1);
    //f->SetParLimits(0, 1.0, 2.0);
    //f->SetParLimits(1, 0.01,0.5);
    //f->SetParLimits(2, 0.001, 0.05);
    //f->SetParLimits(3, 0, 10);
    
    g->GetXaxis()->SetRangeUser(0,distance_cut);
    g->Draw("p.");
    g->Fit(f,"R");
    par_file << str_Energy << " " << f->GetParameter(0) << " " << f->GetParameter(1) << " " << f->GetParameter(2) << " " << f->GetParameter(3) << " " << f->GetParameter(4) << " " << f->GetParameter(5) << " " << f->GetParameter(6) << " " << f->GetParameter(7) << " " << f->GetParameter(8) << endl;
    //cout << f->GetParameter(0) << ", " << f->GetParameter(1) << ", " << f->GetParameter(2) << ", " << f->GetParameter(3) << ", " << f->GetParameter(4) << ", " << f->GetParameter(5) << ", " << f->GetParameter(6) << ", " << f->GetParameter(7) << ", " << f->GetParameter(8) << endl;*/
    par_file.close();
    return 0;
}
