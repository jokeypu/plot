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
            Int_t n = 3;
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
    Double_t shower_Digi(Double_t distance,Double_t angle){
        if ( angle > 90 && angle <= 180 ) angle = 180 - angle;
        else if ( angle > 45 && angle <= 90 ) angle = 90 - angle;
        angle *= TMath::DegToRad();
        //Double_t p[9] = {L_1, L_2, L_3, p0, p1, p2, p3, p4, p5};
        Double_t p[9] = {1.57344, 0.943614, 1.20918, 0.139749, 4.94142, 1.18306, 0.494592, 3.37466, 1.89294};
        SetPar(1,p[0]);
        Double_t T1 = block_Int(distance,angle,p[4])/p[4];
        SetPar(1,p[1]);
        Double_t T2 = p[5]*block_Int(distance,angle,p[6])/p[6];
        SetPar(1,p[2]);
        Double_t T3 = p[7]*block_Int(distance,angle,p[8])/p[8];
        return p[3]*(T1+T2+T3)/3;
    }
}Shower_Function;

Double_t FABC(Double_t x, Double_t ShowerEnergy){
    Double_t A = 0.0466804*TMath::Log(ShowerEnergy)+0.199612;
    Double_t p2 = 0.280941*exp(-0.56037*ShowerEnergy)+1.42325;
    Double_t c1 = 0.0477413*exp(-1.89313*ShowerEnergy)+0.00457459;
    Double_t c2 = 0.0162451*ShowerEnergy+0.0483701;
    p2 *= (1-exp(-A*pow(x,3)));
    c2 *= (1-exp(-A*pow(x,3)));
    return exp(-p2*x)+c1*exp(-c2*x);
}

int Fit_DigiEnergy_test5(std::string dir_name, Int_t NO_Angle, Double_t Energy){
    string title = "doc/"+ dir_name +"_R.txt";
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
        //g->SetPoint(N,distance,TMath::Log(energy));
        g->SetPoint(N,distance,energy);
        N++;
    }
    
    TF1* f1=new TF1("f1","0.8*FABC(x,1.0)",0,distance_cut);
    f1->SetLineColor(kGreen);
    TF1* f2=new TF1("f1","Shower_Function.shower_Digi(x,20)",0,distance_cut);
    f2->SetLineColor(kBlack);
    //TF1* f=new TF1("f1","[2]*FABC(x,[0],[1],1,[3],[4],[5])",0,distance_cut);
    g->Draw("AP.");
    f1->Draw("SAME");
    f2->Draw("SAME");
    //g->Fit(f,"R");
    
    
    par_file.close();
    return 0;
}
