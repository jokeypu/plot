Double_t ShowerEnergy;
void SetPar(Double_t E){
    ShowerEnergy = E;
}

Double_t FABC(Double_t x){
    Double_t A = 0.04585*TMath::Log(ShowerEnergy)+0.1983;
    Double_t p2 = 0.2675*exp(-0.4915*ShowerEnergy)+1.426;
    Double_t c1 = 0.04491*exp(-1.465*ShowerEnergy)+0.004277;
    Double_t c2 = 0.01389*ShowerEnergy+0.02699;
    p2 *= (1-exp(-A*pow(x,3)));
    c2 *= (1-exp(-A*pow(x,3)));
    return exp(-p2*x)+c1*exp(-c2*x);
}

int Fit_DigiEnergy_test4(std::string dir_name, Int_t NO_Angle, Double_t Energy){
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
    
    SetPar(Energy);
    //TF1* f=new TF1("f1","TMath::Log(FABC(x,[0],[1],[2],[3],[4]))",0,distance_cut);
    TF1* f=new TF1("f1","[0]*FABC(x)",0,distance_cut);
    //f->SetParameters(0.119652, 0.78, 1.8, 0.0198265, 0.088);
    g->Draw("AP.");
    g->Fit(f,"R");
    f->Draw("SAME");

    //TString picture_name= "doc/A"+str_NO_Angle+"_FitPicture/A"+str_NO_Angle+"_E"+str_Energy+"_FitPar.png";
    //c1->Print(picture_name);
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
