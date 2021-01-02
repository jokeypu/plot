Double_t FABC(Double_t x,Double_t A, Double_t p1, Double_t p2, Double_t c1, Double_t c2){
    p2 *= (1-exp(-A*pow(x,3)));
    c2 *= 4*(1-exp(-A*pow(x,3)));
    return p1*exp(-p2*x)+c1*exp(-c2*x);
}

int Fit_DigiEnergy_cp(std::string dir_name, const char title[30], Int_t NO_Angle, Double_t Energy){
    //title[20] = "doc/"+ dir_name +"_R.txt"
    ostringstream out1,out2;
    out1 << NO_Angle;
    out2 << fixed << setprecision(1) << Energy;
    string str_NO_Angle = out1.str(), str_Energy = out2.str();
    std::string out_name = "doc/A"+str_NO_Angle+"_par_cp.txt";
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
        g->SetPoint(N,distance,TMath::Log(energy));
        //g->SetPoint(N,distance,energy);
        N++;
    }
    
    TF1* f=new TF1("f1","TMath::Log(FABC(x,[0],[1],[2],[3],[4]))",0,distance_cut);
    //TF1* f=new TF1("f1","[2]*FABC(x,[0],[1],1,[3],[4],[5])",0,distance_cut);
    f->SetParameters(0.18, 0.8, 1.6, 0.0087, 0.055);
    g->Draw("AP.");
    g->Fit(f,"R");
    
    if (f->GetParameter(2)>f->GetParameter(4))
    par_file << str_Energy << " " << f->GetParameter(0) << " " << f->GetParameter(1) << " " << f->GetParameter(2) << " " << f->GetParameter(3) << " " << f->GetParameter(4) << endl;
    else par_file << str_Energy << " " << f->GetParameter(0) << " " << f->GetParameter(3) << " " << f->GetParameter(4) << " " << f->GetParameter(1) << " " << f->GetParameter(2) << endl;

    //c1->SetLogy();
    TString picture_name= "doc/A"+str_NO_Angle+"_FitPicture_cp/A"+str_NO_Angle+"_E"+str_Energy+"_FitPar_cp.png";
    c1->Print(picture_name);
    
    par_file.close();
    return 0;
}
