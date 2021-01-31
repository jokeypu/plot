/*Double_t FABC(Double_t x,Double_t A, Double_t B, Double_t C, Double_t p1, Double_t p2){
    //p2 *= (1.0-(1.0-A)*exp(-pow(x/B,C)));
    p2 *= (1.0-A*exp(-pow(x/B,C)));
    //c2 *= 4*(1-exp(-A*pow(x,3)));
    return p1*exp(-p2*x);
}*/
const Double_t X0 = 0.89;
const Double_t RM = 2.00;
Double_t FABC(Double_t t,Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4){
    Double_t xi = t - p2*t*exp(-pow(X0*t/p3/RM,p4));
    return p0*exp(-p1*xi*X0/RM);
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
    
    TCanvas* c2=new TCanvas("PANDA2","test2",800,600);
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
    gStyle->SetTitleOffset(1.2,"xyz");
    
    std::string file_name(title);
    //std::string file_name = "doc/WorkData_1Gamma_A7_E1.0_OR_R.txt";
    std::ifstream in_file;
    in_file.open(file_name,std::ios::in);
    
    TGraph *g = new TGraph();
    g->SetMarkerStyle(7);
    g->SetMarkerColorAlpha(kAzure+3, 0.5);
    g->GetYaxis()->SetTitle("E_{truth}   [GeV]");
    g->GetXaxis()->SetTitle("t   [X_{0}]");
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    
    string str;
    Int_t N = 0;
    Double_t distance_cut = 3.5/X0;
    
    TH2D* h = new TH2D("Hist","h",200,0,distance_cut,200,0,(Energy));
    h->SetMarkerStyle(7);
    h->SetMarkerColorAlpha(kAzure+3, 0.5);
    h->GetYaxis()->SetTitle("E_{truth}   [GeV]");
    h->GetXaxis()->SetTitle("t   [X_{0}]");
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->GetZaxis()->CenterTitle();
    
    while (std::getline(in_file, str)) {
        std::stringstream strStream(str);
        float distance, angle, energy;
        strStream >> distance >> angle >> energy;
        //if (angle>10 || angle<0) continue;
        if (distance > distance_cut) continue;
        //g->SetPoint(N,distance,energy);
        h->Fill(distance/X0,(energy));
        //g->SetPoint(N,distance,energy);
        N++;
    }
    //cout << N << endl;
    
    for (int i = 1; i < 201; i++){
        for (int j = 1; j < 201; j++){
            if (N < 30000) {if (h->GetBinContent(i,j)<2) h->SetBinContent(i,j,0);}
            else if (h->GetBinContent(i,j)<=2) h->SetBinContent(i,j,0);
            //if (h->GetBinContent(i,j) != 0) h->SetBinContent(i, j, (int)(6*TMath::Log(h->GetBinContent(i,j))));
        }
    }
    
    //TF1* f=new TF1("f1","TMath::Log(FABC(x,[0],[1],[2],[3],[4]))",0,distance_cut);
    TF1* f=new TF1("f1","(FABC(x,[0],[1],[2],[3],[4]))",0,distance_cut);
    f->SetParameters(Energy, 2.5, 0.9, 0.7, 3);
    //f->SetParLimits(0, 0.5, 1);
    //f->SetParLimits(2, 0, 1);
    //f->SetParLimits(2, 1.01, 25);

    //f->SetParameters(1, 0.6, 5, 3.37, 1.45);
    c1->cd();
    //c1->SetLogx();
    //g->Draw("AP.");
    //g->Fit(f,"R");
    h->Draw("PCOLZ");
    h->Fit(f,"R");
    par_file << str_Energy << " " << f->GetParameter(0) << " " << f->GetParameter(1) << " " << f->GetParameter(2) << " " << f->GetParameter(3) << " " << f->GetParameter(4) << endl;
    //if (f->GetParameter(2)>f->GetParameter(4))
    //par_file << str_Energy << " " << f->GetParameter(0) << " " << f->GetParameter(1) << " " << f->GetParameter(2) << " " << f->GetParameter(3) << " " << f->GetParameter(4) << endl;
    //else par_file << str_Energy << " " << f->GetParameter(0) << " " << f->GetParameter(3) << " " << f->GetParameter(4) << " " << f->GetParameter(1) << " " << f->GetParameter(2) << endl;

    //c1->SetLogx();
    //c1->SetLogy();
    TString picture_name= "doc/A"+str_NO_Angle+"_FitPicture_cp/A"+str_NO_Angle+"_E"+str_Energy+"_FitPar_cp.png";
    c1->Print(picture_name);
    
    TGraph *g_Error = new TGraph();
    g_Error->SetMarkerStyle(7);
    g_Error->SetMarkerColorAlpha(kAzure+3, 0.5);
    g_Error->GetYaxis()->SetTitle("E_{func}-E_{truth}   [GeV]");
    g_Error->GetXaxis()->SetTitle("t   [X_{0}]");
    g_Error->GetXaxis()->CenterTitle();
    g_Error->GetYaxis()->CenterTitle();
    
    TH2D* h_Error = new TH2D("Hist_Error","h_Error",200,0,distance_cut,200,-Energy,Energy);
    h_Error->SetMarkerStyle(7);
    h_Error->SetMarkerColorAlpha(kAzure+3, 0.5);
    h_Error->GetYaxis()->SetTitle("E_{func}-E_{truth}   [GeV]");
    h_Error->GetXaxis()->SetTitle("t   [X_{0}]");
    h_Error->GetXaxis()->CenterTitle();
    h_Error->GetYaxis()->CenterTitle();
    h_Error->GetZaxis()->CenterTitle();
    
    in_file.clear();
    in_file.seekg(0, ios::beg);
    N = 0;
    while (std::getline(in_file, str)) {
        std::stringstream strStream(str);
        float distance, angle, energy;
        strStream >> distance >> angle >> energy;
        if (distance > distance_cut) continue;
        //g_Error->SetPoint(N,distance,(FABC(distance,f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3),f->GetParameter(4))-energy));
        h_Error->Fill(distance/X0,(FABC(distance/X0,f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3),f->GetParameter(4))-energy));
        N++;
    }
    c2->cd();
    //g_Error->Draw("AP.");
    h_Error->Draw("PCOLZ");
    TString picture_name_error= "doc/A"+str_NO_Angle+"_FitPicture_cp/Error_A"+str_NO_Angle+"_E"+str_Energy+"_FitPar_cp.png";
    c2->Print(picture_name_error);
    
    in_file.close();
    par_file.close();
    return 0;
}
