int nn(){
    TF1 *f = new TF1("f","(1-[0]*exp(-pow(x/[1],2)))*x",0,5);
    f->SetParameters(0.2,2,2);
    f->SetParLimits(0,0,1);
    f->Draw();
    return 0;
}

Double_t FABC(Double_t x,Double_t A, Double_t B, Double_t C, Double_t p1, Double_t p2){
    //return p1*exp(-p2*(1.0-A*exp(-pow(x/B,C)))*x);
    return p1*exp(-p2*(1.0-A*exp(-pow(x/B,C)))*x);
    //return p1*exp(-p2*x)*exp(p2*A*exp(-pow(x/B,C))*x);
    //return p1*exp(-p2*x)*exp(A*exp(-pow(x/B,C))*x);
    //return p1*exp(-p2*pow(x/B,C));
}

int z_func_test(){
    TCanvas* c1=new TCanvas("PANDA1","test1",800,600);
    TCanvas* c2=new TCanvas("PANDA2","test2",800,600);
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
    
    //std::string file_name("doc/WorkData_1Gamma_A3_E6.0_OR_R.txt");
    //std::string file_name("doc/WorkData_1Gamma_A3_E6.0_OR.txt");
    std::string file_name("doc/WorkData_1Gamma_A7_E1.0_OR.txt");
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
    Double_t distance_cut = 3;
    while (std::getline(in_file, str)) {
        std::stringstream strStream(str);
        float distance, angle, energy;
        strStream >> distance >> angle >> energy;
        //if (angle>10 || angle<0) continue;
        if (distance > distance_cut) continue;
        g->SetPoint(N,distance,(energy));
        N++;
    }
    
    TF1* f=new TF1("f1","(FABC(x,[0],[1],[2],[3],[4]))",0,distance_cut);
    f->SetParameters(0.9, 1.5, 3, 3.37, 1.45);
    f->SetParLimits(0, 0.5, 1);
    f->SetParLimits(1, 1, 3);
    
    c1->cd();
    g->Draw("AP.");
    g->Fit(f,"R");
    
    in_file.clear();
    in_file.seekg(0, ios::beg);
    TGraph *g_Error = new TGraph();
    g_Error->SetMarkerStyle(7);
    g_Error->SetMarkerColorAlpha(kAzure+3, 0.5);
    g_Error->GetYaxis()->SetTitle("E_{func}-E_{truth}");
    g_Error->GetXaxis()->SetTitle("distance");
    g_Error->GetXaxis()->CenterTitle();
    g_Error->GetYaxis()->CenterTitle();
    N = 0;
    while (std::getline(in_file, str)) {
        std::stringstream strStream(str);
        float distance, angle, energy;
        strStream >> distance >> angle >> energy;
        if (distance > distance_cut) continue;
        g_Error->SetPoint(N,distance,(FABC(distance,f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3),f->GetParameter(4))-energy));
        N++;
    }
    c2->cd();
    g_Error->Draw("AP.");
    
    //nn();
    
    in_file.close();
    return 0;
}
