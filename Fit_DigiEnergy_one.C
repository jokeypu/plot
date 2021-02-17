const Double_t RM = 2.00;
Double_t FABC(Double_t r,Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4){
    Double_t xi = r - p2*r*exp(-pow(r/p3/RM,p4));
    return p0*exp(-p1*xi/RM);
}

int Fit_DigiEnergy_one(std::string dir_name, const char title[30], Int_t NO_Angle, Double_t Energy){
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
    //gStyle->SetPalette(1);
    
    c1->SetLeftMargin(0.13);
    c1->SetRightMargin(0.13);
    c1->SetTopMargin(-0.13);
    c1->SetBottomMargin(0.13);
    
    c2->SetLeftMargin(0.13);
    c2->SetRightMargin(0.13);
    c2->SetTopMargin(-0.13);
    c2->SetBottomMargin(0.13);
    
    std::string file_name(title);
    std::ifstream in_file;
    in_file.open(file_name,std::ios::in);
    
    TGraph *g = new TGraph();
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kYellow-4);
    g->GetYaxis()->SetTitle("E_{truth}   [GeV]");
    g->GetXaxis()->SetTitle("t   [X_{0}]");
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    
    string str;
    Int_t N = 0;
    Double_t distance_cut = 3.5;
    Int_t binx = 200, biny = 200;
    
    TH2D* h = new TH2D("Hist", "h", binx, 0, distance_cut, biny, 0, 1);
    h->SetMarkerStyle(7);
    h->SetMarkerColorAlpha(kAzure+3, 1);
    h->GetYaxis()->SetTitle("E_{digi}/E_{Shower}   [GeV]");
    h->GetXaxis()->SetTitle("r   [cm]");
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->GetZaxis()->CenterTitle();
    
    while (std::getline(in_file, str)) {
        std::stringstream strStream(str);
        float distance, angle, energy;
        strStream >> distance >> angle >> energy;
        if (distance > distance_cut) continue;
        h->Fill(distance,energy/Energy);
        N++;
    }
    
    for (int i = 1; i < 201; i++){
        for (int j = 1; j < 201; j++){
            if (N < 30000) {if (h->GetBinContent(i,j)<2) h->SetBinContent(i,j,0);}
            else if (h->GetBinContent(i,j)<2) h->SetBinContent(i,j,0);
        }
    }
    
    int cunt = 0;
    double d_s = 0.2;
    double  step = (distance_cut-d_s)/15;
    std::map<int, TH1D*> h_map;
    for (double  i = d_s; i < distance_cut; i+=step){
        in_file.clear();
        in_file.seekg(0, ios::beg);
        double dc_min = i, dc_max = dc_min+step;
        Double_t nx = (dc_min + dc_max)/2.0;
        h_map[cunt] = new TH1D(Form("t%f_to_t%f",i,i+step), Form("t%f_to_t%f",i,i+step),200,0,1);
        while (std::getline(in_file, str)) {
            std::stringstream strStream(str);
            float distance, angle, energy;
            strStream >> distance >> angle >> energy;
            if (distance > dc_min && distance < dc_max) h_map[cunt]->Fill(energy/Energy);
        }
        TF1* f_temp = new TF1("f_temp", "[0]*TMath::Gaus(x,[1],[2])",0,1);
        f_temp->SetParameters(1000, 0.5, 1);
        //h_map[cunt]->Fit(f_temp,"R");
        if ( f_temp->GetParameter(1) > 1) continue;
        //g->SetPoint(cunt, nx, f_temp->GetParameter(1));
        g->SetPoint(cunt, nx, h_map[cunt]->GetMean());
        cunt++;
    }
    
    TF1* f=new TF1("f1","(FABC(x,[0],[1],[2],[3],[4]))",0,distance_cut);
    f->SetParameters(0.8, 2.5, 0.9, 1.4, 3);
    f->SetLineWidth(3);
    f->SetLineColor(kRed);
    
    c1->cd();
    h->Draw("PCOLZ");
    g->Draw("Psame");
    g->Fit(f,"R");
    par_file << str_Energy << " " << f->GetParameter(0) << " " << f->GetParameter(1) << " " << f->GetParameter(2) << " " << f->GetParameter(3) << " " << f->GetParameter(4) << endl;
    
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
        h_Error->Fill(distance,(FABC(distance,f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3),f->GetParameter(4))-energy/Energy));
        N++;
    }
    c2->cd();
    h_map[5]->Draw();
    //h_Error->Draw("PCOLZ");
    TString picture_name_error= "doc/A"+str_NO_Angle+"_FitPicture_cp/Error_A"+str_NO_Angle+"_E"+str_Energy+"_FitPar_cp.png";
    c2->Print(picture_name_error);
    
    in_file.close();
    par_file.close();
    return 0;
}
