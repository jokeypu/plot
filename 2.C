int Fit_DigiEnergy_test(){
    std::string in_name("doc/DigiEnergy_R.txt");
    std::ifstream File_in;
    File_in.open(in_name, std::ios::in);
    
    std::string out_name1("doc/p1.txt");
    std::ofstream p1;
    p1.open(out_name1,std::ios::out);
    std::string out_name2("doc/p2.txt");
    std::ofstream p2;
    p2.open(out_name2,std::ios::out);
    
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
    
    TGraph2D *g = new TGraph2D("doc/DigiEnergy_R.txt","%lg %lg %lg");
    g->SetMarkerStyle(7);
    g->SetMarkerColorAlpha(kAzure+3, 0.5);
    g->GetZaxis()->SetTitle("E_{digi}");
    g->GetXaxis()->SetTitle("distance");
    g->GetYaxis()->SetTitle("angle");
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    g->GetZaxis()->CenterTitle();
    
    Double_t distance_cut = 4;
    Double_t *x = g->GetX();
    for (int i = 0; i < g->GetN(); i++) if (x[i] > distance_cut) g->RemovePoint(i);
    
    TF1* f=new TF1("f2","[1]*x+[0]",0,45);
    TGraph *g1 = new TGraph();
    TGraph *gp1 = new TGraph();
    gp1->SetMarkerStyle(22);
    gp1->SetMarkerColorAlpha(kRed, 0.5);
    gp1->GetXaxis()->SetTitle("distance");
    gp1->GetYaxis()->SetTitle("p1");
    gp1->GetXaxis()->CenterTitle();
    gp1->GetYaxis()->CenterTitle();
    TGraph *gp2 = new TGraph();
    gp2->SetMarkerStyle(33);
    gp2->SetMarkerColorAlpha(kAzure+3, 0.5);
    gp2->GetXaxis()->SetTitle("distance");
    gp2->GetYaxis()->SetTitle("p2");
    gp2->GetXaxis()->CenterTitle();
    gp2->GetYaxis()->CenterTitle();
    
    std::string str;
    int cout(-1),n(0);
    float distance_test(-1);
    while (std::getline(File_in, str)) {
        std::stringstream strStream(str);
        float distance, angle, energy;
        strStream >> distance >> angle >> energy;
        if ( distance >= distance_cut ) continue;
        if (distance != distance_test) {
            if (cout >= 0 ) {
                g1->Fit(f,"R");
                g1->RemovePoint();
                //gp1->SetPoint(cout,distance,(f->GetParameter(0))/exp(-1.5*(distance-1.415)*(distance-1.415)));
                //gp2->SetPoint(cout,distance,(f->GetParameter(1))/exp(-distance));
                gp1->SetPoint(cout,distance_test,(f->GetParameter(0)));
                gp2->SetPoint(cout,distance_test,(f->GetParameter(1)));
                p1 << distance_test << " " << f->GetParameter(0) << endl;
                p2 << distance_test << " " << f->GetParameter(1) << endl;
            }
            distance_test = distance;
            cout++;
            n = 0;
        }
        g1->SetPoint(n,angle,energy);
        n++;
    }
    File_in.close();
    p1.close();
    p2.close();
    
    g->GetXaxis()->SetRangeUser(0,distance_cut);
    //g->Draw("p.");
    c1->cd();
    gp1->Draw("ap.");
    c2->cd();
    gp2->Draw("ap.");
    
    //TF1* fp1=new TF1("fp1","[0]*exp(-[1]*(x-1.5))",1.5,3.5);
    //TF1* fp2=new TF1("fp2","([0]+[1]*x+[2]*x*x+[3]*x*x*x)*exp(-[3]*x)",0,3.5);
    //TF1* fp2=new TF1("fp2","[0]*exp(-[1]*x)+[2]*exp(-[3]*x*x)+[4]*exp(-[5]*x*x*x)+[6]*exp(-[7]*x*x*x*x)",0.2,3.5);
    
    //gp1->Fit(fp1,"R");
    //gp2->Fit(fp2,"R");
    
    
    
    
    return 0;
}
