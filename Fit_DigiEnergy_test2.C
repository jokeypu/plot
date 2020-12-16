int Fit_DigiEnergy_test2(){
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
    //TCanvas* c2=new TCanvas("PANDA2","test2",800,600);
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
    
    TGraph *g = new TGraph("doc/DigiEnergy_R.txt","1.3 %lg %lg");
    g->SetMarkerStyle(33);
    g->SetMarkerColorAlpha(kRed-3, 0.5);
    g->GetXaxis()->SetTitle("E_{digi}");
    g->GetYaxis()->SetTitle("angle");
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    
    TF1* ff=new TF1("ff","[0]+[1]*pow(sin(0.03536776*x),2)+[2]*pow(sin(0.03536776*x),1)",0,90);
    //TF1* ff=new TF1("ff","[0]*pow(x-45,2)+[1]",0,90);
    ff->SetLineWidth(2);
    ff->SetLineColor(kBlack);
    
    g->Draw("ap.");
    g->Fit(ff,"R");
    
    return 0;
}
