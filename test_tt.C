Double_t FABC(Double_t x, Double_t ShowerEnergy){
    Double_t A = 0.0466804*TMath::Log(ShowerEnergy)+0.199612;
    Double_t p2 = 0.280941*exp(-0.56037*ShowerEnergy)+1.42325;
    Double_t c1 = 0.0477413*exp(-1.89313*ShowerEnergy)+0.00457459;
    Double_t c2 = 0.0162451*ShowerEnergy+0.0483701;
    p2 *= (1-exp(-A*pow(x,3)));
    c2 *= (1-exp(-A*pow(x,3)));
    return exp(-p2*x)+c1*exp(-c2*x);
}

int test_tt(){
    TCanvas* c1=new TCanvas("PANDA1","c1",800,600);
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
    gStyle->SetNdivisions(510,"xyz");
    gStyle->SetTitleFont(42,"xyz");
    gStyle->SetTitleColor(1,"xyz");
    gStyle->SetTitleSize(0.05,"xyz");
    gStyle->SetTitleOffset(1.0,"xyz");
    
    TF1* f1=new TF1("f1","(FABC(x,0.5))",0,8);
    TF1* f2=new TF1("f2","(FABC(x,1))",0,8);
    TF1* f3=new TF1("f3","(FABC(x,1.5))",0,8);
    TF1* f4=new TF1("f4","(FABC(x,2))",0,8);
    TF1* f5=new TF1("f5","(FABC(x,2.5))",0,8);
    TF1* f6=new TF1("f6","(FABC(x,3))",0,8);
    TF1* f7=new TF1("f7","(FABC(x,3.5))",0,8);
    TF1* f8=new TF1("f8","(FABC(x,4))",0,8);
    TF1* f9=new TF1("f9","(FABC(x,4.5))",0,8);
    TF1* f10=new TF1("f10","(FABC(x,5))",0,8);
    TF1* f11=new TF1("f11","(FABC(x,5.5))",0,8);
    TF1* f12=new TF1("f12","(FABC(x,6))",0,8);
    
    f12->SetLineColor(12);
    f11->SetLineColor(11);
    f10->SetLineColor(10);
    f9->SetLineColor(9);
    f8->SetLineColor(8);
    f7->SetLineColor(7);
    f6->SetLineColor(6);
    f5->SetLineColor(5);
    f4->SetLineColor(4);
    f3->SetLineColor(3);
    f2->SetLineColor(2);
    f1->SetLineColor(1);
    
    c1->SetLogy();
    f12->Draw();
    f12->GetXaxis()->CenterTitle();
    f12->GetYaxis()->CenterTitle();
    f12->GetXaxis()->SetTitle("distance");
    f12->GetYaxis()->SetTitle("Energy");
    f2->Draw("same");
    f3->Draw("same");
    f4->Draw("same");
    f5->Draw("same");
    f6->Draw("same");
    f7->Draw("same");
    f8->Draw("same");
    f9->Draw("same");
    f10->Draw("same");
    f11->Draw("same");
    f1->Draw("same");
    
    
return 0;
}
