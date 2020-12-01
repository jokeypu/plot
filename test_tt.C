int test_tt(){
TGraph *g = new TGraph();
    g->SetMarkerStyle(7);
    g->SetMarkerColorAlpha(kAzure+3, 0.5);
//    g->GetZaxis()->SetTitle("E_{digi}");
    g->GetXaxis()->SetTitle("distance");
    g->GetYaxis()->SetTitle("angle");
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
  //  g->GetZaxis()->CenterTitle();
    
    g->SetPoint(0,0.05
    ,0.6503786433587659
);
    g->SetPoint(1,0.046
    ,0.627248939406038
);
    g->SetPoint(2,0.038963557474332396
    ,0.5837674119066867
);
    g->SetPoint(3,0.032
    ,0.5340168621074897
);
    g->SetPoint(4,0.023481622900254368
    ,0.4636849672014407
);
    g->SetPoint(5,0.014311482483855885
    ,0.3674018096083358
);
    g->SetPoint(6,0.009198382897534801
    ,0.2985092947751189
);
    g->SetPoint(7,0.002605946436549501
    ,0.16588725303228416
);
    g->SetPoint(8,8.26526487213163E-4
    ,0.10638811265053175
);
    g->SetPoint(9,1.1884684630644338E-4
    ,0.033455666611466706
);
    g->SetPoint(10,0.024
    ,0.4681438302811466
);
    g->SetPoint(11,0.008629079417502522
    ,0.28720859014697403
);
    g->Draw("AP.");
    TF1 *f=new TF1("f","[0]*pow(x,0.46)",0,0.07);
    g->Fit(f,"R");
    f->Draw("same");
return 0;
}
