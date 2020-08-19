int test1(){

double theta=3.5,phi=4.5;

TGraph* g=new TGraph();
g->SetPoint(0,theta+0.1,phi);
		g->SetPoint(1,theta+2,phi);
                g->SetMarkerStyle(21);
                
                
                g->SetMarkerColorAlpha(kRed, 0.3);
                g->SetMarkerSize(8);
                g->Draw("AP");



//		TGraph* g1=new TGraph();
//g1->SetPoint(0,theta+1,phi);
  //              g1->SetMarkerStyle(21);


    //            g1->SetMarkerColorAlpha(kBlue, 0.3);
      //          g1->SetMarkerSize(4);
        //        g1->Draw("APSAME");

TMarker *m = new TMarker (theta+1.5, phi, 21);
m->Draw("SAME");
//delete m;
//delete *m;
TMarker *m1 = new TMarker (theta+1, phi, 21);
m1->SetY(phi);
m1->Draw("SAME");

return 0;
}
