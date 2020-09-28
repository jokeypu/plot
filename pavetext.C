
TCanvas *pavetext(){
   TCanvas *c = new TCanvas("c");
   TPaveText *pt = new TPaveText(.05,.1,.95,.8);
 
   TText *tt0 = pt->AddText("Before");
   TText *t1 = pt->AddText("mean:	1.02170e+00");
   TText *t2 = pt->AddText("sigma:	2.30096e-02");
   pt->AddLine(.0,.5,1.,.5);
   TText *tt1 = pt->AddText("After");
   TText *t3 = pt->AddText("mean:	1.02125e+00");
   TText *t4 = pt->AddText("sigma:	1.62561e-02");
 
   tt0->SetTextColor(kBlue);
   t1->SetTextColor(kBlue);
   t2->SetTextColor(kBlue);
   tt1->SetTextColor(kRed);
   t3->SetTextColor(kRed);
   t4->SetTextColor(kRed);

   t1->SetTextSize(0.08);
   t2->SetTextSize(0.08);
   t3->SetTextSize(0.08);
   t4->SetTextSize(0.08);
 
   pt->Draw();
 
   return c;
}
