int NewGaus(){
	TF1 *f1 = new TF1("f1","exp(-1.51*x)+0.005*exp(-0.075*x)",0,8);
	TF1 *f2 = new TF1("f2","1-exp(-0.24*x*x*x)",0,8);
	f1->SetLineColor(kRed);
	f2->SetLineColor(kBlue);
	f1->Draw();
	f2->Draw("same");
return 0;
}
