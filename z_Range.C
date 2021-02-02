string print(double min, double max, const double N = 2){
    double dev = max - min;
    double mean = (max+min)/2.0;
    //cout << mean - N*dev << ", " << mean + N*dev << endl;
    
    ostringstream out1,out2;
    out1 << mean - N*dev;
    out2 << mean + N*dev;
    string str_out1 = out1.str(), str_out2 = out2.str();
    return str_out1 + ", " + str_out2;
}

void e1(){
    const double N = 3;
    const double N_m = 1.3;

	double p11_min = -1.1;double p11_max = 0.2;
	double p12_min = 2;double p12_max = 10;
	double p13_min = 2.75;double p13_max = 2.95;


	double p21_min = -1.2;double p21_max = -0.4;
	double p22_min = 6;double p22_max = 12;
	double p23_min = 0.885;double p23_max = 0.91;

	double p31_min = 0.12;double p31_max = 0.3;
	double p32_min = 4;double p32_max = 7;
	double p33_min = 0.72;double p33_max = 0.78;

	double p41_min = -4.2;double p41_max = -3;
	double p42_min = 0.7;double p42_max = 1.4;
	double p43_min = 5.6;double p43_max = 7;
    cout << endl;

	cout << "g11->GetYaxis()->SetRangeUser(" << print(p11_min , p11_max, N) << ");" << endl;
	cout << "g12->GetYaxis()->SetRangeUser(" << print(p12_min , p12_max, N) << ");" << endl;
	cout << "g13->GetYaxis()->SetRangeUser(" << print(p13_min , p13_max, N_m) << ");" << endl;
	cout << "//----------------" << endl;

	cout << "g21->GetYaxis()->SetRangeUser(" << print(p21_min , p21_max, N) << ");" << endl;
	cout << "g22->GetYaxis()->SetRangeUser(" << print(p22_min , p22_max, N) << ");" << endl;
	cout << "g23->GetYaxis()->SetRangeUser(" << print(p23_min , p23_max, N_m) << ");" << endl;
	cout << "//----------------" << endl;

	cout << "g31->GetYaxis()->SetRangeUser(" << print(p31_min , p31_max, N) << ");" << endl;
	cout << "g32->GetYaxis()->SetRangeUser(" << print(p32_min , p32_max, N) << ");" << endl;
	cout << "g33->GetYaxis()->SetRangeUser(" << print(p33_min , p33_max, N_m) << ");" << endl;
	cout << "//----------------" << endl;

	cout << "g41->GetYaxis()->SetRangeUser(" << print(p41_min , p41_max, N) << ");" << endl;
	cout << "g42->GetYaxis()->SetRangeUser(" << print(p42_min , p42_max, N) << ");" << endl;
	cout << "g43->GetYaxis()->SetRangeUser(" << print(p43_min , p43_max, N_m) << ");" << endl;
    
    cout << endl;
}

void e2(const double N = 1.3){
    double p1_min = 0;
    double p1_max = 1;
    print(p1_min , p1_max, N);
}
int z_Range(){
    e1();
    return 0;
}
