void print(double min, double max, const double N = 2);
int test3(){
	//const double N = 2.0;
	const double N = 1.3;

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

	print(p11_min , p11_max, N);
	print(p12_min , p12_max, N);
	print(p13_min , p13_max, N);
	cout << "----------------" << endl;

	print(p21_min , p21_max, N);
	print(p22_min , p22_max, N);
	print(p23_min , p23_max, N);
	cout << "----------------" << endl;

	print(p31_min , p31_max, N);
	print(p32_min , p32_max, N);
	print(p33_min , p33_max, N);
	cout << "----------------" << endl;

	print(p41_min , p41_max, N);
	print(p42_min , p42_max, N);
	print(p43_min , p43_max, N);
	return 0;
}
void print(double min, double max, const double N){
	double dev = max - min;
	double mean = (max+min)/2.0;
	cout << mean - N*dev << ", " << mean + N*dev << endl;
}
