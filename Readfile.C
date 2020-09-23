int Readfile(){
	ifstream readfile;
    	readfile.open("doc/Fit_par_1.txt", ios::in);
	string line;
	int cunt(0), N(0);
	while (readfile >> line) {
		if ((cunt % 3) == 2){
		cout  << line << ", ";
		N++;
		}
		cunt++;
	}
	cout << "N:" << N << endl;
	readfile.close();
	return 0;
}
