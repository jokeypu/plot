int test(){
	std::map<Int_t, std::map<Int_t, Double_t> > shower;
	shower[1882202][7577]=3;
	cout<< shower[1882202][7577] << endl;
	int id =1882202;
	int md =7577;
	if (shower.find(id) == shower.end() || shower[id].find(md) == shower[id].end() )
	shower[id][md] = 3;
	else shower[id][md] += 3;
	cout << shower[id][md] << endl;
	return 0;

}
