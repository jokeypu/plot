void mode_comb_TO_2(int NCut = 5);
void mode_comb_TO_3(int NCut = 5);
void mode_comb_TO_4(int NCut = 5);
int sniff_comb(){
    mode_comb_TO_2(5);
    //mode_comb_TO_3(5);
    //mode_comb_TO_4(5);
	return 0;
}

void mode_comb_TO_2(int NCut){
    string str1, str2;
    std::ofstream out_file;
    out_file.open("sniff_comb_2.txt",std::ios::out);
    ifstream file1;
    file1.open("sniff_dir.txt", ios::in);
    
    ifstream file2;
    file2.open("sniff_dir.txt", ios::in);
    
    int maxEvtNo = 0;
    while (getline(file1,str1)) maxEvtNo++;
    file1.clear();
    file1.seekg(0, ios::beg);
    
    int ievt = 0;
    while (getline(file1,str1)) {
        file2.clear();
        file2.seekg(0, ios::beg);
        if (ievt%(maxEvtNo/100)==0) cout << 100 * (int)ievt/maxEvtNo << "%" << endl;
        while (getline(file2,str2)) {
            string new_str = str1+str2;
            if (new_str.length() < NCut) continue;
            out_file << new_str << endl;
        }
        ievt++;
    }
    
    out_file.close();
    file1.close();
    file2.close();
}

void mode_comb_TO_3(int NCut){
    string str1, str2, str3;
    std::ofstream out_file;
    out_file.open("sniff_comb_3.txt",std::ios::out);
    ifstream file1;
    file1.open("sniff_dir.txt", ios::in);
    
    ifstream file2;
    file2.open("sniff_dir.txt", ios::in);
    
    ifstream file3;
    file3.open("sniff_dir.txt", ios::in);
    
    int maxEvtNo = 0;
    while (getline(file1,str1)) maxEvtNo++;
    file1.clear();
    file1.seekg(0, ios::beg);
    
    int ievt = 0;
    while (getline(file1,str1)) {
        file2.clear();
        file2.seekg(0, ios::beg);
        if (ievt%(maxEvtNo/100)==0) cout << 100 * (int)ievt/maxEvtNo << "%" << endl;
        while (getline(file2,str2)) {
            file3.clear();
            file3.seekg(0, ios::beg);
            while (getline(file3,str3)) {
                string new_str = str1+str2+str3;
                if (new_str.length() < NCut) continue;
                out_file << new_str << endl;
            }
        }
        ievt++;
    }
    
    out_file.close();
    file1.close();
    file2.close();
    file3.close();
}

void mode_comb_TO_4(int NCut){
    string str1, str2, str3, str4;
    std::ofstream out_file;
    out_file.open("sniff_comb_4.txt",std::ios::out);
    ifstream file1;
    file1.open("sniff_dir.txt", ios::in);
    
    ifstream file2;
    file2.open("sniff_dir.txt", ios::in);
    
    ifstream file3;
    file3.open("sniff_dir.txt", ios::in);
    
    ifstream file4;
    file4.open("sniff_dir.txt", ios::in);
    
    int maxEvtNo = 0;
    while (getline(file1,str1)) maxEvtNo++;
    file1.clear();
    file1.seekg(0, ios::beg);
    
    int ievt = 0;
    while (getline(file1,str1)) {
        file2.clear();
        file2.seekg(0, ios::beg);
        if (ievt%(maxEvtNo/100)==0) cout << 100 * (int)ievt/maxEvtNo << "%" << endl;
        while (getline(file2,str2)) {
            file3.clear();
            file3.seekg(0, ios::beg);
            while (getline(file3,str3)) {
                file4.clear();
                file4.seekg(0, ios::beg);
                while (getline(file4,str4)) {
                    string new_str = str1+str2+str3+str4;
                    if (new_str.length() < NCut) continue;
                    out_file << new_str << endl;
                }
            }
        }
        ievt++;
    }
    
    out_file.close();
    file1.close();
    file2.close();
    file3.close();
    file4.close();
}

