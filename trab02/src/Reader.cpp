#include "Reader.h"
#include <fstream>
#include <iostream>
using namespace std;

vector<vector<double>> Reader::GetData(){

	vector<vector<double>> data;

	fstream file;

    file.open(fname, fstream::in);
    char line[100];

    file.getline(line, 100);
    file.getline(line, 100);

    double val;
    file >> val;
    int j = 0;
    vector <double> vet;

    while(!file.eof())
    {
    	vet.emplace_back(val);
    	for (int i = 0; i < cols-1; i++)
        {
    		file >> val;
    		vet.emplace_back(val);
    	}
    	data.emplace_back(vet);
    	vet.clear();

    	//read first element of line -> here because it sets the eof flag when end of file is reached
    	file >> val;
    	++j;

    }

    file.close();

    return data;
}