#ifndef __MP2_H__
#define __MP2_H__
#include <iostream>
#include <vector>
using namespace std;

struct cell{
    float center_coo[3] = {0, 0, 0}; // cm
    float normal[3]; // unitary vector
    float area; //cm^2
    float power; //W
};

struct source{
    float center_coo[3] = {0, 0, 0};
    float power = 100;
};

void makegrid(cell**, int, int);
vector<vector<cell>> makegrid(int &, int &);
double irradiancia(double phi, double a, double r);
double GetR(source src, cell c);
double GetCos(source src, cell c);
double irradiancia(source &src, cell &cl);
double potencia(source src, cell c);
void FillGrid(vector<vector<cell>> &grid, source src, double comprimento, double largura);
bool sort_cell(cell a, cell b);
void sorting_vcells(vector<cell> vcell);
vector<cell> transform_dvectovec(vector<vector<cell>> vcell);
double total_power(vector<vector<cell>> grid);

#endif