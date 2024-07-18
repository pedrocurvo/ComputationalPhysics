#include <algorithm>
#include <cmath>
#include "MP2s.h"
using namespace std;

void makegrid(cell** grid, int x, int y){
    for(int i=0; i < x; i++){
        for(int j=0; j < y; j++){
            cell c;
            grid[i][j] = c;
        }
    }
}

vector<vector<cell>> makegrid(int &rows, int &cols){
    vector<vector<cell>> grid;
    for(int i=0; i < rows; i++){
        vector<cell> row;
        for(int j=0; j < cols; j++){
            cell c;
            row.push_back(c);
        }
        grid.push_back(row);
    }
    return grid;
}

double GetR(source src, cell c){
    return sqrt(pow(src.center_coo[0] - c.center_coo[0], 2) + pow(src.center_coo[1] - c.center_coo[1], 2) + pow(src.center_coo[2] - c.center_coo[2], 2));
}
double GetCos(source src, cell c){
    double norm_vector = sqrt(pow(c.normal[0], 2) + pow(c.normal[1], 2) + pow(c.normal[2], 2));
    double norm_r = GetR(src, c);
    double vector_r[3] = {c.center_coo[0] - src.center_coo[0], c.center_coo[1] - src.center_coo[1], c.center_coo[2] - src.center_coo[2]};
    double inner  = c.normal[0] * vector_r[0] + c.normal[1]* vector_r[1] + c.normal[2]* vector_r[2];
    return abs(inner) / norm_r / norm_vector;
}

double irradiancia(source &src, cell &cl){
    return src.power * GetCos(src, cl) / (4 * M_PI * pow(GetR(src, cl), 2));
}

double potencia(source src, cell c){
   return irradiancia(src, c) * c.area;
}


void FillGrid(vector<vector<cell>> &grid, source src, double comprimento, double largura){
    int ncelly = grid.size();
    int ncellx = grid[0].size();
    double c_cell = comprimento / ncellx;
    double l_cell = largura / ncelly;
    double area1 = comprimento * largura / (ncellx * ncelly);
    for(int i = 0; i < ncelly; i++){
        for(int j = 0; j < ncellx; j++){
            grid[i][j].center_coo[0] = 0.5 * c_cell + c_cell * j;
            grid[i][j].center_coo[1] = 0.5 * l_cell + l_cell * i;
            grid[i][j].center_coo[2] = 0;
            grid[i][j].area = area1;
            grid[i][j].normal[0] = 0;
            grid[i][j].normal[1] = 0;
            grid[i][j].normal[2] = 1;
            grid[i][j].power = potencia(src, grid[i][j]);
        }
    }
}

double total_power(vector<vector<cell>> grid){
    double value = 0;
    for(int i = 0; i < grid.size(); i++){
        for(int j = 0; j < grid[0].size(); j++){
            value += grid[i][j].power;
        }
    }
    return value / (grid.size() * grid[0].size());
}

bool sort_cell(cell a, cell b){
    return a.power < b.power;
}

void sorting_vcells(vector<cell> vcell){
    sort(vcell.begin(), vcell.end(), sort_cell);
}

vector<cell> transform_dvectovec(vector<vector<cell>> vcell){
    vector<cell> res;
    for(int i = 0; i < vcell.size(); i++){
        for(int j = 0; j < vcell[0].size(); j++){
            res.push_back(vcell[i][j]);
        }
    }
    return res;
}