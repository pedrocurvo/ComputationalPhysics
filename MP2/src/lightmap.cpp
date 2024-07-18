#include <lightmap.h>
#include "MP2s.h"
#include <iostream>
#include <stdexcept>
using namespace std;

lightmap::lightmap(int ncellx, int ncelly, double comprimento, double largura, source src){
    GRID = makegrid(ncellx, ncelly);
    FillGrid(GRID, src, comprimento, largura);
    srce = src;
    ncx = ncellx;
    ncy = ncelly;
    com = comprimento;
    lar = largura;
}

pair<float, float> lightmap::GetCellCoo(int index_x, int index_y){
    pair<float, float> res;
    res.first = GRID[index_x][index_y].center_coo[0];
    res.first = GRID[index_x][index_y].center_coo[1];
    return res;
}
pair <int, int> lightmap::GetCellIndex(float x, float y){
    pair<int, int> res;
    res.first = (int)(x / abs(GRID[0][1].center_coo[0] - GRID[0][0].center_coo[0]));
    res.second = (int)(y / abs(GRID[1][0].center_coo[0] - GRID[0][0].center_coo[0]));
    if(x < 0 || y < 0 || res.first > GRID[0].size() || res.second > GRID.size()){
        throw logic_error("Value out of lightmap");
    }

}

double lightmap::GetCellPower(int index_x, int index_y){
    return GRID[index_x][index_y].power;
}

int lightmap::GetCellNx(){
    return ncx;
}

int lightmap::GetCellNy(){
    return ncy;
}

const cell &lightmap::GetMaxCell(){
    vector<cell> res = transform_dvectovec(GRID);
    sorting_vcells(res);
    return res[res.size() - 1];
}

vector<vector<cell>>& lightmap::GetCells(){
    return GRID;
}