#ifndef __LIGHTMAP_H__
#define __LIGHTMAP_H__
#include <iostream>
#include <MP2s.h>
#include <vector>
using namespace std;

class lightmap{
    public:
        lightmap(int ncellx, int ncelly, double comprimento, double largura, source src); //number of cells along x and y
        lightmap(const vector<float>& vx, const vector<float>& vy);
        pair <int, int> GetCellIndex(float x, float y); // return cell indices
        pair<float, float> GetCellCoo(int index_x, int index_y); // return cell center_coo
        double GetCellPower(int index_x, int index_y); // return cell center coo
        double GetCellPower(float x, float y); //return cell power Watts

        int GetCellNx(); // get number of cells along x
        int GetCellNy(); // get number of cells along y

        const cell& GetMaxCell(); // get cell with max power
        vector<vector<cell>>& GetCells(); //return cells grid

        private:
            vector<vector<cell>> GRID;
            source srce;
            int ncx;
            int ncy;
            double com;
            double lar;
};

#endif