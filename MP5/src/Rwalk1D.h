#ifndef __RWALK1D__
#define __RWALK1D__

#include <vector>
#include <map>
#include <algorithm>
#include <TCanvas.h>
#include <iostream>
using namespace std;

class Rwalk1D {
    
    static int seed;

    public:
        Rwalk1D(int n = 1, // N = number of particles
                double x = 0.0, // x = x(0)
                double pl = 0.5, double pr = 0.5, // probabilities Left Right
                double dT = 1, double dX = 1 // time and space steps
                );
        ~Rwalk1D() = default; // default destructor

        // particle simulation
        void Run(int nsteps); // number of steps
        void Draw(TCanvas* c); // draw trajectory

        // getters
        const std::vector<double>& GetTrajectory(int n = 1); // particle number
        double GetTimeStep();
        double GetSpaceStep();
        double GetNumber();

    private:
        double x0; // init coo
        int N; // number of particles
        double pL, pR; // probabilities (left, same, righ)
        double dt, dx; // steps (time, space)
        std::map <int, std::vector<double> > mT; // trajectories
};

#endif