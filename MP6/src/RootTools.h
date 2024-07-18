#ifndef __ROOTTOOLS_H__
#define __ROOTTOOLS_H__

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMultiGraph.h>
#include <TH1.h>
#include <TH2.h>

using namespace std;

class RootTools {
    public:
        RootTools();
        ~RootTools() = default;
        static void TDrawVector(double, double, double);
        static void TDrawVector(double, double);
        static void CanvasUpdater(TCanvas*, bool = false, string = "");
        static void CanvasUpdater(TCanvas, bool = false, string = "");

};

#endif