#include <iostream>
#include "lightmap.h"
#include <fstream>
#include <MP2s.h>
// include class TCanvas from ROOT library #include "TRootCanvas.h"
#include "TH2F.h" // histogram 2D
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TApplication.h" 
#include "TGraph.h" 
#include "TSystem.h"
#include "TF1.h"

using namespace std;

int main(){
    //cell grid[ncellx][ncelly];
    //cell* ptr = &grid[0][0];
    auto lambda = [&](cell** grid, int x, int y){
        for(int i=0; i < x; i++){
        for(int j=0; j < y; j++){
            cell c;
            grid[i][j] = c;
        }
    }
    };
    source src;
    src.center_coo[0]= 100;
    src.center_coo[1]= 150;
    src.center_coo[2]= 100;
    src.power = 100;
    // makegrid(&ptr, ncellx, ncelly);
    lightmap L(2000, 3000, 200, 300, src);
    int nx = 200, ny=300;
    double size_x=200, size_y=300;
    auto h2 = new TH2F("h2", "Mapa de Luz; [cm]; [cm]", nx, 0, size_x, ny, 0, size_y);
    // now we have to fill every histogram cell with the calculated power

    for ( auto& v: L.GetCells() ) {
        for ( auto& c: v) {
            h2->Fill(c.center_coo[0], c.center_coo[1], c.power);
        }
    }
    // Draw
    // - we need to instatiate TApplication to have a graphics display
    // - produce a canvas (tela gráfica) where graphics objects will be placed
    // - draw histogram: check the many options you have available;
    //   here we choose a colored gradient representation
    // - save plot to file
    // - Run application
    TApplication app("app", nullptr, nullptr);
    auto c = new TCanvas("canvas", "lightmap canvas", 0, 0, 1300, 900); // size 800x800
    h2->SetStats(false); // apaga caixa lateral
    h2->Draw("COLZ");
    c->Update(); // update display canvas
    c->SaveAs("lightmap.pdf"); // save graphics in pdf format (eps, png, ...)
    // this gives you control of graphics window buttons
    c->SaveAs("lightmap.png");
    TCanvas *rc = (TCanvas *)c->GetCanvasImp();
    //rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    // without the line that follows ( very last line of your program), nothing is displayed
    app.Run();


    std::ofstream f("lightmap.pgm", std::ios::binary | std::ios::out);
    int maxColorValue = 255;
    // grid of cells [200][300]
    int width=200, height=300;
    f << "P2\n" << width << " " << height << "\n" << maxColorValue << "\n";
    // convert power to grayscale
    // - assign maximal power to 255
    // - assign zero power to 0
    int scale = (100/200) * maxColorValue;
    f << scale << " ";
    f.close();

    

    


    return 0;
}