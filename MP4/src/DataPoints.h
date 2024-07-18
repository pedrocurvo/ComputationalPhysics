#ifndef __DATAPOINTS_H__
#define __DATAPOINTS_H__


#include <iostream>
#include <vector>
#include "TCanvas.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TSystem.h"
#include <initializer_list>
using namespace std;

class DataPoints { public:
// constructors, destructor
	DataPoints() = default; //default constructor (nothing to be done?)
	DataPoints(int N, double* x, double* y); // build DataPoints from C-arrays of x and y values DataPoints(const std::vector< std::pair<double,double> >&);
	DataPoints(const vector< double> & x, const vector< double> & y);
	DataPoints(const vector< pair<double, double>> &);
	DataPoints(initializer_list< pair<double, double>> L) : P(L) {;}

	
	~DataPoints();
	// getters
	const vector<pair<double,double> >& GetPoints() const;
	void GetGraph(TGraph&);
	// draw points using ROOT object TGraph
	virtual void Draw();
	// friend functions (optional)
	friend ostream& operator<< (ostream&, const DataPoints&);


	protected:
	static TCanvas* c;
	vector<pair<double,double> > P; // points
};

#endif