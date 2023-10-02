#include <iostream>
#include <armadillo>
#include <conio.h>
#include "NRLoadFlowAnalysis.h"
#include <complex>
#include <cmath>
using namespace arma;
using namespace std;
int main() {
	
	NRLoadFlowAnalysis::NewtonRaphsonLoadFlow();
	
	_getch();
	return 0;

	}