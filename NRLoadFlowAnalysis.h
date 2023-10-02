#pragma once

#include<iostream>
#include <stdio.h>
#include <string>
#include <list>
#include <armadillo>
//#include "rtdb/CimDba.h"
//#include <Common\IrDateTime.h>
using namespace std;
using namespace arma;

#define ACCURACY			    0.001

namespace NRLoadFlowAnalysis
{
	int NewtonRaphsonLoadFlow();
	void SetRequiredData(void);
	Mat<double> GetBusInfo();
	int GetNumberofBuses();
	Mat<cx_double> GetFTMAT();
	Mat<cx_double> GetLineInfo();
	Mat<double> GetVoltage();
	Mat<double> SetVoltageAngle();
	Row<double> SetScheduledActivePower();
	Row<double> SetScheduledReactivePower();
	complex<double> SetBaseApparentPower();
	void CreateYbus();
	Mat<double> CreateJacob();
	Mat<double> CreateJ11();
	Mat<double> CreateJ12();
	Mat<double> CreateJ21();
	Mat<double> CreateJ22();
	void CalculatePQI();
	void UpdateDeltaVectores();
}
