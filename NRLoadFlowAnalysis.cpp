
# include <iostream>
# include <armadillo>
# include <conio.h>
# include <string>
# include <stdio.h>
# include <list>
# include "NRLoadFlowAnalysis.h"

using namespace std;
using namespace arma;

//Global Variables

Mat<double> busData;
Mat<cx_double> FTMAT;
Mat<double> linePLimit;
Mat<double> busQLimit;
Mat<double> tap;
Row<double> Psch;
Row<double> Qsch;
complex<double> SB;
complex<double> reall = 1;
complex<double> image = 1i;
int n;
Mat<cx_double> lineData;
Mat<cx_double> Ybus;
Mat<double> Y, teta, V, VA;
Row<double> I, P, Q;
Mat<double> J11, J12, J21, J22, Jacob;
Row<double> deltaP, deltaQ, delPQ,delDV;
list<double> delPQList, delDVList;
Row<double> VControlBuses = zeros<rowvec>(10);



int NRLoadFlowAnalysis::NewtonRaphsonLoadFlow()
{
	Row<double> test;
	bool flag = true;
	int count;
	SetRequiredData();
	CreateYbus();
	CalculatePQI();
	while (flag)
	{
		UpdateDeltaVectores();
		Jacob = CreateJacob();
 		delDV =trans((inv(Jacob)*trans(delPQ)));
		// Update Voltage Magnitude and Angle
		count = 0;
		for (int i = 1; i < n; i++)
		{
			VA(i) = VA(i) + delDV(count);
			count++;
		}
		for (int i = 1; i < n; i++)
		{
			if (any(vectorise(VControlBuses) == i) == 0)
			{
				V(i) = V(i) + delDV(count);
				count++;
			}
		}
      
    // Checking if the accuracy is in the desired range.
		flag = false;
	  for (int i = 0; i < delPQ.size(); i++)
		{
			if (abs(delPQ(i) > ACCURACY))
				flag = true; break;
		}
		
	  CalculatePQI();
	  
	}
	cout << "Voltage Mag(pu) = " << V << endl;
	cout << "Voltage ang(rad) = " << VA << endl;
	cout << "Active Power(mw)= " << P*100 << endl;
	cout << "Reactive Power(mw)= " << Q*100 << endl;
	
	return 0;
}

void NRLoadFlowAnalysis::SetRequiredData()
{
	busData = GetBusInfo();
	FTMAT = GetFTMAT();
	lineData = GetLineInfo();
	n = GetNumberofBuses();
	V = GetVoltage();
	VA = SetVoltageAngle();
	Psch = SetScheduledActivePower();
	Qsch = SetScheduledReactivePower();
	SB = SetBaseApparentPower();
}

Mat<double> NRLoadFlowAnalysis::GetBusInfo()
// busData = [Bus_Number  Vmag  Vang  Pg  Qg   Pl   Ql  Qmax  Qmin   BusType]  
{	    busData << 1 << 1.060  << 0 << 0 << 0  << 0  << 0 << endr
	         	<< 2 << 1 << 0 << 40<< 30<< 20 << 10 << endr
		        << 3 << 1 << 0 << 0 << 0 << 45 << 15 << endr
		        << 4 << 1 << 0 << 0 << 0 << 60 << 10 << endr
		        << 5 << 1 << 0 << 0 << 0 << 60 << 10 << endr;
	return busData;
}

Mat<cx_double> NRLoadFlowAnalysis::GetFTMAT()
// FTMAT = [From   To   R   X]
{	FTMAT << 1 << 2 << 0.02 << 0.06 << endr
		  << 1 << 3 << 0.08 << 0.24 << endr
		  << 2 << 3 << .06  << 0.18 << endr
		  << 2 << 4 << 0.06 << 0.18 << endr
		  << 2 << 5 << 0.04 << 0.12 << endr
		  << 3 << 4 << 0.01 << 0.03 << endr
		  << 4 << 5 << 0.08 << 0.24 << endr;
	return FTMAT;
}

Mat<cx_double> NRLoadFlowAnalysis::GetLineInfo()
{// lineData = [lineNumber    from   to    impedance]                                                                       

	  Mat<cx_double> FTMAT;
	  FTMAT = NRLoadFlowAnalysis::GetFTMAT();
	  Mat<cx_double> lineData(FTMAT.n_rows, 4);
	  Mat<cx_double> lineDataOrg(FTMAT.n_rows, 4);
	lineData.col(0) = linspace(1, FTMAT.n_rows, FTMAT.n_rows)* reall;
	lineData.cols(1, 2) = FTMAT.cols(0, 1);
	lineData.col(3) = 1 / (FTMAT.col(2) + FTMAT.col(3)*image);
  //lineData.col(4) = tap settings; 
	lineDataOrg = lineData;
	return lineData;
}

Mat<double> NRLoadFlowAnalysis::GetVoltage()
{
	Mat<double> v;
	Mat<double> busData;
	busData = NRLoadFlowAnalysis::GetBusInfo();
	v = trans(busData.col(1));
	return v;
}

Mat<double> NRLoadFlowAnalysis::SetVoltageAngle()
{
	int n;
	n = NRLoadFlowAnalysis::GetNumberofBuses();
	Mat<double> va = zeros(1, n);
	return va;
}

int NRLoadFlowAnalysis::GetNumberofBuses()
{
	int n;
	Mat<double> busData;
	busData = GetBusInfo();
	n = busData.n_rows;
	return n;
}

Row<double> NRLoadFlowAnalysis::SetScheduledActivePower()
{
	Psch = zeros<rowvec>(n);
	SB = NRLoadFlowAnalysis::SetBaseApparentPower();
	busData = NRLoadFlowAnalysis::GetBusInfo();
	Psch =real( (trans(busData.col(3)) - trans(busData.col(5))) / SB);
	return Psch;
}
Row<double> NRLoadFlowAnalysis::SetScheduledReactivePower()
{
	Qsch = zeros<rowvec>(n);
	SB = NRLoadFlowAnalysis::SetBaseApparentPower();
	busData = NRLoadFlowAnalysis::GetBusInfo();
	Qsch = real((trans(busData.col(4)) - trans(busData.col(6))) / SB);
	return Qsch;
}
complex<double> NRLoadFlowAnalysis::SetBaseApparentPower()
{
	return 100;
}

void NRLoadFlowAnalysis::UpdateDeltaVectores()
{   
	list<double>::iterator it, it1;
	int count = 0;
	deltaP = zeros<rowvec>(n);
	deltaQ = zeros<rowvec>(n);
    deltaP = Psch - P;
	deltaQ = Qsch - Q;
  // Creating Delta(P,Q,delta,V) matrix
	delPQList.clear();
	delDVList.clear();
		for (int i = 1; i < n; i++)
		{
			delPQList.push_back(deltaP(i));
			delDVList.push_back(VA(i));
		}
		for (int i = 1; i < n; i++)
		{
			if (any(vectorise(VControlBuses) == i) == 0)
			{
				delPQList.push_back(deltaQ(i));
				delDVList.push_back(V(i));
			}
			
		}
	delPQ = zeros<rowvec>(delPQList.size());
	delDV = zeros<rowvec>(delDVList.size());
	
	for (it = delPQList.begin(), it1 = delDVList.begin(); it != delPQList.end(), it1 != delDVList.end(); ++it, ++it1)
	{
		delPQ(count) = *it;
		delDV(count) = *it1;
		count++;
	}
}


Mat<double> NRLoadFlowAnalysis::CreateJacob()
{
   	J11 = CreateJ11();
	J12 = CreateJ12();
	J21 = CreateJ21();
	J22 = CreateJ22(); 
	Jacob = join_vert((join_horiz(J11, J12)), join_horiz(J21, J22));
	return Jacob;
}


Mat<double> NRLoadFlowAnalysis::CreateJ22()
{
	
	double temp = 0;
	double temp1 = 0;
	J22 = zeros(n - 1, n - 1);
	
	for (int i = 1; i < n; i++)
	{
		for (int k = 1; k < n; k++)
		{
			if (k == i)
			{
				temp = 0;
				for (int j = 0; j < n; j++)
				{
					if (j != i)
					{
						temp1 = V(j)*Y(i, j)*sin(teta(i, j) + VA(j) - VA(i));
					}
					else
					{
						temp1 = 2 * V(i)*Y(i, j)*sin(teta(i, j) + VA(j) - VA(i));
					}
					temp = temp - temp1;
			}
				J22(i - 1, k - 1) = temp;
		 }
			else { J22(i - 1, k - 1) = -V(i)*Y(i, k)*sin(teta(i, k) + VA(k) - VA(i)); }
		}
	}
		return J22;
}


Mat<double> NRLoadFlowAnalysis::CreateJ21()
{
	double temp;
	J21 = zeros(n - 1, n - 1);
	for (int i = 1; i < n; i++)
	{
		for (int k = 1; k < n; k++)
		{
			if (k == i)
			{
				temp = 0;
				for (int j = 0; j < n; j++)
				{
					if (j != i)
					{
						temp = temp + V(i)*V(j)*Y(i, j)*cos(teta(i, j) + VA(j) - VA(i));
					}
				}
				J21(i - 1, k - 1) = temp;
			}

			else
			{
				J21(i - 1, k - 1) = -V(i)*V(k)*Y(i, k)*cos(teta(i, k) + VA(k) - VA(i));
			}

		}
	}

	return J21;
}


Mat<double> NRLoadFlowAnalysis::CreateJ12()
{
	double temp;
	double temp1;
	J12 = zeros(n - 1, n - 1);
	for (int i = 1; i < n; i++)
	{
		for (int k = 1; k < n; k++)
		{

			if (k == i)
			{
				temp = 0;
				temp1 = 0;
				for (int j = 0; j < n; j++)
				{
					if (i != j)
					{
						temp1 = V(j)*Y(i, j)*cos(teta(i, j) + VA(j) - VA(i));
					}
					else
					{
						temp1 = 2 * V(i)*Y(i, j)*cos(teta(i, j) + VA(j) - VA(i));
					}
					temp = temp + temp1;
				}

				J12(i - 1, k - 1) = temp;
			}
			else
			{
				J12(i - 1, k - 1) = V(i)*Y(i, k)*cos(teta(i, k) + VA(k) - VA(i));
			}
    	}
	}
	return J12;
}

Mat<double> NRLoadFlowAnalysis::CreateJ11()
{
	J11 = zeros(n-1, n-1);
	for (int i = 1; i < n; i++)
	{
		for (int k = 1; k < n; k++)
		{
			if (k == i)
			{
				double temp = 0;
				for (int j = 0; j < n; j++)
				{
					if (i != j)
						temp = temp + V(i)*V(j)*Y(i, j)*sin(teta(i, j) + VA(j) - VA(i));
					    					
				}
				J11(i - 1, k - 1) = temp;
			}
			else
			{
				J11(i - 1, k - 1) = -V(i)*V(k)*Y(i, k)*sin(teta(i, k) + VA(k) - VA(i));
			}

		}

	}
	return J11;
}

void NRLoadFlowAnalysis::CalculatePQI()
{
	I = zeros<rowvec>(n);
	P = zeros<rowvec>(n);
	Q = zeros<rowvec>(n);
	for (int i = 0; i < n; i++)
	{   
		
		for (int j = 0; j < n; j++)
		{
			
			I(i) = I(i) + Y(i, j)*V(j);
			P(i) = P(i) + V(i)*V(j)*Y(i, j)*cos(teta(i, j) + VA(j) - VA(i));
			Q(i) = Q(i) - V(i)*V(j)*Y(i, j)*sin(teta(i, j) + VA(j) - VA(i));
		}	
	}
}

void NRLoadFlowAnalysis::CreateYbus()
{
	SetRequiredData();
	Ybus = zeros(n, n) * 1i;
	for(int i = 0; i < lineData.n_rows; i++)
	{	complex<double> from, to;
		from = lineData(i, 1);
		to = lineData(i, 2);
		Ybus(real(from)-1,real(to)-1) = -(lineData(i, 3));
		
		Ybus((real(to)-1), (real(from)-1)) = Ybus((real(from)-1),( real(to)-1));
	}
	for (int k = 0; k < n; k++)
	{	Ybus(k, k) = 0;
		for(int lineNum = 0; lineNum < lineData.n_rows; lineNum++)
		{
			if ((real(lineData(lineNum, 1)) == (k + 1)) || (real(lineData(lineNum, 2)) == (k + 1)))
				Ybus(k, k) = Ybus(k, k) + lineData(lineNum, 3);
		}
	}
	Y = abs(Ybus);
	teta = arg(Ybus);
	}




