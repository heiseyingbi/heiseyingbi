#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <ilcplex/ilocplex.h>
#include <time.h>
#include <string>
#include <list>
#include <algorithm>
#include "randomGenerator.h"
#include "functions.h"
using namespace std;

ILOSTLBEGIN

#define inf 0x7fffffff

extern int nbCities;//number of sites
extern double totoalhardeningBudget;// the available hardening budget
extern vector<double> dj, pj;//the demand and failure probability of each site
extern int maxCities;
extern vector<double> fj;// the fixed cost to setup a regular facility at each site
extern vector<double> bj;// the hardening cost for a site which equals to sj+rj*pj
extern vector<vector<double>> cijP;//the transporation cost for unit demand
extern vector<vector<double>> cijB;

extern vector<int> feasiblexJU;//location decisions 
extern vector<int> feasiblexJR;
extern vector<vector<int>> feasibleSolution;// assignment decisions
extern vector<double> assignmentCost;// the assignment cost

extern double cplex_time;
extern double cplex_Result;//record the cplex result

//record the optimal solution
extern vector<int> optimalxJU;//optimal upper bound solution
extern vector<int> optimalxJR;
extern vector<vector<int>> optimalSolution;

// decision variables
extern vector<int> xJU;
extern vector<int> xJR;
extern vector<vector<int>> zIJ;
extern vector<vector<vector<int>>> yIKJ;

//parameters related to the lowerbound solution
extern vector<double> aJ;
extern vector<double> lJ;
extern vector<vector<double>> bIJ;
extern vector<vector<vector<double>>> rIKJ;
extern vector<vector<int>> minKrIKJ;// record the k=argmin rIKJ for j and i, the k!=j
extern vector<double> vJU;// the value when an unreliable facility is built in j
extern vector<double> vJR;// the value when an reliable facility is built in j
extern vector<double> namedaI;
extern vector<vector<double>> miuIJ;
extern vector<double> namedaIDirection;
extern vector<vector<double>> miuIJDirection;

extern double LR_time;
extern double LR_Result;
extern double LR_OptimalityGAP;//gap between the upper bound and the lower bound
extern double LR_PerformanceGAP;//LR_PerformanceGAP=(double)((LR_Result-cplex_Result)/cplex_Result);


//the lagrangian relaxation algorithm to solve the reliable location problem with hardening budget
//the relaxation problem is solved by using the cplex directly
double PostOptimize(int instanceNum, double B, int bRatio, double & initialsolution, double & Localsearchsolution)// the lr method with local search
{
	///////////////////////////////////// load the data///////////////////////////////
	int i, j, k;
	//string outputFilename="twoLayerRFLP_LR_MyModel";
	char numa[10];
	itoa(instanceNum, numa, 10);
	string outputFilename = "Logfile\\";
	outputFilename += numa;
	char numa2[10];
	itoa(bRatio, numa2, 10);
	outputFilename += "-";
	outputFilename += numa2;
	outputFilename += "-RFLPwithBudget_LR1";
	stringstream ss;
	ss << nbCities;
	outputFilename += "-";
	outputFilename += ss.str();
	ofstream fout(outputFilename);
	/////////////////////////////////////////data structures//////////////////////////////
	//record the initial solution value
	double solutionValueBeforeLS;// the initial greedy solution
	double solutionValueAfterLS;// the solution after the local search
	double recordInitialSolution;// initial solution value found by the greedy and local search algorithm
	double improvementRatioLS;//improvementRationLS=(solutionValueBeforeLS-solutionValueAfterLS)/solutionValueBeforeLS;

	double minElmValue;
	int indexMinEle;
	int indexMinEle_vJR;// record the index of the site j*=argmin VjR, this site should be opened if there is no reliable facility in the system
	int iPriamry, iBackup;
	//----------------------LR parameters
	double lowerBound, upperBound;
	double minUpperBound = inf, maxLowerBound = -inf;
	int iterNum = 0;
	int iterNumLim = 10000;// the maximal iteration number
	if (instanceNum == 4)
		iterNumLim = 2000;
	int lastUpdateUpbound = 0;//
	int lastUpdateLowerbound = 0;
	int noIncreaseStepNum = 0;
	int noIncreaseStepNumLim = 30;


	//double elpso = 0.001;//gap between the upper bound and the lower bound
	double elpso = 0.00001;//gap between the upper bound and the lower bound

	double stepSize = 0;
	//Revised
	double stepSizeLim = 0.00001;//0.005
	double smoothParameter = 2;
	double normOfSubgradient = 0;

	time_t t_start, t_end1, t_end2;
	t_start = time(NULL);

	///////////////////////// generate the initial solution//////////////////////////////////
	minElmValue = inf;
	indexMinEle = -1;
	for (j = 0; j < nbCities; j++)
	{
		if (bj[j] > B)
			continue;
		double centerDis = fj[j];
		for (i = 0; i < nbCities; i++)
			centerDis += cijP[i][j] * dj[i];
		if (centerDis < minElmValue)
		{
			minElmValue = centerDis;
			indexMinEle = j;// select site to build an reliable facility
		}
	}
	if (indexMinEle < 0)//test
	{
		cout << "Erro: indexMinEle<0, no facility can be hardened... " << endl;
		getchar();
		return 0;
	}
	feasiblexJR[indexMinEle] = 1;// setup the reliable facility
	upperBound = fj[indexMinEle];//the fixed cost for the reliable facility
	//initialization
	//assign all the customers to the reliable facility£¬ obatain an initial solution
	for (i = 0; i < nbCities; i++)
	{
		feasibleSolution[i][0] = indexMinEle;
		feasibleSolution[i][1] = -1;
		assignmentCost[i] = cijP[i][indexMinEle] * dj[i];
		upperBound += assignmentCost[i];
	}

	solutionValueBeforeLS = upperBound;

	// local search: add an unreliable faciliy or  a reliable facility k, 0->1 0-2
	local_search(nbCities, upperBound, feasiblexJU, feasiblexJR, feasibleSolution);

	//the solution after local search
	solutionValueAfterLS = upperBound;

	initialsolution = solutionValueBeforeLS;
	Localsearchsolution = solutionValueAfterLS;
	return 0;
}
