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
double RFLPBudget_MyModel_LR1(int instanceNum, double B, int bRatio)// the lr method with local search
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
	int iterNumLim = 2000;// the maximal iteration number
	if (instanceNum == 4)
		iterNumLim = 2000;
	int lastUpdateUpbound = 0;//
	int lastUpdateLowerbound = 0;
	int noIncreaseStepNum = 0;
	int noIncreaseStepNumLim = 30;


	//double elpso = 0.001;//gap between the upper bound and the lower bound
	double elpso = 0.0001;//gap between the upper bound and the lower bound

	double stepSize = 0;
	//Revised
	double stepSizeLim = 0.0001;//0.005
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
		//cout << "Erro: indexMinEle<0, no facility can be hardened... " << endl;
		getchar();
		return 0;
	}
	feasiblexJR[indexMinEle] = 1;// setup the reliable facility
	upperBound = fj[indexMinEle];//the fixed cost for the reliable facility
	//initialization
	//assign all the customers to the reliable facility， obatain an initial solution
	for (i = 0; i < nbCities; i++)
	{
		feasibleSolution[i][0] = indexMinEle;
		feasibleSolution[i][1] = -1;
		assignmentCost[i] = cijP[i][indexMinEle] * dj[i];
		upperBound += assignmentCost[i];
	}

	solutionValueBeforeLS = upperBound;

	if (bRatio == 1)
		// local search: add an unreliable faciliy or  a reliable facility k, 0->1 0-2
		local_search(nbCities, upperBound, feasiblexJU, feasiblexJR, feasibleSolution);

	//the solution after local search
	solutionValueAfterLS = upperBound;


	// to avoid the case when the lagrangian relaxtion can not improve the initial solution
	for (j = 0; j < nbCities; j++)
	{
		optimalxJR[j] = feasiblexJR[j];
		optimalxJU[j] = feasiblexJU[j];
	}
	for (i = 0; i < nbCities; i++)
	{
		optimalSolution[i][0] = feasibleSolution[i][0];
		optimalSolution[i][1] = feasibleSolution[i][1];
	}

	//get the time
	t_end1 = time(NULL);

	//the ratio of the solution improvement by local search
	improvementRatioLS = (solutionValueBeforeLS - solutionValueAfterLS) / solutionValueBeforeLS;
	recordInitialSolution = solutionValueAfterLS;//record the initial solution 

	//the value of initial solution, minimum upper bound 
	minUpperBound = upperBound;
	//cout << "InitialSolution=" << minUpperBound << endl;
	fout << "InitialSolution=" << minUpperBound << endl;

	
	//////////////////////////////////////// the subgradient algrithm/////////////////////////////
	//------------given the initial lagrangian multipliers
	// initialize the lagrangian multipliers
	for (i = 0; i < nbCities; i++)
	{
		namedaI[i] = 0;
		//namedaI[i]=dj[i]*cijP[i][0];
		/*for(j=1;j<nbCities;j++)
			if(namedaI[i]<dj[i]*cijP[i][j])
				namedaI[i]=dj[i]*cijP[i][j];*/
				/*double aa=0;
				for(j=1;j<nbCities;j++)
					aa+=cijP[i][j]+cijB[i][j];
				namedaI[i]=aa/nbCities;*/
	}
	/*double a=fjU[0];
	for(j=1;j<nbCities;j++)
		if(a<fjU[j])
			a=fjU[j];*/
	for (i = 0; i < nbCities; i++)
	{
		for (j = 0; j < nbCities; j++)
			//miuIJ[i][j]=5*(fjR[j]+fjU[j])+dj[i];
		//	miuIJ[i][j]=5*(fj[j]+fj[j])/nbCities+dj[i];
			//miuIJ[i][j]=a/nbCities;	
			miuIJ[i][j] = 0;
	}



	IloEnv env;
	IloModel model(env);// model
	IloCplex cplex(env);
	cplex.extract(model);
	cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	cplex.setError(env.getNullStream());//cplex
	double exceedB;
	//----------------------the subgradient iteration procedure
	while (true)
	{
		///////////////////////solve the relaxation problem and compute the lowerbound////////////////////////
		//update the coefficient 
		for (j = 0; j < nbCities; j++)
		{
			aJ[j] = fj[j];
			for (i = 0; i < nbCities; i++)
				aJ[j] -= miuIJ[i][j];
		}
		for (i = 0; i < nbCities; i++)
			for (j = 0; j < nbCities; j++)
			{
				bIJ[i][j] = dj[i] * cijP[i][j] - namedaI[i];
			}
		for (i = 0; i < nbCities; i++)
			for (k = 0; k < nbCities; k++)
				for (j = 0; j < nbCities; j++)
				{
					if (j == k)
						continue;
					rIKJ[i][k][j] = dj[i] * cijP[i][k] * (1 - pj[k]) + dj[i] * cijB[i][j] * pj[k] + miuIJ[i][k] - namedaI[i];
				}


		//-----------------compute the lowerbound solution
		//compute the VJU and VJR for each facility j 
		for (j = 0; j < nbCities; j++)
		{
			xJU[j] = xJR[j] = 0;
			vJU[j] = aJ[j];// compute the vJU
			// compute the vJR
			vJR[j] = fj[j];
			for (i = 0; i < nbCities; i++)
			{
				if (bIJ[i][j] < 0)
					vJR[j] += bIJ[i][j];
				zIJ[i][j] = 0;//initialize the zIJ
				minKrIKJ[i][j] = -1;//initialize the minKrIKJ
			}
			for (i = 0; i < nbCities; i++)
			{
				minElmValue = inf;
				indexMinEle = -1;
				for (k = 0; k < nbCities; k++)
				{
					if (k == j)
						continue;
					if (rIKJ[i][k][j] < minElmValue)
					{
						minElmValue = rIKJ[i][k][j];
						indexMinEle = k;// record the k=argmin rIKJ for the i and j
					}
					yIKJ[i][k][j] = 0;//initialize the yIKJ
				}
				/*if(indexMinEle<0)
				{
					//cout<<"erro: indexMinEle="<<indexMinEle<<endl;
					break;
				}else*/
				if (rIKJ[i][indexMinEle][j] < 0)//? do I need to add if indexMinEle!=-1?
				{
					vJR[j] += rIKJ[i][indexMinEle][j];
					minKrIKJ[i][j] = indexMinEle;// record the indexMinEle=argmin rIKJ for the i and j
					// which is used for computing yIKJ[i][k][j]
				}
			}
		}// for j
		//solve the reduced problem
		int subProblemNum = 0;
		try {
			//vJU[j]=aJ[j];// compute the vJU
			//vJR[j]=fj[j];
			subProblemNum++;
			IloBoolVarArray xjr(env, nbCities);//knapsack decision variables
			IloBoolVarArray xju(env, nbCities);//knapsack decision variables

			IloExpr obj(env), budgetConstraint(env);
			for (IloInt j = 0; j < nbCities; j++)//objective
			{
				obj += vJU[j] * xju[j];
				obj += vJR[j] * xjr[j];
				budgetConstraint += bj[j] * xjr[j];
			}
			IloObjective lagrange_obj(env);
			lagrange_obj = IloAdd(model, IloMinimize(env, obj));
			obj.end();
			//model.add(IloMinimize(env, obj));
			//constraints
			IloRangeArray rangeConstraints(env);
			for (IloInt j = 0; j < nbCities; j++)
				rangeConstraints.add(xju[j] + xjr[j] <= 1);
			rangeConstraints.add(IloSum(xjr) >= 1);
			rangeConstraints.add(budgetConstraint <= B);
			budgetConstraint.end();
			model.add(rangeConstraints);

			cplex.solve();

			IloNum tolerance = cplex.getParam(IloCplex::EpInt);
			for (IloInt j = 0; j < nbCities; j++)
			{
				if (cplex.getValue(xjr[j]) > 1 - tolerance)
				{
					xJR[j] = 1;
					xJU[j] = 0;
					// compute the value for other decision variables
					for (i = 0; i < nbCities; i++)//setup reliable facility
					{
						if (bIJ[i][j] < 0)
							zIJ[i][j] = 1;
					}
					for (i = 0; i < nbCities; i++)
					{
						if (minKrIKJ[i][j] != -1)
							yIKJ[i][minKrIKJ[i][j]][j] = 1;
					}// for i
				}
				else if (cplex.getValue(xju[j]) > 1 - tolerance)
				{
					xJR[j] = 0;
					xJU[j] = 1;
				}
				else
				{
					xJR[j] = 0;
					xJU[j] = 0;
				}
			}
			model.remove(lagrange_obj);
			model.remove(rangeConstraints);
		}
		catch (IloException &e) {
			env.out() << "ERROR: " << e << endl;
			getchar();
		}
		catch (...) {
			env.out() << "Unknown exception" << endl;
			getchar();
		}

		//compute the xJU[j] and xJR[j] for each facility j
		lowerBound = 0;
		for (i = 0; i < nbCities; i++)
			lowerBound += namedaI[i];
		for (j = 0; j < nbCities; j++)
			lowerBound += xJR[j] * vJR[j] + xJU[j] * vJU[j];




		///////////////////////get a upperbound solution////////////////////////
		// compute the upperBound according to the xJU[j] and xJR[j]
		// get the upper bound based on the lower bound solution
		upperBound = 0;
		for (j = 0; j < nbCities; j++)// location cost
		{
			feasiblexJU[j] = xJU[j];// record the feasible solution
			feasiblexJR[j] = xJR[j];
			if (feasiblexJU[j] == 1 && feasiblexJR[j] == 1)//test
			{
				//cout << "erro: feasiblexJU[j]==true && feasiblexJR[j]==true" << endl;
				return 0;
			}
			if (feasiblexJU[j] == 1 || feasiblexJR[j] == 1)
				upperBound += fj[j];//fixed location cost
		}
		
		// notice : I don't add the fixed cost to the upperBound, so all the reliable facility is opened 
		for (i = 0; i < nbCities; i++)//assignment cost 
		{
			iPriamry = -1;
			iBackup = -1;
			minElmValue = inf;

			for (j = 0; j < nbCities; j++)// find the j =argmin cijB[i][j] 
				if (feasiblexJR[j] == 1 && cijB[i][j] < minElmValue)
				{
					iBackup = j;
					minElmValue = cijB[i][j];
				}
			minElmValue = inf;
			for (k = 0; k < nbCities; k++)// find the min i->k->iBackup
				if (feasiblexJU[k] == 1 && cijP[i][k] * (1 - pj[k]) + cijB[i][iBackup] * pj[k] < minElmValue)
				{
					minElmValue = cijP[i][k] * (1 - pj[k]) + cijB[i][iBackup] * pj[k];
					iPriamry = k;
				}

			for (j = 0; j < nbCities; j++) // find the min i->j
				if (feasiblexJR[j] == 1 && cijP[i][j] < minElmValue)
				{
					minElmValue = cijP[i][j];
					iPriamry = j;//  
					iBackup = -1;//no backup facility,since its primary facility is reliable
				}

			feasibleSolution[i][0] = iPriamry;// record the feasible solution
			feasibleSolution[i][1] = iBackup;
			assignmentCost[i] = dj[i] * minElmValue;
			upperBound += assignmentCost[i];
		}// for each customer i

		//Local Search
		// 原来的最好可行解，当前迭代（heuristic构建）
		//if (bRatio == 1)
		//{
			if (minUpperBound > upperBound)
			{
				lastUpdateUpbound = iterNum;
				minUpperBound = upperBound;
				solutionValueAfterLS = minUpperBound;
				local_search(nbCities, solutionValueAfterLS, feasiblexJU, feasiblexJR, feasibleSolution);
				if (solutionValueAfterLS < minUpperBound)
				{
					minUpperBound = solutionValueAfterLS;
					//	//cout<<"the minUpperBound has been improved by the local search---"<<endl;
					//	fout<<"the minUpperBound has been improved by the local search---"<<endl;
				}
				//record the optimal solution
				for (j = 0; j < nbCities; j++)
				{
					optimalxJR[j] = feasiblexJR[j];
					optimalxJU[j] = feasiblexJU[j];
				}
				for (i = 0; i < nbCities; i++)
				{
					optimalSolution[i][0] = feasibleSolution[i][0];
					optimalSolution[i][1] = feasibleSolution[i][1];
				}

			}// if the upper bound is updated
		//}

		

		if (maxLowerBound < lowerBound)//is the lower bound is improved
		{
			lastUpdateLowerbound = iterNum;
			maxLowerBound = lowerBound;
			////cout<<"the maxLowerBound= "<<maxLowerBound<<" at "<<iterNum+1<<"th iteration"<<endl;
			noIncreaseStepNum = 0;
		}
		else
			noIncreaseStepNum++;

		// stop condition 2
		////cout<<"(minUpperBound-maxLowerBound)="<<minUpperBound-maxLowerBound<<endl;
		if ((minUpperBound - maxLowerBound) / minUpperBound < elpso)
		{
			//cout << "(minUpperBound-maxLowerBound)/minUpperBound=" << (minUpperBound - maxLowerBound) / minUpperBound << " < " << elpso << endl;
			fout << "(minUpperBound-maxLowerBound)/minUpperBound=" << (minUpperBound - maxLowerBound) / minUpperBound << " < " << elpso << endl;
			break;
		}

		//update the stepsize
		if (noIncreaseStepNum >= noIncreaseStepNumLim)
		{
			smoothParameter = smoothParameter / 2;
			noIncreaseStepNum = 0;
		}
		///////////////////////update the lagrangian multipliers////////////////////////
		//compute the direction
		for (i = 0; i < nbCities; i++)
		{
			double d = 1;
			for (k = 0; k < nbCities; k++)
				for (j = 0; j < nbCities; j++)
					if (j != k)
						d -= yIKJ[i][k][j];
			for (j = 0; j < nbCities; j++)
				d -= zIJ[i][j];
			namedaIDirection[i] = 0.3*namedaIDirection[i] + d;
		}
		for (i = 0; i < nbCities; i++)
			for (j = 0; j < nbCities; j++)
			{
				double d = 0;
				for (k = 0; k < nbCities; k++)
					if (k != j)
						d += yIKJ[i][j][k];
				d -= xJU[j];
				miuIJDirection[i][j] = 0.3*miuIJDirection[i][j] + d;
			}
		//compute the stepSize
		stepSize = smoothParameter * (minUpperBound - lowerBound);
		normOfSubgradient = 0;
		for (i = 0; i < nbCities; i++)
			normOfSubgradient += namedaIDirection[i] * namedaIDirection[i];
		for (i = 0; i < nbCities; i++)
			for (j = 0; j < nbCities; j++)
				normOfSubgradient += miuIJDirection[i][j] * miuIJDirection[i][j];
		if (normOfSubgradient < 0.0001)
		{
			break;
			//cout << "the normOfSubgradient is small engough,<0.0001" << endl;
			fout << "the normOfSubgradient is small engough,<0.0001" << endl;
		}
		stepSize = (double)(stepSize / normOfSubgradient);
		if (stepSize < stepSizeLim)
		{
			//cout << "stepSize<stepSizeLim" << endl;
			fout << "stepSize<stepSizeLim" << endl;
			break;
		}
		//update the mulitpliers
		for (i = 0; i < nbCities; i++)
		{
			namedaI[i] = namedaI[i] + stepSize * namedaIDirection[i];
			if (namedaI[i] < 0)
				namedaI[i] = 0;
		}
		for (i = 0; i < nbCities; i++)
			for (j = 0; j < nbCities; j++)
			{
				miuIJ[i][j] = miuIJ[i][j] + stepSize * miuIJDirection[i][j];
				if (miuIJ[i][j] < 0)// do I need it?
					miuIJ[i][j] = 0;
			}

		// the standard subgradient algorithm
		/*stepSize=smoothParameter*(minUpperBound-lowerBound);
		normOfSubgradient=0;
		for(i=0;i<nbCities;i++)
		{
			double d=1;
			for(k=0;k<nbCities;k++)
				for(j=0;j<nbCities;j++)
					if(j!=k)
						d-=yIKJ[i][k][j];
			for(j=0;j<nbCities;j++)
				d-=zIJ[i][j];
			normOfSubgradient+=d*d;
		}
		for(i=0;i<nbCities;i++)
			for(j=0;j<nbCities;j++)
			{
				double a=0;
				for(k=0;k<nbCities;k++)
					if(k!=j)
					{
						a+=yIKJ[i][j][k];
					}
				a-=xJU[j];
				normOfSubgradient+=a*a;
			}
		if(normOfSubgradient<0.0001)
		{
			break;
			//cout<<"the normOfSubgradient is small engough,<0.0001"<<endl;
		}
		stepSize=(double)(stepSize/normOfSubgradient);
		if(stepSize<stepSizeLim)
		{
			//cout<<"stepSize<stepSizeLim"<<endl;
			break;
		}
		for(i=0;i<nbCities;i++)
		{
			double d=1;
			for(k=0;k<nbCities;k++)
				for(j=0;j<nbCities;j++)
					if(j!=k)
						d-=yIKJ[i][k][j];
			for(j=0;j<nbCities;j++)
				d-=zIJ[i][j];
			namedaI[i]=namedaI[i]+stepSize*d;
			if(namedaI[i]<0)
				namedaI[i]=0;// if I use the NRL paper's subgradient method, this can not be deleted
		}
		for(i=0;i<nbCities;i++)
			for(j=0;j<nbCities;j++)
			{
				double d=0;
				for(k=0;k<nbCities;k++)
					if(k!=j)
					{
						d+=yIKJ[i][j][k];
					}
				d-=xJU[j];
				miuIJ[i][j]=miuIJ[i][j]+stepSize*d;
				if(miuIJ[i][j]<0)// do I need it?
					miuIJ[i][j]=0;
			}*/
		iterNum++;
		//cout << minUpperBound << "," << maxLowerBound << endl;
		fout << minUpperBound << "," << maxLowerBound << endl;
		if (iterNum >= iterNumLim)
		{
			//cout << "The LR has already reach the iteration limit" << iterNumLim << ", stop" << endl;
			fout << "The LR has already reach the iteration limit" << iterNumLim << ", stop" << endl;
			break;
		}
	}//while
	env.end();
	t_end2 = time(NULL);
	LR_time = difftime(t_end2, t_start);
	/////////////////////////////////////////// print the results/////////////////////////////////
	double reliableFacilityFixedCost = 0;
	double facilityFixedCost = 0;
	double unreliableFacilityFixedCost = 0;
	double expectedTransportationCost = 0;
	double improvementRatioLR = (double)(recordInitialSolution - minUpperBound) / recordInitialSolution;
	LR_Result = minUpperBound;
	LR_OptimalityGAP = (minUpperBound - maxLowerBound) / minUpperBound;

	//cout << endl;
	fout << endl;
	//cout << "the upperBound found by the initial random solution is " << recordInitialSolution << endl;
	fout << "the upperBound found by the initial random solution is " << recordInitialSolution << endl;

	//cout << "the ratio of local search improvement is " << improvementRatioLS << endl << endl;
	fout << "the ratio of local search improvement is " << improvementRatioLS << endl << endl;

	//cout << "the ratio of LR improvement is " << improvementRatioLR << endl;
	fout << "the ratio of LR improvement is " << improvementRatioLR << endl;

	//cout << "the maximal lowerBound found by the Lagrangian Relaxation is " << maxLowerBound << endl;
	fout << "the maximal lowerBound found by the Lagrangian Relaxation is " << maxLowerBound << endl;

	//cout << "the upperBound found by the Lagrangian Relaxation is " << minUpperBound << endl;
	fout << "the upperBound found by the Lagrangian Relaxation is " << minUpperBound << endl;

	//cout << "the LR gap is " << LR_OptimalityGAP << endl;
	fout << "the LR gap is " << LR_OptimalityGAP << endl;

	//cout << "the time for initial solution generation is " << difftime(t_end1, t_start) << endl;
	fout << "the time for initial solution generation is " << difftime(t_end1, t_start) << endl;

	//cout << "the total cup time is " << LR_time << endl;
	fout << "the total cup time is " << LR_time << endl;

	//cout << "Unreliable facilities:";
	fout << "Unreliable facilities:";
	double checkBudget = 0;
	double checkAssignmentCost = 0;

	int numOfhardeningFacilities = 0;
	int numOfregularFacilities = 0;

	for (i = 0; i < nbCities; i++)
		if (optimalxJU[i] == 1)
		{
			numOfregularFacilities++;
			unreliableFacilityFixedCost += fj[i];
			//cout << i << ",";
			fout << i << ",";
		}
	//cout << endl;
	fout << endl;

	//cout << "Reliable facilities:";
	fout << "Reliable facilities:";
	for (i = 0; i < nbCities; i++)
		if (optimalxJR[i] == 1)
		{
			numOfhardeningFacilities++;
			reliableFacilityFixedCost += fj[i];
			//cout << i << ",";
			fout << i << ",";
			checkBudget += bj[i];
			//checkFixedCost+=fj[i];
		}
	//cout << endl;
	fout << endl;

	for (i = 0; i < nbCities; i++)
		if (optimalSolution[i][1] == -1)
		{
			//cout << i << "->" << optimalSolution[i][0];
			fout << i << "->" << optimalSolution[i][0];
			if (optimalxJR[optimalSolution[i][0]] == 1)
			{
				//cout << "(reliable)" << endl;
				fout << "(reliable)" << endl;
			}
			else
			{
				//cout << "(erro)" << endl;//test
				fout << "(erro)" << endl;//test
			}
			checkAssignmentCost += dj[i] * cijP[i][optimalSolution[i][0]];
		}
		else
		{
			//cout << i << "->" << optimalSolution[i][0];
			fout << i << "->" << optimalSolution[i][0];
			if (optimalxJU[optimalSolution[i][0]] == 1)
			{
				//cout << "(unreliable)";
				fout << "(unreliable)";
			}
			else
			{
				//cout << "(erro)";
				fout << "(erro)";
			}
			//cout << "->" << optimalSolution[i][1];
			fout << "->" << optimalSolution[i][1];
			if (optimalxJR[optimalSolution[i][1]] == 1)
			{
				//cout << "(reliable)" << endl;
				fout << "(reliable)" << endl;
			}
			else
			{
				//cout << "(erro)" << endl;
				fout << "(erro)" << endl;
			}
			int theK = optimalSolution[i][0];
			int theJ = optimalSolution[i][1];
			checkAssignmentCost += dj[i] * cijP[i][theK] * (1 - pj[theK]) + dj[i] * cijB[i][theJ] * pj[theK];
		}
	//cout << endl;
	fout << endl;
	expectedTransportationCost = checkAssignmentCost;
	//cout << "Available budget=" << B << ",Used budget is " << checkBudget << endl;
	fout << "Available budget=" << B << ",Used budget is " << checkBudget << endl;
	facilityFixedCost = reliableFacilityFixedCost + unreliableFacilityFixedCost;
	//cout << "reliableFacilityFixedCost=" << reliableFacilityFixedCost << endl;
	fout << "reliableFacilityFixedCost=" << reliableFacilityFixedCost << endl;
	//cout << "unreliableFacilityFixedCost=" << unreliableFacilityFixedCost << endl;
	fout << "unreliableFacilityFixedCost=" << unreliableFacilityFixedCost << endl;
	//cout << "facilityFixedCost=" << facilityFixedCost << endl;
	fout << "facilityFixedCost=" << facilityFixedCost << endl;
	//cout << "expectedTransportationCost=" << expectedTransportationCost << endl;
	fout << "expectedTransportationCost=" << expectedTransportationCost << endl;
	fout << endl;
	//cout << endl;

	//cout << "numOfopenedFacilities=" << numOfhardeningFacilities + numOfregularFacilities << endl;
	fout << "numOfopenedFacilities=" << numOfhardeningFacilities + numOfregularFacilities << endl;

	//cout << "numOfhardeningFacilities=" << numOfhardeningFacilities << endl;
	fout << "numOfhardeningFacilities=" << numOfhardeningFacilities << endl;

	//cout << "numOfregularFacilities=" << numOfregularFacilities << endl;
	fout << "numOfregularFacilities=" << numOfregularFacilities << endl;

	fout << endl;
	//cout << endl;


	//cout << "the difference beteween the objective and  LR_solution is " << abs(LR_Result - checkAssignmentCost - facilityFixedCost) << endl;
	fout << "the difference beteween the objective and  LR_solution is " << abs(LR_Result - checkAssignmentCost - facilityFixedCost) << endl;
	//cout << endl;
	//cout << "lastUpdateUpperbound= " << lastUpdateUpbound << endl;
	//cout << "lastUpdateLowerbound= " << lastUpdateLowerbound << endl;
	//	file.close();
	fout.close();
	return maxLowerBound;
}