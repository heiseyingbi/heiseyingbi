
#include <ilcplex/ilocplex.h>
#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <time.h>
#include <string>
#include <list>
#include <algorithm>
#include "randomGenerator.h"
#include "functions.h"
using namespace std;

ILOSTLBEGIN

extern int nbCities;//number of sites
extern double totoalhardeningBudget;// the available hardening budget
extern vector<double> dj, pj;//the demand and failure probability of each site
extern int maxCities;
extern vector<double> fj;// the fixed cost to setup a regular facility at each site
extern vector<double> bj;// the hardening cost for a site which equals to sj+rj*pj
extern vector<vector<double>> cijP;//the transporation cost for unit demand
extern vector<vector<double>> cijB;

extern double cplex_time;
extern double cplex_Result;//record the cplex result

//this funciton is used to do the managerial analysis
int RFLPBudget_MyModel_Cplex_MA(int instanceNum, double B, int bRatio)
{
	IloEnv env;
	try {
		///////////////// read data file ////////////////////////////////
		//loadData(instanceNum);
		 //build the model		 
		int  npCities = 30;
		int potentialSite[30] = { 0,1,4,8,10,14,15,21,22,24,30,32,37,39,43,47,50,53,54,55,57,65,68,75,77,80,87,90,92,95 };

		IloBoolVarArray  XjU(env, npCities);
		IloBoolVarArray  XjR(env, npCities);
		IloArray<IloNumVarArray>  Zij(env);
		for (IloInt i = 0; i < nbCities; i++)
			Zij.add(IloNumVarArray(env, npCities, 0.0, 1.0, ILOBOOL));
		IloArray<IloArray<IloNumVarArray>> Yikj(env);
		for (IloInt i = 0; i < nbCities; i++)
		{
			IloArray<IloNumVarArray> item(env);
			for (IloInt k = 0; k < npCities; k++)
				item.add(IloNumVarArray(env, npCities, 0.0, 1.0, ILOBOOL));
			Yikj.add(item);
		}
		IloModel model(env);
		//the objective function
		IloExpr oFixedCost(env), oTravelCost(env);
		for (IloInt j = 0; j < npCities; j++)
			oFixedCost += fj[potentialSite[j]] * XjR[j];
		for (IloInt j = 0; j < npCities; j++)
			oFixedCost += fj[potentialSite[j]] * XjU[j];
		for (IloInt i = 0; i < nbCities; i++)
			for (IloInt j = 0; j < npCities; j++)
				oTravelCost += cijP[i][potentialSite[j]] * dj[i] * Zij[i][j];
		for (IloInt i = 0; i < nbCities; i++)
			for (IloInt k = 0; k < npCities; k++)
				for (IloInt j = 0; j < npCities; j++)
				{
					if (j == k)
						continue;
					oTravelCost += cijP[i][potentialSite[k]] * dj[i] * (1 - pj[potentialSite[k]])*Yikj[i][k][j];
					oTravelCost += cijB[i][potentialSite[j]] * dj[i] * pj[potentialSite[k]] * Yikj[i][k][j];
				}
		model.add(IloMinimize(env, oFixedCost + oTravelCost));
		oFixedCost.end();
		oTravelCost.end();
		// the constraints
		 // we can at most build an reliable or unreliable facility at one place
		for (IloInt j = 0; j < npCities; j++)
		{
			model.add(XjU[j] + XjR[j] <= 1);
		}
		// each customer should be assinged to a reliable facility or two layer facilities
		for (IloInt i = 0; i < nbCities; i++)
		{
			IloExpr expr_primary(env);
			for (IloInt k = 0; k < npCities; k++)
				for (IloInt j = 0; j < npCities; j++)
				{
					if (j == k)
						continue;
					expr_primary += Yikj[i][k][j];
				}
			for (IloInt j = 0; j < npCities; j++)
				expr_primary += Zij[i][j];
			model.add(expr_primary == 1);
			expr_primary.end();
		}
		// the facilities in the one layer assignment should be reliable
		for (IloInt i = 0; i < nbCities; i++)
			for (IloInt j = 0; j < npCities; j++)
			{
				model.add(Zij[i][j] <= XjR[j]);
			}
		// if there are two layer assignment,then the primary one should be unreliable
		for (IloInt i = 0; i < nbCities; i++)
			for (IloInt k = 0; k < npCities; k++)
			{
				IloExpr expr_primaryUnreliable(env);
				for (IloInt j = 0; j < npCities; j++)
				{
					if (j == k)
						continue;
					expr_primaryUnreliable += Yikj[i][k][j];
				}
				model.add(expr_primaryUnreliable <= XjU[k]);
				expr_primaryUnreliable.end();
			}
		// if there are two layer assignment,then the backup facility should be reliable
		for (IloInt i = 0; i < nbCities; i++)
			for (IloInt j = 0; j < npCities; j++)
			{
				IloExpr expr_backup(env);
				for (IloInt k = 0; k < npCities; k++)
				{
					if (j == k)
						continue;
					expr_backup += Yikj[i][k][j];
				}
				model.add(expr_backup <= XjR[j]);
				expr_backup.end();
			}
		// the should be at least one reliable facility in the system
		IloExpr expr_reliableNum(env);
		for (IloInt j = 0; j < npCities; j++)
			expr_reliableNum += XjR[j];
		model.add(expr_reliableNum >= 1);
		expr_reliableNum.end();
		//model.add(IloSum(XjR) >= 1 );
		// the cost for reliable facility can not exceed the hardening budget
		IloExpr expr_knapsack(env);
		for (IloInt j = 0; j < npCities; j++)
			expr_knapsack += XjR[j] * bj[potentialSite[j]];
		model.add(expr_knapsack <= B);
		expr_knapsack.end();

		// Optimize
		IloCplex cplex(model);
		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());

		cout << "MyModel CPLEX solving..." << endl;
		time_t t_start, t_end;
		t_start = time(NULL);
		cplex.solve();
		t_end = time(NULL);
		cplex_time = difftime(t_end, t_start);
		//////output the result///////
		char numa[10];
		itoa(instanceNum, numa, 10);
		string outputFilename = numa;
		char numa2[10];
		itoa(bRatio, numa2, 10);
		outputFilename += "-";
		outputFilename += numa2;
		outputFilename += "-RFLPwithBudget_Cplex";
		stringstream ss;
		ss << nbCities;
		outputFilename += "-";
		outputFilename += ss.str();
		ofstream  fout(outputFilename);

		cout << "the cup time is " << difftime(t_end, t_start) << endl;
		fout << "the cup time is " << difftime(t_end, t_start) << endl;

		if (cplex.getStatus() == IloAlgorithm::Infeasible)
			env.out() << "No Solution" << endl;
		// Print results
		double reliableFacilityFixedCost = 0;
		double unreliableFacilityFixedCost = 0;
		double facilityFixedCost = 0;
		double expectedTransportationCost = 0;

		cout << "Cost=" << cplex.getObjValue() << endl;
		fout << "Cost=" << cplex.getObjValue() << endl;
		cplex_Result = cplex.getObjValue();

		double checkTotalCost = 0;// to test if the total cost corresponding to the location and assignment
		double cplexTotalCost = cplex.getObjValue();
		//decisions and the cplex objective are the same, if not output erro
		int numOfhardeningFacilities = 0;
		int numOfregularFacilities = 0;

		cout << "Unreliable facilities:" << endl;
		fout << "Unreliable facilities:" << endl;
		IloNum tolerance = cplex.getParam(IloCplex::EpInt);
		//potentialSite
		for (IloInt j = 0; j < npCities; j++)
			if (cplex.getValue(XjU[j]) > 1 - tolerance)
			{
				numOfregularFacilities++;
				unreliableFacilityFixedCost += fj[potentialSite[j]];
				checkTotalCost += fj[potentialSite[j]];
				cout << potentialSite[j] << " ";
				fout << potentialSite[j] << " ";
			}
		cout << endl;
		fout << endl;
		double usedBudget = 0;
		cout << "Reliable facilities:" << endl;
		fout << "Reliable facilities:" << endl;
		for (IloInt j = 0; j < npCities; j++)
			if (cplex.getValue(XjR[j]) > 1 - tolerance)
			{
				numOfhardeningFacilities++;
				reliableFacilityFixedCost += fj[potentialSite[j]];
				checkTotalCost += fj[potentialSite[j]];
				cout << potentialSite[j] << " ";
				fout << potentialSite[j] << " ";
				usedBudget += bj[potentialSite[j]];
			}
		cout << endl;
		fout << endl;
		cout << "the customer assignments are " << endl;
		fout << "the customer assignments are " << endl;
		for (IloInt i = 0; i < nbCities; i++)
			for (IloInt j = 0; j < npCities; j++)
			{
				if (cplex.getValue(Zij[i][j]) > 1 - tolerance)
				{
					cout << i << "->" << potentialSite[j] << "(reliable)" << endl;
					fout << i << "->" << potentialSite[j] << "(reliable)" << endl;
					checkTotalCost += cijP[i][potentialSite[j]] * dj[i];
				}
				else
				{
					for (IloInt k = 0; k < npCities; k++)
					{
						if (k == j)
							continue;
						if (cplex.getValue(Yikj[i][j][k]) > 1 - tolerance)
						{
							cout << i << "->" << potentialSite[j] << "(unreliable)" << "->" << potentialSite[k] << "(reliable)" << endl;
							fout << i << "->" << potentialSite[j] << "(unreliable)" << "->" << potentialSite[k] << "(reliable)" << endl;
							checkTotalCost += cijP[i][potentialSite[j]] * dj[i] * (1 - pj[potentialSite[j]]) + cijB[i][potentialSite[k]] * dj[i] * pj[potentialSite[j]];
						}
					}
				}
			}
		cout << endl;
		fout << endl;
		expectedTransportationCost = cplex.getObjValue() - reliableFacilityFixedCost - unreliableFacilityFixedCost;
		facilityFixedCost = reliableFacilityFixedCost + unreliableFacilityFixedCost;
		cout << "reliableFacilityFixedCost=" << reliableFacilityFixedCost << endl;
		fout << "reliableFacilityFixedCost=" << reliableFacilityFixedCost << endl;
		cout << "unreliableFacilityFixedCost=" << unreliableFacilityFixedCost << endl;
		fout << "unreliableFacilityFixedCost=" << unreliableFacilityFixedCost << endl;
		cout << "facilityFixedCost=" << facilityFixedCost << endl;
		fout << "facilityFixedCost=" << facilityFixedCost << endl;
		cout << "expectedTransportationCost=" << expectedTransportationCost << endl;
		fout << "expectedTransportationCost=" << expectedTransportationCost << endl;

		fout << endl;
		cout << endl;

		cout << "numOfopenedFacilities=" << numOfhardeningFacilities + numOfregularFacilities << endl;
		fout << "numOfopenedFacilities=" << numOfhardeningFacilities + numOfregularFacilities << endl;

		cout << "numOfhardeningFacilities=" << numOfhardeningFacilities << endl;
		fout << "numOfhardeningFacilities=" << numOfhardeningFacilities << endl;

		cout << "numOfregularFacilities=" << numOfregularFacilities << endl;
		fout << "numOfregularFacilities=" << numOfregularFacilities << endl;

		fout << endl;
		cout << endl;
		cout << "total budget=" << B << ", usedBudget=" << usedBudget << endl;
		fout << "total budget=" << B << ", usedBudget=" << usedBudget << endl;
		cout << "TEST: the difference between the checkTotalCost and cpelx objective is" << abs(cplexTotalCost - checkTotalCost) << endl;
		fout << "TEST: the difference between the checkTotalCost and cpelx objective is" << abs(cplexTotalCost - checkTotalCost) << endl;
		fout.close();
	}
	catch (IloException &e) {
		env.out() << "ERROR: " << e << endl;
		getchar();
	}
	catch (...) {
		env.out() << "Unknown exception" << endl;
		getchar();
	}
	env.end();
	return 0;
}