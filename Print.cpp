
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

extern double cplex_time;
extern double cplex_Result;//record the cplex result

extern double LR_time;
extern double LR_Result;
extern double LR_OptimalityGAP;//gap between the upper bound and the lower bound
extern double LR_PerformanceGAP;//LR_PerformanceGAP=(double)((LR_Result-cplex_Result)/cplex_Result);


int printResults(int instanceNum, int bRatio)// print the statistics for the 3 different alogrithm
{
	/*cout<<" cplex_time="<<cplex_time<<endl;
	cout<<" LR_time="<<LR_time<<endl<<endl;
	cout<<" cplex_Result="<<cplex_Result<<endl;
	cout<<" LR_Result="<<LR_Result<<endl<<endl;
	LR_PerformanceGAP=(double)((LR_Result-cplex_Result)/cplex_Result);
	cout<<" LR_PerformanceGAP="<<LR_PerformanceGAP<<endl;
	cout<<" LR_OptimalityGAP="<<LR_OptimalityGAP<<endl<<endl;*/

	//string outputFilename="twoLayerRFLP_Statistic";
	char numa[10];
	itoa(instanceNum, numa, 10);
	string outputFilename = numa;
	char numa2[10];
	itoa(bRatio, numa2, 10);
	outputFilename += "-";
	outputFilename += numa2;
	outputFilename += "-RFLPwithBudget_Statistic";

	//twoLayerRFLP_LR_NoLocalSearch_MyModel
	ofstream  fout(outputFilename, ios::app);

	fout << "=====================================" << endl;
	fout << " cplex_time=" << cplex_time << endl;
	fout << " LR_time=" << LR_time << endl << endl;

	fout << " cplex_Result=" << cplex_Result << endl;
	fout << " LR_Result=" << LR_Result << endl << endl;
	LR_PerformanceGAP = (double)((LR_Result - cplex_Result) / cplex_Result);
	fout << " LR_PerformanceGAP=" << LR_PerformanceGAP << endl;
	fout << " LR_OptimalityGAP=" << LR_OptimalityGAP << endl << endl;

	fout.close();
	return 1;
}

//this function is used to print the results for the 3 algorithms
int printGAP(int instanceNum, int bRatio)// print the statistics for the 3 different alogrithm
{
	double lr1, lr2;
	double LR_1_PerformanceGAP, LR_2_PerformanceGAP;
	string filename = "150cplex.txt";
	ifstream file(filename);
	if (!file)
	{
		cerr << "ERROR: could not open file '" << filename << "' for reading" << endl;
		throw(-1);
	}
	//for(int i=31;i<instanceNum;i++)
	for (int j = 41; j <= instanceNum; j++)
		file >> cplex_Result;
	file.close();

	string filename1 = "150lr.txt";
	ifstream filelr(filename1);
	if (!filelr)
	{
		cerr << "ERROR: could not open file '" << filename1 << "' for reading" << endl;
		throw(-1);
	}
	//for(int i=31;i<instanceNum;i++)
	for (int j = 41; j <= instanceNum; j++)
	{
		filelr >> lr1;
		filelr >> lr2;
	}
	filelr.close();

	cout << " cplex_Result=" << cplex_Result << endl;
	cout << " LR_1_Result=" << lr1 << endl;
	cout << " LR_2_Result=" << lr2 << endl;
	LR_1_PerformanceGAP = (double)((lr1 - cplex_Result) / cplex_Result);
	cout << " LR_1_PerformanceGAP=" << LR_1_PerformanceGAP << endl;
	LR_2_PerformanceGAP = (double)((lr2 - cplex_Result) / cplex_Result);
	cout << " LR_2_PerformanceGAP=" << LR_2_PerformanceGAP << endl;
	cout << endl;

	char numa[10];
	itoa(instanceNum, numa, 10);
	string outputFilename = numa;
	char numa2[10];
	itoa(bRatio, numa2, 10);
	outputFilename += "-";
	outputFilename += numa2;
	outputFilename += "-RFLPwithBudget_performanceGAP";

	//twoLayerRFLP_LR_NoLocalSearch_MyModel
	ofstream  fout(outputFilename, ios::app);
	fout << " cplex_Result=" << cplex_Result << endl;
	fout << " LR_1_Result=" << lr1 << endl;
	fout << " LR_2_Result=" << lr2 << endl;
	LR_1_PerformanceGAP = (double)((lr1 - cplex_Result) / cplex_Result);
	fout << " LR_1_PerformanceGAP=" << LR_1_PerformanceGAP << endl;
	LR_2_PerformanceGAP = (double)((lr2 - cplex_Result) / cplex_Result);
	fout << " LR_2_PerformanceGAP=" << LR_2_PerformanceGAP << endl;
	fout << endl;

	fout.close();
	return 1;
}



int printResultswithCplex(int instanceNum, int bRatio)// print the statistics for the 3 different alogrithm
{
	//instanceNum
//	string filename;
//	if(instanceNum<=20)
//		 filename="49cplex.txt";
//	else if(instanceNum<=40)
//		 filename="88cplex.txt";
//	else if(instanceNum<=60)
//		 filename="150cplex.txt";
//	ifstream file(filename);
//	if (!file) 
//	{
//		 cerr << "ERROR: could not open file '" << filename << "' for reading" << endl;
//				  throw(-1);
//	}
//	//for(int i=31;i<instanceNum;i++)
//	if(instanceNum<=20)
//		for(int j=1;j<=instanceNum;j++)
//			file >> cplex_Result;
//	else if(instanceNum<=40)
//		for(int j=1;j<=instanceNum-20;j++)
//			file >> cplex_Result;
//	else if(instanceNum<=60)
//		for(int j=1;j<=instanceNum-40;j++)
//			file >> cplex_Result;
//	else if(instanceNum<=80)
//		for(int j=1;j<=instanceNum-60;j++)
//			file >> cplex_Result;
//	file.close();
//
//	string filename1;
//	if(instanceNum<=20)
//		filename1="49cplextime.txt";
//	else if(instanceNum<=40)
//		filename1="88cplextime.txt";
//	else if(instanceNum<=60)
//		filename1="150cplextime.txt";
//
////	string filename1="88cplextime.txt";
//	ifstream filetime(filename1);
//	if (!filetime) 
//	{
//		 cerr << "ERROR: could not open file '" << filename1 << "' for reading" << endl;
//				  throw(-1);
//	}
//	//for(int i=31;i<instanceNum;i++)
//		for(int j=1;j<=instanceNum;j++)
//			filetime >> cplex_time;
//	filetime.close();
//
	/*cout<<" cplex_time="<<cplex_time<<endl;
	cout<<" LR_time="<<LR_time<<endl<<endl;
	cout<<" cplex_Result="<<cplex_Result<<endl;
	cout<<" LR_Result="<<LR_Result<<endl<<endl;
	LR_PerformanceGAP=(double)((LR_Result-cplex_Result)/cplex_Result);
	cout<<" LR_PerformanceGAP="<<LR_PerformanceGAP<<endl;
	cout<<" LR_OptimalityGAP="<<LR_OptimalityGAP<<endl<<endl;
*/
//string outputFilename="twoLayerRFLP_Statistic";
	char numa[10];
	itoa(instanceNum, numa, 10);
	string outputFilename = numa;
	char numa2[10];
	itoa(bRatio, numa2, 10);
	outputFilename += "-";
	outputFilename += numa2;
	outputFilename += "-RFLPwithBudget_Statistic";

	//twoLayerRFLP_LR_NoLocalSearch_MyModel
	ofstream  fout(outputFilename, ios::app);

	fout << "=====================================" << endl;
	fout << " cplex_time=" << cplex_time << endl;
	fout << " LR_time=" << LR_time << endl << endl;
	fout << " cplex_Result=" << cplex_Result << endl;
	fout << " LR_Result=" << LR_Result << endl << endl;
	LR_PerformanceGAP = (double)((LR_Result - cplex_Result) / cplex_Result);
	fout << " LR_PerformanceGAP=" << LR_PerformanceGAP << endl;
	fout << " LR_OptimalityGAP=" << LR_OptimalityGAP << endl << endl;
	fout.close();

	return 1;
}