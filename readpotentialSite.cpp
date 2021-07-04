#include <iostream>
#include <fstream>
#include <stdlib.h> 
#include <math.h>
#include <vector>
#include <time.h>
#include <string.h>
#include <list>
#include <algorithm>
#include "functions.h"
using namespace std;

int readpotentialSite()// this function is to read data from the 250 node data set
{
	const char* filename = "potentialSite.txt";//test
	ifstream file(filename);
	if (!file)
	{
		cerr << "ERROR: could not open file '" << filename << "' for reading" << endl;
		throw(-1);
	}
	int nbCities = 30;
	int a;
	cout << "[";
	for (int i = 0; i < nbCities - 1; i++)
	{
		file >> a;
		cout << a << ",";
	}
	file >> a;
	cout << a << "]";
	getchar();
	file.close();
	return 1;
}