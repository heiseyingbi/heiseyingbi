
#include "randomGenerator.h"
using namespace std;

double mean(double x[], double phi[], int samplenumber) // discrete case (x, alpha)
{
	int k;
	double expect;

	expect = x[1] * phi[1];
	for (k = 2; k < samplenumber; k++) {
		expect = expect + x[k] * (phi[k] - phi[k - 1]);
	}
	expect = expect + x[samplenumber] * (1 - phi[samplenumber - 1]);
	return expect;

}

double optimistic(double x[], double phi[], int samplenumber, double alpha)
{
	double v;
	int i;

	for (i = 1; i <= samplenumber; i++) {
		if (phi[i] <= 1 - alpha) v = x[i];
		else break;
	}
	return v;
}

double pessimistic(double x[], double phi[], int samplenumber, double alpha)
{
	double v;
	int i;

	for (i = 1; i <= samplenumber; i++) {
		if (phi[i] <= alpha) v = x[i];
		else break;
	}
	return v;
}

double measure(double x[], double phi[], int samplenumber, char relation, double a)
{
	double v;
	int    i;

	if (relation == '<')
		for (i = 1; i <= samplenumber; i++) {
			if (x[i] <= a) v = phi[i];
			else return v;
		}

	if (relation == '>')
		for (i = samplenumber; i >= 1; i--) {
			if (x[i] >= a) v = 1 - phi[i];
			else return v;
		}
}



/////////////////////////////////////////////////////////////////////////////////////////
/// Uncertain Set ///////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

double rectangular(double x, double a, double b)
{
	if (a < x&&x < b) return 1;
	return 0;
}

double triangular(double x, double a, double b, double c)
{
	if (a < x&&x < b) return (x - a) / (b - a);
	if (b <= x && x < c) return (x - c) / (b - c);
	return 0;
}

double trapezoidal(double x, double a, double b, double c, double d)
{
	if (a < x&&x < b) return (x - a) / (b - a);
	if (b <= x && x <= c) return 1;
	if (c < x&&x < d) return (x - d) / (c - d);
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////



void sortminmax(double *p, int min, int max)
{
	int i, j;
	double w;
	for (i = min; i < max; i++)
		for (j = i + 1; j <= max; j++)
			if (*(p + i) > *(p + j)) {
				w = *(p + i);
				*(p + i) = *(p + j);
				*(p + j) = w;
			}
}

void sortmaxmin(double *p, int min, int max)
{
	int i, j;
	double w;
	for (i = min; i < max; i++)
		for (j = i + 1; j <= max; j++)
			if (*(p + i) < *(p + j)) {
				w = *(p + i);
				*(p + i) = *(p + j);
				*(p + j) = w;
			}
}

double findmaxn(double *p, int min, int max, int n)
{
	int i, j;
	double w;
	for (i = min; i <= n; i++)
		for (j = i + 1; j <= max; j++)
			if (*(p + i) < *(p + j)) {
				w = *(p + j);
				*(p + j) = *(p + i);
				*(p + i) = w;
			}
	return *(p + n);
}

double findminn(double *p, int min, int max, int n)
{
	int i, j;
	double w;
	for (i = min; i <= n; i++)
		for (j = i + 1; j <= max; j++)
			if (*(p + i) > *(p + j)) {
				w = *(p + j);
				*(p + j) = *(p + i);
				*(p + i) = w;
			}
	return *(p + n);
}

double myu(double a, double b) // Uniform Distribution
{
	double y;
	if (a > b) {
		printf("\nThe first parameter should be less than the second!");
		exit(1);
	}
	y = (double)rand() / (RAND_MAX);
	return (a + (b - a)*y);
}

double myexp(double beta) // Exponential Distribution
{
	double u;
	do {
		u = myu(0, 1);
	} while (u <= 0 || u >= 1);
	return (-1 * beta*log(u));
}

double myn(double mu, double sigma2) // Normal Distribution
{
	double mu1, mu2, z;
	do {
		mu1 = myu(0, 1);
		mu2 = myu(0, 1);
	} while (mu1 <= 0 || mu1 >= 1);
	z = sqrt(-2 * log(mu1))*sin(2 * 3.14159*mu2);
	return (mu + sqrt(sigma2)*z);
}

int mybe(double p) // Bernoulli Distribution
{
	if (p > 1 || p < 0) {
		printf("\nParameter is out of bound!\n");
		exit(1);
	}
	if (myu(0, 1) <= p) return 1;
	else return 0;
}


int mybn(int n, double p)
{
	int i, y = 0;
	for (i = 1; i <= n; i++)
		y = y + mybe(p);
	return y;
}


double myc(double alpha, double beta)
{
	double u;
	do {
		u = myu(0, 1);
	} while (u <= 0);
	return (alpha - beta / tan(3.14159*u));
}


double myemp(double *a, int n) // Empirical Distribution
{
	int i, m;
	double mu;
	for (i = 1; i < n; i++)
	{
		if (a[i - 1] > a[i])
		{
			printf("\nThe observations should be ordered from small to large!\n");
			exit(1);
		}
	}
	mu = myu(0, 1);
	m = (int)((n - 1)*mu + 1);
	return (a[m - 1] + ((n - 1)*mu - m + 1)*(a[m] - a[m - 1]));
}


double mychi(int k)
{
	int i;
	double z, y = 0;
	if (k < 1)
	{
		printf("\nThe parameter should be a integer larger than zero!\n");
		exit(1);
	}
	for (i = 1; i <= k; i++)
	{
		z = myn(0, 1);
		y = y + z * z;
	}
	return y;
}


double myf(int k1, int k2)
{
	double y1, y2;
	y1 = mychi(k1);
	for (;;) {
		y2 = mychi(k2);
		if (y2 != 0) break;
	}
	return((y1 / k1) / (y2 / k2));
}


double mys(int k)
{
	double y, z;
	z = myn(0, 1);
	y = mychi(k);
	return (z / sqrt(y / k));
}

double myer(int k, double beta)
{
	int i;
	double y;
	y = 0;
	if (k <= 0)
	{
		printf("\nThe first parameter should be a integer larger than zero!\n");
		exit(1);
	}
	for (i = 1; i <= k; i++)
		y = y + myexp(beta);
	return y;
}


double myg(double alpha, double beta)
{
	double x, v;
	x = 0;
	if (alpha <= 0 || beta <= 0)
	{
		printf("\nParameter should all larger than zero!\n");
		exit(1);
	}
mark:
	v = myexp(1);
	x = x + v;
	if (alpha == 1) {
		x = beta * x;
		return x;
	}
	alpha = alpha - 1;
	goto mark;
}


double myb(int alpha, int beta)
{
	double y1, y2;
	y1 = myg(alpha, 1);
	y2 = myg(beta, 1);
	return (y1 / (y1 + y2));
}


double myw(double alpha, double beta)
{
	double v;
	if (alpha <= 0 || beta <= 0)
	{
		printf("\nParameter should all larger than zero!\n");
		exit(1);
	}
	v = myexp(1);
	return (beta*pow(v, (1 / alpha)));
}


int myge(double p)
{
	int x;
	double r;
	if (p >= 1 || p <= 0)
	{
		printf("\nParameter is out of bound!\n");
		exit(1);
	}
	do {
		r = myu(0, 1);
	} while (r <= 0 || r >= 1);
	x = (int)(log(r) / log(1 - p));
	return x;
}


int mynb(int k, double p)
{
	int i, y;
	if (k <= 0)
	{
		printf("\nThe first parameter should be a integer larger than zero!\n");
		exit(1);
	}
	y = 0;
	for (i = 1; i <= k; i++)
		y = y + myge(p);
	return y;
}


double myl(double a, double b)
{
	double mu;
	do {
		mu = myu(0, 1);
	} while (mu <= 0 || mu >= 1);
	return(a - b * log(1 / mu - 1));
}

int myp(double lamd) // Poisson Distribution
{
	int x;
	double b, u;
	if (lamd <= 0)
	{
		printf("\nThe parameter should be larger than zero!\n");
		exit(1);
	}
	x = 0;
	b = 1;
	do
	{
		u = myu(0, 1);
		b = b * u;
		x = x + 1;
	} while (b >= exp(-lamd));
	return(x - 1);
}

double mynormal(double x, double c, double w)
{
	return(exp(-(x - c)*(x - c) / (w*w)));
}

double mylogn(double mu, double sigma2) // Lognormal Distribution
{
	double z;
	z = myn(mu, sigma2);
	return (exp(z));
}

double myt(double a, double b, double m) // Triangular Distribution
{
	double c, u, y;
	if (a >= m || m >= b) {
		printf("\nThe first parameter should be the smallest!");
		printf("\nAnd the second one should be the largest!\n");
		exit(1);
	}
	c = (m - a) / (b - a);
	u = myu(0, 1);
	if (u < c)
		y = sqrt(c*u);
	else
		y = 1 - sqrt((1 - c)*(1 - u));
	return(a + (b - a)*y);
}

double randommean(double x[], int samplenumber)
{
	int i;
	double Expect = 0;

	for (i = 1; i <= samplenumber; i++)
		Expect = Expect + x[i];

	return Expect / samplenumber;
}
