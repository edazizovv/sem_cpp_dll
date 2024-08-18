// SEM.cpp : Defines the exported functions for the DLL.
#include "pch.h"
#include <iostream>
#include <boost/math/special_functions/polygamma.hpp>
#include <math.h>
#include <nlopt.hpp>
#include <iomanip>
#include <vector>

#include<gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "SEM.h"

double one_get_cpp(double teth, void* others)
{
	double k;
	double ki_sum;
	double ki_sumln;
	double poly;
	double funa;

	// repacking 

	double* dated = static_cast<double*>(others);
	int n = static_cast<int>(*dated);

	// computations

	ki_sum = 0;
	ki_sumln = 0;
	for (int i = 1; i <= n; i++)
	{
		ki_sum += *(dated + i);
		ki_sumln += log(*(dated + i));
	}

	k = ki_sum / (n * teth);
	poly = boost::math::polygamma(0, k);
	funa = -(n * poly) - (n * log(teth)) + ki_sumln;

	return abs(funa);
}

double targetta(const std::vector<double>& x, std::vector<double>& grad, void* my_func_data)
{
	double resa;
	resa = one_get_cpp(x[0], my_func_data);
	return resa;
}

double* one_optima(int data_len, double* data, double low_prec)
{
	low_prec = 1e-10;

	double data_sum;
	double teth_gained;
	double ka_gained;

	int over = data_len + 1;

	double* dat = new double[over];

	*dat = (double)data_len;
	for (int i = 0; i < data_len; i++)
	{
		*(dat + i + 1) = *(data + i);
	}

	std::vector<double> x(1);
	x[0] = 1.0;

	nlopt::opt opt(nlopt::LN_NELDERMEAD, 1);
	std::vector<double> lb(1);
	lb[0] = low_prec;
	opt.set_lower_bounds(lb);
	opt.set_min_objective(targetta, dat);

	double minf;

	try
	{
		nlopt::result result = opt.optimize(x, minf);
	}
	catch (std::exception & e)
	{
		x[0] = -999999;
	}


	teth_gained = x[0];

	data_sum = 0.0;
	for (int i = 0; i < data_len; i++)
	{
		data_sum += *(data + i);
	}

	ka_gained = data_sum / ((double)data_len * teth_gained);
	double rei = 0;

	double* res = new double[3]{ teth_gained, ka_gained, rei };

	return res;
}

double* givemedata(double* donor, int len)
{
	double* result = new double[len];
	for (int i = 0; i < len; i++)
	{
		result[i] = donor[i];
	}
	return result;
}

double* sem(int niter, int k, double* data, int n)
{

	// Initialisation

	double** G = new double* [n];
	double** y = new double* [n];
	double* P = new double[k];
	int* V = new int[k];
	double* tethas_gained = new double[k];
	double* ks_gained = new double[k];
	double** clusters = new double* [k];
	int* clusters_n = new int[k];
	double* res;
	int recounter;
	double numerator;
	double denominator;

	//		here we initialise a random number generator
	const gsl_rng* r = gsl_rng_alloc(gsl_rng_taus);

	//		sampling size
	size_t K = k;

	//		number of items
	unsigned int N = 1;

	//		where to save samples
	//unsigned int* neine;

	for (int j = 0; j < n; j++)
	{
		G[j] = new double[k];
		y[j] = new double[k];

		for (int i = 0; i < k; i++)
		{
			G[j][i] = 1.0 / (double)k;
		}

	}


	//iterations
	for (int t = 0; t < niter; t++)
	{
		//std::cout << std::endl;
		//std::cout << "Iteration: " << t << std::endl;
		//std::cout << std::endl;

		// S-step

		for (int j = 0; j < n; j++)
		{

			const double* p = givemedata(G[j], k);

			/*
			if (j == 0)
			{
				std::cout << std::endl;
				for (int dd = 0; dd < k; dd++)
				{
					std::cout << p[dd] << "\t";
				}
				std::cout << std::endl;
				char a;
				std::cin >> a;
			}
			*/
			// temp array to save sample
			unsigned int* neine = new unsigned int[k];

			// generate
			gsl_ran_multinomial(r, K, N, p, neine);

			// copy data
			//std::cout << std::endl;
			for (int i = 0; i < k; i++)
			{
				y[j][i] = neine[i];
				//std::cout << neine[i] << "\t";
			}
			//std::cout << std::endl;

			delete[] neine;
			delete[] p;
		}

		// clusters
		for (int i = 0; i < k; i++)
		{
			V[i] = 0;
		}
		for (int j = 0; j < n; j++)
		{
			for (int i = 0; i < k; i++)
			{
				if (y[j][i] == 1)
				{
					V[i]++;
				}
			}
		}

		// clusters data
		for (int i = 0; i < k; i++)
		{
			clusters[i] = new double[V[i]];
		}


		for (int i = 0; i < k; i++)
		{
			recounter = 0;
			for (int j = 0; j < n; j++)
			{
				if (y[j][i] == 1)
				{
					clusters[i][recounter] = data[j];
					recounter++;
				}
			}
		}

		// M-step

		// gaining
		for (int i = 0; i < k; i++)
		{
			res = one_optima(V[i], clusters[i], 0);
			tethas_gained[i] = res[0];
			ks_gained[i] = res[1];
			// outer niggers from alien gay space
			//std::cout << tethas_gained[i] << "\t";
		}
		//std::cout << std::endl;
		//char pudge;
		//std::cin >> pudge;

		// del clusters
		for (int i = 0; i < k; i++)
		{
			delete[] clusters[i];
		}

		// WUT?
		for (int i = 0; i < k; i++)
		{
			P[i] = (double)V[i] / (double)n;
			//std::cout << P[i] << "\t";
		}
		//std::cout << std::endl;

		// E-step
		for (int j = 0; j < n; j++)
		{
			for (int i = 0; i < k; i++)
			{
				numerator = P[i] * gsl_ran_gamma_pdf(data[j], tethas_gained[i], ks_gained[i]);
				denominator = 0;
				for (int l = 0; l < k; l++)
				{
					denominator = denominator + (P[l] * gsl_ran_gamma_pdf(data[j], tethas_gained[i], ks_gained[i]));
				}
				G[j][i] = numerator / denominator;
			}
		}
	}

	// return

	double* rs = new double[((k * 2) + 2)];

	for (int i = 0; i < k; i++)
	{
		rs[i * 2] = tethas_gained[i];
		rs[i * 2 + 1] = ks_gained[i];
	}

	for (int r = 0; r < k; r++)
	{
		rs[(r + k * 2)] = P[r];
	}

	return rs;
	//return data;

}
