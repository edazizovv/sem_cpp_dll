// SEM-ext.h - Contains declarations of SEM-functions
// see https://docs.microsoft.com/en-us/cpp/build/walkthrough-creating-and-using-a-dynamic-link-library-cpp?view=vs-2019
#pragma once

#ifdef SEM_EXPORTS
#define SEM_API __declspec(dllexport)
#else
#define SEM_API __declspec(dllimport)
#endif

//extern "C" SEM_API int test();
extern "C" SEM_API double* one_optima(int data_len, double* data, double low_prec);
extern "C" SEM_API double* sem(int niter, int k, double* data, int n);
//extern "C" SEM_API double sem(int niter, int k, double* data, int n);