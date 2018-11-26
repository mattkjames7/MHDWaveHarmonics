#ifndef __ArgSort_h__
#define __ArgSort_h__
#include <stdio.h>
#endif
using namespace std;
extern "C" {
	void SortArrayInt(int n, int *x, int *xs);
	void SortArrayFlt(int n, float *x, float *xs);
	void SortArrayDbl(int n, double *x, double *xs);
	void ArgSortInt(int n, int *x, int *ind);
	void ArgSortFlt(int n, float *x, int *ind);
	void ArgSortDbl(int n, double *x, int *ind);
}
void SortArray(int n, int *x, int *xs);
void SortArray(int n, float *x, float *xs);
void SortArray(int n, double *x, double *xs);
void ArgSort(int n, int *x, int *ind);
void ArgSort(int n, float *x, int *ind);
void ArgSort(int n, double *x, int *ind);
void Swap(int *a, int *b);
void Swap(float *a, float *b);
void Swap(double *a, double *b);
void QuickSort(int *x, int *ind, int low, int high);
void QuickSort(float *x, int *ind, int low, int high);
void QuickSort(double *x, int *ind, int low, int high);
int Partition(int *x, int *ind, int low, int high);
int Partition(float *x, int *ind, int low, int high);
int Partition(double *x, int *ind, int low, int high);
