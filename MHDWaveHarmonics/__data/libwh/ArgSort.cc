#include "ArgSort.h"

/* Some "C" wrappers to allow python to call the functions easily*/

void SortArrayInt(int n, int *x, int *xs) {
	SortArray(n,x,xs);
}
void SortArrayFlt(int n, float *x, float *xs) {
	SortArray(n,x,xs);
}
void SortArrayDbl(int n, double *x, double *xs) {
	SortArray(n,x,xs);
}
void ArgSortInt(int n, int *x, int *ind) {
	ArgSort(n,x,ind);
}
void ArgSortFlt(int n, float *x, int *ind) {
	ArgSort(n,x,ind);
}
void ArgSortDbl(int n, double *x, int *ind) {
	ArgSort(n,x,ind);
}


/*sorting algorithm based on the quick sort algorithm on wikipedia*/
void SortArray(int n, int *x, int *xs) {
	int *ind = new int[n];
	ArgSort(n,x,ind);
	int i;
	for (i=0;i<n;i++) {
		xs[i] = x[ind[i]];
	}
	delete ind;
}

void SortArray(int n, float *x, float *xs) {
	int *ind = new int[n];
	ArgSort(n,x,ind);
	int i;
	for (i=0;i<n;i++) {
		xs[i] = x[ind[i]];
	}
	delete ind;
}

void SortArray(int n, double *x, double *xs) {
	int *ind = new int[n];
	ArgSort(n,x,ind);
	int i;
	for (i=0;i<n;i++) {
		xs[i] = x[ind[i]];
	}
	delete ind;
}

void ArgSort(int n, int *x, int *ind) {
	/* copy array so we don't affect the original array*/
	int i;
	for (i=0;i<n;i++) {
		ind[i] = i;
	}
	/* now start the sorting process */
	QuickSort(x,ind,0,n-1);
}

void ArgSort(int n, float *x, int *ind) {
	/* copy array so we don't affect the original array*/
	int i;
	for (i=0;i<n;i++) {
		ind[i] = i;
	}
	/* now start the sorting process */
	QuickSort(x,ind,0,n-1);
}

void ArgSort(int n, double *x, int *ind) {
	/* copy array so we don't affect the original array*/
	int i;
	for (i=0;i<n;i++) {
		ind[i] = i;
	}
	/* now start the sorting process */
	QuickSort(x,ind,0,n-1);
}

void Swap(int *a, int *b){
	int tmp;
	tmp = *a;
	*a = *b;
	*b = tmp;
	return;
}

void Swap(float *a, float *b){
	float tmp;
	tmp = *a;
	*a = *b;
	*b = tmp;
	return;
}

void Swap(double *a, double *b){
	double tmp;
	tmp = *a;
	*a = *b;
	*b = tmp;
	return;
}

void QuickSort(int *x, int *ind, int low, int high) {
	if (low < high) {
		/* find partition index */
		int pi = Partition(x,ind,low,high);
		
		QuickSort(x,ind,low,pi-1);
		QuickSort(x,ind,pi+1,high);
	}
}

void QuickSort(float *x, int *ind, int low, int high) {
	if (low < high) {
		/* find partition index */
		int pi = Partition(x,ind,low,high);
		
		QuickSort(x,ind,low,pi-1);
		QuickSort(x,ind,pi+1,high);
	}
}

void QuickSort(double *x, int *ind, int low, int high) {
	if (low < high) {
		/* find partition index */
		int pi = Partition(x,ind,low,high);
		
		QuickSort(x,ind,low,pi-1);
		QuickSort(x,ind,pi+1,high);
	}
}

int Partition(int *x, int *ind, int low, int high) {
	//pivot point is last element of array
	int pivot = x[ind[high]];
	int i = low - 1;
	int j;
	for (j=low;j<high;j++) {
		if (x[ind[j]] <= pivot) {
			i++;
			Swap(&ind[i],&ind[j]);
		}
	}
	Swap(&ind[i+1],&ind[high]);
	return (i+1);
}

int Partition(float *x, int *ind, int low, int high) {
	//pivot point is last element of array
	float pivot = x[ind[high]];
	int i = low - 1;
	int j;
	for (j=low;j<high;j++) {
		if (x[ind[j]] <= pivot) {
			i++;
			Swap(&ind[i],&ind[j]);
		}
	}
	Swap(&ind[i+1],&ind[high]);
	return (i+1);
}

int Partition(double *x, int *ind, int low, int high) {
	//pivot point is last element of array
	double pivot = x[ind[high]];
	int i = low - 1;
	int j;
	for (j=low;j<high;j++) {
		if (x[ind[j]] <= pivot) {
			i++;
			Swap(&ind[i],&ind[j]);
		}
	}
	Swap(&ind[i+1],&ind[high]);
	return (i+1);
}
