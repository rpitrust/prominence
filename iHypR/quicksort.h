/*
 * File:   quicksort.h
 * Author: xiaohui
 *
 * Created on March 20, 2011, 2:50 PM
 */

#ifndef QUICKSORT_H
#define	QUICKSORT_H

#include "basic.h"

// randomization switch for quick sort
// 1 = implement randomized quick sort
// 0 = normal quick sort
#define RANDMIZED	1

// partition of the quick sort alg.
int partition(Aux A[], int p, int r);

// randomized quick sort partition routine
int randomized_partition(Aux A[], int p, int r);

// main entry of quick sort
void randomized_quicksort(Aux A[], int p, int r);

#endif	/* QUICKSORT_H */

