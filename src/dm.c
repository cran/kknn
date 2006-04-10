/*****************************************************************/
/*
 *  Copyright (C)2004 Klaus Schliep
 *               
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

 
#include <R.h>
#include <Rmath.h>
#include <R_ext/Memory.h> 
#include <R_ext/Print.h>  
#include <R_ext/Utils.h>


void dm(double *learn, double *valid, int *n, int *m, int *p, double *dm, int *cl, int *k, double *mink, double *weights){
int i, j, l, t; 
double tmp, *dvec;
int *cvec;
cvec = (int *) R_alloc((*n), sizeof(int));
dvec = (double *) R_alloc((*n), sizeof(double));

for(j=0;j<(*m);j++){
	for(i=0;i<*n;i++){
		tmp=0.0;
		for(l=0;l<*p;l++){
			tmp+=pow(fabs(learn[i+l*n[0]]-valid[j+l*m[0]]),*mink)* weights[l];
			}
		dvec[i]=pow(tmp,(1.0/(*mink)));
		cvec[i]=i;     
		}
	rsort_with_index(dvec, cvec, *n);
	for(t=0;t<*k;t++){
		cl[j+t * *m]=cvec[t];
		dm[j+t * *m]=dvec[t];
		}
	}
}

