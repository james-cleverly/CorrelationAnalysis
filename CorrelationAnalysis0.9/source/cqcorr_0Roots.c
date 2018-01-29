
//    cqcorr_0Roots.c  Root finder for correlation analysis,
//                     Newton's globally convergent root finder
//                     (Palatella et al. 2014 Boundary-Layer Meteorology
//                       153:327-337. DOI: 10.1007/s10546-014-9947-x)
//    Copyright (C) 2018  Dr James Cleverly
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

/* Driver for routine newt, Numerical Recipes
   Proprietary libraries not included and must be purchased from http://numerical.recipes */

#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

void funcva(int n,float x[],float f[])
{
    float sigmac,sigmaq,wue,fluxco2,fluxh2o,corrcq;
    sigmac=0.922605;
    sigmaq=0.008662;
    fluxh2o=0.002123;
    corrcq=-0.763302;
    
    // input here
    wue=-1*7.272259;
    
    fluxco2=-1*0.238311;
    // input here

    float qtilda,ctilda1;
    qtilda=SQR(x[1])*(-1+sqrt(1-(1/SQR(x[1]))*(1-SQR(sigmaq*wue/x[2]))));
    ctilda1=SQR(x[1])*(-1+sqrt(1-(1/SQR(x[1]))*(1-SQR(sigmac/x[2]))));
    
    float sigmacr1,wqewcr1;
    sigmacr1=ctilda1*(x[2]/x[1]);
    wqewcr1=qtilda/(ctilda1*wue);
    
    int j;
    j=n;
    f[1]=((fluxco2/fluxh2o)*((1+qtilda)/(1+ctilda1)))-wue;
    f[2]=(1/(sigmac*sigmaq))*((SQR(x[2])/wue)+(x[1]*x[2]*sigmacr1*((1/wue)+wqewcr1))+(SQR(sigmacr1)*wqewcr1))-corrcq;
}

void funcvb(int n,float x[],float f[])
{
    float sigmac,sigmaq,wue,fluxco2,fluxh2o,corrcq;
    sigmac=0.922605;
    sigmaq=0.008662;
    fluxh2o=0.002123;
    corrcq=-0.763302;
    
    // input here
    wue=-1*7.272259;
    
    fluxco2=-1*0.238311;
    // input here
    
    float qtilda,ctilda2;
    qtilda=SQR(x[1])*(-1+sqrt(1-(1/SQR(x[1]))*(1-SQR(sigmaq*wue/x[2]))));
    ctilda2=SQR(x[1])*(-1-sqrt(1-(1/SQR(x[1]))*(1-SQR(sigmac/x[2]))));
    
    float sigmacr2,wqewcr2;
    sigmacr2=ctilda2*(x[2]/x[1]);
    wqewcr2=qtilda/(ctilda2*wue);
    
    int j;
    j=n;
    f[1]=((fluxco2/fluxh2o)*((1+qtilda)/(1+ctilda2)))-wue;
    f[2]=(1/(sigmac*sigmaq))*((SQR(x[2])/wue)+(x[1]*x[2]*sigmacr2*((1/wue)+wqewcr2))+(SQR(sigmacr2)*wqewcr2))-corrcq;
}

#define N 2

int main(void)
{
    int i,check;
	float *x,*f,j,k,l,m;
	
    for (j=0;j<=10;j++) {
        for (k=0;k<=10;k++) {
            l=(j/10)-1;
            m=k;
            x=vector(1,N);
            f=vector(1,N);
            x[1]=l;
            x[2]=m;
            printf("\n%7s %6.1f %5.1f\n","+root",l,m);
            newt(x,N,&check,funcva);
            funcva(N,x,f);
            if (check) printf("Convergence problems.\n");
            if (f[1]<0.001) {
                if (f[1]>-0.001) {
                    if (f[2]<0.001) {
                        if (f[2]>-0.001) {
                            if (x[1]<0) {
                                if (x[1]>-1) {
                                    if (x[2]>0) {
                                        printf("%7s %3s %12s\n","Index","x","f");
                                        for (i=1;i<=N;i++) printf("%5d %12.6f %12.6f\n",i,x[i],f[i]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            free_vector(f,1,N);
            free_vector(x,1,N);
        }
    }
    for (j=0;j<=10;j++) {
        for (k=0;k<=10;k++) {
            l=(j/10)-1;
            m=k;
            x=vector(1,N);
            f=vector(1,N);
            x[1]=l;
            x[2]=m;
            printf("\n%7s %6.1f %5.1f\n","-root",l,m);
            newt(x,N,&check,funcvb);
            funcvb(N,x,f);
            if (check) printf("Convergence problems.\n");
            if (f[1]<0.001) {
                if (f[1]>-0.001) {
                    if (f[2]<0.001) {
                        if (f[2]>-0.001) {
                            if (x[1]<0) {
                                if (x[1]>-1) {
                                    if (x[2]>0) {
                                        printf("%7s %3s %12s\n","Index","x","f");
                                        for (i=1;i<=N;i++) printf("%5d %12.6f %12.6f\n",i,x[i],f[i]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            free_vector(f,1,N);
            free_vector(x,1,N);
        }
    }

}
#undef NRANSI
