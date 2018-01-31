
//    cqcorr_1Fluxes.c  Computes carbon and water component fluxes from correlation analysis
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

/* Calculates fluxes from correlation analysis output
   Proprietary libraries not included and must be purchased from http://numerical.recipes */


#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"


#define N 2

int main(void)
{
    float sigmac,sigmaq,wue,fluxco2,fluxh2o,corrcq,corrcpcr,sigmacp,PD_add,CE,fluxco2_hp,fluxco2_nom;
    sigmac=1.822587;
    sigmaq=0.026931;
    fluxh2o=0.007336;
    fluxco2=0.062225;
    corrcq=-0.340076;
    
    fluxco2_hp=0.038211;
    
    // input here
    wue=-1*7.272259;
    
    fluxco2_nom=-1*0.238311;
    corrcpcr=-1*0.947267;
    sigmacp=0.133835;
    // input here
    
    float qtilda,ctilda1,Tf,GPP,ER;
    PD_add=fluxco2_hp-fluxco2_nom;
    qtilda=SQR(corrcpcr)*(-1+sqrt(1-(1/SQR(corrcpcr))*(1-SQR(sigmaq*wue/sigmacp))));
    ctilda1=SQR(corrcpcr)*(-1+sqrt(1-(1/SQR(corrcpcr))*(1-SQR(sigmac/sigmacp))));
    Tf=1/(1+qtilda);
    CE=(fluxco2*(1-(1/(1+ctilda1))))+PD_add;
    ER=CE-PD_add;
    GPP=fluxco2-CE;
    printf("\n%12s %12s %12s %12s %12s %12s %12s\n","qtilda","ctilda","Tf","GPP","ER","PD_add","CE");
    printf("%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n\n",qtilda,ctilda1,Tf,GPP,ER,PD_add,CE);
    return 0;
}
#undef NRANSI
