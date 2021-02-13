//
//  apdbifurcationmain.cpp
//
//  Main file for generating APD bifurcation diagrams.
//  .txt output files are organized as:
//  1st column:         pacing cycle length (PCL)
//  2nd column:         multiplication factor of Ito (called "itofac")
//  Remaining columns:  beat-to-beat series of action potential durations (APDs)
//  The number of remaining columns depends on the macros BEATS and REMOVEBEATS
//  where #remaining columns = (BEATS - REMOVEBEATS)
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>

#ifndef PRECTYPE
#define PRECTYPE double
#endif

#define TYPECELL TTCellIto
#define TYPECELLSTRING "TT"
#include "Cells/TTCellIto.cpp"


#ifndef numvar1
#define numvar1 1
#endif

#ifndef numvar2
#define numvar2 1
#endif

#ifndef numpcl
#define numpcl 1
#endif

#define NCELLS (numvar1*numvar2*numpcl)

#ifndef BEATS
#define BEATS 20
#endif

#ifndef REMOVEBEATS
#define REMOVEBEATS 10
#endif

#ifndef stimt
#define stimt 100.0L
#endif

#ifndef VARIABLE1
#define VARIABLE1 ibarcafac
#endif

#ifndef VARIABLE2
#define VARIABLE2 nacafac
#endif

#ifndef VARIABLE3
#define VARIABLE3 iksfac
#endif

#ifndef VARIABLE4
#define VARIABLE4 itofac
#endif

#include "APDBifurcation.cpp"

int main(int argc, char *argv[])
{
    const double var1 = atof(argv[1]);
    const double var2 = atof(argv[2]);
    const double var3 = atof(argv[3]);
    const double var4 = atof(argv[4]);
    const double pcl = atof(argv[5]);
    const long double dt = atof(argv[6]);
    
    double _icalfac = 0.8;
    double _nacafac = 5.0;
    double _ikrfac = 0.0;
    double _iksfac = 0.5;
    
    double _itofac = 0.0;
    double _iskfac = 0.0;
    
    APDBifurcation<TYPECELL, NCELLS, BEATS>* h_cells;
    
    h_cells = new APDBifurcation<TYPECELL, NCELLS, BEATS>();
    
    FILE *ap;
    char fileSpecap[100];
    snprintf(fileSpecap, 100, "%sapPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    ap = fopen(fileSpecap, "w");
    
    FILE *allbifs;
    char fileSpec1[100];
    snprintf(fileSpec1, 100, "%sbifsPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    allbifs = fopen(fileSpec1, "w");
    
    FILE *caibifs;
    char fileSpec2[100];
    snprintf(fileSpec2, 100, "%scaiPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    caibifs = fopen(fileSpec2, "w");
    
    FILE *casrbifs;
    char fileSpec3[100];
    snprintf(fileSpec3, 100, "%scasrPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    casrbifs = fopen(fileSpec3, "w");
    
    FILE *kibifs;
    char fileSpec4[100];
    snprintf(fileSpec4, 100, "%skiPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    kibifs = fopen(fileSpec4, "w");
    
    FILE *naibifs;
    char fileSpec5[100];
    snprintf(fileSpec5, 100, "%snaiPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    naibifs = fopen(fileSpec5, "w");
    
    FILE *cpeakbifs;
    char fileSpec6[100];
    snprintf(fileSpec6, 100, "%scpeakPCL_%g_VAR1_%g_VAR2_%g_VAR3_%g_VAR4_%g.txt", TYPECELLSTRING, pcl, var1, var2, var3, var4);
    cpeakbifs = fopen(fileSpec6, "w");\
    
    int index = 0;
    h_cells->Cells.itofac[index] = _itofac;
    h_cells->Cells.iskfac[index] = _iskfac;
      //h_cells->Cells.nai[index] = _nai;
      //h_cells->Cells.ki[index] = _ki;
    h_cells->Cells.nacafac[index] = _nacafac;

    h_cells->Cells.ikrfac[index] = _ikrfac;
    h_cells->Cells.iksfac[index] = _iksfac;
    h_cells->Cells.ibarcafac[index] = _icalfac;
                
    h_cells->Cells.VARIABLE1[index] = var1;
    h_cells->Cells.VARIABLE2[index] = var2;
    h_cells->Cells.VARIABLE3[index] = var3;
    h_cells->Cells.VARIABLE4[index] = var4;

    printf("%d\t%g\n",index,h_cells->Cells.vcfac[index]);

    h_cells->pcls[index] = pcl;
    
    printf("Byte Size of Cells on Device: %lu\n", sizeof(APDBifurcation<TYPECELL, NCELLS, BEATS>) );
    //Now run the program
    
    long double t = -stimt;
    h_cells->dobif(dt, t, ap);

    index = 0;
    fprintf(allbifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
    for (int l = REMOVEBEATS; l < BEATS; l++) {
        fprintf(allbifs, "\t%g", h_cells->apds[BEATS*(index) + l]);
    }
    fprintf(allbifs, "\n");
          
    fprintf(caibifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
    fprintf(casrbifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
    fprintf(kibifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
    fprintf(naibifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
    fprintf(cpeakbifs, "%g\t%g\t%g\t%g\t%g", pcl, var1, var2, var3, var4);
            
    for (int l = 2*REMOVEBEATS; l < 2*BEATS; l++) {
        fprintf(caibifs, "\t%g", h_cells->cais[2*BEATS*(index) + l]);
        fprintf(casrbifs, "\t%g", h_cells->casrs[2*BEATS*(index) + l]);
        fprintf(kibifs, "\t%g", h_cells->kis[2*BEATS*(index) + l]);
        fprintf(naibifs, "\t%g", h_cells->nais[2*BEATS*(index) + l]);
    }
    for (int l = REMOVEBEATS; l < BEATS; l++) {
        fprintf(cpeakbifs, "\t%g", h_cells->cpeaks[BEATS*(index) + l]);
    }
    fprintf(caibifs, "\n");
    fprintf(casrbifs, "\n");
    fprintf(kibifs, "\n");
    fprintf(naibifs, "\n");
    fprintf(cpeakbifs, "\n");

    delete h_cells;
    fclose(allbifs);
    fclose(caibifs);
    fclose(casrbifs);
    fclose(kibifs);
    fclose(naibifs);
    fclose(cpeakbifs);

    return 0;
}


