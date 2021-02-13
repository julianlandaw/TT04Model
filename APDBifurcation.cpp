//
//  APDBifurcation.cpp
//
//  Structure for performing analysis of action potential durations (APDs)
//  given the pacing cycle length (PCL).
//
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//
#ifndef PRECTYPE
#define PRECTYPE double
#endif

#include "APDBifurcation.h"

#ifndef APDBifurcation_cpp
#define APDBifurcation_cpp

#ifndef threshold
#define threshold -75.0
#endif

#ifndef stimulus
#define stimulus -80.0
#endif

#ifndef stimduration
#define stimduration 0.5
#endif

#ifndef DV_MAX
#define DV_MAX 0.1
#endif

#ifndef ADAPTIVE
#define ADAPTIVE 10
#endif

template <template<int> class typecell, int ncells, int beats>
APDBifurcation<typecell, ncells, beats>::APDBifurcation() {
    int i;
    for (i = 0; i < ncells*beats; i++) {
        apds[i] = 0.1;
    }
    for (i = 0; i < 2*ncells*beats; i++) {
        kis[i] = 0;
        nais[i] = 0;
        cais[i] = 0;
        casrs[i] = 0;
    }
    for (i = 0; i < ncells*beats; i++) {
        cpeaks[i] = 0;
    }

    for (i = 0; i < ncells; i++) {
        pcls[i] = 0.0;
        stimtime[i] = 0.0;
        startapds[i] = 0.0;
        st[i] = 0.0;
        curbeat[i] = -1;
    }
    donebifurcation = 0;
}

template <template<int> class typecell, int ncells, int beats>
void APDBifurcation<typecell, ncells, beats>::iterate(const int id, long double dt, long double t) {
    if (id < ncells && curbeat[id] < beats) {
        PRECTYPE vold = Cells.v[id];
        st[id] = (t > -dt/4.0 && t > stimtime[id] - dt/4.0 && t < stimtime[id] + stimduration - dt/4.0) ? stimulus : 0.0;
        Cells.stepdt(id, (PRECTYPE)dt, st[id]);

        PRECTYPE vnew = Cells.v[id];
        if (vnew >= threshold) {
            if (Cells.cai[id] > cpeaks[beats*id + curbeat[id] + 1]) {
                cpeaks[beats*id + curbeat[id] + 1] = Cells.cai[id];
            }
        }
        if (vold < threshold && vnew >= threshold) {
            startapds[id] = t;
            if (curbeat[id] >= -1 && curbeat[id] < beats - 1) {
                kis[2*(beats*id + curbeat[id] + 1)] = Cells.ki[id];
                nais[2*(beats*id + curbeat[id] + 1)] = Cells.nai[id];
                cais[2*(beats*id + curbeat[id] + 1)] = Cells.cai[id];
                casrs[2*(beats*id + curbeat[id] + 1)] = Cells.casr[id];
            }
        }
        else if (vold >= threshold && vnew < threshold) {
            curbeat[id] += 1;
            stimtime[id] += pcls[id];
            while (stimtime[id] < t - dt/4.0) {
               stimtime[id] += pcls[id];
            }
            kis[2*(beats*id + curbeat[id]) + 1] = Cells.ki[id];
            nais[2*(beats*id + curbeat[id]) + 1] = Cells.nai[id];
            cais[2*(beats*id + curbeat[id]) + 1] = Cells.cai[id];
            casrs[2*(beats*id + curbeat[id]) + 1] = Cells.casr[id];
            
            apds[beats*id + curbeat[id]] = t - startapds[id];
            if (curbeat[id] == beats - 1) {
                donebifurcation = donebifurcation + 1;
                printf("%d\n",donebifurcation);
            }
            
        }
    }
}

template <template<int> class typecell, int ncells, int beats>
void APDBifurcation<typecell, ncells, beats>::iterateall(long double dt, long double t) {
    for (int i = 0; i < ncells; i++) {
        iterate(i, dt, t);
    }
}

template <template<int> class typecell, int ncells, int beats>
void APDBifurcation<typecell, ncells, beats>::dobif(long double dt, double startt, FILE* ap) {
    long double t = startt;
    long double t_save = 0;
    double dum;
    long double t_saveap = REMOVEBEATS*pcls[0] - 100;
    while (donebifurcation < ncells) {
        iterateall(dt, t);
        t = t + dt;
        if ((t > stimtime[0] - dt/4.0 && t < stimtime[0] + 3*dt/4.0)) {
            printf("%g\t%d\t%g\t%g\t%g\t%g\n",(double)t, curbeat[0],Cells.v[0],Cells.cai[0],Cells.nai[0], apds[curbeat[0]]);
            t_save = t_save + pcls[0];
        }
        if (t > t_saveap - dt/4.0) {
            fprintf(ap,"%g\t%g\t%g\t%g\n",(double)t - REMOVEBEATS*pcls[0], Cells.v[0],Cells.cai[0], Cells.nai[0]);
            t_saveap = t_saveap + 1;
        }
    }
}
#endif
