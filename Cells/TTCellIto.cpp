//
//  TTCellIto.cpp
//
//  Implementation of the ten Tusscher model, with UCLA Ito
//  Created by Julian Landaw on 12/25/16.
//  Copyright © 2016 Julian Landaw. All rights reserved.
//
#include <math.h>
#include <fstream>
#include "TTCellIto.h"

#ifndef TTCellIto_cpp
#define TTCellIto_cpp

#define ko 5.4
#define cao 2.0
#define nao 140.0
#define vc 0.016404
#define vsr 0.001094
#define bufc 0.15
#define kbufc 0.001
#define bufsr 10.0
#define kbufsr 0.3
#define taufca 2.0
#define taug 2.0
#define vmaxup 0.000425
#define kup 0.00025
#define capacitance 0.185
#define gitodv 0.2
#define gkr 0.096
#define gks 0.245
#define pkna 0.03
#define gk1 5.405
#define gna 14.838
#define gbna 0.00029
#define kmk 1.0
#define kmna 40.0
#define knak 1.362
#define gcal 0.000175
#define gbca 0.000592
#define knaca 1000.0
#define kmnai 87.5
#define kmca 1.38
#define ksat 0.1
#define gamma 0.35
#define alpha 2.5
#define gpca 0.825
#define kpca 0.0005
#define gpk 0.0146
#define arel 0.016464
#define brelsq 0.0625
#define crel 0.008232
#define vleak 0.00008

#define R 8314.472
#define frdy 96485.3415
#define temp 310.0
#define zna 1.0
#define zk 1.0
#define zca 2.0
#define frt 0.037433890822745
//F/(RT)
#define rtf 26.713760659695648
//(RT)/F

//SK-current variables
#define gsk 0.005
//#define skh 0.0006
//#define skn 4

#ifndef DV_MAX
#define DV_MAX 0.1
#endif

#ifndef ADAPTIVE
#define ADAPTIVE 10
#endif

#ifndef stimulus
#define stimulus -52.0
#endif 

#ifndef stimduration
#define stimduration 1.0
#endif

template <int ncells>
TTCellIto<ncells>::TTCellIto()
{
    for (int i = 0; i < ncells; i++) {
        v[i] = -86.2;
        nai[i] = 11.6;
        ki[i] = 138.3;
        cai[i] = 0.0002;
        casr[i] = 0.2;
        m[i] = 0.0;
        h[i] = 0.75;
        j[i] = 0.75;
        xr1[i] = 0.0;
        xr2[i] = 1.0;
        xs[i] = 0.0;
        d[i] = 0.0;
        f[i] = 1.0;
        fca[i] = 1.0;
        g[i] = 1.0;
        zdv[i] = 0.0;
        ydv[i] = 1.0;
        diffcurrent[i] = 0.0;
        itofac[i] = 0.0;
        nacafac[i] = 1.0;
        nakfac[i] = 1.0;
        ibarcafac[i] = 1.0;
        ikrfac[i] = 1.0;
        iksfac[i] = 1.0;
	ik1fac[i] = 1.0;
        
        // SK Variables
        iskfac[i] = 0.0;
        skh[i] = 0.008;
        skn[i] = 4;
        
        #ifndef UCLAito
        rinfshift[i] = 0.0;
        sinfshift[i] = 0.0;
        taurshift[i] = 0.0;
        tausshift[i] = 0.0;
        #endif
        
        // Na and K clamp variables;
        #ifdef clampnai
        naiclamped[i] = true;
        #else
        naiclamped[i] = false;
        #endif
        #ifdef clampki
        kiclamped[i] = true;
        #else
        kiclamped[i] = false;
        #endif
        
        vcfac[i] = 1.0;
    }
}

template <int ncells>
bool TTCellIto<ncells>::iterate(const int id, double dt, double st, double dv_max) {
    double dm, dh, dj, dd, df, dfca, dzdv, dydv, dxs, dxr1, dxr2, dg;
    double dcai, dcasr, dnai, dki, dv;
    double ina, ical, ito, iks, isk, ikr, ik1, inaca, inak, ipca, ipk, ibna, ibca;
    
    ina = comp_ina(id, dt, dm, dh, dj);         // fast sodium current
    ical = comp_ical(id, dt, dd, df, dfca);     // window current
    ito = comp_ito(id, dt, dzdv, dydv);         // UCLA Ito current
    isk = comp_isk(id);                         // SK current
    iks = comp_iks(id, dt, dxs);                // Slow rectifying K current
    ikr = comp_ikr(id, dt, dxr1, dxr2);         // Rapid rectifying K current
    ik1 = comp_ik1(id);                         // IK current
    inaca = comp_inaca(id);                     // Sodium-Ca exchanger
    inak = comp_inak(id);                       // Na-K pump
    ipk = comp_ipk(id);                         // Background K current
    ipca = comp_ipca(id);                       // Background Ca current
    ibca = comp_ibca(id);                       // ... Another background Ca current
    ibna = comp_ibna(id);                       // Background Na current
    dv = (diffcurrent[id] - (ikr + iks + isk + ik1 + ito + ina + ibna + ical + ibca + inak + inaca + ipca + ipk + st))*dt;
    
    // Comment out if one would like to avoid using a dynamic step size that checks for dv being too big
    if (dv_max > 0 && dv*dv > dv_max*dv_max) {return false;}
    
    comp_calcdyn(id, dt, ical, ibca, ipca, inaca, dg, dcasr, dcai);

    if (!naiclamped[id]) {
        dnai = -(ina + ibna + 3.0*inak + 3.0*inaca)/((vcfac[id]*vc)*frdy)*capacitance;
        nai[id] += dnai*dt;
    }
    if (!kiclamped[id]) {
        dki = -(ik1 + ito + isk + ikr + iks - 2.0*inak + ipk + st - diffcurrent[id])/((vcfac[id]*vc)*frdy)*capacitance;
        ki[id] += dki*dt;
    }
    
    v[id] += dv;
    m[id] += dm*dt;
    h[id] += dh*dt;
    j[id] += dj*dt;
    d[id] += dd*dt;
    f[id] += df*dt;
    fca[id] += dfca*dt;
    zdv[id] += dzdv*dt;
    ydv[id] += dydv*dt;
    xs[id] += dxs*dt;
    xr1[id] += dxr1*dt;
    xr2[id] += dxr2*dt;
    g[id] += dg*dt;
    cai[id] += dcai*dt;
    casr[id] += dcasr*dt;
    //nai[id] += dnai*dt;
    //ki[id] += dki*dt;
    
    return true;
}

template <int ncells>
void TTCellIto<ncells>::stepdt (const int id, double dt, double st) {
    if (id > -1 && id < ncells) {
        bool success = iterate(id, dt, st, DV_MAX);
        if (!success) {
            for (int i = 0; i < ADAPTIVE; i++) {
                iterate(id, dt/ADAPTIVE, st, -1.0);
            }
        }
    }
}

#ifndef gnaped
#define gnaped 0
#endif

template <int ncells>
double TTCellIto<ncells>::comp_ina (int id, double dt, double& dm, double& dh, double& dj) //Fast Sodium Current
{
    double ena, mtau, htau, jtau, mss, hss, ina;
    
    ena = (rtf/zna)*log(nao/nai[id]);
    
    if (v[id] < -40.0) {
        htau = 1.0/(0.057*exp(-(v[id]+80)/6.8) + 2.7*exp(0.079*v[id]) + 3.1e5*exp(0.3485*v[id])); //1.0/(ah + bh);
        jtau = 1.0/((-25428.0*exp(0.2444*v[id]) - 6.948e-6*exp(-0.04391*v[id]))*(v[id]+37.78)/(1.0 + exp(0.311*(v[id]+79.23))) + 0.02424*exp(-0.01052*v[id])/(1.0 + exp(-0.1378*(v[id]+40.14)))); // 1.0/(aj + bj);
    }
    else {
        htau = (0.13*(1.0 + exp(-(v[id]+10.66)/11.1)))/0.77;
        jtau = (1.0 + exp(-0.1*(v[id]+32.0)))/(0.6*exp(0.057*v[id]));
    }
    
    mtau = (1.0/(1.0 + exp((-60.0-v[id])/5.0)))*(0.1/(1.0 + exp((v[id]+35.0)/5.0)) + 0.1/(1.0 + exp((v[id]-50.0)/200.0)));
    
    mss = 1.0/((1.0 + exp((-56.86-v[id])/9.03))*(1.0 + exp((-56.86-v[id])/9.03)));
    hss = 1.0/((1.0 + exp((v[id]+71.55)/7.43))*(1.0 + exp((v[id]+71.55)/7.43)));
    
    hss = gnaped + (1 - gnaped)*hss;
    //jss = 0.04 + (1 - 0.04)*jss;
    //jss = hss;
    
    ina = gna*m[id]*m[id]*m[id]*h[id]*j[id]*(v[id]-ena);
    
    dm = (mss-(mss-m[id])*exp(-dt/mtau) - m[id])/dt;
    dh = (hss-(hss-h[id])*exp(-dt/htau) - h[id])/dt;
    dj = (hss-(hss-j[id])*exp(-dt/jtau) - j[id])/dt;
    
    return ina;
}

template <int ncells>
double TTCellIto<ncells>::comp_ical (int id, double dt, double& dd, double& df, double& dfca) // L-type Calcium Current
{
#ifdef xiaodong    
    double taud, tauf, dss, fss, fcass, ical, vfrt;
    dss = 1.0/(1.0+exp(-(v[id]+10.0)/6.24)); 
    taud = dss*1.0/(6.24*0.035);
    if (fabs(v[id] + 10.0) > 0.001) {
        taud = dss*(1.0-exp(-(v[id]+10.0)/6.24))/(0.035*(v[id]+10.0)); 
    }

    fss = (1.0/(1.0+exp((v[id]+32.0)/8.0)))+(0.6/(1.0+exp((50.0-v[id])/20.0))); 
    tauf = 1.0/(0.0197*exp(-(0.0337*(v[id]+10.0)*0.0337*(v[id]+10.0))) + 0.02);
    
    dd = (dss-(dss-d[id])*exp(-dt/taud) - d[id])/dt;
    df = (fss-(fss-f[id])*exp(-dt/tauf) - f[id])/dt;
    dfca = 0.0;
    
    fcass = 1.0;
    if (cai[id] > 0) {
        fcass = 1.0/(1.0 + cai[id]/0.0006);
    }
    vfrt = v[id]*zca*frt;
    
    if (fabs(vfrt) < 1e-3) {
        ical = ibarcafac[id]*gcal*d[id]*f[id]*fcass*zca*frdy*(cai[id]*exp(vfrt) - 0.341*cao)/(1.0 + vfrt/2.0 + vfrt*vfrt/6.0 + vfrt*vfrt*vfrt/24.0);
    }
    else {
        ical = ibarcafac[id]*gcal*d[id]*f[id]*fcass*zca*frdy*vfrt*(cai[id]*exp(vfrt) - 0.341*cao)/(exp(vfrt) - 1.0);
    }
#else    
    
    double taud, tauf, dss, fss, fcass, ical, vfrt;
    
    taud = (1.4/(1.0 + exp((-35.0-v[id])/13.0)) + 0.25)*(1.4/(1.0 + exp((v[id]+5.0)/5.0))) + 1.0/(1.0 + exp((50.0-v[id])/20.0));
    tauf = 1125.0*exp(-(v[id]+27.0)*(v[id]+27.0)/240.0) + 165.0/(1.0 + exp((25-v[id])/10.0)) + 80.0;
    dss = 1.0/(1.0 + exp((-5.0-v[id])/7.5));
    fss = 1.0/(1.0 + exp((v[id]+20.0)/7.0));
    fcass = (1.0/(1.0 + pow(cai[id]/0.000325,8.0)) + 0.1/(1.0 + exp((cai[id]-0.0005)/0.0001)) + 0.2/(1.0 + exp((cai[id]-0.00075)/0.0008)) + 0.23)/1.46;
    
    
    vfrt = v[id]*zca*frt;
    
    if (fabs(vfrt) < 1e-3) {
        ical = ibarcafac[id]*gcal*d[id]*f[id]*fca[id]*zca*frdy*(cai[id]*exp(vfrt) - 0.341*cao)/(1.0 + vfrt/2.0 + vfrt*vfrt/6.0 + vfrt*vfrt*vfrt/24.0);
    }
    else {
        ical = ibarcafac[id]*gcal*d[id]*f[id]*fca[id]*zca*frdy*vfrt*(cai[id]*exp(vfrt) - 0.341*cao)/(exp(vfrt) - 1.0);
    }
    
    dd = (dss-(dss-d[id])*exp(-dt/taud) - d[id])/dt;
    df = (fss-(fss-f[id])*exp(-dt/tauf) - f[id])/dt;
    dfca = (fcass > fca[id] && v[id] > -60.0) ? 0.0 : (fcass - (fcass - fca[id])*exp(-dt/taufca) - fca[id])/dt;
#endif    
    return ical;
}

template <int ncells>
double TTCellIto<ncells>::comp_ito (int id, double dt, double& dzdv, double& dydv)
{
    #ifdef UCLAito
    double ek, zssdv, tauzdv, yssdv, tauydv, ito;
    
    ek = (rtf/zk)*log(ko/ki[id]);
    
    zssdv = 1.0/(1.0 + exp(-(v[id] + 3.0)/15.0));
    tauzdv = 3.5*exp(-(v[id]/30.0)*(v[id]/30.0))+1.5;
    
    yssdv = 1.0/(1.0+exp((v[id] +33.5)/10.0));
    tauydv = 20.0/(1.0+exp((v[id] +33.5)/10.0)) + 20.0;
    
    ito = itofac[id]*gitodv*zdv[id]*ydv[id]*(v[id]-ek);
    
    dzdv = (zssdv-(zssdv-zdv[id])*exp(-dt/tauzdv) - zdv[id])/dt;
    dydv = (yssdv-(yssdv-ydv[id])*exp(-dt/tauydv) - ydv[id])/dt;
    
    return ito;
    #else
    #ifdef ENDO
    double Gto=0.073;
    #else
    double Gto = 0.294;
    #endif
    
    double R_INF, S_INF, TAU_R, TAU_S, ek, ito;
    
    ek = (rtf/zk)*log(ko/ki[id]);
    
    #ifdef EPI
    R_INF=1./(1.+exp((20-(v[id] - rinfshift[id]))/6.));
    S_INF=1./(1.+exp(((v[id] - sinfshift[id])+20)/5.));
    TAU_R=9.5*exp(-((v[id] - taurshift[id])+40.)*((v[id] - taurshift[id])+40.)/1800.)+0.8;
    TAU_S=85.*exp(-((v[id] - tausshift[id])+45.)*((v[id] - tausshift[id])+45.)/320.)+5./(1.+exp(((v[id] - tausshift[id])-20.)/5.))+3.;
    #endif
    #ifdef ENDO
    R_INF=1./(1.+exp((20-(v[id] - rinfshift[id]))/6.));
    S_INF=1./(1.+exp(((v[id] - sinfshift[id])+28)/5.));
    TAU_R=9.5*exp(-((v[id] - taurshift[id])+40.)*((v[id] - taurshift[id])+40.)/1800.)+0.8;
    TAU_S=1000.*exp(-((v[id] - tausshift[id])+67)*((v[id] - tausshift[id])+67)/1000.)+8.;
    #endif
    #ifdef MCELL
    R_INF=1./(1.+exp((20-(v[id] - rinfshift[id]))/6.));
    S_INF=1./(1.+exp(((v[id] - sinfshift[id])+20)/5.));
    TAU_R=9.5*exp(-((v[id] - taurshift[id])+40.)*((v[id] - taurshift[id])+40.)/1800.)+0.8;
    TAU_S=85.*exp(-((v[id] - tausshift[id])+45.)*((v[id] - tausshift[id])+45.)/320.)+5./(1.+exp(((v[id] - tausshift[id])-20.)/5.))+3.;
    #endif
    
    ito = itofac[id]*Gto*zdv[id]*ydv[id]*(v[id]-ek);
    
    dzdv = (R_INF - (R_INF - zdv[id])*exp(-dt/TAU_R) - zdv[id])/dt;
    dydv = (S_INF - (S_INF - ydv[id])*exp(-dt/TAU_S) - ydv[id])/dt;
    
    return ito;
    #endif
}

template <int ncells>
double TTCellIto<ncells>::comp_isk (int id)
{
    double isk, ek, z;
    z = 1.0/(1.0 + pow(skh[id]/cai[id],skn[id]));
    ek = (rtf/zk)*log(ko/ki[id]);
    
    isk = iskfac[id]*gsk*z*(v[id] - ek);
    
    return isk;
}

template <int ncells>
double TTCellIto<ncells>::comp_iks (int id, double dt, double& dxs)
{
    double eks, xsss, tauxs, iks;
    
    eks = rtf*log((ko + pkna*nao)/(ki[id] + pkna*nai[id]));
    xsss = 1.0/(1.0 + exp((-5.0-v[id])/14.0));
    tauxs = 1100.0/(sqrt(1.0 + exp((-10.0-v[id])/6.0))*(1.0 + exp((v[id]-60.0)/20.0)));
    
    iks = iksfac[id]*gks*xs[id]*xs[id]*(v[id]-eks);
    
    dxs = (xsss - (xsss - xs[id])*exp(-dt/tauxs) - xs[id])/dt;
    
    return iks;
}

template <int ncells>
double TTCellIto<ncells>::comp_ikr (int id, double dt, double& dxr1, double& dxr2)
{
    double ek, axr1, bxr1, axr2, bxr2, tauxr1, tauxr2, xr1ss, xr2ss, ikr;
    
    ek = (rtf/zk)*log(ko/ki[id]);
    
    axr1 = 450.0/(1.0 + exp((-45.0-v[id])/10.0));
    bxr1 = 6.0/(1.0 + exp((v[id]+30.0)/11.5));
    
    axr2 = 3.0/(1.0 + exp((-60.0-v[id])/20.0));
    bxr2 = 1.12/(1.0 + exp((v[id]-60.0)/20.0));
    
    tauxr1 = axr1*bxr1;
    tauxr2 = axr2*bxr2;
    
    xr1ss = 1.0/(1.0 + exp((-26.0-v[id])/7.0));
    xr2ss = 1.0/(1.0 + exp((v[id]+88.0)/24.0));
    
    ikr = ikrfac[id]*gkr*sqrt(ko/5.4)*xr1[id]*xr2[id]*(v[id] - ek);
    
    dxr1 = (xr1ss - (xr1ss - xr1[id])*exp(-dt/tauxr1) - xr1[id])/dt;
    dxr2 = (xr2ss - (xr2ss - xr2[id])*exp(-dt/tauxr2) - xr2[id])/dt;
    
    return ikr;
}

template <int ncells>
double TTCellIto<ncells>::comp_ik1 (int id)
{
    double ek, ak1, bk1, xk1ss, ik1;
    ek = (rtf/zk)*log(ko/ki[id]);
    ak1 = 0.1/(1.0 + exp(0.06*(v[id] - ek - 200.0)));
    bk1 = 3.0*exp(0.0002*(v[id] - ek + 100.0)) + exp(0.1*(v[id] - ek - 10.0))/(1.0 + exp(-0.5*(v[id] - ek)));
    xk1ss = ak1/(ak1 + bk1);
    ik1 = ik1fac[id]*gk1*sqrt(ko/5.4)*xk1ss*(v[id] - ek);
    return ik1;
}

template <int ncells>
double TTCellIto<ncells>::comp_inaca (int id)
{
    double inaca;
    inaca = nacafac[id]*knaca*(exp(gamma*v[id]*frt)*nai[id]*nai[id]*nai[id]*cao - exp((gamma-1.0)*v[id]*frt)*nao*nao*nao*cai[id]*alpha)/((kmnai*kmnai*kmnai + nao*nao*nao)*(kmca + cao)*(1.0 + ksat*exp((gamma-1.0)*v[id]*frt)));
    
    return inaca;
}

template <int ncells>
double TTCellIto<ncells>::comp_inak (int id)
{
    double inak;
    inak = nakfac[id]*knak*ko*nai[id]/((ko+kmk)*(nai[id]+kmna)*(1.0 + 0.1245*exp(-0.1*v[id]*frt) + 0.0353*exp(-v[id]*frt)));
    
    return inak;
}

template <int ncells>
double TTCellIto<ncells>::comp_ipca (int id)
{
    double ipca;
    ipca = gpca*cai[id]/(kpca + cai[id]);
    return ipca;
}

template <int ncells>
double TTCellIto<ncells>::comp_ipk (int id)
{
    double ipk, ek;
    ek = (rtf/zk)*log(ko/ki[id]);
    ipk = gpk*(v[id]-ek)/(1.0 + exp((25.0-v[id])/5.98));
    return ipk;
}

template <int ncells>
double TTCellIto<ncells>::comp_ibna (int id)
{
    double ena, ibna;
    ena = (rtf/zna)*log(nao/nai[id]);
    ibna = gbna*(v[id] - ena);
    return ibna;
}

template <int ncells>
double TTCellIto<ncells>::comp_ibca (int id)
{
    double eca, ibca;
    eca = (rtf/zca)*log(cao/cai[id]);
    ibca = gbca*(v[id] - eca);
    return ibca;
}

template <int ncells>
void TTCellIto<ncells>::comp_calcdyn (int id, double dt, double ical, double ibca, double ipca, double inaca, double &dg, double &dcasr, double &dcai) // MAKE SURE THIS IS COMPUTED AFTER ALL OTHER CURRENTS
{
    double gss, icac, irel, ileak, iup, icasr, casrbuf, bjsr, cjsr, cabuf, bc, cc;
    
    gss = (cai[id] > 0.00035) ? 1.0/(1.0 + pow(cai[id]/0.00035,16.0)) : 1.0/(1.0 + pow(cai[id]/0.00035,6.0));
    icac = -(ical + ibca + ipca - 2.0*inaca)/(zca*(vcfac[id]*vc)*frdy)*capacitance;
    //icac = -(ical + ibca + ipca - 2.0*inaca)/(zca*vc*frdy)*capacitance; // NO CHANGE IN Ca VOLUME
    irel = (arel*casr[id]*casr[id]/(brelsq + casr[id]*casr[id]) + crel)*d[id]*g[id];
    
    dg = (g[id] < gss && v[id] > -60.0) ? 0.0 : (gss - (gss - g[id])*exp(-dt/taug) - g[id])/dt;
    ileak = vleak*(casr[id] - cai[id]);
    iup = vmaxup/(1.0 + (kup/cai[id])*(kup/cai[id]));
    icasr = iup - irel - ileak;
    casrbuf = bufsr*casr[id]/(casr[id]+kbufsr);
    //dcasr = dt*((vcfac[id]*vc)/vsr)*icasr;
    dcasr = dt*(vc/vsr)*icasr; //NO CHANGE IN Ca VOLUME
    bjsr = bufsr - casrbuf - casr[id] - dcasr + kbufsr; //
    cjsr = kbufsr*(casrbuf+dcasr+casr[id]);
    if (cjsr < 0.0) {
        bjsr = bufsr - casrbuf - casr[id]*exp(dcasr/casr[id]) + kbufsr;
        cjsr = kbufsr*(casrbuf + casr[id]*exp(dcasr/casr[id]));
    }
    dcasr = ((sqrt(bjsr*bjsr + 4.0*cjsr)-bjsr)/2.0 - casr[id])/dt;
    cabuf = bufc*cai[id]/(cai[id]+kbufc);
    dcai = dt*(icac - icasr);
    bc = bufc - cabuf - cai[id] - dcai + kbufc;
    cc = kbufc*(cabuf+dcai+cai[id]);
    if (cc < 0.0) {
        bc = bufc - cabuf - cai[id]*exp(dcai/cai[id]) + kbufc;
        cc = kbufc*(cabuf + cai[id]*exp(dcai/cai[id]));
    }
    dcai = ((sqrt(bc*bc+4.0*cc)-bc)/2.0 - cai[id])/dt;
}

template <int ncells>
void TTCellIto<ncells>::setcell (int id, TTCellIto<1>* newcell)
{
    v[id] = newcell->v[0];
    nai[id] = newcell->nai[0];
    ki[id] = newcell->ki[0];
    cai[id] = newcell->cai[0];
    casr[id] = newcell->casr[0];
    m[id] = newcell->m[0];
    h[id] = newcell->h[0];
    j[id] = newcell->j[0];
    xr1[id] = newcell->xr1[0];
    xr2[id] = newcell->xr2[0];
    xs[id] = newcell->xs[0];
    d[id] = newcell->d[0];
    f[id] = newcell->f[0];
    fca[id] = newcell->fca[0];
    g[id] = newcell->g[0];
    zdv[id] = newcell->zdv[0];
    ydv[id] = newcell->ydv[0];
    diffcurrent[id] = newcell->diffcurrent[0];
    itofac[id] = newcell->itofac[0];
    ibarcafac[id] = newcell->ibarcafac[0]; //
    nacafac[id] = newcell->nacafac[0];
    nakfac[id] = newcell->nakfac[0];
    ikrfac[id] = newcell->ikrfac[0];
    iksfac[id] = newcell->iksfac[0];
    ik1fac[id] = newcell->ik1fac[0];
    
    iskfac[id] = newcell->iskfac[0];
    skh[id] = newcell->skh[0];
    skn[id] = newcell->skn[0];
    #ifndef UCLAito
    rinfshift[id] = newcell->rinfshift[0];
    sinfshift[id] = newcell->sinfshift[0];
    taurshift[id] = newcell->taurshift[0];
    tausshift[id] = newcell->tausshift[0];
    #endif
    
    naiclamped[id] = newcell->naiclamped[0];
    kiclamped[id] = newcell->kiclamped[0];
    
    vcfac[id] = newcell->vcfac[0];
}

template <int ncells>
void TTCellIto<ncells>::getcell (int id, TTCellIto<1>* newcell)
{
    newcell->v[0] = v[id];
    newcell->nai[0] = nai[id];
    newcell->ki[0] = ki[id];
    newcell->cai[0] = cai[id];
    newcell->casr[0] = casr[id];
    newcell->m[0] = m[id];
    newcell->h[0] = h[id];
    newcell->j[0] = j[id];
    newcell->xr1[0] = xr1[id];
    newcell->xr2[0] = xr2[id];
    newcell->xs[0] = xs[id];
    newcell->d[0] = d[id];
    newcell->f[0] = f[id];
    newcell->fca[0] = fca[id];
    newcell->g[0] = g[id];
    newcell->zdv[0] = zdv[id];
    newcell->ydv[0] = ydv[id];
    newcell->diffcurrent[0] = diffcurrent[id];
    newcell->itofac[0] = itofac[id];
    newcell->nacafac[0] = nacafac[id];
    newcell->nakfac[0] = nakfac[id];
    newcell->ikrfac[0] = ikrfac[id];
    newcell->iksfac[0] = iksfac[id];
    newcell->ik1fac[0] = ik1fac[id];
    newcell->ibarcafac[0] = ibarcafac[id];
    
    newcell->iskfac[0] = iskfac[id];
    newcell->skh[0] = skh[id];
    newcell->skn[0] = skn[id];
    
    #ifndef UCLAito
    newcell->rinfshift[0] = rinfshift[id];
    newcell->sinfshift[0] = sinfshift[id];
    newcell->taurshift[0] = taurshift[id];
    newcell->tausshift[0] = tausshift[id];
    #endif
    
    newcell->naiclamped[0] = naiclamped[id];
    newcell->kiclamped[0] = kiclamped[id];
    
    newcell->vcfac[0] = vcfac[id];
}

template <int ncells>
void TTCellIto<ncells>::saveconditions(FILE* file, int id, bool header, double t) {
    if (header) {
        fprintf(file,"t\tv\tm\th\tj\td\tf\tfca\tzdv\tydv\txs\txr1\txr2\tg\tcai\tcasr\tnai\tki\n");
    }
    fprintf(file,"%g\t",t);
    fprintf(file,"%.12f\t",v[id]);
    fprintf(file,"%.12f\t",m[id]);
    fprintf(file,"%.12f\t",h[id]);
    fprintf(file,"%.12f\t",j[id]);
    fprintf(file,"%.12f\t",d[id]);
    fprintf(file,"%.12f\t",f[id]);
    fprintf(file,"%.12f\t",fca[id]);
    fprintf(file,"%.12f\t",zdv[id]);
    fprintf(file,"%.12f\t",ydv[id]);
    fprintf(file,"%.12f\t",xs[id]);
    fprintf(file,"%.12f\t",xr1[id]);
    fprintf(file,"%.12f\t",xr2[id]);
    fprintf(file,"%.12f\t",g[id]);
    fprintf(file,"%.12f\t",cai[id]);
    fprintf(file,"%.12f\t",casr[id]);
    fprintf(file,"%.12f\t",nai[id]);
    fprintf(file,"%.12f\n",ki[id]);
}

template <int ncells>
void TTCellIto<ncells>::write(std::fstream& file) {
    file.write((char*)&this, sizeof(this) );
}

template <int ncells>
void TTCellIto<ncells>::read(std::fstream& file) {
    file.read((char*)&this, sizeof(this) );
}
/*
#undef ko 
#undef cao 
#undef nao 
#undef vc 
#undef vsr 
#undef bufc 
#undef kbufc 
#undef bufsr 
#undef kbufsr 
#undef taufca 
#undef taug 
#undef vmaxup 
#undef kup 
#undef capacitance 
#undef gitodv 
#undef gkr 
#undef gks 
#undef pkna 
#undef gk1 
#undef gna 
#undef gbna 
#undef kmk 
#undef kmna 
#undef knak 
#undef gcal 
#undef gbca 
#undef knaca 
#undef kmnai 
#undef kmca 
#undef ksat 
#undef gamma 
#undef alpha 
#undef gpca 
#undef kpca 
#undef gpk 
#undef arel 
#undef brelsq 
#undef crel 
#undef vleak 

#undef R 
#undef frdy 
#undef temp 
#undef zna 
#undef zk 
#undef zca
#undef frt 
//F/(RT)
#undef rtf
//(RT)/F
#undef gsk
*/

#endif // TTCellIto_cpp

