/*==============================================================================
 MODULE: mgpt_diffeqs.c			[mgpt]
 Written by: Mario A. Rodriguez-Meza
 Starting date:	January 2018
 Purpose:
 Language: C
 Use:
 Routines and functions:
 External modules, routines and headers:
 Comments and notes:
 Info: Mario A. Rodriguez-Meza
 Depto. de Fisica, ININ
 Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
 e-mail: marioalberto.rodriguez@inin.gob.mx
 http://www.inin.gob.mx
 
 Major revisions:
 Copyright: (c) 2005-2018 Mar.  All Rights Reserved
 ================================================================================
 Legal matters:
 The author does not warrant that the program and routines it contains
 listed below are free from error or suitable for particular applications,
 and he disclaims all liability from any consequences arising from their	use.
 ==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"

local void integration(double ystart[], int nvar, double x1, double x2, double eps, double h1,
                       double hmin, int *nok, int *nbad, int maxnsteps,
                       void (*derivsin)(double, double [], double []),
                       void (*jacobnin)(double, double [], double [], double **, int));

local void derivsFirstOrder(double x,double y[],double dydx[]);
local void derivsFirstOrder_LCDM(double x,double y[],double dydx[]);

local void derivsSecondOrder(double x,double y[],double dydx[]);
local void derivsThirdOrder(double eta,double y[],double dydx[]);

local void jacobnFirstOrder(double x, double y[], double dfdx[], double **dfdy, int n);
local void jacobnFirstOrder_LCDM(double x, double y[], double dfdx[], double **dfdy, int n);
local void jacobnSecondOrder(double x, double y[], double dfdx[], double **dfdy, int n);
local void jacobnThirdOrder(double x, double y[], double dfdx[], double **dfdy, int n);

local void derivsFirstOrder(double x,double y[],double dydx[])
{
    nrhs++;

    dydx[1] = y[2];
    dydx[2]=f1(x)*mu(x,gd.kf)*y[1]-f2(x)*y[2];
}

local void jacobnFirstOrder(double x, double y[], double dfdx[], double **dfdy, int n)
{
    int i;
    
    for (i=1;i<=n;i++) dfdx[i]=0.0;
    dfdy[1][2]=1.0;
    dfdy[2][1]=f1(x)*mu(x,gd.kf);
    dfdy[2][2] = -f2(x);
}

local void derivsFirstOrder_LCDM(double x,double y[],double dydx[])
{
    nrhs++;

    dydx[1] = y[2];
    dydx[2]=f1(x)*y[1]-f2(x)*y[2];
}

local void jacobnFirstOrder_LCDM(double x, double y[], double dfdx[], double **dfdy, int n)
{
    int i;
    
    for (i=1;i<=n;i++) dfdx[i]=0.0;
    dfdy[1][2]=1.0;
    dfdy[2][1]=f1(x);
    dfdy[2][2] = -f2(x);
}

local void derivsSecondOrder(double x,double y[],double dydx[])
{
    nrhs++;
    
    dydx[1] = y[2];
    dydx[2]=f1(x)*mu(x,gd.k1)*y[1]-f2(x)*y[2];
//
    dydx[3] = y[4];
    dydx[4]=f1(x)*mu(x,gd.k2)*y[3]-f2(x)*y[4];
//
    dydx[5] = y[6];
    dydx[6]=f1(x)*mu(x,gd.kf)*y[5]-f2(x)*y[6]
    + sourceA(x,gd.kf,gd.k1,gd.k2)*y[1]*y[3];
//
    dydx[7] = y[8];
    dydx[8]=f1(x)*mu(x,gd.kf)*y[7]-f2(x)*y[8]
    + sourceb(x,gd.kf,gd.k1,gd.k2)*y[1]*y[3];
}

local void jacobnSecondOrder(double x, double y[], double dfdx[], double **dfdy, int n)
{
    int i;
    
    for (i=1;i<=n;i++) dfdx[i]=0.0;
    dfdy[1][2]=1.0;
    dfdy[2][1]=f1(x)*mu(x,gd.k1);
    dfdy[2][2] = -f2(x);
    dfdy[3][4] = 1.0;
    dfdy[4][3] = f1(x)*mu(x,gd.k2);
    dfdy[4][4] = -f2(x);
    dfdy[5][6] = 1.0;
    dfdy[6][1] = sourceA(x,gd.kf,gd.k1,gd.k2)*y[3];
    dfdy[6][3] = sourceA(x,gd.kf,gd.k1,gd.k2)*y[1];
    dfdy[6][5] = f1(x)*mu(x,gd.kf);
    dfdy[6][6] = -f2(x);
    dfdy[7][8] = 1.0;
    dfdy[8][1] = sourceb(x,gd.kf,gd.k1,gd.k2)*y[3];
    dfdy[8][3] = sourceb(x,gd.kf,gd.k1,gd.k2)*y[1];
    dfdy[8][7] = f1(x)*mu(x,gd.kf);
    dfdy[8][8] = -f2(x);
}

local void derivsThirdOrder(double eta,double y[],double dydx[])
{
    nrhs++;
    real Dpk, Dpp, D2f, D2mf, kplusp, kpluspm;

    kplusp = kpp(gd.x,gd.k,gd.p);
    kpluspm = kpp(-gd.x,gd.k,gd.p);
    Dpk  = y[1];
    Dpp  = y[3];
    D2f  = y[5];
    D2mf = y[7];

    dydx[1] = y[2];
    dydx[2] = f1(eta)*mu(eta,gd.k)*y[1]-f2(eta)*y[2];
//
    dydx[3] = y[4];
    dydx[4] = f1(eta)*mu(eta,gd.p)*y[3]-f2(eta)*y[4];
//
    dydx[5] = y[6];
    dydx[6] = f1(eta)*mu(eta,kplusp)*y[5]-f2(eta)*y[6] //;
    + SD2(eta,gd.x,gd.k,gd.p)*y[1]*y[3];
//
    dydx[7] = y[8];
    dydx[8] = f1(eta)*mu(eta,kpluspm)*y[7]-f2(eta)*y[8] //;
    + SD2(eta,-gd.x,gd.k,gd.p)*y[1]*y[3];
//
    dydx[9] = y[10];
    dydx[10] = f1(eta)*mu(eta,gd.k)*y[9]-f2(eta)*y[10]
    + S3I( eta,gd.x,gd.k,gd.p,Dpk,Dpp,D2f,D2mf)
    + S3II(eta,gd.x,gd.k,gd.p,Dpk,Dpp,D2f,D2mf)
    + S3FL(eta,gd.x,gd.k,gd.p,Dpk,Dpp,D2f,D2mf)
    + S3dI(eta,gd.x,gd.k,gd.p,Dpk,Dpp,D2f,D2mf);
}

local void jacobnThirdOrder(double x, double y[], double dfdx[], double **dfdy, int n)
{
    int i;
    real Dpk, Dpp, D2f, D2mf, kplusp, kpluspm;
    
    kplusp = kpp(gd.x,gd.k,gd.p);
    kpluspm = kpp(-gd.x,gd.k,gd.p);

    for (i=1;i<=n;i++) dfdx[i]=0.0;
    dfdy[1][2]=1.0;
    dfdy[2][1]=f1(x)*mu(x,gd.k);
    dfdy[2][2] = -f2(x);
    dfdy[3][4] = 1.0;
    dfdy[4][3] = f1(x)*mu(x,gd.p);
    dfdy[4][4] = -f2(x);
    dfdy[5][6] = 1.0;
    dfdy[6][1] = SD2(x,gd.x,gd.k,gd.p)*y[3];
    dfdy[6][3] = SD2(x,gd.x,gd.k,gd.p)*y[1];
    dfdy[6][5] = f1(x)*mu(x,kplusp);
    dfdy[6][6] = -f2(x);
    dfdy[7][8] = 1.0;
    dfdy[8][1] = SD2(x,-gd.x,gd.k,gd.p)*y[3];
    dfdy[8][3] = SD2(x,-gd.x,gd.k,gd.p)*y[1];
    dfdy[8][7] = f1(x)*mu(x,kpluspm);
    dfdy[8][8] = -f2(x);
    dfdy[9][10] = 1.0;
    dfdy[10][9] = f1(x)*mu(x,gd.k);
    dfdy[10][10] = -f2(x);
}

global real DpFunction(real k)
{
    int nbad,nok;
    double *ystart;
    real Dptmp;

    gd.kf = k;
    
    ystart=dvector(1,NEQS1Order);
    xp=dvector(1,200);
    yp=dmatrix(1,NEQS1Order,1,200);

    ystart[1]=rexp(gd.xnow);
    ystart[2]=rexp(gd.xnow);

    nrhs=0;
    kmax=100;
    dxsav=(gd.xstop-gd.xnow)/20.0;
//
    integration(ystart,NEQS1Order,gd.xnow,gd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                 &nok,&nbad,cmd.maxnsteps,derivsFirstOrder,jacobnFirstOrder);

    Dptmp = yp[1][kount];
    
    free_dmatrix(yp,1,NEQS1Order,1,200);
    free_dvector(xp,1,200);
    free_dvector(ystart,1,NEQS1Order);

    return (Dptmp);
}

global real DpFunction_LCDM(real k)
{
    int nbad,nok;
    double *ystart;
    real Dptmp;
    
    gd.kf = k;
    
    ystart=dvector(1,NEQS1Order);
    xp=dvector(1,200);
    yp=dmatrix(1,NEQS1Order,1,200);
    
    ystart[1]=rexp(gd.xnow);
    ystart[2]=rexp(gd.xnow);
    
    nrhs=0;
    kmax=100;
    dxsav=(gd.xstop-gd.xnow)/20.0;
//
    integration(ystart,NEQS1Order,gd.xnow,gd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                &nok,&nbad,cmd.maxnsteps,derivsFirstOrder_LCDM,jacobnFirstOrder_LCDM);
    
    Dptmp = yp[1][kount];
    
    free_dmatrix(yp,1,NEQS1Order,1,200);
    free_dvector(xp,1,200);
    free_dvector(ystart,1,NEQS1Order);
    
    return (Dptmp);
}


global global_D2v2_ptr DsSecondOrder_func(real kf, real k1, real k2)
{
    int nbad,nok;
    double *ystart;
    global_D2v2_ptr ptmp;
//
    gd.kf = kf;
    gd.k1 = k1;
    gd.k2 = k2;
//
    ystart=dvector(1,NEQS2Orderv2);
    xp=dvector(1,200);
    yp=dmatrix(1,NEQS2Orderv2,1,200);
    
    ystart[1]=rexp(gd.xnow);
    ystart[2]=rexp(gd.xnow);
    ystart[3]=rexp(gd.xnow);
    ystart[4]=rexp(gd.xnow);
    ystart[5]=3.0*rexp(2.0*gd.xnow)/7.0;
    ystart[6]=6.0*rexp(2.0*gd.xnow)/7.0;
    ystart[7]=3.0*rexp(2.0*gd.xnow)/7.0;
    ystart[8]=6.0*rexp(2.0*gd.xnow)/7.0;
    
    nrhs=0;
    kmax=100;
    dxsav=(gd.xstop-gd.xnow)/20.0;

    integration(ystart,NEQS2Orderv2,gd.xnow,gd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                &nok,&nbad,cmd.maxnsteps,derivsSecondOrder,jacobnSecondOrder);
    
    ptmp = (global_D2v2_ptr) allocate(1 * sizeof(global_D2v2));

    etaD2v2(ptmp) = xp[kount];
    Dpk1D2v2(ptmp) = yp[1][kount];
    Dpk2D2v2(ptmp) = yp[3][kount];
    DA2D2(ptmp) = yp[5][kount];
    DB2D2(ptmp) = yp[7][kount];
//
    free_dmatrix(yp,1,NEQS2Orderv2,1,200);
    free_dvector(xp,1,200);
    free_dvector(ystart,1,NEQS2Orderv2);
    
    return ptmp;
}

global global_D3v2_ptr DsThirdOrder_func(real x, real k, real p)
{
    int nbad,nok;
    double *ystart;
    global_D3v2_ptr ptmp;

    gd.x = x;
    gd.k = k;
    gd.p = p;
//
    ystart=dvector(1,NEQS3Orderv2);
    xp=dvector(1,200);
    yp=dmatrix(1,NEQS3Orderv2,1,200);
//
    ystart[1]=rexp(gd.xnow);
    ystart[2]=rexp(gd.xnow);
    ystart[3]=rexp(gd.xnow);
    ystart[4]=rexp(gd.xnow);
//
    ystart[5]=(3.0/7.0)*rexp(2.0*gd.xnow)*(1.0-rsqr(x));
    ystart[6]=(6.0/7.0)*rexp(2.0*gd.xnow)*(1.0-rsqr(x));
    ystart[7]=(3.0/7.0)*rexp(2.0*gd.xnow)*(1.0-rsqr(x));
    ystart[8]=(6.0/7.0)*rexp(2.0*gd.xnow)*(1.0-rsqr(x));
//
    ystart[9]=(5.0/(7.0*9.0))*rexp(3.0*gd.xnow)*rsqr(1.0-rsqr(x))
        *(
            1.0/(1.0 + rsqr(p/k) + 2.0 * (p/k) * x)
          + 1.0/(1.0 + rsqr(p/k) - 2.0 * (p/k) * x)
        );
//
    ystart[10]=(15.0/(7.0*9.0))*rexp(3.0*gd.xnow)*rsqr(1.0-rsqr(x))
        *(
            1.0/(1.0 + rsqr(p/k) + 2.0 * (p/k) * x)
          + 1.0/(1.0 + rsqr(p/k) - 2.0 * (p/k) * x)
          );
//
    nrhs=0;
    kmax=100;
    dxsav=(gd.xstop-gd.xnow)/20.0;
    integration(ystart,NEQS3Orderv2,gd.xnow,gd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                &nok,&nbad,cmd.maxnsteps,derivsThirdOrder,jacobnThirdOrder);
// CHOOSE MANUALLY THE SOLVER:
//    odeint(ystart,NEQS3Orderv2,gd.xnow,gd.xstop,cmd.eps,gd.dx,cmd.dxmin,
//           &nok,&nbad,cmd.maxnsteps,derivsThirdOrder,jacobnThirdOrder,bsstep);


    ptmp = (global_D3v2_ptr) allocate(1 * sizeof(global_D3v2));

    etaD3v2(ptmp) = xp[kount];
    DpkD3v2(ptmp) = yp[1][kount];
    DppD3v2(ptmp) = yp[3][kount];
    D2fD3v2(ptmp) = yp[5][kount];
    D2mfD3v2(ptmp) = yp[7][kount];
    D3symmD3v2(ptmp) = yp[9][kount];
//
    free_dmatrix(yp,1,NEQS3Orderv2,1,200);
    free_dvector(xp,1,200);
    free_dvector(ystart,1,NEQS3Orderv2);
    
    return ptmp;
}

#define BSSTEP          0
#define NULLMETHOD      1
#define RKQS            2
#define STIFBS          3
#define STIFF          4

local void integration(double ystart[], int nvar, double x1, double x2, double eps, double h1,
                       double hmin, int *nok, int *nbad, int maxnsteps,
                       void (*derivsin)(double, double [], double []),
                       void (*jacobnin)(double, double [], double [], double **, int))
{
    switch (gd.method_int) {
        case BSSTEP:
            odeint(ystart,nvar,gd.xnow,gd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                   nok,nbad,maxnsteps,derivsin,jacobnin,bsstep);
            break;
//
        case RKQS:
            odeint(ystart,nvar,gd.xnow,gd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                   nok,nbad,maxnsteps,derivsin,jacobnin,rkqs);
            break;
//
        case STIFBS:
            odeint(ystart,nvar,gd.xnow,gd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                   nok,nbad,maxnsteps,derivsin,jacobnin,stifbs);
            break;
//
        case STIFF:
            odeint(ystart,nvar,gd.xnow,gd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                   nok,nbad,maxnsteps,derivsin,jacobnin,stiff);
            break;
//
        case NULLMETHOD:
            odeint(ystart,nvar,gd.xnow,gd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                   nok,nbad,maxnsteps,derivsin,jacobnin,bsstep);
            break;
//
        default:
            odeint(ystart,nvar,gd.xnow,gd.xstop,cmd.eps,gd.dx,cmd.dxmin,
                   nok,nbad,maxnsteps,derivsin,jacobnin,bsstep);
            break;
    }
}

void integration_method_string_to_int(string method_str,int *method_int)
{
    *method_int=-1;
    if (strcmp(method_str,"bsstep") == 0) {
        *method_int = BSSTEP;
        strcpy(gd.integration_method_comment, "bsstep integration method");
    }
//
    if (strcmp(method_str,"rkqs") == 0) {
        *method_int = RKQS;
        strcpy(gd.integration_method_comment, "rkqs integration method");
    }
//
    if (strcmp(method_str,"stifbs") == 0) {
        *method_int = STIFBS;
        strcpy(gd.integration_method_comment, "stifbs integration method");
    }
//
    if (strcmp(method_str,"stiff") == 0) {
        *method_int = STIFF;
        strcpy(gd.integration_method_comment, "stiff integration method");
    }
//
    if (strnull(method_str)) {
        *method_int = NULLMETHOD;
        strcpy(gd.integration_method_comment,
               "null integration method ... running deafult (bsstep)");
        fprintf(stdout,"\n\tintegration: default integration method (bsstep)...\n");
    }
//
    if (*method_int == -1) {
        *method_int = BSSTEP;
        strcpy(gd.integration_method_comment,
               "Unknown integration method ... running deafult (bsstep)");
        fprintf(stdout,"\n\tintegration: Unknown method... %s ",cmd.integration_method);
        fprintf(stdout,
                "\n\trunning default integration method (bsstep)...\n");
    }
}

#undef BSSTEP
#undef RKQS
#undef STIFBS
#undef STIFF
#undef NULLMETHOD






