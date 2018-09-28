

#include "../general/stdinc.h"
#include "../math/mathfns.h"
//#include "../general/constant.h"
#include "../math/diffeqs.h"
#include "../math/numrec.h"

void rk4(double y[], double dydx[], int n, double x, double h, double yout[],
               void (*derivs)(double, double [], double []))
{
    int i;
    double xh,hh,h6,*dym,*dyt,*yt;
    
    dym=dvector(1,n);
    dyt=dvector(1,n);
    yt=dvector(1,n);
    hh=h*0.5;
    h6=h/6.0;
    xh=x+hh;
    for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i];
    (*derivs)(xh,yt,dyt);
    for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
    (*derivs)(xh,yt,dym);
    for (i=1;i<=n;i++) {
        yt[i]=y[i]+h*dym[i];
        dym[i] += dyt[i];
    }
    (*derivs)(x+h,yt,dyt);
    for (i=1;i<=n;i++)
        yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
    free_dvector(yt,1,n);
    free_dvector(dyt,1,n);
    free_dvector(dym,1,n);
}



// ODEINT
//
//#include <math.h>
//#define NRANSI
//#include "nrutil_nransi.h"
#define MAXSTP 10000
#define TINY 1.0e-30

//local int kmax,kount;
//local double *xp,**yp,dxsav;


void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
            double hmin, int *nok, int *nbad, int maxnsteps,
            void (*derivsin)(double, double [], double []),
            void (*jacobnin)(double, double [], double [], double **, int),
            void (*rkqsin)(double [], double [], int, double *, double, double, double [],
                         double *, double *, void (*)(double, double [], double []),
                           void (*)(double, double [], double [], double **, int)))
{
    int nstp,i;
    double xsav,x,hnext,hdid,h;
    double *yscal,*y,*dydx;
    
    yscal=dvector(1,nvar);
    y=dvector(1,nvar);
    dydx=dvector(1,nvar);
    x=x1;
    h=SIGN(h1,x2-x1);
    *nok = (*nbad) = kount = 0;
    for (i=1;i<=nvar;i++) y[i]=ystart[i];
    if (kmax > 0) xsav=x-dxsav*2.0;
//    for (nstp=1;nstp<=MAXSTP;nstp++) {
//    for (nstp=1;nstp<=cmd.maxnsteps;nstp++) {
    for (nstp=1;nstp<=maxnsteps;nstp++) {
        (*derivsin)(x,y,dydx);
        for (i=1;i<=nvar;i++)
            yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
        if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
            xp[++kount]=x;
            for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
            xsav=x;
        }
        if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
//        (*rkqsin)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivsSecondOrder);
        (*rkqsin)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivsin,jacobnin);
        if (hdid == h) ++(*nok); else ++(*nbad);
        if ((x-x2)*(x2-x1) >= 0.0) {
            for (i=1;i<=nvar;i++) ystart[i]=y[i];
            if (kmax) {
                xp[++kount]=x;
                for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
            }
            free_dvector(dydx,1,nvar);
            free_dvector(y,1,nvar);
            free_dvector(yscal,1,nvar);
            return;
        }
        if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");
        h=hnext;
    }
    nrerror("Too many steps in routine odeint");
}
#undef MAXSTP
#undef TINY
//#undef NRANSI

// rkqs
//
//#include <math.h>
//#define NRANSI
//#include "nrutil_nransi.h"
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
          double yscal[], double *hdid, double *hnext,
          void (*derivs)(double, double [], double []),
          void (*jacobnin)(double, double [], double [], double **, int))
{
    void rkck(double y[], double dydx[], int n, double x, double h,
              double yout[], double yerr[], void (*derivs)(double, double [], double []));
    int i;
    double errmax,h,htemp,xnew,*yerr,*ytemp;
    
    yerr=dvector(1,n);
    ytemp=dvector(1,n);
    h=htry;
    for (;;) {
        rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
        errmax=0.0;
        for (i=1;i<=n;i++) errmax=DMAX(errmax,fabs(yerr[i]/yscal[i]));
        errmax /= eps;
        if (errmax <= 1.0) break;
        htemp=SAFETY*h*rpow(errmax,PSHRNK);
        h=(h >= 0.0 ? DMAX(htemp,0.1*h) : FMIN(htemp,0.1*h));
        xnew=(*x)+h;
        if (xnew == *x) nrerror("stepsize underflow in rkqs");
    }
    if (errmax > ERRCON) *hnext=SAFETY*h*rpow(errmax,PGROW);
    else *hnext=5.0*h;
    *x += (*hdid=h);
    for (i=1;i<=n;i++) y[i]=ytemp[i];
    free_dvector(ytemp,1,n);
    free_dvector(yerr,1,n);
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
//#undef NRANSI


// rkck
//
//#define NRANSI
//#include "nrutil_nransi.h"

void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
          double yerr[], void (*derivs)(double, double [], double []))
{
    int i;
    static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
    double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
    double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;
    
    ak2=dvector(1,n);
    ak3=dvector(1,n);
    ak4=dvector(1,n);
    ak5=dvector(1,n);
    ak6=dvector(1,n);
    ytemp=dvector(1,n);
    for (i=1;i<=n;i++)
        ytemp[i]=y[i]+b21*h*dydx[i];
    (*derivs)(x+a2*h,ytemp,ak2);
    for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
    (*derivs)(x+a3*h,ytemp,ak3);
    for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
    (*derivs)(x+a4*h,ytemp,ak4);
    for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
    (*derivs)(x+a5*h,ytemp,ak5);
    for (i=1;i<=n;i++)
        ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
    (*derivs)(x+a6*h,ytemp,ak6);
    for (i=1;i<=n;i++)
        yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
    for (i=1;i<=n;i++)
        yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
    free_dvector(ytemp,1,n);
    free_dvector(ak6,1,n);
    free_dvector(ak5,1,n);
    free_dvector(ak4,1,n);
    free_dvector(ak3,1,n);
    free_dvector(ak2,1,n);
}
//#undef NRANSI


//
//#include <math.h>
//#define NRANSI
//#include "nrutil.h"
#define KMAXX 8
#define IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define TINY 1.0e-30
#define SCALMX 0.1

double **d,*x;

void bsstep(double y[], double dydx[], int nv, double *xx, double htry, double eps,
            double yscal[], double *hdid, double *hnext,
            void (*derivs)(double, double [], double []),
            void (*jacobnin)(double, double [], double [], double **, int))
{
    void mmid(double y[], double dydx[], int nvar, double xs, double htot,
              int nstep, double yout[], void (*derivs)(double, double[], double[]));
    void pzextr(int iest, double xest, double yest[], double yz[], double dy[],
                int nv);
    int i,iq,k,kk,km;
    static int first=1,kmax,kopt;
    static double epsold = -1.0,xnew;
    double eps1,errmax,fact,h,red,scale,work,wrkmin,xest;
    double *err,*yerr,*ysav,*yseq;
    static double a[IMAXX+1];
    static double alf[KMAXX+1][KMAXX+1];
    static int nseq[IMAXX+1]={0,2,4,6,8,10,12,14,16,18};
    int reduct,exitflag=0;
    
    d=dmatrix(1,nv,1,KMAXX);
    err=dvector(1,KMAXX);
    x=dvector(1,KMAXX);
    yerr=dvector(1,nv);
    ysav=dvector(1,nv);
    yseq=dvector(1,nv);
    if (eps != epsold) {
        *hnext = xnew = -1.0e29;
        eps1=SAFE1*eps;
        a[1]=nseq[1]+1;
        for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
        for (iq=2;iq<=KMAXX;iq++) {
            for (k=1;k<iq;k++)
                alf[k][iq]=rpow(eps1,(a[k+1]-a[iq+1])/
                               ((a[iq+1]-a[1]+1.0)*(2*k+1)));
        }
        epsold=eps;
        for (kopt=2;kopt<KMAXX;kopt++)
            if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
        kmax=kopt;
    }
    h=htry;
    for (i=1;i<=nv;i++) ysav[i]=y[i];
    if (*xx != xnew || h != (*hnext)) {
        first=1;
        kopt=kmax;
    }
    reduct=0;
    for (;;) {
        for (k=1;k<=kmax;k++) {
            xnew=(*xx)+h;
            if (xnew == (*xx)) nrerror("step size underflow in bsstep");
            mmid(ysav,dydx,nv,*xx,h,nseq[k],yseq,derivs);
            xest=DSQR(h/nseq[k]);
            pzextr(k,xest,yseq,y,yerr,nv);
            if (k != 1) {
                errmax=TINY;
                for (i=1;i<=nv;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
                errmax /= eps;
                km=k-1;
                err[km]=rpow(errmax/SAFE1,1.0/(2*km+1));
            }
            if (k != 1 && (k >= kopt-1 || first)) {
                if (errmax < 1.0) {
                    exitflag=1;
                    break;
                }
                if (k == kmax || k == kopt+1) {
                    red=SAFE2/err[km];
                    break;
                }
                else if (k == kopt && alf[kopt-1][kopt] < err[km]) {
                    red=1.0/err[km];
                    break;
                }
                else if (kopt == kmax && alf[km][kmax-1] < err[km]) {
                    red=alf[km][kmax-1]*SAFE2/err[km];
                    break;
                }
                else if (alf[km][kopt] < err[km]) {
                    red=alf[km][kopt-1]/err[km];
                    break;
                }
            }
        }
        if (exitflag) break;
        red=FMIN(red,REDMIN);
        red=FMAX(red,REDMAX);
        h *= red;
        reduct=1;
    }
    *xx=xnew;
    *hdid=h;
    first=0;
    wrkmin=1.0e35;
    for (kk=1;kk<=km;kk++) {
        fact=FMAX(err[kk],SCALMX);
        work=fact*a[kk+1];
        if (work < wrkmin) {
            scale=fact;
            wrkmin=work;
            kopt=kk+1;
        }
    }
    *hnext=h/scale;
    if (kopt >= k && kopt != kmax && !reduct) {
        fact=FMAX(scale/alf[kopt-1][kopt],SCALMX);
        if (a[kopt+1]*fact <= wrkmin) {
            *hnext=h/fact;
            kopt++;
        }
    }
    free_dvector(yseq,1,nv);
    free_dvector(ysav,1,nv);
    free_dvector(yerr,1,nv);
    free_dvector(x,1,KMAXX);
    free_dvector(err,1,KMAXX);
    free_dmatrix(d,1,nv,1,KMAXX);
}
#undef KMAXX
#undef IMAXX
#undef SAFE1
#undef SAFE2
#undef REDMAX
#undef REDMIN
#undef TINY
#undef SCALMX
//#undef NRANSI


// mmid
//
//#define NRANSI
//#include "nrutil_nransi.h"

void mmid(double y[], double dydx[], int nvar, double xs, double htot, int nstep,
          double yout[], void (*derivs)(double, double[], double[]))
{
    int n,i;
    double x,swap,h2,h,*ym,*yn;
    
    ym=dvector(1,nvar);
    yn=dvector(1,nvar);
    h=htot/nstep;
    for (i=1;i<=nvar;i++) {
        ym[i]=y[i];
        yn[i]=y[i]+h*dydx[i];
    }
    x=xs+h;
    (*derivs)(x,yn,yout);
    h2=2.0*h;
    for (n=2;n<=nstep;n++) {
        for (i=1;i<=nvar;i++) {
            swap=ym[i]+h2*yout[i];
            ym[i]=yn[i];
            yn[i]=swap;
        }
        x += h;
        (*derivs)(x,yn,yout);
    }
    for (i=1;i<=nvar;i++)
        yout[i]=0.5*(ym[i]+yn[i]+h*yout[i]);
    free_dvector(yn,1,nvar);
    free_dvector(ym,1,nvar);
}
//#undef NRANSI

// pzextr
//
//#define NRANSI
//#include "nrutil_nransi.h"

//extern float **d,*x;

void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv)
{
    int k1,j;
    double q,f2,f1,delta,*c;
    
    c=dvector(1,nv);
    x[iest]=xest;
    for (j=1;j<=nv;j++) dy[j]=yz[j]=yest[j];
    if (iest == 1) {
        for (j=1;j<=nv;j++) d[j][1]=yest[j];
    } else {
        for (j=1;j<=nv;j++) c[j]=yest[j];
        for (k1=1;k1<iest;k1++) {
            delta=1.0/(x[iest-k1]-xest);
            f1=xest*delta;
            f2=x[iest-k1]*delta;
            for (j=1;j<=nv;j++) {
                q=d[j][k1];
                d[j][k1]=dy[j];
                delta=c[j]-q;
                dy[j]=f1*delta;
                c[j]=f2*delta;
                yz[j] += dy[j];
            }
        }
        for (j=1;j<=nv;j++) d[j][iest]=dy[j];
    }
    free_dvector(c,1,nv);
}
//#undef NRANSI


//#define NRANSI
//#include "nrutil.h"

void simpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
           double xs, double htot, int nstep, double yout[],
           void (*derivs)(double, double [], double []))
{
    void lubksb(double **a, int n, int *indx, double b[]);
    void ludcmp(double **a, int n, int *indx, double *d);
    int i,j,nn,*indx;
    double d,h,x,**a,*del,*ytemp;
    
    indx=ivector(1,n);
    a=dmatrix(1,n,1,n);
    del=dvector(1,n);
    ytemp=dvector(1,n);
    h=htot/nstep;
    for (i=1;i<=n;i++) {
        for (j=1;j<=n;j++) a[i][j] = -h*dfdy[i][j];
        ++a[i][i];
    }
    ludcmp(a,n,indx,&d);
    for (i=1;i<=n;i++)
        yout[i]=h*(dydx[i]+h*dfdx[i]);
    lubksb(a,n,indx,yout);
    for (i=1;i<=n;i++)
        ytemp[i]=y[i]+(del[i]=yout[i]);
    x=xs+h;
    (*derivs)(x,ytemp,yout);
    for (nn=2;nn<=nstep;nn++) {
        for (i=1;i<=n;i++)
            yout[i]=h*yout[i]-del[i];
        lubksb(a,n,indx,yout);
        for (i=1;i<=n;i++)
            ytemp[i] += (del[i] += 2.0*yout[i]);
        x += h;
        (*derivs)(x,ytemp,yout);
    }
    for (i=1;i<=n;i++)
        yout[i]=h*yout[i]-del[i];
    lubksb(a,n,indx,yout);
    for (i=1;i<=n;i++)
        yout[i] += ytemp[i];
    free_dvector(ytemp,1,n);
    free_dvector(del,1,n);
    free_dmatrix(a,1,n,1,n);
    free_ivector(indx,1,n);
}
//#undef NRANSI

/*
void jacobn(double x, double y[], double dfdx[], double **dfdy, int n)
{
    int i;
    
    for (i=1;i<=n;i++) dfdx[i]=0.0;
    dfdy[1][1] = -0.013-1000.0*y[3];
    dfdy[1][2]=0.0;
    dfdy[1][3] = -1000.0*y[1];
    dfdy[2][1]=0.0;
    dfdy[2][2] = -2500.0*y[3];
    dfdy[2][3] = -2500.0*y[2];
    dfdy[3][1] = -0.013-1000.0*y[3];
    dfdy[3][2] = -2500.0*y[3];
    dfdy[3][3] = -1000.0*y[1]-2500.0*y[2];
}
*/

/*
void derivs(double x, double y[], double dydx[])
{
    dydx[1] = -0.013*y[1]-1000.0*y[1]*y[3];
    dydx[2] = -2500.0*y[2]*y[3];
    dydx[3] = -0.013*y[1]-1000.0*y[1]*y[3]-2500.0*y[2]*y[3];
}
*/


void lubksb(double **a, int n, int *indx, double b[])
{
    int i,ii=0,ip,j;
    double sum;
    
    for (i=1;i<=n;i++) {
        ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if (ii)
            for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
        else if (sum) ii=i;
        b[i]=sum;
    }
    for (i=n;i>=1;i--) {
        sum=b[i];
        for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
        b[i]=sum/a[i][i];
    }
}



//#include <math.h>
//#define NRANSI
//#include "nrutil.h"
#define TINY 1.0e-20

void ludcmp(double **a, int n, int *indx, double *d)
{
    int i,imax,j,k;
    double big,dum,sum,temp;
    double *vv;
    
    vv=dvector(1,n);
    *d=1.0;
    for (i=1;i<=n;i++) {
        big=0.0;
        for (j=1;j<=n;j++)
            if ((temp=fabs(a[i][j])) > big) big=temp;
        if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
        vv[i]=1.0/big;
    }
    for (j=1;j<=n;j++) {
        for (i=1;i<j;i++) {
            sum=a[i][j];
            for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        big=0.0;
        for (i=j;i<=n;i++) {
            sum=a[i][j];
            for (k=1;k<j;k++)
                sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
            if ( (dum=vv[i]*fabs(sum)) >= big) {
                big=dum;
                imax=i;
            }
        }
        if (j != imax) {
            for (k=1;k<=n;k++) {
                dum=a[imax][k];
                a[imax][k]=a[j][k];
                a[j][k]=dum;
            }
            *d = -(*d);
            vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (a[j][j] == 0.0) a[j][j]=TINY;
        if (j != n) {
            dum=1.0/(a[j][j]);
            for (i=j+1;i<=n;i++) a[i][j] *= dum;
        }
    }
    free_dvector(vv,1,n);
}
#undef TINY
//#undef NRANSI


//#include <math.h>
//#define NRANSI
//#include "nrutil.h"
#define KMAXX 7
#define IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define TINY 1.0e-30
#define SCALMX 0.1

double **d,*x;

void stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
            double yscal[], double *hdid, double *hnext,
            void (*derivs)(double, double [], double []),
            void (*jacobn)(double, double [], double [], double **, int))
{
//    void jacobn(double x, double y[], double dfdx[], double **dfdy, int n);
    void simpr(double y[], double dydx[], double dfdx[], double **dfdy,
               int n, double xs, double htot, int nstep, double yout[],
               void (*derivs)(double, double [], double []));
    void pzextr(int iest, double xest, double yest[], double yz[], double dy[],
                int nv);
    int i,iq,k,kk,km;
    static int first=1,kmax,kopt,nvold = -1;
    static double epsold = -1.0,xnew;
    double eps1,errmax,fact,h,red,scale,work,wrkmin,xest;
    double *dfdx,**dfdy,*err,*yerr,*ysav,*yseq;
    static double a[IMAXX+1];
    static double alf[KMAXX+1][KMAXX+1];
    static int nseq[IMAXX+1]={0,2,6,10,14,22,34,50,70};
    int reduct,exitflag=0;
    
    d=dmatrix(1,nv,1,KMAXX);
    dfdx=dvector(1,nv);
    dfdy=dmatrix(1,nv,1,nv);
    err=dvector(1,KMAXX);
    x=dvector(1,KMAXX);
    yerr=dvector(1,nv);
    ysav=dvector(1,nv);
    yseq=dvector(1,nv);
    if(eps != epsold || nv != nvold) {
        *hnext = xnew = -1.0e29;
        eps1=SAFE1*eps;
        a[1]=nseq[1]+1;
        for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
        for (iq=2;iq<=KMAXX;iq++) {
            for (k=1;k<iq;k++)
                alf[k][iq]=pow(eps1,((a[k+1]-a[iq+1])/
                                     ((a[iq+1]-a[1]+1.0)*(2*k+1))));
        }
        epsold=eps;
        nvold=nv;
        a[1] += nv;
        for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
        for (kopt=2;kopt<KMAXX;kopt++)
            if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
        kmax=kopt;
    }
    h=htry;
    for (i=1;i<=nv;i++) ysav[i]=y[i];
    jacobn(*xx,y,dfdx,dfdy,nv);
    if (*xx != xnew || h != (*hnext)) {
        first=1;
        kopt=kmax;
    }
    reduct=0;
    for (;;) {
        for (k=1;k<=kmax;k++) {
            xnew=(*xx)+h;
            if (xnew == (*xx)) nrerror("step size underflow in stifbs");
            simpr(ysav,dydx,dfdx,dfdy,nv,*xx,h,nseq[k],yseq,derivs);
            xest=SQR(h/nseq[k]);
            pzextr(k,xest,yseq,y,yerr,nv);
            if (k != 1) {
                errmax=TINY;
                for (i=1;i<=nv;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
                errmax /= eps;
                km=k-1;
                err[km]=pow(errmax/SAFE1,1.0/(2*km+1));
            }
            if (k != 1 && (k >= kopt-1 || first)) {
                if (errmax < 1.0) {
                    exitflag=1;
                    break;
                }
                if (k == kmax || k == kopt+1) {
                    red=SAFE2/err[km];
                    break;
                }
                else if (k == kopt && alf[kopt-1][kopt] < err[km]) {
                    red=1.0/err[km];
                    break;
                }
                else if (kopt == kmax && alf[km][kmax-1] < err[km]) {
                    red=alf[km][kmax-1]*SAFE2/err[km];
                    break;
                }
                else if (alf[km][kopt] < err[km]) {
                    red=alf[km][kopt-1]/err[km];
                    break;
                }
            }
        }
        if (exitflag) break;
        red=FMIN(red,REDMIN);
        red=FMAX(red,REDMAX);
        h *= red;
        reduct=1;
    }
    *xx=xnew;
    *hdid=h;
    first=0;
    wrkmin=1.0e35;
    for (kk=1;kk<=km;kk++) {
        fact=FMAX(err[kk],SCALMX);
        work=fact*a[kk+1];
        if (work < wrkmin) {
            scale=fact;
            wrkmin=work;
            kopt=kk+1;
        }
    }
    *hnext=h/scale;
    if (kopt >= k && kopt != kmax && !reduct) {
        fact=FMAX(scale/alf[kopt-1][kopt],SCALMX);
        if (a[kopt+1]*fact <= wrkmin) {
            *hnext=h/fact;
            kopt++;
        }
    }
    free_dvector(yseq,1,nv);
    free_dvector(ysav,1,nv);
    free_dvector(yerr,1,nv);
    free_dvector(x,1,KMAXX);
    free_dvector(err,1,KMAXX);
    free_dmatrix(dfdy,1,nv,1,nv);
    free_dvector(dfdx,1,nv);
    free_dmatrix(d,1,nv,1,KMAXX);
}
#undef KMAXX
#undef IMAXX
#undef SAFE1
#undef SAFE2
#undef REDMAX
#undef REDMIN
#undef TINY
#undef SCALMX
//#undef NRANSI



//#include <math.h>
//#define NRANSI
//#include "nrutil.h"
#define SAFETY 0.9
#define GROW 1.5
#define PGROW -0.25
#define SHRNK 0.5
#define PSHRNK (-1.0/3.0)
#define ERRCON 0.1296
#define MAXTRY 40
#define GAM (1.0/2.0)
#define A21 2.0
#define A31 (48.0/25.0)
#define A32 (6.0/25.0)
#define C21 -8.0
#define C31 (372.0/25.0)
#define C32 (12.0/5.0)
#define C41 (-112.0/125.0)
#define C42 (-54.0/125.0)
#define C43 (-2.0/5.0)
#define B1 (19.0/9.0)
#define B2 (1.0/2.0)
#define B3 (25.0/108.0)
#define B4 (125.0/108.0)
#define E1 (17.0/54.0)
#define E2 (7.0/36.0)
#define E3 0.0
#define E4 (125.0/108.0)
#define C1X (1.0/2.0)
#define C2X (-3.0/2.0)
#define C3X (121.0/50.0)
#define C4X (29.0/250.0)
#define A2X 1.0
#define A3X (3.0/5.0)

void stiff(double y[], double dydx[], int n, double *x, double htry, double eps,
           double yscal[], double *hdid, double *hnext,
           void (*derivs)(double, double [], double []),
           void (*jacobn)(double, double [], double [], double **, int))
{
//    void jacobn(double x, double y[], double dfdx[], double **dfdy, int n);
    void lubksb(double **a, int n, int *indx, double b[]);
    void ludcmp(double **a, int n, int *indx, double *d);
    int i,j,jtry,*indx;
    double d,errmax,h,xsav,**a,*dfdx,**dfdy,*dysav,*err;
    double *g1,*g2,*g3,*g4,*ysav;
    
    indx=ivector(1,n);
    a=dmatrix(1,n,1,n);
    dfdx=dvector(1,n);
    dfdy=dmatrix(1,n,1,n);
    dysav=dvector(1,n);
    err=dvector(1,n);
    g1=dvector(1,n);
    g2=dvector(1,n);
    g3=dvector(1,n);
    g4=dvector(1,n);
    ysav=dvector(1,n);
    xsav=(*x);
    for (i=1;i<=n;i++) {
        ysav[i]=y[i];
        dysav[i]=dydx[i];
    }
    jacobn(xsav,ysav,dfdx,dfdy,n);
    h=htry;
    for (jtry=1;jtry<=MAXTRY;jtry++) {
        for (i=1;i<=n;i++) {
            for (j=1;j<=n;j++) a[i][j] = -dfdy[i][j];
            a[i][i] += 1.0/(GAM*h);
        }
        ludcmp(a,n,indx,&d);
        for (i=1;i<=n;i++)
            g1[i]=dysav[i]+h*C1X*dfdx[i];
        lubksb(a,n,indx,g1);
        for (i=1;i<=n;i++)
            y[i]=ysav[i]+A21*g1[i];
        *x=xsav+A2X*h;
        (*derivs)(*x,y,dydx);
        for (i=1;i<=n;i++)
            g2[i]=dydx[i]+h*C2X*dfdx[i]+C21*g1[i]/h;
        lubksb(a,n,indx,g2);
        for (i=1;i<=n;i++)
            y[i]=ysav[i]+A31*g1[i]+A32*g2[i];
        *x=xsav+A3X*h;
        (*derivs)(*x,y,dydx);
        for (i=1;i<=n;i++)
            g3[i]=dydx[i]+h*C3X*dfdx[i]+(C31*g1[i]+C32*g2[i])/h;
        lubksb(a,n,indx,g3);
        for (i=1;i<=n;i++)
            g4[i]=dydx[i]+h*C4X*dfdx[i]+(C41*g1[i]+C42*g2[i]+C43*g3[i])/h;
        lubksb(a,n,indx,g4);
        for (i=1;i<=n;i++) {
            y[i]=ysav[i]+B1*g1[i]+B2*g2[i]+B3*g3[i]+B4*g4[i];
            err[i]=E1*g1[i]+E2*g2[i]+E3*g3[i]+E4*g4[i];
        }
        *x=xsav+h;
        if (*x == xsav) nrerror("stepsize not significant in stiff");
        errmax=0.0;
        for (i=1;i<=n;i++) errmax=FMAX(errmax,fabs(err[i]/yscal[i]));
        errmax /= eps;
        if (errmax <= 1.0) {
            *hdid=h;
            *hnext=(errmax > ERRCON ? SAFETY*h*pow(errmax,PGROW) : GROW*h);
            free_dvector(ysav,1,n);
            free_dvector(g4,1,n);
            free_dvector(g3,1,n);
            free_dvector(g2,1,n);
            free_dvector(g1,1,n);
            free_dvector(err,1,n);
            free_dvector(dysav,1,n);
            free_dmatrix(dfdy,1,n,1,n);
            free_dvector(dfdx,1,n);
            free_dmatrix(a,1,n,1,n);
            free_ivector(indx,1,n);
            return;
        } else {
            *hnext=SAFETY*h*pow(errmax,PSHRNK);
            h=(h >= 0.0 ? FMAX(*hnext,SHRNK*h) : FMIN(*hnext,SHRNK*h));
        }
    }
    nrerror("exceeded MAXTRY in stiff");
}
#undef SAFETY
#undef GROW
#undef PGROW
#undef SHRNK
#undef PSHRNK
#undef ERRCON
#undef MAXTRY
#undef GAM
#undef A21
#undef A31
#undef A32
#undef C21
#undef C31
#undef C32
#undef C41
#undef C42
#undef C43
#undef B1
#undef B2
#undef B3
#undef B4
#undef E1
#undef E2
#undef E3
#undef E4
#undef C1X
#undef C2X
#undef C3X
#undef C4X
#undef A2X
#undef A3X
//#undef NRANSI

