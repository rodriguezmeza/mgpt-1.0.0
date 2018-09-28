
#ifndef _diffeqs_h
#define _diffeqs_h

#include "../general/stdinc.h"

void rk4(double y[], double dydx[], int n, double x, double h, double yout[],
               void (*derivs)(double, double [], double []));


//void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
//            double hmin, int *nok, int *nbad, int maxnsteps,
//            void (*derivsin)(double, double [], double []),
//            void (*rkqsin)(double [], double [], int, double *, double, double, double [],
//                           double *, double *, void (*)(double, double [], double [])));
void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
            double hmin, int *nok, int *nbad, int maxnsteps,
            void (*derivsin)(double, double [], double []),
            void (*jacobnin)(double, double [], double [], double **, int),
            void (*rkqsin)(double [], double [], int, double *, double, double, double [],
                           double *, double *, void (*)(double, double [], double []),
                           void (*)(double, double [], double [], double **, int)));


//void bsstep(double y[], double dydx[], int nv, double *xx, double htry, double eps,
//            double yscal[], double *hdid, double *hnext,
//            void (*derivs)(double, double [], double []));
void bsstep(double y[], double dydx[], int nv, double *xx, double htry, double eps,
            double yscal[], double *hdid, double *hnext,
            void (*derivs)(double, double [], double []),
            void (*jacobnin)(double, double [], double [], double **, int));

//void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
//          double yscal[], double *hdid, double *hnext,
//          void (*derivs)(double, double [], double []));
void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
          double yscal[], double *hdid, double *hnext,
          void (*derivs)(double, double [], double []),
          void (*jacobnin)(double, double [], double [], double **, int));

//void stifbs(double y[], double dydx[], int nv, double *xx,
//            double htry, double eps, double yscal[], double *hdid, double *hnext,
//            void (*derivs)(double, double [], double []));
void stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
            double yscal[], double *hdid, double *hnext,
            void (*derivs)(double, double [], double []),
            void (*jacobn)(double, double [], double [], double **, int));

//void stiff(double y[], double dydx[], int n, double *x,
//           double htry, double eps, double yscal[], double *hdid, double *hnext,
//           void (*derivs)(double, double [], double []));
void stiff(double y[], double dydx[], int n, double *x, double htry, double eps,
           double yscal[], double *hdid, double *hnext,
           void (*derivs)(double, double [], double []),
           void (*jacobn)(double, double [], double [], double **, int));

void lubksb(double **a, int n, int *indx, double b[]);
void ludcmp(double **a, int n, int *indx, double *d);
void jacobn(double x, double y[], double dfdx[], double **dfdy, int n);
void pzextr(int iest, double xest, double yest[], double yz[], double dy[],
            int nv);
void simpr(double y[], double dydx[], double dfdx[], double **dfdy,
           int n, double xs, double htot, int nstep, double yout[],
           void (*derivs)(double, double [], double []));


double dxsav,*xp,**yp;
int kmax,kount;
int nrhs;


#endif // ! _diffeqs_h


