/*==============================================================================
 MODULE: mglpt_postprocess.c			[mgpt]
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
 http://www.inin.gob.mx/
 
 Major revisions:
 Copyright: (c) 2005-2018 Mario A. Rodriguez-Meza.  All Rights Reserved
 ================================================================================
 Legal matters:
 The author does not warrant that the program and routines it contains
 listed below are free from error or suitable for particular applications,
 and he disclaims all liability from any consequences arising from their	use.
 ==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"

#define EPSQ 1.0e-6

local real Q12_function(real eta, real ki);

local  real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[]);

local real sigma2L_function_int(real y);
local real sigma2L_function_ver2(void);
local real sigma2L_function(real ki);

local real PSLF(real k);
local real Q1F(real k);
local real Q2F(real k);
local real Q3F(real k);
local real Q5F(real k);
local real QIF(real k);
local real R1F(real k);
local real R2F(real k);
local real R1plus2F(real k);
local real RIF(real k);

local real tildeV(real k);
local real tildeT(real k);
local real xL(real k, real q);
local real xloop(real k, real q);
local real yL(real k, real q);
local real yloop(real k, real q);
local real x10(real k, real q);
local real y10(real k, real q);

local void InputQsRsTable(void);
local void InputQsRsTable_old(void);

// q functions
local real XLT_function(real qi);
global_qfunctions qfunctions(real qi);


local void postprocess_string_to_int(string, int *);
local void correlation_processing(void);
local void prd_processing(void);
//local void biasterms_processing(void);

// Qs and Rs table structure
typedef struct _pointQsRsTable {
    real k;             //1
    real Q1;            //2
    real Q2;            //3
    real Q3;            //4
    real Q5;            //5
    real Q7;            //6
    real Q8;            //7
    real Q9;            //8
    real Q11;           //9
    real Q12;           //10
    real Q13;           //11
    real QI;            //12
    real R1;            //13
    real R2;            //14
    real R1plus2;       //15
    real RI;            //16
    real Dpk;           //17
    real PSMGL;         //18
} pointQsRsTable, *pointQsRsTableptr;

local int nQsRsTable;
local pointQsRsTableptr QsRstab;

#define kQsRs(x)    (((pointQsRsTableptr) (x))->k)
#define Q1QsRs(x)    (((pointQsRsTableptr) (x))->Q1)
#define Q2QsRs(x)    (((pointQsRsTableptr) (x))->Q2)
#define Q3QsRs(x)    (((pointQsRsTableptr) (x))->Q3)
#define Q5QsRs(x)    (((pointQsRsTableptr) (x))->Q5)
#define Q7QsRs(x)    (((pointQsRsTableptr) (x))->Q7)
#define Q8QsRs(x)    (((pointQsRsTableptr) (x))->Q8)
#define Q9QsRs(x)    (((pointQsRsTableptr) (x))->Q9)
#define Q11QsRs(x)    (((pointQsRsTableptr) (x))->Q11)
#define Q12QsRs(x)    (((pointQsRsTableptr) (x))->Q12)
#define Q13QsRs(x)    (((pointQsRsTableptr) (x))->Q13)
#define QIQsRs(x)    (((pointQsRsTableptr) (x))->QI)
#define R1QsRs(x)    (((pointQsRsTableptr) (x))->R1)
#define R2QsRs(x)    (((pointQsRsTableptr) (x))->R2)
#define R1plus2QsRs(x)    (((pointQsRsTableptr) (x))->R1plus2)
#define RIQsRs(x)    (((pointQsRsTableptr) (x))->RI)
#define DpkQsRs(x)    (((pointQsRsTableptr) (x))->Dpk)
#define PSMGLQsRs(x)    (((pointQsRsTableptr) (x))->PSMGL)

local pointQsRsTableptr PQsRstab;


#define CORRELATION                 1
#define BIASTERMS                   3


global void PostProcessing(void)
{
    int process_int;

    fprintf(stdout,"\n\nStarting postprocessing...\n");
    
    postprocess_string_to_int(cmd.options, &process_int);
    switch (process_int){

//        case BIASTERMS: biasterms_processing(); break;
        case CORRELATION: correlation_processing(); break;

        default: error("\nUnknown postprocessing type %s\n\n",cmd.options);
    }
}


local void postprocess_string_to_int(string process_str,int *process_int)
{
    *process_int = -1;
    if (strcmp(process_str,"biasterms") == 0)               *process_int=BIASTERMS;
    if (strcmp(process_str,"correlation") == 0)             *process_int=CORRELATION;
}

#undef BIAS
#undef CORRELATION
#undef PRD
#undef BIASTERMS

//#define fpfnameQsRs "kfunctions.dat"
#define fpfnameak "ak.dat"
#define fpfnamekernel "kernel.dat"
//#define fpfnamebiasterms "biasterms.dat"

local real *kTab;
local real *Q1T;
local real *Q1T2;
local real *Q2T;
local real *Q2T2;
local real *Q3T;
local real *Q3T2;
local real *Q5T;
local real *Q5T2;
local real *QIT;
local real *QIT2;
local real *R1T;
local real *R1T2;
local real *R2T;
local real *R2T2;
local real *R1plus2T;
local real *R1plus2T2;
local real *RIT;
local real *RIT2;

local void correlation_processing(void)
{
    stream outstr;
    pointQsRsTableptr p;
    int i;
    real PSPTk, P11k, P22k, sigma2L=1.0;
    real k, q=1.0;

    fprintf(stdout,"\n\nCorrelation processing...\n");
    InputQsRsTable();
    
    kTab = dvector(1,nQsRsTable);
    Q1T = dvector(1,nQsRsTable);
    Q1T2 = dvector(1,nQsRsTable);
    Q2T = dvector(1,nQsRsTable);
    Q2T2 = dvector(1,nQsRsTable);
    Q5T = dvector(1,nQsRsTable);
    Q5T2 = dvector(1,nQsRsTable);
    QIT = dvector(1,nQsRsTable);
    QIT2 = dvector(1,nQsRsTable);
    R1T = dvector(1,nQsRsTable);
    R1T2 = dvector(1,nQsRsTable);
    R2T = dvector(1,nQsRsTable);
    R2T2 = dvector(1,nQsRsTable);
    R1plus2T = dvector(1,nQsRsTable);
    R1plus2T2 = dvector(1,nQsRsTable);
    RIT = dvector(1,nQsRsTable);
    RIT2 = dvector(1,nQsRsTable);

    i=1;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        kTab[i] = kQsRs(p);
        Q1T[i] = Q1QsRs(p);
        Q2T[i] = Q2QsRs(p);
        Q5T[i] = Q5QsRs(p);
        QIT[i] = QIQsRs(p);
        R1T[i] = R1QsRs(p);
        R2T[i] = R2QsRs(p);
        R1plus2T[i] = R1plus2QsRs(p);
        RIT[i] = RIQsRs(p);
        i++;
    }

    spline(kTab,Q1T,nQsRsTable,1.0e30,1.0e30,Q1T2);
    spline(kTab,Q2T,nQsRsTable,1.0e30,1.0e30,Q2T2);
    spline(kTab,Q5T,nQsRsTable,1.0e30,1.0e30,Q5T2);
    spline(kTab,QIT,nQsRsTable,1.0e30,1.0e30,QIT2);
    spline(kTab,R1T,nQsRsTable,1.0e30,1.0e30,R1T2);
    spline(kTab,R2T,nQsRsTable,1.0e30,1.0e30,R2T2);
    spline(kTab,R1plus2T,nQsRsTable,1.0e30,1.0e30,R1plus2T2);
    spline(kTab,RIT,nQsRsTable,1.0e30,1.0e30,RIT2);

// tildeV   = -(3./35.) (QI[k] - 3. Q2[k] + 2. RI[k] - 6. R2[k])
// tildeT   =-(9./14.)  (QI[k] + 2. Q2[k] + 2. RI[k] + 4. R2[k])
// xL       = PSL[k] (1./3. - j1[k q]/(k q))
// xloop    = (9./98. Q1[k] + 10./21. R1[k]) (1./3. - j1[k q]/(k q))
// yL       = PSL[k] j2[k q]
// yloop    = (9./98. Q1[k] + 10./21. R1[k]) j2[k q]
// x10      = 1./14. (2. RI[k] - 2. R2[k] + 3. RI[k] j0[k q]
//              - 3 (RI[k] + 2. R2[k] + 2. R1plus2[k] + 2. Q5[k]) j1[k q]/(k q))
//
// y10      = -3./14. (RI[k] + 2. R2[k] + 2. R1plus2[k] + 2. Q5[k])
//              * (j0[k q] - 3 j1[k q]/(k q))
//
    fprintf(stdout,"\n\nWriting kernel functions to file %s...\n",fpfnamekernel);
    outstr = stropen(fpfnamekernel,"w!");
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        
        P11k = (10.0/21.0)*R1QsRs(p) + (6.0/7.0)*R2QsRs(p)
        + sigma2L*rsqr(kQsRs(p))*Q3QsRs(p);
        P22k = (9.0/98.0)*Q1QsRs(p) + (3.0/7.0)*Q2QsRs(p) + (1.0/2.0)*Q3QsRs(p);
        PSPTk = PSMGLQsRs(p) + P22k + P11k;
        
        k = kQsRs(p);
        fprintf(outstr,"%g %g %g %g %g %g %g %g %g %g %g %g\n",
                kQsRs(p),
                tildeV(k),tildeT(k),
                xL(k,q),xloop(k,q),
                yL(k,q),yloop(k,q),
                x10(k,q),y10(k,q),
                Q12_function(0., k),sigma2L_function(k),
                sigma2L_function_ver2()
                );
    }
    fclose(outstr);

    free_dvector(RIT2,1,nQsRsTable);
    free_dvector(RIT,1,nQsRsTable);
    free_dvector(R1plus2T2,1,nQsRsTable);
    free_dvector(R1plus2T,1,nQsRsTable);
    free_dvector(R2T2,1,nQsRsTable);
    free_dvector(R2T,1,nQsRsTable);
    free_dvector(R1T2,1,nQsRsTable);
    free_dvector(R1T,1,nQsRsTable);
    free_dvector(QIT2,1,nQsRsTable);
    free_dvector(QIT,1,nQsRsTable);
    free_dvector(Q5T2,1,nQsRsTable);
    free_dvector(Q5T,1,nQsRsTable);
    free_dvector(Q2T2,1,nQsRsTable);
    free_dvector(Q2T,1,nQsRsTable);
    free_dvector(Q1T2,1,nQsRsTable);
    free_dvector(Q1T,1,nQsRsTable);
    free_dvector(kTab,1,nQsRsTable);
    free(PQsRstab);
}

global void biasterms_processing(void)
{
// P22 = 9./98. Q1T[[ii, 2]] + 3./7. Q2T[[ii, 2]] + 1./2. Q3T[[ii, 2]]
// P13 = 10./21. R1T[[ii, 2]] + 6./7. Q2T[[ii, 2]] - sigma2L kT[[ii]] kT[[ii]] PSLT[[ii, 2]]
// a10 = 10./21. R1T[[ii, 2]] + 6./7. R1plus2T[[ii, 2]] + 6./7. R2T[[ii, 2]]
//       + 6./7. Q5T[[ii, 2]] + 2 Q7T[[ii, 2]] + 2. (1. - sigma2L kT[[ii]] kT[[ii]]) PSLT[[ii, 2]]
// a01 = Q9T[[ii, 2]] + 3./7. Q8T[[ii, 2]]
// a20 = Q11T[[ii, 2]] + (1. - sigma2L kT[[ii]] kT[[ii]] ) PSLT[[ii, 2]]
//       + Q9T[[ii, 2]] + 6./7. R1plus2T[[ii, 2]]
// a11 = 2. Q12T[[ii, 2]]
// a02 = 1./2. Q13T[[ii, 2]]
//
    stream outstr;
    pointQsRsTableptr p;
    real sigma2L;
    real P22, P13;
    real a10, a01, a20, a11, a02;
    real k, PSLT;
    real Ploop;
    real a02Off;

    fprintf(stdout,"\n\nBias terms processing...\n");
    InputQsRsTable();
    
    sigma2L = sigma2L_function_ver2();
    fprintf(stdout,"\n\nsigma2L = %g\n",sigma2L);
    
// a02 offset
    a02Off = (1./2.)*Q13QsRs(PQsRstab + 2);

    fprintf(stdout,"\n\nWriting bias terms and the power spectrum to file %s...\n",
            gd.fpfnameSPTPowerSpectrum);
    outstr = stropen(gd.fpfnameSPTPowerSpectrum,"w!");
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {

        k = kQsRs(p);
        PSLT = PSMGLQsRs(p);
        P22 = (9./98.)*Q1QsRs(p) + (3./7.)*Q2QsRs(p) + (1./2.)*Q3QsRs(p);
        P13 = (10./21.)*R1QsRs(p) + (6./7.)*Q2QsRs(p) - sigma2L*k*k*PSLT;
        a10 = (10./21.)*R1QsRs(p) + (6./7.)*R1plus2QsRs(p) + (6./7.)*R2QsRs(p)
                + (6./7.)*Q5QsRs(p) + 2.*Q7QsRs(p) + 2.*(1. - sigma2L*k*k)*PSLT;
        a01 = Q9QsRs(p) + (3./7.)*Q8QsRs(p);
        a20 = Q11QsRs(p) + (1. - sigma2L*k*k)*PSLT
              + Q9QsRs(p) + (6./7.)*R1plus2QsRs(p);
        a11 = 2.*Q12QsRs(p);
        a02 = (1./2.)*Q13QsRs(p) - a02Off;
        Ploop = PSLT + P22 + P13;
        fprintf(outstr,"%g %g %g %g %g %g %g %g %g %g\n",
                k, PSLT,
                P22, rabs(P13), a10, rabs(a01), a20, a11, rabs(a02), Ploop
                );
    }
    fclose(outstr);
    //
    free(PQsRstab);
}

local void InputQsRsTable(void)
{
    stream outstr;
    pointQsRsTableptr p;
    int i;
//
    fprintf(gd.outlog,"\n\nReading solution Q1 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 2, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputQsRsTable: nQsRsTable = %d is absurd\n\n", nQsRsTable);

    PQsRstab = (pointQsRsTableptr) allocate(nQsRsTable * sizeof(pointQsRsTable));
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        kQsRs(p) = inout_xval[i];
        Q1QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Q2 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 3, &nQsRsTable);

    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q2QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Q3 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 4, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q3QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Q5 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 5, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q5QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Q7 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 6, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q7QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Q8 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 7, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q8QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Q9 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 8, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q9QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Q11 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 9, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q11QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Q12 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 10, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q12QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Q13 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 11, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        Q13QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution QI from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 12, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        QIQsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution R1 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 13, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        R1QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution R2 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 14, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        R2QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution R1plus2 from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 15, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        R1plus2QsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution RI from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 16, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        RIQsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution Dpk from file %s...",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 17, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        DpkQsRs(p) = inout_yval[i];
        ++i;
    }
//
    fprintf(gd.outlog,"\nReading solution PSMGL from file %s...\n",gd.fpfnamekfun);
    inout_InputData(gd.fpfnamekfun, 1, 18, &nQsRsTable);
    
    if (nQsRsTable < 1)
        error("\n\nInputSolsTable: nSolsTable = %d is absurd\n\n", nQsRsTable);
    
    i = 0;
    for (p=PQsRstab; p<PQsRstab+nQsRsTable; p++) {
        PSMGLQsRs(p) = inout_yval[i];
        ++i;
    }
//
}

//#undef fpfnameQsRs
#undef fpfnameak
#undef fpfnamekernel
//#undef fpfnamebiasterms


local real tildeV(real k)
{
// -(3./35.) (QI[k] - 3. Q2[k] + 2. RI[k] - 6. R2[k])
    real func;
    
    func= -(3./35.)*( QIF(k) - 3.*Q2F(k) + 2.* RIF(k) - 6.*R2F(k) );
    
    return (func);
}

local real tildeT(real k)
{
// -(9./14.)  (QI[k] + 2. Q2[k] + 2. RI[k] + 4. R2[k])
    real func;
    
    func= -(9./14.)*(QIF(k) + 2.*Q2F(k) + 2.*RIF(k) + 4.*R2F(k) );
    
    return (func);
}

local real xL(real k, real q)
{
// PSL[k] (1./3. - j1[k q]/(k q))
    real func;
    
    func= PSLF(k)*(1./3. - rj1Bessel(k*q)/(k*q) );

    return (func);
}

local real xloop(real k, real q)
{
// (9./98. Q1[k] + 10./21. R1[k]) (1./3. - j1[k q]/(k q))
    real func;
    
    func= ( (9./98.)*Q1F(k) + (10./21.)*R1F(k) )*(1./3. - rj1Bessel(k*q)/(k*q) );
    
    return (func);
}

local real yL(real k, real q)
{
// PSL[k] j2[k q]
    real func;
    
    func= PSLF(k)*rj2Bessel(k*q);
    
    return (func);
}

local real yloop(real k, real q)
{
// (9./98. Q1[k] + 10./21. R1[k]) j2[k q]
    real func;
    
    func= ( (9./98.)*Q1F(k) + (10./21.)*R1F(k) )*rj2Bessel(k*q);
    
    return (func);
}

local real x10(real k, real q)
{
// 1./14. (2. RI[k] - 2. R2[k] + 3. RI[k] j0[k q]
//    -3 (RI[k] + 2. R2[k] + 2. R1plus2[k] + 2. Q5[k]) j1[k q]/(k q))
    real func;
    
    func= (1./14.)*(
                    2.*RIF(k) - 2.*R2F(k) + 3.*RIF(k)*rj0Bessel(k*q)
                    -3.*( RIF(k) + 2.*R2F(k) + 2.*R1plus2F(k) + 2.*Q5F(k) )
                    *rj1Bessel(k*q)/(k*q)
                    );
    
    return (func);
}

local real y10(real k, real q)
{
// -3./14. (RI[k] + 2. R2[k] + 2. R1plus2[k] + 2. Q5[k]) (j0[k q]
//    -3 j1[k q]/(k q))
    real func;
    
    func= (-3./14.)*(
                     RIF(k) + 2.*R2F(k) + 2.*R1plus2F(k) + 2.*Q5F(k)
                     )
                    *(rj0Bessel(k*q) - 3.*rj1Bessel(k*q)/(k*q));
    
    return (func);
}

local real PSLF(real k)
{
    real func;
    func = psInterpolation_nr(k, kPS, pPS, nPSLT);
    return (func);
}

local real Q1F(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, Q1T, nQsRsTable, Q1T2);
    return (func);
}

local real Q2F(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, Q2T, nQsRsTable, Q2T2);
    return (func);
}

local real Q3F(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, Q3T, nQsRsTable, Q3T2);
    return (func);
}

local real Q5F(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, Q5T, nQsRsTable, Q5T2);
    return (func);
}

local real QIF(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, QIT, nQsRsTable, QIT2);
    return (func);
}

local real R1F(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, R1T, nQsRsTable, R1T2);
    return (func);
}

local real R2F(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, R2T, nQsRsTable, R2T2);
    return (func);
}

local real R1plus2F(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, R1plus2T, nQsRsTable, R1plus2T2);
    return (func);
}

local real RIF(real k)
{
    real func;
    func = Interpolation_nr(k, kTab, RIT, nQsRsTable, RIT2);
    return (func);
}

local  real Interpolation_nr(real k, double kPS[], double pPS[], int nPS, double pPS2[])
{
    real psftmp;

    if ( k < kPS[1] || k > kPS[nPS] )
        fprintf(gd.outlog,"\n\nInterpolation_nr: warning! :: k is out of range... %g\n",k);

    splint(kPS,pPS,pPS2,nPS,k,&psftmp);
    
    return (psftmp);
}

#define abskmq      (1.0+rsqr(rr)-2.0*rr*xv)

// Q12
local real Q12_function(real eta, real ki)
{
    int i, j;

    real kmin, kmax;
    real PSLA, PSLB;
    real kk, rr, deltar;
    real xv, w, k2, psl;

    real Q12p, Q12aA, Q12aB, KQ12;
    
    pointPSTableptr p;

    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    fprintf(gd.outlog,"nPSLT, kmin and kmax :: %d %g %g\n",nPSLT,kmin,kmax);

    Q12p = 0.0;
    Q12aA = 0.0;
    Q12aB = 0.0;
    
    PSLA = 0.0;
    p = PSLCDMtab;
//
    for (i=1; i<nPSTable; i++) {
        kk = kPos(p+i);
        rr = kk/ki;
        PSLB = psInterpolation_nr(ki*rr, kPS, pPS, nPSLT);
//
        for (j=1; j<=nGL(pGL); j++) {
            xv = xGL(pGL)[j];
            w = wGL(pGL)[j];
            psl = psInterpolation_nr(ki*rsqrt(abskmq), kPS, pPS, nPSLT);
            KQ12 = rr*xv;
            Q12aB += w*KQ12*psl;
        }
//
        deltar = (kPos(p+i)-kPos(p+i-1))/ki;
        Q12p += deltar*(Q12aA*PSLA + Q12aB*PSLB)/2.0;
        PSLA = PSLB;
        Q12aB = 0.0;
    }
    
    Q12p *= (rpow(ki,3)/FOURPI2);

    return Q12p;
}

// q functions
global_qfunctions qfunctions(real qi)
{
    global_qfunctions_ptr qfuntmp;
    pointPSTableptr p;
    int i;
    real kmin, kmax;
    real kk;
    real deltak;

    qfuntmp = (global_qfunctions_ptr) allocate(1 * sizeof(global_qfunctions));

    real XLp, XLA, XLB;
    real YLp, YLA, YLB;

    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    fprintf(gd.outlog,"nPSLT, kmin and kmax :: %d %g %g\n",nPSLT,kmin,kmax);
    
    p = PSLCDMtab;
//
    XLp = 0.;
    XLA = xL(kmin, qi);
    YLp = 0.;
    YLA = yL(kmin, qi);
//
    for (i=1; i<nPSTable; i++) {
        kk = kPos(p+i);
        deltak = (kPos(p+i) - kPos(p+i-1));
//
        XLB = xL(kk, qi);
        XLp = XLp + (XLA + XLB)*deltak/2.0;
        XLA = XLB;
//
        YLB = yL(kk, qi);
        YLp = YLp + (YLA + YLB)*deltak/2.0;
        YLA = YLB;

    }
    XLp /= PI2;
    YLp /= PI2;
    
    qqfun(qfuntmp) = qi;
    XLqfun(qfuntmp) = XLp;
    YLqfun(qfuntmp) = YLp;

    return *qfuntmp;
}


// XLT
local real XLT_function(real qi)
{
    int i;
    real kmin, kmax;
    real kk;
    real XLp, XLA, XLB;
    real deltak;
    pointPSTableptr p;
    
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    fprintf(gd.outlog,"nPSLT, kmin and kmax :: %d %g %g\n",nPSLT,kmin,kmax);
    
    p = PSLCDMtab;
    //
    XLp = 0.;
    XLA = xL(kmin, qi);
    //
    for (i=1; i<nPSTable; i++) {
        kk = kPos(p+i);
        XLB = xL(kk, qi);
        XLp = XLp + (XLA + XLB)*deltak/2.0;
        XLA = XLB;
    }
    XLp /= PI2;
    
    return XLp;
}


// sigma2L
local real sigma2L_function(real ki)
{
    int i;
    
    real kmin, kmax;
    real PSLA, PSLB;
    real deltar, kk;

    real sigma2L;
    
    pointPSTableptr p;
    
    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    fprintf(gd.outlog,"nPSLT, kmin and kmax :: %d %g %g\n",nPSLT,kmin,kmax);
    
    sigma2L = 0.0;
    
    PSLA = 0.0;
    p = PSLCDMtab;
//
    for (i=1; i<nPSTable; i++) {
        kk = kPos(p+i);
        PSLB = psInterpolation_nr(kk, kPS, pPS, nPSLT);
        deltar = (kPos(p+i)-kPos(p+i-1))/kk;
        sigma2L += deltar*(PSLA + PSLB)/2.0;
        PSLA = PSLB;
    }

    sigma2L *= (1.0/SIXPI2);
    
    return sigma2L;
}

local real sigma2L_function_int(real y)
{
    real p;
    real PSL;
    
    p = rpow(10.0,y);

    PSL = psInterpolation_nr(p, kPS, pPS, nPSLT);

    return p*PSL;
}

local real sigma2L_function_ver2(void)
{
    real result;
    real kmin, kmax;
    real ymin, ymax;

    kmin = kPos(PSLT+1);
    kmax = kPos(PSLT+nPSLT-1);
    ymin = rlog10(kmin);
    ymax = rlog10(kmax);

    result= (1.0/SIXPI2)*rlog(10.0)
    *qromo(sigma2L_function_int,ymin,ymax,midpnt,EPSQ);
    
    return result;

}

#undef EPSQ
