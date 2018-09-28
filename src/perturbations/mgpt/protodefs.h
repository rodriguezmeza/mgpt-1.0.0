/*==============================================================================
 HEADER: protodefs.h				[mglpt]
 Written by: Mario A. Rodriguez-Meza
 Starting date: January 2018
 Purpose: Definitions of global prototypes
 Language: C
 Use: '#include "protodefs.h"
 Use in routines and functions:
 External headers: None
 Comments and notes:
 Info: Mario A. Rodriguez-Meza
 Depto. de Fisica, ININ
 Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
 e-mail: marioalberto.rodriguez@inin.gob.mx
 http://www.astro.inin.mx/mar
 
 Major revisions:
 Copyright: (c) 2005-2018 Mar.  All Rights Reserved
 ================================================================================
 Legal matters:
 The author does not warrant that the program and routines it contains
 listed below are free from error or suitable for particular applications,
 and he disclaims all liability from any consequences arising from their	use.
 ==============================================================================*/

#ifndef _protodefs_h
#define _protodefs_h

void integration_method_string_to_int(string,int *);
void quadraturemethod_string_to_int(string,int *);

//void output(void);

void MainLoop(void);
void StartRun(string, string, string, string);
void StartOutput(void);
void EndRun(void);

// Postprocessing
global void PostProcessing(void);
global void biasterms_processing(void);
global void qfunctions_processing(void);
global void CLPT_correlation_processing(void);

// MGLPT DIFFEQS
global real DpFunction(real k);
global real DpFunction_LCDM(real k);
global global_D2v2_ptr DsSecondOrder_func(real kf, real k1, real k2);
global global_D3v2_ptr DsThirdOrder_func(real x, real k, real p);

// MGLPT QUADS
global_QRs QsRs_functions_driver(real eta, real ki);
global_QRs QsRs_functions_driver_LCDM(real eta, real ki);

// I/O directories:
global void setFilesDirs_log(void);
global void setFilesDirs(void);

// Interpolation of the power spectrum
global real psLCDMf(real k);
global  real psInterpolation(real k, pointPSTableptr PSLtab, int nPSL);
global  real psInterpolation_nr(real k, double kPS[], double pPS[], int nPS);

#endif // ! _protodefs_h
