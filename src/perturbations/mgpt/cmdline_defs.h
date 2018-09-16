/*==============================================================================
 HEADER: cmdline_defs.h		[mgpt]
 Written by: Mario A. Rodriguez-Meza
 Starting date: January 2018
 Purpose: Definitions for importing arguments from the command line
 Language: C
 Use: '#include "cmdline_defs.h"
 Use in routines and functions: (main)
 External headers: stdinc.h
 Comments and notes:
 Info: Mario A. Rodriguez-Meza
 Depto. de Fisica, ININ
 Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
 e-mail: marioalberto.rodriguez@inin.gob.mx
 http://www.inin.gob.mx/
 
 Major revisions:
 Copyright: (c) 2005-2018 Mar.  All Rights Reserved
 ================================================================================
 Legal matters:
 The author does not warrant that the program and routines it contains 
 listed below are free from error or suitable for particular applications, 
 and he disclaims all liability from any consequences arising from their use.
 ==============================================================================*/

#ifndef _cmdline_defs_h
#define _cmdline_defs_h

#define HEAD1	"NagBody"
#define HEAD2	"mgpt Code for the Modified Gravity-Perturbation Theory."
#define HEAD3	"..."

string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",                   ";Parameter input file. Overwritten by what follows",
//
// Power spectrum table:
    "fnamePS=psLCDM.in",		    ";Filename with power spectrum table",
    "kmin=1e-4",                    ";kmin to analyse from the power spectrum table",
    "kmax=100",                      ";kmax to analyse from the power spectrum table",
    "Nk=300",                       ";Total number of k´s analyse from the power spectrum table",":nk",
//
// Modified gravity model parameters:
    "mgModel=HS",                   ";Modified gravity model to study, default f(R) Hu-Sawicki", ":mgm",
    "suffixModel=",               ";Suffix model to add to output filenames", ":suffix",
    "nHS=1",                        ";Hu-Sawicki index",
    "fR0=1.0e-5",                   ";Hu-Sawicki f_R0",
    "screening=1.0",                ";Hu-Sawicki screening", ":sc",
// DGP:
    "eps_DGP=-1.0",                  ";DGP parameter, use with mgModel=DGP",":epsdgp",
    "rc_DGP=1.0",                       ";DGP parameter, use with mgModel=DGP",":rcdgp",
//
    "model_paramfile=fofRHS.in",	";If mgmodel is not the default, give its parameter file name", ":mpf",
//
// Background cosmology:
    "Om=0.281",                     ";Omega matter value (z=0)",":om",
    "OL= 1 - Om",                     ";Omega Lambda value (z=0)",":ol",
    "h=0.697",                      ";Hubble parameter value (z=0)",
//
// Differential equations evolution parameters:
    "etaini=-4.0",                    ";Initial eta value :: Log[1/(1 + zini)]",
    "deta=2/5",                     ";deta integration step",
    "detamin=0.",                   ";Min eta integration step size",
    "eps=1.0e-4",                   ";Error parameter",
    "zout=0.0",                    ";redshift value to stop integration :: eta = exp[-zout] - 1",
    "maxnsteps=10000",              ";Maximum number of integration steps", ":maxn",
    "integration_method=bsstep",	";Integration method to use", ":im",
//
// Quadrature parameters:
    "quadratureMethod=trapezoid",      ";Quadrature method to use", ":quadm",
    "nquadSteps=100",               ";Number of k´s from the power spectrum table to integrate (trapezoid)",":nquad",
    "ngausslegpoints=10",           ";Number of Gauss-Legendre of integration points", ":nglpts",
    "epsquad=1.0e-5",               ";Quadrature tolerance error parameter (open Romberg method)",
//
// Post processing parameters:
    "postprocessing=false",			";Post processing options", ":pp",
    "options=",                     ";Various control options", ":opt",
//
    "Version=1.0.0",                ";Mario A. Rodríguez-Meza 2005-2018",
    NULL,
};

#endif // ! _cmdline_defs_h
