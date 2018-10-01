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
#define HEAD2	"mgpt code for modified gravity perturbation theory."
#define HEAD3	"..."

string defv[] = {  ";"HEAD1": " HEAD2 "\n\t " HEAD3,
    "paramfile=",                   ";Parameter input file. Overwritten by what follows",
//
// Power spectrum table:
    "fnamePS=psLCDM.in",		    ";Filename with power spectrum table",
    "kmin=1e-3",                    ";kmin to analyse from the power spectrum table",
    "kmax=100",                     ";kmax to analyse from the power spectrum table",
    "Nk=250",                       ";Total number of k´s analyse from the power spectrum table",":nk",
// Power spectrum table interpolation and extrapolation parameters:
    "kminT=1e-5",               ";kmin of extended (power spectrum table)",
    "kmaxT=400",                 ";kmax of extended (power spectrum table)",
    "Nkext=800",        ";Total number of extended k´s (power spectrum table)",
    "NkL=50",          ";NLower limit of extended k´s (power spectrum table)",
    "NkU=50",           ";NUper limit of extended k´s (power spectrum table)",
    "NPTL=5",        ";Total number of k to linear fit left side (power spectrum table)",
    "NPTR=12",        ";Total number of k to linear fit right side (power spectrum table)",
//
// CLPT correlation functions table:
    "rmin=50",                      ";rmin of the range for CLPT correlation functions table",
    "rmax=130",                     ";rmax of the range for CLPT correlation functions table",
    "Nr=100",                       ";Total number of r´s to analyse for the CLPT correlation functions table",":nr",
// Modified gravity model parameters:
    "mgModel=LCDM",                 ";Modified gravity model to study, default f(R) Hu-Sawicki", ":mgm",
    "suffixModel=",                 ";Suffix model to add to output filenames", ":suffix",
    "fR0=1.0e-5",                   ";Hu-Sawicky f_R0",
    "screening=1.0",                ";Hu-Sawicky screening", ":sc",
// DGP:
    "epsDGP=-1.0",                  ";DGP parameter, use with mgModel=DGP",":epsdgp",
    "rcDGP=1.0",                    ";DGP parameter, use with mgModel=DGP",":rcdgp",
//
    "modelParamfile=",	            ";If mgmodel=USER, you may give its parameter file name", ":mpf",
//
// Background cosmology:
    "Om=0.281",                     ";Omega matter value (z=0)",":om",
    "OL= 1 - Om",                   ";Omega Lambda value (z=0)",":ol",
    "h=0.697",                      ";Hubble parameter value (z=0)",
//
// Differential equations evolution parameters:
    "etaini=-4.0",                  ";Initial conformal time value :: Log[1/(1 + zini)]",
    "deta=2/5",                     ";Conformal time integration step",
    "detamin=0.",                   ";Min conformal time integration step size",
    "epsSolver=1.0e-4",             ";Differential equations solver tolerance error parameter",":epssolver",
    "zout=0.0",                     ";Output redshift value",
    "maxnsteps=10000",              ";Maximum number of integration steps", ":maxn",
    "solverMethod=bsstep",	        ";Integration method to solve differential equations", ":solver",
//
// Quadrature parameters:
    "quadratureMethod=trapezoid",   ";Quadrature method to use", ":quadm",
    "nquadSteps=200",               ";Number of k´s from the power spectrum table to integrate (trapezoid)",":nquad",
    "ngausslegpoints=10",           ";Number of Gauss-Legendre of integration points", ":nglpts",
    "epsquad=1.0e-5",               ";Quadrature tolerance error parameter (Romberg method: romberg)",
//
// Post processing parameters:
    "postprocessing=false",			";Post processing options", ":pp",
    "options=",                     ";Various control options", ":opt",
//
    "Version=1.0.0",                ";Mario A. Rodríguez-Meza/Alejandro Aviles 2018",
    NULL,
};

#endif // ! _cmdline_defs_h
