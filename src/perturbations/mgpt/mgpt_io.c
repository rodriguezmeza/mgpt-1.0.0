/*==============================================================================
 MODULE: mgpt_io.c		[mgpt]
 Written by: Mario A. Rodriguez-Meza
 Starting date:	January 2018
 Purpose: Routines to drive input and output data
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


void StartOutput(void)
{

    fprintf(stdout,"\n  \t -- %s --\n", gd.model_comment);
//
    fprintf(stdout,"  \t -- %s --\n", gd.integration_method_comment);
    fprintf(stdout,"  \t -- %s --\n\n", gd.quadraturemethod_comment);

    fprintf(gd.outlog,"\n%8s%8s%8s", "maxnsteps", "eta0", "deta");
    fprintf(gd.outlog,"%8s\n","etaout");
    fprintf(gd.outlog,"%8d%8.2f%8.4f%8.4f\n",cmd.maxnsteps,cmd.x,gd.dx,gd.xstop);
    if (! strnull(cmd.options))
        fprintf(stdout,"\n\toptions: %s\n", cmd.options);

}

// I/O directories:
global void setFilesDirs_log(void)
{
    char buf[200];
    
    sprintf(gd.tmpDir,"tmp");
    
    sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi",gd.tmpDir,gd.tmpDir);
    system(buf);
    
    sprintf(gd.logfilePath,"%s/mgpt%s.log",gd.tmpDir,cmd.suffixModel);
}

global void setFilesDirs(void)
{
    char buf[200];
    
    sprintf(gd.clptDir,"Output");
    sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi",gd.clptDir,gd.clptDir);
    fprintf(gd.outlog,"system: %s\n",buf);
    system(buf);
    
    sprintf(gd.inputDir,"Input");
    sprintf(gd.fnamePS,"%s/%s",gd.inputDir,cmd.fnamePS);
    
    sprintf(gd.fpfnamekfun,"Output/kfunctions%s.dat",cmd.suffixModel);
    sprintf(gd.fpfnameSPTPowerSpectrum,"SPTPowerSpectrum%s.dat",cmd.suffixModel);
    sprintf(gd.fpfnameqfunctions,"Output/qfunctions%s.dat",cmd.suffixModel);
    sprintf(gd.fpfnameclptfunctions,"CorrelationFunction%s.dat",cmd.suffixModel);
}

void EndRun(void)
{
    char   buf[200];
	FILE *fd;

	fclose(gd.outlog);
    
    printf("\nFinal CPU time : %g\n\n", second() - gd.cpuinit);
}



