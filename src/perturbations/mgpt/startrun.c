/*==============================================================================
 MODULE: startrun.c				[mgpt]
 Written by: Mario A. Rodriguez-Meza
 Starting date: January 2018
 Purpose: routines to initialize the main code
 Language: C
 Use: 'StartRun();'
 Routines and functions:
 Modules, routines and external headers:
 Coments and notes:
 Info: Mario A. Rodriguez-Meza
 Depto. de Fisica, ININ
 Apdo. Postal 18-1027 Mexico D.F. 11801 Mexico
 e-mail: marioalberto.rodriguez@inin.gob.mx
 http://www.inin.gob.mx/
 
 Mayor revisions:
 Copyright: (c) 2005-2018 Mar.  All Rights Reserved
 ================================================================================
 Legal matters:
 The author does not warrant that the program and routines it contains
 listed below are free from error or suitable for particular applications,
 and he disclaims all liability from any consequences arising from their	use.
 ==============================================================================*/

#include "globaldefs.h"
#include "protodefs.h"

//local void ReadMGModelParameterFile(char *);
//local void PrintMGModelParameterFile(char *);

local void ReadParameterFile(char *);
local void PrintParameterFile(char *);

local void startrun_parameterfile(void);
local void startrun_cmdline(void);
local void ReadParametersCmdline(void);
local void startrun_Common(void);
local void startrun_ParamStat(void);
local void CheckParameters(void);

//local void setFilesDirs_log(void);
//local void setFilesDirs(void);

local void InputPSTable(void);
local void PSLTable(void);
local void GaussLegendrePoints(void);

void StartRun(string head0, string head1, string head2, string head3)
{
    real aTime;
    aTime = cputime();

    gd.headline0 = head0; gd.headline1 = head1;
    gd.headline2 = head2; gd.headline3 = head3;
    printf("\n%s\n%s: %s\n\t %s\n",
		gd.headline0, gd.headline1, gd.headline2, gd.headline3);

//	gd.stopflag = 0;

    cmd.paramfile = GetParam("paramfile");
    if (!strnull(cmd.paramfile))
		startrun_parameterfile();
	else
		startrun_cmdline();

    StartOutput();

    fprintf(stdout,"\nTime to StartRun %g\n\n",cputime()-aTime);
    fflush(stdout);
}

local void startrun_parameterfile(void)
{
	ReadParameterFile(cmd.paramfile);
    startrun_ParamStat();
	startrun_Common();
	PrintParameterFile(cmd.paramfile);
}

#define parameter_null	"parameters_null-mgpt"

local void startrun_cmdline(void)
{
	ReadParametersCmdline();
	startrun_Common();
	PrintParameterFile(parameter_null);
}

local void ReadParametersCmdline(void)
{
// Modified gravity model parameters:
    cmd.mgmodel = GetParam("mgModel");
    cmd.model_paramfile = GetParam("modelParamfile");
    cmd.suffixModel = GetParam("suffixModel");
//    cmd.nHS = GetiParam("nHS");
    cmd.fR0 = GetdParam("fR0");
    cmd.screening = GetdParam("screening");
//
// DGP:
    cmd.eps_DGP = GetdParam("epsDGP");
    cmd.rc_DGP = GetdParam("rcDGP");
//
// Power spectrum table:
    cmd.fnamePS = GetParam("fnamePS");
    cmd.kmin = GetdParam("kmin");
    cmd.kmax = GetdParam("kmax");
    cmd.Nk = GetiParam("Nk");
// Power spectrum table interpolation and extrapolation parameters:
    cmd.kminT = GetdParam("kminT");
    cmd.kmaxT = GetdParam("kmaxT");
    cmd.Nkext = GetiParam("Nkext");
    cmd.NkL = GetiParam("NkL");
    cmd.NkU = GetiParam("NkU");
    cmd.NPTL = GetiParam("NPTL");
    cmd.NPTR = GetiParam("NPTR");
//
// CLPT correlation functions table:
    cmd.rmin = GetdParam("rmin");
    cmd.rmax = GetdParam("rmax");
    cmd.Nr = GetiParam("Nr");
// Background cosmology:
    cmd.om = GetdParam("Om");
    cmd.olstr = GetParam("OL");
    cmd.h = GetdParam("h");
//
// Differential equations evolution parameters:
    cmd.x = GetdParam("etaini");
    cmd.dxstr = GetParam("deta");
    cmd.dxmin = GetdParam("detamin");
    cmd.eps = GetdParam("eps");
    cmd.xstop = GetdParam("zout");
    cmd.maxnsteps = GetiParam("maxnsteps");
	cmd.integration_method = GetParam("solverMethod");
//
// Quadrature parameters:
    cmd.quadratureMethod = GetParam("quadratureMethod");
    cmd.nquadSteps = GetiParam("nquadSteps");
    cmd.ngausslegpoints = GetiParam("ngausslegpoints");
    cmd.epsquad = GetdParam("epsquad");
// Post processing parameters:
    cmd.postprocessing = GetbParam("postprocessing");
    cmd.options = GetParam("options");
}

#undef parameter_null


local void startrun_Common(void)
{
	real dx1, dx2;
    char *ep;

    setFilesDirs_log();
    strcpy(gd.mode,"w");
	if(!(gd.outlog=fopen(gd.logfilePath, gd.mode)))
		error("\nstart_Common: error opening file '%s' \n",gd.logfilePath);

// Modified gravity model parameters:
//

// Background cosmology:
//
    ep = strchr(cmd.olstr, '-');
    if (ep == NULL) {
        gd.ol = GetdParam("OL");
        fprintf(gd.outlog,"\nOL string input without '-' :: %s\n",cmd.olstr);
        fprintf(gd.outlog,"\nOLambda, Om and sum  : %g %g %g\n",gd.ol, cmd.om, gd.ol+cmd.om);
    } else {
        ep = strchr(cmd.olstr, '1');
        if (ep == NULL)
            error("\nstart_Common: OL not in the format '1 - Om' \n");
        else {
            ep = strchr(cmd.olstr, 'O');
            if (ep == NULL)
                error("\nstart_Common: OL not in the format '1 - Om' \n");
            else
                if (*(ep+1) == 'm') {
                    fprintf(gd.outlog,"\nFound Om\n");
                    gd.ol = 1. - cmd.om;
                    fprintf(gd.outlog,"\nOLambda and Om : %g %g\n",gd.ol, cmd.om);
                } else {
                    fprintf(gd.outlog,"\nNot found Om\n");
                    error("\nstart_Common: OL not in the format (1 - Om) \n");
                }
        }
        
    }
//

// Differential equations evolution parameters:
    gd.dx = (sscanf(cmd.dxstr, "%lf/%lf", &dx1, &dx2) == 2 ?
				dx1/dx2 : atof(cmd.dxstr));
    if ( dx2 == 0. )
        error("\n\nstartrun_Common: dx : dx2 must be finite\n");

    CheckParameters();
    GaussLegendrePoints();
    quadraturemethod_string_to_int(cmd.quadratureMethod, &gd.quadmethod_int);
    integration_method_string_to_int(cmd.integration_method, &gd.method_int);
    gd.xnow = cmd.x;
    gd.xout = gd.xnow;
    gd.xoutinfo = gd.xnow;

// eta = Log[1/(1 + z)] -> exp[eta] = 1/(1 + z) -> (1 + z) = exp[-eta] -> z = exp[-eta] - 1
    gd.xstop = rlog(1.0/(1.0+cmd.xstop));

    set_model();
    setFilesDirs();

    if (!strnull(cmd.fnamePS)) {
        InputPSTable();
        PSLTable();
    }
    if (cmd.nquadSteps > nPSLT)
        error("CheckParameters: nquadSteps > nPSLT\n");

/*
    if (!strnull(cmd.model_paramfile)
        && !(strcmp(cmd.mgmodel,"HS") == 0)
        && !(strcmp(cmd.mgmodel,"DGP") == 0)
        && !(strcmp(cmd.mgmodel,"USERMODEL") == 0)
        && !(strcmp(cmd.mgmodel,"LCDM") == 0)
        ) {
        fprintf(stdout,"\n\nNot default model, using parameter file: %s\n",cmd.model_paramfile);
        ReadMGModelParameterFile(cmd.model_paramfile);
        PrintMGModelParameterFile(cmd.model_paramfile);
    }
*/
#define usrmodel_parameter_null    "usrmodel_parameters_null"

    if ( (strcmp(cmd.mgmodel,"USER") == 0) || (strcmp(cmd.mgmodel,"user") == 0))
        if (!strnull(cmd.model_paramfile)) {
            fprintf(stdout,"\n\nUser model :: using parameter file: %s\n",cmd.model_paramfile);
            ReadMGModelParameterFile(cmd.model_paramfile);
            PrintMGModelParameterFile(cmd.model_paramfile);
        } else
            PrintMGModelParameterFile(usrmodel_parameter_null);
#undef usrmodel_parameter_null
}

local void startrun_ParamStat(void)
{
	real dx1, dx2;

// Power spectrum table:
    if (GetParamStat("fnamePS") & ARGPARAM)
        cmd.fnamePS = GetParam("fnamePS");
    if (GetParamStat("kmin") & ARGPARAM)
        cmd.kmin = GetdParam("kmin");
    if (GetParamStat("kmax") & ARGPARAM)
        cmd.kmax = GetdParam("kmax");
    if (GetParamStat("Nk") & ARGPARAM)
        cmd.Nk = GetiParam("Nk");
// Power spectrum table interpolation and extrapolation parameters:
    if (GetParamStat("kminT") & ARGPARAM)
        cmd.kminT = GetdParam("kminT");
    if (GetParamStat("kmaxT") & ARGPARAM)
        cmd.kmaxT = GetdParam("kmaxT");
    if (GetParamStat("Nkext") & ARGPARAM)
        cmd.Nkext = GetiParam("Nkext");
    if (GetParamStat("NkL") & ARGPARAM)
        cmd.NkL = GetiParam("NkL");
    if (GetParamStat("NkU") & ARGPARAM)
        cmd.NkU = GetiParam("NkU");
    if (GetParamStat("NPTL") & ARGPARAM)
        cmd.NPTL = GetiParam("NPTL");
    if (GetParamStat("NPTR") & ARGPARAM)
        cmd.NPTR = GetiParam("NPTR");
//
// CLPT correlation functions table:
    if (GetParamStat("rmin") & ARGPARAM)
        cmd.rmin = GetdParam("rmin");
    if (GetParamStat("rmax") & ARGPARAM)
        cmd.rmax = GetdParam("rmax");
    if (GetParamStat("Nr") & ARGPARAM)
        cmd.Nr = GetiParam("Nr");

// Modified gravity model parameters:
    if (GetParamStat("mgModel") & ARGPARAM)
        cmd.mgmodel = GetParam("mgModel");
    if (GetParamStat("suffixModel") & ARGPARAM)
        cmd.suffixModel = GetParam("suffixModel");
//    if (GetParamStat("nHS") & ARGPARAM)
//        cmd.nHS = GetiParam("nHS");
    if (GetParamStat("fR0") & ARGPARAM)
        cmd.fR0 = GetdParam("fR0");
    if (GetParamStat("screening") & ARGPARAM)
        cmd.rc_DGP = GetdParam("screening");
// DGP:
    if (GetParamStat("epsDGP") & ARGPARAM)
        cmd.eps_DGP = GetdParam("epsDGP");
    if (GetParamStat("rcDGP") & ARGPARAM)
        cmd.rc_DGP = GetdParam("rcDGP");

// Background cosmology:
    if (GetParamStat("Om") & ARGPARAM)
        cmd.om = GetdParam("Om");
    if (GetParamStat("OL") & ARGPARAM)
        cmd.olstr = GetParam("OL");
    if (GetParamStat("h") & ARGPARAM)
        cmd.h = GetdParam("h");

// Differential equations evolution parameters:
	if (GetParamStat("etaini") & ARGPARAM)
		cmd.x = GetdParam("etaini");
	if (GetParamStat("deta") & ARGPARAM) {
		cmd.dxstr = GetParam("deta");
		gd.dx = (sscanf(cmd.dxstr, "%lf/%lf", &dx1, &dx2) == 2 ?
                    dx1/dx2 : atof(cmd.dxstr));
		if ( dx2 == 0. )
			error("\n\nstartrun_ParamStat: deta : deta2 must be finite\n");
	}
	if (GetParamStat("zout") & ARGPARAM)
		cmd.xstop = GetdParam("zout");

	if (GetParamStat("maxnsteps") & ARGPARAM) 
		cmd.maxnsteps = GetiParam("maxnsteps");

	if (GetParamStat("solverMethod") & ARGPARAM) {
		cmd.integration_method = GetParam("solverMethod");
		fprintf(gd.outlog,"\n\nrunning instead %s integration method ...\n",
				cmd.integration_method);
	}

// Quadrature parameters:
    if (GetParamStat("quadratureMethod") & ARGPARAM) {
        cmd.quadratureMethod = GetParam("quadratureMethod");
        fprintf(gd.outlog,"\n\nrunning instead %s quadrature method ...\n",
                cmd.quadratureMethod);
    }
    if (GetParamStat("nquadSteps") & ARGPARAM)
        cmd.nquadSteps = GetiParam("nquadSteps");
    if (GetParamStat("ngausslegpoints") & ARGPARAM)
        cmd.ngausslegpoints = GetiParam("ngausslegpoints");
    if (GetParamStat("epsquad") & ARGPARAM)
        cmd.epsquad = GetdParam("epsquad");

// Post processing parameters:
    if (GetParamStat("postprocessing") & ARGPARAM)
        cmd.postprocessing = GetbParam("postprocessing");
    
    if (GetParamStat("options") & ARGPARAM)
        cmd.options = GetParam("options");
    
}

local void CheckParameters(void)
{
// Modified gravity model parameters:

//
// Power spectrum table:
    if (strnull(cmd.fnamePS))
        error("CheckParameters: You should give a power spectrum filename\n");
    if (cmd.kmin < 0.0)
        error("CheckParameters: absurd value for kmin\n");
    if (cmd.kmax < 0.0)
        error("CheckParameters: absurd value for kmax\n");
    if (cmd.kmin > cmd.kmax)
        error("CheckParameters: kmin can not be greater than kmax\n");
    if (cmd.Nk < 0)
        error("CheckParameters: absurd value for Nk\n");
// Power spectrum table interpolation and extrapolation parameters:
    if (cmd.kminT < 0.0)
        error("CheckParameters: absurd value for kminT\n");
    if (cmd.kmaxT < 0.0)
        error("CheckParameters: absurd value for kmaxT\n");
    if (cmd.kminT > cmd.kmaxT)
        error("CheckParameters: kminT can not be greater than kmaxT\n");
    if (cmd.Nkext < 0)
        error("CheckParameters: absurd value for Nkext\n");
    if (cmd.NkL < 0)
        error("CheckParameters: absurd value for NkL\n");
    if (cmd.NkU < 0)
        error("CheckParameters: absurd value for NkU\n");
    if (cmd.NPTL < 0)
        error("CheckParameters: absurd value for NPTL\n");
    if (cmd.NPTR < 0)
        error("CheckParameters: absurd value for NPTR\n");
//
// CLPT correlation functions table:
    if (cmd.rmin < 0.0)
        error("CheckParameters: absurd value for rmin\n");
    if (cmd.rmax < 0.0)
        error("CheckParameters: absurd value for rmax\n");
    if (cmd.rmin > cmd.rmax)
        error("CheckParameters: rmin can not be greater than rmax\n");
    if (cmd.Nr < 0)
        error("CheckParameters: absurd value for Nr\n");
// Background cosmology:
    if (cmd.om > 1.0 || cmd.om < 0.0)
        error("CheckParameters: absurd value for om\n");
    if ( gd.ol < 0. )
        error("\n\nstartrun_ParamStat: OL (=%g) : must be positive\n",gd.ol);
    if (cmd.h < 0.0)
        error("CheckParameters: absurd value for h\n");

// Differential equations evolution parameters:
    if (gd.dx == 0)
        error("CheckParameters: absurd value for deta\n");
    if(cmd.x == rlog(1.0/(1.0+cmd.xstop)) )
        error("\n\nstartrun_Common: etaini and etastop=exp(-zout)-1 must be different\n");
    if (cmd.maxnsteps < 1)
        error("CheckParameters: absurd value for maxnsteps\n");

// Quadrature parameters:
    if (cmd.nquadSteps <= 1)
        error("CheckParameters: absurd value for nquadSteps\n");
    if (cmd.ngausslegpoints <= 1)
        error("CheckParameters: absurd value for ngausslegpoints\n");
    if (cmd.epsquad >= 1.0e-1 || cmd.epsquad <= 0)
        error("CheckParameters: absurd value for epsquad\n");
}

/*
// I/O directories:
local void setFilesDirs_log(void)
{
    char buf[200];

    sprintf(gd.tmpDir,"tmp");

    sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi",gd.tmpDir,gd.tmpDir);
    system(buf);

    sprintf(gd.logfilePath,"%s/mgpt%s.log",gd.tmpDir,cmd.suffixModel);
}

local void setFilesDirs(void)
{
    char buf[200];

    sprintf(gd.clptDir,"CLPT");
    sprintf(buf,"if [ ! -d %s ]; then mkdir %s; fi",gd.clptDir,gd.clptDir);
    fprintf(gd.outlog,"system: %s\n",buf);
    system(buf);

    sprintf(gd.inputDir,"Input");
    sprintf(gd.fnamePS,"%s/%s",gd.inputDir,cmd.fnamePS);

    sprintf(gd.fpfnamekfun,"CLPT/kfunctions%s.dat",cmd.suffixModel);
    sprintf(gd.fpfnameSPTPowerSpectrum,"SPTPowerSpectrum%s.dat",cmd.suffixModel);
    sprintf(gd.fpfnameqfunctions,"CLPT/qfunctions%s.dat",cmd.suffixModel);
    sprintf(gd.fpfnameclptfunctions,"CorrelationFunction%s.dat",cmd.suffixModel);
}
*/

local void ReadParameterFile(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define BOOLEAN 4
#define MAXTAGS 300

  FILE *fd,*fdout;

  char buf[200],buf1[200],buf2[200],buf3[200];
  int  i,j,nt;
  int  id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int  errorFlag=0;

  nt=0;

// Power spectrum table:
    SPName(cmd.fnamePS,"fnamePS",100);
    RPName(cmd.kmin,"kmin");
    RPName(cmd.kmax,"kmax");
    IPName(cmd.Nk,"Nk");
// Power spectrum table interpolation and extrapolation parameters:
    RPName(cmd.kminT,"kminT");
    RPName(cmd.kmaxT,"kmaxT");
    IPName(cmd.Nkext,"Nkext");
    IPName(cmd.NkL,"NkL");
    IPName(cmd.NkU,"NkU");
    IPName(cmd.NPTL,"NPTL");
    IPName(cmd.NPTR,"NPTR");
//
// CLPT correlation functions table:
    RPName(cmd.rmin,"rmin");
    RPName(cmd.rmax,"rmax");
    IPName(cmd.Nr,"Nr");
// Modified gravity model parameters:
    SPName(cmd.mgmodel,"mgModel",100);
    SPName(cmd.model_paramfile,"modelParamfile",100);
    SPName(cmd.suffixModel,"suffixModel",100);
//    IPName(cmd.nHS,"nHS");
    RPName(cmd.fR0,"fR0");
    RPName(cmd.screening,"screening");
//
// DGP:
    RPName(cmd.eps_DGP,"epsDGP");
    RPName(cmd.rc_DGP,"rcDGP");
//
// Background cosmology:
    RPName(cmd.om,"Om");
    SPName(cmd.olstr,"OL",100);
    RPName(cmd.h,"h");
//
// Differential equations evolution parameters:
    RPName(cmd.x,"etaini");
	SPName(cmd.dxstr,"deta",100);
    RPName(cmd.dxmin,"detamin");
    RPName(cmd.eps,"eps");
	RPName(cmd.xstop,"zout");
    IPName(cmd.maxnsteps,"maxnsteps");
	SPName(cmd.integration_method,"solverMethod",100);
//
// Quadrature parameters:
    SPName(cmd.quadratureMethod,"quadratureMethod",100);
    IPName(cmd.nquadSteps,"nquadSteps");
    IPName(cmd.ngausslegpoints,"ngausslegpoints");
    RPName(cmd.epsquad,"epsquad");
//
// Post processing parameters:
    BPName(cmd.postprocessing,"postprocessing");
    SPName(cmd.options,"options",100);
//

	if((fd=fopen(fname,"r"))) {
		while(!feof(fd)) {
			fgets(buf,200,fd);
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<1)
                continue;
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
                *buf2='\0';
            if(buf1[0]=='%')
				continue;
            for(i=0,j=-1;i<nt;i++)
                if(strcmp(buf1,tag[i])==0) {
                    j=i;
					tag[i][0]=0;
					break;
				}
			if(j>=0) {
                switch(id[j]) {
					case DOUBLE:
						*((double*)addr[j])=atof(buf2); 
						break;
					case STRING:
						strcpy(addr[j],buf2);
						break;
					case INT:
						*((int*)addr[j])=atoi(buf2);
						break;
					case BOOLEAN:
						if (strchr("tTyY1", *buf2) != NULL) {          
							*((bool*)addr[j])=TRUE;
                        } else 
                            if (strchr("fFnN0", *buf2) != NULL)  {
                                *((bool*)addr[j])=FALSE;
                            } else {
                                error("getbparam: %s=%s not bool\n",buf1,buf2);
                            }
						break;
                }
            } else {
                fprintf(stdout, "Error in file %s: Tag '%s' %s.\n",
					fname, buf1, "not allowed or multiple defined");
                errorFlag=1;
            }
        }
        fclose(fd);
    } else {
        fprintf(stdout,"Parameter file %s not found.\n", fname);
        errorFlag=1;
        exit(1); 
    }
  
    for(i=0;i<nt;i++) {
        if(*tag[i]) {
            fprintf(stdout,
                "Error. I miss a value for tag '%s' in parameter file '%s'.\n",
                tag[i],fname);
            exit(0);
        }
    }
#undef DOUBLE 
#undef STRING 
#undef INT 
#undef BOOLEAN
#undef MAXTAGS
}

#define FMTT	"%-35s%s\n"
#define FMTI	"%-35s%d\n"
#define FMTR	"%-35s%g\n"

local void PrintParameterFile(char *fname)
{
    FILE *fdout;
    char buf[200];

    sprintf(buf,"%s/%s%s%s",gd.tmpDir,fname,cmd.suffixModel,"-usedvalues");
    if(!(fdout=fopen(buf,"w"))) {
        fprintf(stdout,"error opening file '%s' \n",buf);
        exit(0);
    } else {
        fprintf(fdout,"%s\n",
		"%-------------------------------------------------------------------");
        fprintf(fdout,"%s %s\n","% Parameter input file for:",gd.headline0);
        fprintf(fdout,"%s\n","%");
        fprintf(fdout,"%s %s: %s\n%s\t    %s\n","%",gd.headline1,gd.headline2,"%",
            gd.headline3);
        fprintf(fdout,"%s\n%s\n",
		"%-------------------------------------------------------------------",
		"%");
// Power spectrum table:
        fprintf(fdout,FMTT,"fnamePS",cmd.fnamePS);
        fprintf(fdout,FMTR,"kmin",cmd.kmin);
        fprintf(fdout,FMTR,"kmax",cmd.kmax);
        fprintf(fdout,FMTI,"Nk",cmd.Nk);
// Power spectrum table interpolation and extrapolation parameters:
        fprintf(fdout,FMTR,"kminT",cmd.kminT);
        fprintf(fdout,FMTR,"kmaxT",cmd.kmaxT);
        fprintf(fdout,FMTI,"Nkext",cmd.Nkext);
        fprintf(fdout,FMTI,"NkL",cmd.NkL);
        fprintf(fdout,FMTI,"NkU",cmd.NkU);
        fprintf(fdout,FMTI,"NPTL",cmd.NPTL);
        fprintf(fdout,FMTI,"NPTR",cmd.NPTR);
//
// CLPT correlation functions table:
        fprintf(fdout,FMTR,"rmin",cmd.rmin);
        fprintf(fdout,FMTR,"rmax",cmd.rmax);
        fprintf(fdout,FMTI,"Nr",cmd.Nr);
// Modified gravity model parameters:
        fprintf(fdout,FMTT,"mgModel",cmd.mgmodel);
        fprintf(fdout,FMTT,"modelParamfile",cmd.model_paramfile);
        fprintf(fdout,FMTT,"suffixModel",cmd.suffixModel);
//        fprintf(fdout,FMTI,"nHS",cmd.nHS);
        fprintf(fdout,FMTR,"fR0",cmd.fR0);
        fprintf(fdout,FMTR,"screening",cmd.screening);
// DGP:
        fprintf(fdout,FMTR,"epsDGP",cmd.eps_DGP);
        fprintf(fdout,FMTR,"rcDGP",cmd.rc_DGP);
//
// Background cosmology:
        fprintf(fdout,FMTR,"Om",cmd.om);
        fprintf(fdout,FMTT,"OL",cmd.olstr);
        fprintf(fdout,FMTR,"h",cmd.h);
//
// Differential equations evolution parameters:
        fprintf(fdout,FMTR,"etaini",cmd.x);
        fprintf(fdout,FMTR,"detamin",cmd.dxmin);
        fprintf(fdout,FMTR,"eps",cmd.eps);
        fprintf(fdout,FMTT,"deta",cmd.dxstr);
        fprintf(fdout,FMTR,"zout",cmd.xstop);
        fprintf(fdout,FMTI,"maxnsteps",cmd.maxnsteps);
        fprintf(fdout,FMTT,"solverMethod",cmd.integration_method);
//
// Quadrature parameters:
        fprintf(fdout,FMTT,"quadratureMethod",cmd.quadratureMethod);
        fprintf(fdout,FMTI,"nquadSteps",cmd.nquadSteps);
        fprintf(fdout,FMTI,"ngausslegpoints",cmd.ngausslegpoints);
        fprintf(fdout,FMTR,"epsquad",cmd.epsquad);
// Post processing parameters:
        fprintf(fdout,FMTT,"postprocessing",cmd.postprocessing ? "true" : "false");
        fprintf(fdout,FMTT,"options",cmd.options);
        fprintf(fdout,"\n\n");
    }
    fclose(fdout);
}

#undef FMTT
#undef FMTI
#undef FMTR


/*
//=============================================================
// Begin: Modified gravity model reading and writing parameters

local void ReadMGModelParameterFile(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define BOOLEAN 4
#define MAXTAGS 300
    
    FILE *fd,*fdout;
    
    char buf[200],buf1[200],buf2[200],buf3[200];
    int  i,j,nt;
    int  id[MAXTAGS];
    void *addr[MAXTAGS];
    char tag[MAXTAGS][50];
    int  errorFlag=0;
    
    nt=0;

// Modified gravity model parameters:
    IPName(cmd.nHS,"nHS");
    RPName(cmd.fR0,"fR0");
    RPName(cmd.screening,"screening");
//
    if((fd=fopen(fname,"r"))) {
        while(!feof(fd)) {
            fgets(buf,200,fd);
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<1)
                continue;
            if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
                *buf2='\0';
            if(buf1[0]=='%')
                continue;
            for(i=0,j=-1;i<nt;i++)
                if(strcmp(buf1,tag[i])==0) {
                    j=i;
                    tag[i][0]=0;
                    break;
                }
            if(j>=0) {
                switch(id[j]) {
                    case DOUBLE:
                        *((double*)addr[j])=atof(buf2);
                        break;
                    case STRING:
                        strcpy(addr[j],buf2);
                        break;
                    case INT:
                        *((int*)addr[j])=atoi(buf2);
                        break;
                    case BOOLEAN:
                        if (strchr("tTyY1", *buf2) != NULL) {
                            *((bool*)addr[j])=TRUE;
                        } else
                            if (strchr("fFnN0", *buf2) != NULL)  {
                                *((bool*)addr[j])=FALSE;
                            } else {
                                error("getbparam: %s=%s not bool\n",buf1,buf2);
                            }
                        break;
                }
            } else {
                fprintf(stdout, "Error in file %s: Tag '%s' %s.\n",
                        fname, buf1, "not allowed or multiple defined");
                errorFlag=1;
            }
        }
        fclose(fd);
    } else {
        fprintf(stdout,"Parameter file %s not found.\n", fname);
        errorFlag=1;
        exit(1);
    }
    
    for(i=0;i<nt;i++) {
        if(*tag[i]) {
            fprintf(stdout,
                    "Error. I miss a value for tag '%s' in parameter file '%s'.\n",
                    tag[i],fname);
            exit(0);
        }
    }
#undef DOUBLE
#undef STRING
#undef INT
#undef BOOLEAN
#undef MAXTAGS
}

#define FMTT	"%-35s%s\n"
#define FMTI	"%-35s%d\n"
#define FMTR	"%-35s%g\n"

local void PrintMGModelParameterFile(char *fname)
{
    FILE *fdout;
    char buf[200];
    
    sprintf(buf,"%s%s%s",fname,cmd.suffixModel,"-usedvalues");
    if(!(fdout=fopen(buf,"w"))) {
        fprintf(stdout,"error opening file '%s' \n",buf);
        exit(0);
    } else {
        fprintf(fdout,"%s\n",
                "%-------------------------------------------------------------------");
        fprintf(fdout,"%s %s\n","% Modified gravity parameter model input file for:",gd.headline0);
        fprintf(fdout,"%s\n","%");
        fprintf(fdout,"%s %s: %s\n%s\t    %s\n","%",gd.headline1,gd.headline2,"%",
                gd.headline3);
        fprintf(fdout,"%s\n%s\n",
                "%-------------------------------------------------------------------",
                "%");
//
// Modified gravity model parameters:
        fprintf(fdout,FMTI,"nHS",cmd.nHS);
        fprintf(fdout,FMTR,"fR0",cmd.fR0);
        fprintf(fdout,FMTR,"screening",cmd.screening);
//
        fprintf(fdout,"\n\n");
    }
    fclose(fdout);
}

#undef FMTT
#undef FMTI
#undef FMTR

// End: Modified gravity model reading and writing parameters
//=============================================================
*/


//#define NPT 50
#define SPREAD 1.0
local void InputPSTable(void)
{
    stream outstr;
    pointPSTableptr p, plog, pn;
    pointPSTableptr PSLCDMtabtmp;
    int nPSTabletmp;
    int i;
    real dk, kval, PSval, kmin, kmax;
    int mwt;
    double al,bl,chi2,q,siga,sigb;
    double *xL,*yL,*sigL;
    double *xR,*yR,*sigR;
    double au, bu;
    real *kPStmp;
    real *pPStmp;
    real *pPS2tmp;
    char namebuf[256];
    real kminext, kmaxext, dktmp, kmn, kmx;
//    int Nkext=600, NkL=50, NkU=50;
//    real kminT=1.0e-6, kmaxT=350.0;


    fprintf(gd.outlog,"\n\nReading power spectrum from file %s...\n",gd.fnamePS);
    inout_InputData(gd.fnamePS, 1, 2, &nPSTabletmp);

    if (nPSTabletmp < 1)
        error("\n\nInputPSTable: nPSTable = %d is absurd\n\n", nPSTabletmp);

    PSLCDMtabtmp = (pointPSTableptr) allocate(nPSTabletmp * sizeof(pointPSTable));
    
    fprintf(gd.outlog,"nPSTable : %d\n", nPSTabletmp);

    i = 0;
    for (p=PSLCDMtabtmp; p<PSLCDMtabtmp+nPSTabletmp; p++) {
        kPos(p) = inout_xval[i];
        PS(p) = inout_yval[i];
        ++i;
    }

    fprintf(gd.outlog,"\n\nCreating log power spectrum...\n");

    PSLCDMLogtab = (pointPSTableptr) allocate(nPSTabletmp * sizeof(pointPSTable));
    plog = PSLCDMLogtab;
    nPSLogT=0;
    for (p=PSLCDMtabtmp; p<PSLCDMtabtmp+nPSTabletmp; p++) {
        kPos(plog) = rlog10(kPos(p));
        PS(plog) = rlog10(PS(p));
        plog++;
        nPSLogT++;
    }

    fprintf(gd.outlog,"\nTotal numbers in Log PS: %d %d\n",nPSLogT,plog-PSLCDMLogtab);
    fprintf(gd.outlog,"Total numbers in Normal PS: %d %d\n\n",nPSTabletmp,p-PSLCDMtabtmp);

    fprintf(gd.outlog,"\n\nLinear fit (a + b x) to log-log power spectrum at minset and maxset...\n");
    fprintf(gd.outlog,"\n\nLinear fit to log-log power spectrum at minset and maxset...\n");

    xL=dvector(1,cmd.NPTL);
    yL=dvector(1,cmd.NPTL);
    sigL=dvector(1,cmd.NPTL);
    xR=dvector(1,cmd.NPTR);
    yR=dvector(1,cmd.NPTR);
    sigR=dvector(1,cmd.NPTR);

// Lower part of the PS
    fprintf(gd.outlog,"\nAt lower part of the spectrum...\n");

    plog = PSLCDMLogtab;
    for (i=1;i<=cmd.NPTL;i++) {
        xL[i]=kPos(plog);
        yL[i]=PS(plog);
        sigL[i]=SPREAD;
        plog++;
    }
    for (mwt=0;mwt<=1;mwt++) {
        fit(xL,yL,cmd.NPTL,sigL,mwt,&al,&bl,&siga,&sigb,&chi2,&q);
        if (mwt == 0)
            fprintf(gd.outlog,"\nIgnoring standard deviations\n");
        else
            fprintf(gd.outlog,"\nIncluding standard deviations\n");
        fprintf(gd.outlog,"%12s %9.6f %18s %9.6f \n",
               "a  =  ",al,"uncertainty:",siga);
        fprintf(gd.outlog,"%12s %9.6f %18s %9.6f \n",
               "b  =  ",bl,"uncertainty:",sigb);
        fprintf(gd.outlog,"%19s %14.6f \n","chi-squared: ",chi2);
        fprintf(gd.outlog,"%23s %10.6f \n","goodness-of-fit: ",q);
    }

// Upper part of the PS
    fprintf(gd.outlog,"\nAt upper part of the spectrum...\n");

    plog = PSLCDMLogtab+nPSLogT-1;
    for (i=1;i<=cmd.NPTR;i++) {
        xR[i]=kPos(plog);
        yR[i]=PS(plog);
        sigR[i]=SPREAD;
        plog--;
    }
    for (mwt=0;mwt<=1;mwt++) {
        fit(xR,yR,cmd.NPTR,sigR,mwt,&au,&bu,&siga,&sigb,&chi2,&q);
        if (mwt == 0)
            fprintf(gd.outlog,"\nIgnoring standard deviations\n");
        else
            fprintf(gd.outlog,"\nIncluding standard deviations\n");
        fprintf(gd.outlog,"%12s %9.6f %18s %9.6f \n",
               "a  =  ",au,"uncertainty:",siga);
        fprintf(gd.outlog,"%12s %9.6f %18s %9.6f \n",
               "b  =  ",bu,"uncertainty:",sigb);
        fprintf(gd.outlog,"%19s %14.6f \n","chi-squared: ",chi2);
        fprintf(gd.outlog,"%23s %10.6f \n","goodness-of-fit: ",q);
    }

// Extending power spectrum
    kPStmp = dvector(1,nPSLogT);
    pPStmp = dvector(1,nPSLogT);
    pPS2 = dvector(1,nPSLogT);
    pn = PSLCDMLogtab;
    i=1;
    for (pn = PSLCDMLogtab; pn<PSLCDMLogtab+nPSLogT; pn++) {
        kPStmp[i] = kPos(pn);
        pPStmp[i] = PS(pn);
        i++;
    }
//
    spline(kPStmp,pPStmp,nPSLogT,1.0e30,1.0e30,pPS2);

    kmin = kPos(PSLCDMtabtmp);
    kmax = kPos(PSLCDMtabtmp+nPSTabletmp-1);
    fprintf(gd.outlog,"\nkmin, kmax of the given power spectrum (with %d values): %g %g",
            kmin, kmax, nPSTabletmp);
    dktmp = (rlog10(kmax) - rlog10(kmin))/((real)(nPSTabletmp - 1));
    kminext = rpow(10.0, rlog10(kmin)-((real)cmd.NkL)*dktmp);
    kmaxext = rpow(10.0, rlog10(kmax)+((real)cmd.NkU)*dktmp);
//
    fprintf(gd.outlog,"\n\nNkL, NkU: %d %d",cmd.NkL, cmd.NkU);
    cmd.NkL = rlog10(kmin/cmd.kminT)/dktmp;
    cmd.NkU = rlog10(cmd.kmaxT/kmax)/dktmp;
    fprintf(gd.outlog,"\nNkL, NkU targets: %d %d\n",cmd.NkL,cmd.NkU);
    kminext = rpow(10.0, rlog10(kmin)-((real)cmd.NkL)*dktmp);
    kmaxext = rpow(10.0, rlog10(kmax)+((real)cmd.NkU)*dktmp);
//
    
    
    fprintf(gd.outlog,"\nkmin, kmax of the extended power spectrum (first try): %g %g",
            kminext, kmaxext);

    kmn = MIN(kminext,cmd.kmin);
    kmx = MAX(kmaxext,cmd.kmax);
    fprintf(gd.outlog,"\nkmin, kmax of the extended power spectrum (second try): %g %g\n",
            kmn, kmx);

    nPSTable = cmd.Nkext;
    fprintf(gd.outlog,"\n\nCreating new PSTable with %d values\n",nPSTable);
    PSLCDMtab = (pointPSTableptr) allocate(nPSTable * sizeof(pointPSTable));
    dk = (rlog10(kmx) - rlog10(kmn))/((real)(nPSTable - 1));
    p = PSLCDMtab;
    
    for (i=1; i<=nPSTable; i++) {
        kval = rlog10(kmn) + dk*((real)(i - 1));
        if (rpow(10.0,kval) > kmin && rpow(10.0,kval) < kmax)
            PSval = psInterpolation_nr(kval, kPStmp, pPStmp, nPSLogT);
        else
            if (rpow(10.0,kval) < kmin)
                PSval = al + bl*kval;
            else
                if (rpow(10.0,kval) > kmax)
                    PSval = au + bu*kval;
                else
                    error("\n\nError: InputPSTable :: kmin, kmax, kval: %g %g %g\n\n", kmin, kmax, rpow(10.0,kval));

        kPos(p) = rpow(10.0,kval);
        PS(p) = rpow(10.0,PSval);
        p++;
    }


    sprintf(namebuf,"%s/%s%s_%s",gd.tmpDir,cmd.fnamePS,cmd.suffixModel,"ext.dat");
    outstr = stropen(namebuf,"w!");
    for (p=PSLCDMtab; p<PSLCDMtab+nPSTable; p++) {
        fprintf(outstr,"%g %g\n",
                kPos(p),PS(p));
    }
    fclose(outstr);

    gd.kminPSext = kmn;
    gd.kmaxPSext = kmx;
    fprintf(gd.outlog,"\nkmin, kmax of the extended power spectrum (final values): %g %g\n",
            gd.kminPSext, gd.kmaxPSext);
    fprintf(gd.outlog,"\nkmin, kmax of the extended power spectrum (PSExt table values): %g %g\n",
            kPos(PSLCDMtab), kPos(PSLCDMtab+nPSTable-1));
//

    free_dvector(pPS2,1,nPSLogT);
    free_dvector(pPStmp,1,nPSLogT);
    free_dvector(kPStmp,1,nPSLogT);

    free_dvector(sigR,1,cmd.NPTR);
    free_dvector(yR,1,cmd.NPTR);
    free_dvector(xR,1,cmd.NPTR);
    free_dvector(sigL,1,cmd.NPTL);
    free_dvector(yL,1,cmd.NPTL);
    free_dvector(xL,1,cmd.NPTL);
}
//#undef NPT
#undef SPREAD


local void PSLTable(void)
{
    char namebuf[256];
    stream outstr;
    real kmin, Dpkmin, Dpk;
    pointPSTableptr p, pn;
    int i;
//
    real xstoptmp, Dp0, Dpzout, fac;

    xstoptmp = gd.xstop;
    gd.xstop = 0.;
//    Dp0 = DpFunction(0.); // LCDM
    Dp0 = DpFunction_LCDM(0.); // LCDM
    fprintf(gd.outlog,"\n\n Dp(0) = %g\n",Dp0);
    gd.xstop = xstoptmp;
    Dpzout = DpFunction(0.);
    fprintf(gd.outlog,"\n\n Dp(%g) = %g\n",cmd.xstop,Dpzout);
    
    fac = rsqr(Dpzout/Dp0);
    for (p = PSLCDMtab; p<PSLCDMtab+nPSTable; p++) {
        PS(p) *= fac;
    }

    sprintf(namebuf,"%s/%s%s_%s",gd.tmpDir,cmd.fnamePS,cmd.suffixModel,"ext2.dat");
    outstr = stropen(namebuf,"w!");
    for (p=PSLCDMtab; p<PSLCDMtab+nPSTable; p++) {
        fprintf(outstr,"%g %g\n",
                kPos(p),PS(p));
    }
    fclose(outstr);

//

    kmin = kPos(PSLCDMtab);
    Dpkmin = DpFunction(kmin);
    fprintf(gd.outlog,"\n\n Dpkmin = %g\n",Dpkmin);

    PSLT = (pointPSTableptr) allocate(nPSTable * sizeof(pointPSTable));
    nPSLT = 0;
    pn = PSLT;
    for (p = PSLCDMtab; p<PSLCDMtab+nPSTable; p++) {
        kPos(pn) = kPos(p);
        Dpk = DpFunction(kPos(p));
        PS(pn) = rsqr(Dpk/Dpkmin)*PS(p);
        pn++;
        nPSLT++;
    }

    kPS = dvector(1,nPSLT);
    pPS = dvector(1,nPSLT);
    pPS2 = dvector(1,nPSLT);
    pn = PSLT;
    i=1;
    for (pn = PSLT; pn<PSLT+nPSLT; pn++) {
        kPS[i] = kPos(pn);
        pPS[i] = PS(pn);
        i++;
    }
//
    spline(kPS,pPS,nPSLT,1.0e30,1.0e30,pPS2);

    sprintf(namebuf,"%s/%s%s%s",gd.clptDir,"PSL",cmd.suffixModel,".dat");
    outstr = stropen(namebuf,"w!");
    for (i=1; i<=nPSLT; i++) {
        fprintf(outstr,"%g %g\n",
                kPS[i],pPS[i]);
    }
    fclose(outstr);

//    fprintf(stdout,"\nPSLCDMTable: %d, %g, %g",nPSTable, kPos(PSLCDMtab), kPos(PSLCDMtab+nPSTable-1));
//    fprintf(stdout,"\nPSLTable: %d, %g, %g",nPSLT, kPos(PSLT), kPos(PSLT+nPSLT-1));

}

local void GaussLegendrePoints(void)
{
    double dx,ss,xval;
    real x1=-1.0, x2=1.0;
    
    int i;
    double xx=0.0;
    
    pGL = (global_GL_ptr) allocate(sizeof(global_GL));

    nGL(pGL) = cmd.ngausslegpoints;
    xGL(pGL)=dvector(1,cmd.ngausslegpoints);
    wGL(pGL)=dvector(1,cmd.ngausslegpoints);
    x1GL(pGL) = -1.0;
    x2GL(pGL) = 1.0;

    fprintf(gd.outlog,"\nComputation of GLs...");
    gauleg(x1GL(pGL),x2GL(pGL),xGL(pGL),wGL(pGL),nGL(pGL));

    fprintf(gd.outlog," end of Gauss-Legendre computing (StartRun)\n\n");

}

//
// Power spectrum interpolation routines:
//
global  real psLCDMf(real k)
{
    pointPSTableptr p, pf, pi;
    int jl, ju, jm;
    real dk, psftmp;
    bool ascnd;
    
    pi = PSLCDMtab;
    pf = PSLCDMtab+nPSTable-1;
    
    if ( k < kPos(pi) || k > kPos(pf) || nPSTable < 2 )
        error("\n\npsLCDMf: k is out of range or nPSTable is wrong... %g %g %g\n",
              k,kPos(pi),kPos(pf));
    
    ascnd = (kPos(pf) >= kPos(pi));
    
    jl=0;
    ju=nPSTable-1;
    while (ju-jl > 1) {
        jm = (ju+jl) >> 1;
        if (k >= kPos(pi+jm) == ascnd)
            jl=jm;
        else
            ju=jm;
    }
    
    p = PSLCDMtab + jl;
    dk = kPos(p+1)-kPos(p);
    psftmp = PS(p)+(PS(p+1)-PS(p))*(k-kPos(p))/dk;
    return (psftmp);
}

global  real psInterpolation(real k, pointPSTableptr PSLtab, int nPSL)
{
    pointPSTableptr p, pf, pi;
    int jl, ju, jm;
    real dk, psftmp;
    bool ascnd;
    
    pi = PSLtab;
    pf = PSLtab+nPSL-1;
    
    if ( nPSL < 2 )
        error("\n\npsInterpolation: nPSL is wrong... %g\n",nPSL);
    
    if ( k > kPos(pf) ) {
        dk = kPos(pf) - kPos(pf-1);
        if ( dk > kPos(pf) - k ) {
            fprintf(stdout,"\n\npsInterpolation: warning!... extrapolating...\n",k);
            psftmp = PS(pf)+(PS(pf)-PS(pf-1))*(k-kPos(pf))/dk;
            return (psftmp);
        } else
            error("\n\npsInterpolation: k is out of range... %g\n",k);
    }
    
    if ( k < kPos(pi) )
        error("\n\npsInterpolation: k is out of range or nPSL is wrong... %g\n",k);
    
    ascnd = (kPos(pf) >= kPos(pi));
    
    jl=0;
    ju=nPSL-1;
    while (ju-jl > 1) {
        jm = (ju+jl) >> 1;
        if (k >= kPos(pi+jm) == ascnd)
            jl=jm;
        else
            ju=jm;
    }
    
    p = PSLtab + jl;
    dk = kPos(p+1)-kPos(p);
    psftmp = PS(p)+(PS(p+1)-PS(p))*(k-kPos(p))/dk;
    fprintf(gd.outlog,"%g %g %d PSL points...\n",k, psftmp, nPSL);
    fflush(gd.outlog);
    
    return (psftmp);
}

global  real psInterpolation_nr(real k, double kPS[], double pPS[], int nPS)
{
    pointPSTableptr pf, pi;
    real psftmp;
    real dps;
    real kf, ki;
    
    pi = PSLCDMtab;
    pf = PSLCDMtab+nPSTable-1;
    ki = kPS[1];
    kf = kPS[nPS];
    
    //
    //    if ( pow(10.0,k) < ki )
    //        fprintf(gd.outlog,"\n\npsInterpolation_nr: warning! :: k is out of range (%g, %g)... %g\n",
    //                ki, kf, pow(10.0,k));
    //    if ( pow(10.0,k) > kf )
    //        fprintf(gd.outlog,"\n\npsInterpolation_nr: warning! :: k is out of range (%g, %g)... %g\n",
    //                ki, kf, pow(10.0,k) );
    
    splint(kPS,pPS,pPS2,nPS,k,&psftmp);
    
    return (psftmp);
}



