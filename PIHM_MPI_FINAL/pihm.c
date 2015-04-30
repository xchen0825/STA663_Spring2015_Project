/*******************************************************************************
 * File        : pihm.c                                                        *
 * Version     : Nov, 2007 (2.0)                                               *
 * Developer of PIHM2.0:        Mukesh Kumar (muk139@psu.edu)                  *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0............................*
 * a) All modifications in physical process representations  in this version   *
 *    are listed as header in f.c and is_sm_et.c.     			       *
 * b) All addition/modifications in variable and structure definition/declarat-*
 *    -ion are listed as header in read_alloc.c and initialize.c	       *
 * c) 3 new input files have been added for geology, landcover and calibration *
 *    data								       *
 * d) Ported to Sundials 2.1.0                                                 *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * PIHM is an integrated finite volume based hydrologic model. It simulates    * 
 * channel routing, overland flow, groundwater flow, macropore based infiltra- *
 * tion and stormflow, throughfall, evaporation from overlandflow-subsurface-  *
 * canopy, transpiration and  snowmelt by full coupling of processes.          * 
 * It uses semi-discrete finite volume approach to discretize PDEs into ODEs,  * 
 * and henceforth solving the global system of ODEs using CVODE. Global ODEs   *
 * are created in f.c. Any modifications in the process equations has to be    *
 * performed in f.c
 *                                                                             *
 *-----------------------------------------------------------------------------*
 * For questions or comments, please contact                                   *
 *      --> Mukesh Kumar (mk176@duke.edu)                                      *
 * This code is free for research purpose only.                                *
 * Please provide relevant references if you use this code in your research work*
 *-----------------------------------------------------------------------------*
 *									       *
 * DEVELOPMENT RELATED REFERENCES:					       *
 * PIHM2.0:								       *
 *	a) Kumar, M., 2008, "Development and Implementation of a Multiscale,   *
 *	Multiprocess Hydrologic Model". PhD Thesis, Penn State University      *
 *	b) Kumar, M, G.Bhatt & C.Duffy, "Coupling of Data and Processes in     *
 *	Mesoscale Watershed", Advances in Water Resources (submitted)          *
 * PIHM1.0:								       *
 *	a) Qu, Y., 2005, "An Integrated hydrologic model for multiproces       *
 *	simulation using semi-discrete finite volume approach".PhD Thesis, PSU *
 *	b) Qu, Y. & C. Duffy, 2007, "A semidiscrete finite volume formulation  *
 *	for multiprocess watershed simulation". Water Resources Research       *
 *-----------------------------------------------------------------------------*
 * LICENSE: 
 *******************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* SUNDIAL Header Files */
#include "sundialstypes.h"   /* realtype, integertype, booleantype defination */
#include "cvode.h"           /* CVODE header file                             */
#include "cvspgmr.h"         /* CVSPGMR linear header file                    */
#include "smalldense.h"      /* use generic DENSE linear solver for "small"   */
//#include "nvector_serial.h"  /* contains the definition of type N_Vector      */		//xchen_20150328
#include "sundialsmath.h"    /* contains UnitRoundoff, RSqrt, SQR functions   */
#include "cvdense.h"         /* CVDENSE header file                           */
#include "dense.h"           /* generic dense solver header file              */
#include "nvector_parallel.h"  /* definition of type N_Vector and macro       */                //xchen_20150328
                               /* NV_DATA_P                             */
#include "mpi.h"             /* MPI constants and types                     */                  //xchen_20150328
#include "pihm.h"            /* Data Model and Variable Declarations     */

/* Function Declarations */
void initialize(char *, Model_Data, Control_Data *, N_Vector);
void is_sm_et(realtype, realtype, Model_Data);	
/* Function to calculate right hand side of ODE systems */
void f(realtype, N_Vector, N_Vector, void *);
void read_alloc(char *, Model_Data, Control_Data *);	/* Variable definition */
void update(realtype, Model_Data);	 
//void PrintData(FILE **,Control_Data *, Model_Data, N_Vector, realtype);

/* Main Function */
int main(int argc, char *argv[])
	{  
  	Model_Data mData;               /* Model Data                */
  	Control_Data cData;             /* Solver Control Data       */
 	N_Vector CV_Y,atol;             /* State Variables Vector    */
	void *cvode_mem;                /* Model Data Pointer        */
  	int flag;                       /* flag to test return value */
  	FILE *Ofile[25];           	/* Output file     */
	char *ofn[25];
	FILE *res_q1_file;		/* Output streamflow at observation site*/	//xchen_20150328
	realtype locQ1;									//xchen_20150328
	char *QFile;									//xchen_20150411
	FILE *iproj;			/* Project File */
  	int npes;			/* number of processors	     */ 	//xchen_20150328
	long int N,local_N;             /* Problem size              */
  	int i,j,k;                      /* loop index                */
  	realtype t,expl_t;    			/* simulation time           */
  	realtype NextPtr, StepSize;     /* stress period & step size */
  	clock_t start, end_r, end_s;    /* system clock at points    */
  	realtype cputime_r, cputime_s;  /* for duration in realtype  */
	char *filename;
        char tmpLname[25][12]={".NetP",".IS",".surf",".unsat",".GW",".et0",".et1",".et2",".P",".Rech",".unsatT",".RechT",".stage",".subriv"};
	
        CV_Y = NULL;
        cvode_mem = NULL;
        mData = NULL;
       
	MPI_Status status;							//xchen_20150328_start
	/* Get processor number, total number of pe's, and my_pe. */
	MPI_Init(&argc, &argv);               /********************/
/*	if(err == MPI_SUCCESS)
		printf("Success in initializing MPI!");*/			//xchen_20150407
	comm = MPI_COMM_WORLD;                /********************/
	MPI_Comm_size(comm, &npes);           /********************/
        MPI_Comm_rank(comm, &my_pe);          /********************/		//xchen_20150328_end
//	printf("Processor %d of %d: PIHM MPI!\n", my_pe, npes);			//xchen_20150409
	
	/* Project Input Name */
	if(argc!=2)
		{
		iproj=fopen("projectName.txt","r");
		if(iproj==NULL)
			{
			printf("\t\nUsage ./pihm project_name");
			printf("\t\n         OR              ");
			printf("\t\nUsage ./pihm, and have a file in the current directory named projectName.txt with the project name in it");
			exit(0);
			}
		else
			{
			filename = (char *)malloc(15*sizeof(char));
			fscanf(iproj,"%s",filename);
			}
		}
	else
		{
  		/* get user specified file name in command line */
    		filename = (char *)malloc(strlen(argv[1])*sizeof(char));
		strcpy(filename,argv[1]);
		}

  	/* allocate memory for model data structure */
  	mData = (Model_Data)malloc(sizeof *mData);
 	if(my_pe==0)								//xchen_20150405
		{ 
  		printf("\n ...  PIHM 2.0 MPI is starting ... \n");
 		}
    	
/*	mData->comm = comm;		//xchen_20150410_start
  	mData->npes = npes;
  	mData->my_pe = my_pe;		//xchen_20150410_end
*/
 	/* read in 9 input files with "filename" as prefix */
	read_alloc(filename, mData, &cData); 

#ifdef SUB_SURF_RIV
	#ifdef LAYER2
  		/* problem size */
        	mData->totele=3*mData->NumEle;
		N = mData->totele + 2*mData->NumRiv;
		local_N = 3*NumEleInProc[my_pe] + 2*NumRivInProc[my_pe];                //xchen_20150328
		atol = N_VNew_Parallel(comm,local_N,N);					//xchen_20150406
	#else
  		/* problem size */
                mData->totele=4*mData->NumEle;
		N = mData->totele + 2*mData->NumRiv;
                local_N = 4*NumEleInProc[my_pe] + 2*NumRivInProc[my_pe];                //xchen_20150328
                atol = N_VNew_Parallel(comm,local_N,N); 				//xchen_20150406
		for(j=0;j<NumEleInProc[my_pe];j++)
			{
			NV_Ith_P(atol,j+3*NumEleInProc[my_pe])=cData.abstol[1];
			}
	#endif 	
	for(j=0;j<NumEleInProc[my_pe];j++)
		{
		NV_Ith_P(atol,j)=cData.abstol[0];
		NV_Ith_P(atol,j+NumEleInProc[my_pe])=cData.abstol[1];
		NV_Ith_P(atol,j+2*NumEleInProc[my_pe])=cData.abstol[2];
		}
	for(j=0;j<NumRivInProc[my_pe];j++)
		{
		NV_Ith_P(atol,j+4*NumEleInProc[my_pe])=cData.abstol[0];
		NV_Ith_P(atol,j+4*NumEleInProc[my_pe]+NumRivInProc[my_pe])=cData.abstol[2];
		}
#elif SURF_RIV
	/* problem size*/
	mData->totele=mData->NumEle;
	N = mData->totele+mData->NumRiv;
     	local_N = NumEleInProc[my_pe] + NumRivInProc[my_pe];                //xchen_20150328
	atol = N_VNew_Parallel(comm,local_N,N);                                 //xchen_20150406
        for(j=0;j<NumEleInProc[my_pe];j++)
                {
                NV_Ith_P(atol,j)=cData.abstol[0];
		}	
        for(j=0;j<NumRivInProc[my_pe];j++)
                {
                NV_Ith_P(atol,j+NumEleInProc[my_pe])=cData.abstol[0];
		}
#endif

	mData->DummyY=(realtype *)malloc(N*sizeof(realtype));
  	if(mData->NumProc!=npes)									//xchen_20150328_start
		{
		printf("\n Number of partitions and processors assigned do not match !!");
		exit(1);
		}
 
	/* initial state variable depending on machine*/
	CV_Y = N_VNew_Parallel(comm, local_N, N);
												
//	/* initial state variable depending on machine*/
//  	CV_Y = N_VNew_Serial(N);									//xchen_20150328_end
  
  	/* initialize mode data structure */
  	initialize(filename, mData, &cData, CV_Y);
 
	if(my_pe == 0) 
  		printf("\nSolving ODE system ... \n");
  
  	/* allocate memory for solver */
  	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  	if(cvode_mem == NULL) {printf("CVodeMalloc failed. \n"); return(1);}
  
  	flag = CVodeSetFdata(cvode_mem, mData);  
  	flag = CVodeSetInitStep(cvode_mem,cData.InitStep);
  	flag = CVodeSetStabLimDet(cvode_mem,TRUE);  
  	flag = CVodeSetMaxStep(cvode_mem,cData.MaxStep); 
	flag = CVodeSetMaxNumSteps(cvode_mem, 4000);
  	flag = CVodeMalloc(cvode_mem, f, cData.StartTime, CV_Y, CV_SV, cData.reltol, atol);  
  	flag = CVSpgmr(cvode_mem, PREC_NONE, 0);
  	flag = CVSpgmrSetGSType(cvode_mem, MODIFIED_GS);
 
	if(my_pe==0)                                            //xchen_20150411_start
        	{
        	QFile = (char *)malloc((strlen(filename)+2)*sizeof(char));
        	strcpy(QFile, filename);
        	res_q1_file = fopen(strcat(QFile, "1.q"), "w");
        	}						//xchen_20150411_end  
  	/* set start time */
  	t=expl_t = cData.StartTime;
  	start = clock();
	
/*	if(my_pe==0)								//xchen_20150328
		{
        	for(i=0;i<11;i++)
                	{
                	sprintf(tmpLname[i+cData.totFiles-11],".rivFlx%d",i);
                	}*/
        	/* Open Output Files */
/*        	for(i=0;i<cData.totFiles;i++)
                	{
                	ofn[i] = (char *)malloc((strlen(filename)+12)*sizeof(char));
                	strcpy(ofn[i], filename);
                	Ofile[i]=fopen(strcat(ofn[i], tmpLname[i]),"w");
                	}
  		}*/
	/* start solver in loops */
  	for(i=0; i<cData.NumSteps; i++)
  		{
    		/* inner loops to next output points with ET step size control */
    		while(expl_t < cData.Tout[i+1])
    			{
      			if (expl_t + cData.ETStep >= cData.Tout[i+1])
      				{
        			NextPtr = cData.Tout[i+1];
      				}
      			else
      				{
        			NextPtr = expl_t + cData.ETStep;
      				}
      			StepSize = NextPtr - expl_t; 
     			 
      			/* calculate Interception Storage */
			is_sm_et(expl_t, StepSize, mData);
			expl_t=expl_t+StepSize;
			}
		NextPtr=cData.Tout[i+1];
		if(my_pe==0)                                            //xchen_20150328_start
                	{
			printf("\n Tsteps = %f ",t);
			fflush(stdout);
			} 						//xchen_20150328_end
     		flag = CVode(cvode_mem, NextPtr, CV_Y, &t, CV_NORMAL);  
                while(flag==CV_TOO_MUCH_WORK)
                	{
                        flag = CVode(cvode_mem, NextPtr, CV_Y, &t, CV_NORMAL);
                        printf("\nCV TOO MUCH WORK in PROC %d", my_pe);			//xchen_20150405
			}
		if(my_pe==0)									//xchen_20150328_start
			{
			if(mData->Riv[0].proc-1!=0)
				{
		              	MPI_Recv(&locQ1,1,PVEC_REAL_MPI_TYPE,mData->Riv[0].proc-1,1,comm,&status);
//               		printf("\n locQ = %f was received from %d by %d",locQ,RivOut,my_pe);
//                		fflush(stdout);
			        fprintf(res_q1_file, "\n%f\t%lf", t, locQ1);
              			fflush(res_q1_file);
         			}
         		else
         			{
                 		locQ1=mData->FluxRiv[0][1];
                 		fprintf(res_q1_file, "\n%f\t%lf", t, locQ1);
                 		fflush(res_q1_file);
         			}
		  	}
 		else
    			{
         		if(my_pe==mData->Riv[0].proc-1)
         			{
                 		locQ1=mData->FluxRiv[0][1];
                 		MPI_Send(&locQ1,1,PVEC_REAL_MPI_TYPE,0, 1, comm);
         			}
			}									//xchen_20150328_end
//		PrintData(Ofile,&cData,mData, CV_Y,t);						//xchen_20150328
/*		fprintf(Ofile[10],"\n%lf\t%lf\t",NextPtr,t);
		for(j=0;j<mData->NumRiv;j++)
			{
			fprintf(Ofile[10],"%lf\t",NV_Ith_S(CV_Y,j+3*mData->NumEle));	
			}
		fflush(Ofile[10]);
*/		update(t,mData);
	}
//   	fflush(stdout);
//        getchar();
        /* close output files */								//xchen_20150328_start
 	if(my_pe==0)
   		{
         	fclose(res_q1_file);
   		}										//xchen_20150328_end
	/* Free memory */
  	N_VDestroy_Parallel(CV_Y);								//xchen_20150328
	/* Free integrator memory */
  	CVodeFree(cvode_mem);
  	free(mData);
	
	MPI_Finalize();										//xchen_20140328
  	return 0;
	}

