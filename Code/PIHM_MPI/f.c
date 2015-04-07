/*******************************************************************************
 *-----------------------------------------------------------------------------*
 * File        : f.c   (PIHM v.2.0)                                            *
 * Function    : Model Kernel: Building ODE system for each physical process   *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *-----------------------------------------------------------------------------*
 * Developer of PIHM v.2.0:  Mukesh Kumar (muk139@psu.edu)		       *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *-----------------------------------------------------------------------------*
 * NOTE: f.c has gone a massive revamp (essentially rewritten) since PIHM v.1.0*
 *                                                                             *
 *...........MODFICATIONS/ADDITIONS incorporated in f.c (PIHM v.2.0)...........*
 * a) Surface Flow: 							       *
 *	--> Correction of diffusion wave approximation (calculation of dh/ds)  *
 *              i. Calculation of dh/ds performed using planar slope connecting*
 *                 neighboring centroids				       *
 *              ii.Reflection of elements at boundaries and rivers for dh/ds   *
 *		   calculation
 *	--> Correction of kinematic wave approximation (dh/ds calculation based*
 *	    on elevation only instead of head				       *
 *	--> Correction of gradient for cases with steep change in topography   *
 * b) Subsurface Flow:							       *
 *	--> Addition of macropore phenomena				       *
 *	--> Addition of rectangular cell beneath a river element	       *
 *	--> Implementation of two layered subsurface model(sat/unsat) based on *
 *	Richard's eqn							       *
 *	--> Incorporation of Vertical and Horizontal Anisotropy                *
 *	--> Use of geologic data					       *
 * c) River Flow:							       *
 *	--> Correction of kinematic and diff. wave approximation of SV eqn     *
 *	--> Incorporation of flexible river shapes			       *
 *	--> Separate incorporation of leakage and lateral flow		       *
 *	--> Correction of bank overland flow for extreme cases		       *
 *	--> Addition of aquifer cells below river elements		       *
 * c) Surface/Subsurface Coupling:					       *
 *	--> Implementation of First order coupling through (in/ex)filtration   *
 *		based on head continuity 				       *
 * d) Evaporation:							       *
 *	--> Incorporation of ET from ground/subsurface/vegetation	       *
 *	--> Incorporation of landcover properties for calculation of each ET   *
 *	    component							       *
 * e) Computational:							       *
 *	--> Use of temporary state variables in calculation. Note: Never change* 
 *		core state variables					       *
 * f) Miscellaneous (other advantages realtive to PIHM1.0): No maximum         *
 *    constraint on gw level. Accordingly, no numerical constraints on subsur- *
 *    face flux terms.Faster Implementation. Led to first large scale model    *
 *    application.
 *-----------------------------------------------------------------------------*
 *									       *
 *-----------------------------------------------------------------------------*
 * For questions or comments, please contact 				       *
 *	--> Mukesh Kumar (mk176@duke.edu)				       *
 * This code is free for research purpose only.		       		       *
 * Please provide relevant references if you use this code in your research work*
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * REFERENCES:                                                                 *
 * PIHM2.0:                                                                    *
 *      a) Kumar, M., 2008, "Development and Implementation of a Multiscale,   *
 *      Multiprocess Hydrologic Model". PhD Thesis, Penn State University      *
 *      b) Kumar, M, G.Bhatt & C.Duffy, "Coupling of Data and Processes in     *
 *      a Mesoscale Watershed", Advances in Water Resources (submitted)        *
 * PIHM1.0:                                                                    *
 *      a) Qu, Y., 2005, "An Integrated hydrologic model for multiproces       *
 *      simulation using semi-discrete finite volume approach".PhD Thesis, PSU *
 *      b) Qu, Y. & C. Duffy, 2007, "A semidiscrete finite volume formulation  *
 *      for multiprocess watershed simulation". Water Resources Research       *
 *******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include "nvector_serial.h"								//xchen_20150329_start
#include "nvector_parallel.h"	
#include "mpi.h"									//xchen_20150329_end
#include "sundialstypes.h"   

#include "f_functions.h"
#define MAXP MaxInProc
#define MAX4P 10*MaxInProc		/*EleSeq: 4Ele	RivSeq: 2*Riv  RivEleSeq: 4*Riv  EleRivSeq: 6*Ele  IN TOTAL: 10Ele+6Riv*/	//xchen_20150405

MPI_Status status;									//xchen_20150329
int locI;										//xchen_20150329
void f(realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS)
	{
  	int i,j,k,l,m,ieleBC,inabr,inabr_left,inabr_right,inabr_plus_2ele,ileft,iright,i_plus_ele,i_plus_2ele,i_plus_3ele,iright_plus_2ele, ileft_plus_2ele, i_plus_totele, i_plus_totele1riv;
  	int nly;									//nly:Total number of layers (surface + subsurface) 	xchen_20140402
	realtype tmp;									//temporaty value   	xchen_20150402
	realtype Delta, Gamma,tmpVar,tmpVar1;
  	realtype Rn, T, Vel, RH, VP,P,LAI,zero_dh,cnpy_h,rl,r_a,r_s,alpha_r,f_r,eta_s,beta_s,Rmax;
  	realtype tempVar1;
	realtype Avg_Sf,Distance;
   	realtype Cwr,Perem, Perem_down,Avg_Rough,Avg_Perem,Avg_Y_Riv,Wid,Wid_down,Avg_Wid;
   	realtype nabrAqDepth,AquiferDepth, Deficit,elemSatn,satKfunc,effKI,effK_unsat,TotalY_Ele,TotalY_Ele_down,TotalY_unsat;
  	realtype *Y, *DY,*DummyY;						//xchen_20150403 DummyY, DummyDY store the local value in a global manner (sizes equal to the global sizes, and only assign value of grids in local machine)  Y and DY is the local Y value and only have # of elements = NumEleInProc + NumRivInProc 
  	realtype comm_buff[MAX4P];							//xchen_20150329
	Model_Data MD;
  	Y = NV_DATA_P(CV_Y);								//xchen_20150329
  	DY = NV_DATA_P(CV_Ydot);							//xchen_20150329
  	MD = (Model_Data) DS;

	/* Initialization of temporary state variables */ 
       	MD->totriv=2*MD->NumRiv;//NumRivInProc[my_pe];
/*	#ifdef SURF_RIV
        MD->totriv=MD->NumRiv;//NumRivInProc[my_pe];
	#endif		*/								//xchen_20150402	SURF_RIV is not considered in MPI since it is super fast
	#ifdef LAYER2
	MD->totele = 3*MD->NumEle;//3*NumEleInProc[my_pe];
	nly = 3;			//Total number of layers (surface + subsurface)
	#elif LAYER3
        MD->totele = 4*MD->NumEle;//4*NumEleInProc[my_pe];
        nly = 4;
	#endif

	/*Dummy and DummyDY are both local variable, but with the size same as global variable*/	
	DummyY=calloc((MD->totele+MD->totriv),sizeof(realtype));			//xchen_20150402
//    	DummyDY=calloc((MD->totele+MD->totriv),sizeof(realtype));			//xchen_20150402

	/*Initialization to 0 		xchen_20150405*/
	for(i=0; i<nly*NumEleInProc[my_pe]+2*NumRivInProc[my_pe]; i++)
                {
                DY[i]=0;
                }
	for(i=0;i<MD->NumRiv;i++)
                {
                MD->FluxRiv[i][0]=0;
                MD->FluxRiv[i][10]=0;
                }
	
/* Initialization states variables: Element then River  */  
/*****************************************************************************************************************/	
	locI=0;										//xchen_201504029_start
	for(m=0; m<EleInProcCounter[my_pe]; m++)
  		{
		i= EleInProc[my_pe][m]-1;
		tmp = 0;
		DummyY[i]=(Y[locI]>=0)?Y[locI]:0;
		for(j=1;j<nly;j++)
			{
			DummyY[i+j*MD->NumEle]=(Y[locI+j*NumEleInProc[my_pe]]>=0)?Y[locI+j*NumEleInProc[my_pe]]:0;
			tmp = tmp + Y[locI+j*NumEleInProc[my_pe]];
			}
		if(tmp>MD->Ele[i].zmax-MD->Ele[i].zmin)
			{
		#ifdef LAYER2
			if((Y[locI+NumEleInProc[my_pe]]>0)&&(Y[locI+NumEleInProc[my_pe]]<MD->Ele[i].zmax-MD->Ele[i].zmin)&&(Y[locI+2*NumEleInProc[my_pe]]>0)&&(Y[locI+2*NumEleInProc[my_pe]]<MD->Ele[i].zmax-MD->Ele[i].zmin))
				{
				DummyY[i+MD->NumEle]=1.0*((MD->Ele[i].zmax-MD->Ele[i].zmin)-Y[locI+2*NumEleInProc[my_pe]]);
				}	
			else if(Y[locI+NumEleInProc[my_pe]]<0)
				{
				//DummyY[i+MD->NumEle]=0;
				DummyY[i+2*MD->NumEle]=MD->Ele[i].zmax - MD->Ele[i].zmin;
				}
			else
				{
				DummyY[i+MD->NumEle]=MD->Ele[i].zmax - MD->Ele[i].zmin;	//do not understand this line (xchen_20150402)
				}
		#elif LAYER3	
			if((Y[locI+NumEleInProc[my_pe]]>0)&&(Y[locI+NumEleInProc[my_pe]]<MD->Ele[i].zmax-MD->Ele[i].zmin)&&(Y[locI+2*NumEleInProc[my_pe]]>0)&&(Y[locI+2*NumEleInProc[my_pe]]<MD->Ele[i].zmax-MD->Ele[i].zmin)&&(Y[locI+3*NumEleInProc[my_pe]]>0)&&(Y[locI+3*NumEleInProc[my_pe]]<MD->Ele[i].zmax-MD->Ele[i].zmin))
				{
				///////////////////To be added/////////////			//xchen_20150405
				}
		#endif
			}
		locI++;
		}
	locI=0;
    	for(m=0; m<RivInProcCounter[my_pe]; m++)
    		{
		i=RivInProc[my_pe][m]-1;
		if(MD->Riv[i].proc-1==my_pe)
			{
			DummyY[i+nly*MD->NumEle]=(Y[locI+nly*NumEleInProc[my_pe]]>=0)?Y[locI+nly*NumEleInProc[my_pe]]:0;
			locI++;
			////////DO not understant this part, why here checking MD->Riv[i].proc-1==my_pe? why assigning element???/////
			}
		else
			{
			DummyY[i+nly*MD->NumEle]=0;
			}
		}
                
		
/* Communication  */  
/* For land element:	communicated all the neighboring cells
   For river element:	communicated all the down stream river element*/		//xchen_20150404
/*****************************************************************************************************************/
	for(i=0;i<MD->NumProc;i++)
		{
		if(my_pe!=i)
			{
			if(my_pe>i)
				{
				for(j=0;j<EleSeqCounter[my_pe][i];j++)
					{
					for(k=0;k<nly;k++)
						{
						comm_buff[j+k*EleSeqCounter[my_pe][i]]=DummyY[EleNabrSequence[my_pe][i][j]+k*MD->NumEle];
						}
					}
				for(j=0;j<RivSeqCounter[my_pe][i];j++)
					{
					comm_buff[j+nly*EleSeqCounter[my_pe][i]]=DummyY[nly*MD->NumEle+RivNabrSequence[my_pe][i][j]];
					comm_buff[j+RivSeqCounter[my_pe][i]+nly*EleSeqCounter[my_pe][i]]=DummyY[nly*MD->NumEle+MD->NumRiv+RivNabrSequence[my_pe][i][j]];
					}
				for(j=0;j<RivEleSeqCounter[my_pe][i];j++)
					{
					comm_buff[j+2*RivSeqCounter[my_pe][i]+nly*EleSeqCounter[my_pe][i]]=DummyY[RivEleNabrSequence[my_pe][i][j]];
					//xchen_20150404 flow exchange between surface element to river 
					comm_buff[j+RivEleSeqCounter[my_pe][i]+2*RivSeqCounter[my_pe][i]+nly*EleSeqCounter[my_pe][i]]=DummyY[2*MD->NumEle+RivEleNabrSequence[my_pe][i][j]];					
					//xchen_20150404 flow exchange between groundwater cells to river
					}
				for(j=0;j<EleRivSeqCounter[my_pe][i];j++)
                                        {
                                        comm_buff[j+2*RivEleSeqCounter[my_pe][i]+2*RivSeqCounter[my_pe][i]+nly*EleSeqCounter[my_pe][i]]=DummyY[nly*MD->NumEle+EleRivNabrSequence[my_pe][i][j]];                                  
					//xchen_20150404 flow exchange between river layer1 and surface element
                                        comm_buff[j+EleRivSeqCounter[my_pe][i]+2*RivEleSeqCounter[my_pe][i]+2*RivSeqCounter[my_pe][i]+nly*EleSeqCounter[my_pe][i]]=DummyY[nly*MD->NumEle+MD->NumRiv+EleRivNabrSequence[my_pe][i][j]];             
					//xchen_20150404 flow exchange between river layer2 and groundwater element 
                                        }
				if(EleSeqCounter[my_pe][i]+RivSeqCounter[my_pe][i]+RivEleSeqCounter[my_pe][i]+EleRivSeqCounter[my_pe][i]>0)
					{
					MPI_Send(&comm_buff[0],2*EleRivSeqCounter[my_pe][i]+2*RivEleSeqCounter[my_pe][i]+2*RivSeqCounter[my_pe][i]+nly*EleSeqCounter[my_pe][i], PVEC_REAL_MPI_TYPE,i,0, comm);
/*                     			printf("\nSent %d data from proc %d to proc %d",RivSeqCounter[my_pe][i]+2*EleSeqCounter[my_pe][i],my_pe,i);
 	                                fflush(stdout); */
					MPI_Recv(&comm_buff[0],2*EleRivSeqCounter[i][my_pe]+2*RivEleSeqCounter[i][my_pe]+2*RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe],PVEC_REAL_MPI_TYPE,i,0,comm,&status);
/*                      		printf("\nReceived %d data from proc %d by proc %d",RivSeqCounter[i][my_pe]+2*EleSeqCounter[i][my_pe],i,my_pe);
                                  	fflush(stdout); */
					}
				for(j=0;j<EleSeqCounter[i][my_pe];j++)
					{
					for(k=0;k<nly;k++)
						{
						DummyY[EleNabrSequence[i][my_pe][j]+k*MD->NumEle]=comm_buff[j+k*EleSeqCounter[i][my_pe]]<0?0:comm_buff[j+k*EleSeqCounter[i][my_pe]];
/*                     				printf("\n %d %f",EleNabrSequence[i][my_pe][j],comm_buff[j],comm_buff[j+EleSeqCounter[my_pe][i]]);
                              			fflush(stdout); */	
						}
					}
				for(j=0;j<RivSeqCounter[i][my_pe];j++)
					{
					DummyY[nly*MD->NumEle+RivNabrSequence[i][my_pe][j]]=comm_buff[j+nly*EleSeqCounter[i][my_pe]]<0?0:comm_buff[j+nly*EleSeqCounter[i][my_pe]];
					DummyY[nly*MD->NumEle+MD->NumRiv+RivNabrSequence[i][my_pe][j]]=comm_buff[j+RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]]<0?0:comm_buff[j+RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]];
 //                                	printf("\nTo Recevie from %d by %d val %f",i,my_pe,DummyY[3*MD->NumEle+RivNabrSequence[i][my_pe][j]]);
 //                                	fflush(stdout);
					}
				for(j=0;j<RivEleSeqCounter[i][my_pe];j++)
					{
					DummyY[RivEleNabrSequence[i][my_pe][j]]=comm_buff[j+2*RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]]<0?0:comm_buff[j+2*RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]];
					DummyY[2*MD->NumEle+RivEleNabrSequence[i][my_pe][j]]=comm_buff[j+RivEleSeqCounter[i][my_pe]+2*RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]]<0?0:comm_buff[j+RivEleSeqCounter[i][my_pe]+2*RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]];		
					}
				for(j=0;j<EleRivSeqCounter[i][my_pe];j++)
                                        {   
                                        DummyY[nly*MD->NumEle+EleRivNabrSequence[i][my_pe][j]]=comm_buff[j+2*RivEleSeqCounter[i][my_pe]+2*RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]]<0?0:comm_buff[j+2*RivEleSeqCounter[i][my_pe]+2*RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]];
					DummyY[nly*MD->NumEle+MD->NumRiv+EleRivNabrSequence[i][my_pe][j]]=comm_buff[j+EleRivSeqCounter[i][my_pe]+2*RivEleSeqCounter[i][my_pe]+2*RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]]<0?0:comm_buff[j+EleRivSeqCounter[i][my_pe]+2*RivEleSeqCounter[i][my_pe]+2*RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]];                       
                                        }	
				}
			else
				{
				if(EleSeqCounter[i][my_pe]+RivSeqCounter[i][my_pe]+RivEleSeqCounter[i][my_pe]+EleRivSeqCounter[i][my_pe]>0)
                        		{
                        		MPI_Recv(&comm_buff[0],2*EleRivSeqCounter[i][my_pe]+2*RivEleSeqCounter[i][my_pe]+2*RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe],PVEC_REAL_MPI_TYPE,i,0,comm,&status);
                /*      		printf("\nReceived %d data from proc %d by proc %d",RivSeqCounter[i][my_pe]+2*EleSeqCounter[i][my_pe],i,my_pe);        
                        		fflush(stdout); */
                        		}    
                        	for(j=0;j<EleSeqCounter[i][my_pe];j++)
                        		{
				    	for(k=0;k<nly;k++)
						{
						DummyY[EleNabrSequence[i][my_pe][j]+k*MD->NumEle]=comm_buff[j+k*EleSeqCounter[i][my_pe]]<0?0:comm_buff[j+k*EleSeqCounter[i][my_pe]];                                
                        			}
					}
                        	for(j=0;j<RivSeqCounter[i][my_pe];j++)
                        		{    
                                	DummyY[nly*MD->NumEle+RivNabrSequence[i][my_pe][j]]=comm_buff[j+nly*EleSeqCounter[i][my_pe]]<0?0:comm_buff[j+nly*EleSeqCounter[i][my_pe]];                               
                                	DummyY[nly*MD->NumEle+MD->NumRiv+RivNabrSequence[i][my_pe][j]]=comm_buff[j+RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]]<0?0:comm_buff[j+nly*EleSeqCounter[i][my_pe]];                               
//                                	printf("\nTo Recevie from %d by %d val %f",i,my_pe,DummyY[3*MD->NumEle+RivNabrSequence[i][my_pe][j]]);
//                                	fflush(stdout);
                        		}
                        	for(j=0;j<RivEleSeqCounter[i][my_pe];j++)
                        		{    
                                	DummyY[RivEleNabrSequence[i][my_pe][j]]=comm_buff[j+2*RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]]<0?0:comm_buff[j+2*RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]];					
					DummyY[2*MD->NumEle+RivEleNabrSequence[i][my_pe][j]]=comm_buff[j+RivEleSeqCounter[i][my_pe]+2*RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]]<0?0:comm_buff[j+RivEleSeqCounter[i][my_pe]+2*RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]];		
                        		}
				for(j=0;j<EleRivSeqCounter[i][my_pe];j++)
                                        {   
					DummyY[nly*MD->NumEle+EleRivNabrSequence[i][my_pe][j]]=comm_buff[j+2*RivEleSeqCounter[i][my_pe]+2*RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]]<0?0:comm_buff[j+2*RivEleSeqCounter[i][my_pe]+2*RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]];
                                        DummyY[nly*MD->NumEle+MD->NumRiv+EleRivNabrSequence[i][my_pe][j]]=comm_buff[j+EleRivSeqCounter[i][my_pe]+2*RivEleSeqCounter[i][my_pe]+2*RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]]<0?0:comm_buff[j+EleRivSeqCounter[i][my_pe]+2*RivEleSeqCounter[i][my_pe]+2*RivSeqCounter[i][my_pe]+nly*EleSeqCounter[i][my_pe]];
                                        }	
                        	for(j=0;j<EleSeqCounter[my_pe][i];j++)
                        		{    
					for(k=0;k<nly;k++)
						{
                                		comm_buff[j+k*EleSeqCounter[my_pe][i]]=DummyY[EleNabrSequence[my_pe][i][j]+k*MD->NumEle]; 
                       				}
					}             
                        	for(j=0;j<RivSeqCounter[my_pe][i];j++)
                        		{    
                                	comm_buff[j+nly*EleSeqCounter[my_pe][i]]=DummyY[nly*MD->NumEle+RivNabrSequence[my_pe][i][j]];
                                	comm_buff[j+RivSeqCounter[my_pe][i]+nly*EleSeqCounter[my_pe][i]]=DummyY[nly*MD->NumEle+MD->NumRiv+RivNabrSequence[my_pe][i][j]];
//                                	printf("\nTo Send from %d to %d val %f",my_pe,i,comm_buff[j+3*EleSeqCounter[my_pe][i]]);
//                                	fflush(stdout);
					}
                        	for(j=0;j<RivEleSeqCounter[my_pe][i];j++)
                        		{
                                	comm_buff[j+2*RivSeqCounter[my_pe][i]+nly*EleSeqCounter[my_pe][i]]=DummyY[RivEleNabrSequence[my_pe][i][j]];	
                                	comm_buff[j+RivEleSeqCounter[my_pe][i]+2*RivSeqCounter[my_pe][i]+nly*EleSeqCounter[my_pe][i]]=DummyY[2*MD->NumEle+RivEleNabrSequence[my_pe][i][j]];								
                        		}
				for(j=0;j<EleRivSeqCounter[my_pe][i];j++)
                                        {   
                                        comm_buff[j+2*RivEleSeqCounter[my_pe][i]+2*RivSeqCounter[my_pe][i]+nly*EleSeqCounter[my_pe][i]]=DummyY[nly*MD->NumEle+EleRivNabrSequence[my_pe][i][j]];        
                                        //xchen_20150404 flow exchange between river layer1 and surface element
                                        comm_buff[j+EleRivSeqCounter[my_pe][i]+2*RivEleSeqCounter[my_pe][i]+2*RivSeqCounter[my_pe][i]+nly*EleSeqCounter[my_pe][i]]=DummyY[nly*MD->NumEle+MD->NumRiv+EleRivNabrSequence[my_pe][i][j]];        
                                        //xchen_20150404 flow exchange between river layer2 and groundwater element 
                                        } 
                        	if(EleSeqCounter[my_pe][i]+RivSeqCounter[my_pe][i]+RivEleSeqCounter[my_pe][i]+EleRivSeqCounter[my_pe][i]>0)
                        		{
                        		MPI_Send(&comm_buff[0],2*EleRivSeqCounter[my_pe][i]+2*RivEleSeqCounter[my_pe][i]+2*RivSeqCounter[my_pe][i]+nly*EleSeqCounter[my_pe][i], PVEC_REAL_MPI_TYPE,i, 0, comm);
                /*      		printf("\nSent %d data from proc %d to proc %d",RivSeqCounter[my_pe][i]+2*EleSeqCounter[my_pe][i],my_pe,i);            
                        		fflush(stdout); */
                        		}
				}
			}	
		}
/*****************************************************************************************************************/
//  printf("\n Just AFTER communication %d",my_pe);								//xchen_20150402_end
	locI=0;
	for(m=0; m<EleInProcCounter[my_pe]; m++)
  		{
    		i=EleInProc[my_pe][m]-1;
		if(MD->Ele[i].procNo-1==my_pe)
			{
   	#ifdef DIFFUSION
			for(j=0;j<3;j++)
				{
                		inabr=MD->Ele[i].nabr[j]-1;
                		ieleBC=-(MD->Ele[i].BC[j]/4)-1;
			      	MD->Ele[i].surfH[j]=(inabr>-1)?((-ieleBC>0)?(MD->Ele[inabr].zmax+DummyY[inabr]):((DummyY[ieleBC+MD->totele]>MD->Riv[ieleBC].depth)?MD->Riv[ieleBC].zmin+DummyY[ieleBC+MD->totele]:MD->Riv[ieleBC].zmax)):((MD->Ele[i].BC[j]!=1)?(MD->Ele[i].zmax+DummyY[i]):Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j])-1], t));
                        	}
                	MD->Ele[i].dhBYdx=-1*(MD->Ele[i].surfY[2]*(MD->Ele[i].surfH[1]-MD->Ele[i].surfH[0])+MD->Ele[i].surfY[1]*(MD->Ele[i].surfH[0]-MD->Ele[i].surfH[2])+MD->Ele[i].surfY[0]*(MD->Ele[i].surfH[2]-MD->Ele[i].surfH[1]))/(MD->Ele[i].surfX[2]*(MD->Ele[i].surfY[1]-MD->Ele[i].surfY[0])+MD->Ele[i].surfX[1]*(MD->Ele[i].surfY[0]-MD->Ele[i].surfY[2])+MD->Ele[i].surfX[0]*(MD->Ele[i].surfY[2]-MD->Ele[i].surfY[1]));
                	MD->Ele[i].dhBYdy=-1*(MD->Ele[i].surfX[2]*(MD->Ele[i].surfH[1]-MD->Ele[i].surfH[0])+MD->Ele[i].surfX[1]*(MD->Ele[i].surfH[0]-MD->Ele[i].surfH[2])+MD->Ele[i].surfX[0]*(MD->Ele[i].surfH[2]-MD->Ele[i].surfH[1]))/(MD->Ele[i].surfY[2]*(MD->Ele[i].surfX[1]-MD->Ele[i].surfX[0])+MD->Ele[i].surfY[1]*(MD->Ele[i].surfX[0]-MD->Ele[i].surfX[2])+MD->Ele[i].surfY[0]*(MD->Ele[i].surfX[2]-MD->Ele[i].surfX[1]));
   	#endif
	/* Lateral Flux Calculation between Triangular elements Follows  */
        	i_plus_ele=i+MD->NumEle;
        	i_plus_2ele=i+2*MD->NumEle;
		AquiferDepth=(MD->Ele[i].zmax-MD->Ele[i].zmin);
    		for(j=0; j<3; j++)
    			{
      			if(MD->Ele[i].nabr[j] > 0)
      				{
               			inabr=MD->Ele[i].nabr[j]-1;   
			#ifdef SUB_SURF_RIV
                                if(MD->Ele[i].Calc[j]>0)
                                MD->FluxSub[i][j] = - MD->FluxSub[inabr][MD->Ele[i].Calc[j]-1];
                                else
                                {
                        	inabr_plus_2ele=inabr+2*MD->NumEle;
			        /***************************************************************************/
			        /* Subsurface Lateral Flux Calculation between Triangular elements Follows */
			        /***************************************************************************/
                    		GradCalc(DummyY[inabr_plus_2ele], MD->Ele[inabr].zmin, DummyY[i_plus_2ele], MD->Ele[i].zmin, MD->Ele[i].Dist[j],0,0,0,0,ELE_ELE);
        			/* take care of macropore effect */
                    		nabrAqDepth=(MD->Ele[inabr].zmax-MD->Ele[inabr].zmin);
                    		avgKH(MD->Ele[i].Macropore,DummyY[i_plus_2ele],AquiferDepth,MD->Ele[i].macD,MD->Ele[i].macKsatH,MD->Ele[i].vAreaF,MD->Ele[i].KsatH,MD->Ele[inabr].Macropore,DummyY[inabr_plus_2ele],nabrAqDepth,MD->Ele[inabr].macD,MD->Ele[inabr].macKsatH,MD->Ele[inabr].vAreaF,MD->Ele[inabr].KsatH);
        			/* groundwater flow modeled by Darcy's law */
        			MD->FluxSub[i][j] = Avg_Ksat*Grad_Y*Avg_Y*MD->Ele[i].edge[j];
                               }
			#endif
		#ifndef NO_UNSAT
			    	/***************************************************************************/
			    	/* Surface Lateral Flux Calculation between Triangular elements Follows    */
			    	/***************************************************************************/
                        if(MD->Ele[i].Calc[j]>0)
                                MD->FluxSurf[i][j] = - MD->FluxSurf[inabr][MD->Ele[i].Calc[j]-1];
                         else
                        	{
                 	#ifdef DIFFUSION
                    		GradCalc(DummyY[inabr], MD->Ele[inabr].zmax, DummyY[i], MD->Ele[i].zmax, MD->Ele[i].Dist[j],0,0,0,0,ELE_ELE);
                    		Avg_Sf=sqrt(pow(MD->Ele[i].dhBYdx,2)+pow(MD->Ele[i].dhBYdy,2));
                    		Avg_Sf=(Avg_Sf>EPS_5)?Avg_Sf:EPS_5;
                    	#else
                    		GradCalc(DummyY[inabr], MD->Ele[inabr].zmax-DummyY[inabr], DummyY[i], MD->Ele[i].zmax-DummyY[i], MD->Ele[i].Dist[j],0,0,0,0,ELE_ELE);
                    		Avg_Sf=(fabs(Grad_Y)>EPS_5)?Grad_Y:EPS_5;
                    	#endif
                    		/* Weighting needed */
        			Avg_Rough = 0.5*(MD->Ele[i].Rough + MD->Ele[inabr].Rough);
        			CrossA = Avg_Y*MD->Ele[i].edge[j];
                    		/* INCLUDE CROSSADOWN */
				OverlandFlow(MD->FluxSurf,i,j, Avg_Y,Grad_Y,Avg_Sf,CrossA,Avg_Rough);
	                	}
            	#endif
      				}
			/************************************************/
			/* Boundary condition Flux Calculations Follows */
			/************************************************/
      			else
      				{
        			/*  No flow (natural) boundary condition is default */
        			if(MD->Ele[i].BC[j] == 0)
        				{
          				MD->FluxSurf[i][j] = 0;
          				MD->FluxSub[i][j] = 0;
        				}
        			else if(MD->Ele[i].BC[j] == 1)	/* Note: ideally different boundary conditions need to be incorporated	for surf and subsurf respectively */
					/* Note: the formulation assumes only dirichlet TS right now */
        				{
          				MD->FluxSurf[i][j] = 0;	/* Note the assumption here is no flow for surface*/ 
                        		effK=effKH(MD->Ele[i].Macropore,DummyY[i_plus_2ele],AquiferDepth,MD->Ele[i].macD,MD->Ele[i].macKsatH,MD->Ele[i].vAreaF,MD->Ele[i].KsatH);
           				Avg_Ksat = effK;
                        		GradCalc(Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j])-1], t)- MD->Ele[i].zmin, MD->Ele[i].zmin, DummyY[i_plus_2ele], MD->Ele[i].zmin, MD->Ele[i].Dist[j],0,0,0,0,ELE_ELE);
           				MD->FluxSub[i][j] = Avg_Ksat*Grad_Y*Avg_Y*MD->Ele[i].edge[j];
          				}
          			else       /* Neumann BC (Note: MD->Ele[i].BC[j] value have to be = 2+(index of neumann boundary TS)*/
          				{
            				MD->FluxSurf[i][j] = Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j])-1], t);
            				MD->FluxSub[i][j] = Interpolation(&MD->TSD_EleBC[(-MD->Ele[i].BC[j])-1], t);
          				}
      				}
    			}
#ifndef NO_UNSAT
#ifdef SUB_SURF_RIV		
		/* Calculation of ET0, ET1, ViR, Recharge, RechargeI is performed in fluxCalc_Ele() */
		fluxCalc_Ele(MD,i,t,DummyY);							//xchen_20150405	add Dummy into the function 
   	#ifdef LAYER3
                i_plus_3ele=i+3*MD->NumEle;							//xchen_20150403_start Change al DY[i] into DY[locI]
                DY[locI+3*NumEleInProc[my_pe]] = MD->RechargeI[i]-MD->Recharge[i];
        #endif
		DY[locI] = MD->EleNetPrep[i] - MD->EleViR[i]-((DummyY[i]<EPS_3)?0:MD->EleET[i][2]);
		DY[locI+NumEleInProc[my_pe]] = MD->EleViR[i]-MD->RechargeI[i]-((DummyY[i]<EPS_3)?MD->EleET[i][2]:0);
		DY[locI+2*NumEleInProc[my_pe]]= MD->Recharge[i];
		if(DummyY[i_plus_2ele]>AquiferDepth-MD->Ele[i].RzD)
			{
			DY[locI+2*NumEleInProc[my_pe]]=DY[locI+2*NumEleInProc[my_pe]]-MD->EleET[i][1];
			}
		else 
			{
		#ifdef LAYER2
			DY[locI+NumEleInProc[my_pe]] = DY[locI+NumEleInProc[my_pe]]-MD->EleET[i][1];
		#endif 
/*	xc 20130529  	The following block is to incorporate the second unsaturated zone into participating ET1 calculation. 
			Purpose is to extract ET1 from the 2 unsaturated zones based on water availability.				*/		
#ifdef LAYER3
			if(s1<=EPS_4)
				{
				if (s2<=EPS_4)
					{
					}
				else
					{
					DY[locI+3*NumEleInProc[my_pe]]=DY[locI+3*NumEleInProc[my_pe]]-MD->EleET[i][1];
					}
				}
			else
				{
				if (s2<=EPS_4)
					{
					DY[locI+NumEleInProc[my_pe]] = DY[locI+NumEleInProc[my_pe]]-MD->EleET[i][1];
					}
				else
					{
					tmpVar1=s1/(s1+s2);
                        		DY[locI+NumEleInProc[my_pe]] = DY[locI+NumEleInProc[my_pe]]-tmpVar1*MD->EleET[i][1];
                        		DY[locI+3*NumEleInProc[my_pe]] = DY[locI+3*NumEleInProc[my_pe]]-(1-tmpVar1)*MD->EleET[i][1];
					}
				}
		#endif
			}
#elif SURF_RIV
		DY[locI] = MD->EleNetPrep[i];
#endif
#endif
  			locI++;
			}
		} 
	/* Lateral Flux Calculation between River-River and River-Triangular elements Follows */  
	//for both river in this processor & element in this processor with neighboring river in other processor(s) xchen_20150405_start
	/* river in this processor */
	for(m=0; m<RivInProcCounter[my_pe]; m++)
		{
		i=RivInProc[my_pe][m]-1;
		if(MD->Riv[i].proc-1==my_pe)
			{
			fluxCalc_Riv(MD,i,t,DummyY);
			}
		}
	/* element in this processor with neighboring river in other processor(s)  */
	for(m=0; m<EleInProcCounter[my_pe]; m++)
        	{   
        	k=EleInProc[my_pe][m]-1;
        	if(MD->Ele[k].procNo-1==my_pe)
                	{  
                	for(i=0;i<MD->NumProc;i++)
                        	{   
                        	if(my_pe!=i && EleRivSeqCounter[my_pe][i]>0)
                                	{   
                                	for(l=0;l<EleRivSeqCounter[my_pe][i];l++)
                                        	{   
                                        	j=EleRivNabrSequence[my_pe][i][l];
						fluxCalc_Riv_light(MD,j,k,t,DummyY); 
						}   
                                	}   
                        	}   
                	}   
        	}
	//for both river in this processor & element in this processor with neighboring river in other processor(s) xchen_20150405_end		

	locI = 0;
	for(m=0; m<EleInProcCounter[my_pe]; m++)
		{
		i=EleInProc[my_pe][m]-1;
        	if(MD->Ele[i].procNo-1==my_pe)
			{	
			i=EleInProc[my_pe][m]-1;   
        		i_plus_ele=i+MD->NumEle;
        		i_plus_2ele=i+2*MD->NumEle;
     			for(j=0; j<3; j++)
      				{
        			DY[locI] =  DY[locI] - MD->FluxSurf[i][j]/MD->Ele[i].area;
			#ifdef SUB_SURF_RIV
        			MD->FluxSource[i] =(MD->Ele[i].source>0)?Interpolation(&MD->TSD_Source[MD->Ele[i].source - 1], t):0;           //xchen_20141231
                        	if (DummyY[i_plus_2ele]>0)                                          //xchen_20141021
                        		{
                        		DY[locI+2*NumEleInProc[my_pe]] = DY[locI+2*NumEleInProc[my_pe]] - MD->FluxSub[i][j]/MD->Ele[i].area - MD->FluxSource[i];      //xchen_20141006
                        		}
			#endif
      				}
			DY[locI]=DY[locI]/(UNIT_C);
		#ifdef NO_UNSAT
			DY[locI+2*NumEleInProc[my_pe]] = DY[locI+2*NumEleInProc[my_pe]] + MD->ElePrep[i];
		#endif
		#ifdef SUB_SURF_RIV
	      		DY[locI+NumEleInProc[my_pe]] = DY[locI+NumEleInProc[my_pe]]/(MD->Ele[i].Porosity*UNIT_C);
      			DY[locI+2*NumEleInProc[my_pe]] = DY[locI+2*NumEleInProc[my_pe]]/(MD->Ele[i].Porosity*UNIT_C);
		#ifdef LAYER3
        		i_plus_3ele=i+3*MD->NumEle;
      			DY[locI+3*NumEleInProc[my_pe]] = DY[locI+3*NumEleInProc[my_pe]]/(MD->Ele[i].Porosity*UNIT_C);
		#endif
		#endif
    			locI++;
			}
		}

	locI=0;
	for(m=0; m<RivInProcCounter[my_pe]; m++)
    		{
		i=RivInProc[my_pe][m]-1;
		if(MD->Riv[i].proc-1==my_pe)
			{
        		i_plus_ele=i+MD->NumEle;
        		i_plus_totele=i+MD->totele;
        		i_plus_totele1riv=i+MD->totele+MD->NumRiv;
			for(j=0;j<=6;j++)
				{
				/* Note the limitation due to d(v)/dt=a*dy/dt+y*da/dt for CS other than rectangle */
				DY[locI+nly*NumEleInProc[my_pe]] = DY[locI+nly*NumEleInProc[my_pe]]-MD->FluxRiv[i][j]/(MD->Riv[i].Length*CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->Riv[i].depth,MD->Riv[i].coeff,3));
				}
			DY[locI+nly*NumEleInProc[my_pe]] = DY[locI+nly*NumEleInProc[my_pe]]/(UNIT_C);
	#ifdef SUB_SURF_RIV
			DY[locI+nly*NumEleInProc[my_pe]+NumRivInProc[my_pe]] = DY[locI+nly*NumEleInProc[my_pe]+NumRivInProc[my_pe]] -MD->FluxRiv[i][7] -MD->FluxRiv[i][8]-MD->FluxRiv[i][9] -MD->FluxRiv[i][10]+MD->FluxRiv[i][6];
      			DY[locI+nly*NumEleInProc[my_pe]+NumRivInProc[my_pe]] = DY[locI+nly*NumEleInProc[my_pe]+NumRivInProc[my_pe]]/(MD->Ele[i_plus_ele].Porosity*MD->Riv[i].Length*CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->Riv[i].depth,MD->Riv[i].coeff,3)*UNIT_C);
	#endif
    			locI++;
			}
//	printf("\nf  %lf, %lf, %lf, %lf, %lf, %lf",MD->FluxRiv[0][0],MD->FluxRiv[1][0],MD->FluxRiv[0][1],MD->FluxRiv[1][1],MD->FluxRiv[1][2],MD->FluxRiv[1][3]);
		}
	free(DummyY);
	}
