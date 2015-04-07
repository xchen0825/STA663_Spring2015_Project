/*******************************************************************************
 * File        : print.c	                                               *
 * Version     : Nov, 2007 (2.0)                                               *
 * Function    : print out model results output files                          *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *                                                                             *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0............................*
 * a) This file is downgraded from Version 1.0, as no ancillary results are    *
 *    being output			                                       *
 * b) Only state variables and flux to/in/accross river and its bed are being  *
 *    output							               *
 * c) Addition of Average Function to output average variables at regular time *
 *    intervals								       *
 *******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cvode.h" 
#include "cvdense.h"
#include "print_functions.h"
/*Temporal average of State vectors */
void avgResults_NV(int fileCounter,FILE **fpin,realtype *varToPrint,N_Vector tmpNV,Model_Data tmpDS,int tmpIntv,realtype tmpt,int tmpCounter,int modIntv)
        {
        int j,tmpNumObj;
	tmpNumObj=(fileCounter+tmpCounter)<12?tmpDS->NumEle:tmpDS->NumRiv;
	switch(fileCounter+tmpCounter)
		{
		case 0:
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->EleNetPrep[j];
				}
                        break;
		case 1:
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->EleIS[j];
				}
                        break;
		case 2:
	#ifdef SUB_SURF_RIV
		case 3:
		case 4:
	#endif
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+NV_Ith_S(tmpNV,j+(fileCounter-2)*tmpDS->NumEle);
				}
                        break;
#ifdef SUB_SURF_RIV
		case 5:
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+tmpDS->EleET[j][0];
				}
                        break;
		case 6:
			for(j=0;j<tmpNumObj;j++) 
				{
				forcingInitialize(tmpDS,j,tmpt);
				et1Calc(tmpDS,j,tmpNV);
				varToPrint[j]=varToPrint[j]+tmpDS->EleET[j][1];
				}
                        break;
		case 7:
			for(j=0;j<tmpNumObj;j++) 
				{
				forcingInitialize(tmpDS,j,tmpt);
				et2Calc(tmpDS,j,tmpNV);
				varToPrint[j]=varToPrint[j]+tmpDS->EleET[j][2];
				}
                        break;
		case 8:
			for(j=0;j<tmpNumObj;j++) 
				{
                                varToPrint[j]=varToPrint[j]+tmpDS->ElePrep[j];
				}
                        break;
		case 9:
			for(j=0;j<tmpNumObj;j++) 
				{
				RechFluxCalc(tmpDS,j,tmpNV);
				varToPrint[j]=varToPrint[j]+tmpDS->Recharge[j];
				}
                        break;
	#ifdef LAYER3
		case 10:	
			for(j=0;j<tmpNumObj;j++) 
				{
/*			if(j==70)
				{
					printf("test");
				}
*/				varToPrint[j]=varToPrint[j]+NV_Ith_S(tmpNV,j+3*tmpDS->NumEle);
				}
                        break;
		case 11:
			for(j=0;j<tmpNumObj;j++) 
				{
				RechIFluxCalc(tmpDS,j,tmpNV);
				varToPrint[j]=varToPrint[j]+tmpDS->RechargeI[j];
				}
                        break;
	#endif
#endif
		case 12:
			tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+NV_Ith_S(tmpNV,j+tmpDS->totele);
				}
                        break;
	#ifdef SUB_SURF_RIV
		case 13:
			tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				varToPrint[j]=varToPrint[j]+NV_Ith_S(tmpNV,j+tmpDS->totele+tmpDS->NumRiv);
				}
                        break;
	#endif
                case 14:
		case 15:
			tmpNumObj=tmpDS->NumRiv;
	               	for(j=0;j<tmpNumObj;j++)
        		       {
		               tmpDS->FluxRiv[j][0]=0;
              		       }	
			for(j=0;j<tmpNumObj;j++) 
				{
				Riv_Outflow_Inflow(tmpDS,j,tmpt,tmpNV);
				varToPrint[j]=varToPrint[j]+tmpDS->FluxRiv[j][fileCounter+tmpCounter-14];
				}
                        break;
		case 16:
			tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				Riv_Left_Surf_Oflow(tmpDS,j,tmpNV);
				varToPrint[j]=varToPrint[j]+tmpDS->FluxRiv[j][fileCounter+tmpCounter-14];
				}
                        break;
		case 17:
			tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				Riv_Right_Surf_Oflow(tmpDS,j,tmpNV);
				varToPrint[j]=varToPrint[j]+tmpDS->FluxRiv[j][fileCounter+tmpCounter-14];
				}
                        break;
	#ifdef SUB_SURF_RIV
		case 18:
			tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				Riv_Left_Surf_Subflow(tmpDS,j,tmpNV);
				varToPrint[j]=varToPrint[j]+tmpDS->FluxRiv[j][fileCounter+tmpCounter-14];
				}
                        break;
		case 19:
			tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				Riv_Right_Surf_Subflow(tmpDS,j,tmpNV);
				varToPrint[j]=varToPrint[j]+tmpDS->FluxRiv[j][fileCounter+tmpCounter-14];
				}
                        break;
		case 20:
			tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				Riv_Bot_flow(tmpDS,j,tmpNV);
				varToPrint[j]=varToPrint[j]+tmpDS->FluxRiv[j][fileCounter+tmpCounter-14];
				}
                        break;
		case 21:
			tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				Riv_Left_Sub_Subflow(tmpDS,j,tmpNV);
				varToPrint[j]=varToPrint[j]+tmpDS->FluxRiv[j][fileCounter+tmpCounter-14];
				}
                        break;
		case 22:
			tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				Riv_Right_Sub_Subflow(tmpDS,j,tmpNV);
				varToPrint[j]=varToPrint[j]+tmpDS->FluxRiv[j][fileCounter+tmpCounter-14];
				}
                        break;
		case 23:
		case 24:
		       tmpNumObj=tmpDS->NumRiv;
                       for(j=0;j<tmpNumObj;j++)
                               {
                       tmpDS->FluxRiv[j][10]=0;
                               }

                 	tmpNumObj=tmpDS->NumRiv;
			for(j=0;j<tmpNumObj;j++) 
				{
				Riv_Sub_Outflow_Inflow(tmpDS,j,tmpNV);
				varToPrint[j]=varToPrint[j]+tmpDS->FluxRiv[j][fileCounter+tmpCounter-14];
				}
                        break;
	#endif
                default:
                        break;
        	}	
       	if(((int)tmpt%tmpIntv)==0)      
               	{
                fprintf(fpin[fileCounter],"%lf\t",tmpt);
       	        for(j=0;j<tmpNumObj;j++)
               	        {               
                       	fprintf(fpin[fileCounter],"%lf\t",varToPrint[j]/(tmpIntv/modIntv));
                        varToPrint[j]=0; 
       	                }
               	fprintf(fpin[fileCounter],"\n");     
                fflush(fpin[fileCounter]);           
       	        } 
	} 
/* print individual states */
void PrintData(FILE **outp,Control_Data *cD, Model_Data DS, N_Vector CV_Y, realtype t)
	{
	int k,m=0;
	for(k=0;k<cD->NumFilesToPrint;k++)
		{
/*	#ifdef LAYER2
		if(cD->FileNoToPrint[k]>9)
			{
			m=2;
			}
	#endif
*/		avgResults_NV(cD->FileNoToPrint[k],outp,DS->PrintVar[cD->FileNoToPrint[k]],CV_Y,DS,cD->fileInt[cD->FileNoToPrint[k]],t,m,cD->b);
		}
	}
  
