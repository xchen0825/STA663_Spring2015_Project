/*********************************************************************************
 * File        : initialize.c                                                    *
 * Function    : initialization of elemental attributes using relational database*
 * Version     : Nov, 2007 (2.0)                                                 *
 *-------------------------------------------------------------------------------*
 *                                                                               *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0..............................*
 * a) Correction of edge calculation term                                        *
 * b) Initialization of new variables for new process, shape representations  and*
 *    calibration (e.g. ET, Infiltration, Macropore, Stormflow, Element beneath  *
 *    river, river shapes, river bed property, thresholds for root zone,         *
 *    infiltration and macropore depths, land cover attributes etc)              *
 *--------------------------------------------------------------------------------*
 * For questions or comments, please contact                                      *
 *      --> Mukesh Kumar (mk176@duke.edu)                                         *
 * This code is free for research purpose only.                                   *
 * Please provide relevant references if you use this code in your research work  *
 *--------------------------------------------------------------------------------*
 *********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "sundialstypes.h"
#include "nvector_serial.h"
#include "pihm.h"  

void updateZmax(void *mData){
	int i, j;
	Model_Data MD;
	MD = (Model_Data) mData;
	for(i=0; i<MD->NumEle; i++){
		MD->Ele[i].zmax = (1.0/3.0)*(MD->Node[MD->Ele[i].node[0] - 1].zmax + MD->Node[MD->Ele[i].node[1] - 1].zmax + MD->Node[MD->Ele[i].node[2] - 1].zmax);
	}
    
	for(i=0; i<MD->NumRiv; i++){
		MD->Riv[i].zmax = (1.0/2.0)*(MD->Node[MD->Riv[i].FromNode - 1].zmax + MD->Node[MD->Riv[i].ToNode - 1].zmax);
	}
    
}

void updateZmin(void *mData){
	int i, j;
	Model_Data MD;
	MD = (Model_Data) mData;
	for(i=0; i<MD->NumNode; i++)
		MD->Node[i].zmin = MD->Node[i].zmax - MD->Node[i].delZ;
	for(i=0; i<MD->NumEle; i++){
		MD->Ele[i].zmin = (1.0/3.0)*(MD->Node[MD->Ele[i].node[0] - 1].zmin + MD->Node[MD->Ele[i].node[1] - 1].zmin + MD->Node[MD->Ele[i].node[2] - 1].zmin);
	}
    
	for(i=0; i<MD->NumRiv; i++){
		MD->Riv[i].zmin = MD->Riv[i].zmax - MD->Riv[i].depth; // (1.0/2.0)*(MD->Node[MD->Riv[i].FromNode - 1].zmin + MD->Node[MD->Riv[i].ToNode - 1].zmin);
	}
	
	for(i=0; i<MD->NumRiv; i++){
		MD->Ele[i + MD->NumEle].zmax = MD->Riv[i].zmin;
		if(MD->Riv[i].LeftEle > 0 && MD->Riv[i].RightEle > 0) MD->Ele[i + MD->NumEle].zmin = MD->Riv[i].zmax - (0.5 * (MD->Ele[MD->Riv[i].LeftEle - 1].zmax + MD->Ele[MD->Riv[i].RightEle - 1].zmax) - 0.5 * (MD->Ele[MD->Riv[i].LeftEle - 1].zmin + MD->Ele[MD->Riv[i].RightEle - 1].zmin));
		if(MD->Riv[i].LeftEle == 0 ) if(MD->Riv[i].RightEle > 0) MD->Ele[i + MD->NumEle].zmin = MD->Riv[i].zmax - (MD->Ele[MD->Riv[i].RightEle - 1].zmax -  MD->Ele[MD->Riv[i].RightEle - 1].zmin); else MD->Ele[i + MD->NumEle].zmin = MD->Riv[i].zmax;
		if(MD->Riv[i].RightEle == 0) if(MD->Riv[i].LeftEle > 0 ) MD->Ele[i + MD->NumEle].zmin = MD->Riv[i].zmax - (MD->Ele[MD->Riv[i].LeftEle - 1].zmax  -  MD->Ele[MD->Riv[i].LeftEle - 1].zmin ); else MD->Ele[i + MD->NumEle].zmin = MD->Riv[i].zmax;
	}
}


void calcDir(void *mData){
	int i, j, k;
	Model_Data MD;
	MD = (Model_Data) mData;
	printf("\ncalcDir...");
	for(i=0; i<MD->NumEle; i++){
		for(j=0; j<3; j++){
			MD->Ele[i].eleFlowDir[j] = 0;
            if(MD->Ele[i].BC[j] > -4){
                if(MD->Ele[i].nabr[j] > 0){
                    if(MD->Ele[i].zmax - MD->Ele[MD->Ele[i].nabr[j]-1].zmax >= -0.000001)
						MD->Ele[i].eleFlowDir[j] = 1;
                }
            }
            else
                if(MD->Ele[i].zmax - MD->Riv[-(MD->Ele[i].BC[j]/4)-1].zmax >= -0.000001)
                    MD->Ele[i].eleFlowDir[j] = 1;
		}
	}
	printf("Done!\n");
}


void idSinks(void *mData, int flag){
	int i, j, k;
	int updated;
	int out, sinks, not_a_sink;
	Model_Data MD;
    MD = (Model_Data) mData;
	printf("idSinks...");
	for(i=0; i<MD->NumEle; i++){
		MD->Ele[i].SINK = 0;
		if(MD->Ele[i].eleFlowDir[0] + MD->Ele[i].eleFlowDir[1] + MD->Ele[i].eleFlowDir[2] < 1)
			MD->Ele[i].SINK = 1;
	}
    
    
	sinks = 0;
	for(i=0; i<MD->NumRiv; i++){
        if(MD->Riv[i].LeftEle > 0) if(MD->Ele[MD->Riv[i].LeftEle-1].SINK == 1)
            if(flag == 0) printf("\nSTR %d LEle %d ", i+1, MD->Riv[i].LeftEle);
        if(MD->Riv[i].RightEle > 0) if(MD->Ele[MD->Riv[i].RightEle-1].SINK == 1)
            if(flag == 0) printf("\nSTR %d REle %d ", i+1, MD->Riv[i].RightEle);
        if(MD->Riv[i].LeftEle > 0 && MD->Riv[i].RightEle > 0) if(MD->Ele[MD->Riv[i].LeftEle-1].SINK == 1 || MD->Ele[MD->Riv[i].RightEle-1].SINK == 1)
            sinks++;
	}
	if(flag == 0) printf("\nStreams with Sinks = %d ", sinks);
	
	if(flag == 1){
        updated = 1;
        while(updated == 1){
            updated = 0;
            for(i=0; i<MD->NumEle; i++){
                //IF ONLY OUTLETS ARE TO THE SINKS => ELEMENT SHOULD BE CONSIDERED SINK STILL !!
                sinks = 0; out = 0; not_a_sink = 0;
                for(j=0; j<3; j++){
                    if(MD->Ele[i].nabr[j] > 0 && MD->Ele[i].BC[j] > -4){
						sinks += MD->Ele[MD->Ele[i].nabr[j] - 1].SINK;
						out   += MD->Ele[i].eleFlowDir[j];
                    }
                    
                    if(MD->Ele[i].BC[j] <= -4 && MD->Ele[i].eleFlowDir[j] == 1)
                        not_a_sink = 1;
                }
                if(not_a_sink == 0 && out > 0 && sinks > 0 && sinks == out && MD->Ele[i].SINK == 0){
                    MD->Ele[i].SINK = 1;
                    updated = 1;
                }
            }
        }
	}
	printf("Done!\n");
    
}


void updateZmaxRiv(void *mData){
	int i, j;
	Model_Data MD;
    MD = (Model_Data) mData;
    
	for(i=0; i<MD->NumRiv; i++){
		MD->Riv[i].zmax = (1.0/2.0)*(MD->Node[MD->Riv[i].FromNode - 1].zmax + MD->Node[MD->Riv[i].ToNode - 1].zmax);
		if(MD->Riv[i].LeftEle > 0) MD->Ele[MD->Riv[i].LeftEle-1].zmax = (1.0/3.0)*(MD->Node[MD->Ele[MD->Riv[i].LeftEle-1].node[0] - 1].zmax + MD->Node[MD->Ele[MD->Riv[i].LeftEle-1].node[1] - 1].zmax + MD->Node[MD->Ele[MD->Riv[i].LeftEle-1].node[2] - 1].zmax);
		if(MD->Riv[i].RightEle > 0) MD->Ele[MD->Riv[i].RightEle-1].zmax = (1.0/3.0)*(MD->Node[MD->Ele[MD->Riv[i].RightEle-1].node[0] - 1].zmax + MD->Node[MD->Ele[MD->Riv[i].RightEle-1].node[1] - 1].zmax + MD->Node[MD->Ele[MD->Riv[i].RightEle-1].node[2] - 1].zmax);
	}
    
}

void calcDirRiv(void *mData){
	int i, j, k;
	Model_Data MD;
    MD = (Model_Data) mData;
	printf("\ncalcDirRiv...");
	for(i=0; i<MD->NumRiv; i++){
		for(j=0; j<3; j++){
		    if(MD->Riv[i].LeftEle > 0){
                MD->Ele[MD->Riv[i].LeftEle-1].eleFlowDir[j] = 0;
                if(MD->Ele[MD->Riv[i].LeftEle-1].BC[j]>-4){
                    if(MD->Ele[MD->Riv[i].LeftEle-1].nabr[j] > 0)
                        if(MD->Ele[MD->Riv[i].LeftEle-1].zmax - MD->Ele[MD->Ele[MD->Riv[i].LeftEle-1].nabr[j]-1].zmax >= -0.000001)
                            MD->Ele[MD->Riv[i].LeftEle-1].eleFlowDir[j] = 1;
                }
                else
                    if(MD->Ele[MD->Riv[i].LeftEle-1].zmax - MD->Riv[-(MD->Ele[MD->Riv[i].LeftEle-1].BC[j]/4)-1].zmax >= -0.000001)
						MD->Ele[MD->Riv[i].LeftEle-1].eleFlowDir[j] = 1;
		    }
		    if(MD->Riv[i].RightEle > 0){
                MD->Ele[MD->Riv[i].RightEle-1].eleFlowDir[j] = 0;
                if(MD->Ele[MD->Riv[i].RightEle-1].BC[j]>-4){
                    if(MD->Ele[MD->Riv[i].RightEle-1].nabr[j] > 0)
                        if(MD->Ele[MD->Riv[i].RightEle-1].zmax - MD->Ele[MD->Ele[MD->Riv[i].RightEle-1].nabr[j]-1].zmax >= -0.000001)
                            MD->Ele[MD->Riv[i].RightEle-1].eleFlowDir[j] = 1;
                }
                else
                    if(MD->Ele[MD->Riv[i].RightEle-1].zmax - MD->Riv[-(MD->Ele[MD->Riv[i].RightEle-1].BC[j]/4)-1].zmax >= -0.000001)
						MD->Ele[MD->Riv[i].RightEle-1].eleFlowDir[j] = 1;
                
            }
            
		}
	}

	printf("Done!\n");
}


realtype sinkdepth(void *mData, int i){
	int j;
	realtype mindepth = 99999;
	realtype tmp = -999999.99;
	Model_Data MD;
    MD = (Model_Data) mData;
	for(j=0; j<3; j++){
		if(MD->Ele[i].nabr[j] == 0)
			continue;
		else if(MD->Ele[i].BC[j] <= -4)
			tmp = MD->Riv[-(MD->Ele[i].BC[j]/4)-1].zmax - MD->Ele[i].zmax > tmp ? MD->Riv[-(MD->Ele[i].BC[j]/4)-1].zmax - MD->Ele[i].zmax : tmp;
		else
			tmp = MD->Ele[MD->Ele[i].nabr[j]-1].zmax - MD->Ele[i].zmax > tmp ? MD->Ele[MD->Ele[i].nabr[j]-1].zmax - MD->Ele[i].zmax : tmp;
	}
	return tmp;
}

void initialize(char *filename, Model_Data DS, Control_Data *CS, N_Vector CV_Y)
	{
  	int i,j,k,tmpBool,BoolBR,BoolR=0;
  	Model_Data MD;
        MD = (Model_Data) DS;
	int *eleFlag, updateFlag, tmp, tmpE, tmpR;
	int updateSink, J, K, Ktmp, before, after;
	realtype del;
	FILE *newelemesh; char newelemeshfilename[200];
	realtype a_x, a_y, b_x, b_y, c_x, c_y,distX,distY;
  	realtype a_zmin, a_zmax, b_zmin, b_zmax, c_zmin, c_zmax; 
  	realtype tempvalue1,tempvalue2,tempvalue3,tempvalue4;
  	FILE *init_file;
  	char *fn;
  	realtype *zmin_cor;
	realtype locMult,locKMult;
  	zmin_cor=(realtype *)malloc(DS->NumEle*sizeof(realtype));
  
  	printf("\nInitializing data structure ... ");

  	/* allocate memory storage to flux terms */
  	DS->FluxSurf = (realtype **)malloc(DS->NumEle*sizeof(realtype));
  	DS->FluxSub = (realtype **)malloc(DS->NumEle*sizeof(realtype));
  	DS->FluxRiv = (realtype **)malloc(DS->NumRiv*sizeof(realtype));
  	DS->EleET = (realtype **)malloc(DS->NumEle*sizeof(realtype));
  	DS->ElePrep = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleViR = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->Recharge = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->RechargeI = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleIS = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleISmax = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleISsnowmax = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleSnow = (realtype *)malloc(DS->NumEle*sizeof(realtype));  
  	DS->EleSnowGrnd = (realtype *)malloc(DS->NumEle*sizeof(realtype));  
  	DS->EleSnowCanopy = (realtype *)malloc(DS->NumEle*sizeof(realtype));  
  	DS->EleTF = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleETloss = (realtype *)malloc(DS->NumEle*sizeof(realtype));
  	DS->EleNetPrep = (realtype *)malloc(DS->NumEle*sizeof(realtype));
	DS->FluxSource = (realtype *)malloc(DS->NumEle*sizeof(realtype));               //xchen_20141006
//        DS->Calc = (int *)malloc(DS->NumEle*sizeof(int));
	 for(i=0;i<DS->NumEle;i++)
		{
    		DS->FluxSurf[i] = (realtype *)malloc(3*sizeof(realtype));
    		DS->FluxSub[i] = (realtype *)malloc(3*sizeof(realtype));
    		DS->EleET[i] = (realtype *)malloc(4*sizeof(realtype));
//                DS->Calc[i] = (int *)malloc(3*sizeof(int));
    		a_x = DS->Node[DS->Ele[i].node[0]-1].x;
    		b_x = DS->Node[DS->Ele[i].node[1]-1].x;
    		c_x = DS->Node[DS->Ele[i].node[2]-1].x;
    		a_y = DS->Node[DS->Ele[i].node[0]-1].y;
    		b_y = DS->Node[DS->Ele[i].node[1]-1].y;
    		c_y = DS->Node[DS->Ele[i].node[2]-1].y;

    		a_zmin = DS->Node[DS->Ele[i].node[0]-1].zmin;
    		b_zmin = DS->Node[DS->Ele[i].node[1]-1].zmin;
    		c_zmin = DS->Node[DS->Ele[i].node[2]-1].zmin;
    		a_zmax = DS->Node[DS->Ele[i].node[0]-1].zmax;
    		b_zmax = DS->Node[DS->Ele[i].node[1]-1].zmax;
    		c_zmax = DS->Node[DS->Ele[i].node[2]-1].zmax;
    
    		DS->Ele[i].area = 0.5*fabs((b_x - a_x)*(c_y - a_y) - (b_y - a_y)*(c_x - a_x));
	//	printf("\n%lf,",DS->Ele[i].area);
         	DS->Ele[i].zmax = (a_zmax + b_zmax + c_zmax)/3.0;
    		DS->Ele[i].zmin = (a_zmin + b_zmin + c_zmin)/3.0; 
    		DS->Ele[i].edge[0] = pow((b_x - c_x), 2) + pow((b_y - c_y), 2);
    		DS->Ele[i].edge[1] = pow((c_x - a_x), 2) + pow((c_y - a_y), 2);
    		DS->Ele[i].edge[2] = pow((a_x - b_x), 2) + pow((a_y - b_y), 2);

    		/* calculate centroid of triangle */
    		DS->Ele[i].x = (a_x + b_x + c_x)/3.0;
    		DS->Ele[i].y = (a_y + b_y + c_y)/3.0;
    
    
    		/* calculate circumcenter of triangle */
  		/*  DS->Ele[i].x = a_x - ((b_y - a_y)*DS->Ele[i].edge[2] - (c_y - a_y)*DS->Ele[i].edge[0])/(4*DS->Ele[i].area);
    		DS->Ele[i].y = a_y + ((b_x - a_x)*DS->Ele[i].edge[2] - (c_x - a_x)*DS->Ele[i].edge[0])/(4*DS->Ele[i].area); 
    		*/
    		DS->Ele[i].edge[0] = sqrt(DS->Ele[i].edge[0]);
    		DS->Ele[i].edge[1] = sqrt(DS->Ele[i].edge[1]);
    		DS->Ele[i].edge[2] = sqrt(DS->Ele[i].edge[2]);
    		DS->Ele[i].KsatV = CS->Cal.KsatV*DS->Geol[(DS->Ele[i].geol-1)].KsatV;
    		DS->Ele[i].infKsatV = CS->Cal.infKsatV*DS->Soil[(DS->Ele[i].soil-1)].KsatV;
		DS->Ele[i].macD=CS->Cal.macD*DS->Geol[DS->Ele[i].geol-1].macD;
		/* Note above porosity statement should be replaced by geologic porosity (in comments below) if the data is available */
		DS->Ele[i].Porosity = CS->Cal.Porosity*(DS->Soil[(DS->Ele[i].soil-1)].ThetaS - DS->Soil[(DS->Ele[i].soil-1)].ThetaR);
//		DS->Ele[i].Porosity = CS->Cal.Porosity*(DS->Geol[(DS->Ele[i].geol-1)].ThetaS - DS->Geol[(DS->Ele[i].geol-1)].ThetaR);
 		if((DS->Ele[i].Porosity>1)&&(DS->Ele[i].Porosity<=0))
			{
			printf("Warning: Porosity value out of bounds");
			getchar();
			} 
		DS->Ele[i].Alpha = CS->Cal.Alpha*DS->Soil[(DS->Ele[i].soil-1)].Alpha;
    		DS->Ele[i].Beta = CS->Cal.Beta*DS->Soil[(DS->Ele[i].soil-1)].Beta; 
		/* Note above van genuchten statement should be replaced by geologic parameters (in comments below) if the data is available */
//		DS->Ele[i].Alpha = CS->Cal.Alpha*DS->Geol[(DS->Ele[i].geol-1)].Alpha;
//    		DS->Ele[i].Beta = CS->Cal.Beta*DS->Geol[(DS->Ele[i].geol-1)].Beta; 
    		DS->Ele[i].hAreaF = CS->Cal.hAreaF*DS->Soil[(DS->Ele[i].soil-1)].hAreaF; 
    		DS->Ele[i].vAreaF = CS->Cal.vAreaF*DS->Geol[(DS->Ele[i].geol-1)].vAreaF; 
		DS->Ele[i].KsatH = CS->Cal.KsatH*DS->Geol[(DS->Ele[i].geol-1)].KsatH;
//		DS->Ele[i].macKsatH = CS->Cal.KsatH*DS->Geol[(DS->Ele[i].geol-1)].KsatH;//Xing_03312013_remove the effect of macropore line 124 and 126
  		DS->Ele[i].macKsatV = CS->Cal.macKsatV*DS->Soil[(DS->Ele[i].soil-1)].macKsatV; 
  //  		DS->Ele[i].macKsatV = CS->Cal.KsatV*DS->Geol[(DS->Ele[i].geol-1)].KsatV;
    		DS->Ele[i].macKsatH = CS->Cal.macKsatH*DS->Geol[(DS->Ele[i].geol-1)].macKsatH;
		DS->Ele[i].macD=CS->Cal.macD*DS->Geol[DS->Ele[i].geol-1].macD; 
    		DS->Ele[i].infD=CS->Cal.infD*DS->Soil[DS->Ele[i].soil-1].infD;
    
    		DS->Ele[i].RzD=CS->Cal.RzD*DS->LandC[DS->Ele[i].LC-1].RzD;  
    		DS->Ele[i].LAImax = DS->LandC[DS->Ele[i].LC-1].LAImax;
    		DS->Ele[i].Rmin = DS->LandC[DS->Ele[i].LC-1].Rmin;
    		DS->Ele[i].Rs_ref = DS->LandC[DS->Ele[i].LC-1].Rs_ref;
    		DS->Ele[i].Albedo = CS->Cal.Albedo*DS->LandC[DS->Ele[i].LC-1].Albedo;
		if(DS->Ele[i].Albedo>1)
			{
                        printf("Warning: Albedo out of bounds");
                        getchar();
			}
    		DS->Ele[i].VegFrac = CS->Cal.VegFrac*DS->LandC[DS->Ele[i].LC-1].VegFrac;                
    		DS->Ele[i].Rough = CS->Cal.Rough*DS->LandC[DS->Ele[i].LC-1].Rough;

		DS->Ele[i].windH=DS->windH[DS->Ele[i].WindVel-1];
		}
	//fflush(stdout);
	//getchar();
        
  	for(i=0; i<DS->NumRiv; i++)
  		{
    		DS->FluxRiv[i] = (realtype *)malloc(11*sizeof(realtype));
		for(j=0;j<3;j++)
			{
			/* Note: Strategy to use BC < -4 for river identification */
			if(DS->Ele[DS->Riv[i].LeftEle-1].nabr[j]==DS->Riv[i].RightEle)
				{
				DS->Ele[DS->Riv[i].LeftEle-1].BC[j]=-4*(i+1);
				}
			if(DS->Ele[DS->Riv[i].RightEle-1].nabr[j]==DS->Riv[i].LeftEle)
				{
				DS->Ele[DS->Riv[i].RightEle-1].BC[j]=-4*(i+1);
				}
			}
    		DS->Riv[i].x = (DS->Node[DS->Riv[i].FromNode-1].x + DS->Node[DS->Riv[i].ToNode-1].x)/2;
    		DS->Riv[i].y = (DS->Node[DS->Riv[i].FromNode-1].y + DS->Node[DS->Riv[i].ToNode-1].y)/2;
    		DS->Riv[i].zmax = (DS->Node[DS->Riv[i].FromNode-1].zmax + DS->Node[DS->Riv[i].ToNode-1].zmax)/2;
    		DS->Riv[i].depth = CS->Cal.rivDepth*DS->Riv_Shape[DS->Riv[i].shape-1].depth;
		DS->Riv[i].coeff=CS->Cal.rivShapeCoeff*DS->Riv_Shape[DS->Riv[i].shape - 1].coeff;
    		DS->Riv[i].zmin = DS->Riv[i].zmax - DS->Riv[i].depth;  
    		DS->Riv[i].Length = sqrt(pow(DS->Node[DS->Riv[i].FromNode-1].x -DS->Node[DS->Riv[i].ToNode-1].x, 2) + pow(DS->Node[DS->Riv[i].FromNode-1].y - DS->Node[DS->Riv[i].ToNode-1].y, 2));
		DS->Riv[i].KsatH=CS->Cal.rivKsatH*DS->Riv_Mat[DS->Riv[i].material-1].KsatH;
		DS->Riv[i].KsatV=CS->Cal.rivKsatV*DS->Riv_Mat[DS->Riv[i].material-1].KsatV;
		DS->Riv[i].bedThick=CS->Cal.rivbedThick*DS->Riv_Mat[DS->Riv[i].material-1].bedThick;
		DS->Riv[i].Rough=CS->Cal.rivRough*DS->Riv_Mat[DS->Riv[i].material - 1].Rough;
		/* Initialization for rectangular cells beneath river */
		/* Note: Ideally this data should be read from the decomposition itself */
		/* but it is not supported right now in PIHMgis (Bhatt, G and Kumar, M; 2007) */
		DS->Ele[i+DS->NumEle].zmax=DS->Riv[i].zmin;
		DS->Ele[i+DS->NumEle].zmin=DS->Riv[i].zmax-(0.5*(DS->Ele[DS->Riv[i].LeftEle-1].zmax+DS->Ele[DS->Riv[i].RightEle-1].zmax)-0.5*(DS->Ele[DS->Riv[i].LeftEle-1].zmin+DS->Ele[DS->Riv[i].RightEle-1].zmin));
//		DS->Ele[i+DS->NumEle].zmin=DS->Riv[i].zmax-40;
		DS->Ele[i+DS->NumEle].macD=0.5*(DS->Ele[DS->Riv[i].LeftEle-1].macD+DS->Ele[DS->Riv[i].RightEle-1].macD)>DS->Riv[i].depth?0.5*(DS->Ele[DS->Riv[i].LeftEle-1].macD+DS->Ele[DS->Riv[i].RightEle-1].macD)-DS->Riv[i].depth:0;	
		DS->Ele[i+DS->NumEle].macKsatH=0.5*(DS->Ele[DS->Riv[i].LeftEle-1].macKsatH+DS->Ele[DS->Riv[i].RightEle-1].macKsatH);
		DS->Ele[i+DS->NumEle].vAreaF=0.5*(DS->Ele[DS->Riv[i].LeftEle-1].vAreaF+DS->Ele[DS->Riv[i].RightEle-1].vAreaF);
		DS->Ele[i+DS->NumEle].KsatH=0.5*(DS->Ele[DS->Riv[i].LeftEle-1].KsatH+DS->Ele[DS->Riv[i].RightEle-1].KsatH);
		DS->Ele[i+DS->NumEle].Porosity=0.5*(DS->Ele[DS->Riv[i].LeftEle-1].Porosity+DS->Ele[DS->Riv[i].RightEle-1].Porosity);
 		}
    for(i=0;i<DS->NumRiv;i++)
        {
        if(DS->Riv[i].down > 0)
           {
           DS->Riv[i].Dist[0] = 0.5*(DS->Riv[i].Length+DS->Riv[DS->Riv[i].down - 1].Length);
           }
        else
           {
           DS->Riv[i].Dist[0] = sqrt(pow(DS->Riv[i].x - DS->Node[DS->Riv[i].ToNode-1].x, 2) + pow(DS->Riv[i].y - DS->Node[DS->Riv[i].ToNode-1].y, 2));
           }
        if(DS->Riv[i].LeftEle > 0)
           {
           DS->Riv[i].Dist[1] = sqrt(pow((DS->Riv[i].x - DS->Ele[DS->Riv[i].LeftEle - 1].x), 2) + pow((DS->Riv[i].y - DS->Ele[DS->Riv[i].LeftEle - 1].y), 2));
           for(j=0; j < 3; j++)
      		   {
       		   if(DS->Ele[DS->Riv[i].LeftEle - 1].nabr[j] == DS->Riv[i].RightEle)
        				{
                        DS->Riv[i].lrEdge[0]=j;
          				break;
        				}
               }
            }
        if(DS->Riv[i].RightEle > 0)
           {
           DS->Riv[i].Dist[2] = sqrt(pow((DS->Riv[i].x - DS->Ele[DS->Riv[i].RightEle - 1].x), 2) + pow((DS->Riv[i].y - DS->Ele[DS->Riv[i].RightEle - 1].y), 2));
           for(j=0; j < 3; j++)
      		   {
               if(DS->Ele[DS->Riv[i].RightEle - 1].nabr[j] == DS->Riv[i].LeftEle)
        				{
                        DS->Riv[i].lrEdge[1]=j;
          				break;
        				}
               }
           }
        }
	for(i=0;i<DS->NumPrep;i++)
		{
		for(j=0; j<DS->TSD_Prep[i].length; j++)
			{
				DS->TSD_Prep[i].TS[j][1]=CS->Cal.Prep*DS->TSD_Prep[i].TS[j][1];
			}
		}
        for(i=0;i<DS->NumTemp;i++)
                {
                for(j=0; j<DS->TSD_Temp[i].length; j++)
                        {
                        DS->TSD_Temp[i].TS[j][1]=CS->Cal.Temp*DS->TSD_Temp[i].TS[j][1];
                        }
                }
	CS->FileNoToPrint=(int *)malloc(sizeof(int)*CS->totFiles);
	DS->PrintVar=(realtype **)malloc(sizeof(realtype **)*CS->totFiles);
	CS->NumFilesToPrint=0;
        for(i=0;i<CS->totFiles;i++)
                {
		if(CS->fileID[i]==1)
			{
			CS->FileNoToPrint[CS->NumFilesToPrint++]=i;
			}
		if(i<CS->totFiles-13)
			{
                        DS->PrintVar[i]=(realtype *)calloc(DS->NumEle,sizeof(realtype));
			}
		else
			{
                        DS->PrintVar[i]=(realtype *)calloc(DS->NumRiv,sizeof(realtype));
			}
                }
	/* Debugging artifacts in data created due to coarser resolution of model elements */
        if (CS->Debug == 1) {
            
            printf("\n ************************************************************");
            printf("\n **************** RESOLVE TOPOGRAPHIC ERROR ***************** \n");
            printf(" ************************************************************\n");
            
            
            for(i=0; i<MD->NumNode; i++) //delZ stores active sub-surface layer thickness
                MD->Node[i].delZ = MD->Node[i].zmax - MD->Node[i].zmin;
            
            // CHECK IF INDIVIDUAL STREAM SEG HAS CORRECT TOPOGRAPHY
            printf("\nChecking Individual Stream Topology... \n");
            tmp = 0;
            for(i=0; i<MD->NumRiv; i++){
                printf("\nRiv\t%d\tFromTo\t%lf", i+1, MD->Node[MD->Riv[i].FromNode-1].zmax - MD->Node[MD->Riv[i].ToNode-1].zmax);
                if(MD->Node[MD->Riv[i].FromNode-1].zmax - MD->Node[MD->Riv[i].ToNode-1].zmax < 0){
                    printf("\nRiv Seg # %d", i+1); tmp++;
                }
            }
            if(tmp > 0) { printf("\nAbove streams should be fixed\n"); getchar(); }
            else printf("\nDone!\n");
            //getchar();
            
            // CHECK IF UP-DOWN STREAM SEG HAS CORRECT TOPOGRAPHY
            printf("\nChecking Up-Down Stream Pair Topology... \n");
            tmp = 0;
            for(i=0; i<MD->NumRiv; i++)
                if(MD->Riv[i].down > 0){
                    printf("\nRiv\t%d\tSlope\t%lf", i+1, MD->Riv[i].zmax - MD->Riv[MD->Riv[i].down - 1].zmax);
                    if(MD->Riv[i].zmax - MD->Riv[MD->Riv[i].down - 1].zmax < 0){
                        printf("\n Riv Segs %d %d", i+1, MD->Riv[i].down); tmp++;
                    }
                } 
            if(tmp > 0) { printf("\nAbove up-down streams should be fixed\n"); getchar(); }
            else printf("\nDone!\n");
            //getchar();
            
            // FLAG NODES IF THEY ARE ON THE STREAMS :: FROM HERE ONWARMD ELEVATION OF THOSE NODES SHOULD NOT BE CHANGED
            for(i=0; i<MD->NumNode; i++)
                MD->Node[i].onRiver = 0;
            for(i=0; i<MD->NumRiv; i++){
                MD->Node[MD->Riv[i].FromNode-1].onRiver = 1;
                MD->Node[MD->Riv[i].ToNode  -1].onRiver = 1;
            }
            
            // Initialize every triangle as NO-sink
            for(i=0; i<MD->NumEle; i++)
                MD->Ele[i].SINK = 0;
            
            calcDir(MD);
            idSinks(MD, 0);
            printf("\nSinks Next to River Segs? ");
            for(i=0; i<MD->NumRiv; i++){
                if(MD->Riv[i].LeftEle > 0) if(MD->Ele[MD->Riv[i].LeftEle-1].SINK == 1)
                    printf(" #L %d ", i+1);
                if(MD->Riv[i].RightEle > 0) if(MD->Ele[MD->Riv[i].RightEle-1].SINK == 1)
                    printf(" #R %d ", i+1);
            }
            //getchar();
            
            calcDir(MD);
            idSinks(MD, 0);
            
            updateSink = 1;
            while(updateSink == 1){
                updateSink = 0; printf("\n\n\nTrying to Fix Sinks Sharing Edge with Stream Reach...\n"); //getchar();
                for(i=0; i<MD->NumRiv; i++){
                    //printf("%d ", i+1);
                    if(MD->Riv[i].LeftEle > 0) if(MD->Ele[MD->Riv[i].LeftEle-1].SINK == 1){
                        if(MD->Ele[MD->Riv[i].LeftEle - 1].zmax - MD->Riv[i].zmax < -0.000001){ //this is necessary to make sure that the element is sink with respect to this particular stream segment
                            //printf(" L %d ", i+1);
                            //Find the node that is not on the Riv AND
                            //update that with the average of the rest two nodes
                            for(j=0; j<3; j++){ // FIND THE NODE OF TRIANGLE THAT IS NOT ON THS STREAM SEGMENT
                                if(! (MD->Ele[MD->Riv[i].LeftEle-1].node[j] == MD->Riv[i].FromNode || MD->Ele[MD->Riv[i].LeftEle-1].node[j] == MD->Riv[i].ToNode) )
                                    break;
                            }
                            //printf("\n>%d %d %d", i, j, MD->Ele[MD->Riv[i].LeftEle-1].BC[j]); getchar();
                            if(MD->Node[MD->Ele[MD->Riv[i].LeftEle-1].node[j]-1].onRiver == 0){ //ENSURE THAT PREVIOUSLY SELECTED NODE IS NOT ON ANY OTHER RIVER SEG
                                
                                printf("\n\nRiv %d LeftEle %d Sink = %d >> ", i+1, MD->Riv[i].LeftEle, MD->Ele[MD->Riv[i].LeftEle-1].SINK);
                                
                                MD->Node[MD->Ele[MD->Riv[i].LeftEle-1].node[j]-1].zmax =  (0.5)*(MD->Node[MD->Riv[i].FromNode-1].zmax + MD->Node[MD->Riv[i].ToNode-1].zmax) + 0.000001;
                                updateSink = 1;
                                MD->Ele[MD->Riv[i].LeftEle-1].SINK = 0;
                                
                                updateZmaxRiv(MD);
                                calcDirRiv(MD); //
                                idSinks(MD, 0); //
                                printf(" >> Now LeftEle %d Sink = %d", i+1, MD->Ele[MD->Riv[i].LeftEle-1].SINK);
                            }
                        }
                    }
                    
                    if(MD->Riv[i].RightEle > 0) if(MD->Ele[MD->Riv[i].RightEle-1].SINK == 1){
                        //if(MD->Riv[i].zmax > MD->Ele[MD->Riv[i].RightEle - 1].zmax){
                        if(MD->Ele[MD->Riv[i].RightEle - 1].zmax - MD->Riv[i].zmax < -0.000001){
                            //printf(" L %d ", i+1);
                            //Find the node that is not on the Riv AND
                            //update that with the average of the rest two nodes
                            for(j=0; j<3; j++){ // FIND THE NODE OF TRIANGLE THAT IS NOT ON THS STREAM SEGMENT
                                if(! (MD->Ele[MD->Riv[i].RightEle-1].node[j] == MD->Riv[i].FromNode || MD->Ele[MD->Riv[i].RightEle-1].node[j] == MD->Riv[i].ToNode) )
                                    break;
                            }
                            if(MD->Node[MD->Ele[MD->Riv[i].RightEle-1].node[j]-1].onRiver == 0){ //ENSURE THAT PREVIOUSLY SELECTED NODE IS NOT ON ANY OTHER RIVER SEG
                                
                                printf("\n\nRiv %d RightEle %d Sink = %d >> ", i+1, MD->Riv[i].RightEle, MD->Ele[MD->Riv[i].RightEle-1].SINK);
                                MD->Node[MD->Ele[MD->Riv[i].RightEle-1].node[j]-1].zmax =  (0.5)*(MD->Node[MD->Riv[i].FromNode-1].zmax + MD->Node[MD->Riv[i].ToNode-1].zmax) + 0.000001;
                                updateSink = 1;
                                MD->Ele[MD->Riv[i].RightEle-1].SINK = 0;
                                
                                updateZmaxRiv(MD);
                                calcDirRiv(MD); //
                                idSinks(MD, 0); //
                                printf(" >> Now RightEle %d Sink = %d", i+1, MD->Ele[MD->Riv[i].LeftEle-1].SINK);
                            }
                        }
                    }
                }
            }
            printf("\n\nTrying to Fix Sinks Sharing Edge with Stream Reach... Done!\n"); //getchar();
            updateZmax(MD);
            
            
            printf("\nSinks Sharing Edge with River... Remaining? ");
            for(i=0; i<MD->NumRiv; i++){
                if(MD->Riv[i].LeftEle > 0) if(MD->Ele[MD->Riv[i].LeftEle-1].SINK == 1)
                    printf(" $L %d ", i+1);
                if(MD->Riv[i].RightEle > 0) if(MD->Ele[MD->Riv[i].RightEle-1].SINK == 1)
                    printf(" $R %d ", i+1);
            }
            //getchar();
            
            
            updateSink = 1;
            while(updateSink == 1){
                updateSink = 0;
                
                printf("\n\nTRYING TO FIX 'ELE' SINKS...\n");
                after = 0;
                for(i=0; i<MD->NumEle; i++)
                    if(MD->Ele[i].SINK){
                        printf(" %d ", i+1);
                        after++;
                    }
                printf("\nAbsolute Sinks = %d\n", after); //getchar();
                
                // Call Function to calculate outlets from a triangle
                calcDir(MD);
                // Call Function to identify sinks
                idSinks(MD, 1);
                
                after = 0;
                for(i=0; i<MD->NumEle; i++)
                    if(MD->Ele[i].SINK){
                        printf(" %d ", i+1);
                        after++;
                    }
                printf("\nRelative Sinks = %d\n", after); //getchar();
                
                for(i=0; i<MD->NumEle; i++){
                    if(MD->Ele[i].SINK==1){ // Check if Sink
                        del = 999999; J = 4; K = 4;
                        for(j=0; j<3; j++){ // Find 
                            if(MD->Ele[i].nabr[j] > 0 && MD->Ele[i].BC[j] > -4){
                                // Find which nabr (0 to 2) of ele J is I
                                for(k=0; k<3; k++){
                                    if(MD->Ele[MD->Ele[i].nabr[j]-1].nabr[k]-1 == i){
                                        Ktmp = k; break;
                                    }
                                }
                                if(MD->Ele[MD->Ele[i].nabr[j]-1].SINK != 1 && MD->Node[MD->Ele[i].node[j]-1].onRiver == 0){
                                    if(del > MD->Node[MD->Ele[MD->Ele[i].nabr[j]-1].node[Ktmp]-1].zmax - MD->Node[MD->Ele[i].node[j]-1].zmax){
                                        del = MD->Node[MD->Ele[MD->Ele[i].nabr[j]-1].node[Ktmp]-1].zmax - MD->Node[MD->Ele[i].node[j]-1].zmax;
                                        J = j;
                                        K = Ktmp;
                                    }
                                }
                            }
                        }
                        
                        if(J < 4  && K < 4) if(MD->Node[MD->Ele[i].node[J]-1].zmax != MD->Node[MD->Ele[MD->Ele[i].nabr[J]-1].node[K]-1].zmax){
                            MD->Node[MD->Ele[i].node[J]-1].zmax = MD->Node[MD->Ele[MD->Ele[i].nabr[J]-1].node[K]-1].zmax;
                            MD->Ele[i].SINK = 0;
                            updateSink = 1; //break;
                        }
                        
                    }
                }
                // Update Zmax of Triangles & Rivers
                updateZmax(MD);
            }
            
            
            
            printf("\n\nTrying to Fix Ele Sinks... DONE!\n\n");
            
            updateZmax(MD);
            updateZmin(MD);
            sprintf(newelemeshfilename, "%s.mesh.topo", filename);
            newelemesh=fopen(newelemeshfilename, "w");//newelemesh=fopen("mesh.new","w");
            fprintf(newelemesh,"%d\t%d\n", MD->NumEle, MD->NumNode);
            for (i=0; i<MD->NumEle; i++)
            { 
                fprintf(newelemesh, "%d\t", MD->Ele[i].index);
                fprintf(newelemesh, "%d\t%d\t%d\t", MD->Ele[i].node[0], MD->Ele[i].node[1], MD->Ele[i].node[2]);
                fprintf(newelemesh, "%d\t%d\t%d\n", MD->Ele[i].nabr[0], MD->Ele[i].nabr[1], MD->Ele[i].nabr[2]);
            }
            
        	// write in nodes information 
        	for (i=0; i<MD->NumNode; i++)
            {
                fprintf(newelemesh, "%d\t", MD->Node[i].index);
                fprintf(newelemesh, "%lf\t%lf\t", MD->Node[i].x, MD->Node[i].y);
                fprintf(newelemesh, "%lf\t%lf\n", MD->Node[i].zmin,MD->Node[i].zmax);
            }	
            fclose(newelemesh);
            
            
            
            
            for (i = 0; i < MD->NumRiv; i++) {
                if (MD->Riv[i].down > 0) {
                    if (MD->Riv[i].zmin < MD->Riv[MD->Riv[i].down - 1].zmin) {
                        BoolR = 1;
                        printf("\n Riv %d is lower than downstream Riv %d by %lf", i + 1, MD->Riv[i].down, MD->Riv[i].zmin - MD->Riv[MD->Riv[i].down - 1].zmin);
                    }
                }
            }
            if (BoolR == 1) {
                printf("\n\tRiver elevation correction needed");
                getchar();
            }
        }
  	for(i=0; i<DS->NumEle; i++)
  		{
                a_x = DS->Node[DS->Ele[i].node[0]-1].x;
                b_x = DS->Node[DS->Ele[i].node[1]-1].x;
                c_x = DS->Node[DS->Ele[i].node[2]-1].x;
                a_y = DS->Node[DS->Ele[i].node[0]-1].y;
                b_y = DS->Node[DS->Ele[i].node[1]-1].y;
                c_y = DS->Node[DS->Ele[i].node[2]-1].y;
		for(j=0;j<3;j++)
			{
			/* Note: Assumption here is that the forumulation is circumcenter based */
			switch(j)
				{
				case 0:
                			distX=(DS->Ele[i].x-0.5*(b_x+c_x));
                			distY=(DS->Ele[i].y-0.5*(b_y+c_y));
					break;
                                case 1:
                			distX=(DS->Ele[i].x-0.5*(c_x+a_x));
                			distY=(DS->Ele[i].y-0.5*(c_y+a_y));
					break;
                                case 2:
                			distX=(DS->Ele[i].x-0.5*(a_x+b_x));
                			distY=(DS->Ele[i].y-0.5*(a_y+b_y));
					break;
				}
			DS->Ele[i].surfH[j]=(DS->Ele[i].nabr[j]>0)?(DS->Ele[i].BC[j]>-4?(DS->Ele[DS->Ele[i].nabr[j]-1].zmax):DS->Riv[-(DS->Ele[i].BC[j]/4)-1].zmax):DS->Ele[i].BC[j]<=-4?DS->Riv[-(DS->Ele[i].BC[j]/4)-1].zmax:(DS->Ele[i].zmax); 
			DS->Ele[i].surfX[j]=(DS->Ele[i].nabr[j]>0)?(DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].x:DS->Riv[-(DS->Ele[i].BC[j]/4)-1].x):(DS->Ele[i].x-2*distX);          
			DS->Ele[i].surfY[j]=DS->Ele[i].nabr[j]>0?(DS->Ele[i].BC[j]>-4?DS->Ele[DS->Ele[i].nabr[j]-1].y:DS->Riv[-(DS->Ele[i].BC[j]/4)-1].y):(DS->Ele[i].y-2*distY);
                        DS->Ele[i].Calc[j]=0;
            if(DS->Ele[i].nabr[j] > 0)
               {
               DS->Ele[i].Dist[j] = sqrt(pow((DS->Ele[i].x - DS->Ele[DS->Ele[i].nabr[j] - 1].x), 2) + pow((DS->Ele[i].y - DS->Ele[DS->Ele[i].nabr[j] - 1].y), 2));
                       if(DS->Ele[i].index>DS->Ele[i].nabr[j])
                       {
                          for(k=0;k<3;k++)
                          if(DS->Ele[i].index==DS->Ele[DS->Ele[i].nabr[j]-1].nabr[k])
                          DS->Ele[i].Calc[j]=k+1;
                       }
               }       
            else
               {
               DS->Ele[i].Dist[j] = sqrt(pow(DS->Ele[i].edge[0]*DS->Ele[i].edge[1]*DS->Ele[i].edge[2]/(4*DS->Ele[i].area), 2) - pow(DS->Ele[i].edge[j]/2, 2));
               }
                	}
		DS->Ele[i].dhBYdx=-(DS->Ele[i].surfY[2]*(DS->Ele[i].surfH[1]-DS->Ele[i].surfH[0])+DS->Ele[i].surfY[1]*(DS->Ele[i].surfH[0]-DS->Ele[i].surfH[2])+DS->Ele[i].surfY[0]*(DS->Ele[i].surfH[2]-DS->Ele[i].surfH[1]))/(DS->Ele[i].surfX[2]*(DS->Ele[i].surfY[1]-DS->Ele[i].surfY[0])+DS->Ele[i].surfX[1]*(DS->Ele[i].surfY[0]-DS->Ele[i].surfY[2])+DS->Ele[i].surfX[0]*(DS->Ele[i].surfY[2]-DS->Ele[i].surfY[1]));  
		DS->Ele[i].dhBYdy=-(DS->Ele[i].surfX[2]*(DS->Ele[i].surfH[1]-DS->Ele[i].surfH[0])+DS->Ele[i].surfX[1]*(DS->Ele[i].surfH[0]-DS->Ele[i].surfH[2])+DS->Ele[i].surfX[0]*(DS->Ele[i].surfH[2]-DS->Ele[i].surfH[1]))/(DS->Ele[i].surfY[2]*(DS->Ele[i].surfX[1]-DS->Ele[i].surfX[0])+DS->Ele[i].surfY[1]*(DS->Ele[i].surfX[0]-DS->Ele[i].surfX[2])+DS->Ele[i].surfY[0]*(DS->Ele[i].surfX[2]-DS->Ele[i].surfX[1]));

  		}
  	/* initialize state variable */
  	/* relax case */
	if (CS->init_type == 0)
  		{
    		for(i=0; i<DS->NumEle; i++)
    			{
      			DS->EleIS[i] = 0;
			DS->EleSnow[i]=0;
			/* Note Two components can be separately read too */
			DS->EleSnowGrnd[i]=(1-DS->Ele[i].VegFrac)* DS->EleSnow[i];
			DS->EleSnowCanopy[i]=DS->Ele[i].VegFrac* DS->EleSnow[i];
      			NV_Ith_S(CV_Y, i) = 0;
      			NV_Ith_S(CV_Y, i + DS->NumEle) = 0;
      			NV_Ith_S(CV_Y, i + 2*DS->NumEle) = DS->Ele[i].zmax - DS->Ele[i].zmin -0.1;
			}
    		for(i=0; i<DS->NumRiv; i++)
    			{
      			NV_Ith_S(CV_Y, i + 4*DS->NumEle) = 0;
			/* Note once the element beneath river is incorporated in decomposition and .mesh file, initialization should be perfomed based on the location data instead of average of neighbor properties */
			NV_Ith_S(CV_Y, i + 4*DS->NumEle+DS->NumRiv)=(DS->Ele[i+DS->NumEle].zmax - DS->Ele[i+DS->NumEle].zmin) -0.1;
    			}
  		}
  	/* data initialization mode */  
  	else if (CS->init_type == 1)
  		{
    			for(i=0; i<DS->NumEle; i++)
    				{
      				DS->EleIS[i] = DS->Ele_IC[i].interception;
				DS->EleSnow[i]=DS->Ele_IC[i].snow;
				/* Note Two components can be separately read too */
				DS->EleSnowGrnd[i]=(1-DS->Ele[i].VegFrac)* DS->EleSnow[i];
				DS->EleSnowCanopy[i]=DS->Ele[i].VegFrac* DS->EleSnow[i];
      				NV_Ith_S(CV_Y, i) = DS->Ele_IC[i].surf;
		#ifdef SUB_SURF_RIV
      				NV_Ith_S(CV_Y, i + DS->NumEle) = DS->Ele_IC[i].unsat+0.1;
				NV_Ith_S(CV_Y, i + 2*DS->NumEle) = DS->Ele_IC[i].sat; 
			#ifdef LAYER3
                                NV_Ith_S(CV_Y, i + 3*DS->NumEle)=(DS->Ele[i].zmax-DS->Ele[i].zmin-DS->Ele[i].infD-NV_Ith_S(CV_Y, i + 2*DS->NumEle))>0?0.75*(DS->Ele[i].zmax-DS->Ele[i].zmin-DS->Ele[i].infD-NV_Ith_S(CV_Y, i + 2*DS->NumEle)):0.1; 
			#endif
				}
                        for(i=0; i<DS->NumRiv; i++)
                                {
                                NV_Ith_S(CV_Y, i + DS->totele) = DS->Riv_IC[DS->Riv[i].IC-1].value;
				/* Note once the element beneath river is incorporated in decomposition and .mesh file, initialization should be perfomed based on the location data instead of average of neighbor properties */
                                NV_Ith_S(CV_Y, i + DS->totele+DS->NumRiv) =(DS->Ele[i+DS->NumEle].zmax - DS->Ele[i+DS->NumEle].zmin)-0.1;
                #endif
                                }

  		}  
  	/* hot start mode */
  	else
  		{
    		fn = (char *)malloc((strlen(filename)+5)*sizeof(char));
    		strcpy(fn, filename);
    		init_file = fopen(strcat(fn, ".init"), "r");
  
    		if(init_file == NULL)
    			{
      			printf("\n  Fatal Error: %s.init is in use or does not exist!\n", filename);
      			exit(1);
    			}
    		else
    			{
      			for(i=0; i<DS->NumEle; i++)
      				{
        			fscanf(init_file, "%lf %lf %lf %lf %lf %lf", &DS->EleIS[i],&DS->EleSnow[i],&tempvalue1,&tempvalue2,&tempvalue3,&tempvalue4);
                                DS->EleSnowGrnd[i]=0*(1-DS->Ele[i].VegFrac)* DS->EleSnow[i];
                                DS->EleSnowCanopy[i]=0*DS->Ele[i].VegFrac* DS->EleSnow[i];
				NV_Ith_S(CV_Y, i)=tempvalue1;
		#ifdef SUB_SURF_RIV
				NV_Ith_S(CV_Y, i + DS->NumEle)=tempvalue2;
				NV_Ith_S(CV_Y, i + 2*DS->NumEle)=tempvalue3;
			#ifdef LAYER3	
//				NV_Ith_S(CV_Y, i + 3*DS->NumEle)=(DS->Ele[i].zmax-DS->Ele[i].zmin-DS->Ele[i].infD-NV_Ith_S(CV_Y, i + 2*DS->NumEle))>0?0.5*(DS->Ele[i].zmax-DS->Ele[i].zmin-DS->Ele[i].infD-NV_Ith_S(CV_Y, i + 2*DS->NumEle)):0.1;
				NV_Ith_S(CV_Y, i + 3*DS->NumEle)=tempvalue4;
			#endif
      				}  
      			for(i=0; i<DS->NumRiv; i++)
      				{
        			fscanf(init_file, "%lf %lf",&tempvalue1,&tempvalue2);
        			NV_Ith_S(CV_Y, i + DS->totele) = tempvalue1;
				NV_Ith_S(CV_Y, i + DS->totele+DS->NumRiv) = (DS->Ele[i+DS->NumEle].zmax - DS->Ele[i+DS->NumEle].zmin);
//				NV_Ith_S(CV_Y, i + DS->totele+DS->NumRiv) = tempvalue2; //(DS->Ele[i+DS->NumEle].zmax - DS->Ele[i+DS->NumEle].zmin);
		#endif
      				} 
    			}
    		fclose(init_file); 
  		}

        printf("done.\n");
	}


