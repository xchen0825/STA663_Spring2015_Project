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
#include "nvector_serial.h"
#include "sundialstypes.h"   

#include "f_functions.h"

void f(realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS)
	{
  	int i, j,ieleBC,inabr,inabr_left,inabr_right,inabr_plus_2ele,ileft,iright,i_plus_ele,i_plus_2ele,i_plus_3ele,iright_plus_2ele, ileft_plus_2ele, i_plus_totele, i_plus_totele1riv;
  	realtype Delta, Gamma,tmpVar,tmpVar1;
  	realtype Rn, T, Vel, RH, VP,P,LAI,zero_dh,cnpy_h,rl,r_a,r_s,alpha_r,f_r,eta_s,beta_s,Rmax;
  	realtype tempVar1;
	realtype Avg_Sf,Distance;
   	realtype Cwr,Perem, Perem_down,Avg_Rough,Avg_Perem,Avg_Y_Riv,Wid,Wid_down,Avg_Wid;
   	realtype nabrAqDepth,AquiferDepth, Deficit,elemSatn,satKfunc,effKI,effK_unsat,TotalY_Ele,TotalY_Ele_down,TotalY_unsat;
  	realtype *Y, *DY;
  	Model_Data MD;
  	Y = NV_DATA_S(CV_Y);
  	DY = NV_DATA_S(CV_Ydot);
  	MD = (Model_Data) DS;

	/* Initialization of temporary state variables */
       	MD->totriv=2*MD->NumRiv;
        #ifdef SURF_RIV
        MD->totriv=MD->NumRiv;
        #endif
	for(i=0; i<MD->totele+MD->totriv; i++)
  		{
		MD->DummyY[i]=(Y[i]>=0)?Y[i]:0;
		DY[i]=0;
  		if(i<MD->NumRiv)
			{
			MD->FluxRiv[i][0]=0;
			MD->FluxRiv[i][10]=0;
			}
   	#ifdef DIFFUSION
		if(i<MD->NumEle)
			{
			for(j=0;j<3;j++)
				{
                		inabr=MD->Ele[i].nabr[j]-1;
                		ieleBC=-(MD->Ele[i].BC[j]/4)-1;
			      	MD->Ele[i].surfH[j]=(inabr>-1)?((-ieleBC>0)?(MD->Ele[inabr].zmax+MD->DummyY[inabr]):((MD->DummyY[ieleBC+MD->totele]>MD->Riv[ieleBC].depth)?MD->Riv[ieleBC].zmin+MD->DummyY[ieleBC+MD->totele]:MD->Riv[ieleBC].zmax)):((MD->Ele[i].BC[j]!=1)?(MD->Ele[i].zmax+MD->DummyY[i]):Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j])-1], t));
			   //  MD->Ele[i].surfH[j]=0;
                        	}
                	MD->Ele[i].dhBYdx=-1*(MD->Ele[i].surfY[2]*(MD->Ele[i].surfH[1]-MD->Ele[i].surfH[0])+MD->Ele[i].surfY[1]*(MD->Ele[i].surfH[0]-MD->Ele[i].surfH[2])+MD->Ele[i].surfY[0]*(MD->Ele[i].surfH[2]-MD->Ele[i].surfH[1]))/(MD->Ele[i].surfX[2]*(MD->Ele[i].surfY[1]-MD->Ele[i].surfY[0])+MD->Ele[i].surfX[1]*(MD->Ele[i].surfY[0]-MD->Ele[i].surfY[2])+MD->Ele[i].surfX[0]*(MD->Ele[i].surfY[2]-MD->Ele[i].surfY[1]));
                	MD->Ele[i].dhBYdy=-1*(MD->Ele[i].surfX[2]*(MD->Ele[i].surfH[1]-MD->Ele[i].surfH[0])+MD->Ele[i].surfX[1]*(MD->Ele[i].surfH[0]-MD->Ele[i].surfH[2])+MD->Ele[i].surfX[0]*(MD->Ele[i].surfH[2]-MD->Ele[i].surfH[1]))/(MD->Ele[i].surfY[2]*(MD->Ele[i].surfX[1]-MD->Ele[i].surfX[0])+MD->Ele[i].surfY[1]*(MD->Ele[i].surfX[0]-MD->Ele[i].surfX[2])+MD->Ele[i].surfY[0]*(MD->Ele[i].surfX[2]-MD->Ele[i].surfX[1]));
			}
   	#endif
  		}	
	/* Lateral Flux Calculation between Triangular elements Follows  */
	for(i=0; i<MD->NumEle; i++)
  		{
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
                    		GradCalc(MD->DummyY[inabr_plus_2ele], MD->Ele[inabr].zmin, MD->DummyY[i_plus_2ele], MD->Ele[i].zmin, MD->Ele[i].Dist[j],0,0,0,0,ELE_ELE);
        			/* take care of macropore effect */
                    		nabrAqDepth=(MD->Ele[inabr].zmax-MD->Ele[inabr].zmin);
                    		avgKH(MD->Ele[i].Macropore,MD->DummyY[i_plus_2ele],AquiferDepth,MD->Ele[i].macD,MD->Ele[i].macKsatH,MD->Ele[i].vAreaF,MD->Ele[i].KsatH,MD->Ele[inabr].Macropore,MD->DummyY[inabr_plus_2ele],nabrAqDepth,MD->Ele[inabr].macD,MD->Ele[inabr].macKsatH,MD->Ele[inabr].vAreaF,MD->Ele[inabr].KsatH);
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
                    		GradCalc(MD->DummyY[inabr], MD->Ele[inabr].zmax, MD->DummyY[i], MD->Ele[i].zmax, MD->Ele[i].Dist[j],0,0,0,0,ELE_ELE);
                    		Avg_Sf=sqrt(pow(MD->Ele[i].dhBYdx,2)+pow(MD->Ele[i].dhBYdy,2));
                    		Avg_Sf=(Avg_Sf>EPS_5)?Avg_Sf:EPS_5;
                    	#else
                    		GradCalc(MD->DummyY[inabr], MD->Ele[inabr].zmax-MD->DummyY[inabr], MD->DummyY[i], MD->Ele[i].zmax-MD->DummyY[i], MD->Ele[i].Dist[j],0,0,0,0,ELE_ELE);
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
                        		effK=effKH(MD->Ele[i].Macropore,MD->DummyY[i_plus_2ele],AquiferDepth,MD->Ele[i].macD,MD->Ele[i].macKsatH,MD->Ele[i].vAreaF,MD->Ele[i].KsatH);
           				Avg_Ksat = effK;
                        		GradCalc(Interpolation(&MD->TSD_EleBC[(MD->Ele[i].BC[j])-1], t)- MD->Ele[i].zmin, MD->Ele[i].zmin, MD->DummyY[i_plus_2ele], MD->Ele[i].zmin, MD->Ele[i].Dist[j],0,0,0,0,ELE_ELE);
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
		fluxCalc_Ele(MD,i,t); 
   	#ifdef LAYER3
                i_plus_3ele=i+3*MD->NumEle;
                DY[i_plus_3ele] = MD->RechargeI[i]-MD->Recharge[i];
        #endif
		DY[i] = MD->EleNetPrep[i] - MD->EleViR[i]-((MD->DummyY[i]<EPS_3)?0:MD->EleET[i][2]);
		DY[i_plus_ele] = MD->EleViR[i]-MD->RechargeI[i]-((MD->DummyY[i]<EPS_3)?MD->EleET[i][2]:0);
		DY[i_plus_2ele]= MD->Recharge[i];
		if(MD->DummyY[i_plus_2ele]>AquiferDepth-MD->Ele[i].RzD)
			{
			DY[i_plus_2ele]=DY[i_plus_2ele]-MD->EleET[i][1];
			}
		else 
			{
		#ifdef LAYER2
			DY[i_plus_ele] = DY[i_plus_ele]-MD->EleET[i][1];
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
					DY[i_plus_3ele]=DY[i_plus_3ele]-MD->EleET[i][1];
					}
				}
			else
				{
				if (s2<=EPS_4)
					{
					DY[i_plus_ele] = DY[i_plus_ele]-MD->EleET[i][1];
					}
				else
					{
					tmpVar1=s1/(s1+s2);
                        		DY[i_plus_ele] = DY[i_plus_ele]-tmpVar1*MD->EleET[i][1];
                        		DY[i_plus_3ele] = DY[i_plus_3ele]-(1-tmpVar1)*MD->EleET[i][1];
					}
				}
		#endif
			}
#elif SURF_RIV
		DY[i] = MD->EleNetPrep[i];
#endif
#endif
  		} 
	/* Lateral Flux Calculation between River-River and River-Triangular elements Follows */   
	for(i=0; i<MD->NumRiv; i++)
  		{
		fluxCalc_Riv(MD,i,t);
  		}
	for(i=0; i<MD->NumEle; i++)
    		{   
        	i_plus_ele=i+MD->NumEle;
        	i_plus_2ele=i+2*MD->NumEle;
     		for(j=0; j<3; j++)
      			{
        		DY[i] =  DY[i] - MD->FluxSurf[i][j]/MD->Ele[i].area;
		#ifdef SUB_SURF_RIV
        		MD->FluxSource[i] =(MD->Ele[i].source>0)?Interpolation(&MD->TSD_Source[MD->Ele[i].source - 1], t):0;           //xchen_20141231
                        if (MD->DummyY[i_plus_2ele]>0)                                          //xchen_20141021
                        {
                        DY[i_plus_2ele] = DY[i_plus_2ele] - MD->FluxSub[i][j]/MD->Ele[i].area - MD->FluxSource[i];      //xchen_20141006
                        }
		#endif
      			}
		DY[i]=DY[i]/(UNIT_C);
	#ifdef NO_UNSAT
		DY[i_plus_2ele] = DY[i_plus_2ele] + MD->ElePrep[i];
	#endif
#ifdef SUB_SURF_RIV
	      	DY[i_plus_ele] = DY[i_plus_ele]/(MD->Ele[i].Porosity*UNIT_C);
      		DY[i_plus_2ele] = DY[i_plus_2ele]/(MD->Ele[i].Porosity*UNIT_C);
	#ifdef LAYER3
        	i_plus_3ele=i+3*MD->NumEle;
      		DY[i_plus_3ele] = DY[i_plus_3ele]/(MD->Ele[i].Porosity*UNIT_C);
	#endif
#endif
    		}
   	for(i=0; i<MD->NumRiv; i++)
    		{
        	i_plus_ele=i+MD->NumEle;
        	i_plus_totele=i+MD->totele;
        	i_plus_totele1riv=i+MD->totele+MD->NumRiv;
		for(j=0;j<=6;j++)
			{
			/* Note the limitation due to d(v)/dt=a*dy/dt+y*da/dt for CS other than rectangle */
			DY[i_plus_totele] = DY[i_plus_totele]-MD->FluxRiv[i][j]/(MD->Riv[i].Length*CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->Riv[i].depth,MD->Riv[i].coeff,3));
			}
		DY[i_plus_totele] = DY[i_plus_totele]/(UNIT_C);
	#ifdef SUB_SURF_RIV
		DY[i_plus_totele1riv] = DY[i_plus_totele1riv] -MD->FluxRiv[i][7] -MD->FluxRiv[i][8]-MD->FluxRiv[i][9] -MD->FluxRiv[i][10]+MD->FluxRiv[i][6];
      		DY[i_plus_totele1riv] = DY[i_plus_totele1riv]/(MD->Ele[i_plus_ele].Porosity*MD->Riv[i].Length*CS_AreaOrPerem(MD->Riv_Shape[MD->Riv[i].shape - 1].interpOrd,MD->Riv[i].depth,MD->Riv[i].coeff,3)*UNIT_C);
	#endif
    		}
//	printf("\nf  %lf, %lf, %lf, %lf, %lf, %lf",MD->FluxRiv[0][0],MD->FluxRiv[1][0],MD->FluxRiv[0][1],MD->FluxRiv[1][1],MD->FluxRiv[1][2],MD->FluxRiv[1][3]);
	}

