#include "pihm.h"
realtype returnVal(realtype rArea, realtype rPerem, realtype eqWid,realtype ap_Bool)
        {
        if(ap_Bool==1)
                {
                return rArea;
                }
        else if(ap_Bool==2)
                {
                return rPerem;
                }
        else
                {
                return eqWid;
                }
        }
realtype CS_AreaOrPerem(int rivOrder, realtype rivDepth, realtype rivCoeff, realtype a_pBool)
        {
        realtype rivArea, rivPerem, eq_Wid;
        switch(rivOrder)
                {
                case 1:
                        rivArea = rivDepth*rivCoeff;
                        rivPerem= 2.0*rivDepth+rivCoeff;
                        eq_Wid=rivCoeff;
                        return returnVal(rivArea, rivPerem, eq_Wid, a_pBool);
                case 2:
                        rivArea = pow(rivDepth,2)/rivCoeff;
                        rivPerem = 2.0*rivDepth*pow(1+pow(rivCoeff,2),0.5)/rivCoeff;
                        eq_Wid=2.0*pow(rivDepth+EPS_2,1/(rivOrder-1))/pow(rivCoeff,1/(rivOrder-1));
                        return returnVal(rivArea, rivPerem, eq_Wid, a_pBool);
                case 3:
                        rivArea = 4*pow(rivDepth,1.5)/(3*pow(rivCoeff,0.5));
                        rivPerem =(pow(rivDepth*(1+4*rivCoeff*rivDepth)/rivCoeff,0.5))+(log(2*pow(rivCoeff*rivDepth,0.5)+pow(1+4*rivCoeff*rivDepth,0.5))/(2*rivCoeff));
                        eq_Wid=2.0*pow(rivDepth+EPS_2,1/(rivOrder-1))/pow(rivCoeff,1/(rivOrder-1));
                        return returnVal(rivArea, rivPerem, eq_Wid, a_pBool);
                case 4:
                        rivArea = 3*pow(rivDepth,4.0/3.0)/(2*pow(rivCoeff,1.0/3.0));
                        rivPerem = 2*((pow(rivDepth*(1+9*pow(rivCoeff,2.0/3.0)*rivDepth),0.5)/3)+(log(3*pow(rivCoeff,1.0/3.0)*pow(rivDepth,0.5)+pow(1+9*pow(rivCoeff,2.0/3.0)*rivDepth,0.5))/(9*pow(rivCoeff,1.0/3.0))));
                        eq_Wid=2.0*pow(rivDepth+EPS_2,1/(rivOrder-1))/pow(rivCoeff,1/(rivOrder-1));
                        return returnVal(rivArea, rivPerem, eq_Wid, a_pBool);
                default:
                        printf("\n Relevant Values entered are wrong");
                        printf("\n Depth: %lf\tCoeff: %lf\tOrder: %d\t",rivDepth, rivCoeff, rivOrder);
                        return 0;
                }
        }
OverlandFlow(realtype **flux, int loci, int locj, realtype avg_y, realtype grad_y, realtype avg_sf, realtype crossA, realtype avg_rough)
        {
        flux[loci][locj] = crossA*pow(avg_y, 2.0/3.0)*grad_y/(sqrt(fabs(avg_sf))*avg_rough);
//      flux[loci][locj] = (grad_y>0?1:-1)*crossA*pow(avg_y, 2.0/3.0)*sqrt(fabs(grad_y))/(avg_rough);
        }
OLFeleToriv(realtype eleYtot,realtype EleZ,realtype cwr,realtype rivZmax,realtype rivYtot,realtype **fluxriv,int loci,int locj,realtype length)
        {
        realtype threshEle;
	if(fabs(eleYtot-EleZ)<EPS_3)
		{
		fluxriv[loci][locj]=0;
		}
	else
	{
        if(rivZmax < EleZ)
                {
                threshEle = EleZ;
                }
        else
                {
                threshEle = rivZmax;
                }
        if (rivYtot > eleYtot)
                {
                if (eleYtot > threshEle)
                        {
                        fluxriv[loci][locj] = cwr*2.0*sqrt(2*GRAV*UNIT_C*UNIT_C)*length*sqrt(rivYtot-eleYtot)*(rivYtot-threshEle)/3.0;
                        }
                else
                        {
                        if(threshEle<rivYtot)
                                {
                                fluxriv[loci][locj] = cwr*2.0*sqrt(2*GRAV*UNIT_C*UNIT_C)*length*sqrt(rivYtot-threshEle)*(rivYtot-threshEle)/3.0;
                                }
                        else
                                {
                                fluxriv[loci][locj]=0.0;
                                }
                        }
                }
        else
                {
                if (rivYtot > threshEle)
                        {
                        fluxriv[loci][locj] = -cwr*2.0*sqrt(2*GRAV*UNIT_C*UNIT_C)*length*sqrt(eleYtot - rivYtot)*(eleYtot - threshEle)/3.0;
                        }
                else
                        {
                        if(threshEle<eleYtot)
                                {
                                fluxriv[loci][locj] = -cwr*2.0*sqrt(2*GRAV*UNIT_C*UNIT_C)*length*sqrt(eleYtot - threshEle)*(eleYtot - threshEle)/3.0;
                                }
                        else
                                {
                                fluxriv[loci][locj]=0.0;
                                }
                        }
                }
	}
        }
/*
realtype avgY(realtype zi,realtype zinabr,realtype yi,realtype yinabr)
        {
        if(zinabr>zi)
                {
                if(zinabr>zi+yi)
                        {
                        return yinabr/2;
                        }
                else
                        {
                        return (yi+zi-zinabr+yinabr)/2;
                        }
                }
        else
                {
                if(zi>zinabr+yinabr)
                        {
                        return yi/2;
                        }
                else
                        {
                        return (yi+yinabr+zinabr-zi)/2;
                        }
                }
        }
*/
realtype avgY(realtype diff, realtype yi, realtype yinabr)
        {
        if(diff>0)
                {
                if(yi>EPS_3)
                        {
//                      return 0.5*(yi+yinabr);
//                      return ((yinabr>yi)?0:1.0*yi);  /* Note the if-else TRUE case can be possible only for Kinematic case */
                        return yi;
                        }
                else
                        {
                        return 0;
                        }
                }
        else
                {
                if(yinabr>EPS_3)
                        {
//                      return 0.5*(yi+yinabr);
//                      return ((yi>yinabr)?0:1.0*yinabr);  /* Note the if-else TRUE case can be possible only for Kinematic case */
                        return yinabr;
                        }
                else
                        {
                        return 0;
                        }
                }
        }
realtype effKV(realtype ksatFunc,realtype gradY,realtype macKV,realtype KV,realtype areaF)
        {
        if(ksatFunc>=0.98)
                {
                return (macKV*areaF+KV*(1-areaF)*ksatFunc);
                }
        else
                {
                if(fabs(gradY)*ksatFunc*KV<=1*KV*ksatFunc)
                        {
                        return KV*ksatFunc;
                        }
                else
                        {
                        if(fabs(gradY)*ksatFunc*KV<(macKV*areaF+KV*(1-areaF)*ksatFunc))
                                {
                                return (macKV*areaF*ksatFunc+KV*(1-areaF)*ksatFunc);
                                }
                        else
                                {
                                return (macKV*areaF+KV*(1-areaF)*ksatFunc);
                                }
                        }
                }
        }
realtype effKH(int mp,realtype tmpY, realtype aqDepth, realtype MacD, realtype MacKsatH, realtype areaF, realtype ksatH)
        {
                if(mp==1)
                        {
                        if(tmpY>aqDepth-MacD)
                                {
                                if(tmpY>aqDepth)
                                        {
                                        return  (MacKsatH*MacD*areaF+ksatH*(aqDepth-MacD*areaF))/aqDepth;
                                        }
                                else
                                        {
                                        return  (MacKsatH*(tmpY-(aqDepth-MacD))*areaF+ksatH*(aqDepth-MacD+(tmpY-(aqDepth-MacD))*(1-areaF)))/tmpY;
                                        }
                                }
                        else
                                {
                                return ksatH;
                                }
                        }
                else
                        {
                                return ksatH;
                        }
        }
realtype Interpolation(TSD *Data, realtype t)
        {
        int i, success;
        realtype result;
        i=Data->iCounter;
        success = 0;
        t=t/(UNIT_C);
        while(i<Data->length && t>Data->TS[i][0])
                {
                i++;
                }
        if(i==0)
                {
                /* t is smaller than the 1st node */
                result = Data->TS[i][1];
                }
        else if(i >= Data->length)
                {
                result = Data->TS[i-1][1];
                }
        else
                {
                result = ((Data->TS[i][0]-t)*Data->TS[i-1][1] + (t-Data->TS[i-1][0])*Data->TS[i][1])/(Data->TS[i][0]-Data->TS[i-1][0]);
                success = 1;
                }
        if(success == 0)
                {
                /*
                printf("\nWarning:  Extrapolation is used ...\n");
                */
                }
        return result;
        }
realtype GradCalc(realtype yNabr, realtype zNabr, realtype y, realtype z, realtype dist, int  interp, int interpNabr, realtype coeff, realtype coeffNabr, int booltype)
    	{
        realtype perem, peremNabr,zriv,totalY,totalY_nabr;
        totalY = y+z;
        totalY_nabr = yNabr+zNabr;
        Dif_Y=totalY - totalY_nabr;
        switch(booltype)
        	{
           	case ELE_ELE:
                Avg_Y=avgY(Dif_Y,y,yNabr);
                break;
           	case RIVSUB_RIVSUB:
                Avg_Y=avgY(Dif_Y,y,yNabr);
                break;
           	case RIV_ELESUB:
        	/* This is head at river edge representation */
//      	Avg_Y_Sub = ((MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].LeftEle-1].zmax-MD->Ele[MD->Riv[i].LeftEle-1].zmin)+MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle])>MD->Riv[i].zmin)?MD->Riv[i].zmin-(MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].LeftEle-1].zmax-MD->Ele[MD->Riv[i].LeftEle-1].zmin)):MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle];
        	/* This is head in neighboring cell represention */
                Avg_Y = zNabr>z?yNabr:((zNabr+yNabr)>z?(zNabr+yNabr-z):0);
                Avg_Y=avgY(Dif_Y,y,Avg_Y);
                break;
           	case RIVSUB_ELESUB:
		zriv=coeff;
        	/* This is head at river edge representation */
//      	Avg_Y_Sub = ((MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].LeftEle-1].zmax-MD->Ele[MD->Riv[i].LeftEle-1].zmin)+MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle])>MD->Riv[i].zmin)?MD->Riv[i].zmin-(MD->Riv[i].zmax-(MD->Ele[MD->Riv[i].LeftEle-1].zmax-MD->Ele[MD->Riv[i].LeftEle-1].zmin)):MD->DummyY[MD->Riv[i].LeftEle-1 + 2*MD->NumEle];
        	/* This is head in neighboring cell represention */
                Avg_Y = zNabr>zriv?0:((zNabr+yNabr)>zriv?(zriv-zNabr):yNabr);
                Avg_Y=avgY(Dif_Y,y,Avg_Y);
                break;
           	case RIVSURF_RIVSURF:
                perem = CS_AreaOrPerem(interp,y,coeff,2);
                peremNabr = CS_AreaOrPerem(interpNabr,yNabr,coeff,2);
                CrossA = CS_AreaOrPerem(interp,y,coeff,1);
                CrossANabr = CS_AreaOrPerem(interpNabr,yNabr,coeffNabr,1);
                y=perem*(CrossA/(perem*perem));
                yNabr=peremNabr*(CrossANabr/(peremNabr*peremNabr));
                Avg_Y=avgY(Dif_Y,y,yNabr);
//                Avg_Y=((perem+peremNabr)==0)?0:(CrossA+CrossANabr)/(perem+peremNabr);
		AvgCrossA=(Dif_Y>0)?CrossA:CrossANabr;
                break;
           	default:
                return 0;
		}
        Grad_Y = Dif_Y/dist;
	}

realtype avgKH(int macp,realtype y,realtype aqD,realtype macD,realtype macKsatH,realtype vAreaF,realtype KsatH,int macpnabr,realtype ynabr,realtype aqDnabr,realtype macDnabr,realtype macKsatHnabr,realtype vAreaFnabr,realtype KsatHnabr)
	{
        effK=effKH(macp,y,aqD,macD,macKsatH,vAreaF,KsatH);
        effKnabr=effKH(macpnabr,ynabr,aqDnabr,macDnabr,macKsatHnabr,vAreaFnabr,KsatHnabr);
        /* It should be weighted average. However, there is an ambiguity about distance used */
        Avg_Ksat=0.5*(effK+effKnabr);
    	}

void conditionals(Model_Data DS, Control_Data *CS)
	{
	}
void fluxCalc_Ele(Model_Data DS, int i,realtype t,realtype *DummyY)
	{
		realtype S,Rn,G,T,Vel,RH,VP,LAI,rl,P,Gamma,Delta,r_a,Rmax,alpha_r,f_r,eta_s,elemSatn,beta_s,r_s,TotalY_unsat,effK_unsat,satKfunc,tmpVar,effKI,AquiferDepth,pet,mf;
		realtype elemSatn1,elemSatn2,beta_s1,beta_s2,r_s1,r_s2;
		int i_plus_ele,i_plus_2ele,i_plus_3ele;	
                i_plus_ele=i+DS->NumEle;
                i_plus_2ele=i+2*DS->NumEle;
                AquiferDepth=(DS->Ele[i].zmax-DS->Ele[i].zmin);
                /**************************************************************************************************/
                /* Evaporation Module: [2] is ET from OVLF/SUBF, [1] is Transpiration, [0] is ET loss from canopy */
                /**************************************************************************************************/
                /* Physical Unit Dependent. Change this */
		//DS->ElePrep[i] = Interpolation(&DS->TSD_Prep[DS->Ele[i].prep-1], t);
		Rn = Interpolation(&DS->TSD_Rn[DS->Ele[i].Rn-1], t);
        //      G = Interpolation(&DS->TSD_G[DS->Ele[i].G-1], t);
                T = Interpolation(&DS->TSD_Temp[DS->Ele[i].temp-1], t);
                Vel = Interpolation(&DS->TSD_WindVel[DS->Ele[i].WindVel-1], t);
                RH = Interpolation(&DS->TSD_Humidity[DS->Ele[i].humidity-1], t);
                VP = Interpolation(&DS->TSD_Pressure[DS->Ele[i].pressure-1], t);
                LAI = Interpolation(&DS->TSD_LAI[DS->Ele[i].LC-1], t);
                rl=Interpolation(&DS->TSD_RL[DS->Ele[i].LC-1], t);
		mf=Interpolation(&DS->TSD_MeltF[DS->Ele[i].meltF-1], t);
/*
                Actual Expression is as follows
                P = 101.325*pow(10,3)*pow((293-0.0065*DS->Ele[i].zmax)/293,5.26);
                Gamma = P*1.0035*0.92/(0.622*2441);
                However we use the following expression to avoid repetetive computations
*/
//              Write the reference here
                P = 101325*pow((1.0-0.00002218*DS->Ele[i].zmax),5.26);
                Gamma = P*0.000608;
                Delta = 2503*pow(10,3)*exp(17.27*T/(T+237.3))/(pow(237.3 + T, 2));
//              Write the reference here
                r_a = log(DS->Ele[i].windH/rl)*log(10*DS->Ele[i].windH/rl)/(Vel*0.16);
//		S=9.81*DS->Ele[i].VegFrac/pow(Vel/86400.0,2);
//		r_a_new=(1+0.5*((S<10)?S:10))*log(DS->Ele[i].windH/rl)*log(10*DS->Ele[i].windH/rl)/(Vel*0.16);		
		//r_a_new=0.004/(1+0.5*((S<10)?S:10));
		pet=Rn*(1-DS->Ele[i].Albedo)*Delta+(1.2*1003.5*((VP/RH)-VP)/r_a);
                DS->EleET[i][2] = DS->pcCal.Et2*(pet)/(2441000000.0*(Delta+Gamma));
		//pet=Rn*(1-DS->Ele[i].Albedo)*Delta+(1.2*1003.5*((VP/RH)-VP)/r_a);
                //DS->EleET[i][2] = DS->pcCal.Et2*(LAI/DS->Ele[i].LAImax)*(1-DS->Ele[i].VegFrac)*(pet)/(2441000000.0*(Delta+Gamma));

		if(LAI>0.0)
                        {
/*
                        Actual Expression is as follows
                        Rmax = 5000.0/(60*UNIT_C);
                        However we use the following expression to avoid repetetive computations
*/
                        Rmax = 83.334/(UNIT_C);         /* Unit day_per_m */
                        f_r= 1.1*Rn*(1-exp(-LAI))/(DS->Ele[i].Rs_ref*LAI);
                        alpha_r= (1+f_r)/(1+(DS->Ele[i].Rmin/Rmax));
                        eta_s= 1- 0.0016*(pow((24.85-T),2));
                        if(AquiferDepth-DummyY[i_plus_2ele]<DS->Ele[i].RzD)
                                {
                                elemSatn=1.0;
				beta_s=elemSatn;
                       		r_s=((DS->Ele[i].Rmin*alpha_r/(beta_s*LAI*pow(eta_s,4)))> Rmax)?Rmax:(DS->Ele[i].Rmin*alpha_r/(beta_s*LAI*pow(eta_s,4)));
                        	DS->EleET[i][1] = DS->pcCal.Et1*((1-pow(((DS->EleIS[i]+DS->EleSnowCanopy[i]<0)?0:(DS->EleIS[i]+DS->EleSnowCanopy[i]))/(DS->EleISmax[i]+DS->EleISsnowmax[i]),2.0/3))*(pet)/(2441000000.0*(Delta+Gamma*(1+r_s/r_a))));
                                }
                         else
                                {
                	 #ifdef LAYER2        
		        	elemSatn =DummyY[i_plus_ele]/DS->Ele[i].infD;
				beta_s= ((elemSatn<EPS_4)?EPS_4:((elemSatn>1.0)?1:elemSatn));
                       		r_s=((DS->Ele[i].Rmin*alpha_r/(beta_s*LAI*pow(eta_s,4)))> Rmax)?Rmax:(DS->Ele[i].Rmin*alpha_r/(beta_s*LAI*pow(eta_s,4)));
                        	DS->EleET[i][1] = ((beta_s<=EPS_4)?0:1.0)*DS->pcCal.Et1*((1-pow(((DS->EleIS[i]+DS->EleSnowCanopy[i]<0)?0:(DS->EleIS[i]+DS->EleSnowCanopy[i]))/(DS->EleISmax[i]+DS->EleISsnowmax[i]),2.0/3))*(pet)/(2441000000.0*(Delta+Gamma*(1+r_s/r_a))));
			#endif
                       	#ifdef LAYER3
                		i_plus_3ele=i+3*DS->NumEle;
				elemSatn1 =DummyY[i_plus_ele]/DS->Ele[i].infD;
                                elemSatn2 =DummyY[i_plus_3ele]/(((AquiferDepth-DS->Ele[i].infD-DummyY[i_plus_2ele])<DS->Ele[i].infD)?DS->Ele[i].infD:(AquiferDepth-DS->Ele[i].infD-DummyY[i_plus_2ele]));
                                beta_s1= ((elemSatn1<EPS_4)?EPS_4:((elemSatn1>1.0)?1:elemSatn1));
                                beta_s2= ((elemSatn2<EPS_4)?EPS_4:((elemSatn2>1.0)?1:elemSatn2));
                                s1 = beta_s1;
                                s2 = beta_s2;
                                r_s1=((DS->Ele[i].Rmin*alpha_r/(beta_s1*LAI*pow(eta_s,4)))> Rmax)?Rmax:(DS->Ele[i].Rmin*alpha_r/(beta_s1*LAI*pow(eta_s,4)));
                                r_s2=((DS->Ele[i].Rmin*alpha_r/(beta_s2*LAI*pow(eta_s,4)))> Rmax)?Rmax:(DS->Ele[i].Rmin*alpha_r/(beta_s2*LAI*pow(eta_s,4)));
                                r_s=(r_s1+r_s2)/2;
                                DS->EleET[i][1]=(beta_s1<=EPS_4)?((beta_s2<=EPS_4)?0:DS->pcCal.Et1*((1-pow(((DS->EleIS[i]+DS->EleSnowCanopy[i]<0)?0:(DS->EleIS[i]+DS->EleSnowCanopy[i]))/(DS->EleISmax[i]+DS->EleISsnowmax[i]),2.0/3))*(pet)/(2441000000.0*(Delta+Gamma*(1+r_s2/r_a))))):(beta_s2<=EPS_4)?DS->pcCal.Et1*((1-pow(((DS->EleIS[i]+DS->EleSnowCanopy[i]<0)?0:(DS->EleIS[i]+DS->EleSnowCanopy[i]))/(DS->EleISmax[i]+DS->EleISsnowmax[i]),2.0/3))*(pet)/(2441000000.0*(Delta+Gamma*(1+r_s1/r_a)))):DS->pcCal.Et1*((1-pow(((DS->EleIS[i]+DS->EleSnowCanopy[i]<0)?0:(DS->EleIS[i]+DS->EleSnowCanopy[i]))/(DS->EleISmax[i]+DS->EleISsnowmax[i]),2.0/3))*(pet)/(2441000000.0*(Delta+Gamma*(1+r_s/r_a))));
			#endif 
	        		}
                       	} 
                else
                        {
                        DS->EleET[i][1] =0.0;
                        }
                /* Note: Assumption is OVL flow depth less than EPS_3 is immobile water */
                /* Calculation of saturation of top layer */
                elemSatn=DummyY[i_plus_ele]/DS->Ele[i].infD;
		elemSatn=((elemSatn<EPS_4)?EPS_4:((elemSatn>1.0)?1:elemSatn));
                DS->EleET[i][2]=(DummyY[i]<EPS_3)?((elemSatn>EPS_4)?elemSatn*DS->EleET[i][2]:0):DS->EleET[i][2];
               
		Avg_Y=-(pow(pow(1/elemSatn,DS->Ele[i].Beta/(DS->Ele[i].Beta-1))-1,1/DS->Ele[i].Beta)/DS->Ele[i].Alpha);
                TotalY_unsat=((elemSatn>=0.5)?(DummyY[i_plus_ele]-0.5*DS->Ele[i].infD):0)+((Avg_Y<MINpsi)?MINpsi:Avg_Y)+DS->Ele[i].zmin+AquiferDepth-0.5*DS->Ele[i].infD;                /* Calculation of gradient between overland flow and top unsaturated zone */
                Grad_Y=(DummyY[i]+DS->Ele[i].zmax-TotalY_unsat)/(0.5*DS->Ele[i].infD);
                Grad_Y= (Grad_Y>0)?((DummyY[i]<EPS_3)?0:Grad_Y):((DummyY[i_plus_ele]>EPS_3)?Grad_Y:0);
                /* The effective conductivity for infiltration is calculated based on saturated condition */
/*
                Following are the expressions being evaluated in the next line
                elemSatn=1.0;
                satKfunc=pow(elemSatn,0.5)*pow(-1+pow(1-pow(elemSatn,DS->Ele[i].Beta/(DS->Ele[i].Beta-1)),(DS->Ele[i].Beta-1)/DS->Ele[i].Beta),2);
                satKfunc=1.0;
                effK_unsat=(DS->Ele[i].Macropore==1)?effKV(satKfunc,Grad_Y,DS->Ele[i].macKsatV,DS->Ele[i].infKsatV,DS->Ele[i].hAreaF):DS->Ele[i].infKsatV;
*/
                effK_unsat=(DS->Ele[i].Macropore==1)? DS->Ele[i].macKsatV*DS->Ele[i].hAreaF+(1-DS->Ele[i].hAreaF)*DS->Ele[i].infKsatV:DS->Ele[i].infKsatV;
                /* Infiltration/exfiltration rate calculated */
		DS->EleViR[i] = ((DS->EleNetPrep[i]<=(effK_unsat*Grad_Y)) && (DummyY[i]<=EPS_3))?DS->EleNetPrep[i]:(effK_unsat*Grad_Y);
//		DS->EleViR[i] = (effK_unsat)*Grad_Y;
		/* Calculation of effective conductivity of top layer */
                elemSatn=DummyY[i_plus_ele]/DS->Ele[i].infD;
                elemSatn=(elemSatn<EPS_4)?EPS_4:(elemSatn>1?1:elemSatn);
 //             satKfunc=pow(elemSatn,0.5)*pow(-1+pow(1-pow(elemSatn,DS->Ele[i].Beta/(DS->Ele[i].Beta-1)),(DS->Ele[i].Beta-1)/DS->Ele[i].Beta),2);
                effK_unsat=(DS->Ele[i].Macropore==1)?effKV(elemSatn,Grad_Y,DS->Ele[i].macKsatV,DS->Ele[i].infKsatV,DS->Ele[i].hAreaF):DS->Ele[i].infKsatV*elemSatn;        
//              effK_unsat=(DS->Ele[i].Macropore==1)?effKV(satKfunc,Grad_Y,DS->Ele[i].macKsatV,DS->Ele[i].infKsatV,DS->Ele[i].hAreaF):DS->Ele[i].infKsatV*satKfunc;        
	#ifdef LAYER3
                /* Calculation of saturation of intermediate layer */
                tmpVar=((AquiferDepth-DS->Ele[i].infD-DummyY[i_plus_2ele])<DS->Ele[i].infD)?DS->Ele[i].infD:(AquiferDepth-DS->Ele[i].infD-DummyY[i_plus_2ele]);
                i_plus_3ele=i+3*DS->NumEle;					//xchen_20141231
		elemSatn=DummyY[i_plus_3ele]/tmpVar;
                elemSatn=(elemSatn<EPS_4)?EPS_4:(elemSatn>1?1:elemSatn);
                Avg_Y=-(pow(pow(1/elemSatn,DS->Ele[i].Beta/(DS->Ele[i].Beta-1))-1,1/DS->Ele[i].Beta)/DS->Ele[i].Alpha);                
		TotalY=((Avg_Y<MINpsi)?MINpsi:Avg_Y)+DS->Ele[i].zmin-0.5*tmpVar+AquiferDepth-DS->Ele[i].infD+(((elemSatn>0.5)&&(tmpVar>DS->Ele[i].infD))?(DummyY[i_plus_3ele]-0.5*(AquiferDepth-DS->Ele[i].infD-DummyY[i_plus_2ele])):0)+((tmpVar==DS->Ele[i].infD)?((DummyY[i_plus_3ele]-0.5*tmpVar)>0?(DummyY[i_plus_3ele]-0.5*tmpVar):0.0*tmpVar):0);
                /* Calculation of gradient between top and intermediate layer */
                Grad_Y=(TotalY_unsat-TotalY)/(0.5*(DS->Ele[i].infD+tmpVar));
		
		/* Calculation of effective conductivity of the intermediate layer */
//                satKfunc=pow(elemSatn,0.5)*pow(-1+pow(1-pow(elemSatn,DS->Ele[i].Beta/(DS->Ele[i].Beta-1)),(DS->Ele[i].Beta-1)/DS->Ele[i].Beta),2);
                effKI=(DS->Ele[i].Macropore==1)?effKV(elemSatn,Grad_Y,DS->Ele[i].macKsatV,DS->Ele[i].KsatV,DS->Ele[i].hAreaF):DS->Ele[i].KsatV*elemSatn;
  //              effKI=(DS->Ele[i].Macropore==1)?effKV(satKfunc,Grad_Y,DS->Ele[i].macKsatV,DS->Ele[i].KsatV,DS->Ele[i].hAreaF):DS->Ele[i].KsatV*satKfunc;
                /* Calculation of exchange flux between top and intermediate layer */
                /* Using upwind method */
                DS->RechargeI[i] =(Grad_Y>0)?((DummyY[i_plus_ele]>EPS_3)?effK_unsat*Grad_Y:0):((DummyY[i_plus_3ele]>EPS_3)?effKI*Grad_Y:0);
                /* Calculation of gradient between intermediate and bottom layer */
                Grad_Y=(TotalY-(DummyY[i_plus_2ele]+DS->Ele[i].zmin))/(0.5*(AquiferDepth-DS->Ele[i].infD));
                /* Calculation of effective conductivity of the bottom layer */
                effKI=(DS->Ele[i].Macropore==1)?effKV(elemSatn,Grad_Y,DS->Ele[i].macKsatV,DS->Ele[i].KsatV,DS->Ele[i].hAreaF):DS->Ele[i].KsatV*elemSatn;
                tmpVar=(DummyY[i_plus_2ele]>AquiferDepth-DS->Ele[i].macD)?effKI:DS->Ele[i].KsatV*elemSatn;
//                tmpVar=(DS->DummyY[i_plus_2ele]>AquiferDepth-DS->Ele[i].macD)?effKI:DS->Ele[i].KsatV*satKfunc;
/*
                Following statements will be rewritted to avoid repetetive computations
                effKI calculates top layer effective conductivity
                satKfunc=1.0;
                effKI=(DS->Ele[i].Macropore==1)?effKV(satKfunc,Grad_Y,DS->Ele[i].macKsatV,DS->Ele[i].KsatV,DS->Ele[i].hAreaF):DS->Ele[i].KsatV*satKfunc;
                effK=(DS->Ele[i].Macropore==1)?((DS->DummyY[i_plus_2ele]>AquiferDepth-DS->Ele[i].macD)?effKI:DS->Ele[i].KsatV*satKfunc):DS->Ele[i].KsatV*satKfunc;
*/
                effKI=(DS->Ele[i].Macropore==1)? DS->Ele[i].macKsatV*DS->Ele[i].hAreaF+(1-DS->Ele[i].hAreaF)*DS->Ele[i].KsatV:DS->Ele[i].KsatV;
                effK=(DS->Ele[i].Macropore==1)?((DummyY[i_plus_2ele]>AquiferDepth-DS->Ele[i].macD)?effKI:DS->Ele[i].KsatV):DS->Ele[i].KsatV;
                /* Calculation of exchange flux between top and bottom layer */
                /* Using upwind method */
                DS->Recharge[i] =(Grad_Y>0)?((DummyY[i_plus_3ele]>EPS_3)?tmpVar*Grad_Y:0):((DummyY[i_plus_2ele]>EPS_3)?effK*Grad_Y:0);
        #else
/*                Effective conductivity of the top layer. The conductivity should be infKsatV assuming that this is the conductivity of the top layer. That is why the line below is commented

*/
//              effK_unsat=(DS->Ele[i].Macropore==1)?effKV(satKfunc,Grad_Y,DS->Ele[i].macKsatV,DS->Ele[i].KsatV,DS->Ele[i].hAreaF):DS->Ele[i].KsatV*satKfunc;
                elemSatn=DummyY[i_plus_2ele]/(AquiferDepth-DS->Ele[i].infD);
                elemSatn=(elemSatn<EPS_4)?EPS_4:(elemSatn>1.0?1.0:elemSatn);
//                satKfunc=pow(elemSatn,0.5)*pow(-1+pow(1-pow(elemSatn,DS->Ele[i].Beta/(DS->Ele[i].Beta-1)),(DS->Ele[i].Beta-1)/DS->Ele[i].Beta),2);
                Avg_Y=-(pow(pow(1/elemSatn,DS->Ele[i].Beta/(DS->Ele[i].Beta-1))-1,1/DS->Ele[i].Beta)/DS->Ele[i].Alpha);                
		TotalY=((Avg_Y<MINpsi)?MINpsi:Avg_Y)+DS->Ele[i].zmin+DummyY[i_plus_2ele]+(DummyY[i_plus_2ele]>(AquiferDepth-DS->Ele[i].infD)?0:0.5*(AquiferDepth-DS->Ele[i].infD-DummyY[i_plus_2ele]));
                /* Calculation of gradient between intermediate and bottom layer */                
		Grad_Y=(TotalY_unsat-TotalY)/(0.5*DS->Ele[i].infD);                
		tmpVar=(DS->Ele[i].Macropore==1)?effKV(elemSatn,Grad_Y,DS->Ele[i].macKsatV,DS->Ele[i].KsatV,DS->Ele[i].hAreaF):DS->Ele[i].KsatV*elemSatn;
//		tmpVar=(DS->Ele[i].Macropore==1)?effKV(satKfunc,Grad_Y,DS->Ele[i].macKsatV,DS->Ele[i].KsatV,DS->Ele[i].hAreaF):DS->Ele[i].KsatV*satKfunc;
                effK=(DS->Ele[i].Macropore==1)?((DummyY[i_plus_2ele]>AquiferDepth-DS->Ele[i].macD)?tmpVar:DS->Ele[i].KsatV*elemSatn):DS->Ele[i].KsatV*elemSatn;
//               effK=(DS->Ele[i].Macropore==1)?((DS->DummyY[i_plus_2ele]>AquiferDepth-DS->Ele[i].macD)?tmpVar:DS->Ele[i].KsatV*satKfunc):DS->Ele[i].KsatV*satKfunc;
                /* Calculation of exchange flux between top and bottom layer */
                /* Using upwind method */
                DS->Recharge[i] =(Grad_Y>0)?((DummyY[i_plus_ele]>EPS_3)?effK_unsat*Grad_Y:0):((DummyY[i_plus_2ele]>EPS_3)?effK*Grad_Y:0);
                DS->RechargeI[i]=DS->Recharge[i];
	#endif
	}
void fluxCalc_Riv(Model_Data DS, int i,realtype t,realtype *DummyY)
	{
		int i_plus_ele,i_plus_totele,i_plus_totele1riv,ileft,iright,ileft_plus_2ele,iright_plus_2ele, inabr,inabr_left,inabr_right;
		realtype tmpVar,Wid, Wid_down, Avg_Wid, Avg_Rough, Avg_Sf, Perem, Avg_Perem, nabrAqDepth;
		i_plus_ele=i+DS->NumEle;
                i_plus_totele=i+DS->totele;
                i_plus_totele1riv=i+DS->totele+DS->NumRiv;
                ileft=DS->Riv[i].LeftEle-1;
                iright=DS->Riv[i].RightEle-1;
                ileft_plus_2ele=ileft+2*DS->NumEle;
                iright_plus_2ele=iright+2*DS->NumEle;
                TotalY = DummyY[i_plus_totele] + DS->Riv[i].zmin;
                if(DS->Riv[i].down > 0)
                        {
                        inabr=DS->Riv[i].down - 1;
                        inabr_left=DS->Riv[inabr].LeftEle-1;
                        inabr_right=DS->Riv[inabr].RightEle-1;
		#ifdef SUB_SURF_RIV
                        /************************************************************************/
                        /* Lateral Flux Calculation between Element Beneath River (EBR) and EBR */
                        /************************************************************************/
                        GradCalc(DummyY[inabr+DS->totele+DS->NumRiv],DS->Ele[inabr+DS->NumEle].zmin,DummyY[i_plus_totele1riv],DS->Ele[i_plus_ele].zmin,DS->Riv[i].Dist[0],0,0,0,0,RIVSUB_RIVSUB);
                        Wid = CS_AreaOrPerem(DS->Riv_Shape[DS->Riv[i].shape - 1].interpOrd,DS->Riv[i].depth,DS->Riv[i].coeff,3);
                        Wid_down = CS_AreaOrPerem(DS->Riv_Shape[DS->Riv[inabr].shape - 1].interpOrd,DS->Riv[inabr].depth,DS->Riv[inabr].coeff,3);
                        /* take care of macropore effect */
                        avgKH(DS->Ele[ileft].Macropore,DummyY[ileft_plus_2ele],DS->Ele[ileft].zmax-DS->Ele[ileft].zmin,DS->Ele[ileft].macD,DS->Ele[ileft].macKsatH,DS->Ele[ileft].vAreaF,DS->Ele[ileft].KsatH,DS->Ele[iright].Macropore,DummyY[iright_plus_2ele],DS->Ele[iright].zmax-DS->Ele[iright].zmin,DS->Ele[iright].macD,DS->Ele[iright].macKsatH,DS->Ele[iright].vAreaF,DS->Ele[iright].KsatH);
                        tmpVar=Avg_Ksat;
                        avgKH(DS->Ele[inabr_left].Macropore,DummyY[inabr_left+2*DS->NumEle],DS->Ele[inabr_left].zmax-DS->Ele[inabr_left].zmin,DS->Ele[inabr_left].macD,DS->Ele[inabr_left].macKsatH,DS->Ele[inabr_left].vAreaF,DS->Ele[inabr_left].KsatH,DS->Ele[inabr_right].Macropore,DummyY[inabr_right+2*DS->NumEle],DS->Ele[inabr_right].zmax-DS->Ele[inabr_right].zmin,DS->Ele[inabr_right].macD,DS->Ele[inabr_right].macKsatH,DS->Ele[inabr_right].vAreaF,DS->Ele[inabr_right].KsatH);
                        Avg_Ksat=0.5*(Avg_Ksat+tmpVar);
                        Avg_Wid=0.5*(Wid+Wid_down);
                        /* groundwater flow modeled by Darcy's law */
                        DS->FluxRiv[i][9] = Avg_Ksat*Grad_Y*Avg_Y*Avg_Wid;
                        /* accumulate to get in-flow for down segments: [10] for inflow, [9] for outflow */
                        DS->FluxRiv[inabr][10] = DS->FluxRiv[inabr][10] - DS->FluxRiv[i][9];
		#endif
                        /****************************************************************/
                        /* Lateral Flux Calculation between River-River element Follows */
                        /****************************************************************/
                #ifdef DIFFUSION
                        GradCalc(DummyY[inabr + DS->totele],DS->Riv[inabr].zmin,DummyY[i_plus_totele],DS->Riv[i].zmin,DS->Riv[i].Dist[0],DS->Riv_Shape[DS->Riv[i].shape - 1].interpOrd,DS->Riv_Shape[DS->Riv[inabr].shape - 1].interpOrd,DS->Riv[i].coeff,DS->Riv[inabr].coeff,RIVSURF_RIVSURF);
                #else
                        GradCalc(DummyY[inabr + DS->totele],DS->Riv[inabr].zmin-DummyY[inabr + DS->totele],DummyY[i_plus_totele],DS->Riv[i].zmin-DummyY[i_plus_totele],DS->Riv[i].Dist[0],DS->Riv_Shape[DS->Riv[i].shape - 1].interpOrd,DS->Riv_Shape[DS->Riv[inabr].shape - 1].interpOrd,DS->Riv[i].coeff,DS->Riv[inabr].coeff,RIVSURF_RIVSURF);
                #endif
//                      AvgCrossA=0.5*(CrossA+CrossANabr);
                        Avg_Rough = (DS->Riv[i].Rough + DS->Riv[inabr].Rough)/2.0;
//                      Avg_Sf = (Grad_Y>0)?Grad_Y:EPS;
                        Avg_Sf = (fabs(Grad_Y)>EPS_5)?Grad_Y:EPS_5;
                        /* INCLUDE CROSSADOWN */
                        OverlandFlow(DS->FluxRiv,i,1, Avg_Y,Grad_Y,Avg_Sf,AvgCrossA,Avg_Rough);
                        /* accumulate to get in-flow for down segments: [0] for inflow, [1] for outflow */
                        DS->FluxRiv[inabr][0] = DS->FluxRiv[inabr][0] - DS->FluxRiv[i][1];
                        }
                else
                        {
                        switch(DS->Riv[i].down)
                                {
                                case -1:
                                GradCalc(Interpolation(&DS->TSD_Riv[(DS->Riv[i].BC)-1], t),DS->Node[DS->Riv[i].ToNode-1].zmax -DS->Riv[i].depth,DummyY[i_plus_totele],DS->Riv[i].zmin,DS->Riv[i].Dist[0],DS->Riv_Shape[DS->Riv[i].shape - 1].interpOrd,DS->Riv_Shape[DS->Riv[i].shape - 1].interpOrd,DS->Riv[i].coeff,DS->Riv[i].coeff,RIVSURF_RIVSURF);
                                /* Dirichlet boundary condition */
                                /* Note: do i need to change else part here for diff wave */
                                Avg_Sf = (DS->RivMode==1)?Grad_Y:Grad_Y;;
                                Avg_Rough = DS->Riv_Mat[DS->Riv[i].material-1].Rough;
                                OverlandFlow(DS->FluxRiv,i,1, Avg_Y,Grad_Y,Avg_Sf,CrossA,Avg_Rough);
                                break;
                                case -2:
                                /* Neumann boundary condition */
                                DS->FluxRiv[i][1] = Interpolation(&DS->TSD_Riv[DS->Riv[i].BC-1], t);
                                break;
                                case -3:
                                /* zero-depth-gradient boundary conditions */
//                              Distance = sqrt(pow(DS->Riv[i].x - DS->Node[DS->Riv[i].ToNode-1].x, 2) + pow(DS->Riv[i].y - DS->Node[DS->Riv[i].ToNode-1].y, 2));
                                Perem = CS_AreaOrPerem(DS->Riv_Shape[DS->Riv[i].shape - 1].interpOrd,DummyY[i_plus_totele],DS->Riv[i].coeff,2);
                                Grad_Y = (DS->Riv[i].zmin - (DS->Node[DS->Riv[i].ToNode-1].zmax -DS->Riv[i].depth))/DS->Riv[i].Dist[0];
                                Avg_Rough = DS->Riv_Mat[DS->Riv[i].material-1].Rough;
                                Avg_Perem = Perem;
                                CrossA = CS_AreaOrPerem(DS->Riv_Shape[DS->Riv[i].shape - 1].interpOrd,DummyY[i_plus_totele],DS->Riv[i].coeff,1);
                                DS->FluxRiv[i][1] = sqrt(Grad_Y)*CrossA*((Avg_Perem>0)?pow(CrossA/Avg_Perem,2.0/3.0):0)/Avg_Rough;
                                break;
                                case -4:
                                /* Critical Depth boundary conditions */
                                CrossA = CS_AreaOrPerem(DS->Riv_Shape[DS->Riv[i].shape - 1].interpOrd,DummyY[i_plus_totele],DS->Riv[i].coeff,1);
                                DS->FluxRiv[i][1] = CrossA*sqrt(GRAV*UNIT_C*UNIT_C*DummyY[i_plus_totele1riv]);  /* Note the dependence on physical units */
                                break;
                                default:
                                printf("Fatal Error: River Routing Boundary Condition Type Is Wrong!");
                                exit(1);
                                }
                        /* Note: bdd condition for subsurface element can be changed. Assumption: No flow condition */
                        DS->FluxRiv[i][9]=0;
                        }
                if(DS->Riv[i].LeftEle > 0)
                        {
                        /*****************************************************************************/
                        /* Lateral Surface Flux Calculation between River-Triangular element Follows */
                        /*****************************************************************************/
                        OLFeleToriv(DummyY[ileft]+DS->Ele[ileft].zmax,DS->Ele[ileft].zmax,DS->Riv_Mat[DS->Riv[i].material-1].Cwr, DS->Riv[i].zmax,TotalY,DS->FluxRiv,i,2,DS->Riv[i].Length);
                        DS->FluxSurf[ileft][DS->Riv[i].lrEdge[0]] = -DS->FluxRiv[i][2];
		#ifdef SUB_SURF_RIV
                        /*********************************************************************************/
                        /* Lateral Sub-surface Flux Calculation between River-Triangular element Follows */
                        /*********************************************************************************/
                        GradCalc(DummyY[ileft_plus_2ele], DS->Ele[ileft].zmin,DummyY[i_plus_totele],DS->Riv[i].zmin, DS->Riv[i].Dist[1], 0,0,0,0,RIV_ELESUB);
                        nabrAqDepth=(DS->Ele[ileft].zmax-DS->Ele[ileft].zmin);                        
			avgKH(0,0,0,0,0,0,DS->Riv[i].KsatH,DS->Ele[ileft].Macropore,DummyY[ileft_plus_2ele],nabrAqDepth,DS->Ele[ileft].macD,DS->Ele[ileft].macKsatH,DS->Ele[ileft].vAreaF,DS->Ele[ileft].KsatH);
                        DS->FluxRiv[i][4]=DS->Riv[i].Length*Avg_Ksat*Grad_Y*Avg_Y;
                        /***********************************************************************************/
                        /* Lateral Flux between rectangular element (beneath river) and triangular element */
                        /***********************************************************************************/                        
			GradCalc(DummyY[ileft_plus_2ele], DS->Ele[ileft].zmin,DummyY[i_plus_totele1riv],DS->Ele[i_plus_ele].zmin, DS->Riv[i].Dist[1], 0,0,DS->Riv[i].zmin,0,RIVSUB_ELESUB);                        
			avgKH(DS->Ele[ileft].Macropore,DummyY[ileft_plus_2ele],nabrAqDepth,DS->Ele[ileft].macD,DS->Ele[ileft].macKsatH,DS->Ele[ileft].vAreaF,DS->Ele[ileft].KsatH,DS->Ele[iright].Macropore,DummyY[iright_plus_2ele],DS->Ele[iright].zmax-DS->Ele[iright].zmin,DS->Ele[iright].macD,DS->Ele[iright].macKsatH,DS->Ele[iright].vAreaF,DS->Ele[iright].KsatH);                        
			avgKH(0,0,0,0,0,0,Avg_Ksat,DS->Ele[ileft].Macropore,DummyY[ileft_plus_2ele],nabrAqDepth,DS->Ele[ileft].macD,DS->Ele[ileft].macKsatH,DS->Ele[ileft].vAreaF,DS->Ele[ileft].KsatH);
                        DS->FluxRiv[i][7]=DS->Riv[i].Length*Avg_Ksat*Grad_Y*Avg_Y;
                        /* replace flux term */
                        DS->FluxSub[ileft][DS->Riv[i].lrEdge[0]] =  -DS->FluxRiv[i][4];
                        DS->FluxSub[ileft][DS->Riv[i].lrEdge[0]] = DS->FluxSub[ileft][DS->Riv[i].lrEdge[0]] -DS->FluxRiv[i][7];
		#endif
                        }
                if (DS->Riv[i].RightEle > 0)
                        {
                        /*****************************************************************************/
                        /* Lateral Surface Flux Calculation between River-Triangular element Follows */
                        /*****************************************************************************/ 
                       	OLFeleToriv(DummyY[iright]+DS->Ele[iright].zmax,DS->Ele[iright].zmax,DS->Riv_Mat[DS->Riv[i].material-1].Cwr, DS->Riv[i].zmax,TotalY,DS->FluxRiv,i,3,DS->Riv[i].Length);
                        DS->FluxSurf[iright][DS->Riv[i].lrEdge[1]] = -DS->FluxRiv[i][3];                        
		#ifdef SUB_SURF_RIV
                        /*********************************************************************************/
                        /* Lateral Sub-surface Flux Calculation between River-Triangular element Follows */
                        /*********************************************************************************/
                        GradCalc(DummyY[iright_plus_2ele], DS->Ele[iright].zmin,DummyY[i_plus_totele],DS->Riv[i].zmin, DS->Riv[i].Dist[2], 0,0,0,0,RIV_ELESUB);
                        /* take care of macropore effect */
                        nabrAqDepth=(DS->Ele[iright].zmax-DS->Ele[iright].zmin);                        
			avgKH(0,0,0,0,0,0,DS->Riv[i].KsatH,DS->Ele[iright].Macropore,DummyY[iright_plus_2ele],nabrAqDepth,DS->Ele[iright].macD,DS->Ele[iright].macKsatH,DS->Ele[iright].vAreaF,DS->Ele[iright].KsatH);                        
			DS->FluxRiv[i][5]=DS->Riv[i].Length*Avg_Ksat*Grad_Y*Avg_Y;
                        /***********************************************************************************/
                        /* Lateral Flux between rectangular element (beneath river) and triangular element */
                        /***********************************************************************************/                        
			GradCalc(DummyY[iright_plus_2ele], DS->Ele[iright].zmin,DummyY[i_plus_totele1riv],DS->Ele[i_plus_ele].zmin, DS->Riv[i].Dist[2],0,0,DS->Riv[i].zmin,0,RIVSUB_ELESUB);
                        /* Note this is same as that was calculated for left ele */                        
			avgKH(DS->Ele[ileft].Macropore,DummyY[ileft_plus_2ele],DS->Ele[ileft].zmax-DS->Ele[ileft].zmin,DS->Ele[ileft].macD,DS->Ele[ileft].macKsatH,DS->Ele[ileft].vAreaF,DS->Ele[ileft].KsatH,DS->Ele[iright].Macropore,DummyY[iright_plus_2ele],nabrAqDepth,DS->Ele[iright].macD,DS->Ele[iright].macKsatH,DS->Ele[iright].vAreaF,DS->Ele[iright].KsatH);                        
			avgKH(0,0,0,0,0,0,Avg_Ksat,DS->Ele[iright].Macropore,DummyY[iright_plus_2ele],nabrAqDepth,DS->Ele[iright].macD,DS->Ele[iright].macKsatH,DS->Ele[iright].vAreaF,DS->Ele[iright].KsatH);
                        DS->FluxRiv[i][8]=DS->Riv[i].Length*Avg_Ksat*Grad_Y*Avg_Y;                        
			/* replace flux item */
			DS->FluxSub[iright][DS->Riv[i].lrEdge[1]] =  -DS->FluxRiv[i][5];                        
			DS->FluxSub[iright][DS->Riv[i].lrEdge[1]] = DS->FluxSub[iright][DS->Riv[i].lrEdge[1]] -DS->FluxRiv[i][8];
		#endif
                        }                
	#ifdef SUB_SURF_RIV
		Avg_Wid=CS_AreaOrPerem(DS->Riv_Shape[DS->Riv[i].shape - 1].interpOrd,DummyY[i_plus_totele],DS->Riv[i].coeff,3);                
		Dif_Y=(DS->Riv[i].zmin-(DummyY[i_plus_totele1riv]+DS->Ele[i_plus_ele].zmin))>0?DummyY[i_plus_totele]:DummyY[i_plus_totele]+DS->Riv[i].zmin-(DummyY[i_plus_totele1riv]+DS->Ele[i_plus_ele].zmin);
        	Grad_Y=Dif_Y/DS->Riv[i].bedThick;
        	DS->FluxRiv[i][6]=DS->Riv[i].KsatV*Avg_Wid*DS->Riv[i].Length*Grad_Y;
	#elif SURF_RIV
		DS->FluxRiv[i][6]=0;
	#endif

	}


void fluxCalc_Riv_light(Model_Data DS, int i,int j, realtype t,realtype *DummyY)		
//xchen_20150405 flux calculation between River-Triangular element light version
	{
		int i_plus_ele,i_plus_totele,i_plus_totele1riv,ileft,iright,ileft_plus_2ele,iright_plus_2ele, inabr,inabr_left,inabr_right;
		realtype tmpVar,Wid, Wid_down, Avg_Wid, Avg_Rough, Avg_Sf, Perem, Avg_Perem, nabrAqDepth;
//		i_plus_ele=i+DS->NumEle;
                i_plus_totele=i+DS->totele;
                i_plus_totele1riv=i+DS->totele+DS->NumRiv;
                ileft=DS->Riv[i].LeftEle-1;
                iright=DS->Riv[i].RightEle-1;
                ileft_plus_2ele=ileft+2*DS->NumEle;
                iright_plus_2ele=iright+2*DS->NumEle;
                TotalY = DummyY[i_plus_totele] + DS->Riv[i].zmin;
                if(ileft == j)
                        {
                        /*****************************************************************************/
                        /* Lateral Surface Flux Calculation between River-Triangular element Follows */
                        /*****************************************************************************/
                        OLFeleToriv(DummyY[ileft]+DS->Ele[ileft].zmax,DS->Ele[ileft].zmax,DS->Riv_Mat[DS->Riv[i].material-1].Cwr, DS->Riv[i].zmax,TotalY,DS->FluxRiv,i,2,DS->Riv[i].Length);
                        DS->FluxSurf[ileft][DS->Riv[i].lrEdge[0]] = -DS->FluxRiv[i][2];
		#ifdef SUB_SURF_RIV
                        /*********************************************************************************/
                        /* Lateral Sub-surface Flux Calculation between River-Triangular element Follows */
                        /*********************************************************************************/
                        GradCalc(DummyY[ileft_plus_2ele], DS->Ele[ileft].zmin,DummyY[i_plus_totele],DS->Riv[i].zmin, DS->Riv[i].Dist[1], 0,0,0,0,RIV_ELESUB);
                        nabrAqDepth=(DS->Ele[ileft].zmax-DS->Ele[ileft].zmin);                        
			avgKH(0,0,0,0,0,0,DS->Riv[i].KsatH,DS->Ele[ileft].Macropore,DummyY[ileft_plus_2ele],nabrAqDepth,DS->Ele[ileft].macD,DS->Ele[ileft].macKsatH,DS->Ele[ileft].vAreaF,DS->Ele[ileft].KsatH);
                        DS->FluxRiv[i][4]=DS->Riv[i].Length*Avg_Ksat*Grad_Y*Avg_Y;
                        /***********************************************************************************/
                        /* Lateral Flux between rectangular element (beneath river) and triangular element */
                        /***********************************************************************************/                        
			GradCalc(DummyY[ileft_plus_2ele], DS->Ele[ileft].zmin,DummyY[i_plus_totele1riv],DS->Ele[ileft].zmin, DS->Riv[i].Dist[1], 0,0,DS->Riv[i].zmin,0,RIVSUB_ELESUB);   										//xchen_20150405                     
//			GradCalc(DummyY[ileft_plus_2ele], DS->Ele[ileft].zmin,DummyY[i_plus_totele1riv],DS->Ele[i_plus_ele].zmin, DS->Riv[i].Dist[1], 0,0,DS->Riv[i].zmin,0,RIVSUB_ELESUB);                        
			avgKH(DS->Ele[ileft].Macropore,DummyY[ileft_plus_2ele],nabrAqDepth,DS->Ele[ileft].macD,DS->Ele[ileft].macKsatH,DS->Ele[ileft].vAreaF,DS->Ele[ileft].KsatH,DS->Ele[iright].Macropore,DummyY[iright_plus_2ele],DS->Ele[iright].zmax-DS->Ele[iright].zmin,DS->Ele[iright].macD,DS->Ele[iright].macKsatH,DS->Ele[iright].vAreaF,DS->Ele[iright].KsatH);                        
			avgKH(0,0,0,0,0,0,Avg_Ksat,DS->Ele[ileft].Macropore,DummyY[ileft_plus_2ele],nabrAqDepth,DS->Ele[ileft].macD,DS->Ele[ileft].macKsatH,DS->Ele[ileft].vAreaF,DS->Ele[ileft].KsatH);
                        DS->FluxRiv[i][7]=DS->Riv[i].Length*Avg_Ksat*Grad_Y*Avg_Y;
                        /* replace flux term */
                        DS->FluxSub[ileft][DS->Riv[i].lrEdge[0]] =  -DS->FluxRiv[i][4];
                        DS->FluxSub[ileft][DS->Riv[i].lrEdge[0]] = DS->FluxSub[ileft][DS->Riv[i].lrEdge[0]] -DS->FluxRiv[i][7];
		#endif
                        }
                if (iright == j)
                        {
                        /*****************************************************************************/
                        /* Lateral Surface Flux Calculation between River-Triangular element Follows */
                        /*****************************************************************************/ 
                       	OLFeleToriv(DummyY[iright]+DS->Ele[iright].zmax,DS->Ele[iright].zmax,DS->Riv_Mat[DS->Riv[i].material-1].Cwr, DS->Riv[i].zmax,TotalY,DS->FluxRiv,i,3,DS->Riv[i].Length);
                        DS->FluxSurf[iright][DS->Riv[i].lrEdge[1]] = -DS->FluxRiv[i][3];                        
		#ifdef SUB_SURF_RIV
                        /*********************************************************************************/
                        /* Lateral Sub-surface Flux Calculation between River-Triangular element Follows */
                        /*********************************************************************************/
                        GradCalc(DummyY[iright_plus_2ele], DS->Ele[iright].zmin,DummyY[i_plus_totele],DS->Riv[i].zmin, DS->Riv[i].Dist[2], 0,0,0,0,RIV_ELESUB);
                        /* take care of macropore effect */
                        nabrAqDepth=(DS->Ele[iright].zmax-DS->Ele[iright].zmin);                        
			avgKH(0,0,0,0,0,0,DS->Riv[i].KsatH,DS->Ele[iright].Macropore,DummyY[iright_plus_2ele],nabrAqDepth,DS->Ele[iright].macD,DS->Ele[iright].macKsatH,DS->Ele[iright].vAreaF,DS->Ele[iright].KsatH);                        
			DS->FluxRiv[i][5]=DS->Riv[i].Length*Avg_Ksat*Grad_Y*Avg_Y;
                        /***********************************************************************************/
                        /* Lateral Flux between rectangular element (beneath river) and triangular element */
                        /***********************************************************************************/                        
			GradCalc(DummyY[iright_plus_2ele], DS->Ele[iright].zmin,DummyY[i_plus_totele1riv],DS->Ele[iright].zmin, DS->Riv[i].Dist[2],0,0,DS->Riv[i].zmin,0,RIVSUB_ELESUB);										//xchen_20150405
//			GradCalc(DummyY[iright_plus_2ele], DS->Ele[iright].zmin,DummyY[i_plus_totele1riv],DS->Ele[i_plus_ele].zmin, DS->Riv[i].Dist[2],0,0,DS->Riv[i].zmin,0,RIVSUB_ELESUB);
                        /* Note this is same as that was calculated for left ele */                        
			avgKH(DS->Ele[ileft].Macropore,DummyY[ileft_plus_2ele],DS->Ele[ileft].zmax-DS->Ele[ileft].zmin,DS->Ele[ileft].macD,DS->Ele[ileft].macKsatH,DS->Ele[ileft].vAreaF,DS->Ele[ileft].KsatH,DS->Ele[iright].Macropore,DummyY[iright_plus_2ele],nabrAqDepth,DS->Ele[iright].macD,DS->Ele[iright].macKsatH,DS->Ele[iright].vAreaF,DS->Ele[iright].KsatH);                        
			avgKH(0,0,0,0,0,0,Avg_Ksat,DS->Ele[iright].Macropore,DummyY[iright_plus_2ele],nabrAqDepth,DS->Ele[iright].macD,DS->Ele[iright].macKsatH,DS->Ele[iright].vAreaF,DS->Ele[iright].KsatH);
                        DS->FluxRiv[i][8]=DS->Riv[i].Length*Avg_Ksat*Grad_Y*Avg_Y;                        
			/* replace flux item */
			DS->FluxSub[iright][DS->Riv[i].lrEdge[1]] =  -DS->FluxRiv[i][5];                        
			DS->FluxSub[iright][DS->Riv[i].lrEdge[1]] = DS->FluxSub[iright][DS->Riv[i].lrEdge[1]] -DS->FluxRiv[i][8];
		#endif
                        }                
	#ifdef SUB_SURF_RIV
		Avg_Wid=CS_AreaOrPerem(DS->Riv_Shape[DS->Riv[i].shape - 1].interpOrd,DummyY[i_plus_totele],DS->Riv[i].coeff,3);                
		Dif_Y=(DS->Riv[i].zmin-(DummyY[i_plus_totele1riv]+DS->Ele[i_plus_totele1riv].zmin))>0?DummyY[i_plus_totele]:DummyY[i_plus_totele]+DS->Riv[i].zmin-(DummyY[i_plus_totele1riv]+DS->Ele[i_plus_totele1riv].zmin);						//xchen_20150405
//		Dif_Y=(DS->Riv[i].zmin-(DummyY[i_plus_totele1riv]+DS->Ele[i_plus_ele].zmin))>0?DummyY[i_plus_totele]:DummyY[i_plus_totele]+DS->Riv[i].zmin-(DummyY[i_plus_totele1riv]+DS->Ele[i_plus_ele].zmin);
        	Grad_Y=Dif_Y/DS->Riv[i].bedThick;
        	DS->FluxRiv[i][6]=DS->Riv[i].KsatV*Avg_Wid*DS->Riv[i].Length*Grad_Y;
	#elif SURF_RIV
		DS->FluxRiv[i][6]=0;
	#endif
	}
