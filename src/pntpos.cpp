/*------------------------------------------------------------------------------
* pntpos.c : standard positioning
*
*          Copyright (C) 2007-2018 by T.TAKASU, All rights reserved.
*
* version : $Revision:$ $Date:$
* history : 2010/07/28 1.0  moved from rtkcmn.c
*                           changed api:
*                               pntpos()
*                           deleted api:
*                               pntvel()
*           2011/01/12 1.1  add option to include unhealthy satellite
*                           reject duplicated observation data
*                           changed api: ionocorr()
*           2011/11/08 1.2  enable snr mask for single-mode (rtklib_2.4.1_p3)
*           2012/12/25 1.3  add variable snr mask
*           2014/05/26 1.4  support galileo and beidou
*           2015/03/19 1.5  fix bug on ionosphere correction for GLO and BDS
*           2018/10/10 1.6  support api change of satexclude()
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants -----------------------------------------------------------------*/

#define SQR(x)      ((x)*(x))

#define NX          (4+3)       /* # of estimated parameters */

#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay std (m) */
#define ERR_TROP    3.0         /* tropspheric delay std (m) */
#define ERR_SAAS    0.3         /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5         /* broadcast iono model error factor */
#define ERR_CBIAS   0.3         /* code bias error std (m) */
#define REL_HUMI    0.7         /* relative humidity for saastamoinen model */


extern "C" {
	extern double ionmodel_nequick(double * galionpar, int month, double UTC, double usr_longitude_degree, double usr_latitude_degree, double usr_height_meters,
		double sat_longitude_degree, double sat_latitude_degree, double sat_height_meters);
}

/* pseudorange measurement error variance ------------------------------------*/
static double varerr(const prcopt_t *opt, double el, int sys)
{
    double fact,varr;
    fact=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);
    varr=SQR(opt->err[0])*(SQR(opt->err[1])+SQR(opt->err[2])/sin(el));
    if (opt->ionoopt==IONOOPT_IFLC) varr*=SQR(3.0); /* iono-free */
    return SQR(fact)*varr;
}

/* psendorange with code bias correction -------------------------------------*/
static double prange(const obsd_t *obs, const nav_t *nav, const double *azel,
	int iter, const prcopt_t *opt, double *var, double *tgd1, double *tgd2)
{
	static double gammaGPSL2 = SQR(FREQ1 / FREQ2), gammaGPSL5 = SQR(FREQ1 / FREQ5);
    const double *lam=nav->lam[obs->sat-1];
	double PC, P1, P2, P1_P2, P1_C1, P2_C2, gamma, tgd, bgd_E5aE1, bgd_E5bE1;
    int i=0,j=1,f,sys,prn;
	int fidx[MAXFREQ] = { 0 };
	int galcode;
    *var=0.0;
	prcopt_t opt2 = *opt;
	opt2.ionoopt = IONOOPT_IFLC;

    if (!(sys=satsys(obs->sat,&prn))) return 0.0;
    
	frqidx(opt2, fidx);
    /* L1-L2 for GPS/GLO/QZS, L1-L5 for GAL/SBS */
   // if (NFREQ>=3&&(sys&(SYS_GAL|SYS_SBS))) j=2;

	if (opt->ionoopt == IONOOPT_IFLC) {
		if (fidx[0]<0 || fidx[1]<0) return 0.0;
	}
	else if (lam[fidx[0]] == 0.0) return 0.0;
    

    /* test snr mask */
    if (iter>0) {
		if (testsnr(0, fidx[0], azel[1], obs->SNR[fidx[0]] * 0.25, &opt->snrmask)) {
            trace(4,"snr mask: %s sat=%2d el=%.1f snr=%.1f\n",
				time_str(obs->time, 0), obs->sat, azel[1] * R2D, obs->SNR[fidx[0]] * 0.25);
            return 0.0;
        }
        if (opt->ionoopt==IONOOPT_IFLC) {
			if (testsnr(0, fidx[1], azel[1], obs->SNR[fidx[1]] * 0.25, &opt->snrmask)) return 0.0;
        }
    }
	gamma = opt->ionoopt == IONOOPT_IFLC ? SQR(lam[fidx[1]]) / SQR(lam[fidx[0]]) : 0.0; /* f1^2/f2^2 */
	P1 = obs->P[fidx[0]];
	P2 = opt->ionoopt == IONOOPT_IFLC ? obs->P[fidx[1]] : 0.0;
    P1_P2=nav->cbias[obs->sat-1][0];
    P1_C1=nav->cbias[obs->sat-1][1];
    P2_C2=nav->cbias[obs->sat-1][2];
    
	if (P1 == 0 || (P2 == 0 && opt->ionoopt == IONOOPT_IFLC))return 0.0;

	if (sys == SYS_GPS || sys == SYS_QZS){
		tgd = gettgd(obs->sat, nav, 0,NULL);
		if (obs->code[fidx[0]] == CODE_L1C) P1 += P1_C1; /* C1->P1 */
		if (obs->code[fidx[1]] == CODE_L2C) P2 += P2_C2; /* C2->P2 */
		P1 = P1 - tgd*(fidx[0] == 0 ? 1.0 : (fidx[0] == 1 ? gammaGPSL2 : (fidx[0] == 2 ? gammaGPSL5 : 0)));
		P2 = P2 - tgd*(fidx[1] == 0 ? 1.0 : (fidx[1] == 1 ? gammaGPSL2 : (fidx[1] == 2 ? gammaGPSL5 : 0)));
		if (tgd1 != NULL)*tgd1 = tgd / CLIGHT;
		if (tgd2 != NULL)*tgd2 = tgd*gammaGPSL2 / CLIGHT;
	}
	else if (sys == SYS_CMP){
		P1 = P1 - gettgd(obs->sat, nav, fidx[0], NULL);
		P2 = P2 - gettgd(obs->sat, nav, fidx[1], NULL);
		if (tgd1 != NULL)*tgd1 = gettgd(obs->sat, nav, fidx[0], NULL) / CLIGHT;;
		if (tgd2 != NULL)*tgd2 = gettgd(obs->sat, nav, fidx[1], NULL) / CLIGHT;;
	}
	else if (sys == SYS_GAL){
		//#define FREQ1       1.57542E9           /* L1/E1/B1C  frequency (Hz) */
		//#define FREQ2       1.22760E9           /* L2     frequency (Hz) */
		//#define FREQ5       1.17645E9           /* L5/E5a/B2a frequency (Hz) */
		//#define FREQ6       1.27875E9           /* E6/LEX frequency (Hz) */
		//#define FREQ7       1.20714E9           /* E5b/B2b    frequency (Hz) */
		//#define FREQ8       1.191795E9          /* E5a+b  frequency (Hz) */
		bgd_E5bE1 = gettgd(obs->sat, nav, 1, &galcode);
		bgd_E5aE1 = gettgd(obs->sat, nav, 2, &galcode);
		if (fidx[0] == 0 && galcode!=2){
			P1 = P1 - bgd_E5bE1;/* E1 from I/NAV*/
		}
		else if (fidx[0] == 0 && galcode == 2){
			P1 = P1 - bgd_E5aE1;/* E1 from F/NAV*/
		}

		if (fidx[0] == 1){
			P1 = P1 - bgd_E5bE1*SQR(FREQ1/FREQ7);/* E5b*/
		}

		if (fidx[0] == 2){
			P1 = P1 -bgd_E5aE1*SQR(FREQ1 / FREQ5);/* E5a*/
		}


		if (fidx[1] == 0){
			/* no data */
		}
		else if (fidx[1] == 1){
			P2 = P2 - bgd_E5bE1*SQR(FREQ1 / FREQ7);/* E5b*/
		}
		else if (fidx[1] == 2){
			P2 = P2 - bgd_E5aE1*SQR(FREQ1 / FREQ5);/* E5a*/
		}
		if (tgd1 != NULL)*tgd1 = gettgd(obs->sat, nav, fidx[0], NULL) / CLIGHT;;
		if (tgd2 != NULL)*tgd2 = gettgd(obs->sat, nav, fidx[1], NULL) / CLIGHT;;
	}

	if (opt->ionoopt == IONOOPT_IFLC) { /* dual-frequency */
		/* iono-free combination */
		PC = (gamma*P1 - P2) / (gamma - 1.0);
	}
	else { /* single-frequency */
		PC = P1;
	}

	if (opt->sateph == EPHOPT_SBAS) PC -= P1_C1; /* sbas clock based C1 */

	*var = SQR(ERR_CBIAS);

	return PC;
}

static int selion(char *modelname, gtime_t time, bds_ion_t bds_ion)
{/* selecte specified ion according experience */
	/* GEO first; IGSO second; MEO last */
	int i, j = -1, iGEO = -1, iIGSO = -1, iMEO = -1;
	double ep[6], hour;
	double t_GEO = 24.0, t_IGSO = 24.0, t_MEO = 24.0, dt;
	double tmin = 24.0;
	time2epoch(time, ep);
	hour = ep[3];
	if (strcmp(modelname, "BDSK8") == 0){
		for (i = 0; i < bds_ion.nk8; i++){
			dt = fabs(bds_ion.bdsk8[i].hour - hour);
			if (dt <= tmin) {
				j = i;
				tmin = dt;
			}
		}

		if (j<0) {
			//trace(2, "no broadcast ion: %s \n", time_str(time, 0));
			return -1;
		}

		return j;
	}
	else if (strcmp(modelname, "BDSSH9") == 0){
		for (i = 0; i < bds_ion.nsh9; i++){
			dt = fabs(bds_ion.bdssh9[i].hour - hour);
			if (dt <= tmin) {
				j = i;
				tmin = dt;
			}
		}
		if (j<0) {
			//trace(2, "no broadcast ion: %s \n", time_str(time, 0));
			return -1;
		}

		return j;
	}
}

/* ionospheric correction ------------------------------------------------------
* compute ionospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          int    sat       I   satellite number
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    ionoopt   I   ionospheric correction option (IONOOPT_???)
*          double *ion      O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric delay (L1) variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int ionocorr(const prcopt_t opt, gtime_t time, const nav_t *nav, int sat, const double *pos, const double *satpos,
                    const double *azel, int ionoopt, double *ion, double *var)
{
    ////trace(4,"ionocorr: time=%s opt=%d sat=%2d pos=%.3f %.3f azel=%.3f %.3f\n",
    //      time_str(time,3),ionoopt,sat,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
    //      azel[1]*R2D);
	int stat, i, idx;
	double tow, k;
	gtime_t bdt;
	double ep[6];
	double ionk8[8] = { 0.0 };
	double ionsh9[9] = { 0.0 };
	prcopt_t opt2 = opt;
	opt2.ionoopt = IONOOPT_IFLC;
	int fidx[MAXFREQ] = { 0 };
	double freq = 0;
	frqidx(opt2, fidx);

	freq = nav->lam[sat - 1][fidx[0]] > 0 ? CLIGHT / nav->lam[sat - 1][fidx[0]] : FREQ1_CMP;

	k = FREQ1_CMP*FREQ1_CMP/freq/freq;
    /* broadcast model */
    if (ionoopt==IONOOPT_BRDC) {
        *ion=k*ionmodel(time,nav->ion_gps,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
	
	/* BDSSH9 first, BDSK8 2nd */
	if (ionoopt == IONOOPT_BDSK8){
		bdt = gpst2bdt(time);
		idx = selion("BDSK8", time, nav->ion_bdsk9->bds_ion);
		if (idx >= 0) for (i = 0; i < 8; i++) ionk8[i] = nav->ion_bdsk9->bds_ion.bdsk8[idx].ion[i];
		*ion = k*ionmodel_BDSK8(bdt, ionk8, pos, azel);
		//if (fabs(*ion)>10.0)printf("sat=%d ,ion=%lf\n", sat, *ion);
		*var = SQR(*ion*ERR_BRDCI);
		return 1;
	}
	if (ionoopt == IONOOPT_BDSSH9){
		bdt = gpst2bdt(time);
		//idx = selion("BDSSH9", time, nav->ion_bdsk9->bds_ion);
		//if (idx >= 0) for (i = 0; i < 9; i++) ionsh9[i] = nav->bds_ion.bdssh9[idx].ion[i];

		//double ep[6] = { 2020, 11, 5, 17, 0, 0 };
		//gtime_t time2 = epoch2time(ep);
		////time2epoch(time, ep);

		//double stapos[3] = { -2159945.3149, 4383237.1901, 4085412.3177 },pos2[3];
		//ecef2pos(stapos, pos2);
		//double satpos2[3] = { -34342443.2642, 24456393.0586, 12949.4754 };
		//double ionsh2[9] = {
		//	13.82642069, -1.78004616, 5.17200000,   //BDS1
		//	3.13030940, -4.58961329, 0.32483867,    //BDS2
		//	-0.07802383, 1.24312590, 0.37763627     //BDS3
		//};
			
	    *ion = k*ionmodel_BDSK9(bdt, nav->ion_bdsk9, pos, satpos);
		//*ion = k*ionmodel_BDSK9_2(time2, ionsh2, pos2, satpos2);

		*var = SQR(*ion*ERR_BRDCI);

		return 1;
	}
	/* ionex tec model */
	if (ionoopt == IONOOPT_TEC) {
		*ion = 0.0;
		*var = 0.0;
		iontec(time, nav, pos, azel, 1, ion, var);
		*ion = k*(*ion);
		return 1;
	}
	if (ionoopt == IONOOPT_GALION){
		double ionpar[3],tec=0.0;
		double satepos[3];
		int i;
		double ep[6];


		/*transfer xyz to blh */
		ecef2pos(satpos, satepos);

		for (i = 0; i < 3; i++){
			satepos[i] = i<2?satepos[i]*R2D: satepos[i];
			ionpar[i] = nav->ion_gal[i];
		}
		time2epoch(time, ep);/*
#if 0
		double coef[3] = { 236.831641, -0.39362878, 0.00402826613 };
		for (i = 0; i < 3; i++) {
			ionpar[i] = coef[i];
		}
		ep[0] = 0; ep[1] = 4; ep[3] = 0.0; ep[4] = 0.0; ep[5] = 0.0;
		int month = 1;
		double UTC = 0.0;
		double usr_longitude_degree = 297.66;
		double usr_latitude_degree = 82.49;
		double usr_height_meters = 78.11;

		double sat_longitude_degree = 8.23;
		double sat_latitude_degree = 54.29;
		double sat_height_meters = 20281546.18;
		tec = ionmodel_nequick(ionpar, ep[1], ep[3] + ep[4] / 60.0 + ep[5] / 3600.0,
		usr_longitude_degree, usr_latitude_degree, usr_height_meters, sat_longitude_degree, sat_latitude_degree, sat_height_meters);
		printf("%lf \n", tec);
#endif*/

		tec = ionmodel_nequick(ionpar, ep[1], ep[3] + ep[4] / 60.0 + ep[5] / 3600.0,pos[1] * R2D, pos[0] * R2D, pos[2], satepos[1], satepos[0], satepos[2]);

		/*transfer to B1I */
		/* switch to BDS B1I */
		tec = tec*40.28e16 / FREQ1_CMP / FREQ1_CMP;
		tec = tec < 0 ? 0.0 : tec;

		//printf("%lf \n", tec);

		*ion = k*tec;
		*var = SQR(*ion*ERR_BRDCI);
		return 1;
	}
    *ion=0.0;
    *var=ionoopt==IONOOPT_OFF?SQR(ERR_ION):0.0;
    return 1;
}
/* tropospheric correction -----------------------------------------------------
* compute tropospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    tropopt   I   tropospheric correction option (TROPOPT_???)
*          double *trp      O   tropospheric delay (m)
*          double *var      O   tropospheric delay variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                    const double *azel, int tropopt, double *trp, double *var)
{
    //trace(4,"tropcorr: time=%s opt=%d pos=%.3f %.3f azel=%.3f %.3f\n",
    //      time_str(time,3),tropopt,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
    //      azel[1]*R2D);
    
    /* saastamoinen model */
    if (tropopt==TROPOPT_SAAS||tropopt==TROPOPT_EST||tropopt==TROPOPT_ESTG) {
        *trp=tropmodel(time,pos,azel,0.45);
        *var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
        return 1;
    }
    /* sbas troposphere model */
    if (tropopt==TROPOPT_SBAS) {
        *trp=sbstropcorr(time,pos,azel,var);
        return 1;
    }
    /* no correction */
    *trp=0.0;
    *var=tropopt==TROPOPT_OFF?SQR(ERR_TROP):0.0;
    return 1;
}
/* pseudorange residuals -----------------------------------------------------*/
static int rescode(int iter, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *vare, const int *svh,
                   const nav_t *nav, const double *x, const prcopt_t *opt,
                   double *v, double *H, double *var, double *azel, int *vsat,
				   double *resp, int *ns,int *sat)
{
    double r,dion,dtrp,vmeas,vion,vtrp,rr[3],pos[3],dtr,e[3],P,lam_L1;
    int i,j,nv=0,sys,mask[4]={0};
	char cprn[128];
	double res, tgd1, tgd2, dr;

    trace(3,"resprng : n=%d\n",n);
    
    for (i=0;i<3;i++) rr[i]=x[i];
    dtr=x[3];
	res=tgd1=tgd2=dr=0.0;

    ecef2pos(rr,pos);
    
	for (i = *ns = 0; i < n&&i < MAXOBS; i++) {
		vsat[i] = 0; azel[i * 2] = azel[1 + i * 2] = resp[i] = 0.0;
		if (!(sys = satsys(obs[i].sat, NULL))) continue;
		satno2id(obs[i].sat, cprn);



		/* reject duplicated observation data */
		if (i < n - 1 && i < MAXOBS - 1 && obs[i].sat == obs[i + 1].sat) {
			trace(2, "duplicated observation data %s sat=%2d\n",
				time_str(obs[i].time, 3), obs[i].sat);
			i++;
			continue;
		}

		/* geometric distance/azimuth/elevation angle */
		if ((r = geodist(sys, rs + i * 6, rr, e, &dr)) <= 0.0 ||
			satazel(pos, e, azel + i * 2) < opt->elmin)
		{
			//trace(4, "rescode: sat=%d (%s) removed! no satpos/obs or ele<%lf, elev=%lf\n", obs[i].sat, cprn, opt->elmin,satazel(pos, e, azel + i * 2)*R2D);
			continue;
		}
        
        /* psudorange with code bias correction */
        if ((P=prange(obs+i,nav,azel+i*2,iter,opt,&vmeas,&tgd1,&tgd2))==0.0) continue;
        

        /* excluded satellite? */
        if (satexclude(obs[i].sat,vare[i],svh[i],opt)) continue;
        
        /* ionospheric corrections */
		if (!ionocorr(*opt, obs[i].time, nav, obs[i].sat, pos, rs + i * 6, azel + i * 2,
                      iter>0?opt->ionoopt:IONOOPT_BRDC,&dion,&vion)) continue;
        
        /* tropospheric corrections */
        if (!tropcorr(obs[i].time,nav,pos,azel+i*2,
                      iter>0?opt->tropopt:TROPOPT_SAAS,&dtrp,&vtrp)) {
            continue;
        }
        /* pseudorange residual */
        v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);
        
        /* design matrix */
        for (j=0;j<NX;j++) H[j+nv*NX]=j<3?-e[j]:(j==3?1.0:0.0);
        
        /* time system and receiver bias offset correction */
        if      (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*NX]=1.0; mask[1]=1;}
        else if (sys==SYS_GAL) {v[nv]-=x[5]; H[5+nv*NX]=1.0; mask[2]=1;}
        else if (sys==SYS_CMP) {v[nv]-=x[6]; H[6+nv*NX]=1.0; mask[3]=1;}
        else mask[0]=1;

		sat[obs[i].sat - 1] = 1;
        vsat[i]=1; resp[i]=v[nv]; (*ns)++;
		res = v[nv];

        /* error variance */
		//var[nv++] = (varerr(opt, azel[1 + i * 2], sys) + vare[i] + vmeas + vion + vtrp)  *(exp(-azel[i * 2 + 1] / 15.0 / D2R));
		var[nv++] = (varerr(opt, azel[1 + i * 2], sys) + vare[i] + vmeas + vion + vtrp);
		satno2id(obs[i].sat,cprn);
		if(opt->mode==PMODE_SINGLE)trace(2, "rescode: sat=%2d (%s)azel=%5.1f %4.1f res=%8.2f r=%14.3f dtr=%8.2f dts=%14.3f ion=%7.3f dtrp=%7.3f d_E_Rt=%7.3f tgd1=%19.12E tgd2=%19.12e\n", obs[i].sat, cprn,
			azel[i * 2] * R2D, azel[1 + i * 2] * R2D, res, r, dtr, CLIGHT*dts[i * 2], dion, dtrp, dr, tgd1, tgd2);
    }
    /* constraint to avoid rank-deficient */
    for (i=0;i<4;i++) {
        if (mask[i]) continue;
        v[nv]=0.0;
        for (j=0;j<NX;j++) H[j+nv*NX]=j==i+3?1.0:0.0;
        var[nv++]=0.01;
    }
    return nv;
}
/* validate solution ---------------------------------------------------------*/
static int valsol(const double *azel, const int *vsat, int n,
                  const prcopt_t *opt, const double *v, int nv, int nx,
                  char *msg,double *dop0)
{
    double azels[MAXOBS*2],dop[4],vv;
    int i,ns;
    
    trace(3,"valsol  : n=%d nv=%d\n",n,nv);
    
	if (opt->coordfixed > 1E3) return 1;

	for (i = 0; i<6; i++) { dop0[i] = 0.0; }
    /* chi-square validation of residuals */
    vv=dot(v,v,nv);
    if (nv>nx&&vv>chisqr[nv-nx-1]) {
        sprintf(msg,"chi-square error nv=%d vv=%.1f cs=%.1f",nv,vv,chisqr[nv-nx-1]);
        return 0;
    }
    /* large gdop check */
    for (i=ns=0;i<n;i++) {
        if (!vsat[i]) continue;
        azels[  ns*2]=azel[  i*2];
        azels[1+ns*2]=azel[1+i*2];
        ns++;
    }
    dops(ns,azels,opt->elmin,dop);
    if (dop[0]<=0.0||dop[0]>opt->maxgdop) {
        sprintf(msg,"gdop error nv=%d gdop=%.1f",nv,dop[0]);
        return 0;
    }
	for (i = 0; i<4; i++) { dop0[i] = dop[i]; }
    return 1;
}
/* estimate receiver position ------------------------------------------------*/
static int estpos(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *vare, const int *svh, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel, int *vsat,
                  double *resp, char *msg)
{
    double x[NX]={0},dx[NX],Q[NX*NX],*v,*H,*var,sig;
    int i,j,k,info,stat,nv,ns;
	double presp[MAXSAT], alpha[MAXSAT] = { 0.0};
	double mean, std,extd;
    trace(3,"estpos  : n=%d\n",n);
    
    v=mat(n+4,1); H=mat(NX,n+4); var=mat(n+4,1);
    
	for (i = 0; i<3; i++) x[i] = sol->rr[i];
	if (opt->coordfixed != 0 && norm(opt->ru, 3) > 0.0)
	{
		for (i = 0; i<3; i++) x[i] = opt->ru[i];
		//printf("%lf %lf %lf\n", x[0],x[1],x[2]);
	}
    
	trace(4, "x(%2.2d)=%14.3f %14.3f %14.3f\n", 0, x[0], x[1], x[2]);
    for (i=0;i<MAXITR;i++) {
        
        /* pseudorange residuals */
        nv=rescode(i,obs,n,rs,dts,vare,svh,nav,x,opt,v,H,var,azel,vsat,resp,
			&ns, sol->sat);
        
		//printf("iter=%d\n", i + 1);
        if (nv<NX) {
            sprintf(msg,"lack of valid sats ns=%d",nv);
            break;
        }
		for (j = 0; j < nv-3; j++) {
			presp[j] = v[j];
		}
		//median_t(resp, nres); 
		mean=median3_t(presp, nv-3,&std);


        /* weight by variance */
        for (j=0;j<nv;j++) {
            sig=sqrt(var[j]);
			if (std!=0.0&&std < 5.0&&alpha[j] == 0.0)
			{
				alpha[j] = fabs(presp[j] / std)<3.0 ? 1.0 : fabs(presp[j] / std/2);
			}
			if (alpha[j] == 0.0||j>=nv-3)alpha[j] = 1.0;
			alpha[j] = 1.0;
			sig = sig*alpha[j];
			if (j < nv-3)trace(1, "v(%2.2d)=%14.3f, sigma0=%14.3f std=%14.3f sig=%14.3f\n", i+1, presp[j], sqrt(var[j]), std, sig);



            v[j]/=sig;
            for (k=0;k<NX;k++) H[k+j*NX]/=sig;
        }
        /* least square estimation */
		if ((info = lsq(H, v, NX, nv, dx, Q, opt->coordfixed))) {
            sprintf(msg,"lsq error info=%d",info);
            break;
        }
        for (j=0;j<NX;j++) x[j]+=dx[j];
        
		trace(4, "x(%2.2d)=%14.3f %14.3f %14.3f\n", i + 1, x[0], x[1], x[2]);

        if (norm(dx,NX)<1E-4) {
            sol->type=0;
            sol->time=timeadd(obs[0].time,-x[3]/CLIGHT);
            sol->dtr[0]=x[3]/CLIGHT; /* receiver clock bias (s) */
            sol->dtr[1]=x[4]/CLIGHT; /* glo-gps time offset (s) */
            sol->dtr[2]=x[5]/CLIGHT; /* gal-gps time offset (s) */
            sol->dtr[3]=x[6]/CLIGHT; /* bds-gps time offset (s) */
            for (j=0;j<6;j++) sol->rr[j]=j<3?x[j]:0.0;
            for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*NX];
            sol->qr[3]=(float)Q[1];    /* cov xy */
            sol->qr[4]=(float)Q[2+NX]; /* cov yz */
            sol->qr[5]=(float)Q[2];    /* cov zx */
            sol->ns=(unsigned char)ns;
            sol->age=sol->ratio=0.0;
            
            /* validate solution */
			if ((stat = valsol(azel, vsat, n, opt, v, nv, NX, msg, sol->dop))) {
                sol->stat=opt->sateph==EPHOPT_SBAS?SOLQ_SBAS:SOLQ_SINGLE;
            }
			free(v); free(H); free(var);
            
            return stat;
        }
    }
    if (i>=MAXITR) sprintf(msg,"iteration divergent i=%d",i);
    
	free(v); free(H); free(var);
    
    return 0;
}
/* raim fde (failure detection and exclution) -------------------------------*/
static int raim_fde(const obsd_t *obs, int n, const double *rs,
                    const double *dts, const double *vare, const int *svh,
                    const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                    double *azel, int *vsat, double *resp, char *msg)
{
    obsd_t *obs_e;
    sol_t sol_e={{0}};
    char tstr[32],name[16],msg_e[128];
    double *rs_e,*dts_e,*vare_e,*azel_e,*resp_e,rms_e,rms=100.0;
    int i,j,k,nvsat,stat=0,*svh_e,*vsat_e,sat=0;
    
    trace(3,"raim_fde: %s n=%2d\n",time_str(obs[0].time,0),n);
    
    if (!(obs_e=(obsd_t *)malloc(sizeof(obsd_t)*n))) return 0;
    rs_e = mat(6,n); dts_e = mat(2,n); vare_e=mat(1,n); azel_e=zeros(2,n);
    svh_e=imat(1,n); vsat_e=imat(1,n); resp_e=mat(1,n); 
    
    for (i=0;i<n;i++) {
        
        /* satellite exclution */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            obs_e[k]=obs[j];
            matcpy(rs_e +6*k,rs +6*j,6,1);
            matcpy(dts_e+2*k,dts+2*j,2,1);
            vare_e[k]=vare[j];
            svh_e[k++]=svh[j];
        }
        /* estimate receiver position without a satellite */
        if (!estpos(obs_e,n-1,rs_e,dts_e,vare_e,svh_e,nav,opt,&sol_e,azel_e,
                    vsat_e,resp_e,msg_e)) {
            trace(3,"raim_fde: exsat=%2d (%s)\n",obs[i].sat,msg);
            continue;
        }
        for (j=nvsat=0,rms_e=0.0;j<n-1;j++) {
            if (!vsat_e[j]) continue;
            rms_e+=SQR(resp_e[j]);
            nvsat++;
        }
        if (nvsat<5) {
            trace(3,"raim_fde: exsat=%2d lack of satellites nvsat=%2d\n",
                  obs[i].sat,nvsat);
            continue;
        }
        rms_e=sqrt(rms_e/nvsat);
        
        trace(3,"raim_fde: exsat=%2d rms=%8.3f\n",obs[i].sat,rms_e);
        
        if (rms_e>rms) continue;
        
        /* save result */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            matcpy(azel+2*j,azel_e+2*k,2,1);
            vsat[j]=vsat_e[k];
            resp[j]=resp_e[k++];
        }
        stat=1;
        *sol=sol_e;
        sat=obs[i].sat;
        rms=rms_e;
        vsat[i]=0;
        strcpy(msg,msg_e);
    }
    if (stat) {
        time2str(obs[0].time,tstr,2); satno2id(sat,name);
        trace(2,"%s: %s excluded by raim\n",tstr+11,name);
    }
    free(obs_e);
    free(rs_e ); free(dts_e ); free(vare_e); free(azel_e);
    free(svh_e); free(vsat_e); free(resp_e);
    return stat;
}
/* doppler residuals ---------------------------------------------------------*/
static int resdop(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const nav_t *nav, const double *rr, const double *x,
                  const double *azel, const int *vsat, double *v, double *H)
{
    double lam,rate,pos[3],E[9],a[3],e[3],vs[3],cosel;
    int i,j,nv=0;
	double doppler;
    trace(3,"resdop  : n=%d\n",n);
    
    ecef2pos(rr,pos); xyz2enu(pos,E);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        
        lam=nav->lam[obs[i].sat-1][0];
        
        if (obs[i].D[0]==0.0||lam==0.0||!vsat[i]||norm(rs+3+i*6,3)<=0.0) {
            continue;
        }
        /* line-of-sight vector in ecef */
        cosel=cos(azel[1+i*2]);
        a[0]=sin(azel[i*2])*cosel;
        a[1]=cos(azel[i*2])*cosel;
        a[2]=sin(azel[1+i*2]);
        matmul("TN",3,1,3,1.0,E,a,0.0,e);
        
        /* satellite velocity relative to receiver in ecef */
        for (j=0;j<3;j++) vs[j]=rs[j+3+i*6]-x[j];
        
        /* range rate with earth rotation correction */
        rate=dot(vs,e,3)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
                                      rs[3+i*6]*rr[1]-rs[  i*6]*x[1]);
        
        /* doppler residual */
		/* igmas special */
		if (nav->igmasta==1){ doppler = obs[i].D[0]; }
		else{ doppler = obs[i].D[0]; }
		v[nv] = -lam*doppler - (rate + x[3] - CLIGHT*dts[1 + i * 2]);
		trace(5, "resdopper: %d dopv=%lf\n", obs[i].sat, v[nv]);
        /* design matrix */
        for (j=0;j<4;j++) H[j+nv*4]=j<3?-e[j]:1.0;
        
        nv++;
    }
    return nv;
}
/* estimate receiver velocity ------------------------------------------------*/
static void estvel(const obsd_t *obs, int n, const double *rs, const double *dts,
                   const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                   const double *azel, const int *vsat)
{
    double x[4]={0},dx[4],Q[16],*v,*H;
    int i,j,nv;
    
    trace(3,"estvel  : n=%d\n",n);
    
    v=mat(n,1); H=mat(4,n);
    
    for (i=0;i<MAXITR;i++) {
        
        /* doppler residuals */
        if ((nv=resdop(obs,n,rs,dts,nav,sol->rr,x,azel,vsat,v,H))<4) {
            break;
        }
        /* least square estimation */
        if (lsq(H,v,4,nv,dx,Q,0.0)) break;
        
        for (j=0;j<4;j++) x[j]+=dx[j];
        
        if (norm(dx,4)<1E-6) {
            for (i=0;i<3;i++) sol->rr[i+3]=x[i];
            break;
        }
    }
    free(v); free(H);
}
/* single-point positioning ----------------------------------------------------
* compute receiver position, velocity, clock bias by single-point positioning
* with pseudorange and doppler observables
* args   : obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          prcopt_t *opt    I   processing options
*          sol_t  *sol      IO  solution
*          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output)
*          ssat_t *ssat     IO  satellite status              (NULL: no output)
*          char   *msg      O   error message for error exit
* return : status(1:ok,0:error)
* notes  : assuming sbas-gps, galileo-gps, qzss-gps, compass-gps time offset and
*          receiver bias are negligible (only involving glonass-gps time offset
*          and receiver bias)
*-----------------------------------------------------------------------------*/
extern int pntpos(const obsd_t *obs, int n, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel, ssat_t *ssat,
                  char *msg)
{
    prcopt_t opt_=*opt;
    double *rs,*dts,*var,*azel_,*resp;
    int i,stat,vsat[MAXOBS]={0},svh[MAXOBS];
	double rescode[MAXSAT];
	int nres = 0;
    sol->stat=SOLQ_NONE;
    
    if (n<=0) {strcpy(msg,"no observation data"); return 0;}
    
    trace(3,"pntpos  : tobs=%s n=%d\n",time_str(obs[0].time,3),n);
    
    sol->time=obs[0].time; msg[0]='\0';
    
    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel_=zeros(2,n); resp=mat(1,n);
    
    if (opt_.mode!=PMODE_SINGLE) { /* for precise positioning */
#if 0
        opt_.sateph =EPHOPT_BRDC;
#endif
        opt_.ionoopt=IONOOPT_BRDC;
        opt_.tropopt=TROPOPT_SAAS;
    }
    /* satellite positons, velocities and clocks */
    satposs(sol->time,obs,n,nav,&opt_,opt_.sateph,rs,dts,var,svh);
    
    /* estimate receiver position with pseudorange */
    stat=estpos(obs,n,rs,dts,var,svh,nav,&opt_,sol,azel_,vsat,resp,msg);
    
    /* raim fde */
    if (!stat&&n>=6&&opt->posopt[4]) {
        //stat=raim_fde(obs,n,rs,dts,var,svh,nav,&opt_,sol,azel_,vsat,resp,msg);
    }
    /* estimate receiver velocity with doppler */
    if (stat) estvel(obs,n,rs,dts,nav,&opt_,sol,azel_,vsat);
    
    if (azel) {
        for (i=0;i<n*2;i++) azel[i]=azel_[i];
    }
    if (ssat) {
        for (i=0;i<MAXSAT;i++) {
            ssat[i].vs=0;
            ssat[i].azel[0]=ssat[i].azel[1]=0.0;
            ssat[i].resp[0]=ssat[i].resc[0]=0.0;
            ssat[i].snr[0]=0;
        }
		nres = 0;
		for (i = 0; i<n; i++) {
			//if (!vsat[i]) continue;
			ssat[obs[i].sat - 1].resp[0] = resp[i];
			rescode[nres++] = resp[i];
		}
		//median_t(resp, nres);
		if (isepoch(obs[0].time, "2022/08/20 00:01:44")){
			//int ok;
			//ok = 1;
		}
		median2_t(resp, nres);
        for (i=0;i<n;i++) {
            ssat[obs[i].sat-1].azel[0]=azel_[  i*2];
            ssat[obs[i].sat-1].azel[1]=azel_[1+i*2];
            ssat[obs[i].sat-1].snr[0]=obs[i].SNR[0];
			ssat[obs[i].sat - 1].resp[0] = resp[i];
            if (!vsat[i]) continue;
            ssat[obs[i].sat-1].vs=1;
        }
    }
    free(rs); free(dts); free(var); free(azel_); free(resp);
    return stat;
}
