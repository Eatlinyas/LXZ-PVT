/****************************************************************************************************************
 *   bdgim.cpp : BDGIM model function definition.
 *   BDS Broadcast Iononspheric model(BDGIM) calculation module
 *   Author: Zishen LI, Ningbo WANG
 *   e-mail: lizishen@aircas.ac.cn, wangnb@aircas.ac.cn
 *
 *   Copyright (C) 2021 by Aerospace Information Research Institute, Chinese Academy of Sciences, Beijing, China.
 *   This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
 *   Public License as published by the Free Software Foundation; either version 2 of the License, or (at your
 *   option) any later version.
 *
 *		Create at  2021,	Apr. 23
 ***************************************************************************************************************/
#define _CRT_SECURE_NO_WARNINGS
#include "bdgim.h"

/*********** define some basic consts **************************/
const double LAT_POLE = 80.27;                // latitude of north geomagnetic polar [unit:degree]
const double LON_POLE = -72.58;               // east longitude of north geomagnetic polar [unit:degree]
const double Hion_bdgim = 400000.0;           // heigth of Ionospheric layer [unit:m]
const double EARTH_RADIUS = 6378137.0;     	  // the average radius of earth,for compute the IPP [unit:m]

static double Init_mjd = 0.0;                 // initial mjd every day for determining whether recalculate non-broadcast parameters

/**** BDGIM Periodic Table for Non-Broadcast Coefficient Forecast. Corresponds to degree/order [ 3/0 3/1 3/-1 3/2 ... 5/2 5/-2 ] ****/
const double NonBrdPara_table[NONBRDNUM][TRISERINUM] = {
	{-0.610000,-0.510000, 0.230000,-0.060000, 0.020000, 0.010000, 0.000000,-0.010000,-0.000000, 0.000000, 0.010000,-0.190000,-0.090000,-0.180000, 0.150000, 1.090000, 0.500000,-0.340000, 0.000000,-0.130000, 0.050000,-0.060000, 0.030000,-0.030000, 0.040000},
	{-1.310000,-0.430000,-0.200000,-0.050000,-0.080000,-0.030000,-0.020000, 0.000000,-0.020000, 0.000000, 0.000000,-0.020000, 0.070000, 0.060000,-0.310000,-0.140000,-0.080000,-0.090000,-0.110000, 0.070000, 0.030000, 0.130000,-0.020000, 0.080000,-0.020000},
	{-2.000000, 0.340000,-0.310000, 0.060000,-0.060000, 0.010000,-0.030000, 0.010000, 0.010000, 0.030000, 0.000000, 0.120000, 0.030000,-0.550000, 0.130000,-0.210000,-0.380000,-1.220000,-0.220000,-0.370000, 0.070000,-0.070000, 0.040000,-0.010000,-0.040000},
	{-0.030000,-0.010000, 0.160000, 0.170000,-0.110000,-0.010000,-0.050000,-0.000000, 0.000000, 0.010000, 0.010000,-0.100000, 0.060000,-0.020000, 0.050000, 0.520000, 0.360000, 0.050000, 0.010000, 0.050000, 0.020000, 0.030000,-0.010000, 0.040000,-0.000000},
	{ 0.150000, 0.170000,-0.030000, 0.150000, 0.150000, 0.050000,-0.010000, 0.010000,-0.010000, 0.020000,-0.000000, 0.060000, 0.090000, 0.090000,-0.090000, 0.270000, 0.140000, 0.150000, 0.020000, 0.060000,-0.010000, 0.020000,-0.030000, 0.010000,-0.010000},
	{-0.480000, 0.020000, 0.020000, 0.000000,-0.140000,-0.030000,-0.070000,-0.000000, 0.010000, 0.010000, 0.000000, 0.000000, 0.010000,-0.080000,-0.030000,-0.000000, 0.040000,-0.290000,-0.030000,-0.110000, 0.030000,-0.050000, 0.020000,-0.020000, 0.000000},
	{-0.400000,-0.060000, 0.040000, 0.110000, 0.010000, 0.050000,-0.030000,-0.010000,-0.000000, 0.000000, 0.000000,-0.020000, 0.020000, 0.000000, 0.060000, 0.110000, 0.000000,-0.170000,-0.010000,-0.070000, 0.020000,-0.050000, 0.010000,-0.020000, 0.010000},
	{ 2.280000, 0.300000, 0.180000,-0.050000, 0.010000,-0.030000,-0.010000,-0.010000,-0.020000,-0.020000,-0.000000,-0.080000, 0.000000, 0.860000,-0.360000, 0.170000, 0.250000, 1.580000, 0.490000, 0.460000,-0.040000, 0.010000, 0.040000,-0.040000, 0.070000},
	{-0.160000, 0.440000, 0.340000,-0.160000, 0.040000,-0.010000, 0.020000, 0.000000, 0.000000, 0.000000, 0.000000,-0.020000,-0.040000,-0.180000, 0.080000, 0.230000, 0.170000,-0.060000,-0.030000,-0.000000,-0.010000, 0.000000, 0.000000, 0.000000, 0.000000},
	{-0.210000,-0.280000, 0.450000, 0.020000,-0.140000,-0.000000,-0.010000, 0.000000, 0.000000, 0.000000, 0.000000,-0.070000,-0.020000,-0.050000, 0.050000, 0.350000, 0.270000,-0.150000,-0.020000,-0.040000,-0.010000, 0.000000, 0.000000, 0.000000, 0.000000},
	{-0.100000,-0.310000, 0.190000, 0.110000,-0.050000,-0.080000, 0.030000, 0.000000, 0.000000, 0.000000, 0.000000, 0.010000,-0.010000,-0.070000, 0.060000,-0.050000,-0.030000, 0.000000, 0.010000, 0.010000, 0.020000, 0.000000, 0.000000, 0.000000, 0.000000},
	{-0.130000,-0.170000,-0.250000, 0.040000, 0.080000,-0.040000,-0.100000, 0.000000, 0.000000, 0.000000, 0.000000, 0.030000, 0.010000, 0.040000,-0.020000, 0.020000,-0.030000, 0.130000, 0.020000, 0.070000, 0.030000, 0.000000, 0.000000, 0.000000, 0.000000},
	{ 0.210000, 0.040000,-0.120000, 0.120000, 0.080000,-0.000000, 0.010000, 0.000000, 0.000000, 0.000000, 0.000000, 0.150000,-0.100000, 0.140000,-0.050000,-0.600000,-0.320000, 0.280000, 0.040000, 0.090000, 0.020000, 0.000000, 0.000000, 0.000000, 0.000000},
	{ 0.680000, 0.390000, 0.180000, 0.070000,-0.010000,-0.020000, 0.050000, 0.000000, 0.000000, 0.000000, 0.000000, 0.060000, 0.000000,-0.030000, 0.060000, 0.020000,-0.100000,-0.080000,-0.040000,-0.050000,-0.040000, 0.000000, 0.000000, 0.000000, 0.000000},
	{ 1.060000,-0.120000, 0.400000, 0.020000, 0.010000,-0.030000,-0.010000, 0.000000, 0.000000, 0.000000, 0.000000,-0.050000,-0.010000, 0.370000,-0.200000, 0.010000, 0.200000, 0.620000, 0.160000, 0.150000,-0.040000, 0.000000, 0.000000, 0.000000, 0.000000},
	{ 0.000000, 0.120000,-0.090000,-0.140000, 0.110000, 0.000000, 0.040000, 0.000000, 0.000000, 0.000000, 0.000000,-0.030000, 0.020000,-0.110000, 0.040000, 0.270000, 0.100000,-0.010000,-0.020000,-0.010000,-0.010000, 0.000000, 0.000000, 0.000000, 0.000000},
	{-0.120000,-0.000000, 0.210000,-0.140000,-0.120000,-0.030000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,-0.100000, 0.050000,-0.120000, 0.070000, 0.320000, 0.300000,-0.040000,-0.010000, 0.010000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}
};
/*** degree/order table for BDGIM Non-Broadcast Coefficient  [ 3/0 3/1 3/-1 3/2 ... 5/2 5/-2 ] ***/
const int NonBrdPara_degord_table[NONBRDNUM][2] = {
	{3,0},
	{3,1},
	{3,-1},
	{3,2},
	{3,-2},
	{3,3},
	{3,-3},
	{4,0},
	{4,1},
	{4,-1},
	{4,2},
	{4,-2},
	{5,0},
	{5,1},
	{5,-1},
	{5,2},
	{5,-2}
};
/*** degree/order table for BDGIM Broadcast Coefficient ***/
const int BrdPara_degord_table[BRDPARANUM][2] = {
	{0,0},
	{1,0},
	{1,1},
	{1,-1},
	{2,0},
	{2,1},
	{2,-1},
	{2,2},
	{2,-2}
};

/*****************************************************************************
* Description : Transform earth-fixed latitude/longitude into sun-fixed or/and geomagnetic latitude/longitude
* Parameters  :
*		double mjd		I		the calculate epoch
*		double *lat		I		input latitude (in arc.)
*		double *lon		I		input longitude (in arc)
*		int geomag		I		change into geomagnetic frame
*		int sunframe		I		change into sun-fixed frame
*		double *lat1		O		output latitude (in arc.)
*		double *lon1	O		output longitude (in arc.)
*****************************************************************************/
void EFLSFL(double mjd,double *lat,double *lon,int geomag,int sunframe,double *lat1,double *lon1)
{
	double Bp=LAT_POLE;
	double Lp=LON_POLE;
	double sinlon1 = 0.0, coslon1 = 0.0;

	if(geomag)
	{
		*lat1 = asin(sin(Bp*PI/180)*sin(*lat)+cos(Bp*PI/180)*cos(*lat)*cos(*lon-Lp*PI/180));
		sinlon1= (cos(*lat)*sin((*lon)-Lp*PI/180))/cos(*lat1);
		coslon1 = -(sin(*lat)-sin(Bp*PI/180)*sin(*lat1))/(cos(Bp*PI/180)*cos(*lat1));
		*lon1 = atan2(sinlon1,coslon1);
	}else
	{
		*lat1 = *lat;
		*lon1 = *lon;
	}

	if(sunframe)
	{
		double SUNLON = PI * (1.0 - 2.0 * (mjd - (int)(mjd)));
		double sinlon2 = 0.0, coslon2 = 0.0;
		sinlon2 = sin(SUNLON-Lp*PI/180);
		coslon2 = sin(Bp*PI/180)*cos(SUNLON-Lp*PI/180);
		SUNLON = atan2(sinlon2,coslon2);
		*lon1=*lon1-SUNLON;
		*lon1=atan2(sin(*lon1),cos(*lon1));
	}
}

/*****************************************************************************
* Description : different ionospheric mapping functions 
* Parameters  : int type	I
*					type = 0: none mapping function, return 1.0
*					type = 1: SLM mapping function
*					type = 2: MSLM mapping function by CODE (Schaer)
*					type = 3: Klobuchar mapping function
*				double ipp_elev  I    Ionospheric puncture point elevation
*				double sat_elev  I    satellite elevation
*				double Hion      I    Ionospheric layer height
* return: double IMF    mapping function factor
*****************************************************************************/
double IonMapping(int type, double ipp_elev, double sat_elev, double Hion)
{
	double IMF=0.0;
	double RE = 6378.0;	 // the average radius of earth [unit: km]

	if(Hion<10)
		Hion = 400000;

	switch(type) 
	{
	case 0:
		IMF=1.0;
		break;
	case 1:
		// SLM mapping function
		IMF=1/sin(ipp_elev);
		break;
	case 2:
		// MSLM mapping function by CODE (Schaer)
		IMF=1/sqrt(1-pow(RE*sin(0.9782*(PI/2-sat_elev))/(RE+Hion),2.0));
		break;
	case 3:
		// Klobuchar mapping function
		IMF=1+16*pow((0.53-sat_elev/PI),3.0);
		break;

	default:
		printf("can not find this kind of mapping model! The model type is %-4d\n",type);
		exit(-1);
	}
	return IMF;
}

/*****************************************************************************
* Description : calculate the XYZ, L, B and elevation of the IPP (accurate)
*					source code: BERNESE 5.0
* Parameters  :
*		double *sta		  I		station xyz	[m]			
*		double *sat		  I		satellite xyz [m]
*		double Hion	      I		ionosphere single layer height [m]
*		double *IPPXYZ	  O		xyz of IPP  [m]
*		double *IPP_B	  O		geographic latitude of the IPP [radian]
*		double *IPP_L	  O		geographic longitude  of the IPP [radian]
*		double *IPP_E	  O		elevation of the IPP [radian]
*		double *sat_ele	  O		elevation of satellite-station [radian]
* return: int
*****************************************************************************/
int IPPBLH1(double *sta,double *sat,double Hion,double *IPPXYZ,double *IPP_B,double *IPP_L,double *IPP_E, double* sat_ele)
{
	int n = 3;
	double dot = 0.0;
	while (--n >= 0) dot += sta[n] * sta[n];
	if (dot < 10e-6) return 0;
	n = 3; dot = 0.0;
	while (--n >= 0) dot += sat[n] * sat[n];
	if (dot < 10e-6) return 0;
	
	double DS,R1,R2,R3,zenith,zenith1,alpha;
	double sta_ipp;
	DS = Distance(sta,sat);
	R1 = sqrt(sta[0]*sta[0]+sta[1]*sta[1]+sta[2]*sta[2]);
	R2 = EARTH_RADIUS+Hion;
	R3 = sqrt(sat[0]*sat[0]+sat[1]*sat[1]+sat[2]*sat[2]);
	zenith = PI - acos((R1 * R1 + DS * DS - R3 * R3) / (2 * R1 * DS));
	*sat_ele = PI / 2 - zenith;
	// the zenit distance of the satellite from intersection with ionosphere
	zenith1 = asin(R1/R2*sin(zenith)); 
	//zenith2=asin(R1/R3*sin(zenith));
	alpha = zenith-zenith1;
	sta_ipp = sqrt(R1*R1+R2*R2-2*R1*R2*cos(alpha));
	// computation of coordinate of intersection of line receiver-satellite
	IPPXYZ[0]=sta[0]+sta_ipp*(sat[0]-sta[0])/DS;
	IPPXYZ[1]=sta[1]+sta_ipp*(sat[1]-sta[1])/DS;
	IPPXYZ[2]=sta[2]+sta_ipp*(sat[2]-sta[2])/DS;

	*IPP_L = atan2(IPPXYZ[1],IPPXYZ[0]);
	*IPP_B = atan(IPPXYZ[2]/sqrt(IPPXYZ[1]*IPPXYZ[1]+IPPXYZ[0]*IPPXYZ[0]));
	*IPP_E = PI/2.0-zenith1; 

	return 1;
}

/*****************************************************************************
* Description : calculate latitude, longitude and elevaiton of the IPP (approximate) according to user latitude ,longitude and satellite elevation, azimuth
* Parameters  :
*		double lat_u		I		user latitude  [unit: radian]
*		double lon_u		I		user longitude  [unit: radian]
*		double hion		    I		ionospheric single layer height [in meter]
*		double sat_ele	    I		satellite elevation [in radian]
*		double sat_azimuth  I	    satellite azimuth [in radian]
*		double &ipp_b       O		latitude of the ipp [in radian]
*		double &ipp_l	    O		longitude of the ipp [in radian]
*		double &ipp_e	    O		elevation of  the IPP [in radian]
* return: int
*****************************************************************************/
int IPPBLH2(double lat_u, double lon_u, double hion, double sat_ele, double sat_azimuth, double* ipp_b, double* ipp_l, double* ipp_e)
{
	double phiu = 0.0;
	double temp1 = 0.0, temp2 = 0.0;

	phiu = PI/2 - sat_ele - asin(EARTH_RADIUS*cos(sat_ele)/(EARTH_RADIUS+hion));
	*ipp_b = asin(sin(lat_u)*cos(phiu)+cos(lat_u)*sin(phiu)*cos(sat_azimuth));

	temp1 = sin(phiu)*sin(sat_azimuth)/cos(*ipp_b);
	temp2 = (cos(phiu)-sin(lat_u)*sin(*ipp_b))/(cos(lat_u)*cos(*ipp_b));

	*ipp_l = lon_u + atan2(temp1,temp2);
	*ipp_e = asin(EARTH_RADIUS*cos(sat_ele)/(EARTH_RADIUS+hion));

	return 1;
}

/****** calculate distance between xyz1 and xyz2 ******/
double Distance(double *xyz1,double *xyz2)
{
	return sqrt((xyz1[0]-xyz2[0])*(xyz1[0]-xyz2[0])+(xyz1[1]-xyz2[1])*(xyz1[1]-xyz2[1])+(xyz1[2]-xyz2[2])*(xyz1[2]-xyz2[2]));
}

/*****************************************************************************
* Description : Set the period term of the non-broadcast perdTable for BDGIM model
* Parameters  : 
*              NonBrdIonData* nonBrdData   I    BDGIM Non-Broadcast Ionospheric Parameters
*****************************************************************************/
void SetNonBrdCoefPeriod(NonBrdIonData* nonBrdData)
{
	double year_day = 365.25;       // days
	double solar_cycle = 4028.71;  // days
	double Moon_Month = 27.0;   // days
	double peridday = 0;
	int m = 0;
	int period_num = 0;

	nonBrdData->omiga[period_num++] = 0.0;
	// day period
	nonBrdData->omiga[period_num++] = 2*PI;
	nonBrdData->omiga[period_num++] = 2*PI/0.5;
	nonBrdData->omiga[period_num++] = 2*PI/0.33;
	// semi-month period
	nonBrdData->omiga[period_num++] = 2*PI/14.6;
	// month period
	nonBrdData->omiga[period_num++] = 2*PI/27.0;
	// one-third year period
	nonBrdData->omiga[period_num++] = 2*PI/121.6;
	// semi-year period
	nonBrdData->omiga[period_num++] = 2*PI/182.51;
	// year period
	nonBrdData->omiga[period_num++] = 2*PI/365.25;
	// solar period
	nonBrdData->omiga[period_num++] = 2*PI/4028.71;
	nonBrdData->omiga[period_num++] = 2*PI/2014.35;
	nonBrdData->omiga[period_num++] = 2*PI/1342.90;
	nonBrdData->omiga[period_num++] = 2*PI/1007.18;
}

/*****************************************************************************
* Description : Calculate the non-broadcast BDGIM parameters at the first compute epoch every day
* Parameters : 
*		double mjd		            I		the compute epoch [in mjd]
*       NonBrdIonData* nonBrdData   I       BDGIM Non-Broadcast Ionospheric Parameters
*****************************************************************************/
int CalNonBrdCoef(double mjd, NonBrdIonData* nonBrdData)
{
	if (mjd >= Init_mjd && (mjd - Init_mjd) < 1.0)
		return 1;

	double tmjd = 0.0, dmjd = 0.0, coef=0.0;
	int n=0, igroup = 0, icoef=0, ipar=0;

	dmjd = 2.0 / 24.0;

	for (tmjd = (int)(mjd); tmjd<(int)(mjd) + 1; tmjd = tmjd + dmjd)
	{
	// set the non-broadcast parameter
		for (icoef=0;icoef<NONBRDNUM;icoef++)
		{
			coef = 0.0; ipar = 0;
			for(n=0;n<PERIODNUM;n++)
			{
				if(nonBrdData->omiga[n]==0)
				{
					coef = nonBrdData->perdTable[icoef][ipar++];				
				}else
				{
					coef += nonBrdData->perdTable[icoef][ipar++]* cos(nonBrdData->omiga[n]*(tmjd + dmjd / 2.0));
					coef += nonBrdData->perdTable[icoef][ipar++]* sin(nonBrdData->omiga[n]*(tmjd + dmjd / 2.0));
				}
			}
			nonBrdData->nonBrdCoef[icoef][igroup] = coef;
		}
		igroup++;
	}

	Init_mjd = (int)(mjd);
	return 1;
}

/*****************************************************************************
* Description : find BDGIM non-broadcast coefficients group according to the mjd, and set the non-broadcast parameter
*		for BDGIM  model
* Parameters  :
*		double mjd		            I		the compute epoch (in mjd)
*       NonBrdIonData* nonBrdData   I       BDGIM Non-Broadcast Ionospheric Parameters
* return :    
*		int group		O		the session group of the current epoch
*
*****************************************************************************/
int BrdCoefGroupIndex(double mjd, NonBrdIonData* nonBrdData)
{
	double tmjd = 0.0, dmjd = 0.0;
	int igroup = -1;

	// calculate the non-broadcast parameters of BDGIM model
	if (mjd<Init_mjd || mjd>Init_mjd + 1.0)
		CalNonBrdCoef(mjd, nonBrdData);

	// set the sh coefficient group time interval
	dmjd = 2.0 / 24.0;

	for (tmjd = (int)(Init_mjd); tmjd<(int)(Init_mjd) + 1; tmjd = tmjd + dmjd)
	{
		if (mjd >= tmjd && mjd < tmjd + dmjd)
		{
			igroup++;
			break;
		}
		else
			igroup++;
	}

	if (igroup > 11)
		igroup = -1;

	return igroup;
}

/*****************************************************************************
 * Name        : VTECFromBroadSH
 * Description : Obtains the vertical ionospheric TEC using BDGIM ionospheric mode	
 * Parameters  :
 *      NonBrdIonData* nonBrdData   I               BDGIM Non-Broadcast Ionospheric Parameters
 *      BrdIonData* brdData			I               BDGIM Broadcast Ionospheric Parameters
 *		double mjd					I	[MJD]		The calculate time (Modified Julian Day)
 *		double lat					I	[arc]		The geomagnetic latitude of the Ionospheric Puncture Point (IPP)
 *		double lon                  I	[arc]		The geomagnetic longitude of the Ionospheric Puncture Point (IPP) 
 *		double *vtec			    O   [TECU]	    Ionospheric correction in TECU (in electrons per area unit)
 *****************************************************************************/
int VtecBrdSH(NonBrdIonData* nonBrdData, BrdIonData* brdData,double mjd,double ipp_b,double ipp_l, double* vtec)
{
	double temp =0.0;
	int ipar=0, igroup = 0;
	double ipp_b1=0.0,ipp_l1=0.0;
	double vtec_brd=0.0, vtec_A0=0.0;
	*vtec = 0.0;
	
	// obtains BDGIM model session group according to the mjd
	igroup = BrdCoefGroupIndex(mjd, nonBrdData);

	if (igroup == -1)
		return 0;

	ipp_b1 = ipp_b;
	ipp_l1 = ipp_l;

	// calculate the VTEC computed from the broadcast coefficients
	for(ipar=0;ipar< BRDPARANUM;ipar++)
	{
		if(brdData->degOrd[ipar][1]> brdData->degOrd[ipar][0])
			return 0;
		vtec_brd = vtec_brd + brdData->brdIonCoef[ipar]*(ASLEFU(ipp_b1,ipp_l1, brdData->degOrd[ipar][0], brdData->degOrd[ipar][1]));
	}

	// calculate the VTEC coomputed from the non-broadcast coefficients
	for (ipar=0;ipar<NONBRDNUM;ipar++)
	{
		if(nonBrdData->degOrd[ipar][1]>nonBrdData->degOrd[ipar][0])
			return 0;
		vtec_A0 = vtec_A0 + nonBrdData->nonBrdCoef[ipar][igroup]*(ASLEFU(ipp_b1,ipp_l1,nonBrdData->degOrd[ipar][0],nonBrdData->degOrd[ipar][1]));
	}

	*vtec = vtec_brd + vtec_A0;

	if (brdData->brdIonCoef[0]  > 35.0)
		*vtec = MAX(brdData->brdIonCoef[0] / 10.0, *vtec);
	else if (brdData->brdIonCoef[0] > 20.0)
		*vtec = MAX(brdData->brdIonCoef[0] / 8.0, *vtec);
	else if (brdData->brdIonCoef[0] > 12.0)
		*vtec = MAX(brdData->brdIonCoef[0] / 6.0, *vtec);
	else
		*vtec = MAX(brdData->brdIonCoef[0] / 4.0, *vtec);

	return 1;
}

/*****************************************************************************
* Description : Normalized legendre polynomial
* Parameters  :
*			INN: Degree     IMM: Order:
*****************************************************************************/	
double ASLEFU(double XLAT,double XLON,int INN,int IMM)
{
	int NN,MM;
	int II,LL;
	double XX,PMM,FACT,SOMX2,COEF,PMMP1,result,FACTN,PLL;
	double KDELTA;    
	NN=abs(INN);
	MM=abs(IMM);

	// compute the associated legendre polynomial
	if(INN>=0) 
		XX=sin(XLAT);
	else 
		XX=cos(XLAT);
	PMM=1.0;
	if(MM>0) 
	{
		SOMX2=sqrt((1.0-XX)*(1.0+XX));
		FACT=1.0;
		for(II=1;II<=MM;II++)
		{
			PMM=PMM*FACT*SOMX2;
			FACT=FACT+2.0;
		}
	}

	if(NN==MM) 
		COEF=PMM;
	else
	{
		PMMP1=XX*(2*MM+1)*PMM;
		if(NN==MM+1) 
			COEF=PMMP1;
		else
		{
			for(LL=MM+2;LL<=NN;LL++)
			{
				PLL=(XX*(2*LL-1)*PMMP1-(LL+MM-1)*PMM)/(LL-MM);
				PMM=PMMP1;
				PMMP1=PLL;
			}
			COEF=PLL;
		}
	}

	/* compute the normalization factor and (normalized) coefficient */
	if(MM==0) 
		KDELTA=1.0;
	else 
		KDELTA=0;

	FACTN=sqrt(2.0*(2.0*NN+1.0)/(1.0+KDELTA)* FAKULT(NN-MM)/FAKULT(NN+MM));

	if(IMM>=0) 
		result=FACTN*COEF*cos(MM*XLON);
	else
		result=FACTN*COEF*sin(MM*XLON);
	return result;
}

/* compute the factorial of N */
double FAKULT(int N)
{
	double FAK;
	int I;
	FAK=1.0;
	if(N<=1) return FAK;
	for(I=2;I<=N;I++)  
		FAK=FAK*I;
	return FAK;
}

/*****************************************************************************
* Description : obtains the slant ionospheric delay in B1C using BDGIM ionospheric model
* Parameters  :
*      NonBrdIonData* nonBrdData   I    BDGIM Non-Broadcast Ionospheric Parameters
*      BrdIonData* brdData         I    BDGIM Broadcast Ionospheric Parameters
*	   double mjd			       I	current epoch
*	   double* sta_xyz		       I	station x,y,z
*	   double* sat_xyz		       I	satellite x,y,z
*	   double* brdPara		       I	broadcast ionospheric parameters [const: 9 parameters model]
*	   double* ion_delay	       O	ionospheric delay in B1C [m]
*****************************************************************************/
int IonBdsBrdModel(NonBrdIonData* nonBrdData, BrdIonData* brdData, double mjd, double* sta_xyz, double* sat_xyz, double* brdPara, double* ion_delay) {
	double ipp_xyz[3] = { 0.0 };
	double ipp_b = 0.0, ipp_l = 0.0, ipp_e = 0.0;
	double geomag_b = 0.0, geomag_l = 0.0;
	double sat_ele = 0.0;
	double mf = 0.0, K = 0.0;
	double vtec = 0.0;
	int num, i, j;

	// 1:set basic parameters : period, degree/order, omiga, non-broadcast parameters
	for (num = 0; num < NONBRDNUM; num++) {
		nonBrdData->degOrd[num][0] = NonBrdPara_degord_table[num][0];
		nonBrdData->degOrd[num][1] = NonBrdPara_degord_table[num][1];

		for (j = 0; j < TRISERINUM; j++) {
			nonBrdData->perdTable[num][j] = NonBrdPara_table[num][j];
		}
	}

	for (i = 0; i < BRDPARANUM; i++) {
		brdData->degOrd[i][0] = BrdPara_degord_table[i][0];
		brdData->degOrd[i][1] = BrdPara_degord_table[i][1];

		brdData->brdIonCoef[i] = brdPara[i];
	}
	SetNonBrdCoefPeriod(nonBrdData);

	// 2:calculate IPP information
	IPPBLH1(sta_xyz, sat_xyz, Hion_bdgim, ipp_xyz, &ipp_b, &ipp_l, &ipp_e,&sat_ele);

	// 3:Transform earth-fixed coordinate to sun-fixed and geomagnetic coordinate
	EFLSFL(mjd, &ipp_b, &ipp_l, 1, 1, &geomag_b, &geomag_l);

	// 4:Calcute the vertical ionospheric TEC
	VtecBrdSH(nonBrdData, brdData,mjd, geomag_b, geomag_l, &vtec);

	// 5:Calcute the mapping factor
	mf=IonMapping(2, ipp_e, sat_ele, Hion_bdgim);

	// 6:calculate the Delay conversion factor
	K = 40.3e16 / (1.0*pow(FREQ1_BDS, 2)); 

	// 7:calculate the ionospheric Delay in BDS B1C frequency
	*ion_delay = mf * K * vtec;

	return 1;
}
