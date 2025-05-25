/****************************************************************************************************************
 *   bdgim.h : BDGIM model constants, basic structures and function prototypes.
 *   BDS Broadcast Ionospheric model(BDGIM) calculation module
 *   Author: Zishen LI, Ningbo WANG
 *   e-mail: lizishen@aircas.ac.cn, wangnb@aircas.ac.cn
 *
 *   Copyright (C) 2021 by Aerospace Information Research Institute, Chinese Academy of Sciences, Beijing, China.
 *   This program is free software; you can redistribute it and/or modify it under the terms of the GNU General 
 *   Public License as published by the Free Software Foundation; either version 2 of the License, or (at your 
 *   option) any later version. 
 *
 *	 Create at  2021,	Apr. 23
 ***************************************************************************************************************/

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/*********** define some basic consts **************************/

#define PI	          (4.0*atan(1.0)) 
#define MIN(a,b)   (((a)>(b))?(b):(a))    // minimum between a and b
#define MAX(a,b)  (((a)<(b))?(b):(a))     // maximum between a and b
#define FREQ1_BDS       1575420000.0      // BDS-3 B1C  frequency (Hz) 
#define BRDPARANUM      9                 // Number of broadcast ionospheric parameters
#define PERIODNUM      13                 // Number of forecast period for every non-broadcast parameters
#define TRISERINUM    (PERIODNUM*2-1)     // Trigonometric series number 
#define NONBRDNUM      17                 // Number of non-broadcast one group
#define MAXGROUP       12                 // 12 groups non-broadcast coefficient every day

/*********** BDGIM Non-Broadcast Ionospheric Parameters Struct **************************/
typedef struct{
	int degOrd[NONBRDNUM][2];		              // spheric harmonic function degree and order index	
	double omiga[PERIODNUM];				      // omiga calculated from the period
	double perdTable[NONBRDNUM][PERIODNUM*2-1];	  // the array for storing the non-broadcast parameter period table for perdTable
	double nonBrdCoef[NONBRDNUM][MAXGROUP];	      // the non-broadcast bdgim parameter of the calculate day
} NonBrdIonData;

/*********** BDGIM Broadcast Ionospheric Parameters Struct **************************/
typedef struct{
	double brdIonCoef[BRDPARANUM];		// the body array for storing the bdgim broadcast parameter
	int degOrd[BRDPARANUM][2];		    // spheric harmonic function degree and order index for broadcast parameter
} BrdIonData;

/*-------------------------------- BDGIM model function ----------------------------*/

/* obtains the slant ionospheric delay in B1C using BDGIM ionospheric model */
 int IonBdsBrdModel(NonBrdIonData* nonBrdData, BrdIonData* brdData, double mjd, double* sta_xyz, double* sat_xyz, double* brdPara, double* ion_delay);


/* obtains the vertical ionospheric tec using BDGIM ionospheric mode */
int VtecBrdSH(NonBrdIonData* nonBrdData, BrdIonData* brdData, double mjd, double ipp_b, double ipp_l, double* vtec);

/* transform earth-fixed latitude/longitude into sun-fixed latitude/longitude */
void EFLSFL(double mjd, double* lat, double* lon, int geomag, int sunframe, double* lat1, double* lon1);

/* ionospheric mapping functions */
double IonMapping(int type, double ipp_elev, double sat_elev, double Hion);

/* calculate XYZ, latitude, longitude and elevation of the IPP (accurate) according to user and satellite position */
int IPPBLH1(double* sta, double* sat, double Hion, double* IPPXYZ, double* IPP_B, double* IPP_L, double* IPP_E, double* sat_ele);

/* calculate latitude, longitude and elevation of the IPP (approximate) according to user latitude ,longitude and satellite elevation, azimuth */
int IPPBLH2(double lat_u, double lon_u, double hion, double sat_ele, double sat_azimuth, double* ipp_b, double* ipp_l, double* ipp_e);

int BrdCoefGroupIndex(double mjd, NonBrdIonData* nonBrdData);	// obtains non-broadcast coefficient session group according to the mjd
int CalNonBrdCoef(double mjd, NonBrdIonData* nonBrdData);		// calculate the non-broadcast BDGIM coefficients at the first compute epoch
void SetNonBrdCoefPeriod(NonBrdIonData* nonBrdData);		    // Set the period term of the non-broadcast perdTable for BDGIM model
double ASLEFU(double XLAT,double XLON,int INN,int IMM);		    // Normalized legendre polynomial
double FAKULT(int N);										    // compute the factorial of N
double Distance(double *xyz1,double *xyz2);

