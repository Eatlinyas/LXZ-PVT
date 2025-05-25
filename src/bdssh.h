/***************************************************************************
*   Copyright (C) 2012 by Institute of Geodesy and Geophysics, CAS, Wuhan   *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                             *
*		Create at  2012,	Apr. 8
*			update 1: 2014.04.28	defines the non-broadcast bdssh model coefficients structure
***************************************************************************/
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#ifndef BDSSH_H
#define BDSSH_H

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


#define BRDCOUNT			9

#ifndef MAXPERIODCOUNT
#define MAXPERIODCOUNT  20
#endif



typedef struct {        /* GPS/QZS/GAL broadcast ephemeris type */
	int sat;            /* satellite number */
	int iodi;     /* IODI,reserved */
	int hour;
	double ion[8];      /* Klobuchar 8 parameters */
} bdsk8_t;

typedef struct {        /* GPS/QZS/GAL broadcast ephemeris type */
	int sat;            /* satellite number */
	int iodi;     /* IODI,reserved */
	int hour;
	double ion[14];      /* Klobuchar 8 parameters */
} bdsk14_t;

typedef struct {        /* GPS/QZS/GAL broadcast ephemeris type */
	int sat;            /* satellite number */
	int iodi;     /* IODI,reserved */
	int hour;
	double ion[9];      /* bds sh9 parameters */
} bdssh9_t;

typedef struct {        /* navigation data type */
	int nk8, nk14, nsh9;
	bdsk8_t bdsk8[1024];
	bdssh9_t bdssh9[1024];
}bds_ion_t;

class BDSSH
{
public:
	int BrdIonCoefNum;
	BDSSH(void);
	~BDSSH(void);

	bool uniqion(double *epoch);

	int GetBrdCoefNum(){ return BrdIonCoefNum; }
	double BrdIonCoef[BRDCOUNT][24];			// the body array for storing the bdssh broadcast parameter
	bds_ion_t bds_ion;
	int BrdIonCoefGroup;
private:


};

#endif
