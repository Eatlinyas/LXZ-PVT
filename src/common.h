#pragma once
/**************************************************************
* common.h : declare some functions related to coordinate and time conversion.
***************************************************************/

typedef struct
{
	double mjd;
	double daysec;
} MjdData;

/* calculate MJD from UTC --------------------------------*/
void UTC2MJD(int year, int month, int day, int hour, int min, double second, MjdData* mjdata);