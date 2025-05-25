/************************************************************************************
 *   main.cpp : A simple example of calling BDGIM module.
 *   We support an interface to calculate the ionospheric delay in BDS B1C frequency.
 ************************************************************************************/

#include "bdgim.h"
#include"common.h"

int main()
{
	MjdData mjdData;
	NonBrdIonData nonBrdData;				 // non-broadcast parameter structure
	BrdIonData brdData;                      // broadcast parameter structure
		
	int year=0, month=0, day=0, hour=0, min=0;
	double second = 0.0;
	double mjd = 0.0;
	double sat_xyz[3] = { -34342443.2642, 24456393.0586, 12949.4754 };
	double sta_xyz[3] = { -2159945.3149, 4383237.1901, 4085412.3177 };
	double ion_delay = 0.0;

	memset(&mjdData, 0, sizeof(MjdData));
	memset(&nonBrdData, 0, sizeof(NonBrdIonData));
	memset(&brdData, 0, sizeof(BrdIonData));
	
	year = 2020; month = 11; day = 5; hour = 17; min = 0; second = 0.0;
	UTC2MJD(year, month, day, hour, min, second, &mjdData);
	mjd = mjdData.mjd;

	double brdPara[9] = {
		13.82642069, -1.78004616, 5.17200000,   //BDS1
		3.13030940, -4.58961329, 0.32483867,    //BDS2
		-0.07802383, 1.24312590, 0.37763627     //BDS3
	};

	/* the BDGIM module interface */
	IonBdsBrdModel(&nonBrdData, &brdData, mjd, sta_xyz, sat_xyz, brdPara, &ion_delay);
	printf(" the ionosphere delay in B1C : %7.2lf [m]\n", ion_delay);
	
	return 0;
}
