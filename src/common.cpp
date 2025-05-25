#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include"common.h"

/*****************************************************************************
* Name:  UTC2MJD
* Description :  calculate MJD from UTC 
*	Return mjd
*****************************************************************************/
 void UTC2MJD(int year, int month, int day, int hour, int min, double second, MjdData* mjdata)
{
    double hourn = 0.0;
    double m_julindate = 0.0;

    memset(mjdata, 0, sizeof(MjdData));

    if (year < 80)
        year = year + 2000;
    else if (year >= 80 && year <= 1000)
    {
        year = year + 1900;
    }

    hourn = hour + min / 60.0 + second / 3600.0;
    if (month <= 2)
    {
        year -= 1;
        month += 12;
    }

    m_julindate = (int)(365.25 * year) + (int)(30.6001 * (month + 1)) + day + hourn / 24.0 + 1720981.5;

    mjdata->mjd = m_julindate - 2400000.5;

    mjdata->daysec = hour * 3600.0 + min * 60.0 + second;

    if (mjdata->daysec == 86400.0)
    {
        mjdata->daysec = 0.0;
    }

    return;
}