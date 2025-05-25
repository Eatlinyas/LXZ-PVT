
//#include "bdssh.h"
#include "rtklib.h"

#include "bdgim.h"
#include "common.h"


BDSSH::BDSSH(void)
{
	int i;
	this->BrdIonCoefNum = 0;															// the parameter number for broadcast [9 in default] 
	this->bds_ion.nk8 = 0;
	this->bds_ion.nk14 = 0;
	this->bds_ion.nsh9 = 0;// the total parameter number for SH resolution [26 in default] 
	memset(this->BrdIonCoef, 0, BRDCOUNT*MAXGROUP*sizeof(double));			        // the body array for storing the BDSSH broadcast parameter														// the begin epoch of the compute epoch [in MJD]
}

BDSSH::~BDSSH(void)
{

}


void sort(double *data, int n){
	int i, j;
	double help, tmp;
	for (i = 0; i < n; i++){
		tmp = fabs(data[i]);
		for (j = i + 1; j< n; j++){
			if (tmp > fabs(data[j])){
				help = data[j];
				data[j] = data[i];
				data[i] = help;
			}
		}
	}
}

/* unique ion ----------------------------------------------------------
* unique ion in navigation data and update carrier wave length
* args   : nav_t *nav    IO     navigation data
* return : number of epochs
*-----------------------------------------------------------------------------*/
int  selemajority(double *tmp, int n){
	int i, j, k, l;
	int freq[512], nf = 0;
	if (n <= 0)return -1;

	sort(tmp, n);
	for (i = 0, j = 0; i < n; i++){
		if (tmp[j] == tmp[i])nf++;
		else{
			freq[j] = nf;
			tmp[++j] = tmp[i];
			nf = 1;
		}
	}
	freq[j] = nf;
	k = freq[0];
	l = 0;
	for (i = 0; i <= j; i++){
		if (k<freq[i]) {
			k = freq[i];
			l = i;
		}
	}
	for (i = 0; i < n; i++){
		if (tmp[i] == tmp[l])return i;
	}
	return -1;
}

bool BDSSH::uniqion(double *ep)
{
	int i, j, k;
	int sys, nigmas = 0;
	double tmp[512];
	int n = 0;
	int ii, idx;
	double fjd;
	int year, month, day;

	trace(3, "uniqion:\n");

	if (bds_ion.nsh9 <= 0)return -1;

	this->BrdIonCoefNum = 9;
	for (i = 0; i < bds_ion.nsh9; i++){
		if (bds_ion.bdssh9[i].hour != -1)nigmas++;
	}
	if (nigmas == 0)
	{ /* not igmas format, keeps only the first one */
		this->BrdIonCoefGroup = 1;
		for (j = 0; j < this->BrdIonCoefNum; j++){
			this->BrdIonCoef[j][0] = bds_ion.bdssh9[0].ion[j];
		}
	}
	else{
		this->BrdIonCoefGroup = 24;
		for (i = 0; i < 24; i++){
			n = 0;
			for (j = 0; j < bds_ion.nsh9; j++){
				if (bds_ion.bdssh9[j].hour != i)continue;
				tmp[n] = norm(bds_ion.bdssh9[j].ion, this->BrdIonCoefNum);
				if (tmp[n] == 0.0)continue;
				n++;
			}
			ii = selemajority(tmp, n);
			if (ii < 0){
				for (j = 0; j < this->BrdIonCoefNum; j++){
					this->BrdIonCoef[j][i] = 0.0;
				}
			}
			else
			{
				idx = -1;
				for (j = 0; j < bds_ion.nsh9; j++){
					if (bds_ion.bdssh9[j].hour != i)continue;
					if (tmp[ii] == norm(bds_ion.bdssh9[j].ion, 9)) {
						for (k = 0; k < this->BrdIonCoefNum; k++){
							this->BrdIonCoef[k][i] = bds_ion.bdssh9[j].ion[k];
						}
						break;
					}
				}
				//if (idx == -1){ printf("Error!\n"); }
			}
		}
	}
	return true;
}



/* ionosphere model ------------------------------------------------------------
* compute ionospheric delay by broadcast ionosphere model (klobuchar model - BDS K8)
* args   : gtime_t t        I   time (gpst)
*          double *ion      I   iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3}
*                               broadcasted by BDS K8
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric delay (B1I) (m)
*-----------------------------------------------------------------------------*/

extern double ionmodel_BDSK8(gtime_t t, const double *ion, const double *pos,
	const double *azel)
{
	const double ion_default[] = { /* 2004/1/1 */
		2.7940E-009,  3.7253E-008, -2.9802E-007,  4.1723E-007,
		1.4950E+005, -7.2090E+005,  7.0779E+006, -7.8643E+006
	};

	double tt, f, psi, phi, lam, amp, per, x, k, R = RE_WGS84, H = 375000;
	int week;
	double vtime;
	double sow;
	gtime_t bdt;
	//double tmp[4];



	if (pos[2]<-1E3 || azel[1] <= 0) return 0.0;
	if (norm(ion, 8) <= 0.0){
		//printf("non K8!\n");
		ion = ion_default;
	}

	psi = PI / 2 - azel[1] - asin(R / (R + H)*cos(azel[1]));

	phi = asin(sin(pos[0])*cos(psi) + cos(pos[0])*sin(psi)*cos(azel[0]));
	lam = pos[1] + asin(sin(psi)*sin(azel[0]) / cos(phi));

	/* geomagnetic latitude (semi-circle) */

	/* local time (s) */
	sow = time2bdt(t, &week);
	tt = 43200.0*lam / PI + sow;
	tt -= floor(tt / 86400.0)*86400.0; /* 0<=tt<86400 */

	/* slant factor */
	f = R / (R + H)*cos(azel[1]);
	f = 1.0 / sqrt(1 - f*f);

	/* ionospheric delay */
	phi = fabs(phi / PI);
	amp = ion[0] + phi*(ion[1] + phi*(ion[2] + phi*ion[3]));
	per = ion[4] + phi*(ion[5] + phi*(ion[6] + phi*ion[7]));
	amp = amp<    0.0 ? 0.0 : amp;
	per = per<72000.0 ? 72000.0 : per;
	per = per>172800.0 ? 172800.0 : per;
	x = 2.0*PI*(tt - 50400.0) / per;
	vtime = fabs(x) < 1.5708 ? 5E-9 + amp*cos(x) : 5E-9;
	double varr = azel[1] * R2D<15 ? (0.6 + sin(azel[1])) : 1.0;
	return CLIGHT*f*vtime*varr;
}

/* BDSSH9 for B1I */
extern double ionmodel_BDSK9(gtime_t time, const BDSSH *bdssh, const double *pos, const double *satxyz)
{

	MjdData mjdData;
	NonBrdIonData nonBrdData;				 // non-broadcast parameter structure
	BrdIonData brdData;                      // broadcast parameter structure
	double ep[6], brdPara[9], sta_xyz[3], sat_xyz[3],iondelay=0.0;
	memset(&mjdData, 0, sizeof(MjdData));
	memset(&nonBrdData, 0, sizeof(NonBrdIonData));
	memset(&brdData, 0, sizeof(BrdIonData));
	BDSSH *bdsk9 = (BDSSH*)bdssh;
	int igroup = -1,i,j;
	double ionsh9[9],dt;

	time2epoch(time, ep);
	

	/* select reference ion par and ref hour */
	if (bdsk9->BrdIonCoefGroup == 1){
		igroup = 0;
		for (i = 0; i < bdsk9->BrdIonCoefNum; i++)brdPara[i] = bdsk9->BrdIonCoef[i][0];
	}
	else if (bdsk9->BrdIonCoefGroup == 24){
		dt = 48.0;
		for (i = 0; i<24; i++){
			for (j = 0; j < bdsk9->BrdIonCoefNum; j++)ionsh9[j] = bdsk9->BrdIonCoef[j][i];
			if (norm(ionsh9, bdsk9->BrdIonCoefNum) == 0.0)continue;
			if (fabs(ep[3] - i) < dt){
				dt = fabs(ep[3] - i);
				igroup = i;
				for (j = 0; j < bdsk9->BrdIonCoefNum; j++)brdPara[j] = ionsh9[j];
			}
		}
	}
	if (dt >= 48.0 || igroup<0)return 0.0;


	UTC2MJD((int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5], &mjdData);
	pos2ecef(pos, sta_xyz);
	for (i = 0; i < 3; i++)sat_xyz[i] = satxyz[i];

	IonBdsBrdModel(&nonBrdData, &brdData, mjdData.mjd, sta_xyz, sat_xyz, brdPara, &iondelay);
	//printf(" the ionosphere delay in B1C : %7.2lf [m]\n", iondelay);

	// B1C to B1I 
	double k = FREQ1*FREQ1/ FREQ1_CMP / FREQ1_CMP;
	iondelay = iondelay*k;

	return iondelay<0.0 ? 0.0 : iondelay;
}

/* BDSSH9 for B1I */
extern double ionmodel_BDSK9_2(gtime_t time,  double* brdPara, const double *pos, const double *satxyz)
{

	MjdData mjdData;
	NonBrdIonData nonBrdData;				 // non-broadcast parameter structure
	BrdIonData brdData;                      // broadcast parameter structure
	double ep[6], sta_xyz[3], sat_xyz[3], iondelay = 0.0;
	memset(&mjdData, 0, sizeof(MjdData));
	memset(&nonBrdData, 0, sizeof(NonBrdIonData));
	memset(&brdData, 0, sizeof(BrdIonData));
	int igroup = -1, i, j;
	double ionsh9[9], dt;

	time2epoch(time, ep);


	
	UTC2MJD((int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5], &mjdData);
	pos2ecef(pos, sta_xyz);
	for (int i = 0; i < 3; i++)sat_xyz[i] = satxyz[i];
	IonBdsBrdModel(&nonBrdData, &brdData, mjdData.mjd, sta_xyz, sat_xyz, brdPara, &iondelay);
	printf(" the ionosphere delay in B1C : %7.2lf [m]\n", iondelay);

	// B1C to B1I 
	double k = FREQ1*FREQ1 / FREQ1_CMP / FREQ1_CMP;
	iondelay = iondelay*k;

	return iondelay<0.0 ? 0.0 : iondelay;
}

//#endif