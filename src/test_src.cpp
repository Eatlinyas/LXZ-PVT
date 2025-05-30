#include "rtklib.h"
#include "math.h"
#include "stdio.h"
#include "DBSCAN.h"
#define SQR(x)      ((x)*(x))
// threshold of cycle-slip
#define THRES_SLIP  2.0             

static int edited = 0;
static int IsWriteHeader = 0;
static int IsWriteSNRHeader = 0;
//sys,selectedfrq
static int selectedfrqs[6][7] = { 0 };
double *Az, *El, *Nadir, *AzInSa, *Mp[NFREQ + NEXOBS];
char outpath[1024];
FILE* fpSat;


extern double permutation(int flag, const double *x0, const int nx0, const double ratio)
{
	int i, j, index, nx;
	double help;
	double *x = zeros(nx0, 1);

	for (i = 0, nx = 0; i < nx0; i++){
		if (fabs(x0[i])>1E-6){ x[nx++] = x0[i]; }
	}
	for (i = 0; i < nx; i++){
		for (j = i + 1; j < nx; j++){
			if (flag == 1){
				if (x[i]>x[j]){
					help = x[j];
					x[j] = x[i];
					x[i] = help;
				}
			}
			else {
				if (x[i] < x[j]){
					help = x[i];
					x[i] = x[j];
					x[j] = help;
				}
			}
		}
	}
	index = nx>0 ? (int)round(ratio*nx) - 1 : -1;
	if (index >= 0){
		return x[index];
	}
	else
	{
		return 0.0;
	}
}


extern double median_t(double *data, int n){
	double median, MAD, mean_t = 0.0;
	double vv[MAXSAT];
	int i, nok = 0;

	median = permutation(1, data, n, 0.5);
	for (i = 0; i < n; i++){
		//printf("%d %lf\n", i + 1, data[i]);
		vv[i] = fabs(data[i] - median) / 0.6745;
	}
	MAD = permutation(1, vv, n, 0.5);
	for (i = 0; i<n; i++) {
		vv[i] = vv[i] * 0.6745; /* recover to ddres */
		//if (fabs(v[i] - median)>6 * MAD){
		if (fabs(vv[i])<4 * MAD){
			nok++;
			mean_t = mean_t + (data[i] - mean_t) / nok;
		}
		//else { printf("%d %lf\n", i + 1, data[i]); }
	}
	for (i = 0; i<n; i++) {
		data[i] = data[i] - mean_t;
		data[i] = fabs(data[i])>99999.9999 ? 99999.9999 : data[i];
	}
	return mean_t;
}

extern double median2_t(double *res, int npt){
	vector<point> data;
	double mean_t = 0.0;
	int i, j,n,nok=0;
	for (i = 0; i < npt; i++){
		if (fabs(res[i]) < 1E-9)continue;
		vector<float> xn;
		xn.push_back(res[i]);	
		//printf("%d: %lf\n", i,res[i]);
		point p(xn, nok++);
		data.push_back(p);
	}

	mean_t=DBSCAN(data, 1.5, (int)(nok*0.5));

	for (i = 0; i < npt; i++)
	{
		if (fabs(res[i]) < 1E-9)continue;
		res[i] -= mean_t;
	}
	trace(3, "mean(%d): %lf\n", nok,mean_t);
	return mean_t;
}


extern double median3_t(double *data, int n,double *mad){
	double median, MAD, mean_t = 0.0;
	double vv[MAXSAT];
	int i, nok = 0;

	median = permutation(1, data, n, 0.5);	//每个可用卫星都有伪距值，返回的median为伪距中位数
	for (i = 0; i < n; i++){
		trace(1,"%d %lf\n", i + 1, data[i]);
		vv[i] = fabs(data[i] - median) / 0.6745;
	}
	MAD = permutation(1, vv, n, 0.5);
	for (i = 0; i<n; i++) {
		vv[i] = vv[i] * 0.6745; /* recover to ddres */
		//if (fabs(v[i] - median)>6 * MAD){
		if (fabs(vv[i])<4 * MAD){
			nok++;
			mean_t = mean_t + (data[i] - mean_t) / nok;
		}
		else { trace(1,"%d* %lf\n", i + 1, data[i]); }
	}
	for (i = 0; i<n; i++) {
		data[i] = data[i] - mean_t;
	}
	*mad = MAD;
	return mean_t;
}


int getSatmeaaageHeaderGNSSSYS(int sys, char line[])
{
	int i, nsat;
	char SYS, result[4096] = "", data[128];
	switch (sys)
	{
	case SYS_GPS:nsat = NSATGPS; SYS = 'G'; break;
	case SYS_GLO:nsat = NSATGLO; SYS = 'R'; break;
	case SYS_GAL:nsat = NSATGAL; SYS = 'E'; break;
	case SYS_QZS:nsat = NSATQZS; SYS = 'J'; break;
	case SYS_CMP:nsat = NSATCMP; SYS = 'C'; break;
	case SYS_SBS:nsat = NSATSBS; SYS = 'S'; break;
	default:nsat = 3; SYS = 'N'; break;
	}
	for (i = 0; i<nsat; i++)
	{
		sprintf(data, "%10c%2.2d", SYS, i + 1);
		/*sprintf(data,"%58c%2d",SYS,i+1);*/
		strcat(result, data);
	}
	strcpy(line, result);
	return 1;
	//printf("%s",hehe);
	//printf("hehe\n");
}

extern void outsatres(rtk_t* rtk, int* sat, int ns)
{
	int i, j, k;
	char line[4096];
	//rtk->opt.
	//if(!IsOpen)fpSat=fopen("E:\\learnprogram\\ReBuild_RTKLIB\\option\\SatStatis.txt","w");
	if (!fpSat)fpSat = fopen(strcat(outpath, "line_res"), "w");
	//if(!fpSat_p)fpSat_p=fopen(strcat(outpath,"psu-line_res"),"w");
	if (!IsWriteHeader)
	{
		for (i = 0; i<6; i++)
		{
			line[0] = '\0';
			if (rtk->opt.navsys&(int)pow(2.0, i))
				getSatmeaaageHeaderGNSSSYS((int)pow(2.0, i), line);
			fputs(line, fpSat);
			//fputs(line,fpSat_p);	
		}
		fputs("\n", fpSat);//fputs("\n",fpSat_p);
		//rtk->opt.navsys
		IsWriteHeader++;
	}

	for (i = 1; i <= MAXSAT; i++)
	{
		if (!(rtk->opt.navsys&satsys(i, NULL)))continue;
		for (j = 0; j<ns; j++)
		{
			if (sat[j] == i)break;
		}
		if (j >= ns)
		{
			sprintf(line, "%12d", 0);
			fputs(line, fpSat);//fputs(line,fpSat_p);
		}
		else
		{
			sprintf(line, "%12.4f", rtk->ssat[sat[j] - 1].resc[0]);
			fputs(line, fpSat);
			//sprintf(line,"%12.4f",rtk->ssat[sat[j]-1].resp[0]);
			//fputs(line,fpSat_p);
		}
		/*if(j>=ns)
		for(k=0;k<5;k++)
		{
		sprintf(line,"%12d",0);
		fputs(line,fpSat);fputs(line,fpSat_p);
		}
		else for(k=0;k<5;k++)
		{
		sprintf(line," %11.4f",rtk->ssat[sat[j]-1].resc[k]);
		fputs(line,fpSat);
		sprintf(line," %11.4f",rtk->ssat[sat[j]-1].resp[k]);
		fputs(line,fpSat_p);
		}*/
	}

	fputs("\n", fpSat);//fputs("\n",fpSat_p);
}



extern void outsatres_single(FILE *fpSat_p, rtk_t* rtk, obsd_t *obs, int n)
{
	int i, j;
	char line[2048];
	gtime_t time;
	/*FCB *myfcb;*/
	//rtk->opt.
	//if(!IsOpen)fpSat=fopen("E:\\learnprogram\\ReBuild_RTKLIB\\option\\SatStatis.txt","w");

	if (!IsWriteHeader)
	{
		sprintf(line, "%23s", "%BDT            Obs_Time"); fputs(line, fpSat_p);
		for (i = 0; i<6; i++)
		{
			line[0] = '\0';
			if (rtk->opt.navsys&(int)pow(2.0, i))
				getSatmeaaageHeaderGNSSSYS((int)pow(2.0, i), line);
			fputs(line, fpSat_p);
		}
		fputs("\n", fpSat_p);
		//rtk->opt.navsys
		IsWriteHeader++;
	}
	time = rtk->sol.time;
	if (rtk->tsys == TSYS_CMP){
		time = gpst2bdt(rtk->sol.time);
	}
	time2str(time, line, 3); fputs(line, fpSat_p);
	for (i = 1; i <= MAXSAT; i++)
	{
		if (!(rtk->opt.navsys&satsys(i, NULL)))continue;
		for (j = 0; j<n; j++)
		{
			if (obs[j].sat == i)break;
		}
		if (j >= n ) sprintf(line, "%12d", 999999);
		else if (fabs(rtk->ssat[i - 1].resp[0]) >= 9999.9) sprintf(line, "%12.4f", 9999.9);
		else sprintf(line, "%12.4f", rtk->ssat[i - 1].resp[0]);
		fputs(line, fpSat_p);
	}
	fputs("\n", fpSat_p);
}

int getSatmeaaageHeaderGNSSSYS_SNR(int sys, char line[])
{
	int i, nsat;
	char SYS, result[4096] = "", data[128];
	switch (sys)
	{
	case SYS_GPS:nsat = NSATGPS; SYS = 'G'; break;
	case SYS_GLO:nsat = NSATGLO; SYS = 'R'; break;
	case SYS_GAL:nsat = NSATGAL; SYS = 'E'; break;
	case SYS_QZS:nsat = NSATQZS; SYS = 'J'; break;
	case SYS_CMP:nsat = NSATCMP; SYS = 'C'; break;
	case SYS_SBS:nsat = NSATSBS; SYS = 'S'; break;
	default:nsat = 3; SYS = 'N'; break;
	}
	for (i = 0; i<nsat; i++)
	{
		sprintf(data, "%7c%2.2d", SYS, i + 1);
		/*sprintf(data,"%58c%2d",SYS,i+1);*/
		strcat(result, data);
	}
	strcpy(line, result);
	return 1;
	//printf("%s",hehe);
	//printf("hehe\n");
}

extern void outsatsnr_single(FILE *fpSat_snr, rtk_t* rtk, obsd_t *obs, int n)
{
	int i, j;
	char line[2048];
	gtime_t time;
	/*FCB *myfcb;*/
	//rtk->opt.
	//if(!IsOpen)fpSat=fopen("E:\\learnprogram\\ReBuild_RTKLIB\\option\\SatStatis.txt","w");
	//if (!fpSat_snr)
	//{
	//	/*FILE *fp;
	//	char  fd = fp->fd;
	//	FCB *fcb;
	//	char *filiname = fcb[fd].filiname*/
	//	strcpy(line, outpath);
	//	fpSat_snr = fopen(strcat(line, "psu_snr"), "w");
	//}
	if (!IsWriteSNRHeader)
	{
		sprintf(line, "%23s", "%BDT           Obs_Time"); fputs(line, fpSat_snr);
		for (i = 0; i<6; i++)
		{
			line[0] = '\0';
			if (rtk->opt.navsys&(int)pow(2.0, i))
				getSatmeaaageHeaderGNSSSYS_SNR((int)pow(2.0, i), line);
			fputs(line, fpSat_snr);
		}
		fputs("\n", fpSat_snr);
		//rtk->opt.navsys
		IsWriteSNRHeader++;
	}

	time = rtk->sol.time;
	if (rtk->tsys == TSYS_CMP){
		time = gpst2bdt(rtk->sol.time);
	}
	time2str(time, line, 3); fputs(line, fpSat_snr);
	for (i = 1; i <= MAXSAT; i++)
	{
		if (!(rtk->opt.navsys&satsys(i, NULL)))continue;
		for (j = 0; j<n; j++)
		{
			if (obs[j].sat == i)break;
		}
		if (j >= n) sprintf(line, " %08.1f", 0.0);
		else sprintf(line, " %3.3d%05.1f", (int)(rtk->ssat[i - 1].snr[0] * 0.25 * 10), rtk->ssat[i - 1].azel[1] * R2D);
		fputs(line, fpSat_snr);
	}

	fputs("\n", fpSat_snr);
}

/* discard space characters at tail ------------------------------------------*/
/* 删除“#”及其之后的内容 */
/* 从字符串末尾开始，逐个删除不可见字符，直到遇到可见字符 */
static void chop(char *str)
{
	char *p;
	if ((p = strchr(str, '#'))) *p = '\0'; /* comment */
	for (p = str + strlen(str) - 1; p >= str&&!isgraph((int)*p); p--) *p = '\0';
}

/* load options ----------------------------------------------------------------
* load options from file
* args   : char   *file     I  options file
*          opt_t  *opts     IO options table
*                              (terminated with table[i].name="")
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int loadfiles(const char* file, opt_t* opts, char* infile[], char* outfile)
{
	FILE* fp;
	opt_t* opt;
	char buff[2048], * p;
	int n = 0, ni = 0;

	trace(3, "loadopts: file=%s\n", file);

	if (!(fp = fopen(file, "r"))) {
		trace(1, "loadopts: options file open error (%s)\n", file);
		return 0;
	}
	while (fgets(buff, sizeof(buff), fp)) {
		n++;
		chop(buff);
		if (buff[0] == '\0') continue;
		if (!(p = strstr(buff, "=")))	//查找buff中“=”的位置，p为返回的指针
		{
			fprintf(stderr, "invalid option %s (%s:%d)\n", buff, file, n);
			continue;
		}
		*p++ = '\0';		//将“=”后面的内容清空
		chop(buff);

		if (!strcmp(buff, "-infile")) {
			if (strlen(p) > 0)
			{
				strcpy(infile[ni++], p);
				continue;
			}
		}
		if (!strcmp(buff, "-outfile")) {
			strcpy(outfile, p);
			continue;
		}
	}
	fclose(fp);

	return ni;
}