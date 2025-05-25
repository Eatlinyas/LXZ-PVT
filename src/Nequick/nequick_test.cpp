#include <stdio.h>
#include <windows.h>


extern "C" {
	extern double  ionmodel_nequick(int month, double UT, double usr_longitude_degree, double usr_latitude_degree, double usr_height_meters,
		double sat_longitude_degree, double sat_latitude_degree, double sat_height_meters);
	extern void test(char* string);
}
void main() 
{
	//1 20 307.19 5.25 -25.76 10.94 44.72 20450566.19 381.63445
	int month = 1;
	double UTC = 20.0;
	double usr_longitude_degree = 307.19;
	double usr_latitude_degree = 5.25;
	double usr_height_meters = -25.76;

	double sat_longitude_degree = 10.94;
	double sat_latitude_degree = 44.72;
	double sat_height_meters = 20450566.19;
	
	double tec=0.0;

	tec = ionmodel_nequick(month, UTC, usr_longitude_degree, usr_latitude_degree, usr_height_meters,sat_longitude_degree, sat_latitude_degree, sat_height_meters);


	printf("%lf \n", tec);
	test("OK\n");
	system("pause");

}