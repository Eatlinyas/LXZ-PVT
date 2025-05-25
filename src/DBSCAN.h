#if !defined(_DBSCAN_H_INCLUDED_)  
#define _DBSCAN_H_INCLUDED_  

#include <iostream>  
#include <sstream>  
#include <fstream>  
#include <vector>  
#include <ctime>  
#include <cstdlib>  
#include <limits>  
#include <cmath>  
#include <stack>  
#include <map>  

using namespace std;

enum
{
	pointType_UNDO,
	pointType_NOISE,
	pointType_BORDER,
	pointType_CORE

};

class point{
public:
	float x;
	float y;
	vector<float> xn;
	//vector<float> yn;  
	int cluster;
	int pointType;  //1 noise 2 border 3 core  
	int pts;        //points in MinPts   
	int corePointID;
	vector<int> corepts;
	int  visited;
	void init();
	point();
	point(float a, float b, int c){
		x = a;
		y = b;
		cluster = c;
	};
	point(vector<float>  an, int c){//three or more dimension  
		xn = an;
		cluster = c;
	};
};


float stringToFloat(string i);
vector<point> openFile(const char* dataset);
float squareDistance(point a, point b);
float squareDistanceVect(point a, point b);
double  DBSCAN(vector<point> dataset, float Eps, int MinPts);


#endif  