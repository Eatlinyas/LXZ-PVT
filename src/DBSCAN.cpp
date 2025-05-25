#include "DBSCAN.h"  

int clusterID = 0;

void point::init()
{
	cluster = 0;
	pointType = pointType_UNDO;//pointType_NOISE pointType_UNDO  
	pts = 0;
	visited = 0;
	corePointID = -1;
}


float stringToFloat(string i){
	stringstream sf;
	float score = 0;
	sf << i;
	sf >> score;
	return score;
}

vector<point> openFile(const char* dataset){
	fstream file;
	file.open(dataset, ios::in);
	if (!file)
	{
		cout << "Open File Failed!" << endl;
		vector<point> a;
		return a;
	}
	vector<point> data;

	int i = 1;
	while (!file.eof()){
		string temp;
		file >> temp; //read one row data  
		int split = temp.find(',', 0);

		//point p(stringToFloat(temp.substr(0,split)),stringToFloat(temp.substr(split+1,temp.length()-1)),i++);  

		//three or more dimension data process  
		int numSpit = 0;
		vector <int> splitN1, splitN2;
		numSpit = count(temp.begin(), temp.end(), ',');
		if (numSpit>0){ // 2d or 3d 4d....
			vector<float> xn;
			int m = 0; int cPos = 0;
			while (1){//  
				int splitTemp = temp.find_first_of(',', cPos);
				if (splitTemp != string::npos){
					m++;
					splitN1.push_back(cPos);
					cPos = splitTemp + 1;
					splitN2.push_back(splitTemp);
				}
				else
					break;
			}
			for (int m = 0; m<numSpit; m++){
				xn.push_back(stringToFloat(temp.substr(splitN1[m], splitN2[m])));
			}
			xn.push_back(stringToFloat(temp.substr(splitN2[numSpit - 1] + 1, temp.length() - 1)));
			point p(xn, i++);
			data.push_back(p);
		}
		else
		{
			vector<float> xn;
			int m = 0; int cPos = 0;
			xn.push_back(stringToFloat(temp));
			point p(xn, i++);
			data.push_back(p);
		}
	}
	file.close();
	cout << "successful!" << endl;
	return data;
}

float squareDistance(point a, point b){
	return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}

float squareDistanceVect(point a, point b){
	vector <float> xn1 = a.xn;
	vector <float> yn1 = b.xn;

	float sumSqrt = 0;
	for (int i = 0; i<xn1.size(); i++){
		sumSqrt = sumSqrt + (xn1[i] - yn1[i])*(xn1[i] - yn1[i]);
	}
	return sqrt(sumSqrt);
}

double  DBSCAN(vector<point> dataset, float Eps, int MinPts){
	double mean_t = 0.0;
	int n = 0;
	int len = dataset.size();//数据长度  
	if (len <= 0)return mean_t;

	for (int i = 0; i<len; i++)//参数初始化  
	{
		dataset[i].init();
	}
	clusterID = 0;
	vector<vector <float>> distP2P(len);

	//calculate pts  
	//cout << "calculate pts" << endl;
	for (int i = 0; i<len; i++){
		for (int j = 0; j<len; j++){//i+1  
			float distance = squareDistanceVect(dataset[i], dataset[j]);//squareDistanceVect squareDistance   
			distP2P[i].push_back(distance);//disp for debug  
			if (distance <= Eps){
				dataset[i].pts++;
			}
		}
	}
	//core point   
	//cout << "core point " << endl;
	vector<point> corePoint;
	for (int i = 0; i<len; i++){
		int tempPts = dataset[i].pts;
		if (tempPts >= MinPts) {
			dataset[i].pointType = pointType_CORE;
			dataset[i].corePointID = i;
			corePoint.push_back(dataset[i]);
		}
	}
	//cout << "joint core point" << endl;

	//joint core point  
	int numCorePoint = corePoint.size(); //core point number  

	for (int i = 0; i<numCorePoint; i++){
		for (int j = 0; j<numCorePoint; j++){
			float distTemp = distP2P[corePoint[i].corePointID][corePoint[j].corePointID];//display for debug  
			if (distTemp <= Eps){//squareDistance(corePoint[i],corePoint[j])  
				corePoint[i].corepts.push_back(j);//other point orderID link to core point  
			}
		}
	}
	for (int i = 0; i<numCorePoint; i++){
		stack<point*> ps;
		if (corePoint[i].visited == 1) continue;
		clusterID++;
		corePoint[i].cluster = clusterID; //create a new cluster  
		ps.push(&corePoint[i]);
		point *v;
		while (!ps.empty()){
			v = ps.top();
			v->visited = 1;
			ps.pop();
			for (int j = 0; j<v->corepts.size(); j++){
				if (corePoint[v->corepts[j]].visited == 1) continue;
				corePoint[v->corepts[j]].cluster = corePoint[i].cluster;
				corePoint[v->corepts[j]].visited = 1;
				ps.push(&corePoint[v->corepts[j]]);
			}
		}

	}

	//cout << "border point,joint border point to core point" << endl;
	//border point,joint border point to core point  
	for (int i = 0; i<len; i++){
		for (int j = 0; j<numCorePoint; j++){
			float distTemp = distP2P[i][corePoint[j].corePointID];
			if (distTemp <= Eps) {
				dataset[i].pointType = pointType_BORDER;
				dataset[i].cluster = corePoint[j].cluster;
				break;
			}
		}
	}
	//cout << "output" << endl;
	//output  
	//display    
	//save in .txt format named clustering.txt  
	//fstream clustering;//save .txt  
	//clustering.open("clustering.txt", ios::out);//save .txt  

	char dispInfo[500];
	int dataDim = dataset[0].xn.size();//data dimension  
	//for (int i = 0; i<len; i++){
	//	//%11.4lf,%11.4lf,%11.4lf,%11.4lf  
	//	sprintf(dispInfo, "第%3d个数据：", i + 1);
	//	for (int j = 0; j<dataDim; j++){
	//		char dataSrc[30];
	//		if (j == 0)
	//			sprintf(dataSrc, "%11.4lf", dataset[i].xn[j]);
	//		else
	//			sprintf(dataSrc, ",%11.4lf", dataset[i].xn[j]);
	//		strcat(dispInfo, dataSrc);
	//	}
	//	char dataClust[30];
	//	sprintf(dataClust, ",%4d\n", dataset[i].cluster);
	//	strcat(dispInfo, dataClust);
	//	//trace(5, "%s\n", dispInfo);
	//	//cout << dispInfo;      //display in cmd window  
	//	//clustering << dispInfo;//save the results in .txt format named clustering.txt  

	//}
	//clustering.close();//save .txt  

	for (int i = 0; i < dataset.size(); i++)
	{
		//printf("%d: %d \n", i,dataset[i].cluster);
		if (dataset[i].cluster != 1)continue;
		mean_t += (dataset[i].xn[0] - mean_t) / ++n;
	}
	return mean_t;
}