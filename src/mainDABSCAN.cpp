#if 1  
//代码原始出处：http://blog.csdn.net/k76853/article/details/50440182   
//zhouxianen,20170322  
/*修改：
（1）维度扩展，原来程序只能读取、聚类和显示二维特征数据，扩展后，可以处理n维数据，(2)添加了枚举类型，
（3）定义了init()函数初始化化类成员变量(效率较低)，
（4）领域半径eps判断条件修改为小于等于，这是为保证与matlab版本形式一致，
（5)输出所有数据所属类型，结果输出到控制台显示，同时保存为clutering.txt文件，其中文件的读取和保存格式与原作者一致。
修改后代码运行结果与matlab官网代码（http://cn.mathworks.com/matlabcentral/fileexchange/52905-dbscan-clustering-algorithm）运行结果一致。
*/
/*
输入数据类型：*.txt文件

输入txt文件格式如下(数据维度可以是n维)：
1.400000,0.200000
1.400000,0.200000
1.300000,0.200000

输出结果：控制台输出，同时输出clustering.txt文件
输出txt文件格式如下：
1.4,0.2,1
1.4,0.2,1
1.3,0.2,1
其中第三列为分类标记，数据之间用逗号分开
本文件中附上了两组组测试数据和测试结果（测试数据来源：http://yarpiz.com/255/ypml110-dbscan-clustering）
二维测试数据：dataPoint2.txt
四维测试数据：dataPoint4.txt
输出结果：clustering2.txt,clustering4.txt (eps=0.25,MinPts=3时结果)
*/

#include "DBSCAN.h"  

int main(int argc, char** argv) {
	vector<point> dataset = openFile("matlab2.txt");
	DBSCAN(dataset, 0.02, 6);
	return 0;
}

#endif  