/*
程序说明：
运行本工程的方法有两种：
1. 在工程目录下的 release\\批处理目录下，直接运行17个测试的批处理文件 RUN.bat
2. 打开工程编译运行工程文件，手动输入测试序号 1-17， 当输入 -1 时程序结束输入

对于每个测试，最短路径的解和运行时间都以追加方式写入文件源文件所在目录下的 out.txt文件中
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <time.h>
using namespace std;

// 常量
const int MAX = 60;
const int NumChromos = MAX * 5;
const int MaxNumGen = 5000;				//最大迭代次数为5000代
const int ConvergNum = 200;
const double DIFF = 1e-15;				//收敛判断代数
const double PC = 0.3;				//交叉概率
const double PM = 0.1;				//变异概率

// 结构体
struct City
{
	int id;//城市序号
	double x, y;//城市坐标
};
struct Chromo
{
	double sumOfDis;//染色体的值
	vector<int> vGene;//染色体的基因序列
};

// 容器
vector<City> vCity(MAX);//城市向量
vector<vector<double>> vCityDisMatr(MAX);//城市距离矩阵
vector<Chromo> vChromo(NumChromos + 1);	//第一个染色体只记录上一代遗传的最优解而不参与遗传运算
vector<Chromo> vTempChromo(NumChromos + 1);//用于种群拷贝
vector<double> Pro(NumChromos + 1);//生存概率

// 迭代函数 用于容器 for_each
void init_Group();// 初始化种群并且优化初值
void ChromoDis(Chromo& a);// 计算染色体的回路路径总长

// 打印函数
void print_Chromo(int cc); //打印染色体
void printCity(const City& city);//打印城市
void printResult();//打印最终解
void startMsg();
// 全局函数
double CalCircusDis(const vector<int>& v);//计算回路路径长度
void Evolution(int);//进化函数
bool ReadData();//读入文件数据
double RandProduce(double, double);//产生一定范围的随机数
void Crossover();// 交配
void Mutation();//变异
void Selection();//选择
void CalCityDis();//计算城市距离
void Solve();// GA算法解决TSP问题的总调用函数
void Fitness();

// 全局变量 ----------------
int mode;//输入模式
int nGen = 0;
int NumCities; // 城市数目
double preDis; // 上一代最优解
char tspFile[20]; // tsp文件名
char inputFile[20];//测试文件名

// main --------------------------------
int main()
{
	while (1)
	{
		startMsg();
		printf("程序开始...\n输入测试序号( 1~17 ) / (输入 -1 结束程序)： ");
		scanf("%s", &inputFile);
		if (strcmp(inputFile,"-1")==0)
			break;

		strcat(inputFile, ".txt");
		printf("读入数据...\n");
		if (!ReadData())
		{
			cerr << "---读入城市数据出错!\n";
			exit(1);
		}
		printf("OK-----------------------------\n\n");
		clock_t startTime = clock();

		Solve();

		double tt = static_cast<double>(clock() - startTime) / CLOCKS_PER_SEC;
		printf("\n共使用时间: %.8f 秒\n", tt);

		// 最小路径长和时间写入文件
		FILE *fp;
		fp = fopen("out.txt", "a");
		fprintf(fp, "%.8lf %.8f\n", vChromo[1].sumOfDis, tt);
	}

	return 0;
}

void startMsg()
{
	printf("**************************************************\n");
	printf("***********                          *************\n");
	printf("***********          TSP             *************\n");
	printf("***********                          *************\n");
	printf("**************************************************\n");
}

// 文件读取
bool ReadData()
{
	try
	{
		//----------------------------------------------
		ifstream fin;
		fin.open(inputFile);
		cout << "输入文件类型(1 或 2): ";
		//cin >> mode;
		fin >> mode;
		cout << mode << endl;
		cout << "输入文件名: ";
		//cin >> tspFile;
		fin >> tspFile;
		cout << tspFile << endl;
		fin.close();
		//----------------------------------------------

		fin.open(tspFile);
		fin >> NumCities;
		printf("城市数目： %d\n", NumCities);
		int index = 0;
		vCity.resize(NumCities);
		//----------------------------------------------

		if (mode == 1)
		{
			for (vector<City>::iterator i = vCity.begin(); i != vCity.end(); ++i)
			{
				(*i).id = ++index;
				fin >> (*i).x >> (*i).y;
			}
			fin.close();
			printf("计算城市距离矩阵...\n");
			CalCityDis();
			printf("OK-----------------------------\n\n");
		}
		else
		{
			vCityDisMatr.resize(NumCities);
			vector<vector<double>>::iterator it = vCityDisMatr.begin();
			for (; it != vCityDisMatr.end(); ++it)
				(*it).resize(NumCities);

			for (int i = 0; i < NumCities; ++i)
			{
				vCity[i].id = ++index;
				for (int j = 0; j < NumCities; ++j)
					fin >> vCityDisMatr[i][j];
			}
		}
	}
	catch (...)
	{
		return false;
	}
	return true;
}

void print_Chromo(int cc)
{
	printf("%d-> ", cc);
}

void printCity(const City& city)
{
	printf("City #%d: ( %d, %d )\n", city.id, city.x, city.y);
}

//打印最后结果
void printResult()
{
	printf("最优解回路： ");
	for_each(vChromo[1].vGene.begin(), vChromo[1].vGene.end(), print_Chromo);
	printf("1.\n最短距离= %.8lf\n", vChromo[1].sumOfDis);
}

// 染色体比较函数，用于sort排序
bool cmp_Chromo(const Chromo &a, const Chromo &b)
{
	return a.sumOfDis < b.sumOfDis;
}

// 迭代函数 用于容器 for_each
void ChromoDis(Chromo& a)
{
	a.sumOfDis = CalCircusDis(a.vGene);
}

void init_Chromo(Chromo& a)
{
	a.vGene.resize(NumCities);
	for (int i = 1; i <= NumCities; ++i)
		a.vGene[i - 1] = i;

	// random_shuffle(first, last):　Rearranges the elements in the range [first,last) randomly.
	// 重新打乱染色体基因排布顺序
	random_shuffle(a.vGene.begin(), a.vGene.end());
	a.sumOfDis = CalCircusDis(a.vGene);
}

// 初始化种群并且优化初值
void init_Group()
{
	printf("初始化 %d 个染色体: \n", NumChromos);
	for_each(vChromo.begin() + 1, vChromo.end(), init_Chromo);
	//先淘汰一次，以优化初值
	vChromo[0] = vChromo[1];
	sort(vChromo.begin() + 1, vChromo.end(), cmp_Chromo);
}

// 遗传算法函数
void Solve()
{
	Fitness();
	printf("初始化种群开始...\n");
	init_Group();
	printf("OK-----------------------------\n\n");
	printf("最大迭代次数5000，结束条件为连续 %d 代变化小于 %.15lf\n", ConvergNum, DIFF);
	printf("开始算法...\n");
	int g;
	for (g = 1; g <= MaxNumGen; ++g)
	{
		Selection();//选择函数		
		srand(static_cast<unsigned>(time(NULL)));//设定随机数种子
		Crossover();// 交叉原则
		Mutation();//变异原则
		Evolution(g);//评估进化原则

		if (ConvergNum < nGen)	break;
	}
	printf("进化结束，共遗传 %d 代\n\n", g);
	printResult();
}

// 产生随机数，在 [a,b]范围内
double RandProduce(double a, double b)
{
	double y;
	if (a>b) {
		printf("\nThe first parameter should be less than the second!\n");
		exit(1);
	}
	y = static_cast<double>(rand()) / RAND_MAX;
	return (a + (b - a)*y);
}

// 计算城市距离矩阵
void CalCityDis()
{
	double dx, dy;
	vCityDisMatr.resize(NumCities);
	vector<vector<double>>::iterator it = vCityDisMatr.begin();
	for (; it != vCityDisMatr.end(); ++it)
		(*it).resize(NumCities);

	for (int i = 0; i<NumCities; ++i)
		for (int j = i; j<NumCities; ++j)
			if (i == j) vCityDisMatr[i][j] = 0.0;
			else
			{
				dx = vCity[i].x - vCity[j].x;
				dy = vCity[i].y - vCity[j].y;
				vCityDisMatr[i][j] = sqrt(dx*dx + dy*dy);
				vCityDisMatr[j][i] = vCityDisMatr[i][j];
			}
}

// 计算回路的距离
double CalCircusDis(const vector<int>& Chm)
{
	double sum = 0.0;

	for (int it = 0; it<NumCities - 1; ++it)
		sum += vCityDisMatr[Chm[it] - 1][Chm[it + 1] - 1];

	sum += vCityDisMatr[Chm[NumCities - 1] - 1][Chm[0] - 1];
	return sum;
}

// 求适应函数
void Fitness()
{
	Pro[0] = 0.05;
	double a = 0.05;
	for (int i = 1; i <= NumChromos; ++i)
	{
		a *= 0.95;
		Pro[i] = Pro[i - 1] + a;
	}
}

//  选择算子，轮盘赌
void Selection()
{
	double r;
	int label;

	vTempChromo[0] = vChromo[0];
	for (int i = 1; i <= NumChromos; ++i)
	{
		r = RandProduce(0, Pro[NumChromos]);
		label = 0;
		for (int j = 0; j <= NumChromos; ++j)
		{
			if (r <= Pro[j])
			{
				label = j;
				break;
			}
		}
		vTempChromo[i] = vChromo[label];
	}
	swap(vChromo, vTempChromo);
}

// 变异，对两个随机位置之间的城市随机重排
void Mutation()
{
	vector<int>::iterator it;
	for (int i = 1; i <= NumChromos; ++i)
		if (PM > RandProduce(0, 1))
		{
			int left = static_cast<int>(RandProduce(0, NumCities / 2));
			int right = static_cast<int>(RandProduce(NumCities / 2, NumCities));
			it = vChromo[i].vGene.begin();

			random_shuffle(it + left, it + right);
		}
}

// 进化函数
void Evolution(int i)
{
	//对每个染色体求距离
	for_each(vChromo.begin() + 1, vChromo.end(), ChromoDis);
	//进化，好解(小)在前，上一代的最优解也带入
	sort(vChromo.begin(), vChromo.end(), cmp_Chromo);
	//printf("第 %d 代最短距离 = %.15f\n",i, vChromo[0].sumOfDis);
	preDis = vChromo[0].sumOfDis;
	if (preDis - vChromo[0].sumOfDis < DIFF) ++nGen;
	else nGen = 0;
}

//----------------------------较为优化的交配算法 -------------------------
//贪婪算法Right
int GreedyRight(vector<int> &genes, int index)
{
	bool Find = false;
	vector<int>::iterator it;
	for (it = genes.begin(); it != genes.end(); ++it)
		if (*it == index)
		{
			Find = true;
			break;
		}
	if (Find)
	{
		++it;
		if (it == genes.end())
			it = genes.begin();

		return *it;
	}
	else
		return -1;
}

//贪婪算法Left
int GreedyLeft(vector<int> &genes, int index)
{
	bool Find = false;
	vector<int>::iterator it;
	for (it = genes.begin(); it != genes.end(); ++it)
		if (*it == index)
		{
			Find = true;
			break;
		}

	if (Find)
	{
		if (it == genes.begin())
			return genes.back();
		else
		{
			--it;
			return *it;
		}
	}
	else
		return -1;
}

//贪婪消除
void GreedyErase(vector<int> &genes, int incre)
{
	bool Find = false;
	vector<int>::iterator it;
	for (it = genes.begin(); it != genes.end(); ++it)
		if (*it == incre)
		{
			Find = true;
			break;
		}
	if (Find)
		genes.erase(it);
}

// 交配算法：贪心交叉方式（Greedy Crossover），
// 具体算法可参见 谢胜利,等.求解TSP问题的一种改进的遗传算法
void CrossFunc(int A, int B)
{
	int randCity, cur, next, rightA, rightB, leftA, leftB;
	vector<int>  ChildA, ChildB;
	randCity = static_cast<int>(RandProduce(1, NumCities));
	cur = randCity;
	ChildA.push_back(cur);
	vector<int> ParentA = vChromo[A].vGene;
	vector<int> ParentB = vChromo[B].vGene;
	while (ParentA.size() > 1 && ParentB.size() > 1)
	{
		rightA = GreedyRight(ParentA, cur);
		rightB = GreedyRight(ParentB, cur);

		if (vCityDisMatr[cur - 1][rightA - 1]  < vCityDisMatr[cur - 1][rightB - 1])
		{
			ChildA.push_back(rightA);
			next = rightA;
		}
		else
		{
			ChildA.push_back(rightB);
			next = rightB;
		}
		GreedyErase(ParentA, cur);
		GreedyErase(ParentB, cur);
		cur = next;
	}
	cur = randCity;
	ChildB.push_back(cur);
	ParentA = vChromo[A].vGene;
	ParentB = vChromo[B].vGene;
	while (ParentA.size() > 1 && ParentB.size() > 1)
	{
		leftA = GreedyLeft(ParentA, cur);
		leftB = GreedyLeft(ParentB, cur);

		if (vCityDisMatr[cur - 1][leftA - 1] < vCityDisMatr[cur - 1][leftB - 1])
		{
			ChildB.push_back(leftA);
			next = leftA;
		}
		else
		{
			ChildB.push_back(leftB);
			next = leftB;
		}
		GreedyErase(ParentA, cur);
		GreedyErase(ParentB, cur);
		cur = next;
	}
	swap(vChromo[A].vGene, ChildA);
	swap(vChromo[B].vGene, ChildB);
}

//交配函数
void Crossover()
{
	vector<int> CrossGene;
	double pRand;

	for (int i = 1; i <= NumChromos; i++)
	{
		pRand = static_cast<double>(RandProduce(0, 1));
		if (pRand < PC)
			CrossGene.push_back(i);
	}

	int CN = CrossGene.size();
	if (CN % 2 != 0)
		CrossGene.pop_back();

	CN = CrossGene.size();
	for (int i = 0; i<CN; i += 2)
	{
		int A = CrossGene[i];
		int B = CrossGene[i + 1];

		CrossFunc(A, B);
	}
}
