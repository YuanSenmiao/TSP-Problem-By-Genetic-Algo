/*
����˵����
���б����̵ķ��������֣�
1. �ڹ���Ŀ¼�µ� release\\������Ŀ¼�£�ֱ������17�����Ե��������ļ� RUN.bat
2. �򿪹��̱������й����ļ����ֶ����������� 1-17�� ������ -1 ʱ�����������

����ÿ�����ԣ����·���Ľ������ʱ�䶼��׷�ӷ�ʽд���ļ�Դ�ļ�����Ŀ¼�µ� out.txt�ļ���
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

// ����
const int MAX = 60;
const int NumChromos = MAX * 5;
const int MaxNumGen = 5000;				//����������Ϊ5000��
const int ConvergNum = 200;
const double DIFF = 1e-15;				//�����жϴ���
const double PC = 0.3;				//�������
const double PM = 0.1;				//�������

// �ṹ��
struct City
{
	int id;//�������
	double x, y;//��������
};
struct Chromo
{
	double sumOfDis;//Ⱦɫ���ֵ
	vector<int> vGene;//Ⱦɫ��Ļ�������
};

// ����
vector<City> vCity(MAX);//��������
vector<vector<double>> vCityDisMatr(MAX);//���о������
vector<Chromo> vChromo(NumChromos + 1);	//��һ��Ⱦɫ��ֻ��¼��һ���Ŵ������Ž���������Ŵ�����
vector<Chromo> vTempChromo(NumChromos + 1);//������Ⱥ����
vector<double> Pro(NumChromos + 1);//�������

// �������� �������� for_each
void init_Group();// ��ʼ����Ⱥ�����Ż���ֵ
void ChromoDis(Chromo& a);// ����Ⱦɫ��Ļ�··���ܳ�

// ��ӡ����
void print_Chromo(int cc); //��ӡȾɫ��
void printCity(const City& city);//��ӡ����
void printResult();//��ӡ���ս�
void startMsg();
// ȫ�ֺ���
double CalCircusDis(const vector<int>& v);//�����··������
void Evolution(int);//��������
bool ReadData();//�����ļ�����
double RandProduce(double, double);//����һ����Χ�������
void Crossover();// ����
void Mutation();//����
void Selection();//ѡ��
void CalCityDis();//������о���
void Solve();// GA�㷨���TSP������ܵ��ú���
void Fitness();

// ȫ�ֱ��� ----------------
int mode;//����ģʽ
int nGen = 0;
int NumCities; // ������Ŀ
double preDis; // ��һ�����Ž�
char tspFile[20]; // tsp�ļ���
char inputFile[20];//�����ļ���

// main --------------------------------
int main()
{
	while (1)
	{
		startMsg();
		printf("����ʼ...\n����������( 1~17 ) / (���� -1 ��������)�� ");
		scanf("%s", &inputFile);
		if (strcmp(inputFile,"-1")==0)
			break;

		strcat(inputFile, ".txt");
		printf("��������...\n");
		if (!ReadData())
		{
			cerr << "---����������ݳ���!\n";
			exit(1);
		}
		printf("OK-----------------------------\n\n");
		clock_t startTime = clock();

		Solve();

		double tt = static_cast<double>(clock() - startTime) / CLOCKS_PER_SEC;
		printf("\n��ʹ��ʱ��: %.8f ��\n", tt);

		// ��С·������ʱ��д���ļ�
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

// �ļ���ȡ
bool ReadData()
{
	try
	{
		//----------------------------------------------
		ifstream fin;
		fin.open(inputFile);
		cout << "�����ļ�����(1 �� 2): ";
		//cin >> mode;
		fin >> mode;
		cout << mode << endl;
		cout << "�����ļ���: ";
		//cin >> tspFile;
		fin >> tspFile;
		cout << tspFile << endl;
		fin.close();
		//----------------------------------------------

		fin.open(tspFile);
		fin >> NumCities;
		printf("������Ŀ�� %d\n", NumCities);
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
			printf("������о������...\n");
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

//��ӡ�����
void printResult()
{
	printf("���Ž��·�� ");
	for_each(vChromo[1].vGene.begin(), vChromo[1].vGene.end(), print_Chromo);
	printf("1.\n��̾���= %.8lf\n", vChromo[1].sumOfDis);
}

// Ⱦɫ��ȽϺ���������sort����
bool cmp_Chromo(const Chromo &a, const Chromo &b)
{
	return a.sumOfDis < b.sumOfDis;
}

// �������� �������� for_each
void ChromoDis(Chromo& a)
{
	a.sumOfDis = CalCircusDis(a.vGene);
}

void init_Chromo(Chromo& a)
{
	a.vGene.resize(NumCities);
	for (int i = 1; i <= NumCities; ++i)
		a.vGene[i - 1] = i;

	// random_shuffle(first, last):��Rearranges the elements in the range [first,last) randomly.
	// ���´���Ⱦɫ������Ų�˳��
	random_shuffle(a.vGene.begin(), a.vGene.end());
	a.sumOfDis = CalCircusDis(a.vGene);
}

// ��ʼ����Ⱥ�����Ż���ֵ
void init_Group()
{
	printf("��ʼ�� %d ��Ⱦɫ��: \n", NumChromos);
	for_each(vChromo.begin() + 1, vChromo.end(), init_Chromo);
	//����̭һ�Σ����Ż���ֵ
	vChromo[0] = vChromo[1];
	sort(vChromo.begin() + 1, vChromo.end(), cmp_Chromo);
}

// �Ŵ��㷨����
void Solve()
{
	Fitness();
	printf("��ʼ����Ⱥ��ʼ...\n");
	init_Group();
	printf("OK-----------------------------\n\n");
	printf("����������5000����������Ϊ���� %d ���仯С�� %.15lf\n", ConvergNum, DIFF);
	printf("��ʼ�㷨...\n");
	int g;
	for (g = 1; g <= MaxNumGen; ++g)
	{
		Selection();//ѡ����		
		srand(static_cast<unsigned>(time(NULL)));//�趨���������
		Crossover();// ����ԭ��
		Mutation();//����ԭ��
		Evolution(g);//��������ԭ��

		if (ConvergNum < nGen)	break;
	}
	printf("�������������Ŵ� %d ��\n\n", g);
	printResult();
}

// ������������� [a,b]��Χ��
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

// ������о������
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

// �����·�ľ���
double CalCircusDis(const vector<int>& Chm)
{
	double sum = 0.0;

	for (int it = 0; it<NumCities - 1; ++it)
		sum += vCityDisMatr[Chm[it] - 1][Chm[it + 1] - 1];

	sum += vCityDisMatr[Chm[NumCities - 1] - 1][Chm[0] - 1];
	return sum;
}

// ����Ӧ����
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

//  ѡ�����ӣ����̶�
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

// ���죬���������λ��֮��ĳ����������
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

// ��������
void Evolution(int i)
{
	//��ÿ��Ⱦɫ�������
	for_each(vChromo.begin() + 1, vChromo.end(), ChromoDis);
	//�������ý�(С)��ǰ����һ�������Ž�Ҳ����
	sort(vChromo.begin(), vChromo.end(), cmp_Chromo);
	//printf("�� %d ����̾��� = %.15f\n",i, vChromo[0].sumOfDis);
	preDis = vChromo[0].sumOfDis;
	if (preDis - vChromo[0].sumOfDis < DIFF) ++nGen;
	else nGen = 0;
}

//----------------------------��Ϊ�Ż��Ľ����㷨 -------------------------
//̰���㷨Right
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

//̰���㷨Left
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

//̰������
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

// �����㷨��̰�Ľ��淽ʽ��Greedy Crossover����
// �����㷨�ɲμ� лʤ��,��.���TSP�����һ�ָĽ����Ŵ��㷨
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

//���亯��
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
