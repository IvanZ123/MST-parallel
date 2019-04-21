#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <math.h>
using namespace std;

struct Edge
{
	int from, to;
	double weight;
	//Edge(int from, int to, double weight):from(from), to(to), weight(weight){}
};

//Count the number of lines in a file, https://blog.csdn.net/wangshihui512/article/details/8921929 
//Can use m instead...
int CountLines(char *filename)
{
    ifstream ReadFile;
    int n=0;
    string tmp;
    ReadFile.open(filename,ios::in);//ios::in read only
    if(ReadFile.fail())
    {
        return 0;
    }
    else
    {
        while(getline(ReadFile,tmp,'\n'))
        {
            n++;
        }
        ReadFile.close();
        return n;
    }
}


void Par_Simulate_Priority_CW_using_Binary_Search(int n, int m, Edge* E,int* R )
{
	//int B[n],l[n],h[n],lo[n],hi[n],md[n];
	int* B = new int[n];
	int* l = new int[n];
	int* h = new int[n];
	int* lo = new int[n];
	int* hi = new int[n];
	int* md = new int[n];
	#pragma omp parallel for 
	for(int i = 0; i < n; i++)
	{
		//cout<<"omp initialize l,h"<<endl;
		l[i] = 0;
		h[i] = m-1;
	}
	for(int k = 0; k <1 + log2(m); k++)
	{
		#pragma omp parallel for
		for (int i = 0; i < n; ++i)
		{
			//cout<<"omp initialize B"<<endl;
			B[i] = 0;
			lo[i] = l[i];
			hi[i] = h[i];
		}
		#pragma omp parallel for
		for (int i = 0; i < m; ++i)
		{
			int u = E[i].from;
			if(u>n) continue;
			md[u] = floor((lo[u]+hi[u])/2);
			if(i >= lo[u] && i <= md[u]) B[u] = 1;
		}
		#pragma omp parallel for
		for (int i = 0; i < m; ++i)
		{
			int u = E[i].from;
			if(u>n) continue;
			md[u] = floor((lo[u]+hi[u])/2);
			if(B[u] == 1 && i >= lo[u] && i <= md[u]) h[u] = md[u];
			else if(B[u] == 0 && i > md[u] && i<= hi[u]) l[u] = md[u] + 1;
		}
	}
	#pragma omp parallel for
	for(int i = 0; i < m; ++i)
	{
		int u = E[i].from;
		//cout<<"u = "<<u<<endl;s
		if(u>n) continue;
		if(l[u] == i) 
		{
			R[u] = i;
		}
		//cout<<"Assign R["<<i<<"]"<<R[u]<<endl;
	}
	delete []B;
	delete []l;
	delete []h;
	delete []lo;
	delete []hi;
	delete []md;
}


void Par_Randomized_MST_Priority_CW (int n, int m, Edge *E, int *MST)
{
	int L[n];
	int *C = new int[n];
	int *R = new int[n];
	bool F;
	#pragma omp parallel for
	for(int i =0; i < n; i++)
	{
		L[i] = i;
	}
	F = (m > 0) ? true:false;
	while(F)
	{
		//toss a coin
		#pragma omp parallel for 
		for(int i = 0; i < n; i++)
		{
			C[i] = rand()%2;
			//cout<<"coin "<<C[i]<<endl;
		}
		//Use radix sort to get R
		//R[i] is the position that n+1 vertext first appear in the sorted E
		Par_Simulate_Priority_CW_using_Binary_Search(n, m, E, R);
		//cout<<"R value = "<<R[0]<<endl;
		#pragma omp parallel for 
		for(int i = 0; i < m; i++)
		{
			
			int u = E[i].from;
			int v = E[i].to;
			if(C[u]==0 && C[v] == 1 && R[u] == i && u < n && v < n )
			{
				L[u] = v;
				MST[i] = 1;
			}
		}

		#pragma omp parallel for 
		for(int i = 0; i < m; i++)
		{
			int uu = E[i].from;
			int vv = E[i].to;
			if(uu<n && vv<n && L[uu] != L[vv]){
				E[i].from = L[uu];
				E[i].to = L[vv];
			}
			else
			{
				E[i].from = n+1;
				E[i].to = n+1;
			}

		}
		F = false;
		#pragma omp parallel for 
		for(int i = 0; i < m; i++)
		{
			
			if(E[i].from != E[i].to) {
				F = true;
			}
		}
	}
	delete []R;
	delete []C;
}



int main(){
	char filename[512] = "in1.txt";
	fstream fin;
	string tmp;
	int numL,n,m;
	fin.open(filename,ios::in);
	if(fin.fail()) 
	{
		cout << "File not exit!" <<endl;
		return 0;
	}
	else
	{
		numL = CountLines(filename);
		fin >> n;	
		fin >> m;
		cout << "n="<<n<<" m="<<m<<endl;
		double edge[numL-1][3];
		// vector<vector<double> >E;
		Edge* EG = new Edge[numL-1];
		int i = 0;
		while(!fin.eof())
		{
			fin>>edge[i][0];
			fin>>edge[i][1];
			fin>>edge[i][2];
			EG[i].from = edge[i][0];
			EG[i].to = edge[i][1];
			EG[i].weight = edge[i][2];
			i++;
		}
		fin.close();
		// cout << "edge 0: "<<EG[0].from << " "<<EG[0].to << " "<<EG[0].weight << " "<<endl;
		// cout << "edge 1: "<<EG[1].from << " "<<EG[1].to << " "<<EG[1].weight << " "<<endl;
		// cout << "edge 2: "<<EG[2].from << " "<<EG[2].to << " "<<EG[].weight << " "<<endl;
		int* MST = new int[m];
		#pragma omp parallel for 
		for (int i = 0; i < m; ++i)
		{
			MST[i] = 0;
		}
		for (int i = 0; i < m; ++i)
		{
			cout <<"MST["<<i<<"] = "<<MST[i]<<endl;
		}
		Par_Randomized_MST_Priority_CW(n,m,EG,MST);

		for (int i = 0; i < m; ++i)
		{
			cout <<"MST["<<i<<"] = "<<MST[i]<<" ";
		}

	}

}
