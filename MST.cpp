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
	int* B = new int[n];
	int* l = new int[n];
	int* h = new int[n];
	int* lo = new int[n];
	int* hi = new int[n];
	int* md = new int[n];
	//l[n],h[n],lo[n],hi[n],md[n];
	#pragma omp parallel for 
	for(int i = 0; i < n; i++)
	{
		l[i] = 0;
		h[i] = m-1;
	}
	for(int k = 0; k <1 + log2(m); k++)
	{
		#pragma omp parallel for
		for (int i = 0; i < n; ++i)
		{
			B[i] = 0;
			lo[i] = l[i];
			hi[i] = h[i];
		}
		
		#pragma omp parallel for
		for (int i = 0; i < m; ++i)
		{
			int u = E[i].from;
			if(u>=n) 
			{
				continue;
			}

			md[u] = floor((lo[u]+hi[u])/2);
			
			if(i >= lo[u] && i <= md[u])
			{
				B[u] = 1;
			}
		}
		#pragma omp parallel for
		for (int i = 0; i < m; ++i)
		{
			int u = E[i].from;
			if(u>=n) continue;
			md[u] = floor((lo[u]+hi[u])/2);
			if(B[u] == 1 && i >= lo[u] && i <= md[u]) h[u] = md[u];
			else if(B[u] == 0 && i > md[u] && i<= hi[u]) l[u] = md[u] + 1;
		}
	}
	#pragma omp parallel for
	for(int i = 0; i < m; ++i)
	{
		int u = E[i].from;
		if(u>=n) continue;
		else if(l[u] == i) 
		{
			R[u] = i;
		}
	}
	delete []B;
	delete []l;
	delete []lo;
	delete []h;
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
	

	#pragma omp parallel for
	for(int i =0; i < n; i++)
	{
		R[i] = m;
	}
	F = (m > 0) ? true:false;
	while(F)
	{
		//toss a coin
		#pragma omp parallel for 
		for(int i = 0; i < n; i++)
		{
			C[i] = rand()%2;
		}
		//Use radix sort to get R
		//R[i] is the position that n+1 vertext first appear in the sorted E

		Par_Simulate_Priority_CW_using_Binary_Search(n, m, E, R);

		//cout<<"R is"<<R[0]<<" "<<R[1]<<" "<<R[2]<<" "<<R[3]<<" "<<R[4]<<endl;
		#pragma omp parallel for 
		for(int i = 0; i < m; i++)
		{
			
			int u = E[i].from;
			int v = E[i].to;
			if(C[u]==1 && C[v] == 0 && R[u] == i && u < n && v < n )
			{
				L[v] = u;
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
				E[i].from = n;
				E[i].to = n;
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

bool sortweight( const vector<double>& v1,
				const vector<double>& v2 ) 
		{ 
			return v1[2] < v2[2];
		}
bool sortfrom( const vector<double>& v1,
				const vector<double>& v2 ) 
		{ 
			return v1[0] < v2[0];
		}

int main(){
	char filename[512] = "in.txt";
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
		Edge* backup = new Edge[numL-1];
		int i = 0;
		vector<vector<double>> unsort;
		while(!fin.eof())
		{
			fin>>edge[i][0];
			fin>>edge[i][1];
			fin>>edge[i][2];
			vector<double> tmp;
			tmp.push_back(edge[i][0]);
			tmp.push_back(edge[i][1]);
			tmp.push_back(edge[i][2]);
			unsort.push_back(tmp);
			i++;
		}
		fin.close();
		//https://www.geeksforgeeks.org/sorting-2d-vector-in-c-set-1-by-row-and-column/
		std::sort(unsort.begin(), unsort.end(),sortweight);

		//#pragma omp parallel for 
		for (int i = 0; i < m; ++i)
		{
			int j = i+1;
			EG[i].from = unsort[j][0]-1;
			// cout<<"EG "<<EG[i].from<<endl;
			EG[i].to = unsort[j][1]-1;
			EG[i].weight = unsort[j][2];
			backup[i].from = unsort[j][0];
			backup[i].to = unsort[j][1];
			backup[i].weight = unsort[j][2];
		}
		//cout<<"EG "<<EG[1].from<<endl;

		int* MST = new int[m];
		#pragma omp parallel for 
		for (int i = 0; i < m; ++i)
		{
			MST[i] = 0;
		}

		Par_Randomized_MST_Priority_CW(n,m,EG,MST);

		// for (int i = 0; i < m; ++i)
		// {
		// 	cout <<"MST["<<i<<"] = "<<MST[i];
		// 	cout << " "<<endl;
		// }
		vector<vector<double>> ans;
		int mout = 0;
		double cost = 0;
		for (int i = 0; i < m; ++i)
		{
			if(MST[i]==1)
			{
				vector<double> tmp;
				tmp.push_back(backup[i].from);
				tmp.push_back(backup[i].to);
				tmp.push_back(backup[i].weight);
				ans.push_back(tmp);
				mout++;
				cost += backup[i].weight;
			}
		}
		std:sort(ans.begin(), ans.end(),sortfrom);
		// for (int i = 0; i < 20; ++i)
		// {
		// 	cout <<ans[i][0]<<" "<<ans[i][1]<<" "<<ans[i][2]<<endl;
		// }
		cout << "num edges output: "<< mout<<endl;
		cout << "TOTAL cost output: "<< cost<<endl;

		delete []EG;
		delete []backup;

	}

}
