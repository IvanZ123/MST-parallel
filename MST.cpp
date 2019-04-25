#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <unordered_map>
#include <math.h>
using namespace std;

void bubblesort(vector<double> &num_list, int left, int right){
	bool swapped;
	int n = right - left + 1;
	for(int i = left; i < n - 1; i++){ 
		swapped = false; 
		for (int j = left; j < n - i - 1; j++){ 
			if (num_list[j] > num_list[j+1]){
			 	double temp = num_list[j];
			 	num_list[j + 1] = num_list[j];
			 	num_list[j] = temp;
			   	swapped = true; 
			} 
		} 
		if(!swapped){
			break;
		}      
	}
}

void prefix_sum(vector<double> &num_list){
	int size = num_list.size();
	if(size <= 1){
		return;
	}
	vector<double> left(size/2);

	#pragma omp parallel for
	for(int i = 0; i < size/2; i++){
		left[i] = num_list[2*i] + num_list[2*i + 1];
	}

	prefix_sum(left);

	#pragma omp parallel for
	for(int i = 1; i < size; i++){
		if(i % 2 != 0){
			num_list[i] = left[i/2];
		}
		else{
			num_list[i] = left[(i - 1)/2] + num_list[i];
		}
	}
}

int partition(vector<double> &num_list, int left, int right, int partition_val){
	int partition_index;
	int size = right - left + 1;
	int B[size];
	vector<double> less_than(size);
	vector<double> greater_than(size);

	#pragma omp parallel for
	for(int i = 0; i < size; i++){
		B[i] = num_list[left + i];
		if(B[i] < partition_val){
			less_than[i] = 1;
		}
		else if(B[i] > partition_val){
			greater_than[i] = 1;
		}
	}

	prefix_sum(less_than);
	prefix_sum(greater_than);

	partition_index = left + less_than[size - 1];
	num_list[partition_index] = partition_val;

	#pragma omp parallel for
	for(int i = 0; i < size; i++){
		if(B[i] < partition_val){
			num_list[left + less_than[i] - 1] = B[i];
		}
		else if(B[i] > partition_val){
			num_list[partition_index + greater_than[i]] = B[i];
		}
	}

	return partition_index;
}

void quicksort(vector<double> &num_list, int left, int right, int m){
	int size = right - left + 1;
	//perform bubble sort for small size
	if(size <= m){
		//bubblesort(num_list, left, right);

		sort(num_list.begin() + left, num_list.begin() + right + 1);
	}
	else{
		int random_index = rand() % size; 
		int partition_val = num_list[left + random_index];

		int partition_index = partition(num_list, left, right, partition_val);

		#pragma omp task shared(num_list)
		//left part
		quicksort(num_list, left, partition_index - 1, m);
		//right part
		quicksort(num_list, partition_index + 1, right, m);

		#pragma omp taskwait
	}
	return;
}




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


void Par_Simulate_Priority_CW_using_Binary_Search(long int n,long int m, Edge* E,int* R )
{
	int* B = new int[n];
	int* l = new int[n];
	int* h = new int[n];
	int* lo = new int[n];
	int* hi = new int[n];
	int* md = new int[n];
	#pragma omp parallel for 
	for(int i = 0; i < n; i++)
	{
		l[i] = 0;
		h[i] = m-1;
	}
	for(int k = 0; k <=1 + log2(m); k++)
	{
		#pragma omp parallel for
		for (int i = 0; i < n; ++i)
		{
			B[i] = 0;
			// lo[i] = l[i];
			// hi[i] = h[i];
		}
		
		// #pragma omp parallel for
		// for (int i = 0; i < m; ++i)
		// {
		// 	int u = E[i].from;
		// 	if(u>=n) continue;
		// 	md[u] = floor((lo[u]+hi[u])/2);
			
		// 	if(i >= lo[u] && i <= md[u])
		// 	{
		// 		B[u] = 1;
		// 	}
		// }
		// #pragma omp parallel for
		// for (int i = 0; i < m; ++i)
		// {
		// 	int u = E[i].from;
		// 	if(u>=n) continue;
		// 	md[u] = floor((lo[u]+hi[u])/2);
		// 	if(B[u] == 1 && i >= lo[u] && i <= md[u]) h[u] = md[u];
		// 	else if(B[u] == 0 && i > md[u] && i<= hi[u]) l[u] = md[u] + 1;
		// }
		#pragma omp parallel for
		for (int i = 0; i < m; ++i)
		{
			int u = E[i].from;
			if(u>=n) continue;
			md[u] = floor((l[u]+h[u])/2);
			
			if(i >= l[u] && i <= md[u])
			{
				B[u] = 1;
			}
		}
		#pragma omp parallel for
		for (int i = 0; i < m; ++i)
		{
			int u = E[i].from;
			if(u>=n) continue;
			md[u] = floor((l[u]+h[u])/2);
			if(B[u] == 1 && i >= l[u] && i <= md[u]) h[u] = md[u];
			else if(B[u] == 0 && i > md[u] && i<= h[u]) l[u] = md[u] + 1;
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


void Par_Randomized_MST_Priority_CW (long int n,long int m, Edge *E, int *MST)
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
		// #pragma omp parallel for
		// for(int i =0; i < n; i++)
		// {
		// 	R[i] = m;
		// }	
		//toss a coin
		#pragma omp parallel for 
		for(int i = 0; i < n; i++)
		{
			C[i] = rand()%2;
		}
//Use radix sort to get R
//R[i] is the position that n+1 vertext first appear in the sorted E

		Par_Simulate_Priority_CW_using_Binary_Search(n, m, E, R);


		#pragma omp parallel for 
		for(int i = 0; i < m; i++)
		{
			
			int u = E[i].from;
			int v = E[i].to;
			if(C[u]==0 && C[v] == 1 && R[u] == i && u < n && v < n )
			// if(C[u]==0 && C[v] == 1 && R[u] == i)
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
			// if(L[uu] != L[vv]){
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

bool edgecomp(Edge e1, Edge e2)
{
	return e1.weight<e2.weight;
}

int main(int argc, char *argv[]){

	// char filename[512] = "com-dblp-in.txt";
	char filename[512] = "com-amazon-in.txt";
	// char filename[512] = "in.txt";
	omp_set_num_threads(68);
	fstream fin;
	long numL,n,m;
	fin.open(filename,ios::in);
    double start_time, final_time;
	if(fin.fail()) 
	{
		cout << "File not exit!" <<endl;
		return 0;
	}
	else
	{
		// numL = CountLines(filename);
		numL = CountLines(filename);
		fin >> n;	
		fin >> m;
		cout << "n="<<n<<" m="<<m<<endl;
		
		Edge* EG = new Edge[numL-1];
		Edge* backup = new Edge[numL-1];
		double (*edge)[3] = new double[m][3];

		int i = 0;
		cout<<"EG defined"<<endl;
		while(!fin.eof())
		{
			fin>>edge[i][0];
			fin>>edge[i][1];
			fin>>edge[i][2];
			i++;
		}
		fin.close();

		#pragma omp parallel for
		for (int i = 0; i < m; ++i)
		{
			EG[i].from = edge[i][0]-1;
			EG[i].to = edge[i][1]-1;
			EG[i].weight = edge[i][2];
			backup[i].from = edge[i][0];
			backup[i].to = edge[i][1];
			backup[i].weight = edge[i][2];
		}

		cout<<"data got"<<endl;
		sort(EG,EG+m,&edgecomp);
		sort(backup,backup+m,&edgecomp);
		cout<<"data sorted"<<endl;
		//https://www.geeksforgeeks.org/sorting-2d-vector-in-c-set-1-by-row-and-column/
		// std::sort(unsort.begin(), unsort.end(),sortweight);
		//use quick sort here instead of the sort. We want the sorted vetor<Edge> E as further input
		//quicksort(weight,0,weight.size()-1,32);
		int* MST = new int[m];
		#pragma omp parallel for 
		for (int i = 0; i < m; ++i)
		{
			MST[i] = 0;
		}
		/* starting the clock */
	    start_time = omp_get_wtime();

		cout<<"program start"<<endl;
		Par_Randomized_MST_Priority_CW(n,m,EG,MST);
		cout<<"program end"<<endl;
		final_time = omp_get_wtime() - start_time;

		printf("Total time: %.6f\n",final_time);
	
		vector<vector<double>> ans;
		long mout = 0;
		long double cost = 0;
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
		
		cout << "num edges output: "<< mout<<endl;
		cout << "TOTAL cost output: "<< cost<<endl;
 		ofstream output;
		output.open ("com-amazon-search-out.txt");
		// output.open ("out2.txt");
		output <<mout<<" "<<cost <<"\n";
		for (int i = 0; i < ans.size(); ++i)
		{
			output <<ans[i][0]<<" "<<ans[i][1]<<" "<<ans[i][2]<<"\n";
		}
 		output.close();

		delete [] EG;
		delete [] backup;
		delete [] MST;
		delete [] edge;

	}

}
