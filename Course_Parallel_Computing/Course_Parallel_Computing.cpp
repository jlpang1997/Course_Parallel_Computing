#include <iostream>
#include"lab1_OpenMP.h"
#include<omp.h>
#define N 30
#include<time.h>
void part_sort(int a[], int begin, int end)//begin end 分别是最小最大下标 用最简单的选择排序
{
	for (int i = begin; i <= end; i++)
	{
		int k = i;
		for (int j = i+1; j <= end; j++)
		{
			if (a[j] < a[k])
			{
				k = j;
			}
		}
		if (k != i)
		{
			int tmp = a[i];
			a[i] = a[k];
			a[k] =tmp;
		}

	}
}
using namespace std;
int main()
{

	//PiComputing_Parallel_Region();
	//PiComputing_Parallel_ShareTask();

	//PiComputing_Parallel_Private_Critial();
	//PiComputing_Parallel_Reduction();
	int array[N] = { 15,46,48,93,39,6,72,91,14,
				36,69,40,89,61,97,12,21,54,
				53,97,84,58,32,27,33,72,20 },
		array_tmp[N];
	for(int i=0;i<N;i++)
	{
		array[i] = rand() % 1000000;
		//printf("%d ", array[i]);
	}
	cout << endl;


	int main_elem[NUM_THREADS - 1];//这里就规定了，至少要两个线程，同时，要从每个子数组中选取线程数个元素，所以线程数的平方要小于N
	int partition[NUM_THREADS][NUM_THREADS];//记录每段元素个数
	int part_len = N / NUM_THREADS;
	//int index[NUM_THREADS+1];

	int a[NUM_THREADS*NUM_THREADS];
	//index[0] = 0;
	//index[NUM_THREADS] = N;

	omp_set_num_threads(NUM_THREADS);
#pragma omp parallel
	{
		//局部排序
		int id = omp_get_thread_num();
		int tmp = (id == NUM_THREADS - 1) ? N : (id + 1)*part_len;
		part_sort(array, id*part_len, tmp - 1);

#pragma omp barrier
		//选取样本
		int step = part_len / NUM_THREADS;
		for (int j = 0; j < part_len; j+=step)
		{
			a[id*NUM_THREADS+ j/step] = array[id*part_len+j];
		}
#pragma omp barrier
		//样本排序,只给0线程干
		if(id==0)
			part_sort(a, 0, NUM_THREADS*NUM_THREADS-1);
#pragma omp barrier

		//选取主元
		if(id<NUM_THREADS-1)
			main_elem[id] = a[NUM_THREADS*(id+1)];//最后一个线程不用工作
#pragma omp barrier
		//主元划分
		int start = id * part_len;
		for (int j = 0; j < NUM_THREADS; j++)
		{
			int tmp = (id == NUM_THREADS - 1) ? N : (id + 1)*part_len;
			if (j == NUM_THREADS - 1)
			{
				partition[id][j] = tmp - start;
			}
			else
			{
				int i;
				for (i = start; i < tmp&&array[i] <= main_elem[j]; i++);
				partition[id][j] = i - start;
				start = i;
			}
		}
#pragma omp barrier
		//全局交换
		for (int j = 0; j < NUM_THREADS; j++)
		{
			int start_index_tmp = 0;//找到临时数组的起点
			for (int k = 0; k < j; k++)
			{
				for(int b=0;b<NUM_THREADS;b++)
				{
					start_index_tmp += partition[b][k];
				}
			}
			for (int b = 0; b < id; b++)
			{
				start_index_tmp += partition[b][j];
			}
			int start_index_init = 0;//找到原来数组地起点
			for (int m = 0; m < j; m++)
			{
				start_index_init += partition[id][m];
			}
			int len = partition[id][j];//获得要交换的子子数组的大小
			for (int i = start_index_tmp, k = start_index_init, n = 0; n < len; i++, k++, n++)
			{
				array_tmp[i] = array[id*part_len + k];
			}
		}

//#pragma omp barrier//专门用一个处理器来找到新的分段
//
//		if (id == 0)
//		{
//			index[0] = 0;
//			for (int j = 1; j < NUM_THREADS; j++)
//			{
//				int tmp = 0;
//				for (int idi = 0; idi < NUM_THREADS; idi++)
//				{
//					tmp += partition[idi][j - 1];
//				}
//				index[j] = tmp + index[j - 1];
//			}
//		}
#pragma omp barrier
		//归并排序
		int starti=0, endi=0;
		for (int j=0; j < id; j++)
		{
			for (int idi = 0; idi < NUM_THREADS; idi++)
			{
				starti += partition[idi][j];
			}
		}
		endi = starti;
		for (int idi = 0; idi < NUM_THREADS; idi++)
		{
			endi += partition[idi][id];
		}
			part_sort(array_tmp,starti,endi- 1);
	}
	for (int i = 0; i < N; i++)
	{
		cout << array_tmp[i] << " ";
	}
	cout << endl;
	return 0;
}
