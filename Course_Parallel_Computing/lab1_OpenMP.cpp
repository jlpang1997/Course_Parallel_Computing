#include"lab1_OpenMP.h"
static long num_steps = 100000;
double step;
using namespace std;
void PiComputing_Parallel_Region()
{
	//int i;

	double x, pi, sum[THREAD_NUM_MAX];
	step = 1.0 / (double)num_steps;
	omp_set_num_threads(THREAD_NUM_MAX);
	#pragma omp parallel //这里面定义的全都是线程私有变量
	{
		double x;
		int id;
		id = omp_get_thread_num();
		int i;
		for (i = id, sum[id] = 0.0; i < num_steps; i += THREAD_NUM_MAX)
		{
			x = (i + 0.5)*step;
			sum[id] += 4.0 / (1.0 + x * x);
		}
	}

	int i;
	for (i = 0, pi = 0.0; i < THREAD_NUM_MAX; i++)
	{
		pi += sum[i] * step;

	}
	printf("PiComputing_Parallel_Region:%lf\n", pi);

}
void PiComputing_Parallel_ShareTask()
{
	int i;

	double x, pi, sum[THREAD_NUM_MAX];
	step = 1.0 / (double)num_steps;
	omp_set_num_threads(THREAD_NUM_MAX);
	#pragma omp parallel 
	{
		double x;
		int id;
		id = omp_get_thread_num();
		sum[id] = 0;

		#pragma omp for
		for (i = 0; i < num_steps; i++)
		{
			x = (i + 0.5)*step;
			sum[id] += 4.0 / (1.0 + x * x);
		}


	}
	for (i = 0, pi = 0.0; i < THREAD_NUM_MAX; i++)
	{
		pi += sum[i] * step;

	}
	printf("PiComputing_Parallel_ShareTask:%lf\n", pi);

}

void PiComputing_Parallel_Private_Critial()
{
	int i;
	double x=0.0, pi=0.0, sum=0.0;
	step = 1.0 / (double)num_steps;
	omp_set_num_threads(THREAD_NUM_MAX);
#pragma omp parallel private(x,sum,i)//i也应该作为私有变量，要不然循环迭代的时候，两个就混在一起了
	{
		
		int id;
		id = omp_get_thread_num();

		for (i = id,sum=0.0; i < num_steps; i+=THREAD_NUM_MAX)
		{
			x = (i + 0.5)*step;
			sum+= 4.0 / (1.0 + x * x);
		}
		#pragma omp critical//以下部分顺序执行
		{
			pi += sum * step;
		}
		


	}

	printf("PiComputing_Parallel_Private_Critial:%lf\n", pi);

}

void PiComputing_Parallel_Reduction()
{
	int i;
	double x = 0.0, pi = 0.0, sum = 0.0;
	step = 1.0 / (double)num_steps;
	omp_set_num_threads(THREAD_NUM_MAX);
#pragma omp parallel for reduction(+:sum) private(x)//i也应该作为私有变量，要不然循环迭代的时候，两个就混在一起了
		for (i = 1; i < num_steps; i++)
		{
			x = (i - 0.5)*step;
			sum += 4.0 / (1.0 + x * x);
		}
	pi = sum * step;
	printf("PiComputing_Parallel_Reduction:%lf\n", pi);

}
void part_sort(int a[], int begin, int end)//begin end 分别是最小最大下标 用最简单的选择排序
{
	for (int i = begin; i <= end; i++)
	{
		int k = i;
		for (int j = i + 1; j <= end; j++)
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
			a[k] = tmp;
		}

	}
}
void PSRS_OpenMP(int array[], int n, int thread_num)//n为数组实际容量
{
	if (thread_num*thread_num >= n||thread_num<=1||thread_num>THREAD_NUM_MAX||n>NMAX)
	{
		printf("error\n");
		exit(0);
	}
	int array_tmp[NMAX];
	int main_elem[THREAD_NUM_MAX - 1];//这里就规定了，至少要两个线程，同时，要从每个子数组中选取线程数个元素，所以线程数的平方要小于NMAX
	int partition[THREAD_NUM_MAX][THREAD_NUM_MAX];//记录每段元素个数
	int part_len = n / thread_num;
	int sample[THREAD_NUM_MAX*THREAD_NUM_MAX];
	//int index[THREAD_NUM_MAX+1];
	//index[0] = 0;
	//index[thread_num] = n;
	//omp_set_num_threads(NUM_THREADS);
	omp_set_num_threads(thread_num);
#pragma omp parallel
	{
		//局部排序
		int id = omp_get_thread_num();
		int tmp = (id == thread_num - 1) ? n : (id + 1)*part_len;
		part_sort(array, id*part_len, tmp - 1);

//#pragma omp barrier
		//选取样本
		int step = part_len / thread_num;
		for (int j = 0; j < part_len; j += step)
		{
			sample[id*thread_num + j / step] = array[id*part_len + j];
		}
#pragma omp barrier
		//样本排序,只给0线程干
		if (id == 0)
			part_sort(sample, 0, thread_num*thread_num - 1);
#pragma omp barrier

		//选取主元
		if (id < thread_num - 1)
			main_elem[id] = sample[thread_num*(id + 1)];//最后一个线程不用工作
#pragma omp barrier
		//主元划分
		int start = id * part_len;
		for (int j = 0; j < thread_num; j++)
		{
			int tmp = (id == thread_num - 1) ? n : (id + 1)*part_len;
			if (j == thread_num - 1)
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
		for (int j = 0; j < thread_num; j++)
		{
			int start_index_tmp = 0;//找到临时数组的起点
			for (int k = 0; k < j; k++)
			{
				for (int b = 0; b < thread_num; b++)
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
		//			for (int j = 1; j < thread_num; j++)
		//			{
		//				int tmp = 0;
		//				for (int idi = 0; idi < thread_num; idi++)
		//				{
		//					tmp += partition[idi][j - 1];
		//				}
		//				index[j] = tmp + index[j - 1];
		//			}
		//		}
#pragma omp barrier
		//归并排序
		int starti = 0, endi = 0;
		for (int j = 0; j < id; j++)
		{
			for (int idi = 0; idi < thread_num; idi++)
			{
				starti += partition[idi][j];
			}
		}
		endi = starti;
		for (int idi = 0; idi < thread_num; idi++)
		{
			endi += partition[idi][id];
		}
		part_sort(array_tmp, starti, endi - 1);
		for (int i = starti; i < endi; i++)
		{
			array[i] = array_tmp[i];
		}

	}
	//for (int i = 0; i < n; i++)
	//{
	//	//cout << array_tmp[i] << " ";
	//	
	//}
	//cout << endl;

}
void main_lab1()
{
	//四种pi的计算方法
	//PiComputing_Parallel_Region();
	//PiComputing_Parallel_ShareTask();
	//PiComputing_Parallel_Private_Critial();
	//PiComputing_Parallel_Reduction();
	int array[NMAX];
	for (int i = 0; i < 20; i++)
	{
		array[i] = rand() % 1000000;
	}

	PSRS_OpenMP(array, 20, 3);
	for (int i = 0; i < 20; i++)
		cout << array[i] << endl;

}