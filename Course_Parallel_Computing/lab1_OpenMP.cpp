#include"lab1_OpenMP.h"

static long num_steps = 100000;
double step;
using namespace std;
void PiComputing_Parallel_Region()
{
	//int i;

	double x, pi, sum[NUM_THREADS];
	step = 1.0 / (double)num_steps;
	omp_set_num_threads(NUM_THREADS);
	#pragma omp parallel //这里面定义的全都是线程私有变量
	{
		double x;
		int id;
		id = omp_get_thread_num();
		int i;
		for (i = id, sum[id] = 0.0; i < num_steps; i += NUM_THREADS)
		{
			x = (i + 0.5)*step;
			sum[id] += 4.0 / (1.0 + x * x);
		}
	}

	int i;
	for (i = 0, pi = 0.0; i < NUM_THREADS; i++)
	{
		pi += sum[i] * step;

	}
	printf("PiComputing_Parallel_Region:%lf\n", pi);

}
void PiComputing_Parallel_ShareTask()
{
	int i;

	double x, pi, sum[NUM_THREADS];
	step = 1.0 / (double)num_steps;
	omp_set_num_threads(NUM_THREADS);
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
	for (i = 0, pi = 0.0; i < NUM_THREADS; i++)
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
	omp_set_num_threads(NUM_THREADS);
#pragma omp parallel private(x,sum,i)//i也应该作为私有变量，要不然循环迭代的时候，两个就混在一起了
	{
		
		int id;
		id = omp_get_thread_num();

		for (i = id,sum=0.0; i < num_steps; i+=NUM_THREADS)
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
	omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for reduction(+:sum) private(x)//i也应该作为私有变量，要不然循环迭代的时候，两个就混在一起了
		for (i = 1; i < num_steps; i++)
		{
			x = (i - 0.5)*step;
			sum += 4.0 / (1.0 + x * x);
		}
	pi = sum * step;
	printf("PiComputing_Parallel_Reduction:%lf\n", pi);

}