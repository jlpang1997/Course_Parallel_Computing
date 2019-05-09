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
	#pragma omp parallel //�����涨���ȫ�����߳�˽�б���
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
#pragma omp parallel private(x,sum,i)//iҲӦ����Ϊ˽�б�����Ҫ��Ȼѭ��������ʱ�������ͻ���һ����
	{
		
		int id;
		id = omp_get_thread_num();

		for (i = id,sum=0.0; i < num_steps; i+=THREAD_NUM_MAX)
		{
			x = (i + 0.5)*step;
			sum+= 4.0 / (1.0 + x * x);
		}
		#pragma omp critical//���²���˳��ִ��
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
#pragma omp parallel for reduction(+:sum) private(x)//iҲӦ����Ϊ˽�б�����Ҫ��Ȼѭ��������ʱ�������ͻ���һ����
		for (i = 1; i < num_steps; i++)
		{
			x = (i - 0.5)*step;
			sum += 4.0 / (1.0 + x * x);
		}
	pi = sum * step;
	printf("PiComputing_Parallel_Reduction:%lf\n", pi);

}
void part_sort(int a[], int begin, int end)//begin end �ֱ�����С����±� ����򵥵�ѡ������
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
void PSRS_OpenMP(int array[], int n, int thread_num)//nΪ����ʵ������
{
	if (thread_num*thread_num >= n||thread_num<=1||thread_num>THREAD_NUM_MAX||n>NMAX)
	{
		printf("error\n");
		exit(0);
	}
	int array_tmp[NMAX];
	int main_elem[THREAD_NUM_MAX - 1];//����͹涨�ˣ�����Ҫ�����̣߳�ͬʱ��Ҫ��ÿ����������ѡȡ�߳�����Ԫ�أ������߳�����ƽ��ҪС��NMAX
	int partition[THREAD_NUM_MAX][THREAD_NUM_MAX];//��¼ÿ��Ԫ�ظ���
	int part_len = n / thread_num;
	int sample[THREAD_NUM_MAX*THREAD_NUM_MAX];
	//int index[THREAD_NUM_MAX+1];
	//index[0] = 0;
	//index[thread_num] = n;
	//omp_set_num_threads(NUM_THREADS);
	omp_set_num_threads(thread_num);
#pragma omp parallel
	{
		//�ֲ�����
		int id = omp_get_thread_num();
		int tmp = (id == thread_num - 1) ? n : (id + 1)*part_len;
		part_sort(array, id*part_len, tmp - 1);

//#pragma omp barrier
		//ѡȡ����
		int step = part_len / thread_num;
		for (int j = 0; j < part_len; j += step)
		{
			sample[id*thread_num + j / step] = array[id*part_len + j];
		}
#pragma omp barrier
		//��������,ֻ��0�̸߳�
		if (id == 0)
			part_sort(sample, 0, thread_num*thread_num - 1);
#pragma omp barrier

		//ѡȡ��Ԫ
		if (id < thread_num - 1)
			main_elem[id] = sample[thread_num*(id + 1)];//���һ���̲߳��ù���
#pragma omp barrier
		//��Ԫ����
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
		//ȫ�ֽ���
		for (int j = 0; j < thread_num; j++)
		{
			int start_index_tmp = 0;//�ҵ���ʱ��������
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
			int start_index_init = 0;//�ҵ�ԭ����������
			for (int m = 0; m < j; m++)
			{
				start_index_init += partition[id][m];
			}
			int len = partition[id][j];//���Ҫ��������������Ĵ�С
			for (int i = start_index_tmp, k = start_index_init, n = 0; n < len; i++, k++, n++)
			{
				array_tmp[i] = array[id*part_len + k];
			}
		}

		//#pragma omp barrier//ר����һ�����������ҵ��µķֶ�
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
		//�鲢����
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
	//����pi�ļ��㷽��
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