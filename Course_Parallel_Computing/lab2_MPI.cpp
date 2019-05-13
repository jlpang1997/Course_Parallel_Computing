#include<stdio.h>
#include<string.h>
#include<iostream>
#include<mpi.h>
#include<malloc.h>
#include"lab1_OpenMP.h"
#include"lab2_MPI.h"
#include<iomanip>
struct Meg
{
	int starti;
	int len;
};
using namespace std;
double pi=0.0;
void PiComputing_MPI(int argc, char *argv[])
{
	int id,thread_num;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &thread_num);
	long n = 10000;
	double h = 1.0 / (double)n;
	double sum = 0.0;
	for (int i = id+1; i <= n; i += thread_num)//���̷ֱ߳��󲿷ֺ�
	{
		double x = h * (i - 0.5);
		sum += 4.0 / (1.0 + x * x);
	}
	double mypi = h * sum;
	MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);//��Լ���
	if (id == 0)
	{
		printf("pi:%.20lf\n", pi);
	}
	MPI_Finalize();
}

void PSRS_MPI(int argc, char *argv[], int array[], int n)
{
	int thread_num;
	int id;
	int part_len;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &thread_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	double starttime, endtime;
	starttime = MPI_Wtime();

	//ȫ���ö�̬���飬��󲻻��գ��ɳ����Լ�����
	int *array_tmp = (int *)malloc(sizeof(int)*n);
	int *main_elem = (int *)malloc(sizeof(int)*(thread_num - 1));
	int *sample = (int *)malloc(sizeof(int)*thread_num*thread_num);
	part_len = n / thread_num;

	if (thread_num*thread_num >= n || thread_num <= 1 || thread_num > THREAD_NUM_MAX || n > NMAX)
	{
		if(id==0)
			printf("error\n");
		MPI_Finalize();
		exit(0);
	}
	if (id == 0)
	{
		cout << "ԭʼ���飺" << endl;
		for (int i = 0; i < n; i++)
		{

			//printf("%d %s", array[i], (i >0&&(i+1)%(part_len )==0&&i<part_len*(thread_num-1) )? "\n" : " ");
			cout << array[i] << " ";
		}
		cout << endl;
	}
	//�ֲ�����,Ȼ��ÿһ���̶߳��Ѹ�������õ��ǲ��ֹ㲥�����̣߳����յĽ����array_tmp�������ľֲ��������顣
	int tmp = (id == thread_num - 1) ? n : (id + 1)*part_len;
	Quick_Sort(array, id*part_len, tmp - 1);
	MPI_Barrier(MPI_COMM_WORLD);
	//�����������Ԫ�ظ������߳�������������������ᵼ�»��ֲ����ȣ������޷�ֱ����gather������Ҫ���Ԫ�ظ�����ƫ��������gatherv
	int *revccount = (int *)malloc(sizeof(int)*thread_num);
	int *displs = (int *)malloc(sizeof(int)*thread_num);
	for (int i = 0; i < thread_num; i++)
	{
		revccount[i]= (i == thread_num - 1) ? (n - i * part_len) : part_len;//�̸߳���Ԫ�ظ���
	};
	displs[0] = 0;
	for (int i = 1; i < thread_num; i++)
	{
		displs[i] = displs[i - 1] + revccount[i - 1];//ƫ��������
	}
	//�����ظ����߳�
	MPI_Gatherv(&array[id*part_len], revccount[id], MPI_INT, array_tmp, revccount, displs,MPI_INT,0, MPI_COMM_WORLD);
	free(revccount);
	free(displs);

	//���̲߳�������������ѡȡ��Ԫ��������Ԫ��array�㲥�������߳�
	if (id == 0)
	{
		for (int i = 0; i < n; i++)
			array[i] = array_tmp[i];
		int step = part_len / thread_num;
		for (int idi = 0; idi < thread_num; idi++)
		{
			for (int j = 0; j < part_len; j += step)
			{
				sample[idi*thread_num + j / step] = array_tmp[idi*part_len + j];
			}
			main_elem[idi] = sample[thread_num*(idi + 1)];
		}
		Quick_Sort(sample, 0, thread_num*thread_num - 1);
		for (int idi = 0; idi < thread_num-1; idi++)
		{
			main_elem[idi] = sample[thread_num*(idi + 1)];
		}
	}
	MPI_Bcast(&array[0],NMAX, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(main_elem, thread_num - 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	//��Ԫ����,partition�����ÿһ�������黮�ֳ�����������Ĵ�С
	int* partition = (int *)malloc(sizeof(int)*thread_num);
	int start = id * part_len;
	for (int j = 0; j < thread_num; j++)
	{
		int tmp = (id == thread_num - 1) ? n : (id + 1)*part_len;
		if (j == thread_num - 1)
		{
			partition[j] = tmp - start;
		}
		else
		{
			int i;
			for (i = start; i < tmp&&array[i] <= main_elem[j]; i++);
			partition[j] = i - start;
			start = i;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	//����ÿ������gather��ÿ���̵߳�partition_gather��
	int* partition_gather=(int *)malloc(sizeof(int)*thread_num*thread_num);
	MPI_Allgather(partition, thread_num, MPI_INT, partition_gather, thread_num, MPI_INT, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	//ȫ�ֽ���
	for (int j = 0; j < thread_num; j++)
	{
		int start_index_tmp = 0;//�ҵ���ʱ��������
		for (int k = 0; k < j; k++)
			for (int b = 0; b < thread_num; b++)
				start_index_tmp += partition_gather[b*thread_num+k];
		for (int b = 0; b < id; b++)
			start_index_tmp += partition_gather[b*thread_num + j];
		int start_index_init = 0;//�ҵ�ԭ����������
		for (int m = 0; m < j; m++)
			start_index_init += partition_gather[id*thread_num + m];
		int len = partition_gather[id*thread_num + j];

		for (int i = start_index_tmp, k = start_index_init, n = 0; n < len; i++, k++, n++)
		{
			array_tmp[i] = array[id*part_len + k];//�Ȱ��߳��ڲ���array_tmp��װ��
		}
		//���ڶ���ÿ���߳���˵�����ǰ�ԭ������������黮�ֳ�n���֣�Ȼ����������������n���ֻ��ǲ������ģ�û���������gather��bcast�ķ�����ֻ��һ��һ�鷢�͸����̣߳�Ȼ����װ�ˡ�
		for (int idi = 1; idi < thread_num; idi++)//�����̷߳��͸����߳�
		{
			Meg megsend;
			if (id == idi)
			{
				megsend = { start_index_tmp,len };
				MPI_Send(&megsend, sizeof(megsend), MPI_INT, 0, 1, MPI_COMM_WORLD);
				MPI_Send(&array_tmp[start_index_tmp], len, MPI_INT, 0, 1, MPI_COMM_WORLD);
			}
		}
		if (id == 0)//���߳̽��������̷߳�����array_tmp����װ
		{
			for (int idi = 1; idi < thread_num; idi++)
			{
				Meg megrecv;
				MPI_Recv(&megrecv, sizeof(megrecv), MPI_INT, idi, 1, MPI_COMM_WORLD, &status);//�Ƚ���ͷ����Ϣ��Ȼ�������������
				//cout << megrecv.starti << " " << megrecv.len << endl;
				MPI_Recv(&array_tmp[megrecv.starti], megrecv.len, MPI_INT, idi, 1, MPI_COMM_WORLD, &status);
			}
		}
	}
	MPI_Bcast(array_tmp, n, MPI_INT, 0, MPI_COMM_WORLD);//������̹߳㲥�������̣߳���֤ÿ���̵߳�array_tmp����һ�¡�

	MPI_Barrier(MPI_COMM_WORLD);
	//�鲢����
	int starti = 0, endi = 0;
	for (int j = 0; j < id; j++)
	{
		for (int idi = 0; idi < thread_num; idi++)
		{
			/*starti += partition[idi][j];*/
			starti +=  partition_gather[idi*thread_num + j];
		}
	}
	endi = starti;
	for (int idi = 0; idi < thread_num; idi++)
	{
		/*endi += partition[idi][id];*/
		endi += partition_gather[idi*thread_num + id];
	}
	Quick_Sort(array_tmp, starti, endi - 1);//���ڲ�����

	//���洦��ͬ�ϣ�����ͨ��gatherv��ÿ���źõķŵ����߳��е�array
	int subcount = endi - starti;
	int* count = (int *)malloc(sizeof(int)*thread_num);
	int *disp = (int *)malloc(sizeof(int)*thread_num);
	MPI_Allgather(&subcount, 1, MPI_INT, count, 1, MPI_INT,  MPI_COMM_WORLD);
	disp[0] = 0;
	for (int i = 1; i < thread_num; i++)
	{
		disp[i] = disp[i - 1] + count[i - 1];
	}
	MPI_Gatherv(&array_tmp[starti], count[id], MPI_INT, array, count, disp, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	endtime = MPI_Wtime();
	if (id == 0)
	{
		cout << "��������" << endl;
		for (int i = 0; i < n; i++)
		{
			cout << array[i] <<" ";
		}
		cout << endl;
		printf("runtime:%30.20lf\n", endtime - starttime);
		//cout << "runtime:" <<setw(16)<< endtime - starttime << endl;
	}
	
	MPI_Finalize();
	//free(array_tmp);
	//cout << "array end" << endl;
	
	//cout << "end" << endl;
	//free(sample);
	//cout << "end" << endl;

}

int partition(int a[], int starti, int endi)
{
	int flag = a[endi];
	int i = starti, j = starti;
	for (i = starti; i < endi; i++)//j��Ϊ��ǰС�ڲ��ֵ����һ������һ��
	{
		if (a[i] <= flag)
		{
			swap(a[i], a[j]);
			j++;
		}
	}
	swap(a[j], a[endi]);
	return j;
}
void Quick_Sort(int a[], int starti, int endi)
{
	if (starti >= endi)
		return;
	else
	{
		int q = partition(a, starti, endi);
		Quick_Sort(a, starti, q - 1);
		Quick_Sort(a, q + 1, endi);
	}
}
void main_lab2(int argc,char *argv[])
{
	//PiComputing_MPI( argc, argv);
	int n = 20;
	int array[NMAX];
	//	= {
	//	15,46,48,93,39,6,72,91,14,
	//	36,69,40,89,61,97,12,21,54,
	//	53,97,84,58,32,27,33,72,20
	//};
	for (int i = 0; i < n; i++)
	{
		array[i] = rand() % 1000000;
		//cout << array[i] << endl;
	}
	
	PSRS_MPI(argc,argv,array, n);

	
}