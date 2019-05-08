#include<stdio.h>
#include<string.h>
#include<mpi.h>
double pi=0.0;
void PiComputing_MPI(int argc, char *argv[])
{
	int group_size;
	int my_rank;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &group_size);

	long n = 2000;
	//MPI_Bcast(&n, 1, MPI_LONG, 0, MPI_COMM_WORLD);

	double h = 1.0 / (double)n;
	double sum = 0.0;
	for (int i = my_rank+1; i <= n; i += group_size)
	{
		double x = h * (i - 0.5);
		sum += 4.0 / (1.0 + x * x);
	}
	double mypi = h * sum;
	MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (my_rank == 0)
	{
		printf("pi is appraximately:%.16lf\n", pi);
	}
	MPI_Finalize();


}
void main_lab2(int argc,char *argv[])
{
	PiComputing_MPI( argc, argv);
	//int proceNum;//��������


	//int thisId;//��ǰ����id

	//MPI_Status  status;//״̬��Ϣ

	//char message[1024];//������Ϣ����

	//MPI_Init(&argc, &argv);//��ʼ��mpi����


	//MPI_Comm_rank(MPI_COMM_WORLD, &thisId);//��ȡ��ǰ���̺�

	//if (thisId == 0)
	//{
	//	MPI_Comm_size(MPI_COMM_WORLD, &proceNum);//��ȡ��������

	//	printf("��ǰ����%d������\n", proceNum);
	//}

	//if (thisId != 0) //��ǰ���̲���������
	//{
	//	printf("���ǽ���%d,���ڷ���......", thisId);

	//	sprintf(message, "��������̣���");
	//	// int MPI_Send(void *message, int count, MPI_Datatype datatype, int rank, int tag, MPI_Comm comm);������Ϣ��0���̣������̣�
	//	MPI_Send(message, strlen(message) + 1, MPI_CHAR, 0, 9527, MPI_COMM_WORLD); //MPI_COMM_WORLD:��������������MPI����

	//}
	//else {//thisId != 0�����������̣���ô�������������
	//	for (int i = 1; i < proceNum; i++)//��ʵ��������ʽ�ķ��ͺͽ��ܣ����Կ��Կ�ѭ������ȡmessage��ָ����������Ϣ
	//	{
	//		// int MPI_Recv(void *message, int count, MPI_Datatype datatype, int rank, int tag, MPI_Comm comm, MPI_Status *status);//������Ϣ
	//		MPI_Recv(message, 1024, MPI_CHAR, i, 9527, MPI_COMM_WORLD, &status);


	//		printf("��ã�����%d\n", status.MPI_SOURCE);


	//	}
	//}



	//MPI_Finalize();



	//return 0;
}