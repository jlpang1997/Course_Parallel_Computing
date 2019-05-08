//#include<stdio.h>
//#include<string.h>
//#include<mpi.h>
////#pragma comment(lib,"msmpi.lib")
//
//int main(int argc, char* argv[])
//{
//	int proceNum;//进程数量
//
//
//	int thisId;//当前进程id
//
//	MPI_Status  status;//状态信息
//
//	char message[1024];//发送信息内容
//
//	MPI_Init(&argc, &argv);//初始化mpi进程
//
//
//	MPI_Comm_rank(MPI_COMM_WORLD, &thisId);//获取当前进程号
//
//	if (thisId == 0)
//	{
//		MPI_Comm_size(MPI_COMM_WORLD, &proceNum);//获取进程数量
//
//		printf("当前共有%d个进程\n", proceNum);
//	}
//
//	if (thisId != 0) //当前进程不是主进程
//	{
//		printf("我是进程%d,正在发送......", thisId);
//
//		sprintf(message, "你好主进程！！");
//		// int MPI_Send(void *message, int count, MPI_Datatype datatype, int rank, int tag, MPI_Comm comm);发送信息到0进程（主进程）
//		MPI_Send(message, strlen(message) + 1, MPI_CHAR, 0, 9527, MPI_COMM_WORLD); //MPI_COMM_WORLD:包含程序中所有MPI进程
//
//	}
//	else {//thisId != 0即不是主进程，那么这里就是主进程
//		for (int i = 1; i < proceNum; i++)//其实这是阻塞式的发送和接受，所以可以靠循环来读取message所指缓冲区的信息
//		{
//			// int MPI_Recv(void *message, int count, MPI_Datatype datatype, int rank, int tag, MPI_Comm comm, MPI_Status *status);//接受信息
//			MPI_Recv(message, 1024, MPI_CHAR, i, 9527, MPI_COMM_WORLD, &status);
//
//
//			printf("你好，进程%d\n", status.MPI_SOURCE);
//
//
//		}
//	}
//
//
//
//	MPI_Finalize();
//
//
//
//	return 0;
//}