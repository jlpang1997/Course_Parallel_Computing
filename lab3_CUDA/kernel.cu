
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include<iostream>
#include <stdio.h>
#define CUDACC
using namespace std;
// 矩阵类型，行优先，M(row, col) = *(M.elements + row * M.width + col)
struct Matrix
{
	int width;
	int height;
	int *elements;
};
// 获取矩阵A的(row, col)元素
__device__ int getElement(Matrix *A, int row, int col)
{
	return A->elements[row * A->width + col];
}

// 为矩阵A的(row, col)元素赋值
__device__ void setElement(Matrix *A, int row, int col, float value)
{
	A->elements[row * A->width + col] = value;
}
__device__ void print(int out)
{
	//printf("%d");
}

// 矩阵相加kernel，2-D，每个线程计算一个元素
__global__ void matAddKernel(Matrix *A, Matrix *B, Matrix *C)
{
	int Cvalue = 0;
	int row = threadIdx.y + blockIdx.y * blockDim.y;
	int col = threadIdx.x + blockIdx.x * blockDim.x;
	Cvalue = getElement(A, row, col) + getElement(B, row, col);
	setElement(C, row, col, Cvalue);
}
// 矩阵相乘kernel，2-D，每个线程计算一个元素
__global__ void matMulKernel(Matrix *A, Matrix *B, Matrix *C)
{
	int Cvalue = 0;
	int row = threadIdx.y + blockIdx.y * blockDim.y;//获取该线程所处理的矩阵行号
	int col = threadIdx.x + blockIdx.x * blockDim.x;//获取该线程所处理的矩阵列号

	for (int i = 0; i < A->width; ++i)//普通的矩阵乘法
	{
		int a, b, c;
		a = A->elements[row*A->width + i];
		b = B->elements[i*B->width + col];
		c = a * b;
		Cvalue +=c;
	}
	setElement(C, row, col, Cvalue);//把结果写回到C矩阵
}
//using namespace std;
int main()
{
	int width = 1 << 2;
	int height = 1 << 2;
	Matrix *A, *B, *C, *D;
	

	int height_A = 16;
	int width_A = 4;
	int height_B = width_A;
	int width_B = 4;

	int width_result = width_B;
	int height_result = height_A;

	// 申请托管内存  不用ppt上的拷来拷去的做法 这相当于是共享内存了吧
	cudaMallocManaged((void**)&A, sizeof(Matrix));
	cudaMallocManaged((void**)&B, sizeof(Matrix));
	cudaMallocManaged((void**)&C, sizeof(Matrix));
	cudaMallocManaged((void**)&D, sizeof(Matrix));

	cudaMallocManaged((void**)&A->elements, height_A*width_A*sizeof(int));
	cudaMallocManaged((void**)&B->elements, height_B*width_B * sizeof(int));
	cudaMallocManaged((void**)&C->elements, height_result*width_result * sizeof(int));
	cudaMallocManaged((void**)&D->elements, height_result*width_result * sizeof(int));

	A->height = height_A;
	A->width = width_A;
	B->height = height_B;
	B->width = width_B;
	C->height = height_result;
	C->width = width_result;
	D->height = height_result;
	D->width = width_result;

	for (int i = 0; i < height_A; i++)
	{
		for (int j = 0; j < width_A; j++)
		{
			A->elements[i*width_A + j]=rand()%10;
		}
	}
	for (int i = 0; i < height_B; i++)
	{
		for (int j = 0; j < width_B; j++)
		{
			B->elements[i*width_B + j]=rand()%10;
		}
	}
	dim3 blockSize(2,2);
	dim3 gridSize(width_B/blockSize.x,height_A/blockSize.y);
	// 执行kernel
	matMulKernel <<< gridSize, blockSize >>> (A, B, C);
	 //同步device 保证结果能正确访问

	cudaDeviceSynchronize();
	// 检查执行结果


	std::cout << "A:" << std::endl;
	for (int i = 0; i < height_A; ++i)
	{
		for (int j = 0; j < width_A; j++)
		{
			printf("%3d ", A->elements[i*width_A + j]);
		}
		std::cout << std::endl;
	}
	std::cout << "B:" << std::endl;
	for (int i = 0; i < height_B; ++i)
	{
		for (int j = 0; j < width_B; j++)
		{
			printf("%3d ", B->elements[i*width_B + j]);
		}
		std::cout << std::endl;
	}
	std::cout << "\n\n\nA*B:" << std::endl;
	for (int i = 0; i < height_result; i++)
	{
		for (int j = 0; j < width_result; j++)
		{
			printf("%10d ", C->elements[i*width_result + j]);
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	return 0;
}

/*
	1,学会用Nsight调试
	2，学会一些基本的cuda函数
	3，完成实验2



*/
