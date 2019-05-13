
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include<iostream>
#include <stdio.h>
#define CUDACC
using namespace std;
// �������ͣ������ȣ�M(row, col) = *(M.elements + row * M.width + col)
struct Matrix
{
	int width;
	int height;
	int *elements;
};
// ��ȡ����A��(row, col)Ԫ��
__device__ int getElement(Matrix *A, int row, int col)
{
	return A->elements[row * A->width + col];
}

// Ϊ����A��(row, col)Ԫ�ظ�ֵ
__device__ void setElement(Matrix *A, int row, int col, float value)
{
	A->elements[row * A->width + col] = value;
}
__device__ void print(int out)
{
	//printf("%d");
}

// �������kernel��2-D��ÿ���̼߳���һ��Ԫ��
__global__ void matAddKernel(Matrix *A, Matrix *B, Matrix *C)
{
	int Cvalue = 0;
	int row = threadIdx.y + blockIdx.y * blockDim.y;
	int col = threadIdx.x + blockIdx.x * blockDim.x;
	Cvalue = getElement(A, row, col) + getElement(B, row, col);
	setElement(C, row, col, Cvalue);
}
// �������kernel��2-D��ÿ���̼߳���һ��Ԫ��
__global__ void matMulKernel(Matrix *A, Matrix *B, Matrix *C)
{

	int Cvalue = 0;
	int row = threadIdx.y + blockIdx.y * blockDim.y;//��ȡ���߳�������ľ����к�
	int col = threadIdx.x + blockIdx.x * blockDim.x;//��ȡ���߳�������ľ����к�

	//col = 0;
	//printf("(%d , %d)\n", row, col);
	for (int i = 0; i < A->height; i++)
	{
		for (int j = 0; j < A->width; j++)
		{
			printf("%d ", A->elements[i*A->width + j]);
		}
		printf("\n");
	}
	for (int i = 0; i < B->height; i++)
	{
		for (int j = 0; j < B->width; j++)
		{
			printf("%d ", B->elements[i*A->width + j]);
		}
		printf("\n");
	}
	for (int i = 0; i < A->width; ++i)//��ͨ�ľ���˷�
	{
		int a, b, c;
		//a = getElement(A, row, i);
		a = A->elements[row*A->width + i];
		b = B->elements[i*B->width + col];
		//b = getElement(B, i, col);
		c = a * b;
		Cvalue +=c;
		printf("%d * %d = %d \n", a,b,c);
	}
	//printf("\n");
	setElement(C, row, col, Cvalue);//�ѽ��д�ص�C����
}
//using namespace std;
int main()
{
	int width = 1 << 2;
	int height = 1 << 2;
	Matrix *A, *B, *C, *D;
	// �����й��ڴ�

	int width_A = 4;
	int height_A = 4;
	int height_B = width_A;
	int width_B = 1;

	int width_result = width_B;
	int height_result = height_A;


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
			A->elements[i*height_A + j]=rand()%10;
		}
	}
	for (int i = 0; i < height_B; i++)
	{
		for (int j = 0; j < width_B; j++)
		{
			B->elements[i*height_B + j]=rand()%10;
		}
	}

	//dim3 blockSize(1, 2);
	//dim3 gridSize((width + blockSize.x - 1) / blockSize.x,
	//	(height + blockSize.y - 1) / blockSize.y);
	////std::cout << (width + blockSize.x - 1) / blockSize.x << std::endl;
	//// ִ��kernel
	//matMulKernel <<< gridSize, blockSize >>> (A, B, C);
	 //ͬ��device ��֤�������ȷ����


	dim3 blockSize(1, 2);
	dim3 gridSize((width_result + blockSize.x - 1) / blockSize.x,
		(height_result + blockSize.y - 1) / blockSize.y);
	// ִ��kernel
	matMulKernel <<< gridSize, blockSize >>> (A, B, C);
	 //ͬ��device ��֤�������ȷ����

	cudaDeviceSynchronize();
	// ���ִ�н��


	std::cout << "A:" << std::endl;
	for (int i = 0; i < height_A; ++i)
	{
		for (int j = 0; j < width_A; j++)
		{
			printf("%3d ", A->elements[i*height_A + j]);
		}
		std::cout << std::endl;
	}
	std::cout << "B:" << std::endl;
	for (int i = 0; i < height_B; ++i)
	{
		for (int j = 0; j < width_B; j++)
		{
			printf("%3d ", B->elements[i*height_B + j]);
		}
		std::cout << std::endl;
	}
	std::cout << "A*B:" << std::endl;
	for (int i = 0; i < height_result; ++i)
	{
		for (int j = 0; j < width_result; j++)
		{
			printf("%3d ", C->elements[i*height_result + j]);
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	return 0;
}
