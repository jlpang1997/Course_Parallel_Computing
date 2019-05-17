
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include<iostream>
#include <stdio.h>
#define CUDACC


void MatrixMulti_ShareMem(int *A, int *B, int *C, int n);
__global__ void MatrixMulti_ShareMem_device(int *A, int *B, int *C, int n);
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

	for (int i = 0; i < A->width; ++i)//��ͨ�ľ���˷�
	{
		int a, b, c;
		a = A->elements[row*A->width + i];
		b = B->elements[i*B->width + col];
		c = a * b;
		Cvalue +=c;
	}
	setElement(C, row, col, Cvalue);//�ѽ��д�ص�C����
}
//using namespace std;
int main()
{
	//int width = 1 << 2;
	//int height = 1 << 2;
	//Matrix *A, *B, *C, *D;
	

	int height_A = 4;
	int width_A = 4;
	int height_B = width_A;
	int width_B = 4;

	int width_result = width_B;
	int height_result = height_A;

	//// �����й��ڴ�  ����ppt�ϵĿ�����ȥ������ ���൱���ǹ����ڴ��˰�
	//cudaMallocManaged((void**)&A, sizeof(Matrix));
	//cudaMallocManaged((void**)&B, sizeof(Matrix));
	//cudaMallocManaged((void**)&C, sizeof(Matrix));
	//cudaMallocManaged((void**)&D, sizeof(Matrix));

	//cudaMallocManaged((void**)&A->elements, height_A*width_A*sizeof(int));
	//cudaMallocManaged((void**)&B->elements, height_B*width_B * sizeof(int));
	//cudaMallocManaged((void**)&C->elements, height_result*width_result * sizeof(int));
	//cudaMallocManaged((void**)&D->elements, height_result*width_result * sizeof(int));

	//A->height = height_A;
	//A->width = width_A;
	//B->height = height_B;
	//B->width = width_B;
	//C->height = height_result;
	//C->width = width_result;
	//D->height = height_result;
	//D->width = width_result;

	//for (int i = 0; i < height_A; i++)
	//{
	//	for (int j = 0; j < width_A; j++)
	//	{
	//		A->elements[i*width_A + j]=rand()%10;
	//	}
	//}
	//for (int i = 0; i < height_B; i++)
	//{
	//	for (int j = 0; j < width_B; j++)
	//	{
	//		B->elements[i*width_B + j]=rand()%10;
	//	}
	//}
	//dim3 blockSize(1,2);
	//dim3 gridSize(width_B/blockSize.x,height_A/blockSize.y);
	//// ִ��kernel
	//matMulKernel <<< gridSize, blockSize >>> (A, B, C);
	// //ͬ��device ��֤�������ȷ����

	//cudaDeviceSynchronize();
	//// ���ִ�н��


	//std::cout << "A:" << std::endl;
	//for (int i = 0; i < height_A; ++i)
	//{
	//	for (int j = 0; j < width_A; j++)
	//	{
	//		printf("%3d ", A->elements[i*width_A + j]);
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << "B:" << std::endl;
	//for (int i = 0; i < height_B; ++i)
	//{
	//	for (int j = 0; j < width_B; j++)
	//	{
	//		printf("%3d ", B->elements[i*width_B + j]);
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << "\n\n\nA*B:" << std::endl;
	//for (int i = 0; i < height_result; i++)
	//{
	//	for (int j = 0; j < width_result; j++)
	//	{
	//		printf("%10d ", C->elements[i*width_result + j]);
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;


	int *A = (int *)malloc(sizeof(int)*height_A*width_A);
	int *B = (int*)malloc(sizeof(int)*height_B*width_B);
	int*C = (int *)malloc(sizeof(int)*height_result*width_result);
	for (int i = 0; i < height_A; i++)
	{
		for (int j = 0; j < width_A; j++)
		{
			A[i*width_A + j]=rand()%10;
			printf("%3d ", A[i*width_A + j]);
		}
		printf("\n");
	}
	cout << endl;

	for (int i = 0; i < height_B; i++)
	{
		for (int j = 0; j < width_B; j++)
		{
			B[i*width_B + j]=rand()%10;
			printf("%3d ", B[i*width_B + j]);
		}
		printf("\n");
	}
	cout << endl << endl;
	MatrixMulti_ShareMem(A, B, C, height_A);
	for (int i = 0; i < height_result; i++)
	{
		for (int j = 0; j < width_result; j++)
		{
			printf("%3d ", C[i*width_result + j]);
		}
		printf("\n");
	}
	free(A);
	free(B);
	free(C);
	return 0;
}
void MatrixMulti_ShareMem(int *A, int *B, int *C, int n)
{
	int *cuda_A, *cuda_B, *cuda_C;
	int size = sizeof(int)*n*n;
	cudaMalloc(&cuda_A, size);
	cudaMalloc(&cuda_B, size);
	cudaMalloc(&cuda_C, size);//cuda��ȫ���ڴ�
	
	cudaMemcpy(cuda_A, A, size, cudaMemcpyHostToDevice);
	cudaMemcpy(cuda_B, B, size, cudaMemcpyHostToDevice);
	//����������ŵ�gpu��ȫ���ڴ���

	dim3 blocksize(2, 2);
	dim3 gridsize(n / blocksize.x, n / blocksize.y);



	MatrixMulti_ShareMem_device << <gridsize, blocksize ,blocksize.x*blocksize.y*sizeof(int)>> > (cuda_A, cuda_B, cuda_C, n);//Ҫָ����̬����Ĺ����ڴ��С

	cudaMemcpy(C, cuda_C,size, cudaMemcpyDeviceToHost);
	cudaFree(cuda_A);
	cudaFree(cuda_B);
	cudaFree(cuda_C);
}


__global__ void MatrixMulti_ShareMem_device(int *A, int *B, int *C, int n)
{
	int x_start_A = blockIdx.x*blockDim.x;
	int y_start_A = blockIdx.y*blockDim.y;
	
	//���߳�����block�ĵ�һ��Ԫ�ص�����
	int x_start_B= blockIdx.x*blockDim.x;
	int y_start_B= blockIdx.y*blockDim.y;


	int value = 0;
	for (int i = 0; i < n / blockDim.x; i++)
	{
		x_start_A += i * blockDim.x;
		y_start_B += i * blockDim.y;

		extern __shared__  int matrix_A[],matrix_B[];

		matrix_A[threadIdx.y*blockDim.x + threadIdx.x] = A[(y_start_A+threadIdx.y)*n + x_start_A+threadIdx.x];

		matrix_B[threadIdx.y*blockDim.x + threadIdx.y] = B[(y_start_B + threadIdx.y)*n + x_start_B + threadIdx.x];

		//__syncthreads();

		for (int i = 0; i < blockDim.x; i++)
		{
			value += matrix_A[threadIdx.y*blockDim.x + i] * matrix_B[i*blockDim.x + threadIdx.x];
		}

	}


}


/*
	1,ѧ����Nsight����
	2��ѧ��һЩ������cuda����
	3�����ʵ��2



*/
