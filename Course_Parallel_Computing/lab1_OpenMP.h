#pragma once
#include"stdio.h"
#include<omp.h>
#include<iostream>
#define NMAX 10000//��Ϊ������������
#define THREAD_NUM_MAX 1000
void PiComputing_Parallel_Region();
void PiComputing_Parallel_ShareTask();
void PiComputing_Parallel_Private_Critial();
void PiComputing_Parallel_Reduction();
void PSRS_OpenMP(int array[], int n,int thread_num);
void part_sort(int a[], int begin, int end);
void main_lab1();