#pragma once
#include"stdio.h"
#include<omp.h>
#include<iostream>
#define NUM_THREADS 5
void PiComputing_Parallel_Region();
void PiComputing_Parallel_ShareTask();
void PiComputing_Parallel_Private_Critial();
void PiComputing_Parallel_Reduction();