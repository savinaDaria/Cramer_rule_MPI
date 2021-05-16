//паралельний алгоритм MPI
#include <stdio.h>
#include <random>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <mpi.h>
#include <iostream>

using namespace std;
random_device rd;  //Will be used to obtain a seed for the random number engine
mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
uniform_real_distribution<> distrib(0, 1);

int ProcNum; // Number of available processes
int ProcRank; // Rank of current process

void CopyData(double* Matrix, double* newMatrix, int Size)
{
	for (int i = 0; i < Size; i++)
	{
		newMatrix[i] = Matrix[i];
	}
}
void RandomDataInitialization(double* pMatrix, double* pVectorB, int Size) {
	int i, j; 
	for (i = 0; i < Size; i++) {
		pVectorB[i] = distrib(gen);
		for (j = 0; j < Size; j++)
			pMatrix[i * Size + j] = distrib(gen);
	}
}
void ProcessInitialization(double*& pMatrix, double*& pVectorB, double*& pResultX, double*& pTempMatrix, int*& Indexes, int*& pProcIndexes, double*& pResult, int& Size, int& ColNum) {
	if (ProcRank == 0) {
		do {
			printf("\nEnter size of the initial objects: ");
			scanf_s("%d", &Size);
			printf("\nChosen objects size = %d\n", Size);
			if (Size <= 0)
			{
				printf("\nSize of objects must be greater than 0!\n");
			}
		} while (Size <= 0);
		
	}
	MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	pResultX = new double[Size];
	pResult = new double[Size];
	pTempMatrix = new double[Size * Size];
	pMatrix = new double[Size * Size];
	pVectorB = new double[Size];
	int RestRows; // Number of rows, that haven’t been distributed yet
	RestRows = Size;
	for (int i = 0; i < ProcRank; i++)
		RestRows = RestRows - RestRows / (ProcNum - i);
	ColNum = RestRows / (ProcNum - ProcRank);
	pProcIndexes = new int[ColNum];
	if (ProcRank == 0) {
		Indexes = new int[Size];
		for (int i = 0; i < Size; i++)
			Indexes[i] = i;
		RandomDataInitialization(pMatrix, pVectorB, Size);
	}
	
}
void PrintMatrix(double* pMatrix, int Size) {
	int i, j; // Loop variables
	for (i = 0; i < Size; i++) {
		for (j = 0; j < Size; j++)
			printf("%7.4f ", pMatrix[i * Size + j]);
		printf("\n");
	}
}
void PrintVector(double* pVector, int Size) {
	int i;
	for (i = 0; i < Size; i++)
		printf("%7.4f ", pVector[i]);
	printf("\n");
}
double GaussDeterminantCalculation(double* pMatrix, int Size)
{
	double tmp, determinant = 1;
	int count = 1;
	for (int k = 0; k < Size * Size - Size; k += Size + 1) {
		for (int i = k + Size; i < Size * Size; i += Size) {
			tmp = -pMatrix[i] / pMatrix[k];
			for (int j = i; j < i + Size; j++) {
				for (int h = k; h < k + Size; h++) {

					if (j == h + Size * count) {
						pMatrix[j] += pMatrix[h] * tmp;
					}
				}
			}
			count++;
		}
		count = 1;
	}
	for (int i = 0; i < Size * Size; i += Size + 1) {
		determinant *= pMatrix[i];
	}
	return determinant;
}
void ReplaceMatrixCol(double* pMatrix, double* pVectorB, double* &pTempMatrix, int Size, int pCol)
{
	int i, count = 0;
	for (i = 0; i < Size * Size; i++) {
		if (i == pCol + Size * count) {
			pTempMatrix[i] = pVectorB[count];
			count++;
		}
		else {
			pTempMatrix[i] = pMatrix[i];
		}
	}
}
void MainLoop(double* pMatrix, double* pVectorB, double* pResultX, double* pTempMatrix, int* pProcIndexes, int Size, int ColNum) {
	double determ_main=0;
	if (ProcRank == 0) {
			CopyData(pMatrix, pTempMatrix, Size * Size);
			determ_main = GaussDeterminantCalculation(pTempMatrix, Size);
	}
	MPI_Bcast(&determ_main, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	double tmpd = 0;
	int i;
	if (abs(determ_main)>0.00001 ) {
		for (int j = 0; j < ColNum; j++) {
			i = pProcIndexes[j];
			ReplaceMatrixCol(pMatrix, pVectorB, pTempMatrix, Size, i);
			tmpd= GaussDeterminantCalculation(pTempMatrix, Size);
			pResultX[j] = tmpd / determ_main;
		}
	}	
	else {
		printf("\n The determinant of Matrix = 0. Generate new Matrix \n");
	}
}
void DataDistribution(double* pMatrix, double* pVectorB, int* Indexes, int* pProcIndexes, int Size, int ColNum) {
	MPI_Bcast(pVectorB, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(pMatrix, Size * Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	int* pSendNum; // the number of elements sent to the process
	int* pSendInd; // the index of the first data element sent to the process
	int RestRows = Size; // Number of rows, that haven’t been distributed yet
	pSendInd = new int[ProcNum];
	pSendNum = new int[ProcNum];
	ColNum = Size / ProcNum;
	pSendNum[0] = ColNum;
	pSendInd[0] = 0;
	for (int i = 1; i < ProcNum; i++) {
		RestRows -= ColNum;
		ColNum = RestRows / (ProcNum - i);
		pSendNum[i] = ColNum;
		pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
	}
	MPI_Scatterv(Indexes, pSendNum, pSendInd, MPI_INT, pProcIndexes, pSendNum[ProcRank], MPI_INT, 0, MPI_COMM_WORLD);
	//delete[] pSendNum;
	//delete[] pSendInd;
}
void ResultReplication(double* pResultX, double* pResult, int Size, int ColNum) {
	int i; // Loop variable
	int* pReceiveNum; // Number of elements, that current process sends
	int* pReceiveInd; 
	int RestRows = Size; // Number of rows, that haven’t been distributed yet
	pReceiveNum = new int[ProcNum];
	pReceiveInd = new int[ProcNum];
	pReceiveInd[0] = 0;
	pReceiveNum[0] = Size / ProcNum;
	for (i = 1; i < ProcNum; i++) {
		RestRows -= pReceiveNum[i - 1];
		pReceiveNum[i] = RestRows / (ProcNum - i);
		pReceiveInd[i] = pReceiveInd[i - 1] + pReceiveNum[i - 1];
	}
	MPI_Allgatherv(pResultX, pReceiveNum[ProcRank], MPI_DOUBLE, pResult, pReceiveNum, pReceiveInd, MPI_DOUBLE, MPI_COMM_WORLD);
	delete[] pReceiveNum;
	delete[] pReceiveInd;
}
void ProcessTermination(double* pMatrix, double* pVectorB, double* pTempMatrix, int*& Indexes, int*& pProcIndexes) {
	delete[] pMatrix;
	delete[] pVectorB;
	delete[] pTempMatrix;
	delete[] Indexes;
	delete[] pProcIndexes;
}
void main(int argc, char* argv[]) {
	double* pVectorB, * pResultX, * pResult; // вектор b і результуючий вектор Х
	double* pMatrix, * pTempMatrix; //проміжні матриці
	int* Indexes, * pProcIndexes;//індекси стовпчиків в кожному процесі
	int Size; // розмір матриці
	double start, finish,duration;
	int ColNum; 
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	ProcessInitialization(pMatrix, pVectorB, pResultX, pTempMatrix, Indexes, pProcIndexes, pResult, Size, ColNum);
	DataDistribution(pMatrix, pVectorB, Indexes, pProcIndexes, Size, ColNum);
	//printf("\nRank of current process = %d \n", ProcRank);
	start = MPI_Wtime();
	MainLoop(pMatrix, pVectorB, pResultX, pTempMatrix, pProcIndexes, Size, ColNum);
	ResultReplication(pResultX, pResult, Size, ColNum);
	finish = MPI_Wtime();
	duration = finish - start;
	if (ProcRank == 0) {
		printf("Cramer`s rule + Gauss determinant program\n");
		/*printf("\n Matrix: \n");
		PrintMatrix(pMatrix, Size);
		printf("\n Vector B: \n");
		PrintVector(pVectorB, Size);
		printf("\n Result Vector X: \n");
		PrintVector(pResult, Size);*/
		printf("\n Time of execution: %f\n", duration);
	}
	//ProcessTermination(pMatrix, pVectorB, pTempMatrix, Indexes, pProcIndexes);
	MPI_Finalize();
}


