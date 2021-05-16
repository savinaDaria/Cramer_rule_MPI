//#include <stdio.h>
//#include <random>
//#include <stdlib.h>
//#include <conio.h>
//#include <time.h>
//#include <iostream>
//
//using namespace std;
//random_device rd;  //Will be used to obtain a seed for the random number engine
//mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
//uniform_real_distribution<> distrib(0, 1);
//
//void CopyData(double* Matrix, double* newMatrix, int Size)
//{
//	for (int i = 0; i < Size; i++)
//	{
//		newMatrix[i] = Matrix[i];
//	}
//}
//void RandomDataInitialization(double* pMatrix, double* pVectorB, int Size) {
//	int i, j;
//	for (i = 0; i < Size; i++) {
//		pVectorB[i] = distrib(gen);
//		for (j = 0; j < Size; j++)
//			pMatrix[i * Size + j] = distrib(gen);
//	}
//}
//void ProcessInitialization(double*& pMatrix, double*& pVectorB, double*& pResultX, double*& pTempMatrix, int& Size) {
//	do {
//			printf("\nEnter size of the initial objects: ");
//			scanf_s("%d", &Size);
//			printf("\nChosen objects size = %d\n", Size);
//			if (Size <= 0)
//			{
//				printf("\nSize of objects must be greater than 0!\n");
//			}
//		} while (Size <= 0);
//	pResultX = new double[Size];
//	pTempMatrix = new double[Size * Size];
//	pMatrix = new double[Size * Size];
//	pVectorB = new double[Size];
//	RandomDataInitialization(pMatrix, pVectorB, Size);
//
//}
//void PrintMatrix(double* pMatrix, int Size) {
//	int i, j; // Loop variables
//	for (i = 0; i < Size; i++) {
//		for (j = 0; j < Size; j++)
//			printf("%7.4f ", pMatrix[i * Size + j]);
//		printf("\n");
//	}
//}
//void PrintVector(double* pVector, int Size) {
//	int i;
//	for (i = 0; i < Size; i++)
//		printf("%7.4f ", pVector[i]);
//	printf("\n");
//}
//double GaussDeterminantCalculation(double* pMatrix, int Size)
//{
//	double tmp, determinant = 1;
//	int count = 1;
//	for (int k = 0; k < Size * Size - Size; k += Size + 1) {
//		for (int i = k + Size; i < Size * Size; i += Size) {
//			tmp = -pMatrix[i] / pMatrix[k];
//			for (int j = i; j < i + Size; j++) {
//				for (int h = k; h < k + Size; h++) {
//
//					if (j == h + Size * count) {
//						pMatrix[j] += pMatrix[h] * tmp;
//					}
//				}
//			}
//			count++;
//		}
//		count = 1;
//	}
//	for (int i = 0; i < Size * Size; i += Size + 1) {
//		determinant *= pMatrix[i];
//	}
//	return determinant;
//}
//void ReplaceMatrixCol(double* pMatrix, double* pVectorB, double*& pTempMatrix, int Size, int pCol)
//{
//	int i, count = 0;
//	for (i = 0; i < Size * Size; i++) {
//		if (i == pCol + Size * count) {
//			pTempMatrix[i] = pVectorB[count];
//			count++;
//		}
//		else {
//			pTempMatrix[i] = pMatrix[i];
//		}
//	}
//}
//void MainLoop(double* pMatrix, double* pVectorB, double* pResultX, double* pTempMatrix, int Size) {
//	double determ_main,  tmpd ;
//	CopyData(pMatrix, pTempMatrix, Size * Size);
//	determ_main = GaussDeterminantCalculation(pTempMatrix, Size);
//	if (abs(determ_main) > 0.00001) {
//		for (int j = 0; j < Size; j++) {
//			ReplaceMatrixCol(pMatrix, pVectorB, pTempMatrix, Size, j);
//			tmpd = GaussDeterminantCalculation(pTempMatrix, Size);
//			pResultX[j] = tmpd / determ_main;
//		}
//	}
//	else {
//		printf("\n The determinant of Matrix = 0. Generate new Matrix \n");
//	}
//}
//void ProcessTermination(double* pMatrix, double* pVectorB, double* pTempMatrix) {
//	delete[] pMatrix;
//	delete[] pVectorB;
//	delete[] pTempMatrix;
//}
//void main(int argc, char* argv[]) {
//	double* pVectorB, * pResultX; // вектор b і результуючий вектор Х
//	double* pMatrix, * pTempMatrix; //проміжні матриці
//	int Size; // розмір матриці
//	double start, finish, duration;
//	ProcessInitialization(pMatrix, pVectorB, pResultX, pTempMatrix, Size);
//	start = clock();
//	MainLoop(pMatrix, pVectorB, pResultX, pTempMatrix, Size);
//	finish = clock();
//	duration = finish - start;
//	printf("Cramer`s rule + Gauss determinant program\n");
//	printf("\n Result Vector X: \n");
//	PrintVector(pResultX, Size);
//	printf("\n Time of execution: %f\n", duration);
//	ProcessTermination(pMatrix, pVectorB, pTempMatrix);
//}
