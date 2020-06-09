// ������������ ������ 2 �� ���������� ������
// ��������� ����������� ������ 821702
// �����, ������� �����, ������ ����


#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <thread>
#include <mutex>

using namespace std;

class Matrix
{
	//int p, m, q; ������� ������ A, B, E, G (pxm, mxq, 1xm, pxq) 
	int p, m, q;

	// ����
	int rank;
	
	// ����� ���������� ��������
	int multiplicationTime = 1;
	int additionTime = 1;
	int subtractionTime = 1;
	int comparisonTime = 1;

	// ���������� ������� ��������
	int numberOfMultiplications;
	int numberOfAdditions;
	int numberOfSubtractions;
	int numberOfComparisons;

	// ����� ������ ��������
	int consistentTime;
	int parallelTime;
	int averageTime;

	// ���������
	int acceleration;

	// �������������
	int efficiency;

	// ������� ����� ������
	int averageLength;

	// ��������
	int discrepancy;

	// ���������� ������������ ���������
	int n ;
	
	// ��������� ������������ ���� ����� �������.
	double procSetMultiplication(vector <double> numbers);

	// ��������� ������������ ���� ����� ������� � ������� (1 - x1) * (1 - x2)...(1 - xn).
	double procSetMultiplicationOneMinusNumber(vector <double> numbers);

	// ����� ���������.
	double sup(const double& firstNumber, const double& secondNumber);

	// ����� ������� �� ���� �����
	double minOfTwo(const double& firstNumber, const double& secondNumber);
	
	// ����� �������� �� ���� �����
	double maxOfTwo(const double& firstNumber, const double& secondNumber);

	// ��������� ��������� ����� x � ���������� [fMin;fMax]
	double fRandomNumber(double fMin, double fMax);

	// ���������� ������� �� ���������� ����������
	vector<vector<double>> generateMatrix(int firstSize, int secondSize);

public:

	// �������
	vector <double> e;
	vector <vector <double>> a, b, c, g;
	vector <vector <vector <double>>> d, f;

	Matrix();
	Matrix(const int& p, const int& m, const int& q, const int& n);

	void start();

	void generate();

	void setP(int p);
	void setM(int m);
	void setQ(int q);

	void setA(vector <vector <double>> a);
	void setB(vector <vector <double>> b);
	void setE(vector <double> e);
	void setG(vector <vector <double>> g);

	vector <vector <double>> getC();
	vector <vector <vector <double>>> getD();
	vector <vector <vector <double>>> getF();

	void printNumbersList(vector <double> numbers);
	void printMatrix(vector<vector <double>> matrix);

	void procMatrixC();
	void procMatrixD();
	void procMatrixF();
};


