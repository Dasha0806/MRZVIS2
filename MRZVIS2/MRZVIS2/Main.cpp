// Лабораторная работа 2 по дисциплине МРЗвИС
// Выполнена студентками группы 821702
// БГУИР, Хвесеня Дарья, Дыкова Юлия

#include <iostream>
#include "Matrix.h"

using namespace std;

int main()
{
	setlocale(LC_ALL, "ru");

	int p, m, q, n;
	cout << "Введите размеры матрицы" << endl;
	cin >> p;
	cin >> m;
	cin >> q;
	cout << "Введите количество процессорных элементов" << endl;
	cin >> n;

	cout << endl;

	Matrix abcdef(p, m, q, n);

	abcdef.start();

	cout << "Матрица A:\n";
	abcdef.printMatrix(abcdef.a);
	cout << endl;

	cout << "Матрица B:\n";
	abcdef.printMatrix(abcdef.b);
	cout << endl;

	cout << "Матрица E:\n";
	abcdef.printNumbersList(abcdef.e);
	cout << endl << endl;

	cout << "Матрица G:\n";
	abcdef.printMatrix(abcdef.g);

	cout << "\n---------------------------------\n\n";

	cout << "Матрица С:\n";
	abcdef.printMatrix(abcdef.getC());
	cout << endl;

	cout << "Матрица D:\n";
	for (auto matrix : abcdef.getD()) {
		abcdef.printMatrix(matrix);
		cout << "--\n";
	}
	cout << endl;

	cout << "Матрица F:\n";
	for (auto matrix : abcdef.getF()) {
		abcdef.printMatrix(matrix);
		cout << "--\n";
	}
	cout << endl;

	return 0;
}