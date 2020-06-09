// ������������ ������ 2 �� ���������� ������
// ��������� ����������� ������ 821702
// �����, ������� �����, ������ ����

#include <iostream>
#include "Matrix.h"

using namespace std;

int main()
{
	setlocale(LC_ALL, "ru");

	int p, m, q, n;
	cout << "������� ������� �������" << endl;
	cin >> p;
	cin >> m;
	cin >> q;
	cout << "������� ���������� ������������ ���������" << endl;
	cin >> n;

	cout << endl;

	Matrix abcdef(p, m, q, n);

	abcdef.start();

	cout << "������� A:\n";
	abcdef.printMatrix(abcdef.a);
	cout << endl;

	cout << "������� B:\n";
	abcdef.printMatrix(abcdef.b);
	cout << endl;

	cout << "������� E:\n";
	abcdef.printNumbersList(abcdef.e);
	cout << endl << endl;

	cout << "������� G:\n";
	abcdef.printMatrix(abcdef.g);

	cout << "\n---------------------------------\n\n";

	cout << "������� �:\n";
	abcdef.printMatrix(abcdef.getC());
	cout << endl;

	cout << "������� D:\n";
	for (auto matrix : abcdef.getD()) {
		abcdef.printMatrix(matrix);
		cout << "--\n";
	}
	cout << endl;

	cout << "������� F:\n";
	for (auto matrix : abcdef.getF()) {
		abcdef.printMatrix(matrix);
		cout << "--\n";
	}
	cout << endl;

	return 0;
}