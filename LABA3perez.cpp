#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>


using namespace std;

const unsigned int str = 150;
const unsigned int stolb = 14;
double h = 0.1;
int a = 0, b = 11;


double my_lagrange(int a, int b, double x, double x_arr[], double* z) {
	double sum = 0;
	for (int i = a; i <= b; i++) {
		double l = 1;
		for (int j = a; j <= b; j++)
			if (j != i)
				l *= (x - x_arr[j]) / (x_arr[i] - x_arr[j]);
		sum += z[i] * l;
	}
	return sum;
}



double d1(double m0, double m1) {
	return m1 - m0;

}
double d2(double m0, double m1, double m2) {

	return m2 - 2 * m1 + m0;
}

double d3(double m0, double m1, double m2, double m3) {

	return m3 - 3 * m2 + 3 * m1 - m0;
}

double d4(double m0, double m1, double m2, double m3, double m4) {

	return m4 - 4 * m3 + 6 * m2 - 4 * m1 + m0;
}
double d5(double m0, double m1, double m2, double m3, double m4, double m5) {

	return m5 - 5 * m4 + 10 * m3 - 10 * m2 + 5 * m1 - m0;
}
double d6(double m0, double m1, double m2, double m3, double m4, double m5, double m6) {

	return m6 - 6 * m5 + 15 * m4 - 20 * m3 + 15 * m2 - 6 * m1 + m0;
}

double New1(double* z, double x) {
	double sum = z[0];
	int n = 6;
	double x0 = 1, h = 1;
	double d[6] = { d1(z[0],z[1]),d2(z[0],z[1],z[2]),d3(z[0],z[1],z[2],z[3]),d4(z[0],z[1],z[2],z[3],z[4]),d5(z[0],z[1],z[2],z[3],z[4],z[5]),d6(z[0],z[1],z[2],z[3],z[4],z[5],z[6]) };
	double q = (x - x0) / h;
	sum += q * d[0];
	for (int i = 1; i < n; i++) {

		q *= ((x - x0) / h - i) / (i + 1);
		
		sum += q * d[i];

	}

	return sum;
}
double New2(double* z, double x) {
	double sum = z[6];
	int n = 6;
	double x0 = 7, h = 1;
	double d[6] = { d1(z[5],z[6]),d2(z[4],z[5],z[6]),d3(z[3],z[4],z[5],z[6]),d4(z[2],z[3],z[4],z[5],z[6]),d5(z[1],z[2],z[3],z[4],z[5],z[6]),d6(z[0],z[1],z[2],z[3],z[4],z[5],z[6]) };
	double q = (x - x0) / h;
	sum += q * d[0];
	for (int i = 1; i < n; i++) {

		q *= ((x - x0) / h + i) / (i + 1);
		sum += q * d[i];

	}
	return sum;
}


//Степ полиномом
double Steppolinom(double x, int s, double** koef)
{
	double y = 0.0;
	for (int i = 0; i <= s; i++)  y += koef[i][0] * pow(x, i);
	return y;
}
//обратн матрица
void inversion(double** A, int N)
{
	double temp;
	double** E = new double*[N];
	for (int i = 0; i < N; i++)
		E[i] = new double[N];
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			E[i][j] = 0.0;
			if (i == j)
				E[i][j] = 1.0;
		}
	for (int k = 0; k < N; k++)
	{
		temp = A[k][k];
		for (int j = 0; j < N; j++)
		{
			A[k][j] /= temp;
			E[k][j] /= temp;
		}
		for (int i = k + 1; i < N; i++)
		{
			temp = A[i][k];
			for (int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int k = N - 1; k > 0; k--)
	{
		for (int i = k - 1; i >= 0; i--)
		{
			temp = A[i][k];
			for (int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			A[i][j] = E[i][j];
	for (int i = 0; i < N; i++)
		delete[] E[i];

	delete[] E;
}
//Перемножение матриц. f-столб второй
void multiply(double** A, double** B, double** C, int c, int d, int f)
{
	for (int i = 0; i < c; i++)
	{
		for (int j = 0; j < f; j++)
		{
			for (int r = 0; r < d; r++)
				C[i][j] += A[i][r] * B[r][j];
		}
	}
}

//Печать массива
void print(double** temp, int n, int k)
{
	for (int i = 0; i < n; i++)
	{
		cout << i << "\t |";
		for (int j = 0; j < k; j++) cout << temp[i][j] << "\t |";
		cout << endl;
	}
}

int main() {

	double ip;
	int j;
	string line1;
	double** temp = new double*[150];
	for (int i = 0; i < 150; i++)
		temp[i] = new double[14];

	for (int i = 0; i < 150; i++)
		for (int j = 0; j < 14; j++)
			temp[i][j] = 0;
	int n = 0;
	ifstream file;
	file.open("17_Сочи.csv");
	if (!file.is_open())
		cout << "Error" << endl;
	else {
		while (!file.eof())
		{
			getline(file, line1);
			if (n != 0)
			{
				for (long long unsigned int i = 0; i < line1.size(); i++)
				{
					if (line1[i] == ',') line1[i] = ' ';
				};
				stringstream ss(line1);
				j = 0;
				while ((ss >> ip) && j < 14)
				{
					temp[n - 1][j] = ip;
					j++;
				}
			}
			n++;
		}
		n = 150;
	}


	fstream fin;
	fin.open("17_Сочи.csv", ios::in);

	double** y_arr = new double*[str];
	for (int count = 0; count < str; count++)
		y_arr[count] = new double[stolb];

	string line, word;
	getline(fin, line);
	for (int count_row = 0; count_row < str; count_row++) {
		getline(fin, line);
		stringstream s(line);
		for (int count_column = 0; count_column < stolb; count_column++) {
			getline(s, word, ',');
			y_arr[count_row][count_column] = stof(word);
		}
	}

	ofstream fout;
	fout.open("Lag.txt");

	int c;
	cout << "Viberite stolbik" << endl;
	cin >> c;
	int count9 = 0;
	for (int i = 0; i < 150; i++) {
		if (y_arr[i][c] > 900)
			count9++;
	}

	double* z = new double[150 - count9];

	int p = 0;
	for (int i = 0; i < 150; i++) {
		if (y_arr[i][c] < 999.9) {
			z[p] = y_arr[i][c];
			//cout << z[p] << endl;
			p++;

		}
	}

	double x_arr[str];


	for (int l = 0; l < str; l++) {
		x_arr[l] = l;
	}
	for (double x = a; x <= b; x += h)
	{
		cout << "f(" << x << ") = " << my_lagrange(a, b, x, x_arr, z) << endl;
		fout << x << " " << my_lagrange(a, b, x, x_arr, z) << endl;
	}


	////////////////////////////////////////////////////////////////////////////////////////////////



	ofstream fout54;
	fout54.open("Lag2.txt");
	for (int l = 0; l < str; l++) {
		x_arr[l] = l;
	}
	for (int j = a; j <= b; j++)
	{
		//cout << x_arr[j] << endl;
		fout54 << x_arr[j] << " " << z[j] << endl;
	}



	////////////////////////////////////////////////////////////////////////////////////////////////


	ofstream fout3;
	fout3.open("s.txt");
	for (int x = 0; x < 265; x++) {
		fout3 << x << " " << z[x] << endl;
	}
	cout << endl << endl;
	ofstream fout1;
	fout1.open("New1.txt");
	int a1 = 1, b1 = 6;
	for (double x = a1; x <= b1; x += 0.1) {
		cout << "f(" << x << ") = "
			<< New1(z, x) << endl;
		fout1 << x << " " << New1(z, x) << endl;
	}


	cout << endl << endl;
	ofstream fout2;
	fout2.open("New2.txt");
	int a2 = 1, b2 = 6;
	for (double x = b2; x >= a2; x -= 0.1) {
		cout << "f(" << x << ") = "
			<< New2(z, x) << endl;
		fout2 << x << " " << New2(z, x) << endl;
	}

	//5
	ofstream file2;
	file2.open("sp.txt");
	int s, month = c;
	cout << "vvedite stepen polinoma (do 5): ";
	cin >> s;


	vector<double>  temperature;
	vector<int> stroki;
	int nstrok = 0;
	for (int i = 0; i < 150; i++) if (temp[i][month] != 999.9) { temperature.push_back(temp[i][month]); nstrok++; stroki.push_back(nstrok); }
	cout << "kol-vo prav temp: " << temperature.size() << endl;
	int n1 = stroki.size();//кол-во переменн
	//прав ч ур-я. матрица зн-й y-температур
	double** secondmat = new double*[s + 1];
	for (int i = 0; i < s + 1; i++) secondmat[i] = new double[1];//матрица решений
	for (int i = 0; i < s + 1; i++) secondmat[i][0] = 0;
	for (int i = 0; i < s + 1; i++)
	{
		for (int j = 0; j < n1; j++) secondmat[i][0] += temperature[j] * pow(stroki[j], i);
	}
	//большая матрица сумм произведений х в степенях
	double** firstmat = new double*[s + 1];
	for (int i = 0; i < s + 1; i++) firstmat[i] = new double[s + 1];
	for (int i = 0; i < s + 1; i++)
	{
		for (int j = 0; j < s + 1; j++)
		{
			if (i == 0 && j == 0) firstmat[i][j] = n1;
			else { for (int k = 0; k < n1; k++) firstmat[i][j] += pow(stroki[k], i + j); }
		}
	}
	//искомая матрица коэф-в
	double** koef = new double*[s + 1];
	for (int i = 0; i < s + 1; i++) koef[i] = new double[1];
	for (int i = 0; i < s + 1; i++) koef[i][0] = 0.0;

	inversion(firstmat, s + 1);

	multiply(firstmat, secondmat, koef, s + 1, s + 1, 1);
	cout << "matrich. koef.: " << endl;
	print(koef, s + 1, 1);
	cout << endl;
	int k = stroki[n1 - 1];
	for (double i = 0; i < k; i += 0.1)
		file2 << i << "\t  " << Steppolinom(i, s, koef) << endl;
	file2.close();

	fout.close();
	fout1.close();
	fout2.close();
	fout3.close();
	fout54.close();
	delete[] * y_arr;
	delete[] z;
	return 0;
}
