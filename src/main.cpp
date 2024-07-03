
#include <iostream>
#include <math.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <fftw3.h>
#include <vector>
#include <algorithm>
#include <iterator>
# define M_PI           3.14159265358979323846  /* pi */
using namespace std;
const float FDIS = 1200000;

int main()
{
	///////////////////////////////////////////////////////////////////read first signal
	setlocale(LC_ALL, "Russian");

	vector <float> mas1;
	ifstream firstFILE("psp1cut.bin", ios::binary | ios::in);
	float num;
	while (firstFILE.read(reinterpret_cast<char*>(&num), sizeof(float))) {
		mas1.push_back(num);
	}
	firstFILE.close();

	vector <complex<double>> sgnl1;//ОПОРНЫЙ СИГНАЛ
	for (int i = 0; i < mas1.size(); i++)
	{
		complex<float> Complexvalue(mas1.at(i), mas1.at(i + 1));
		sgnl1.push_back(Complexvalue);
		i++;
	}
	////////////////////////////////////////////////////////////////// read second signal
	vector <float> mas2;

	ifstream FILE2("psp1noised.bin", ios::binary | ios::in);
	float num2;
	while (FILE2.read(reinterpret_cast<char*>(&num2), sizeof(float))) {
		mas2.push_back(num2);
	}
	FILE2.close();

	vector <complex<double>> sgnl2;//ЗАШУМЛЕННАЯ РЕАЛИЗАЦИЯ СИГНАЛА
	for (int i = 0; i < mas2.size(); i++)
	{
		complex<float> Complexvalue(mas2.at(i), mas2.at(i + 1));
		sgnl2.push_back(Complexvalue);
		i++;
	}
	
	/////////////////////////////////////////////////////////////

	int SIZE = sgnl1.size();
	vector<complex<double>> cpysgnl2(sgnl2);//КОПИЯ ЗАШУМЛЕННОГО СИГНАЛА
	float masmax[3][401];

	for (int f0 = -200; f0 <= 200; f0++) //ПОДСТРОЙКА ПО ЧАСТОТЕ
	{
		for (int k = 0; k < sgnl2.size(); k++)//ДОМНОЖЕНИЕ НА КОМПЛЕКСНУЮ ЭКСПОНЕНТУ
		{
			complex<double> e = (cos(2 * M_PI *f0 * k / FDIS), sin(2 * M_PI * f0 * k / FDIS));
			cpysgnl2[k] = sgnl2.at(k) * e;
		}

		vector <complex<double>> Fur1(sgnl1);//СОЗДАНИЕ ВЕКТОРА ДЛЯ БПФ ПЕРВОГО СИГНАЛА
		vector <complex<double>> Fur2(cpysgnl2);//СОЗДАНИЕ ВЕКТОРА ДЛЯ БПФ ВТОРОГО СИГНАЛА
		vector <complex<double>> Fur3(sgnl2);//СОЗДАНИЕ ВЕКТОРА ДЛЯ БПФ СВЕРТКИ ДВУХ СИГНАЛОВ
		vector <complex<double>> Fur4(Fur3);
		//СОЗДАНИЕ ПЛАНОВ ПФ
		fftw_plan p1 = fftw_plan_dft_1d(sgnl1.size(), (fftw_complex*)&sgnl1.at(0), (fftw_complex*)&Fur1.at(0), FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_plan p2 = fftw_plan_dft_1d(sgnl2.size(), (fftw_complex*)&cpysgnl2.at(0), (fftw_complex*)&Fur2.at(0), FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_plan p3 = fftw_plan_dft_1d(sgnl2.size(), (fftw_complex*)&Fur3.at(0), (fftw_complex*)&Fur4.at(0), FFTW_BACKWARD, FFTW_ESTIMATE);

		fftw_execute(p1);//ВЫПОЛНЕНИЕ БПФ ПЕРВОГО СИГНАЛА

		fftw_execute(p2);//ВЫПОЛНЕНИЕ БПФ ВТОРОГО 
		for (int j = 0; j < SIZE; j++) {
			Fur1.at(j) = conj(Fur1.at(j));
		}
	
		for (int i = 0; i < sgnl1.size(); i++) {
			Fur3.at(i) = Fur1.at(i) * Fur2.at(i);
		}
		fftw_execute(p3);//ВЫПОЛНЕНИЕ ОБРАТНОГО БПФ СВЕРТКИ СИГНАЛОВ
		fftw_destroy_plan(p1);
		fftw_destroy_plan(p2);
		fftw_destroy_plan(p3);

		vector <double> ABS(SIZE);
		for (int i = 0; i < SIZE; i++)
		{
			ABS.at(i) = abs(Fur4.at(i));
		}

	
		vector<double> :: iterator MAX;

		MAX = max_element(ABS.begin(), ABS.end());
		masmax[0][f0 + 200] = *MAX;//запись максимальеного значения корреляции
		masmax[1][f0 + 200] = distance(ABS.begin(), MAX);//запись индекса максимального значения корреляции
		masmax[2][f0 + 200] = f0;

	}

	double max_value = masmax[0][0];
	double max_index = 0;
	double Ftune = 0;
	for (int i = 1; i < 401; i++) {
		if (masmax[0][i] > max_value) {
			max_value = masmax[0][i];
			max_index = masmax[1][i];
			Ftune = masmax[2][i];
		}
	}
	//cout << "MAX_VALUE=" << max_value << endl;
	cout << "временная задержка" << max_index << " отсчетов" << endl;
	cout << "частотная отсройка" << -Ftune << " Гц" << endl;


	return 0;
}
