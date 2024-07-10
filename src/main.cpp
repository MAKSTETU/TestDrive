
#include <iostream>
#include <math.h>
#include <cmath>
#include <complex>
#include <string>
#include <fstream>
#include <fftw3.h>
#include <vector>
#include <algorithm>
#include <iterator>
# define M_PI           3.14159265358979323846  /* pi */
using namespace std;
const float FDIS = 1200000;

///////////////////////////////////////////////////////////////////////file reading function
vector<complex<double>> read(string name) {
	vector <float> mass;
	ifstream FILE(name, ios::binary | ios::in);
	float num;
	while (FILE.read(reinterpret_cast<char*>(&num), sizeof(float))) {
		mass.push_back(num);
	}
	FILE.close();
	vector <complex<double>> signal;//ОПОРНЫЙ СИГНАЛ
	for (int i = 0; i < mass.size(); i++)
	{
		complex<float> Complexvalue(mass.at(i), mass.at(i + 1));
		signal.push_back(Complexvalue);
		i++;
	}
	return signal;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Signal
{
public:
	void Signal::culc() {
		vector <complex<double>> signal1 = read("psp1cut.bin");//ОПОРНЫЙ СИГНАЛ
		vector <complex<double>> signal2 = read("psp1noised.bin");//ЗАШУМЛЕННАЯ РЕАЛИЗАЦИЯ СИГНАЛА

		int SIZE = signal1.size();
		vector<complex<double>> cpysgnl2(signal2);//КОПИЯ ЗАШУМЛЕННОГО СИГНАЛА
		float masmax[3][401];
		for (int f0 = -200; f0 <= 200; f0++) //ПОДСТРОЙКА ПО ЧАСТОТЕ
		{
			for (int k = 0; k < SIZE; k++)//ДОМНОЖЕНИЕ НА КОМПЛЕКСНУЮ ЭКСПОНЕНТУ
			{
				complex<double> e = (cos(2 * M_PI * f0 * k / FDIS), sin(2 * M_PI * f0 * k / FDIS));
				cpysgnl2[k] = signal2.at(k) * e;
			}

			vector <complex<double>> Fur1(signal1);//СОЗДАНИЕ ВЕКТОРА ДЛЯ БПФ ПЕРВОГО СИГНАЛА
			vector <complex<double>> Fur2(cpysgnl2);//СОЗДАНИЕ ВЕКТОРА ДЛЯ БПФ ВТОРОГО СИГНАЛА
			vector <complex<double>> Fur3(signal2);//СОЗДАНИЕ ВЕКТОРА ДЛЯ БПФ СВЕРТКИ ДВУХ СИГНАЛОВ
			vector <complex<double>> Fur4(Fur3);
			//СОЗДАНИЕ ПЛАНОВ ПФ
			fftw_plan p1 = fftw_plan_dft_1d(signal1.size(), (fftw_complex*)&signal1.at(0), (fftw_complex*)&Fur1.at(0), FFTW_FORWARD, FFTW_ESTIMATE);
			fftw_plan p2 = fftw_plan_dft_1d(signal2.size(), (fftw_complex*)&cpysgnl2.at(0), (fftw_complex*)&Fur2.at(0), FFTW_FORWARD, FFTW_ESTIMATE);
			fftw_plan p3 = fftw_plan_dft_1d(signal2.size(), (fftw_complex*)&Fur3.at(0), (fftw_complex*)&Fur4.at(0), FFTW_BACKWARD, FFTW_ESTIMATE);

			fftw_execute(p1);//ВЫПОЛНЕНИЕ БПФ ПЕРВОГО СИГНАЛА

			fftw_execute(p2);//ВЫПОЛНЕНИЕ БПФ ВТОРОГО 
			for (int j = 0; j < SIZE; j++) {
				Fur1.at(j) = conj(Fur1.at(j));
			}

			for (int i = 0; i < signal1.size(); i++) {
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


			vector<double> ::iterator MAX;

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
		cout << "временная задержка: " << max_index << " отсчетов" << endl;
		cout << "частотная отсройка: " << -Ftune << " Гц" << endl;

	}
	Signal() {
		
	}
	~Signal() {

		}

private:
	
};



int main()
{
	setlocale(LC_ALL, "Russian");
	Signal test;
	test.culc();

	return 0;
}
