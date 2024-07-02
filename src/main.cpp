
// metra.cpp : Ýòîò ôàéë ñîäåðæèò ôóíêöèþ "main". Çäåñü íà÷èíàåòñÿ è çàêàí÷èâàåòñÿ âûïîëíåíèå ïðîãðàììû.
//
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
constexpr complex<float> _i = { 0.0, 1.0 }; // erunda
const float FDIS = 1200000;

int main()
{
	///////////////////////////////////////////////////////////////////read first signal

	vector <float> mas1;

	ifstream firstFILE("psp1cut.bin", ios::binary | ios::in);
	float num;
	while (firstFILE.read(reinterpret_cast<char*>(&num), sizeof(float))) {
		mas1.push_back(num);
	}
	firstFILE.close();

	vector <complex<float>> sgnl1;//signal 1
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

	vector <complex<float>> sgnl2;//signal 2
	for (int i = 0; i < mas2.size(); i++)
	{
		complex<float> Complexvalue(mas2.at(i), mas2.at(i + 1));
		sgnl2.push_back(Complexvalue);
		i++;
	}
	
	//cout << sgnl2.size() << sgnl2.front() << endl << sgnl2.at(1) << endl << sgnl2.at(2)<<endl<<sgnl2.back();
	/////////////////////////////////////////////////////////////

	int SIZE = sgnl1.size();
	vector<complex<float>> cpysgnl2(sgnl2);//copy of signal 2
	float masmax[3][401];

	for (int f0 = -200; f0 <= 200; f0++) //ÏÎÄÑÒÐÎÉÊÀ ÏÎ ×ÀÑÒÎÒÅ
	{
		for (int k = 0; k < sgnl2.size(); k++)//ÄÎÌÍÎÆÅÍÈÅ ÍÀ ÊÎÌÏËÅÊÑÍÓÞ ÝÊÑÏÎÍÅÍÒÓ
		{
			complex<float> e = (cos(2 * M_PI *f0 * k / FDIS), sin(2 * M_PI * f0 * k / FDIS));
			cpysgnl2[k] = sgnl2.at(k) * e;
		}

		vector <complex<float>> Fur1(sgnl1);//ft vector signal1
		vector <complex<float>> Fur2(cpysgnl2);//ft vector cpysgnl2
		vector <complex<float>> Fur3(sgnl2);//vector S1*S2
		//fft plan
		fftw_plan p1 = fftw_plan_dft_1d(sgnl1.size(), (fftw_complex*)&sgnl1.at(0), (fftw_complex*)&Fur1.at(0), FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_plan p2 = fftw_plan_dft_1d(sgnl2.size(), (fftw_complex*)&cpysgnl2.at(0), (fftw_complex*)&Fur2.at(0), FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_plan p3 = fftw_plan_dft_1d(sgnl2.size(), (fftw_complex*)&Fur3.at(0), (fftw_complex*)&Fur3.at(0), FFTW_BACKWARD, FFTW_ESTIMATE);

		fftw_execute(p1);//fft signal 1

		fftw_execute(p2);//fft copsgnl2

		for (int i = 0; i < sgnl1.size(); i++) {
			Fur3.at(i) = Fur1.at(i) * Fur2.at(i);
		}
		fftw_execute(p3);//fft S1*S2
		fftw_destroy_plan(p1);
		fftw_destroy_plan(p2);
		fftw_destroy_plan(p3);

		vector <float> ABS(SIZE);// vector ABS of fft S1*S2
		for (int i = 0; i < SIZE; i++)
		{
			ABS.at(i) = abs(Fur3.at(i));
		}

	
		vector<float> :: iterator MAX;

		MAX = max_element(ABS.begin(), ABS.end());
		masmax[0][f0 + 200] = *MAX;//max value
		masmax[1][f0 + 200] = distance(ABS.begin(), MAX);//index of max value
		masmax[2][f0 + 200] = f0;

	}

	float max_value = masmax[0][0];
	float max_index = 0;
	float Ftune = 0;
	for (int i = 1; i < 401; i++) {
		if (masmax[0][i] > max_value) {
			max_value = masmax[0][i];
			max_index = masmax[1][i];
			Ftune = masmax[2][i];
		}
	}
	cout << "MAX_VALUE=" << max_value << endl;
	cout << "MAX_INDEX=" << max_index << endl;
	cout << "Ftune=" << Ftune << endl;

	

	return 0;
}
