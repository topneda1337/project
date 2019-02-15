#include <iostream>
#include <fstream>
#include <cmath>
#include "Header.h"
#include <vector>
#include <algorithm>
#include <string>


void mult(std::vector<std::vector<double> > &ans, std::vector<std::vector<double> > &arr1, std::vector<std::vector<double> > &arr2, int n = 8)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			ans[i][j] = 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				ans[i][j] += arr1[i][k] * arr2[k][j];
}

void InvMatr(std::vector<std::vector<double> > &DCT, std::vector<std::vector<double> > &ANS, int n = 8)
{
	std::vector<std::vector<double> > Count(n, std::vector<double>(2 * n, 0.0));
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			Count[i][j] = DCT[i][j];
	for (int i = 0; i < n; i++)
		Count[i][i + n] = 1.0;

	for (int i = 0; i < n; i++)
	{
		int k = i;
		while (k < n && (Count[k][i] < EPS) && (Count[k][i] > -EPS))
			k++;
		if (k == n)
			return;
		std::swap(Count[i], Count[k]);
		double X = Count[i][i];
		for (int j = 0; j < 2 * n; j++)
			Count[i][j] /= X;
		for (int j = i + 1; j < n; j++)
		{
			double P = Count[j][i];
			Count[j][i] = 0;
			for (int t = i + 1; t < 2 * n; t++)
				Count[j][t] -= P * Count[i][t];
		}
	}

	for (int i = n - 1; i > 0; i--)
	{
		for (int j = i - 1; j >= 0; j--)
		{
			double X = Count[j][i];
			for (int k = 0; k < 2 * n; k++)
				Count[j][k] -= X * Count[i][k];
		}
	}

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			ANS[i][j] = Count[i][j + n];
	return;
}

void RevALG(std::vector<std::vector<int> > &in, std::vector<std::vector<int> > &out)
{
	std::vector<std::vector<double> > QUANT(8, std::vector<double>(8, 0.0));

	double qcoeff = 2.0;
	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			QUANT[i][j] = (1 + (double)(1 + i + j) * qcoeff);

	std::vector<std::vector<double> > ARR(8, std::vector<double>(8, 0.0));

	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			ARR[i][j] = (double)in[i][j] * QUANT[i][j];

	std::vector<std::vector<double> > DCT(8, std::vector<double>(8, 0.0));
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			if (i == 0)
				DCT[i][j] = 1.0 / sqrt(8);
			else
				DCT[i][j] = sqrt(2.0 / 8.0) * cos((double)((2 * j + 1) * i) * PI / (double)(2 * 8));
		}
	}

	std::vector<std::vector<double> > INVDCT(8, std::vector<double>(8, 0.0));
	std::vector<std::vector<double> > tmp(8, std::vector<double>(8, 0.0));


	for (int k = 0; k < 8; k++)
		for (int m = k; m < 8; m++)
			std::swap(DCT[k][m], DCT[m][k]);

	InvMatr(DCT, INVDCT);

	mult(tmp, ARR, INVDCT);

	ARR = tmp;

	for (int k = 0; k < 8; k++)
		for (int m = k; m < 8; m++)
			std::swap(DCT[k][m], DCT[m][k]);

	InvMatr(DCT, INVDCT);

	mult(tmp, INVDCT, ARR);

	ARR = tmp;

	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
		{
			out[i][j] = ARR[i][j];
		}
	return;
}



void DECOMPRESS(std::vector<unsigned char> &arr, std::vector<std::vector<RGBQUAD> > &out, int height, int width, long long size0, long long size1, long long size2)
{
	out.resize(height, std::vector<RGBQUAD>(width));

	std::vector<std::vector<YCbCr> > temparr(height, std::vector<YCbCr>(width));

	out.resize(height, std::vector<RGBQUAD>(width));
	std::vector<int> first, second, third;
	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
			third.push_back(arr[i * width + j]);
	size0 = height * width;

	for (int i = size0; i < size0 + size1; i++)
	{
		if (arr[i] >= 128)
		{
			int c = arr[i] - 128;
			for (int j = i + 1; j <= i + c; j++)
				first.push_back(arr[j]);
			i = i + c;
		}
		else
		{
			int c = arr[i];
			for (int j = 0; j < arr[i]; j++)
				first.push_back(arr[i + 1]);
			i++;
		}
	}
	for (int i = size0 + size1; i < size0 + size1 + size2; i++)
	{
		if (arr[i] >= 128)
		{
			int c = arr[i] - 128;
			for (int j = i + 1; j <= i + c; j++)
				second.push_back(arr[j]);
			i = i + c;
		}
		else
		{
			int c = arr[i];
			for (int j = 0; j < arr[i]; j++)
				second.push_back(arr[i + 1]);
			i++;
		}
	}
	for (int i = 0; i < third.size(); i++)
		temparr[i / width][i % width].Y = third[i];
	int tempH = 0, tempW = 0;

	for (int i = 0; i + 64 < first.size(); i += 65)
	{
		int cnttemp = 0;
		std::vector<std::vector<int> > INV(8, std::vector<int>(8, 0.0));
		for (int k = 0; k < 15; k++)
		{
			for (int j = std::max(k - 7, 0); j <= std::min(7, k); j++)
			{
				int x;
				if (cnttemp != 0)
				{
					if (k % 2 == 1)
						INV[j][k - j] = first[i + cnttemp] - 128;
					else
						INV[k - j][j] = first[i + cnttemp] - 128;

				}
				else
				{
					INV[0][0] = (int)first[i] * 256 + (int)first[i + 1];
					INV[0][0] -= 32768;
					cnttemp++;
				}
				cnttemp++;
			}
		}

		std::vector<std::vector<int> > tmp(8, std::vector<int>(8, 0.0));

		tmp = INV;

		RevALG(INV, tmp);

		INV = tmp;

		for (int j = tempH; j < std::min(tempH + 8, height); j++)
			for (int k = tempW; k < std::min(tempW + 8, width); k++)
				temparr[j][k].Cb = INV[j - tempH][k - tempW];



		tempW += 8;
		if (tempW >= width)
		{
			tempW = 0;
			tempH += 8;
		}
	}
	tempH = 0, tempW = 0;
	for (int i = 0; i + 64 < second.size(); i += 65)
	{

		if (tempW == 280)
			int f = 0;
		int cnttemp = 0;
		std::vector<std::vector<int> > INV(8, std::vector<int>(8, 0.0));
		for (int k = 0; k < 15; k++)
		{
			for (int j = std::max(k - 7, 0); j <= std::min(7, k); j++)
			{
				int x;
				if (cnttemp != 0)
				{
					if (k % 2 == 1)
						INV[j][k - j] = second[i + cnttemp] - 128;
					else
						INV[k - j][j] = second[i + cnttemp] - 128;
				}
				else
				{
					INV[0][0] = (int)second[i] * 256 + (int)second[i + 1];
					INV[0][0] -= 32768;
					cnttemp++;
				}
				cnttemp++;
			}
		}

		std::vector<std::vector<int> > tmp(8, std::vector<int>(8, 0.0));

		tmp = INV;

		RevALG(INV, tmp);

		INV = tmp;

		for (int j = tempH; j < std::min(tempH + 8, height); j++)
			for (int k = tempW; k < std::min(tempW + 8, width); k++)
				temparr[j][k].Cr = INV[j - tempH][k - tempW];
		tempW += 8;
		if (tempW >= width)
		{
			tempW = 0;
			tempH += 8;
		}
	}

	std::vector<std::vector<double> > coeffmatr(3, std::vector<double>(3, 0.0));
	std::vector<std::vector<double> > counter(3, std::vector<double>(3, 0.0));

	coeffmatr = { { 0.299, 0.578, 0.114 },{ -0.1678, -0.3313, 0.5 },{ 0.5, -0.4187, -0.0813 } };

	InvMatr(coeffmatr, counter, 3);


	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			double Y = temparr[i][j].Y;
			double Cb = temparr[i][j].Cb - 128;
			double Cr = temparr[i][j].Cr - 128;

			int B = (double)Y + 1.402*(double)Cr;
			int G = (double)Y - 0.34414*(double)Cb - 0.71414*(double)Cr;
			int R = (double)Y + 1.772*(double)Cb;
			R = std::min(R, 255);
			R = std::max(R, 0);
			G = std::min(G, 255);
			G = std::max(G, 0);
			B = std::min(B, 255);
			B = std::max(B, 0);
			out[i][j].rgbRed = R;
			out[i][j].rgbGreen = G;
			out[i][j].rgbBlue = B;
		}
	}

	return;
}


int main(int argc, char *argv[])
{
	std::string fileName = "b";

	std::cout << "Input filename:\n";
	std::cin >> fileName;

	fileName += ".selfjpg";

	std::ifstream fileStream(fileName, std::ifstream::binary);
	if (!fileStream) {
		std::cout << "Error opening file '" << fileName << "'." << std::endl;
		return 0;
	}

	std::vector<unsigned char> ARR;

	unsigned char c;

	short int X1 = 0, Y1 = 0, Z1 = 0;
	fileStream.read((char*)(&X1), sizeof(X1));
	fileStream.read((char*)(&Y1), sizeof(Y1));
	fileStream.read((char*)(&Z1), sizeof(Z1));

	unsigned char X2, XX2;
	unsigned char Y2, YY2;

	fileStream.read((char*)(&X2), sizeof(X2));
	fileStream.read((char*)(&XX2), sizeof(XX2));
	fileStream.read((char*)(&Y2), sizeof(Y2));
	fileStream.read((char*)(&YY2), sizeof(YY2));

	int height = (int)X2 * 256 + (int)XX2;
	int width = (int)Y2 * 256 + (int)YY2;

	long long size0 = 0;
	long long size1 = 0;
	long long size2 = 0;


	fileStream.read((char*)(&size0), sizeof(size0));
	fileStream.read((char*)(&size1), sizeof(size1));
	fileStream.read((char*)(&size2), sizeof(size2));

	for (int i = 0; i < ((int)X1 * 256 * 256 + (int)Y1 * 256 + (int)Z1); i++)
	{
		if (!fileStream.read((char*)(&c), sizeof(c)))
			std::cout << "!!!!!!!!";
		ARR.push_back(c);
	}

	
	std::vector<std::vector<RGBQUAD> > OUTarr;

	std::string OutfileName = "tst2";

	std::cout << "Input Output file name:\n";
	std::cin >> OutfileName;

	OutfileName += ".bmp";

	std::ofstream stream(OutfileName, std::ofstream::binary);

	if (!stream)
	{
		std::cout << "Error creating file '" << fileName << "'." << std::endl;
		return 0;
	}

	DECOMPRESS(ARR, OUTarr, height, width, size0, size1, size2);

	
	unsigned char file[14] = {
		'B','M', // magic
		0,0,0,0, // size in bytes
		0,0, // app data
		0,0, // app data
		40 + 14,0,0,0 // start of data offset
	};
	unsigned char info[40] = {
		40,0,0,0, // info hd size
		0,0,0,0, // width
		0,0,0,0, // heigth
		1,0, // number color planes
		24,0, // bits per pixel
		0,0,0,0, // compression is none
		0,0,0,0, // image bits size
		0x13,0x0B,0,0, // horz resoluition in pixel / m
		0x13,0x0B,0,0, // vert resolutions (0x03C3 = 96 dpi, 0x0B13 = 72 dpi)
		0,0,0,0, // #colors in pallete
		0,0,0,0, // #important colors
	};

	int w = OUTarr[0].size();
	int h = OUTarr.size();

	int padSize = (4 - (w * 3) % 4) % 4;
	int sizeData = w*h * 3 + h*padSize;
	int sizeAll = sizeData + sizeof(file) + sizeof(info);

	file[2] = (unsigned char)(sizeAll);
	file[3] = (unsigned char)(sizeAll >> 8);
	file[4] = (unsigned char)(sizeAll >> 16);
	file[5] = (unsigned char)(sizeAll >> 24);

	info[4] = (unsigned char)(w);
	info[5] = (unsigned char)(w >> 8);
	info[6] = (unsigned char)(w >> 16);
	info[7] = (unsigned char)(w >> 24);

	info[8] = (unsigned char)(h);
	info[9] = (unsigned char)(h >> 8);
	info[10] = (unsigned char)(h >> 16);
	info[11] = (unsigned char)(h >> 24);

	info[20] = (unsigned char)(sizeData);
	info[21] = (unsigned char)(sizeData >> 8);
	info[22] = (unsigned char)(sizeData >> 16);
	info[23] = (unsigned char)(sizeData >> 24);

	stream.write((char*)file, sizeof(file));
	stream.write((char*)info, sizeof(info));

	unsigned char pad[3] = { 0, 0, 0 };

	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			unsigned char pixel[3];
			pixel[0] = OUTarr[i][j].rgbRed;
			pixel[1] = OUTarr[i][j].rgbGreen;
			pixel[2] = OUTarr[i][j].rgbBlue;
	
			stream.write((char*)pixel, 3);
		}
		stream.write((char*)pad, padSize);
	}
	
	
	return 1;
}
