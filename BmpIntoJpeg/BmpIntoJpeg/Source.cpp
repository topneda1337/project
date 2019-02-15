#include <iostream>
#include <fstream>
#include <cmath>
#include "Header.h"
#include <vector>
#include <algorithm>
#include <string>

//https://gist.github.com/ziggi/e15a95b9feac8f59c7c1

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

void ALG(int i, int j, int height, int width, std::vector<std::vector<double> > &DCT, std::vector<std::vector<YCbCr> > &arr, std::vector<int> &compANS, int flagg)
{
	std::vector<std::vector<double> > IMG(8, std::vector<double>(8, 0.0));
	std::vector<std::vector<double> > ANS(8, std::vector<double>(8, 0.0));

	if (j == 280)
		int f = 0;

	if (flagg == 1)
		for (int k = i; k < std::min(i + 8, height); k++)
			for (int m = j; m < std::min(j + 8, width); m++)
				IMG[k - i][m - j] = arr[k][m].Cb;
			
	if(flagg == 0)
		for (int k = i; k < std::min(i + 8, height); k++)
			for (int m = j; m < std::min(j + 8, width); m++)
				IMG[k - i][m - j] = arr[k][m].Cr;

	if (flagg == 2)
		for (int k = i; k < std::min(i + 8, height); k++)
			for (int m = j; m < std::min(j + 8, width); m++)
				IMG[k - i][m - j] = arr[k][m].Y;

	for (int k = 0; k < 8; k++)
		for (int m = k; m < 8; m++)
			std::swap(DCT[k][m], DCT[m][k]);

	mult(ANS, IMG, DCT);


	IMG = ANS;

	for (int k = 0; k < 8; k++)
		for (int m = k; m < 8; m++)
			std::swap(DCT[k][m], DCT[m][k]);

	mult(ANS, DCT, IMG);

	std::vector<std::vector<double> > QUANT(8, std::vector<double>(8, 0.0));

	double qcoeff = 2.0;
	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			QUANT[i][j] = (1 + (double)(1 + i + j) * qcoeff);

	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			ANS[i][j] /= QUANT[i][j];

	std::vector<unsigned char> comparr;
	for (int i = 0; i < 15; i++)
	{
		for (int j = std::max(i - 7, 0); j <= std::min(7, i); j++)
		{
			int x;
			if (i % 2 == 1)
				x = ANS[j][i - j];
			else
				x = ANS[i - j][j];
			if (i == 0 && j == 0)
			{
				x += 32768;
				comparr.push_back((unsigned char)((x / 256)));
				comparr.push_back((unsigned char)((x % 256)));
			}
			else
			{
				comparr.push_back((unsigned char)(x + 128));
				if (x > 127 || x < -128)
					int f = 229 + 7;
			}
		}
	}



	for (int i = 0; i < 65; i++)
	{
		int j = i + 1;
		bool flag = false;
		//1 2 3 4 5 i = 0 j = 5 1 2 3 4 5 5 i = 0 j = 5
		while (j < 65 && comparr[j] != comparr[j - 1])
		{
			flag = true;
			j++;
		}
		if (!flag)
		{
			j = i + 1;
			while (j < 65 && comparr[j] == comparr[j - 1])
				j++;
		}
		if (flag)
		{
			if (j == 65)
				j++;

			unsigned char cnt = j - i - 1;
			unsigned char x = cnt + 128;
			compANS.push_back(x);
			for (int k = i; k < j - 1; k++)
			{
				unsigned char c = comparr[k];
				compANS.push_back(c);
			}
			i = j - 2;
		}
		else
		{
			unsigned char cnt = j - i;
			unsigned char x = cnt;
			compANS.push_back(x);
			unsigned char c = comparr[i];
			compANS.push_back(c);
			i = j - 1;
		}
	}
}

bool COMPRESS(std::vector<std::vector<YCbCr> > &arr, std::vector<std::vector<RGBQUAD> > &in, std::string fileName, int width, int height)
{

	std::ofstream fileStream(fileName, std::ofstream::binary);

	if (!fileStream) {
		std::cout << "Error creating file '" << fileName << "'." << std::endl;
		return 0;
	}


	unsigned char X1 = height / 256, AnX1 = height % 256;
	unsigned char Y1 = width / 256, AnY1 = width % 256;

	std::vector<unsigned char> fullout;

	fullout.push_back(X1);
	fullout.push_back(AnX1);
	fullout.push_back(Y1);
	fullout.push_back(AnY1);


	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			arr[i][j].Y = (unsigned char)(0.299 * (double)in[i][j].rgbRed + 0.578 * (double)in[i][j].rgbGreen + 0.114 * (double)in[i][j].rgbBlue);
			arr[i][j].Cb = (int)(-0.1678 * (double)in[i][j].rgbRed - 0.3313 * (double)in[i][j].rgbGreen + 0.5 * (double)in[i][j].rgbBlue);
			arr[i][j].Cr = (int)(0.5 * (double)in[i][j].rgbRed - 0.4187 * (double)in[i][j].rgbGreen - 0.0813 * (double)in[i][j].rgbBlue);

			arr[i][j].Cb += 128;
			arr[i][j].Cr += 128;
			arr[i][j].Y = std::min(arr[i][j].Y, 255);
			arr[i][j].Y = std::max(arr[i][j].Y, 0);
			arr[i][j].Cb = std::min(arr[i][j].Cb, 255);
			arr[i][j].Cb = std::max(arr[i][j].Cb, 0);
			arr[i][j].Cr = std::min(arr[i][j].Cr, 255);
			arr[i][j].Cr = std::max(arr[i][j].Cr, 0);
			
			fullout.push_back(arr[i][j].Y);
		}
	}
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

	long long size1, size2, size0;

	size0 = width * height;

	for (int i = 0; i < height; i += 8)
	{
		for (int j = 0; j < width; j += 8)
		{
			std::vector<int> compANS;

			ALG(i, j, height, width, DCT, arr, compANS, 1);

			for (int ii = 0; ii < compANS.size(); ii++)
				fullout.push_back(compANS[ii]);
		}
	}

	size1 = fullout.size() - 4 - size0;

	for (int i = 0; i < height; i += 8)
	{
		for (int j = 0; j < width; j += 8)
		{
			std::vector<int> compANS;

			ALG(i, j, height, width, DCT, arr, compANS, 0);

			for (int ii = 0; ii < compANS.size(); ii++)
				fullout.push_back(compANS[ii]);


		}
	}
	size2 = fullout.size() - size1 - 4 - size0;

	short int XXX4 = (fullout.size() - 4) / 256 / 256, XX3 = ((fullout.size() - 4) / 256) % 256, YY3 = (fullout.size() - 4) % 256;

	fileStream.write((char*)(&XXX4), sizeof(XXX4));
	fileStream.write((char*)(&XX3), sizeof(XX3));
	fileStream.write((char*)(&YY3), sizeof(YY3));

	for (int i = 0; i < 4; i++)
		fileStream.write((char*)(&fullout[i]), sizeof(fullout[i]));
	
	fileStream.write((char*)(&size0), sizeof(size0));
	fileStream.write((char*)(&size1), sizeof(size1));
	fileStream.write((char*)(&size2), sizeof(size2));

	for (int i = 4; i < fullout.size(); i++)
		if (!fileStream.write((char*)(&fullout[i]), sizeof(fullout[i])))
			std::cout << "!!!!!!!!!!!";
	return true;
}


int main(int argc, char *argv[])
{
	std::string fileName = "a.bmp";

	std::cout << "Input filename:\n";

	std::cin >> fileName;

	fileName += ".bmp";

	// открываем файл
	std::ifstream fileStream(fileName, std::ifstream::binary);
	if (!fileStream) {
			std::cout << "Error opening file '" << fileName << "'." << std::endl;
		return 0;
	}

	std::vector<unsigned char> allheader;

	// заголовк изображения
	BITMAPFILEHEADER fileHeader;
	read(fileStream, fileHeader.bfType, sizeof(fileHeader.bfType));

	read(fileStream, fileHeader.bfSize, sizeof(fileHeader.bfSize));

	read(fileStream, fileHeader.bfReserved1, sizeof(fileHeader.bfReserved1));

	read(fileStream, fileHeader.bfReserved2, sizeof(fileHeader.bfReserved2));

	read(fileStream, fileHeader.bfOffBits, sizeof(fileHeader.bfOffBits));


	if (fileHeader.bfType != 0x4D42) {
		std::cout << "Error: '" << fileName << "' is not BMP file." << std::endl;
		return 0;
	}



	// информация изображения
	BITMAPINFOHEADER fileInfoHeader;
	read(fileStream, fileInfoHeader.biSize, sizeof(fileInfoHeader.biSize));
	allheader.push_back(fileInfoHeader.biSize);
	// bmp core
	if (fileInfoHeader.biSize >= 12) {
		read(fileStream, fileInfoHeader.biWidth, sizeof(fileInfoHeader.biWidth));

		read(fileStream, fileInfoHeader.biHeight, sizeof(fileInfoHeader.biHeight));

		read(fileStream, fileInfoHeader.biPlanes, sizeof(fileInfoHeader.biPlanes));

		read(fileStream, fileInfoHeader.biBitCount, sizeof(fileInfoHeader.biBitCount));

	}

	// получаем информацию о битности
	int colorsCount = fileInfoHeader.biBitCount >> 3;
	if (colorsCount < 3) {
		colorsCount = 3;
	}

	int bitsOnColor = fileInfoHeader.biBitCount / colorsCount;
	int maskValue = (1 << bitsOnColor) - 1;

	// bmp v1
	if (fileInfoHeader.biSize >= 40) {
		read(fileStream, fileInfoHeader.biCompression, sizeof(fileInfoHeader.biCompression));
		read(fileStream, fileInfoHeader.biSizeImage, sizeof(fileInfoHeader.biSizeImage));
		read(fileStream, fileInfoHeader.biXPelsPerMeter, sizeof(fileInfoHeader.biXPelsPerMeter));
		read(fileStream, fileInfoHeader.biYPelsPerMeter, sizeof(fileInfoHeader.biYPelsPerMeter));
		read(fileStream, fileInfoHeader.biClrUsed, sizeof(fileInfoHeader.biClrUsed));
		read(fileStream, fileInfoHeader.biClrImportant, sizeof(fileInfoHeader.biClrImportant));
	}

	// bmp v2
	fileInfoHeader.biRedMask = 0;
	fileInfoHeader.biGreenMask = 0;
	fileInfoHeader.biBlueMask = 0;

	if (fileInfoHeader.biSize >= 52) {
		read(fileStream, fileInfoHeader.biRedMask, sizeof(fileInfoHeader.biRedMask));
		read(fileStream, fileInfoHeader.biGreenMask, sizeof(fileInfoHeader.biGreenMask));
		read(fileStream, fileInfoHeader.biBlueMask, sizeof(fileInfoHeader.biBlueMask));
	}

	// если маска не задана, то ставим маску по умолчанию
	if (fileInfoHeader.biRedMask == 0 || fileInfoHeader.biGreenMask == 0 || fileInfoHeader.biBlueMask == 0) {
		fileInfoHeader.biRedMask = maskValue << (bitsOnColor * 2);
		fileInfoHeader.biGreenMask = maskValue << bitsOnColor;
		fileInfoHeader.biBlueMask = maskValue;
	}

	// bmp v3
	if (fileInfoHeader.biSize >= 56) {
		read(fileStream, fileInfoHeader.biAlphaMask, sizeof(fileInfoHeader.biAlphaMask));
	}
	else {
		fileInfoHeader.biAlphaMask = maskValue << (bitsOnColor * 3);
	}

	// bmp v4
	if (fileInfoHeader.biSize >= 108) {
		read(fileStream, fileInfoHeader.biCSType, sizeof(fileInfoHeader.biCSType));
		read(fileStream, fileInfoHeader.biEndpoints, sizeof(fileInfoHeader.biEndpoints));
		read(fileStream, fileInfoHeader.biGammaRed, sizeof(fileInfoHeader.biGammaRed));
		read(fileStream, fileInfoHeader.biGammaGreen, sizeof(fileInfoHeader.biGammaGreen));
		read(fileStream, fileInfoHeader.biGammaBlue, sizeof(fileInfoHeader.biGammaBlue));
	}

	// bmp v5
	if (fileInfoHeader.biSize >= 124) {
		read(fileStream, fileInfoHeader.biIntent, sizeof(fileInfoHeader.biIntent));
		read(fileStream, fileInfoHeader.biProfileData, sizeof(fileInfoHeader.biProfileData));
		read(fileStream, fileInfoHeader.biProfileSize, sizeof(fileInfoHeader.biProfileSize));
		read(fileStream, fileInfoHeader.biReserved, sizeof(fileInfoHeader.biReserved));
	}

	// проверка на поддерку этой версии формата
	if (fileInfoHeader.biSize != 12 && fileInfoHeader.biSize != 40 && fileInfoHeader.biSize != 52 &&
		fileInfoHeader.biSize != 56 && fileInfoHeader.biSize != 108 && fileInfoHeader.biSize != 124) {
		std::cout << "Error: Unsupported BMP format." << std::endl;
		return 0;
	}

	if (fileInfoHeader.biBitCount != 16 && fileInfoHeader.biBitCount != 24 && fileInfoHeader.biBitCount != 32) {
		std::cout << "Error: Unsupported BMP bit count." << std::endl;
		return 0;
	}

	if (fileInfoHeader.biCompression != 0 && fileInfoHeader.biCompression != 3) {
		std::cout << "Error: Unsupported BMP compression." << std::endl;
		return 0;
	}

	// rgb info


	//RGBQUAD **rgbInfo = new RGBQUAD*[fileInfoHeader.biHeight];
	std::vector<std::vector<RGBQUAD> > rgbInfo(fileInfoHeader.biHeight);

	for (unsigned int i = 0; i < fileInfoHeader.biHeight; i++) {
		rgbInfo[i].resize(fileInfoHeader.biWidth);
	}

	// определение размера отступа в конце каждой строки
	int linePadding = ((fileInfoHeader.biWidth * (fileInfoHeader.biBitCount / 8)) % 4);


	linePadding = (4 - linePadding) % 4;

	// чтение
	unsigned int bufer;

	for (unsigned int i = 0; i < fileInfoHeader.biHeight; i++) {
		for (unsigned int j = 0; j < fileInfoHeader.biWidth; j++) {
			read(fileStream, bufer, fileInfoHeader.biBitCount / 8);

			rgbInfo[i][j].rgbRed = bitextract(bufer, fileInfoHeader.biRedMask);
			rgbInfo[i][j].rgbGreen = bitextract(bufer, fileInfoHeader.biGreenMask);
			rgbInfo[i][j].rgbBlue = bitextract(bufer, fileInfoHeader.biBlueMask);
			rgbInfo[i][j].rgbReserved = bitextract(bufer, fileInfoHeader.biAlphaMask);
		}
		fileStream.seekg(linePadding, std::ios_base::cur);
	}



	std::vector<std::vector<YCbCr> > rgbInfo1(fileInfoHeader.biHeight);

	for (unsigned int i = 0; i < fileInfoHeader.biHeight; i++) {
		rgbInfo1[i].resize(fileInfoHeader.biWidth);
	}

	std::cout << "Input output file name:\n";

	std::string out;

	std::cin >> out;

	out += ".selfjpg";

	COMPRESS(rgbInfo1, rgbInfo, out, fileInfoHeader.biWidth, fileInfoHeader.biHeight);

	return 1;
}

unsigned char bitextract(const unsigned int byte, const unsigned int mask) {
	if (mask == 0) {
		return 0;
	}

	// определение количества нулевых бит справа от маски
	int
		maskBufer = mask,
		maskPadding = 0;

	while (!(maskBufer & 1)) {
		maskBufer >>= 1;
		maskPadding++;
	}

	// применение маски и смещение
	return (byte & mask) >> maskPadding;
}