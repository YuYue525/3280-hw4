#include "bmp.h"		//	Simple .bmp library
#include<iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>

using namespace std;
#define PI 3.14159265358979
int QuantizationMatrix[8][8] = {
	{3, 5, 7, 9, 11, 13, 15, 17}, 
	{5, 7, 9, 11, 13, 15, 17, 19},
	{7, 9, 11, 13, 15, 17, 19, 21},
	{9, 11,	13,	15,	17,	19,	21,	23},
	{11, 13, 15, 17, 19, 21, 23, 25},
	{13, 15, 17, 19, 21, 23, 25, 27},
	{15, 17, 19, 21, 23, 25, 27, 29},
	{17, 19, 21, 23, 25, 27, 29, 31}
	};


int main(int argc, char** argv)
{
	if (argc != 3)
	{
		cout << "Arguments prompt: dct.exe <img_path> <apply_idct>" << endl;
		return 0;
	}
	string imgPath = argv[1];
	bool need_idct = stoi(argv[2]);

	//! read input image
	Bitmap s_img(imgPath.c_str());
	int rows = s_img.getHeight(), cols = s_img.getWidth();
	cout << "Apply DCT on image ("<<rows<<", "<<cols<< ")." << endl;
	
	//! preprocess by shifting pixel values by 128
	//TODO
    int f[256][256] = {0};
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            unsigned char tmp = 0;
            s_img.getPixel(j, i, tmp);
            f[j][i] = (int)tmp - 128;
        }
    }
	
	//! 2D DCT for every 8x8 block (assume that the input image resolution is fixed to 256)	
	// The quantized coefficients should be stored into 'coeffArray'
	double coeffArray[256][256]={0};
	int blockRow = rows / 8, blockCol = cols / 8;
    
	for (int i = 0; i < blockRow; i++)
	{
		for (int j = 0; j < blockCol; j++)
		{
			int xpos = j*8, ypos = i*8;
			//! apply DCT on block_ij (basic requirement)
			//TODO
            double Fr[8][8] = {0};
            
            for (int v = 0; v < 8; v++)
            {
                for (int u = 0; u < 8; u++)
                {
                    double cu = 1;
                    double sum = 0;
                    if(u == 0)
                        cu = 1.0/sqrt(2.0);
                    for (int x = 0; x < 8; x++)
                    {
                        sum += cos((double)((2*x+1)*u*PI)/16.0) * f[xpos+x][ypos+v];
                    }

                    Fr[u][v] = 0.5 * cu * sum;
                }
            }
            for (int u = 0; u < 8; u++)
            {
                for (int v = 0; v < 8; v++)
                {
                    double cv = 1;
                    if(v == 0)
                        cv = 1.0/sqrt(2.0);
                    double sum = 0;
                    for (int y = 0; y < 8; y++)
                    {
                        sum += cos((double)((2*y+1)*v*PI)/16.0) * Fr[u][y];
                    }
                    coeffArray[xpos+u][ypos+v] = 0.5 * cv * sum;
                }
            }
            
			//! quantize the frequency coefficient of this block (basic requirement)
			//TODO
            for (int u = 0; u < 8; u++)
            {
                for (int v = 0; v < 8; v++)
                {
                    coeffArray[xpos+u][ypos+v] = round((double)coeffArray[xpos+u][ypos+v] / (double)QuantizationMatrix[u][v]);
                }
            }
			
		}
	}
	
	//! output the computed coefficient array
	FILE *fp = fopen("coeffs.txt", "w");
	for (int r = 0; r < rows; r++)
	{
		for (int c = 0; c < cols; c++)
		{
			fprintf(fp, "%3.3lf ", coeffArray[c][r]);
		}
		fprintf(fp, "\n");
	}
	cout << "Quantized coefficients saved!" << endl;
    
	if (need_idct)
	{
		Bitmap reconstructedImg(cols, rows);
		//! apply IDCT on the quantized coefficients (enhancement part)
		for (int i = 0; i < blockRow; i++)
		{
			for (int j = 0; j < blockCol; j++)
			{
				int xpos = j*8, ypos = i*8;
                double Fr_[8][8] = {0};
				//! apply de-quantization on block_ij 
				//TODO
                for (int u = 0; u < 8; u++)
                {
                    for (int v = 0; v < 8; v++)
                    {
                        coeffArray[xpos+u][ypos+v] = (double)coeffArray[xpos+u][ypos+v] * QuantizationMatrix[u][v]);
                    }
                }

				//! apply IDCT on this block
				//TODO
                
                for (int x = 0; x < 8; x++)
                {
                    for (int y = 0; y < 8; y++)
                    {
                        double sum = 0;
                        for (int v = 0; v < 8; v++)
                        {
                            double cv = 1;
                            if(v == 0)
                                cv = 1.0/sqrt(2.0);
                            sum += cv * cos((double)((2*y+1)*v*PI)/16.0) * coeffArray[xpos+x][ypos+v];
                        }
                        Fr_[x][y] = 0.5 * sum;
                    }
                }
                
                for (int y = 0; y < 8; y++)
                {
                    for (int x = 0; x < 8; x++)
                    {
                        double sum = 0;
                        for (int u = 0; u < 8; u++)
                        {
                            double cu = 1;
                            if(u == 0)
                                cu = 1.0/sqrt(2.0);
                            sum += cu * cos((double)((2*x+1)*u*PI)/16.0) * Fr_[u][y];
                        }
                        f[xpos+x][ypos+y] = (int)0.5 * sum;
                    }
                }
				//! shiftting back the pixel value range to 0~255
				//TODO
                for (int y = 0; y < 8; y++)
                {
                    for (int x = 0; x < 8; x++)
                    {
                        f[xpos+x][ypos+y] += 128;
                        if(f[xpos+x][ypos+y]<0)
                            f[xpos+x][ypos+y] = 0;
                        if(f[xpos+x][ypos+y]>255)
                            f[xpos+x][ypos+y]=255;
                        reconstructedImg.setPixel(xpos+x, ypos+y, (unsigned char)f[xpos+x][ypos+y]);
                    }
                }
				
			}
		}
        
        
        
		string savePath = "reconstructedImg.bmp";
		reconstructedImg.save(savePath.c_str());
		cout << "reconstructed image saved!" << endl;
	}
    
	return 0;
}
