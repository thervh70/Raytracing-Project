#ifndef IMAGE_JFDJKDFSLJFDFKSDFDJFDFJSDKSFJSDLF
#define IMAGE_JFDJKDFSLJFDFKSDFDJFDFJSDKSFJSDLF
#ifdef WIN32
#include <windows.h>
#endif
#include "GL/glut.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>

//Image class 
//This class can be used to write your final result to an image. 
//You can open the image using a PPM viewer.

//YOU CAN IGNORE THIS CODE!
 class RGBValue
{
	public:
	RGBValue(float rI=0, float gI=0, float bI=0)
	: r(rI)
	, g(gI)
	, b(bI)
	{
		if (r>1)
			r=1.0;
		if (g>1)
			g=1.0;
		if (b>1)
			b=1.0;

		if (r<0)
			r=0.0;
		if (g<0)
			g=0.0;
		if (b<0)
			b=0.0;
	};
	
	float operator[](int i) const
	{
		switch(i)
		{
			case 0:
				return r;
			case 1:
				return g;
			case 2:
				return b;
			default: 
				return r;
		}
	}
	float & operator[](int i)
	{
		switch(i)
		{
			case 0:
				return r;
			case 1:
				return g;
			case 2:
				return b;
			default: 
				return r;
		}
	}
	float r, b,g;
};





class Image
{
	public:
	Image(int width, int height)
	: _width(width)
	, _height(height)
	{
		_image.resize(3*_width*_height);
	}
	void setPixel(int i, int j, const RGBValue & rgb)
	{
		_image[3*(_width*j+i)]=rgb[0];
		_image[3*(_width*j+i)+1]=rgb[1];
		_image[3*(_width*j+i)+2]=rgb[2];
		
	}
	std::vector<float> _image;
	int _width;
	int _height;

	bool writeImagePPM(const char * filename);
	bool writeImageBMP(const char * filename);
};

bool Image::writeImagePPM(const char * filename)
{
	FILE* file;
	file = fopen(filename, "wb");
	if (!file)
	{
		printf("dump file problem... file\n");
		return false;
	}

	// Print PPM header
	fprintf(file, "P6\n%i %i\n255\n", _width, _height);

	// Print raw pixel data
	std::vector<unsigned char> imageC(_image.size());

	for (unsigned int i = 0; i<_image.size();++i)
		imageC[i] = (unsigned char)(_image[i] * 255.0f);

	int t = fwrite(&(imageC[0]), _width * _height * 3, 1, file);
	if (t != 1)
	{
		printf("Dump file problem... fwrite\n");
		return false;
	}

	fclose(file);
	return true;
}

// WARNING: Width of the file must be dividable by 4!!!
bool Image::writeImageBMP(const char * filename)
{
	FILE* file;
	file = fopen(filename, "wb");
	if (!file)
	{
		printf("dump file problem... file\n");
		return false;
	}

	// Print BMP header
	unsigned char bmpfileheader[14] = { 'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0 };
	unsigned char bmpinfoheader[40] = { 40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0 };
	unsigned char bmppad[3] = { 0,0,0 };

	bmpfileheader[2] = (unsigned char)(_image.size());
	bmpfileheader[3] = (unsigned char)(_image.size() >> 8);
	bmpfileheader[4] = (unsigned char)(_image.size() >> 16);
	bmpfileheader[5] = (unsigned char)(_image.size() >> 24);

	bmpinfoheader[4] = (unsigned char)(_width);
	bmpinfoheader[5] = (unsigned char)(_width >> 8);
	bmpinfoheader[6] = (unsigned char)(_width >> 16);
	bmpinfoheader[7] = (unsigned char)(_width >> 24);
	bmpinfoheader[8] = (unsigned char)(_height);
	bmpinfoheader[9] = (unsigned char)(_height >> 8);
	bmpinfoheader[10] = (unsigned char)(_height >> 16);
	bmpinfoheader[11] = (unsigned char)(_height >> 24);

	fwrite(bmpfileheader, 1, 14, file);
	fwrite(bmpinfoheader, 1, 40, file);

	// Print raw pixel data
	std::vector<unsigned char> imageC(_image.size());

	for (unsigned int i = 0; i<_image.size();++i)
		imageC[i] = (unsigned char)(_image[i] * 255.0f);

	int t = fwrite(&(imageC[0]), _width * _height * 3, 1, file);
	if (t != 1)
	{
		printf("Dump file problem... fwrite\n");
		return false;
	}

	fclose(file);
	return true;
}

#endif