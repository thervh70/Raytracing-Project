#ifndef IMAGE_JFDJKDFSLJFDFKSDFDJFDFJSDKSFJSDLF
#define IMAGE_JFDJKDFSLJFDFKSDFDJFDFJSDKSFJSDLF
#ifdef WIN32
#include <windows.h>
#endif
#include "GL/glut.h"
#include "settings.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fstream>
# define M_PI           3.14159265358979323846  /* pi */
# define M_E            2.71828182845904523536  /* e  */

extern float focusDepth;


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
	Image(void) {
		setSize(WindowSize_X, WindowSize_Y);
		int rowsize = floor((24 * _width + 31) / 32) * 4;
		imageC.resize(rowsize * _height);
	}
	Image(int width, int height)
	: _width(width)
	, _height(height)
	{
		_image.resize(3 * _width*_height);
		_depth.resize(_width*_height);

		int rowsize = floor((24 * _width + 31) / 32) * 4;
		imageC.resize(rowsize * _height);
	}
	void setPixel(int i, int j, const RGBValue & rgb, float depth)
	{
		_image[3*(_width*j+i)]=rgb[0];
		_image[3*(_width*j+i)+1]=rgb[1];
		_image[3*(_width*j+i)+2]=rgb[2];
		_depth[_width*j+i] = depth;
	}
	void setSize(int w, int h) {
		_width = w;
		_height = h;
		_image.resize(3 * w*h);
		_depth.resize(w*h);
	}

	void getDepthRange(float &lowest, float &highest) {
		lowest = std::numeric_limits<float>::max();
		highest = 0;

		for (int i = 0; i < _depth.size(); i++) {
			lowest = std::min(lowest, _depth[i]);
			if(_depth[i] != std::numeric_limits<float>::max())
				highest = std::max(highest, _depth[i]);
		}
	}

	void normailzeDepthRange() {
		float lowest, highest, denom;
		getDepthRange(lowest, highest);
		denom = std::max(std::abs(lowest - focusDepth), std::abs(highest - focusDepth));

		for (int i = 0; i < _depth.size(); i++) {
			if (_depth[i] != std::numeric_limits<float>::max())
				_depth[i] = std::abs(_depth[i] - focusDepth) / denom;
			else
				_depth[i] = 1;
		}
	}

	std::vector<std::vector<float>> gauseanMap(float intesisty) {
		float total = 0.;
		std::vector<std::vector<float>> map;
		map.resize(13);

		if (intesisty >= 1) {
			for (int i = 0; i < 13; i++) {
				std::vector<float> row;
				row.resize(13);

				for (int j = 0; j < 13; j++) {
					row[j] = 0.0f;
				}
				map[i] = row;
			}
			map[6][6] = 1;
			return map;
		}

		for (int i = 0; i < 13; i++) {
			std::vector<float> row;
			row.resize(13);

			for (int j = 0; j < 13; j++) {
				row[j] = 
					1 / (intesisty * depthPower * 2 * M_PI)*
					std::powf(
						M_E, 
						( -(i - 6)*(i - 6) -(j - 6)*(j - 6) ) / 2 * std::powf(intesisty * depthPower, 2)
					);
				total += row[j];
			}
			map[i] = row;
		}

		for (int i = 0; i < 13; i++) {
			for (int j = 0; j < 13; j++) {
				map[i][j] /= total;
			}
		}

		return map;
	}

	void printGauseanMap(float intens) {
		std::vector<std::vector<float>> map = gauseanMap(intens);

		std::cout << "MAP FOR SIGMA = " << intens << std::endl;
		for (int i = 0; i < 13; i++) {

			std::cout << "[ ";
			for (int j = 0; j < 13; j++) {
				std::cout << map[i][j] << ", \t";
			}
			std::cout << " ]" << std::endl;
		}
		std::cout << "_______________________" << std::endl << std::endl;
	}

	int blurrWidth() {
		return _width - 12;
	}

	int blurrHeight() {
		return _height - 12;
	}

	void initBlurrMap() {
		_blurr.resize(3* blurrWidth() * blurrHeight());
	}

	void blurrLine(int y) {
		if (y < 6 || y >= (_height-6)) return;
		std::vector<std::vector<float>> map;
		float r, g, b;
		for (int i = 6; i < (_width-6); i++) {

			r = g = b = 0.0f;
			map = gauseanMap(_depth[y*_width + i]);

			for (int j = 0; j < 13; j++) {
				for (int k = 0; k < 13; k++) {
					r += _image[3 * ((y + j - 6)*_width + i + k - 6)] * map[j][k];
					g += _image[3 * ((y + j - 6)*_width + i + k - 6)+1] * map[j][k];
					b += _image[3 * ((y + j - 6)*_width + i + k - 6)+2] * map[j][k];
				}
			}

			_blurr[3 * ((y - 6)*blurrWidth() + i - 6) ] = r;
			_blurr[3 * ((y - 6)*blurrWidth() + i - 6) + 1] = g;
			_blurr[3 * ((y - 6)*blurrWidth() + i - 6) + 2] = b;
		}
	}

	void blurrAllLines() {
		for (int i = 6; i < (_height - 6); i++) {
			Image::blurrLine(i); 
		}
	}

	std::vector<unsigned char> imageC;
	std::vector<float> _image;
	std::vector<float> _depth;
	std::vector<float> _blurr;
	int _width;
	int _height;

	bool writeImagePPM(const char * filename);
	bool writeImageBMP(const char * filename);
	bool writeDepthBMP(const char * filename);
	bool writeBlurredBMP(const char * filename);
	bool writeImageBMP(const char * filename, int x, int y, int w, int h);
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
	bmpinfoheader[8] = (unsigned char)(-_height);			// Minus height, because standard BMP renders bottom-to-top
	bmpinfoheader[9] = (unsigned char)((-_height) >> 8);
	bmpinfoheader[10] = (unsigned char)((-_height) >> 16);
	bmpinfoheader[11] = (unsigned char)((-_height) >> 24);

	fwrite(bmpfileheader, 1, 14, file);
	fwrite(bmpinfoheader, 1, 40, file);

	// Print raw pixel data
	int rowsize = floor((24 * _width + 31) / 32) * 4;

	for (unsigned int y = 0; y < _height; ++y) {
		unsigned int i = 0;
		for (; i < _width * 3; i += 3) {
			// RGB are flipped to BGR
			imageC[y*rowsize + i] = (unsigned char)(_image[y*_width*3 + i + 2] * 255.0f);
			imageC[y*rowsize + i + 1] = (unsigned char)(_image[y*_width*3 + i + 1] * 255.0f);
			imageC[y*rowsize + i + 2] = (unsigned char)(_image[y*_width*3 + i] * 255.0f);
		}
		for (; i < rowsize; ++i) {
			imageC[y*rowsize + i] = (unsigned char)1;
		}
	}

	int t = fwrite(&(imageC[0]), _width * _height * 3, 1, file);
	if (t != 1)
	{
		printf("Dump file problem... fwrite\n");
		return false;
	}

	fclose(file);
	return true;
}

bool Image::writeDepthBMP(const char * filename)
{
	FILE* file;
	file = fopen(filename, "wb");
	if (!file)
	{
		printf("dump file problem... file\n");
		return false;
	}

	Image::normailzeDepthRange();

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
	bmpinfoheader[8] = (unsigned char)(-_height);			// Minus height, because standard BMP renders bottom-to-top
	bmpinfoheader[9] = (unsigned char)((-_height) >> 8);
	bmpinfoheader[10] = (unsigned char)((-_height) >> 16);
	bmpinfoheader[11] = (unsigned char)((-_height) >> 24);

	fwrite(bmpfileheader, 1, 14, file);
	fwrite(bmpinfoheader, 1, 40, file);

	// Print raw pixel data
	int rowsize = floor((24 * _width + 31) / 32) * 4;

	for (unsigned int y = 0; y < _height; ++y) {
		unsigned int i = 0;
		for (; i < (_width * 3); ++i) {
			if (y*_width + i < _depth.size()) {
				// RGB are flipped to BGR
				imageC[y*rowsize + i * 3] =
					imageC[y*rowsize + i * 3 + 1] =
					imageC[y*rowsize + i * 3 + 2] =
					(unsigned char)((1 - _depth[y*_width + i]) * 255.0f);
			}
		}
		for (; i < rowsize; ++i) {
			imageC[y*rowsize + i] = (unsigned char)1;
		}
	}

	int t = fwrite(&(imageC[0]), _width * _height * 3, 1, file);
	if (t != 1)
	{
		printf("Dump file problem... fwrite\n");
		return false;
	}

	fclose(file);
	return true;
}

bool Image::writeBlurredBMP(const char * filename)
{
	FILE* file;
	file = fopen(filename, "wb");
	if (!file)
	{
		printf("dump file problem... file\n");
		return false;
	}

	Image::normailzeDepthRange();
	Image::initBlurrMap();
	Image::blurrAllLines();

	// Print BMP header
	unsigned char bmpfileheader[14] = { 'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0 };
	unsigned char bmpinfoheader[40] = { 40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0 };
	unsigned char bmppad[3] = { 0,0,0 };


	int s = blurrWidth() * blurrHeight();
	bmpfileheader[2] = (unsigned char)(s);
	bmpfileheader[3] = (unsigned char)(s >> 8);
	bmpfileheader[4] = (unsigned char)(s >> 16);
	bmpfileheader[5] = (unsigned char)(s >> 24);

	bmpinfoheader[4] = (unsigned char)(blurrWidth());
	bmpinfoheader[5] = (unsigned char)(blurrWidth() >> 8);
	bmpinfoheader[6] = (unsigned char)(blurrWidth() >> 16);
	bmpinfoheader[7] = (unsigned char)(blurrWidth() >> 24);

	bmpinfoheader[8] = (unsigned char)(-blurrHeight());			// Minus height, because standard BMP renders bottom-to-top
	bmpinfoheader[9] = (unsigned char)(-blurrHeight() >> 8);
	bmpinfoheader[10] = (unsigned char)(-blurrHeight() >> 16);
	bmpinfoheader[11] = (unsigned char)(-blurrHeight() >> 24);

	fwrite(bmpfileheader, 1, 14, file);
	fwrite(bmpinfoheader, 1, 40, file);

	// Print raw pixel data
	int rowsize = floor((24 * blurrWidth() + 31) / 32) * 4;

	for (unsigned int y = 0; y < blurrHeight(); ++y) {
		unsigned int i = 0;
		for (; i < blurrWidth() * 3; i += 3) {
			// RGB are flipped to BGR
			imageC[y*rowsize + i ]     = (unsigned char)(_blurr[y*blurrWidth() * 3 + i + 2] * 255.0f);
			imageC[y*rowsize + i + 1 ] = (unsigned char)(_blurr[y*blurrWidth() * 3 + i + 1] * 255.0f);
			imageC[y*rowsize + i + 2 ] = (unsigned char)(_blurr[y*blurrWidth() * 3 + i] * 255.0f);
		}
		for (; i < rowsize; ++i) {
			imageC[y*rowsize + i - 3] = (unsigned char)1;
		}
	}

	int t = fwrite(&(imageC[0]), s * 3, 1, file);
	if (t != 1)
	{
		printf("Dump file problem... fwrite\n");
		return false;
	}

	fclose(file);
	return true;
}

bool Image::writeImageBMP(const char* filename, int x, int y, int w, int h) {
	std::fstream s(filename, std::fstream::in | std::fstream::out | std::fstream::binary);
	for (int Y = y; Y < y + h; ++Y) {
		std::vector<char> imageC(w * 3);
		for (int j = x; j < x + w; ++j) {
			// RGB are flipped to BGR
			imageC[(j - x) * 3] = (char)(_image[Y*_width * 3 + j * 3 + 2] * 255.0f);
			imageC[(j - x) * 3 + 1] = (char)(_image[Y*_width * 3 + j * 3 + 1] * 255.0f);
			imageC[(j - x) * 3 + 2] = (char)(_image[Y*_width * 3 + j * 3] * 255.0f);
		}
		const std::vector<char> imageConst(imageC);
		s.seekp(54 + (Y*w + x) * 3, std::ios_base::beg);
		s.write(&(imageConst[0]), w * 3);
	}
	s.close();
	return true;
}

#endif