#include "include\CImg.h";
#include "include\dirent.h";
#include "include\kiss_fft.h";
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>

const double PI = 3.141592653589793238462;

const int VARIANCE_INC_X = 1;
const int VARIANCE_INC_Y = 1;

//const int ITERATIONS = 6;
const int ITERATIONS = 2;
const float IMAGE_SCALE = 0.04166666666666666666666666666667 / 2;
const float SPECTRAL_SCALE = 0.125 / 2; //(400, 500)

using namespace cimg_library;

typedef CImg<float> Image;

std::vector<Image*>* loadImages(char* direc);
void display(Image* image);

Image* toLogPolar(Image* image);
std::vector<Image*>* loadLogPolar(std::vector<Image*>* images);

Image* toSpectral(Image* image);
std::vector<Image*>* loadSpectral(std::vector<Image*>* images);

Image* toSpectralFull(Image* image);
Image* toVarianceFull(Image* image, Image* spectralMean);
float* getOptimalCoordsFull(Image* image, Image* spectralMean);

Image* toVariance(Image* image, Image* spectralMean);
std::vector<Image*>* loadVariance(std::vector<Image*>* images, Image* spectralMean);

float getCovariance(Image* image1, Image* image2);
float* getOptimalCoords(Image* varianceImage);
void printOptimalCoords(Image* varianceImage);

int main()
{
	Image* meanImage = new Image("meanface.jpg");
	toSpectralFull(meanImage);

	std::vector<Image*>* images = loadImages("res/After/");
	
	for(int i=0; i<images->size(); i++)
	{
		Image* img = images->at(i);
		img->rotate(-90);
		img->channel(0);
		img->resize(img->width() * IMAGE_SCALE, img->height() * IMAGE_SCALE);
	}
	meanImage->resize(meanImage->width() * SPECTRAL_SCALE, meanImage->height() * SPECTRAL_SCALE);

	/*
	for(int i=0; i<images->size(); i++)
	{
		Image* img = images->at(i);
		img->rotate(-90);
		img->channel(0);

		Image* smallImage = new Image(img->get_resize(img->width() * IMAGE_SCALE, img->height() * IMAGE_SCALE));
		Image* smallMeanImage = new Image(meanImage->get_resize(meanImage->width() * SPECTRAL_SCALE, meanImage->height() * SPECTRAL_SCALE));

		float* smallCoords = getOptimalCoordsFull(smallImage, smallMeanImage);

		delete smallImage;
		delete smallMeanImage;

		for(int j = 1; j < ITERATIONS; j++)
		{
			const int BUFFER_REGION = 20;

			float prevScale = (int)pow((float)2, j-1);
			float scale = (int)pow((float)2, j);
			if(scale > 1 / IMAGE_SCALE)
			{
				scale = 1 / IMAGE_SCALE;
			}

			smallCoords[0] *= scale / prevScale;
			smallCoords[1] *= scale / prevScale;

			Image* smallImage = new Image(img->get_resize(img->width() * IMAGE_SCALE * scale, img->height() * IMAGE_SCALE * scale));
			Image* smallMeanImage = new Image(meanImage->get_resize(meanImage->width() * SPECTRAL_SCALE * scale, meanImage->height() * SPECTRAL_SCALE * scale));
			smallImage->crop(
				smallCoords[0] - smallMeanImage->width() / 2 - BUFFER_REGION,
				smallCoords[1] - smallMeanImage->height() / 2 + BUFFER_REGION,
				smallCoords[0] + smallMeanImage->width() / 2 - BUFFER_REGION,
				smallCoords[1] + smallMeanImage->height() / 2 + BUFFER_REGION
				);

			smallCoords = getOptimalCoordsFull(smallImage, smallMeanImage);

			std::stringstream ss;
			ss << j;
			smallImage->save(("smallImage" + ss.str() + ".bmp").c_str());
			smallImage->save(("smallMeanImage" + ss.str() + ".bmp").c_str());

			delete smallImage;
			delete smallMeanImage;
		}

		std::cout << smallCoords[0] << " " << smallCoords[1] << std::endl;

		delete smallCoords;

		Image* logPolar = toSpectralFull(img);
	}
	*/

	std::vector<Image*>* logPolarImages = loadLogPolar(images);
	std::vector<Image*>* spectralImages = loadSpectral(logPolarImages);
	std::vector<Image*>* varianceImages = loadVariance(spectralImages, meanImage);

	//CImgDisplay disp(images->at(0)->width(), images->at(0)->height(), "VARIANCE", 0);
	CImgDisplay disp(images->at(0)->width(), images->at(0)->height());

	CImgDisplay imgdisp(*images->at(0), "IMG", 0);
	CImgDisplay meandisp(*meanImage, "MEAN", 0);

	int currentImage = 0;
	bool isArrowLeft = false, isArrowRight = false;
	while(!disp.is_closed())
	{
		if(disp.is_keyARROWLEFT())
		{
			if(!isArrowLeft)
			{
				currentImage = (currentImage + images->size() - 1) % images->size();
				isArrowLeft = true;
			}
		}
		else
		{
			isArrowLeft = false;
		}
		if(disp.is_keyARROWRIGHT())
		{
			if(!isArrowRight)
			{
				currentImage = (currentImage + images->size() + 1) % images->size();
				isArrowRight = true;
			}
		}
		else
		{
			isArrowRight = false;
		}
		//disp.display(*varianceImages->at(currentImage));
		disp.display(*spectralImages->at(currentImage));
	}
	
	images->clear();
	logPolarImages->clear();
	spectralImages->clear();
	varianceImages->clear();
	delete meanImage;
}

std::vector<Image*>* loadImages(char* direc)
{
	std::vector<Image*>* result = new std::vector<Image*>();

	DIR* dir;
	dir = opendir(direc);

	struct dirent *ent;
	if(dir != NULL)
	{
		while((ent = readdir (dir)) != NULL)
		{
			std::string s = std::string(ent->d_name);
			size_t found = s.find(".JPG");
			if(found != std::string::npos)
			{
				Image* image = new Image((direc + std::string(ent->d_name)).c_str());
				result->push_back(image);
				printf(s.c_str());
				std::cout << std::endl;
			}
		}
	}
	else
	{
		printf("ERROR - Cannot load image");
	}
	closedir(dir);

	return result;
}

Image* toLogPolar(Image* image)
{
	float width = image->width();
	float height = image->height();
	Image* result = new Image(width, height);
	cimg_forXY(*result, _x, _y)
	{
		//x = radius (0 to smth)
		//y = degrees (0 to 2PI)
		float r = log((float)_x / 2);
		float a = _y / height * 2 * PI;

		float x = exp(r) * cos(a) + width / 2;
		float y = exp(r) * sin(a) + height / 2;

		(*result)(_x, _y) = (*image)(floor(x), floor(y));
	}
	return result;
}

std::vector<Image*>* loadLogPolar(std::vector<Image*>* images)
{
	std::vector<Image*>* result = new std::vector<Image*>();
	for(int i=0; i<images->size(); i++)
	{
			result->push_back(toLogPolar(images->at(i)));
	}
	return result;
}

Image* toSpectral(Image* image)
{
	float width = image->width();
	float height = image->height();
	Image* result = new Image(width, height);
	
	kiss_fft_cfg cfgF = kiss_fft_alloc(width, 0, NULL, NULL);
	kiss_fft_cfg cfgI = kiss_fft_alloc(width, 1, NULL, NULL);

	kiss_fft_cpx* cx_in;
	kiss_fft_cpx* cx_out;

	cx_in = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx) * width);
	cx_out = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx) * width);

	for(int _y = 0; _y < height; _y++)
	{
		for(int i = 0; i < width; i++)
		{
			cx_in[i].r = (*image)(i, _y);
			cx_in[i].i = 0;
		}

		kiss_fft(cfgF, cx_in, cx_out);

		for(int i = 0; i < width; i++)
		{
			cx_in[i].r = (pow(cx_out[i].r, 2) + pow(cx_out[i].i, 2)) / width;
			cx_in[i].i = 0;
		}

		kiss_fft(cfgI, cx_in, cx_out);

		float r0 = cx_out[0].r;
		r0 = 1;
		for(int i = 0; i < width; i++)
		{
			(*result)(i, _y) = (cx_out[i].r / r0 + 1) / 2 * 255;
		}
	}
	free(cfgF);
	free(cfgI);

	free(cx_in);
	free(cx_out);
	return result;
}

std::vector<Image*>* loadSpectral(std::vector<Image*>* images)
{
	std::vector<Image*>* result = new std::vector<Image*>();
	for(int i=0; i<images->size(); i++)
	{
		result->push_back(toSpectral(images->at(i)));
	}
	return result;
}

Image* toVariance(Image* image, Image* spectralMean)
{
	float width = image->width();
	float height = image->height();

	float swidth = spectralMean->width();
	float sheight = spectralMean->height();

	Image* result = new Image(width, height);
	result->fill(0);
	for(int x = 0; x < width - swidth; x += VARIANCE_INC_X)
	{
		for(int y = 0; y < height - sheight; y += VARIANCE_INC_Y)
		{
			Image* sectionImage = new Image(image->get_crop(x, y, x + swidth - 1, y + sheight - 1));
			Image* spectralImage = toSpectralFull(sectionImage);
			delete sectionImage;

			float covariance = getCovariance(spectralImage, spectralMean);
			delete spectralImage;

			(*result)(x + swidth / 2, y + sheight / 2) = covariance;
		}
	}
	//post process result
	float min = INT_MAX;
	float max = 0;
	for(int x = 0; x < width - swidth; x += VARIANCE_INC_X)
	{
		for(int y = 0; y < height - sheight; y += VARIANCE_INC_Y)
		{
			float covariance = (*result)(x + swidth / 2, y + sheight / 2);
			if(covariance < min)
			{
				min = covariance;
			}
			if(covariance > max)
			{
				max = covariance;
			}
		}
	}
	max -= min;
	for(int x = 0; x < width - swidth; x += VARIANCE_INC_X)
	{
		for(int y = 0; y < height - sheight; y += VARIANCE_INC_Y)
		{
			float value = (*result)(x + swidth / 2, y + sheight / 2);
			value = 255.0f - ((value - min) * 255.0f / max);
			(*result)(x + swidth / 2, y + sheight / 2) = value;
			if(value == 255.0f)
			{
				std::cout << x + swidth / 2 << " " << y + sheight / 2 << std::endl;
			}
		}
	}
	return result;
}
std::vector<Image*>* loadVariance(std::vector<Image*>* images, Image* spectralMean)
{
	std::vector<Image*>* result = new std::vector<Image*>();
	for(unsigned int i=0; i<images->size(); i++)
	{
		result->push_back(toVariance(images->at(i), spectralMean));
	}
	return result;
}

void buildMeanFace() //TEMP
{
	std::vector<Image*>* images = loadImages("res/After/");
	for(unsigned int i=0; i<images->size(); i++)
	{
		images->at(i)->channel(0);
		images->at(i)->crop(720, 168, 1674, 1290);
		images->at(i)->resize(60, 70);
	}
	Image* image = new Image(images->at(0)->width(), images->at(0)->height());
	image->fill(0);
	for(unsigned int i=0; i<images->size(); i++)
	{
		*image += *images->at(i);
	}
	*image /= images->size();
	image->save("meanface.bmp");
	images->clear();
}

Image* toSpectralFull(Image* image)
{
	Image* logPolar = toLogPolar(image);
	Image* spectral = toSpectral(logPolar);
	
	delete logPolar;

	return spectral;
}

float getCovariance(Image* image1, Image* image2)
{
	float variance = 0;
	for(int x = 0; x < image1->width(); x++)
	for(int y = 0; y < image1->height(); y++)
	{
		variance += pow((*image1)(x, y) - (*image2)(x, y), 2);
	}
	variance /= image1->width() * image1->height() - 1;

	return variance;
}

float* getOptimalCoords(Image* varianceImage)
{
	int width = varianceImage->width();
	int height = varianceImage->height();

	float* coords;
	coords = new float[2];
	coords[0] = 0;
	coords[1] = 0;

	float max = 0;
	for(int x = 0; x < width; x++)
	{
		for(int y = 0; y < height; y++)
		{
			float value = (*varianceImage)(x, y);
			if(value > max)
			{
				coords[0] = x;
				coords[1] = y;
				max = value;
			}
		}
	}

	//TODO implement meanshift
	const int MEANSHIFT_RADIUS = 10;
	const float GAUSSIAN_CONSTANT = -1.0f / (2 * pow((float)MEANSHIFT_RADIUS, 2));
	//const float GAUSSIAN_CONSTANT = 1;

	for(int i = 0; i < 5; i++) // 5 = meanshift iterations
	{
		float cx = 0;
		float cy = 0;

		float kernelsum = 0;
		for(int x = 0; x < width; x++)
		for(int y = 0; y < height; y++)
		{
			float norm = pow(x - coords[0], 2) + pow(y - coords[1], 2);
			if(sqrt(norm) <= MEANSHIFT_RADIUS)
			{
				//float kernel = exp(norm * GAUSSIAN_CONSTANT) * (*varianceImage)(x, y);
				float kernel = 1;
				cx += (x - coords[0]) * kernel;
				cy += (y - coords[1]) * kernel;
				kernelsum += kernel;
			}
		}
		cx /= kernelsum;
		cy /= kernelsum;

		coords[0] += cx;
		coords[1] += cy;
	}

	return coords;
}

void printOptimalCoords(Image* varianceImage)
{
	float* coords = getOptimalCoords(varianceImage);
	std::cout << coords[0] << " " << coords[1] << std::endl;
	delete coords;
}

Image* toVarianceFull(Image* image, Image* spectralMean)
{
	Image* spectralImage = toSpectralFull(image);
	Image* varianceImage = toVariance(spectralImage, spectralMean);
	delete spectralImage;
	return varianceImage;
}

float* getOptimalCoordsFull(Image* image, Image* spectralMean)
{
	Image* varianceImage = toVarianceFull(image, spectralMean);
	float* coords = getOptimalCoords(varianceImage);
	delete varianceImage;
	return coords;
}