
#include "utils.h"
#include <vector>

float* loadToGPU(float* data, size_t items);

void getBestSize(int imgsToProcess, int origXSize, int origYSize, int &batchSize, int &xSize, int &ySize);
size_t getFreeMem(int device);

void kernel1(float* imgs, size_t oldX, size_t oldY, int noOfImages, size_t newX, size_t newY,
		float* filter,
		std::complex<float>*& result);
//		float*& result);

void kernel3(float maxShift, size_t noOfImgs,
		const std::complex<float>* imgs, size_t fftXdim, size_t fftYdim,
		float*& result, std::complex<float>*& result2);

