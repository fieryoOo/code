#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
//#include "utils.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <string>

cv::Mat imageRGBA;
cv::Mat imageOut;

uchar4 *d_rgbaImage__;
uchar4 *d_blurImage__;

size_t numRows() { return imageRGBA.rows; }
size_t numCols() { return imageRGBA.cols; }

//return types are void since any internal error will be handled by quitting
//no point in returning error codes...
//returns a pointer to an RGBA version of the input image
//and a pointer to the single channel grey-scale output
//on both the host and device
void preProcess(uchar4 **inputImage, uchar4 **blurImage,
                uchar4 **d_rgbaImage, uchar4 **d_blurImage,
                const std::string &filename) {
  //make sure the context initializes ok
  cudaFree(0);//checkCudaErrors(cudaFree(0));

  cv::Mat image;
  image = cv::imread(filename.c_str(), CV_LOAD_IMAGE_COLOR);
  if (image.empty()) {
    std::cerr << "Couldn't open file: " << filename << std::endl;
    exit(1);
  }

  cv::cvtColor(image, imageRGBA, CV_BGR2RGBA);

  //allocate memory for the output
  imageOut.create(image.rows, image.cols, CV_32FC4); //CV_8UC1

  //This shouldn't ever happen given the way the images are created
  //at least based upon my limited understanding of OpenCV, but better to check
  if (!imageRGBA.isContinuous() || !imageOut.isContinuous()) {
    std::cerr << "Images aren't continuous!! Exiting." << std::endl;
    exit(1);
  }

  *inputImage = (uchar4 *)imageRGBA.ptr<unsigned char>(0);
  *blurImage  = (uchar4 *)imageOut.ptr<unsigned char>(0);

  const size_t numPixels = numRows() * numCols();
  //allocate memory on the device for both input and output
  cudaMalloc(d_rgbaImage, sizeof(uchar4) * numPixels); // checkCudaErrors
  cudaMalloc(d_blurImage, sizeof(uchar4) * numPixels); // checkCudaErrors
  cudaMemset(*d_blurImage, 0, numPixels * sizeof(uchar4)); //make sure no memory is left laying around; checkCudaErrors

  //copy input array to the GPU
  cudaMemcpy(*d_rgbaImage, *inputImage, sizeof(uchar4) * numPixels, cudaMemcpyHostToDevice); // checkCudaErrors

  d_rgbaImage__ = *d_rgbaImage;
  d_blurImage__ = *d_blurImage;
}

void postProcess(const std::string& output_file) {
  const int numPixels = numRows() * numCols();
  //copy the output back to the host
  cudaMemcpy(imageOut.ptr<uchar4>(0), d_blurImage__, sizeof(uchar4) * numPixels, cudaMemcpyDeviceToHost); // checkCudaErrors

  //output the image
  cv::Mat imageOut2;
  imageOut2.create(imageOut.rows, imageOut.cols, CV_8UC4);
  cv::cvtColor(imageOut, imageOut2, CV_RGBA2BGR);
  cv::imwrite(output_file.c_str(), imageOut2);

  //cleanup
  cudaFree(d_rgbaImage__);
  cudaFree(d_blurImage__);
}
