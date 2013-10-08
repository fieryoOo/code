// Homework 1
// Color to Greyscale Conversion

//A common way to represent color images is known as RGBA - the color
//is specified by how much Red, Grean and Blue is in it.
//The 'A' stands for Alpha and is used for transparency, it will be
//ignored in this homework.

//Each channel Red, Blue, Green and Alpha is represented by one byte.
//Since we are using one byte for each color there are 256 different
//possible values for each color.  This means we use 4 bytes per pixel.

//Greyscale images are represented by a single intensity value per pixel
//which is one byte in size.

//To convert an image from color to grayscale one simple method is to
//set the intensity to the average of the RGB channels.  But we will
//use a more sophisticated method that takes into account how the eye 
//perceives color and weights the channels unequally.

//The eye responds most strongly to green followed by red and then blue.
//The NTSC (National Television System Committee) recommends the following
//formula for color to greyscale conversion:

//I = .299f * R + .587f * G + .114f * B

//Notice the trailing f's on the numbers which indicate that they are 
//single precision floating point constants and not double precision
//constants.

//You should fill in the kernel as well as set the block and grid sizes
//so that the entire image is processed.

//#include "image_IO.cu"
//#include "utils.h"
#include <stdio.h>
#include <iostream>
#define blks 32;

__global__
void blur_image(const uchar4* const rgbaImage,
                       uchar4* greyImage,
                       int numRows, int numCols)
{
  //TODO
  //Fill in the kernel to convert from color to greyscale
  //the mapping from components of a uchar4 to RGBA is:
  // .x -> R ; .y -> G ; .z -> B ; .w -> A
  //
  //The output (greyImage) at each pixel should be the result of
  //applying the formula: output = .299f * R + .587f * G + .114f * B;
  //Note: We will be ignoring the alpha channel for this conversion

  //First create a mapping from the 2D block and grid locations
  //to an absolute 2D location in the image, then use that to
  //calculate a 1D offset
    int idx, idxc, npt=0;
    int ix, x = blockIdx.x * 32 + threadIdx.x;
    if( x >= numCols ) return;
    int iy, y = blockIdx.y * 16 + threadIdx.y;
    if( y >= numRows ) return;
    uchar4 rgbain;
    idxc = x + y * numCols;
    for(ix=x-0; ix<x+1; ix++) {
       if( ix<0 || ix>numCols ) continue;
       for(iy=y-0; iy<y+1; iy++) {
          if( iy<0 || iy>=numRows ) continue;
          idx = ix + iy * numCols;
          rgbain = rgbaImage[idx];
          greyImage[idxc].x += rgbain.x;
          greyImage[idxc].y += rgbain.y; 
          greyImage[idxc].z += rgbain.z;
	  npt++;
       }
    }
    //greyImage[idxc].x /= npt;
    //greyImage[idxc].y /= npt; 
    //greyImage[idxc].z /= npt; 
    //greyImage[idxc].w = rgbaImage[idxc].w;
    greyImage[idxc] = rgbaImage[idxc];
}

void your_rgba_to_greyscale(const uchar4 * const h_rgbaImage, uchar4 * const d_rgbaImage,
                            uchar4* d_greyImage, size_t numRows, size_t numCols)
{
  //You must fill in the correct sizes for the blockSize and gridSize
  //currently only one block with one thread is being launched
  const dim3 blockSize(32, 16, 1);  //TODO
  const dim3 gridSize( numCols/32+1, numRows/16+1, 1);  //TODO
  blur_image<<<gridSize, blockSize>>>(d_rgbaImage, d_greyImage, numRows, numCols);
 
  cudaDeviceSynchronize(); //checkCudaErrors(cudaGetLastError());
}
