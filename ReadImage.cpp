/* ReadImage.cpp:
 * 	Contains functionality for inputting the pixel values contained
 * 	within PGM and PPM images into a class of type ImageType.
 * Credit:
 * 	Dr. George Bebis
 * 		The source code was obtained via his website at
 * 		https://www.cse.unr.edu/~bebis/CS479/ .
 * 		It has been slightly modified (renaming functions,
 * 		adding comments, etc) for our intents and purposes.
 * Date:
 * 	15 March 2023
 */


// Libraries
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

using namespace std;

#include "image.h"
#include "rgb.h"


// Functions

/* readImagePGM:
 * 	Inputs the pixel values contained within a PGM image into
 * 	a given storage location.
 * args:
 * 	@fname: Path to file to read pixel values from.
 * 	@image: The location to store the pixel values to.
 * return:
 * 	void
 */
void readImagePGM(char fname[], ImageType& image)
{
 int i, j;
 int N, M, Q;
 unsigned char *charImage;
 char header [100], *ptr;
 ifstream ifp;

 ifp.open(fname, ios::in | ios::binary);

 if (!ifp) {
   cout << "Can't read image: " << fname << endl;
   exit(1);
 }

 // read header

 ifp.getline(header,100,'\n');
 if ( (header[0]!=80) ||    /* 'P' */
      (header[1]!=53) ) {   /* '5' */
      cout << "Image " << fname << " is not PGM" << endl;
      exit(1);
 }

 string head = header;
 //head = head.substr(2);
 strcpy(header, &header[2]);
 while(header[0]=='#')
   ifp.getline(header,100,'\n');

 M=strtol(header,&ptr,0);
 N=atoi(ptr);
 Q=strtol(header,&ptr,0);
 
 ImageType newImage = ImageType(N, M, Q);
 image = newImage;
 
 charImage = (unsigned char *) new unsigned char [M*N];

 ifp.read( reinterpret_cast<char *>(charImage), (M*N)*sizeof(unsigned char));

 if (ifp.fail()) {
   cout << "Image " << fname << " has wrong size" << endl;
   exit(1);
 }

 ifp.close();

 //
 // Convert the unsigned characters to integers
 //

 int val;
 for(i=0; i<N; i++)
   for(j=0; j<M; j++) {
     val = (int)charImage[i*M+j];
     image.setPixelVal(i, j, val);     
   }

 delete [] charImage;

}


/* readImagePPM():
 * 	Inputs the pixel values contained within a PPM image into
 * 	a given storage location.
 * args:
 * 	@fname: Path to file to read pixel values from.
 * 	@image: The location to store the pixel values to.
 * return:
 * 	void
 */
void readImagePPM(char fname[], ImageType& image)
{
 int i, j;
 int N, M, Q;
 unsigned char *charImage;
 char header [100], *ptr;
 ifstream ifp;

 ifp.open(fname, ios::in | ios::binary);
  
 if (!ifp) {
   cout << "Can't read image: " << fname << endl;
   exit(1);
 }

 // read header

 ifp.getline(header,100,'\n');

 if ( (header[0]!=80) ||    /* 'P' */
      (header[1]!=54) ) {   /* '6' */
      cout << "Image " << fname << " is not PPM" << endl;
      exit(1);
 }

 ifp.getline(header,100,'\n');
 while(header[0]=='#')
   ifp.getline(header,100,'\n');

 M=strtol(header,&ptr,0);
 N=atoi(ptr);
 
 ifp.getline(header,100,'\n');
 Q=strtol(header,&ptr,0);

 charImage = (unsigned char *) new unsigned char [3*M*N];

 ifp.read( reinterpret_cast<char *>(charImage), (3*M*N)*sizeof(unsigned char));

 if (ifp.fail()) {
   cout << "Image " << fname << " has wrong size" << endl;
   exit(1);
 }

 ifp.close();
 
 /* Convert the unsigned characters to integers */

 RGB val;

 for(i=0; i < N; i++)
  for(j=0; j < 3*M; j+=3) {
    val.r = (int)charImage[i*3*M+j];
    val.g = (int)charImage[i*3*M+j+1];
    val.b = (int)charImage[i*3*M+j+2];
    image.setPixelVal(i, j/3, val);
  }

delete [] charImage;

}
