#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI  3.1415926535897932384626433832795
#define PI2 6.283185307179586476925286766559

int performWFR(unsigned char *inputFilename, unsigned char *outputFilename, int image_width, int image_height, double sigmax, double wxl,double wxi,double wxh, double sigmay,double wyl, double wyi, double wyh, double thr)
{

    printf("\nEntering C program...\n");

	//read the fringes image which is BYTE format 
    printf("Reading %s \n", inputFilename);
	unsigned char *fringes = (unsigned char *)calloc(image_width * image_height, sizeof(unsigned char));
	FILE *input_file = fopen(inputFilename, "rb");
    if(!input_file)
	{
		printf("\ncan't open file %s\n", inputFilename);
		return 1;
    }
	fread(fringes, sizeof(unsigned char), image_width * image_height, input_file);
	fclose(input_file);
	
	//convert fringes to float data type
	int i;
	float *fringes_float = (float *)malloc(image_width * image_height * sizeof(float));
	for (i=0; i<image_width*image_height; ++i)
	{
		fringes_float[i] = fringes[i];
	}







	//writing the results to the file
    printf("Writing %s \n", outputFilename);
	FILE *output_file = fopen(outputFilename, "wb");
  	fwrite(fringes_float, sizeof(float), image_width * image_height, output_file);
	fclose(output_file);


	//free used memory
	free(fringes);
	free(fringes_float);
	//free(wrapped_phase);
    printf("Leaving C program...\n \n");
	return 0;
}
