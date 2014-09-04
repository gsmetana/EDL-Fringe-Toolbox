/*
This program is written by Dr. Abdulbasit Abid and Dr. Munther Gdeisat
The Genaral Engineering Research Institute (GERI www.ljmu.ac.uk/GERI), Liverpool John Moores University,
James Parsons Building, Byrom Stree, Liverpool L3 3AF, United Kindgodm.
To contact the authors please use the email m.a.gdeisat@ljmu.ac.uk

This program implements the wavelet transform profilometry (WTP) algorithm in C language.

Modified slightly by Gregory Smetana (gsmetana@caltech.edu) to improve ease
of use through a flexible interface with Python

Written on 1st Sept 2008, revised June 2014
Version 1.1

For more information about how to use the software, please see the companion user guide
*/

//#include "stdafx.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI  3.1415926535897932384626433832795
#define PI2 6.283185307179586476925286766559
#define YES 1
#define NO  0
#define COST 1
#define MAX  2
#define MORLET 1
#define PAUL4  2
#define TINY 1E-30

/*** Start of Functions to extend the borders of the image using the LPC method (1)******/
// This code is taken from numerical recipes in C, second edition 1992
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define NR_END 1
#define FREE_ARG char*

float *vector(int nl, int nh)
// allocate a float vector with subscript range v[nl..nh] 
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	return v-nl+NR_END;
}

void free_vector(float *v, int nl, int nh)
// free a float vector allocated with vector() 
{
	free((FREE_ARG) (v+nl-NR_END));
}


void memcof(float data[], int n, int m, float *xms, float d[])
{
	int k,j,i;
	float p=0.0,*wk1,*wk2,*wkm;

	wk1=vector(1,n);
	wk2=vector(1,n);
	wkm=vector(1,m);
	for (j=1;j<=n;j++) p += SQR(data[j]);
	*xms=p/n;
	wk1[1]=data[1];
	wk2[n-1]=data[n];
	for (j=2;j<=n-1;j++) {
		wk1[j]=data[j];
		wk2[j-1]=data[j];
	}
	for (k=1;k<=m;k++) {
		float num=0.0,denom=0.0;
		for (j=1;j<=(n-k);j++) {
			num += wk1[j]*wk2[j];
			denom += SQR(wk1[j])+SQR(wk2[j]);
		}
		d[k]=2.0*num/denom;
		*xms *= (1.0-SQR(d[k]));
		for (i=1;i<=(k-1);i++)
			d[i]=wkm[i]-d[k]*wkm[k-i];
		if (k == m) {
			free_vector(wkm,1,m);
			free_vector(wk2,1,n);
			free_vector(wk1,1,n);
			return;
		}
		for (i=1;i<=k;i++) wkm[i]=d[i];
		for (j=1;j<=(n-k-1);j++) {
			wk1[j] -= wkm[k]*wk2[j];
			wk2[j]=wk2[j+1]-wkm[k]*wk1[j+1];
		}
	}
}

void predic(float data[], int ndata, float d[], int m, float future[], int nfut)
{
    int k,j;
    float sum,discrp,*reg;
    reg=vector(1,m);
    for (j=1;j<=m;j++)
        reg[j]=data[ndata+1-j];
    for (j=1;j<=nfut;j++) {
        discrp=0.0;
        sum=discrp;
        for (k=1;k<=m;k++)
            sum += d[k]*reg[k];
        for (k=m;k>=2;k--)
            reg[k]=reg[k-1];
        future[j]=reg[1]=sum;
    }
    free_vector(reg,1,m);
}


void flipLRrow(float *row, int n)
{
	int i;
	float *temp_row = (float *) calloc(n, sizeof(float) );
	for (i=0; i<n; ++i)
		temp_row[i] = row[n-i-1];

	for (i=0; i<n; ++i)
		row[i] = temp_row[i];

	free(temp_row);
}


//Extend the image using linear predictive algorithm. For more details please see
//chapter 13.6 of numerical recipes in C, second edition 1992
//http://www.nr.com/oldverswitcher.html
float *extend_image_linear_prediction(float *fringes, int *image_width_ptr, int *image_height_ptr, int extend_by_pixels)
{
	int image_width  = *image_width_ptr;
	int image_width_new  = *image_width_ptr + 2 * extend_by_pixels;
	int image_height = *image_height_ptr;
	int i, j;
	float *extended_fringes = (float *) calloc( image_width_new * image_height, sizeof(float) );

	//add zero padding to the right and left of the image (zero padding)
	for (i=0; i < image_height; ++i)
	{
		for (j=extend_by_pixels; j < image_width_new - extend_by_pixels; ++j)
		{
			extended_fringes[i*image_width_new + j] = fringes[i*image_width + j - extend_by_pixels];			
		}
	}

	//predict the right border of the image
	int no_of_coeffs = 100;
	float error = 0; 
	float *coeffs = vector(1,no_of_coeffs);

	for (i=0; i < image_height; ++i)
	{
		memcof(extended_fringes + i*image_width_new + extend_by_pixels - 1, image_width, no_of_coeffs, &error, coeffs);
		predic(extended_fringes + i*image_width_new + extend_by_pixels - 1, image_width, coeffs, no_of_coeffs, 
			extended_fringes + i*image_width_new + extend_by_pixels + image_width - 1, extend_by_pixels);
	}

	//predict the left border of the image
	for (i=0; i < image_height; ++i)
	{
		flipLRrow(extended_fringes + i*image_width_new, image_width_new);
		memcof(extended_fringes + i*image_width_new + extend_by_pixels - 1, image_width, no_of_coeffs, &error, coeffs);
		predic(extended_fringes + i*image_width_new + extend_by_pixels - 1, image_width, coeffs, no_of_coeffs, 
			extended_fringes + i*image_width_new + extend_by_pixels + image_width - 1, extend_by_pixels);
		flipLRrow(extended_fringes + i*image_width_new, image_width_new);

	}

	
	//Update the image width and free allocated memory
	*image_width_ptr = image_width_new;
	free(fringes);
	free_vector(coeffs, 1, no_of_coeffs+1);

	return extended_fringes;
}

float *crop_wrapped_phase(float *wrapped_phase, int *image_width_ptr, int *image_height_ptr, int extend_by_pixels)
{
	int image_width  = *image_width_ptr;
	int image_width_new  = *image_width_ptr - 2 * extend_by_pixels;
	int image_height = *image_height_ptr;
	int crop_by_pixels = extend_by_pixels;
	int i, j;
	float *cropped_wrapped_phase = (float *) calloc( image_width_new * image_height, sizeof(float) );

	//crop the left and right of the image borders
	for (i=0; i < image_height; ++i)
	{
		for (j=crop_by_pixels; j < image_width - crop_by_pixels; ++j)
		{
			cropped_wrapped_phase[i*image_width_new + j - crop_by_pixels] = wrapped_phase[i*image_width + j];			
		}
	}
	
	//Update the image width
	*image_width_ptr = image_width_new;
	free(wrapped_phase);

	return cropped_wrapped_phase;

}
/************end of extending the borders of the image using LPC method (1)**********/

/******* Start of Functions to calculate the ridge (2)********/
struct NODE
{
	double *local_maxima;
	int *indices_local_maxima_scale;
	double *cost;
	int *previous_link_node;
	int no_of_local_maxima;
};

typedef struct NODE NODE;

int ridge_maximum(float *wrapped_phase, double* Wavelet_Result_R, double* Wavelet_Result_C, double *Wavelet_magnitude, int image_width, int no_of_scales)
{
	int scale_matrix_size_int = image_width * no_of_scales;
	double Max;
	int i, j, index1, index2, X_Pos_Max, Y_Pos_Max;

	for (i=0; i<scale_matrix_size_int; ++i)
	{
		Wavelet_magnitude[i] = sqrt(Wavelet_Result_R[i] * Wavelet_Result_R[i] + Wavelet_Result_C[i] * Wavelet_Result_C[i]);
	}

	double temp;
	//find the maximum magnitude in each column in the matrix
	for (i=0; i<image_width; ++i)
	{
		temp = -999999999.;
		for (j=0; j<no_of_scales; ++j)
		{
			index1 = j * image_width + i;
			if(Wavelet_magnitude[index1] > temp)
			{
				temp = Wavelet_magnitude[index1];
				Y_Pos_Max = j;
			}
		}
		index2 = Y_Pos_Max * image_width + i;
		wrapped_phase[i] = atan2(Wavelet_Result_C[index2], Wavelet_Result_R[index2]);
	}
	return 0;
}

int local_maxima(double *Wavelet_magnitude, int image_width, int no_of_scales, NODE *nodes, int node_index, double *diffrentiated_input)
{
	//diffrentiate the incoming signal
	int i, counter = 0;
	double diff;
	int index1, index2;
	for (i=0; i<no_of_scales - 1; ++i)
	{
		index1 = i * image_width;
		index2 = (i + 1) * image_width;
		diff = Wavelet_magnitude[index2] - Wavelet_magnitude[index1];
		if(  (Wavelet_magnitude[index1] != 0) && (Wavelet_magnitude[index2] != 0) && (diff != 0) )
		{
			diffrentiated_input[counter] = abs(diff);
			++counter;
		}
	}
	
	//finding delta
	double delta = 99999999999999.0;
	for (i=0; i<counter; ++i)
		if (diffrentiated_input[i] < delta)	
			delta = diffrentiated_input[i];
	
	//finding the local maxima
	double Max = -999999999999., Min = 999999999999999.9;
	int MaxPos;
	int LookForMax = 1;
	double thisValue;
	int No_of_local_maxima = 0;
	for (i=0; i<no_of_scales; ++i)
	{
		index1 = i * image_width;
		thisValue = Wavelet_magnitude[index1];
		if (thisValue > Max)
		{
			Max = thisValue;
			MaxPos = index1;
		}

		if (thisValue < Min)
			Min = thisValue;

		if (LookForMax)
		{
			if (thisValue <= (Max - delta) )
			{
				nodes[node_index].local_maxima[No_of_local_maxima] = Max;
				nodes[node_index].indices_local_maxima_scale[No_of_local_maxima] = MaxPos;
				No_of_local_maxima++;
				Min = thisValue;
				LookForMax = 0;
			}
		}
		else
		{
			if (thisValue > (Min + delta) )
			{
				Max  = thisValue;
				MaxPos = index1;
				LookForMax = 1;
			}
		}
	}

	nodes[node_index].no_of_local_maxima = No_of_local_maxima;

	//if the local_maxima function fails to find a maximum then consider the global maximum as local maximum
	Max = -999999999999999.9;
	if (No_of_local_maxima == 0)
	{
		for (i=0; i<no_of_scales; ++i)
		{
			index1 = i * image_width;
			if (Wavelet_magnitude[index1] > Max)
			{
				Max = Wavelet_magnitude[index1];
				MaxPos = index1;
			}
		}
		nodes[node_index].local_maxima[No_of_local_maxima] = Max;
		nodes[node_index].indices_local_maxima_scale[No_of_local_maxima] = MaxPos;
		No_of_local_maxima++;
		nodes[node_index].no_of_local_maxima = No_of_local_maxima;
	}

	return 0; 
}

//calculate the phase of a row of a fringe pattern using the Cost Ridge Extraction Algorithm
int ridge_cost(float *wrapped_phase, double* Wavelet_Result_R, double* Wavelet_Result_C, double *Wavelet_magnitude, int image_width, int no_of_scales, NODE *nodes)
{
	int scale_matrix_size_int = image_width * no_of_scales;
	int i, j, current_node, previous_node;;
	double *diffrentiated_input = (double*) calloc(no_of_scales, sizeof(double));
	double *sub_cost = (double*) calloc(no_of_scales, sizeof(double));
	double Min = 999999999999999999.9;
	int MinPos_node, MinPos_scale, previous_link;

	for (i=0; i<scale_matrix_size_int; ++i)
	{
		Wavelet_magnitude[i] = sqrt(Wavelet_Result_R[i] * Wavelet_Result_R[i] + Wavelet_Result_C[i] * Wavelet_Result_C[i]);
	}

	//find the local minima and cost for all the row
	for (i=0; i<image_width; ++i)
	{
		local_maxima(Wavelet_magnitude + i, image_width, no_of_scales, nodes, i, diffrentiated_input);
		//find the cost function for all the nodes
		if (i > 0)
		{
			for (current_node = 0; current_node < nodes[i].no_of_local_maxima; ++current_node)
			{
				for (previous_node = 0; previous_node < nodes[i-1].no_of_local_maxima; ++previous_node)
				{
					sub_cost[previous_node] = nodes[i-1].cost[previous_node] +
						pow( (double) ((nodes[i].indices_local_maxima_scale[current_node] / image_width) - 
						(nodes[i-1].indices_local_maxima_scale[previous_node] / image_width)), 2.0) -
						pow((double)(nodes[i].local_maxima[current_node]), 2.0);
						
					//find the minimum of the subcost vector and assign it to the cost of the node
					if (sub_cost[previous_node] < Min)
					{
						Min = sub_cost[previous_node];
						MinPos_node = previous_node;

					}
				}
				nodes[i].cost[current_node] = Min;
				nodes[i].previous_link_node[current_node] = MinPos_node;
				Min = 999999999999999999.9;
			}
		}
		else
		{
			for (current_node=0; current_node < nodes[i].no_of_local_maxima; ++current_node)
			{
				nodes[i].cost[current_node] = 0;
				nodes[i].previous_link_node[current_node] = 0;
			}
		}//if
	} //for i
	
//find the minimum path through all the nodes (backward)
	int *optimal_path_indexes = (int *) calloc(image_width, sizeof(double));
	//find the minimum cost for the last node
	Min = 999999999999999999.9;
	for (j=0; j<nodes[image_width-1].no_of_local_maxima; ++j)
	{
		if (nodes[image_width-1].cost[j] < Min)
		{
			Min = nodes[image_width-1].cost[j];
			MinPos_scale = nodes[image_width-1].indices_local_maxima_scale[j];
			MinPos_node = j;
		}
	}
	optimal_path_indexes[image_width-1] = MinPos_scale;
	MinPos_node = nodes[image_width-1].previous_link_node[MinPos_node];

	//find the path for the rest of the nodes (right to left)
	for(i=image_width - 2; i > -1 ; --i)
	{
		optimal_path_indexes[i] = nodes[i].indices_local_maxima_scale[MinPos_node];
		MinPos_node = nodes[i].previous_link_node[MinPos_node];
	}

	//find the wrapped phase map for the row
	for (i=0; i<image_width; ++i)
	{
		wrapped_phase[i] = atan2(Wavelet_Result_C[ optimal_path_indexes[i] + i], Wavelet_Result_R[ optimal_path_indexes[i] + i]);
	}

	free(diffrentiated_input);
	free(optimal_path_indexes);
	free(sub_cost);
	return 0;
}


/******* End of Functions to calculate the ridge (2)********/

/**Start of Functions to calculate the WTP in the frequency domain (3)*******/

void remove_DC(float *fringes, int image_width, int image_height)
{
	float mean;
	int i, j;
	for(i=0; i<image_height; ++i)
	{
		//find the mean of each row
		mean = 0;
		for(j=0; j<image_width; ++j)
		{
			mean = mean + fringes[i*image_width + j];
		}

		//remove the DC from each row
		mean = mean/image_width;
		for(j=0; j<image_width; ++j)
		{
			fringes[i*image_width + j] = fringes[i*image_width + j] - mean;
		}
	}
}
//Computing Morlet wavelet in the Frequency domain
int MorletC_FT(double *W_Re, int image_width, double a /*scale*/, double sigma)
{
	int i, n = image_width;
	double *w = W_Re;
	double sqrt_a_n = sqrt(a)/n;
	double nDiv2 = n/2;
	int dx;
	for(dx=0; dx<n; ++dx)
	{
		if(dx <= nDiv2)
			w[dx] = PI2 * dx / n;
		else
			w[dx] = -PI2 * (n-dx) / n;
	}

	double Fb = -0.5/(sigma * sigma);
	double Fc = 1;//theoritiacally it should be set between 0.8 and 0.955 

	for(i=0; i<n; ++i)
	{
		w[i] = w[i]*a - PI2 * Fc;
		w[i] = sqrt_a_n * exp( (double) (Fb * w[i]*w[i]) );
	}

	return 1;
}

// Paul 4. Real part (Frequency domain) 
#define H(w) ( (w>0.0) ? 1.0 : 0.0 ) // Heaviside step function 
int PAUL4_FT(double *W_Re, int image_width, double a /*scale*/)
{
    int i, n = image_width;
	double *w = W_Re;
	double nDiv2 = n/2;
	int dx;
	for(dx=0; dx<n; ++dx)
	{
		if(dx <= nDiv2)
			w[dx] = PI2 * dx / n;
		else
			w[dx] = -PI2 * (n-dx) / n;
	}

    double c = 0.282465006843625;
    double w4;
	for(i=0; i<n; ++i)
	{
		w[i] = w[i]*a;
		w4 = pow(w[i], 4.0);
		w[i] = c * w4 * H(w[i]) * exp( (double) -w[i]);
	}

	return 1;
}

//isign = 1  forward Fourier Transform
//isign = -1 inverse Fourier Transform
void fft(double *re, double *im, unsigned int n, int isign)
{
	unsigned int i, j, k, l, le, le1, ip, n2;
	double wpr, wpi, wr, wi, wtr, wti;
	n2 = n/2;
	j = 1;

	for(i=0; i<n-1; i++) 
	{
		if(i<j) 
		{
			wtr = re[j-1];
			wti = im[j-1];
			re[j-1] = re[i];
			im[j-1] = im[i];
			re[i]   = wtr;
			im[i]   = wti;
		}

		k = n2;
		while(k<j)
		{
			j -= k;
			k >>= 1;
		}
		j += k;
	}//for i

	l=1;
	k=n;
	while(k>>=1)
	{
		le1 = (le=1<<l++) >> 1;
		wtr = PI / (double)le1;
		wpr = cos(wtr); 
		wpi = -isign*sin(wtr);
		wr = 1.0;       
		wi = 0.0;

		for(j=0; j<le1; j++) 
		{
			for(i=j; i<n; i+=le) 
			{
				ip = i + le1;
				wtr    = wr*re[ip] - wi*im[ip];
				wti    = wi*re[ip] + wr*im[ip];
				re[ip] = re[i] - wtr;
				im[ip] = im[i] - wti;
				re[i]  = re[i] + wtr;
				im[i]  = im[i] + wti;
			}

			wr = (wtr=wr)*wpr - wi*wpi;
			wi = wi*wpr + wtr*wpi;
		}
	}//while

	//normalization
	double factor = 1/sqrt( (double) n);
	for(j=0; j<n; j++) 
	{
		re[j] = re[j]*factor;
		im[j] = im[j]*factor;
	}	
}

//isign = 1  forward Fourier Transform
//isign = -1 inverse Fourier Transform

void dft(double *Re,double *Im, int n, int isign, double * Re_temp, double *Im_temp, double *SineLookupTable, double *CosineLookupTable)
{
	int i, k;
	double arg;
	double cosarg, sinarg;
	
	for (i=0; i<n; i++) 
	{
		Re_temp[i] = 0;
		Im_temp[i] = 0;
		for (k=0; k<n; k++) 
		{
			cosarg = CosineLookupTable[i * n + k];
			sinarg = isign * SineLookupTable[i * n + k];
			Re_temp[i] += (Re[k] * cosarg - Im[k] * sinarg);
			Im_temp[i] += (Re[k] * sinarg + Im[k] * cosarg);
		}
	}

	/* Copy the data back */
	for (i=0; i<n; i++) 
	{
		Re[i] = Re_temp[i];
		Im[i] = Im_temp[i];
	}

	//normalization
	double factor = 1/sqrt( (double) n);
	int j;
	for(j=0; j<n; j++) 
	{
		Re[j] = Re[j]*factor;
		Im[j] = Im[j]*factor;
	}	

}

int WFA(float *fringes, float *wrapped_phase, int image_width, int image_height, double starting_scale, double scale_step, double ending_scale, int ridge_alg, double sigma, int extend_by_pixels, int wavelet_type)
{
	int image_size = image_width * image_height;
	int image_width_FFT = 1;
	int i, j, k, order = 1;

	//determine to use FFT or DFT for calculations
	while (image_width_FFT < image_width) 
	{
		image_width_FFT = pow(2.0,(double) order);
		++order;
	}
	order = order - 1;
	int no_of_scales = ((ending_scale - starting_scale)/scale_step + 1);

	//Allocate memory for nodes to be used by the cost function
	struct NODE *nodes = (NODE *) malloc(image_width * sizeof(NODE));
	for(i=0; i<image_width; ++i)
	{
		nodes[i].local_maxima = (double *) calloc(no_of_scales, sizeof(double)); 
		nodes[i].indices_local_maxima_scale = (int *) calloc(no_of_scales, sizeof(int));  
		nodes[i].cost = (double *) calloc(no_of_scales, sizeof(double)); 
		nodes[i].previous_link_node = (int *) calloc(no_of_scales, sizeof(int));  
	}

	//R: real   C: imaginary
	double *fringes_R = (double *) calloc(image_size, sizeof(double));
	double *fringes_C = (double *) calloc(image_size, sizeof(double));	

	//covert one row of the fringe pattern to double
	for(i=0; i<image_size; ++i)
	{
		fringes_R[i] = fringes[i];
	}

	double *Re_temp;
	double *Im_temp;
	double *SineLookupTable ;
	double *CosineLookupTable;

	if(image_width_FFT != image_width)
	{
		//house keeping for the DFT function
		Re_temp = (double*) calloc(image_width, sizeof(double) );
		Im_temp = (double*) calloc(image_width, sizeof(double) );

		//construct a lookup table for the DFT function
		double arg;
		SineLookupTable   = (double*) calloc(image_width * image_width, sizeof(double) );
		CosineLookupTable = (double*) calloc(image_width * image_width, sizeof(double) );
		for (i=0; i<image_width; i++) 
		{
			arg = -PI2 * (double)i / (double)image_width;
			for (k=0; k<image_width; k++) 
			{
				SineLookupTable[i*image_width + k]   = sin(k * arg);
				CosineLookupTable[i*image_width + k] = cos(k * arg);
			}
		}
	}
	
	//calculate the complex FFT or DFT for the fringe pattern row by row
	int increment = 0;
	for(i=0; i<image_height; ++i)
	{
		increment = image_width * i;
		if(image_width_FFT == image_width)
			fft(fringes_R + increment, fringes_C + increment, image_width, 1);
		else
			dft(fringes_R + increment, fringes_C + increment, image_width, 1, Re_temp, Im_temp, SineLookupTable, CosineLookupTable);
	}
	
	//calculate the Morlet coefficients wavelet 
	//--Note that this Morlet wavelet is in the frequency domain
	//constuct a matrix with the size image_width X no_of_scales
	//wavelet coefficients for complex 1D Morlet in the frequency domain
	double *Wavelet_R = (double*) calloc(image_width * no_of_scales, sizeof(double) );
	double *Wavelet_C = (double*) calloc(image_width * no_of_scales, sizeof(double) );
	double scale = starting_scale;
	int dy;
	for (dy=0; dy<no_of_scales; ++dy)
	{
		if (wavelet_type == MORLET)
			MorletC_FT(Wavelet_R + dy * image_width, image_width, scale, sigma);
		else if (wavelet_type == PAUL4)
			PAUL4_FT(Wavelet_R + dy * image_width, image_width, scale);

		scale = scale + scale_step;
	}

	double *Wavelet_Result_R = (double*) calloc(image_width * no_of_scales, sizeof(double) );
	double *Wavelet_Result_C = (double*) calloc(image_width * no_of_scales, sizeof(double) );
	double *Wavelet_magnitude = (double*) calloc(image_width * no_of_scales, sizeof(double) );
	int incrementF, incrementW;
	//process the image row by row using 1D-CWT (Complex Morlet)
	for (i=0; i<image_height; ++i)
	{
		for (dy=0; dy<no_of_scales; ++dy)
		{
			incrementW = dy * image_width;
			incrementF = i * image_width;
			for (j=0; j<image_width; ++j)
			{
				*(Wavelet_Result_R + incrementW + j) = *(Wavelet_R + incrementW + j) * *(fringes_R + incrementF + j);
				*(Wavelet_Result_C + incrementW + j) = *(Wavelet_R + incrementW + j) * *(fringes_C + incrementF + j);
			}
			if(image_width_FFT == image_width)
				fft(Wavelet_Result_R + incrementW, Wavelet_Result_C + incrementW, image_width, -1);
			else
				dft(Wavelet_Result_R + incrementW, Wavelet_Result_C + incrementW, image_width, -1, Re_temp, Im_temp, SineLookupTable, CosineLookupTable);
		}
		//extract wrapped phase using maximum ridge algorithm for the row
		switch (ridge_alg)
		{	
			case MAX: 
				ridge_maximum(wrapped_phase + i * image_width, Wavelet_Result_R, Wavelet_Result_C, Wavelet_magnitude, image_width, no_of_scales);
				break;

			case COST: 
				ridge_cost(wrapped_phase + i * image_width, Wavelet_Result_R, Wavelet_Result_C, Wavelet_magnitude, image_width, no_of_scales, nodes);
				break;

			default:	printf("please specify one of the ridge extraction algorithms\n");
						printf("Input 1 for Maximum ridge extraction algorithm\n");
						printf("Input 2 for Cost Function ridge extraction algorithm\n");
				break;
		}
	}
 
	//Set the value of the wrapped phase to invalid value
	int index_new, index_old;
	int image_width_old = image_width - 2 * extend_by_pixels;
	
	if(image_width_FFT != image_width)
	{
		free(SineLookupTable);
		free(CosineLookupTable);
		free(Re_temp);
		free(Im_temp);
	}
	free(Wavelet_magnitude);
	free(Wavelet_Result_R);
	free(Wavelet_Result_C);
	free(Wavelet_R);
	free(Wavelet_C);
	for(i=0; i<image_width; ++i)
	{
		free(nodes[i].local_maxima);
		free(nodes[i].indices_local_maxima_scale);
		free(nodes[i].cost);
		free(nodes[i].previous_link_node);
	}
	free(nodes);
	free(fringes_R);
	free(fringes_C);

	return 1;
} 
/**End of Functions to calculate the WTP in the frequency domain (3)*****/

/**Start of Functions to calculate the WTP in the time domain (4) *****/
// Complex Morlet
#define FB  2.0 /* bandwidth parameter */
#define FC  1 /* wavelet center frequency and it should be theoritically between 0.8 and 0.95 but 1 is fine practically*/ 

// Complex Morlet. Real part (Time domain) 
double CMORLETreal(double x, double a, double b, double sigma)
{
      double Fb = FB; /* bandwidth parameter */
      double Fc = FC; /* wavelet center frequency */
      double c = 1.0 / sqrt(PI*Fb);
      double x2;
	  double result;
      if ( a == 0.0 ) a = TINY;	  
	  x = (x - b) / a;
      x2 = x * x * sigma * sigma;
      result = c * exp(-x2/Fb) * cos(PI2*Fc*x);
	  return result;
}

// Complex Morlet. Imaginary part (Time domain) 
double CMORLETimag(double x, double a, double b, double sigma)
{
      double Fb = FB; // bandwidth parameter 
      double Fc = FC; // wavelet center frequency 
      double c = 1.0 / sqrt(PI*Fb);
      double x2;
      if ( a == 0.0 ) a = TINY;	  
	  x = (x - b) / a;
      x2 = x * x * sigma * sigma;
      double result = c * exp(-x2/Fb) * sin(PI2*Fc*x);
	  return result;
}

// Paul 4. Real part (Time domain)
double CPAUL4real(double x, double a, double b)
{
	  const double c = 1.078936850151577;
      double x2, x4, x6, x8, x10;
      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      x6 = x4 * x2;
      x8 = x4 * x4;
      x10 = x8 * x2;
      return c * (5.0*x4 - 10.0*x2 + 1.0) /
             (x10 + 5.0*x8 + 10.0*x6 + 10.0*x4 + 5.0*x2 + 1.0);
}

// Paul 4. Imaginary part (Time domain) 
double PAUL4imag(double x, double a, double b)
{
      const double c = 1.078936850151577;
      double x2, x4, x6, x8, x10;

      if ( a == 0.0 ) a = TINY;
      x = (x - b) / a;
      x2 = x * x;
      x4 = x2 * x2;
      x6 = x4 * x2;
      x8 = x4 * x4;
      x10 = x8 * x2;
      return c * x * (x4 - 10.0*x2 + 5.0) /
             (x10 + 5.0*x8 + 10.0*x6 + 10.0*x4 + 5.0*x2 + 1.0);
}

//non optimised continuous wavelet transform calculations in the time domain
int CWT(double *Wavelet_R_Result, double *Wavelet_I_Result, double *fringes_R, int image_width, 
		double starting_scale, double scale_step, double ending_scale, int n_points, double *W_re, double *W_im)
{
	//set some coeffcients required to perform the wavelet calculations
	double ivalp = 2;
	double bstep = 1;
	double d = (double)n_points / 8 ;
	double istep = 1.0 / ivalp;
    double ivlp_amin = ivalp * starting_scale;
    int dx, dy, i; 
	int no_of_scales = ((ending_scale - starting_scale)/scale_step + 1);

	double a_esl = -4;
	double a_esr =  4;
	double a, b, t1, t2, ind, f, factor;
	int R = n_points/2;
	int index;

	//The wavelet matrix has the size image_width X no_of_scales
	for (dy = 0, a = starting_scale; dy < no_of_scales; dy++, a+=scale_step)
    {				
		// calculate ivalp for next a 
		if(dy && ivalp != 1) 
		{
			ivalp = (unsigned int)ceil( ivlp_amin / a );
			istep = 1.0 / (double)ivalp;				
		}

		a_esl = -4 * a;
		a_esr = -a_esl;

        // Translations 
        for (dx = 0, b = 0.0; dx < image_width; dx++, b+=bstep)
        {
			// Compute wavelet boundaries 
			t1 = a_esl + b;    
			t2 = a_esr + b;
			if(t1<0.0)				
				t1=0.0; 
			if(t2>=image_width)		
				t2=(image_width-1);

			// Perform convolution 
			index = dy*image_width + dx;
			Wavelet_R_Result[index] = 0.0;
			Wavelet_I_Result[index] = 0.0;
			for (f = t1; f <= t2; f+=istep) 
			{
				ind = d*(f-b)/a + R;
				Wavelet_R_Result[index] += fringes_R[(unsigned int)f] * W_re[(int)ind];
				Wavelet_I_Result[index] += fringes_R[(unsigned int)f] * W_im[(int)ind];
			}
			factor = sqrt(a);
			Wavelet_R_Result[index] *= 1/(factor * (double)ivalp);
			Wavelet_I_Result[index] *= 1/(factor * (double)ivalp);
		}							
	}  
	return 0;
}

//very optimised continuous wavelet transform calculations in the time domain
int CWTO3(double *Wavelet_R_Result, double *Wavelet_I_Result, double *fringes_R, int image_width, 
		double starting_scale, double scale_step, double ending_scale, int n_points, double *W_re, double *W_im)
{
	//set some coeffcients required to perform the wavelet calculations
	double ivalp = 2;
	double bstep = 1;
	double d = (double)n_points / 8 ;
	double istep = 1.0 / ivalp;
    double ivlp_amin = ivalp * starting_scale;
    int dx, dy, i; 
	int no_of_scales = ((ending_scale - starting_scale)/scale_step + 1);

	double a_esl = -4;
	double a_esr =  4;
	double a, b, t1, t2, ind, f, factor;
	int R = n_points/2;
	int index;

	//The wavelet matrix has the size image_width X no_of_scales
	for (dy = 0, a = starting_scale; dy < no_of_scales; dy++, a+=scale_step)
    {				
		// calculate ivalp for next a 
		if(dy && ivalp != 1) 
		{
			ivalp = (unsigned int)ceil( ivlp_amin / a );
			istep = 1.0 / (double)ivalp;				
		}

		a_esl = -4 * a;
		a_esr = -a_esl;

        //Translations 
        for (dx = 0, b = 0.0; dx < image_width; dx++, b+=bstep)
        {
			// Compute wavelet boundaries 
			t1 = a_esl + b;    
			t2 = a_esr + b;
			if(t1<0.0)				
				t1=0.0; 
			if(t2>=image_width)		
				t2=(image_width-1);

			// Perform convolution 
			index = dy*image_width + dx;
			Wavelet_R_Result[index] = 0.0;
			Wavelet_I_Result[index] = 0.0;
			for (f = t1; f <= t2; f+=istep) 
			{
				ind = d*(f-b)/a + R;
				Wavelet_R_Result[index] += fringes_R[(unsigned int)f] * W_re[(int)ind];
				Wavelet_I_Result[index] += fringes_R[(unsigned int)f] * W_im[(int)ind];
			}
			factor = sqrt(a);
			Wavelet_R_Result[index] *= 1/(factor * (double)ivalp);
			Wavelet_I_Result[index] *= 1/(factor * (double)ivalp);
		}							
	}  
	return 0;
}

//calculate the WTA in the time domain
int WTA(float *fringes, float *wrapped_phase, int image_width, int image_height, double starting_scale, 
		double scale_step, double ending_scale, int ridge_alg, double sigma, int extend_by_pixels, int wavelet_type)
{
	int image_size = image_width * image_height;
	int i, j;
	int no_of_scales = ((ending_scale - starting_scale)/scale_step + 1);

	//Allocate memory for nodes to be used by the cost function
	struct NODE *nodes = (NODE *) malloc(image_width * sizeof(NODE));
	for(i=0; i<image_width; ++i)
	{
		nodes[i].local_maxima = (double *) calloc(no_of_scales, sizeof(double)); 
		nodes[i].indices_local_maxima_scale = (int *) calloc(no_of_scales, sizeof(int));  
		nodes[i].cost = (double *) calloc(no_of_scales, sizeof(double)); 
		nodes[i].previous_link_node = (int *) calloc(no_of_scales, sizeof(int));  
	}
	

	//covert each row of the fringe pattern to double
	double *fringes_R = (double *) calloc(image_size, sizeof(double));
	for(i=0; i<image_size; ++i)
	{
		fringes_R[i] = fringes[i];
	}
	
	//Pre computation required to the Morlet wavelet coefficients  
	//Increase N_points to increase accurucy of the wavelet transform
	int N_points = 60000;
	double wstep = (double) 8 / N_points ;
	int L = -4 / wstep - 1;  
	int R =  4 / wstep + 1;
	double step;
	double *W_re = (double*) calloc(R-L, sizeof(double) );
	double *W_im = (double*) calloc(R-L, sizeof(double) );

	if (wavelet_type == MORLET)
	{
		for(j=0, step = -4; j<(R-L); j++, step+=wstep)
		{			
			W_re[j] = CMORLETreal(step, 1.0, 0.0, sigma);
			W_im[j] = CMORLETimag(step, 1.0, 0.0, sigma);
		}
	}

	if (wavelet_type == PAUL4)
	{
		for(j=0, step = -4; j<(R-L); j++, step+=wstep)
		{			
			W_re[j] = CPAUL4real(step, 1.0, 0.0);
			W_im[j] = PAUL4imag(step, 1.0, 0.0);
		}
	}

	//calculate the complex Morlet wavelet coefficients  
	//--Note that this Morlet wavelet is in the time domain
	//constuct a matrix with the size image_width X no_of_scales
	double *Wavelet_Result_R = (double*) calloc(image_width * no_of_scales, sizeof(double) );
	double *Wavelet_Result_I = (double*) calloc(image_width * no_of_scales, sizeof(double) );
	double *Wavelet_magnitude = (double*) calloc(image_width * no_of_scales, sizeof(double) );
	double scale = starting_scale;
	int increment;
	
	//process the image row by row using 1D-CWT (Complex Morlet)
	printf("\nPlease wait as this may take few minutes to finish\n\n");
	for (i=0; i<image_height; ++i)
	{
		increment = i * image_width;
		//use only one of the two functions
		CWTO3(Wavelet_Result_R, Wavelet_Result_I, fringes_R + increment, image_width, starting_scale, scale_step, ending_scale, N_points, W_re, W_im);
		//CWT(Wavelet_Result_R, Wavelet_Result_I, fringes_R + increment, image_width, starting_scale, scale_step, ending_scale, N_points, W_re, W_im);

		
		//extract wrapped phase using maximum ridge algorithm for the row
		switch (ridge_alg)
		{	
			case MAX: 
				ridge_maximum(wrapped_phase + increment, Wavelet_Result_R, Wavelet_Result_I, Wavelet_magnitude, image_width, no_of_scales);
			break;
				case COST: 
					ridge_cost(wrapped_phase + increment, Wavelet_Result_R, Wavelet_Result_I, Wavelet_magnitude, image_width, no_of_scales, nodes);
			break;
				default:	printf("please specify one of the ridge extraction algorithms\n");
							printf("Input 1 for Maximum ridge extraction algorithm");
							printf("Input 2 for Cost Function ridge extraction algorithm");
			break;
		}
	}

	//Set the value of the wrapped phase to invalid value
	int index_new, index_old;
	int image_width_old = image_width - 2 * extend_by_pixels;
		
	free(Wavelet_magnitude);
	free(Wavelet_Result_R);
	free(Wavelet_Result_I);
	free(W_re);
	free(W_im);
	
	for(i=0; i<image_width; ++i)
	{
		free(nodes[i].local_maxima);
		free(nodes[i].indices_local_maxima_scale);
		free(nodes[i].cost);
		free(nodes[i].previous_link_node);
	}
	free(nodes);
	free(fringes_R);

	return 1;
} 
/**End of Functions to calculate the WTP in the time domain (4) *****/


int performWTP(unsigned char *inputFilename, unsigned char *outputFilename, int image_width, int image_height, int to_extend_image, int extend_by_pixels, int use_FFT, int wavelet_type, int ridge_alg, double starting_scale, double scale_step, double ending_scale, double Morlet_sigma)
{
    printf("\nEntering C program...\n");
	if (to_extend_image == NO)  extend_by_pixels = 0;

	int image_width_FFT = 0, order = 0;
	if (use_FFT == YES && to_extend_image == YES)
	{
		//determine log2(image_width)
		while (image_width_FFT <= image_width) 
		{
			image_width_FFT = pow(2.0,(double) order);
			++order;
		}
		extend_by_pixels = (image_width_FFT - image_width)/2;
	}

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
	//	printf("\n debug %f\n", fringes_float[i]);
	}

	//subtract the image from its mathematical mean
	remove_DC(fringes_float, image_width, image_height);

	if (to_extend_image == YES)
	{
		printf("Extending the image borders using linear prediction,...\n");
		fringes_float = extend_image_linear_prediction(fringes_float, &image_width, &image_height, extend_by_pixels);
	
	}

	float *wrapped_phase = (float *)malloc(image_width * image_height * sizeof(float));

	//use either the WFA or the WTA function.
	if (use_FFT == YES)
	{
		printf("Performing wavelet transform in Frequency domain, please wait...\n");
		WFA(fringes_float, wrapped_phase, image_width, image_height, starting_scale, scale_step, ending_scale, ridge_alg, Morlet_sigma, extend_by_pixels, wavelet_type);
	}
	else
	{
		printf("Performing wavelet transform in time domain, please wait...\n");
		WTA(fringes_float, wrapped_phase, image_width, image_height, starting_scale, scale_step, ending_scale, ridge_alg, Morlet_sigma, extend_by_pixels, wavelet_type);
	}

	//crop the image to return it to the original image dimensions
	if (to_extend_image == YES)
		wrapped_phase = crop_wrapped_phase(wrapped_phase, &image_width, &image_height, extend_by_pixels);

	//writing the results to the file
    printf("Writing %s \n", outputFilename);
	FILE *output_file = fopen(outputFilename, "wb");
  	fwrite(wrapped_phase, sizeof(float), image_width * image_height, output_file);
	fclose(output_file);

	//free used memory
	free(fringes);
	free(fringes_float);
	free(wrapped_phase);
    printf("Leaving C program...\n \n");
	return 0;
}
