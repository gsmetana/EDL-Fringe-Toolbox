#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include <fftw3.h>

#define PI  3.1415926535897932384626433832795
#define PI2 6.283185307179586476925286766559

//*****************************************************************************/
int performWFR(unsigned char *inputFilename, unsigned char *outputFilename, int image_width, int image_height, double sigmax, double wxl,double wxi,double wxh, double sigmay,double wyl, double wyi, double wyh, double thr)
//*****************************************************************************/
{

    printf("\nEntering C program...\n");

    int m, n;
    unsigned char *fringes;
    FILE *input_file;

	int i, j;
	float *f;

    int sx, sy;
    int mm, nn;
    double *w0;
    double complex *w; 

    fftw_complex *fexp;   
    fftw_complex *Ff;

    fftw_complex *wexp;   
    fftw_complex *Fw;   

    fftw_complex *FfFw;   
    fftw_complex *sf;   

    fftw_plan plan_forward;
    fftw_plan plan_backward;

    int step, steps; 
    double wyt, wxt;

    double *gwx, *gwy, *gphase, *gr;

    m = image_height;
    n = image_width;

	//read the fringes image which is BYTE format 
    printf("Reading %s \n", inputFilename);
	fringes = (unsigned char *)calloc(n * m, sizeof(unsigned char));
    input_file = fopen(inputFilename, "rb");

    if(!input_file)
	{
		printf("\ncan't open file %s\n", inputFilename);
		return 1;
    }
	fread(fringes, sizeof(unsigned char), n * m, input_file);
	fclose(input_file);
	
	//convert fringes to float data type
	f = (float *)malloc(n * m * sizeof(float));
	for (i=0; i<n*m; ++i)
	{
		f[i] = fringes[i];
	}

    //half window size along x, by default 3*sigmax; window size: 2*sx+1
    sx=round(3*sigmax); 
    //half window size along y, by default 3*sigmay; window size: 2*sy+1
    sy=round(3*sigmay); 

    mm=m+2*sx;
    nn=n+2*sy;

    fexp = fftw_malloc ( sizeof ( fftw_complex ) * nn * mm );

    fexpand(f, m, n, mm, nn, fexp); 

    Ff = fftw_malloc ( sizeof ( fftw_complex ) * nn * mm );

    plan_forward = fftw_plan_dft_2d ( nn, mm, fexp, Ff, FFTW_FORWARD, 
        FFTW_ESTIMATE );

    fftw_execute ( plan_forward );

    w0 = (double *)malloc( (2*sx+1) * (2*sy+1) * sizeof(double));
    initwindow(w0, sx, sy, sigmax, sigmay);

    steps=floor((wyh-wyl)/wyi)+1;

    w = (complex double*)malloc( (2*sx+1) * (2*sy+1)  * sizeof(double complex ));
    wexp = fftw_malloc ( sizeof ( fftw_complex ) * nn * mm );
    Fw = fftw_malloc ( sizeof ( fftw_complex ) * nn * mm );

    FfFw = fftw_malloc ( sizeof ( fftw_complex ) * nn * mm );
    sf = fftw_malloc ( sizeof ( fftw_complex ) * nn * mm );        

    // allocate to zero
    gwx  = (double *)calloc( m * n,  sizeof(double));
    gwy = (double *)calloc( m * n, sizeof(double));
    gphase = (double *)calloc( m * n,  sizeof(double));
    gr = (double *)calloc( m * n,  sizeof(double));

    for(wyt=wyl; wyt <= wyh; wyt+=wyi)
    {
        step=floor((wyt-wyl)/wyi)+1;
        printf("%f % \n", 100 * (float) step/steps);

        for(wxt=wxl; wxt <= wxh; wxt+=wxi)
        {
            shiftwindow(w0, sx, sy, wxt, wyt, w);

            cfexpand(w, (2*sx+1), (2*sy+1), mm, nn, wexp); 

            plan_forward = fftw_plan_dft_2d ( nn, mm, wexp, Fw, FFTW_FORWARD, FFTW_ESTIMATE );
            fftw_execute ( plan_forward );

            for ( j = 0; j < mm; j++ )
            {
                for ( i = 0; i < nn; i++ )
                {
                    FfFw[i *mm  + j] = Ff[i *mm  + j] * Fw[i *mm  + j];
                }
            }
            plan_backward = fftw_plan_dft_2d ( nn, mm, FfFw, sf, FFTW_BACKWARD, FFTW_ESTIMATE );
            fftw_execute ( plan_backward );

            for ( j = 0; j < m; j++ )
            {
                for ( i = 0; i < n; i++ )
                {
                    // inverse data is scaled by NX*NY
                    sf[ (i+sy) *mm  + (j+sx) ] = sf[(i+sy) *mm  + (j+sx) ] / ( double complex) ( nn * mm);
                    
                    // update where necessary
                    if( cabs( sf[ (i+sy) *mm  + (j+sx) ] ) > gr[i*m + j] )
                    {
                        gr[i*m + j] = cabs( sf[ (i+sy) *mm  + (j+sx) ] );
                        gwx[i*m + j] = wxt;
                        gwy[i*m + j] = wyt;
                        gphase[i*m + j] = carg( sf[ (i+sy) *mm  + (j+sx) ] );
                    }
                }
            }
        }
    }


    // END WFR

	//writing the results to the file
    printf("Writing %s \n", outputFilename);
	FILE *output_file = fopen(outputFilename, "wb");
  	

    fwrite(gphase, sizeof(double), n * m, output_file); 

	fclose(output_file);

	//free used memory
	free(fringes);
	free(f);
    free(w0);

    free(w);

    fftw_free (fexp);
    fftw_free (Ff);
    fftw_free (FfFw);
    fftw_free (sf);

    free(gwx);
    free(gwy);
    free(gphase);
    free(gr);

    fftw_destroy_plan ( plan_forward );
    fftw_destroy_plan ( plan_backward );

	//free(wrapped_phase);
    printf("Leaving C program...\n \n");
	return 0;
}
//*****************************************************************************/
void fexpand(float f[], int m, int n, int mm, int nn, fftw_complex fexp[] )
//*****************************************************************************/
//expand f from [m n] to [mm nn]

{
    int i, j;
    for ( j = 0; j < mm; j++ )
    {
        for ( i = 0; i < nn; i++ )
        {
            if(j < m && i < n)
            {
                fexp[i*mm + j] = f[i*m +j];
            } else{
                fexp[i*mm + j] = 0;
            }
        }
    }
}

//*****************************************************************************/
void cfexpand(double complex f[], int m, int n, int mm, int nn, fftw_complex fexp[] )
//*****************************************************************************/
//expand f to [m n]

{
    int i, j;
    for ( j = 0; j < mm; j++ )
    {
        for ( i = 0; i < nn; i++ )
        {
            if(j < m && i < n)
            {
                fexp[i*mm + j] = f[i*m +j];
            } else{
                fexp[i*mm + j] = 0;
            }
        }
    }
}

//*****************************************************************************/
void initwindow(double w0[], int sx, int sy, double sigmax, double sigmay)
//*****************************************************************************/
// [y x]=meshgrid(-sy:sy,-sx:sx); 
// w0=exp(-x.*x/2/sigmax/sigmax-y.*y/2/sigmay/sigmay); 
// w0=w0/sqrt(sum(sum(w0.*w0))); 
{
    int i,j, x, y;
    double sum, norm;
    sum = 0;
    for ( j = 0; j < (2*sx +1); j++ )
    {
        for ( i = 0; i < (2*sy +1); i++ )
        {
            x = j - sx;
            y = i - sy;
            w0[i * (2*sx + 1) + j] = exp( -0.5*x*x/sigmax/sigmax - 0.5*y*y/sigmay/sigmay);
            sum += w0[i * (2*sx + 1) + j]* w0[i * (2*sx + 1) + j];
        }
    }
    norm = sqrt(sum);

    for ( j = 0; j < (2*sx +1); j++ )
    {
        for ( i = 0; i < (2*sy +1); i++ )
        {
            w0[i * (2*sx + 1) + j] = w0[i * (2*sx + 1) + j] / norm;
        }
    }
}

//*****************************************************************************/
void shiftwindow(double w0[], int sx, int sy, double wxt, double wyt, complex double w[])
//*****************************************************************************/
// w=w0.*exp(j*wxt*x+j*wyt*y);
{
    int i,j, x, y;
    for ( j = 0; j < (2*sx +1); j++ )
    {
        for ( i = 0; i < (2*sy +1); i++ )
        {
            x = j - sx;
            y = i - sy;
            w[i * (2*sx + 1) + j] = w0[i * (2*sx + 1) + j] * cexp(I * wxt*x + I*wyt*y);
        }
    }
}
