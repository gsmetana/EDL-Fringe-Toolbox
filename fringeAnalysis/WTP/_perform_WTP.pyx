cdef extern int performWTP(unsigned char *inputFilename,
                            unsigned char *outputFilename,
                            int image_width, 
                            int image_height, 
                            int to_extend_image, 
                            int extend_by_pixels, 
                            int use_FFT, 
                            int wavelet_type, 
                            int ridge_alg, 
                            double starting_scale, 
                            double scale_step, 
                            double ending_scale, 
                            double Morlet_sigma)

def perform_WTP(unsigned char *inputFilename,
                            unsigned char *outputFilename,
                            int image_width, 
                            int image_height, 
                            int to_extend_image, 
                            int extend_by_pixels, 
                            int use_FFT, 
                            int wavelet_type, 
                            int ridge_alg, 
                            double starting_scale, 
                            double scale_step, 
                            double ending_scale, 
                            double Morlet_sigma):
    performWTP(inputFilename,
                            outputFilename,
                            image_width, 
                            image_height, 
                            to_extend_image, 
                            extend_by_pixels, 
                            use_FFT, 
                            wavelet_type, 
                            ridge_alg, 
                            starting_scale, 
                            scale_step, 
                            ending_scale, 
                            Morlet_sigma)
