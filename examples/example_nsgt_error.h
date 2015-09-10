#ifndef NSGT_EXAMPLE_NSGT_ERROR_H
#define NSGT_EXAMPLE_NSGT_ERROR_H

#include "example_utils.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>

void example_nsgt_error(double *input, size_t len, int samplerate) {
    printf("Test full nsgt\r\n");

    time_point_t init_time = timer_start();

    NsgtScaleDesc scale = {};
    scale.scale_type = NSGT_SCALE_OCT;
    scale.fmin = 50.;
    scale.fmax = 22050.;
    scale.num_bins = 24;
    Nsgt *nsgt = nsgt_create(&scale, samplerate, len, 0);
    NsgtResult *result = nsgt_alloc_result(nsgt);
    // or:
    // NsgtResult *result = (NsgtResult *) malloc(nsgt_result_get_size(nsgt));
    // nsgt_init_result(nsgt, result);

    //todo: use aligned allocation, otherwise fftw is slower since it cannot use SSE
    nsgt_complex_t *backwards = (nsgt_complex_t *) malloc(len * sizeof(nsgt_complex_t));
    nsgt_complex_t *input_complex = (nsgt_complex_t *) malloc(len * sizeof(nsgt_complex_t));

    /* the buffers are cleared by fftw if optimization is enabled*/
    nsgt_forward_set_buffers(nsgt, input_complex, result); //lots of time is spent here, initializing the fftw library
    nsgt_backward_set_buffers(nsgt, result, backwards); //lots of time is spent here, initializing the fftw library

    copy_real_to_complex(input, input_complex, len);

    printf("Initialization in %.3f ms\r\n", timer_end(init_time));

    time_point_t forward_time = timer_start();
    nsgt_forward(nsgt);
    printf("Forward in %.3fms, ", timer_end(forward_time));

    time_point_t backward_time = timer_start();
    nsgt_backward(nsgt);
    printf("Backwards in %.3fms, ", timer_end(backward_time));

    double err = nsgt_util_calc_error(backwards, input_complex, len);
    err = 10 * log10(err); //db
    printf("Error: %.3fdb\r\n\r\n", err);

    free(backwards);
    free(input_complex);

    nsgt_free_result(result);
    nsgt_free(nsgt);
}

#endif //NSGT_EXAMPLE_NSGT_ERROR_H
