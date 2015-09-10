#ifndef NSGT_EXAMPLE_NSGT_STREAM_ERROR_H
#define NSGT_EXAMPLE_NSGT_STREAM_ERROR_H

#include "example_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

void example_stream_error(double *input_real, size_t len, int samplerate) {
    printf("Test nsgt stream (sliCQ)\r\n");
    time_point_t init_time = timer_start();
    int blocksize = 32768;

    NsgtScaleDesc scale = { 0 };
    scale.scale_type = NSGT_SCALE_OCT;
    scale.fmin = 100.;
    scale.fmax = 22050.;
    scale.num_bins = 24;
    int32_t optimize = 0; //NSGT_OPTIMIZE_FFT;
    NsgtStream *nsgt = nsgt_stream_create(&scale, samplerate, (size_t) blocksize, 4096, optimize);

    nsgt_complex_t *f_slice = (nsgt_complex_t *) malloc(blocksize * sizeof(nsgt_complex_t));
    nsgt_complex_t *f_slice_old = (nsgt_complex_t *) malloc(blocksize * sizeof(nsgt_complex_t));
    nsgt_complex_t *out_slice = (nsgt_complex_t *) malloc(blocksize * sizeof(nsgt_complex_t));

    NsgtResult *res = nsgt_alloc_result(nsgt->nsgt);
    nsgt_stream_forward_set_buffers(nsgt, f_slice, res);
    nsgt_stream_backward_set_buffers(nsgt, res, out_slice);

    printf("Initialization in %.3f ms\r\n", timer_end(init_time));

    for (size_t i = 0; i < len / blocksize + 1; i++) {
        //prepare input for nsgt, the last block only contains padding
        memset(f_slice, 0, blocksize * sizeof(nsgt_complex_t));
        int remaining = (int) (len - i * blocksize);
        if (blocksize < remaining)
            remaining = blocksize;
        if (remaining < 0)
            remaining = 0;
        copy_real_to_complex(input_real + i * blocksize, f_slice, (size_t) remaining);

        time_point_t forward_time = timer_start();
        nsgt_stream_forward(nsgt); //2 * hhop in
        printf("Forward in %.3fms, ", timer_end(forward_time));

        time_point_t backwards_time = timer_start();
        nsgt_stream_backward(nsgt); //2 * hhop out
        printf("Backwards in %.3fms, ", timer_end(backwards_time));

        if (i > 0) { // The second backwards output contains the info from the first forward input block and so on
            double err = nsgt_util_calc_error(out_slice, f_slice_old, (size_t) blocksize);
            err = 10 * log10(err); //db
            printf("Error: %.3fdb", err);
        }
        printf("\r\n");

        //store samples from current block, so we can test them in the next iteration
        memcpy(f_slice_old, f_slice, sizeof(nsgt_complex_t) * blocksize);
    }
    printf("\r\n");
    free(f_slice);
    free(f_slice_old);
    free(out_slice);
    nsgt_free_result(res);
    nsgt_stream_free(nsgt);
}

#endif //NSGT_EXAMPLE_NSGT_STREAM_ERROR_H
