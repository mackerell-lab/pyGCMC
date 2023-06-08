/*
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved
        Mingtian Zhao, Abhishek A. Kognole,
        Aoxiang Tao, Alexander D. MacKerell Jr.
    E-mail:
        zhaomt@outerbanks.umaryland.edu
        alex@outerbanks.umaryland.edu
*/


#include <iostream>
#include <cuda_runtime.h>
#include <thrust/device_vector.h>

extern "C" {


    __global__ void vector_add_kernel(const float *a, const float *b, float *c, int N) {
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i < N) {
            c[i] = a[i] + b[i];
        }
    }

    void vector_add_cuda(const float *a, const float *b, float *c, int N) {
        float *d_a, *d_b, *d_c;

        cudaMalloc(&d_a, N * sizeof(float));
        cudaMalloc(&d_b, N * sizeof(float));
        cudaMalloc(&d_c, N * sizeof(float));

        cudaMemcpy(d_a, a, N * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_b, b, N * sizeof(float), cudaMemcpyHostToDevice);

        int blockSize = 256;
        int numBlocks = (N + blockSize - 1) / blockSize;
        vector_add_kernel<<<numBlocks, blockSize>>>(d_a, d_b, d_c, N);

        // Synchronize device
        cudaDeviceSynchronize();

        cudaMemcpy(c, d_c, N * sizeof(float), cudaMemcpyDeviceToHost);
        cudaFree(d_a);
        cudaFree(d_b);
        cudaFree(d_c);
    }

    __global__ void vector_sub_kernel(const float *a, const float *b, float *c, int N) {
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i < N) {
            c[i] = a[i] - b[i];
        }
    }

    void vector_sub_cuda(const float *a, const float *b, float *c, int N) {
        float *d_a, *d_b, *d_c;

        cudaMalloc(&d_a, N * sizeof(float));
        cudaMalloc(&d_b, N * sizeof(float));
        cudaMalloc(&d_c, N * sizeof(float));

        cudaMemcpy(d_a, a, N * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_b, b, N * sizeof(float), cudaMemcpyHostToDevice);

        int blockSize = 256;
        int numBlocks = (N + blockSize - 1) / blockSize;
        vector_sub_kernel<<<numBlocks, blockSize>>>(d_a, d_b, d_c, N);

        // Synchronize device
        cudaDeviceSynchronize();

        cudaMemcpy(c, d_c, N * sizeof(float), cudaMemcpyDeviceToHost);
        cudaFree(d_a);
        cudaFree(d_b);
        cudaFree(d_c);
    }
}

