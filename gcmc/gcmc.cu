/*
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved
        Mingtian Zhao, Alexander D. MacKerell Jr.
    E-mail:
        zhaomt@outerbanks.umaryland.edu
        alex@outerbanks.umaryland.edu
*/

#include "gcmc.h"
#include "gcmc_move.h"

extern "C"{

    // Initialize random number generator states for CUDA threads
    __global__ void setup_rng_states(curandState *states, unsigned long long seed) {
        int global_threadIdx  = blockIdx.x * blockDim.x + threadIdx.x;
        curand_init(seed, global_threadIdx, 0, &states[global_threadIdx]);
    }

    // Main GCMC simulation function
    void runGCMC_cuda(const InfoStruct *info, AtomArray *fragmentInfo, residue *residueInfo, Atom *atomInfo, const float *grid, const float *ff, const int *moveArray){
        
        InfoStruct *Ginfo;
        AtomArray *GfragmentInfo;
        residue *GresidueInfo; 
        Atom *GatomInfo;
        float *Ggrid;
        float *Gff;
        int *GmoveArray;

        cudaMalloc(&Ginfo, sizeof(InfoStruct));
        cudaMalloc(&GfragmentInfo, sizeof(AtomArray)*info->fragTypeNum);
        cudaMalloc(&GresidueInfo, sizeof(residue)*info->totalResNum);
        cudaMalloc(&GatomInfo, sizeof(Atom)*info->totalAtomNum);
        cudaMalloc(&Ggrid, sizeof(float)*info->totalGridNum * 3);
        cudaMalloc(&Gff, sizeof(float)*info->ffXNum*info->ffYNum *2);

        cudaMemcpy(Ginfo, info, sizeof(InfoStruct), cudaMemcpyHostToDevice);
        cudaMemcpy(GfragmentInfo, fragmentInfo, sizeof(AtomArray)*info->fragTypeNum, cudaMemcpyHostToDevice);
        cudaMemcpy(GresidueInfo, residueInfo, sizeof(residue)*info->totalResNum, cudaMemcpyHostToDevice) ;
        cudaMemcpy(GatomInfo, atomInfo, sizeof(Atom)*info->totalAtomNum, cudaMemcpyHostToDevice);
        cudaMemcpy(Ggrid, grid, sizeof(float)*info->totalGridNum * 3, cudaMemcpyHostToDevice);
        cudaMemcpy(Gff, ff, sizeof(float)*info->ffXNum*info->ffYNum *2, cudaMemcpyHostToDevice);

        int maxConf = 0;
        for (int fragType = 0; fragType < info->fragTypeNum; fragType ++ ){
            if (fragmentInfo[fragType].confBias > maxConf){
                maxConf = fragmentInfo[fragType].confBias;
            }
        }

        AtomArray *GTempFrag;
        cudaMalloc(&GTempFrag, sizeof(AtomArray)*maxConf);

        Atom *GTempInfo;
        cudaMalloc(&GTempInfo, sizeof(Atom)*maxConf);

        Atom *TempInfo;
        TempInfo = (Atom *)malloc(sizeof(Atom)*maxConf);

        for (int i = 0;i < maxConf; i++){
            TempInfo[i].type = 0;
        }

        cudaMemcpy(GTempInfo, TempInfo, sizeof(Atom)*maxConf, cudaMemcpyHostToDevice);

        curandState *d_rng_states;
        
        cudaMalloc((void **)&d_rng_states, maxConf * sizeof(curandState) * numThreadsPerBlock);

        srand(info->seed);

        setup_rng_states<<<maxConf, numThreadsPerBlock>>>(d_rng_states, info->seed);

        int step_threshold = info->mcsteps / 20;

        for (int stepi = 0 ; stepi < info->mcsteps; ++stepi){
            int moveFragType = moveArray[stepi] / 4;
            int moveMoveType = moveArray[stepi] % 4;
            int confBias = fragmentInfo[moveFragType].confBias;

            bool accepted = false;
            switch (moveMoveType)
            {
            case 0: // Insert
                accepted = move_add(info, Ginfo,fragmentInfo, GfragmentInfo, GresidueInfo, GatomInfo, Ggrid, Gff, moveFragType, GTempFrag, TempInfo, GTempInfo, d_rng_states);
                break;
            case 1: // Delete
                accepted = move_del(info, Ginfo,fragmentInfo, GfragmentInfo, GresidueInfo, GatomInfo, Ggrid, Gff, moveFragType, GTempFrag, TempInfo, GTempInfo, d_rng_states);
                break;
            case 2: // Translate
                accepted = move_trn(info, Ginfo,fragmentInfo, GfragmentInfo, GresidueInfo, GatomInfo, Ggrid, Gff, moveFragType, GTempFrag, TempInfo, GTempInfo, d_rng_states);
                break;
            case 3: // Rotate
                accepted = move_rot(info, Ginfo,fragmentInfo, GfragmentInfo, GresidueInfo, GatomInfo, Ggrid, Gff, moveFragType, GTempFrag, TempInfo, GTempInfo, d_rng_states);
                break;
            }
        }
        printf("\n");

        cudaDeviceSynchronize();

        cudaMemcpy(fragmentInfo, GfragmentInfo, sizeof(AtomArray)*info->fragTypeNum, cudaMemcpyDeviceToHost);
        cudaMemcpy(residueInfo, GresidueInfo, sizeof(residue)*info->totalResNum, cudaMemcpyDeviceToHost);
        cudaMemcpy(atomInfo, GatomInfo, sizeof(Atom)*info->totalAtomNum, cudaMemcpyDeviceToHost);

        cudaFree(Ginfo);
        cudaFree(GfragmentInfo);
        cudaFree(GresidueInfo);
        cudaFree(GatomInfo);
        cudaFree(Ggrid);
        cudaFree(Gff);
        cudaFree(GTempFrag);
        cudaFree(GTempInfo);
        cudaFree(d_rng_states);

        free(TempInfo);
    }
}


