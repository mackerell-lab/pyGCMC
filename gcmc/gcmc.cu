/*
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved
        Mingtian Zhao, Abhishek A. Kognole,
        Aoxiang Tao, Alexander D. MacKerell Jr.
    E-mail:
        zhaomt@outerbanks.umaryland.edu
        alex@outerbanks.umaryland.edu
*/


#include <cuda_runtime.h>
// #include <unistd.h>
// #include <thrust/device_vector.h>
#include "gcmc.h"

// #include <cstdio>
// #define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
// inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
// {
//    if (code != cudaSuccess) 
//    {
//       fprintf(stderr,"GPU assert: %s %s %d\n", cudaGetErrorString(code), file, line);
//       if (abort) exit(code);
//    }
// }

extern "C"{
    // void runGCMC_cuda(const InfoStruct *info, AtomArray *fragmentInfo, residue *residueInfo, Atom *atomInfo, const float *grid, const float *ff, const int *moveArray){}
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
        // cudaMalloc(&GmoveArray, sizeof(int)*info->mcsteps);

        // printf("cudaMalloc done\n");

        // sleep(360);

        cudaMemcpy(Ginfo, info, sizeof(InfoStruct), cudaMemcpyHostToDevice);
        // printf("cudaMemcpy Ginfo done\n");
        cudaMemcpy(GfragmentInfo, fragmentInfo, sizeof(AtomArray)*info->fragTypeNum, cudaMemcpyHostToDevice);
        // printf("cudaMemcpy GfragmentInfo done\n");
        cudaMemcpy(GresidueInfo, residueInfo, sizeof(residue)*info->totalResNum, cudaMemcpyHostToDevice) ;
        // printf("cudaMemcpy GresidueInfo done\n");
        cudaMemcpy(GatomInfo, atomInfo, sizeof(Atom)*info->totalAtomNum, cudaMemcpyHostToDevice);
        // printf("cudaMemcpy GatomInfo done\n");
        cudaMemcpy(Ggrid, grid, sizeof(float)*info->totalGridNum * 3, cudaMemcpyHostToDevice);
        // printf("cudaMemcpy Ggrid done\n");
        cudaMemcpy(Gff, ff, sizeof(float)*info->ffXNum*info->ffYNum *2, cudaMemcpyHostToDevice);
        // printf("cudaMemcpy Gff done\n");
        // cudaMemcpy(GmoveArray, moveArray, sizeof(int)*info->mcsteps, cudaMemcpyHostToDevice);
        // printf("cudaMemcpy GmoveArray done\n");

        for (int stepi = 0 ; stepi < info->mcsteps; ++stepi){
            // Start MC steps
            int moveFragType = moveArray[stepi] / 4;
            int moveMoveType = moveArray[stepi] % 4;
            int confBias = fragmentInfo[moveFragType].confBias;

            // // perform move
            // bool accepted = false;
            // switch (moveMoveType)
            // {
            // case 0: // Insert
            //     accepted = move_add(Ginfo, GfragmentInfo, GresidueInfo, GatomInfo, Ggrid, Gff, moveFragType, confBias);
            //     break;

            // case 1: // Del
            //     accepted = move_del(frag_index);
            //     break;

            // case 2: // Trn
            //     accepted = move_trans(frag_index);
            //     break;

            // case 3: // Rot
            //     accepted = move_rotate(frag_index);
            //     break;
            // }






        }




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
        // cudaFree(GmoveArray);
        
    }
}

// extern "C" {


//     void runGCMC_cuda(const InfoStruct *info, AtomArray *fragmentInfo, residue *residueInfo, Atom *atomInfo, const float *grid, const float *ff, const int *moveArray){


//         // InfoStruct *Ginfo;
//         // AtomArray *GfragmentInfo;
//         // residue *GresidueInfo; 
//         // Atom *GatomInfo;
//         // float *Ggrid;
//         // float *Gff;
//         // int *GmoveArray;


//         // cudaMalloc(&Ginfo, sizeof(InfoStruct));

//         // cudaMalloc((void**)&Ginfo, sizeof(InfoStruct));
//         // cudaMalloc((void**)&GfragmentInfo, sizeof(AtomArray)*info->fragTypeNum);
//         // cudaMalloc((void**)&GresidueInfo, sizeof(residue)*info->totalResNum);
//         // cudaMalloc((void**)&GatomInfo, sizeof(Atom)*info->totalAtomNum);
//         // cudaMalloc((void**)&Ggrid, sizeof(float)*info->totalGridNum * 3);
//         // cudaMalloc((void**)&Gff, sizeof(float)*info->ffXNum*info->ffYNum *2);
//         // cudaMalloc((void**)&GmoveArray, sizeof(int)*info->mcsteps);

//         // cudaMemcpy(Ginfo, info, sizeof(InfoStruct), cudaMemcpyHostToDevice);
//         // cudaMemcpy(GfragmentInfo, fragmentInfo, sizeof(AtomArray)*info->fragTypeNum, cudaMemcpyHostToDevice);
//         // cudaMemcpy(GresidueInfo, residueInfo, sizeof(residue)*info->totalResNum, cudaMemcpyHostToDevice);
//         // cudaMemcpy(GatomInfo, atomInfo, sizeof(Atom)*info->totalAtomNum, cudaMemcpyHostToDevice);
//         // cudaMemcpy(Ggrid, grid, sizeof(float)*info->totalGridNum * 3, cudaMemcpyHostToDevice);
//         // cudaMemcpy(Gff, ff, sizeof(float)*info->ffXNum*info->ffYNum *2, cudaMemcpyHostToDevice);
//         // cudaMemcpy(GmoveArray, moveArray, sizeof(int)*info->mcsteps, cudaMemcpyHostToDevice);

        
//         // cudaDeviceSynchronize();

//         // sleep(60);

//         // cudaMemcpy(fragmentInfo, GfragmentInfo, sizeof(AtomArray)*info->fragTypeNum, cudaMemcpyDeviceToHost);
//         // cudaMemcpy(residueInfo, GresidueInfo, sizeof(residue)*info->totalResNum, cudaMemcpyDeviceToHost);
//         // cudaMemcpy(atomInfo, GatomInfo, sizeof(Atom)*info->totalAtomNum, cudaMemcpyDeviceToHost);

//         // cudaFree(Ginfo);
//         // cudaFree(GfragmentInfo);
//         // cudaFree(GresidueInfo);
//         // cudaFree(GatomInfo);
//         // cudaFree(Ggrid);
//         // cudaFree(Gff);
//         // cudaFree(GmoveArray);

//     }



// }

