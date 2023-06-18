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
#include <unistd.h>
// #include <thrust/device_vector.h>

struct Atom {
    float position[3];
    float charge;
    int type;
};

struct AtomArray {
    
    char name[4];

    int startRes;

    float muex;
    float conc;
    float confBias;
    float mcTime;
    
    int totalNum;
    int maxNum;

    int num_atoms;
    Atom atoms[20];
};

struct InfoStruct{
    int mcsteps;
    float cutoff;
    float grid_dx;
    float startxyz[3];
    float cryst[3];

    float cavityFactor;
    
    int fragTypeNum;
    
    int totalGridNum;
    int totalResNum;
    int totalAtomNum;
    
    int ffXNum;
    int ffYNum;
};

struct residue{
    float position[3];
    int atomNum;
    int atomStart;
};


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
        cudaMalloc(&GmoveArray, sizeof(int)*info->mcsteps);

        sleep(60);

        cudaMemcpy(Ginfo, info, sizeof(InfoStruct), cudaMemcpyHostToDevice);
        cudaMemcpy(GfragmentInfo, fragmentInfo, sizeof(AtomArray)*info->fragTypeNum, cudaMemcpyHostToDevice);
        cudaMemcpy(GresidueInfo, residueInfo, sizeof(residue)*info->totalResNum, cudaMemcpyHostToDevice);
        cudaMemcpy(GatomInfo, atomInfo, sizeof(Atom)*info->totalAtomNum, cudaMemcpyHostToDevice);
        cudaMemcpy(Ggrid, grid, sizeof(float)*info->totalGridNum * 3, cudaMemcpyHostToDevice);
        cudaMemcpy(Gff, ff, sizeof(float)*info->ffXNum*info->ffYNum *2, cudaMemcpyHostToDevice);
        cudaMemcpy(GmoveArray, moveArray, sizeof(int)*info->mcsteps, cudaMemcpyHostToDevice);

        
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

