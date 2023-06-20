
#include <iostream>
#include <cuda_runtime.h>
#include <curand_kernel.h>



#ifndef GCMC_H_
#define GCMC_H_


#define numThreadsPerBlock 128
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062

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
    int confBias;
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

    int showInfo;

    float cavityFactor;
    
    int fragTypeNum;
    
    int totalGridNum;
    int totalResNum;
    int totalAtomNum;
    
    int ffXNum;
    int ffYNum;

    uint seed;
};

struct residue{
    float position[3];
    int atomNum;
    int atomStart;
};

#endif