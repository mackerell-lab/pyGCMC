/*
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved
        Mingtian Zhao, Alexander D. MacKerell Jr.
    E-mail:
        zhaomt@outerbanks.umaryland.edu
        alex@outerbanks.umaryland.edu
*/

#ifndef GCMC_H_
#define GCMC_H_

#include <iostream>
#include <vector>
#include <unordered_set>
#include <cuda_runtime.h>
#include <curand_kernel.h>

// Constants
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062
#define temperature 300.0
#define OCCUPANCY_MAX 99999999
#define ENERGY_MAX 99999999999999
#define KCAL_TO_KJ 4.184
#define KJ_TO_KCAL 0.239
#define MOLES_TO_MOLECULES 0.0006023
#define MOLECULES_TO_MOLES 1660.539
#define BOLTZMANN 0.0083115
#define beta 1.0 / (BOLTZMANN * temperature)

#define numThreadsPerBlock 128

// Atom structure
struct Atom {
    float position[3];
    float charge;
    int type;
};

// AtomArray structure
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

// InfoStruct structure
struct InfoStruct {
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
    int PME;
    uint seed;
};

// Residue structure
struct Residue {
    float position[3];
    int atomNum;
    int atomStart;
};

#endif // GCMC_H_
