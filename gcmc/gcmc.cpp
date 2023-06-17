/*
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved
        Mingtian Zhao, Abhishek A. Kognole,
        Aoxiang Tao, Alexander D. MacKerell Jr.
    E-mail:
        zhaomt@outerbanks.umaryland.edu
        alex@outerbanks.umaryland.edu
*/


#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

struct Atom {
    float position[3];
    float charge;
    int type;
};

struct AtomArray {
    
    char name[4];

    float muex;
    float conc;
    float confBias;
    float mcTime;
    
    int totalNum;
    int maxNum;

    int num_atoms;
    Atom atoms[20];
};

