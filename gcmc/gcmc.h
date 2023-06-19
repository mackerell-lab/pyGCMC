
#include <iostream>

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

    int showInfo;

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

