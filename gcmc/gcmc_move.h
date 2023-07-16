/*
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved
        Mingtian Zhao, Abhishek A. Kognole,
        Aoxiang Tao, Alexander D. MacKerell Jr.
    E-mail:
        zhaomt@outerbanks.umaryland.edu
        alex@outerbanks.umaryland.edu
*/

// #include <cuda_runtime.h>
// #include <unistd.h>
// #include <thrust/device_vector.h>
#include "gcmc.h"



extern "C"{

        __device__ __host__ inline float Min(float a,float b)
        {
            return !(b<a)?a:b;	
        }

        // Device kernel function
        __device__ inline void rotate_atoms(Atom *atoms, int num_atoms, float axis[3], float angle) {
            // Normalize the axis vector
            float norm = sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
            axis[0] /= norm;
            axis[1] /= norm;
            axis[2] /= norm;

            // Convert the angle to radians and multiply by 2*PI so it rotates between 0 and 2*PI
            // angle *= 2 * M_PI;

            int idx = threadIdx.x;

            // Rotate the atoms with number less than 128

            if (idx < num_atoms) {
                float atom_position[3] = {atoms[idx].position[0], atoms[idx].position[1], atoms[idx].position[2]};

                // Compute the rotation matrix
                float c = cos(angle);
                float s = sin(angle);
                float C = 1 - c;
                float R[3][3];

                R[0][0] = axis[0] * axis[0] * C + c;
                R[0][1] = axis[0] * axis[1] * C - axis[2] * s;
                R[0][2] = axis[0] * axis[2] * C + axis[1] * s;

                R[1][0] = axis[1] * axis[0] * C + axis[2] * s;
                R[1][1] = axis[1] * axis[1] * C + c;
                R[1][2] = axis[1] * axis[2] * C - axis[0] * s;

                R[2][0] = axis[2] * axis[0] * C - axis[1] * s;
                R[2][1] = axis[2] * axis[1] * C + axis[0] * s;
                R[2][2] = axis[2] * axis[2] * C + c;

                // Apply the rotation matrix
                atoms[idx].position[0] = atom_position[0] * R[0][0] + atom_position[1] * R[0][1] + atom_position[2] * R[0][2];
                atoms[idx].position[1] = atom_position[0] * R[1][0] + atom_position[1] * R[1][1] + atom_position[2] * R[1][2];
                atoms[idx].position[2] = atom_position[0] * R[2][0] + atom_position[1] * R[2][1] + atom_position[2] * R[2][2];
            }
        }

        // Device kernel function
        __device__ inline void rotate_atoms_shared(Atom *atoms, int num_atoms, float axis[3], float angle) {

            // Declare shared memory
            __shared__ float sh_axis[3];
            __shared__ float sh_R[3][3];

            int idx = threadIdx.x;

            // Assign axis to shared memory
            if (idx < 3) {
                sh_axis[idx] = axis[idx];
            }

            // Normalize the axis vector
            if (idx == 0) {
                float norm = sqrt(sh_axis[0] * sh_axis[0] + sh_axis[1] * sh_axis[1] + sh_axis[2] * sh_axis[2]);
                sh_axis[0] /= norm;
                sh_axis[1] /= norm;
                sh_axis[2] /= norm;

                // Compute the rotation matrix
                float c = cos(angle);
                float s = sin(angle);
                float C = 1 - c;

                sh_R[0][0] = sh_axis[0] * sh_axis[0] * C + c;
                sh_R[0][1] = sh_axis[0] * sh_axis[1] * C - sh_axis[2] * s;
                sh_R[0][2] = sh_axis[0] * sh_axis[2] * C + sh_axis[1] * s;
 
                sh_R[1][0] = sh_axis[1] * sh_axis[0] * C + sh_axis[2] * s;
                sh_R[1][1] = sh_axis[1] * sh_axis[1] * C + c;
                sh_R[1][2] = sh_axis[1] * sh_axis[2] * C - sh_axis[0] * s;

                sh_R[2][0] = sh_axis[2] * sh_axis[0] * C - sh_axis[1] * s;
                sh_R[2][1] = sh_axis[2] * sh_axis[1] * C + sh_axis[0] * s;
                sh_R[2][2] = sh_axis[2] * sh_axis[2] * C + c;
            }

            __syncthreads(); // Ensure that the shared memory has been initialized before continuing

            // Rotate the atoms with number less than 128
            if (idx < num_atoms) {
                float atom_position[3] = {atoms[idx].position[0], atoms[idx].position[1], atoms[idx].position[2]};
                
                // Apply the rotation matrix
                atoms[idx].position[0] = atom_position[0] * sh_R[0][0] + atom_position[1] * sh_R[0][1] + atom_position[2] * sh_R[0][2];
                atoms[idx].position[1] = atom_position[0] * sh_R[1][0] + atom_position[1] * sh_R[1][1] + atom_position[2] * sh_R[1][2];
                atoms[idx].position[2] = atom_position[0] * sh_R[2][0] + atom_position[1] * sh_R[2][1] + atom_position[2] * sh_R[2][2];
            }
        }

        __device__ void randomFragment(const InfoStruct &SharedInfo, AtomArray &SharedFragmentInfo, Atom *GTempInfo, const float *Ggrid, curandState *rng_states) {

            int tid = threadIdx.x;

            __shared__ float randomR[3];
            __shared__ float randomThi[3];
            __shared__ float randomPhi;
            __shared__ int gridN;

            if (tid < 3){
                randomR[tid] = curand_uniform(rng_states) * SharedInfo.grid_dx;
            }
            if (tid >= 3 && tid < 6){
                randomThi[tid - 3] = curand_uniform(rng_states);
            }
            if (tid == 6){
                randomPhi = curand_uniform(rng_states) * 2 * PI;
            }
            if (tid == 7){
                gridN = curand(rng_states) % SharedInfo.totalGridNum;
            }

            __syncthreads();
            if (tid < 3){
                randomR[tid] += Ggrid[gridN * 3 + tid];
            }

            __syncthreads();

            // if ( threadIdx.x == 0 && blockIdx.x == 636){

            //     printf("Before before 4 rotation:\n");
            //     // printf("randomR: %8.3f%8.3f%8.3f\n", randomR[0], randomR[1], randomR[2]);
            //     // printf("randomThi: %8.3f%8.3f%8.3f\n", randomThi[0], randomThi[1], randomThi[2]);
            //     // printf("randomPhi: %8.3f\n", randomPhi);
            //     // printf("gridN: %d\n", gridN);
            //     for (int i = 0; i < SharedFragmentInfo.num_atoms; i++){
            //         printf("%8.3f%8.3f%8.3f\n", SharedFragmentInfo.atoms[i].position[0], SharedFragmentInfo.atoms[i].position[1], SharedFragmentInfo.atoms[i].position[2]);
            //     }
            // }
            // __syncthreads();


            // if (threadIdx.x == 0 && blockIdx.x == 636){
            //     printf("randomR: %8.3f%8.3f%8.3f\n", randomR[0], randomR[1], randomR[2]);
            //     printf("randomThi: %8.3f%8.3f%8.3f\n", randomThi[0], randomThi[1], randomThi[2]);
            //     printf("randomPhi: %8.3f\n", randomPhi);
            //     printf("gridN: %d\n", gridN);

            //     printf("Before rotation:\n");
            //     for (int i = 0; i < SharedFragmentInfo.num_atoms; i++){
            //         printf("%8.3f%8.3f%8.3f\n", SharedFragmentInfo.atoms[i].position[0], SharedFragmentInfo.atoms[i].position[1], SharedFragmentInfo.atoms[i].position[2]);
            //     }
            // }
            // __syncthreads();
            
            // for (int i=0;i<2000;i++)
            
            // rotate_atoms(SharedFragmentInfo.atoms, SharedFragmentInfo.num_atoms, randomThi, randomPhi);
            rotate_atoms_shared(SharedFragmentInfo.atoms, SharedFragmentInfo.num_atoms, randomThi, randomPhi);

            __syncthreads();



            // if (tid == 0 && blockIdx.x == 636){
            //     printf("randomR: %8.3f%8.3f%8.3f\n", randomR[0], randomR[1], randomR[2]);
            //     printf("randomThi: %8.3f%8.3f%8.3f\n", randomThi[0], randomThi[1], randomThi[2]);
            //     printf("randomPhi: %8.3f\n", randomPhi);
            //     printf("gridN: %d\n", gridN);

            //     printf("After rotation:\n");
            //     for (int i = 0; i < SharedFragmentInfo.num_atoms; i++){
            //         printf("%8.3f%8.3f%8.3f\n", SharedFragmentInfo.atoms[i].position[0], SharedFragmentInfo.atoms[i].position[1], SharedFragmentInfo.atoms[i].position[2]);
            //     }
            // }
            // __syncthreads();
            for (int i= tid ; i < SharedFragmentInfo.num_atoms; i += blockDim.x){
                SharedFragmentInfo.atoms[i].position[0] += randomR[0];
                SharedFragmentInfo.atoms[i].position[1] += randomR[1];
                SharedFragmentInfo.atoms[i].position[2] += randomR[2];
            }

            // if (tid < SharedFragmentInfo.num_atoms){
            //     SharedFragmentInfo.atoms[tid].position[0] += randomR[0];
            //     SharedFragmentInfo.atoms[tid].position[1] += randomR[1];
            //     SharedFragmentInfo.atoms[tid].position[2] += randomR[2];
            // }

            if (tid < 3){
                GTempInfo->position[tid] = randomR[tid];
            }
            if (tid == 4)
                GTempInfo->type = -1;




        }

        __device__ __host__ inline float fast_round(const float a) {
            return a >= 0 ? (int)(a + 0.5f) : (int)(a - 0.5f);
        }

        __device__ __host__ inline float distanceP(const float x[3], const float y[3], const float period[3]){
            float dx = x[0] - y[0];
            float dy = x[1] - y[1];
            float dz = x[2] - y[2];

            dx -= fast_round(dx / period[0]) * period[0];
            dy -= fast_round(dy / period[1]) * period[1];
            dz -= fast_round(dz / period[2]) * period[2];

            return sqrtf(dx * dx + dy * dy + dz * dz);
        }

        
        //If NBFIX entry exists for the pair, then this calculation
        //!V(Lennard-Jones) = 4*Eps,i,j[(sigma,i,j/ri,j)**12 - (sigma,i,j/ri,j)**6]
        // sigma and Eps are the nbfix entries
        // units from force field: sigma (nm), epsilon (kJ)

        __device__ inline float calc_vdw_nbfix (float sigma, float epsilon, float dist_sqrd)
        {
            // convert sigma from nm to Angstroms to match dist
            // sigma *= 10;

            float sigma_sqrd = sigma * sigma * 100;
            float sigma_dist_sqrd = sigma_sqrd / dist_sqrd;
            float E_vdw = 4*epsilon * ( pow(sigma_dist_sqrd, 6) - pow(sigma_dist_sqrd, 3) );

            //cout << "(calc_vdw_nbfix) epsilon: " << epsilon << " sigma_sqrd: " << sigma_sqrd << " dist_sqrd: " << dist_sqrd << " E_vdw: " << E_vdw << endl;

            return E_vdw;
        }

        
        // E_elec = (1/4*pi*eps_0)*(q1*q2/d)
        // units from force field: charge (e)
        // units from this program: distance (A)
        // eps_0 = 8.85e-12 (C^2/J*m)
        // kel = 1/(4*pi*eps_0)
        //     = 8.99e9 (J*m/C^2)
        //     = 1389.3 (kJ/mol * A/e^2)
        //     = 332.05 (kcal/mol * A/e^2)
        // from http://users.mccammon.ucsd.edu/~blu/Research-Handbook/physical-constant.html
        __device__ inline float calc_elec (float charge1, float charge2, float dist)
        {
            // float kel = 331.843 * KCAL_TO_KJ; //1389.3;  what's 331.843?

            // float E_elec = kel * (charge1 * charge2) / dist / eps;
            // float E_elec = 1389.3 * (charge1 * charge2) / dist;
            // E_elec = 0;


            // Differences in the electrostatic energies:
            //  
            // (*) The conversion from charge units to kcal/mol in CHARMM is based 
            // on the value 332.0716 whereas AMBER uses 18.2223**2 or 332.0522173.
            // The actual value is somewhat lower than both these values
            // (~331.843)!  To convert the charges to "CHARMM", they should be
            // multiplied by 1.000058372.  This was not done within this file.
            // [When this is done, the charges are not rounded and therefore
            // non-integral charges for the residues are apparent.]  To get around
            // this problem either the charges can be scaled within CHARMM (which
            // will still lead to non-integral charge) or in versions of CHARMM
            // beyond c25n3, and through the application of the "AMBER" keyword in
            // pref.dat, the AMBER constant can be used.  By default, the "fast"
            // routines cannot be used with the AMBER-style impropers.  In the
            // later versions of CHARMM, the AMBER keyword circumvents this
            // problem.
            // 
            // Ref: https://home.chpc.utah.edu/~cheatham/cornell_rtf


            float E_elec = 1388.431112 * (charge1 * charge2) / dist;
            return E_elec;
        }


        __device__ inline void calcProtEnergy(const InfoStruct SharedInfo, AtomArray &SharedFragmentInfo , AtomArray *GfragmentInfo,
                                    residue *GresidueInfo, Atom *GatomInfo, const float *Gff, Atom *GTempInfo, float *sh_energy){
            
            int tid = threadIdx.x;
            __shared__ int maxResidueNum;
            if (tid == 0)
                maxResidueNum = GfragmentInfo->startRes;
            
            
            __syncthreads();

            for (int resi = tid;resi < maxResidueNum; resi+= numThreadsPerBlock){

                // if (distanceP(GTempInfo->position, GresidueInfo[resi].position, SharedInfo.cryst) > SharedInfo.cutoff) {
                if (distanceP(GTempInfo->position, GresidueInfo[resi].position, SharedInfo.cryst) > SharedInfo.cutoff) {
                    continue;
                }

                int resiStart = GresidueInfo[resi].atomStart;
                int resiEnd = GresidueInfo[resi].atomStart + GresidueInfo[resi].atomNum;
                float resiEnergy = 0;
                for (int atomi = resiStart; atomi < resiEnd; atomi++){
                    // int atomType = GatomInfo[atomi].type;
                    // float atomCharge = GatomInfo[atomi].charge;
                    // float atomEnergy = 0;
                    for (int atomj = 0; atomj < SharedFragmentInfo.num_atoms; atomj++){

                        float distance = distanceP(GatomInfo[atomi].position, SharedFragmentInfo.atoms[atomj].position, SharedInfo.cryst);
                        int typeij = SharedFragmentInfo.atoms[atomj].type * SharedInfo.ffYNum + GatomInfo[atomi].type;
                        // float sigma = Gff[typeij * 2];
                        // float epsilon = Gff[typeij * 2 + 1];
                        resiEnergy += calc_vdw_nbfix(Gff[typeij * 2], Gff[typeij * 2 + 1], distance * distance);
                        resiEnergy += calc_elec(SharedFragmentInfo.atoms[atomj].charge, GatomInfo[atomi].charge, distance);
                        // float sigma = Gff[typeij * 2];
                        // float energy = Gff[atomType * SharedInfo.atomTypeNum + SharedFragmentInfo.atoms[atomj].type] * atomCharge * SharedFragmentInfo.atoms[atomj].charge / distance;
                        // resiEnergy += energy;
                    }
                    // resiEnergy += atomEnergy;
                }
                sh_energy[tid] += resiEnergy;
            }


        }


        __device__ inline void calcFragEnergy(const InfoStruct SharedInfo, AtomArray &SharedFragmentInfo , AtomArray *GfragmentInfo,
                                    residue *GresidueInfo, Atom *GatomInfo, const float *Gff, Atom *GTempInfo, float *sh_energy){
            
            int tid = threadIdx.x;
            __shared__ int startResidueNum;
            __shared__ int endResidueNum;


            for (int fragType = 0; fragType < SharedInfo.fragTypeNum; fragType++){

                if (tid == 0){
                    startResidueNum = GfragmentInfo[fragType].startRes;
                    endResidueNum = GfragmentInfo[fragType].startRes + GfragmentInfo[fragType].totalNum;
                }

                __syncthreads();

                for (int resi = tid + startResidueNum;resi < endResidueNum; resi+= numThreadsPerBlock){

                    if (resi == SharedFragmentInfo.startRes)
                        continue;

                    // if (distanceP(GTempInfo->position, GresidueInfo[resi].position, SharedInfo.cryst) > SharedInfo.cutoff) {
                    if (distanceP(GTempInfo->position, GresidueInfo[resi].position, SharedInfo.cryst) > SharedInfo.cutoff) {
                        continue;
                    }

                    int resiStart = GresidueInfo[resi].atomStart;
                    int resiEnd = GresidueInfo[resi].atomStart + GresidueInfo[resi].atomNum;
                    float resiEnergy = 0;
                    for (int atomi = resiStart; atomi < resiEnd; atomi++){
                        // int atomType = GatomInfo[atomi].type;
                        // float atomCharge = GatomInfo[atomi].charge;
                        // float atomEnergy = 0;
                        for (int atomj = 0; atomj < SharedFragmentInfo.num_atoms; atomj++){

                            float distance = distanceP(GatomInfo[atomi].position, SharedFragmentInfo.atoms[atomj].position, SharedInfo.cryst);
                            int typeij = SharedFragmentInfo.atoms[atomj].type * SharedInfo.ffYNum + GatomInfo[atomi].type;
                            // float sigma = Gff[typeij * 2];
                            // float epsilon = Gff[typeij * 2 + 1];
                            resiEnergy += calc_vdw_nbfix(Gff[typeij * 2], Gff[typeij * 2 + 1], distance * distance);
                            resiEnergy += calc_elec(SharedFragmentInfo.atoms[atomj].charge, GatomInfo[atomi].charge, distance);
                            // float sigma = Gff[typeij * 2];
                            // float energy = Gff[atomType * SharedInfo.atomTypeNum + SharedFragmentInfo.atoms[atomj].type] * atomCharge * SharedFragmentInfo.atoms[atomj].charge / distance;
                            // resiEnergy += energy;
                        }
                        // resiEnergy += atomEnergy;
                    }
                    sh_energy[tid] += resiEnergy;
                }
                
                __syncthreads();
            }
            


        }
        

        __device__ void calcEnergy(const InfoStruct SharedInfo, AtomArray &SharedFragmentInfo , AtomArray *GfragmentInfo,
                                    residue *GresidueInfo, Atom *GatomInfo, const float *Gff, Atom *GTempInfo){

            __shared__ float sh_energy[numThreadsPerBlock];

            int tid = threadIdx.x;

            sh_energy[tid] = 0;

            __syncthreads();

            calcProtEnergy(SharedInfo, SharedFragmentInfo, GfragmentInfo, GresidueInfo, GatomInfo, Gff, GTempInfo, sh_energy);
            
            __syncthreads();

            calcFragEnergy(SharedInfo, SharedFragmentInfo, GfragmentInfo, GresidueInfo, GatomInfo, Gff, GTempInfo, sh_energy);

            __syncthreads();

            for (int s = numThreadsPerBlock / 2; s > 0; s >>= 1) {
                if (tid < s) {
                    sh_energy[tid] += sh_energy[tid + s];
                }
                __syncthreads();
            }

            if (tid == 0){

                GTempInfo->charge = sh_energy[0];
                // GTempInfo->type = -1;
                // if (blockIdx.x == 0){
                //     printf("energy: %f\n", sh_energy[0]);
                // }
            }

        }

        __global__ void Gmove_add(const InfoStruct *Ginfo, AtomArray *GfragmentInfo, residue *GresidueInfo, 
                        Atom *GatomInfo, const float *Ggrid, const float *Gff, const int moveFragType,
                        AtomArray *GTempFrag, Atom *GTempInfo, curandState *d_rng_states) {
                    
            __shared__ InfoStruct SharedInfo;
            __shared__ AtomArray SharedFragmentInfo;

            int threadId = numThreadsPerBlock * blockIdx.x + threadIdx.x;
            curandState *rng_states = &d_rng_states[threadId];


            int tid = threadIdx.x;

            if (threadIdx.x == 0) {
                SharedInfo = Ginfo[0];
            }
            
            if (threadIdx.x == 1) {
                SharedFragmentInfo = GfragmentInfo[moveFragType];
                SharedFragmentInfo.startRes = -1;
                
                // if (blockIdx.x == 636){
                //     printf("The center of the %d fragment is %8.3f%8.3f%8.3f\n", blockIdx.x, GTempInfo[blockIdx.x].position[0], GTempInfo[blockIdx.x].position[1], GTempInfo[blockIdx.x].position[2]);
                //     printf("The energy of the fragment is %f\n", GTempInfo[blockIdx.x].charge);
                //     printf("The type of the fragment is %d\n", GTempInfo[blockIdx.x].type);
                //     for (int i = 0; i < SharedFragmentInfo.num_atoms; i++){
                //         float total = SharedFragmentInfo.atoms[i].position[0] + SharedFragmentInfo.atoms[i].position[1] + SharedFragmentInfo.atoms[i].position[2] + 1;
                //         printf("%8.3f%8.3f%8.3f %8.3f\n",  SharedFragmentInfo.atoms[i].position[0], SharedFragmentInfo.atoms[i].position[1], SharedFragmentInfo.atoms[i].position[2], total);
                //     }
                // }
            }

            __syncthreads();

            randomFragment(SharedInfo, SharedFragmentInfo, &GTempInfo[blockIdx.x], Ggrid, rng_states);

            __syncthreads();

            // if (tid == 0 && blockIdx.x == 636){
            //     printf("The center of the %d fragment is %8.3f%8.3f%8.3f\n", blockIdx.x, GTempInfo[blockIdx.x].position[0], GTempInfo[blockIdx.x].position[1], GTempInfo[blockIdx.x].position[2]);
            //     printf("The energy of the fragment is %f\n", GTempInfo[blockIdx.x].charge);
            //     printf("The type of the fragment is %d\n", GTempInfo[blockIdx.x].type);
            //     for (int i = 0; i < SharedFragmentInfo.num_atoms; i++){
            //         printf("%8.3f%8.3f%8.3f\n", SharedFragmentInfo.atoms[i].position[0], SharedFragmentInfo.atoms[i].position[1], SharedFragmentInfo.atoms[i].position[2]);
            //     }
            // }

            if (tid == 0)
                GTempFrag[blockIdx.x] = SharedFragmentInfo;

            // if (SharedInfo.PME == 0)
                calcEnergy(SharedInfo, SharedFragmentInfo, GfragmentInfo, GresidueInfo, GatomInfo, Gff, &GTempInfo[blockIdx.x]);



            // __syncthreads();

            // if (tid == 0 && abs(GTempInfo[blockIdx.x].charge) > 150 ){
            //     printf("The center of the %d fragment is %8.3f%8.3f%8.3f\n", blockIdx.x, GTempInfo[blockIdx.x].position[0], GTempInfo[blockIdx.x].position[1], GTempInfo[blockIdx.x].position[2]);
            //     printf("The energy of the fragment is %f\n", GTempInfo[blockIdx.x].charge);
            //     printf("The type of the fragment is %d\n", GTempInfo[blockIdx.x].type);
            //     for (int i = 0; i < SharedFragmentInfo.num_atoms; i++){
            //         printf("%8.3f%8.3f%8.3f %8.3f\n", SharedFragmentInfo.atoms[i].position[0], SharedFragmentInfo.atoms[i].position[1], SharedFragmentInfo.atoms[i].position[2]);
            //     }
            // }

            // else
            //     calcEnergyPME(SharedInfo, SharedFragmentInfo, GfragmentInfo, GresidueInfo, GatomInfo, Gff, &GTempInfo[blockIdx.x]);

            // __syncthreads();

            // confBiasAdd(SharedInfo, SharedFragmentInfo, GfragmentInfo, GresidueInfo, GatomInfo, GTempInfo, GTempFrag, rng_states);

            // if (threadIdx.x == 0 && blockIdx.x == 0)

            //     printf("The position of the first atom is %f, %f, %f\n", SharedFragmentInfo.atoms[0].position[0], SharedFragmentInfo.atoms[0].position[1], SharedFragmentInfo.atoms[0].position[2]);







        }

        // __global__ void Gmove_add_check(const InfoStruct *Ginfo, AtomArray *GfragmentInfo, residue *GresidueInfo, 
        //                 Atom *GatomInfo, const int moveFragType,
        //                 AtomArray *GTempFrag, Atom *GTempInfo, curandState *d_rng_states) {
                    

        //     int tid = threadIdx.x;
        //     for (int i = tid ;i < 100;i+= numThreadsPerBlock)

        // }


        // Gupdate<<<1, numThreadsPerBlock>>>(GTempInfo,GfragmentInfo,GresidueInfo, GatomInfo)

        __global__ void GupdateAdd(AtomArray *GfragmentInfo, residue *GresidueInfo, 
                        Atom *GatomInfo,
                        AtomArray *GTempFrag, Atom *GTempInfo, const int moveFragType, const int totalNum) {
            
            int tid = threadIdx.x;
            int bid = blockIdx.x;

            __shared__ int bidStartRes;
            __shared__ int bidStartAtom;
            __shared__ int bidAtomNum;

            if (tid == 0 && bid == 0){
                GfragmentInfo[moveFragType].totalNum = totalNum;

            }

            if (GTempInfo[bid].type == -1)
                return;

            if (tid == 1){

                bidStartRes = GfragmentInfo[moveFragType].startRes + GTempInfo[bid].type;

            }
            __syncthreads();
            if (tid == 0){
                bidAtomNum = GresidueInfo[bidStartRes].atomNum;
                bidStartAtom = GresidueInfo[bidStartRes].atomStart;
                GresidueInfo[bidStartRes].position[0] = GTempInfo[bid].position[0];
                GresidueInfo[bidStartRes].position[1] = GTempInfo[bid].position[1];
                GresidueInfo[bidStartRes].position[2] = GTempInfo[bid].position[2];
            }
            __syncthreads();
            for (int i = tid;i<bidAtomNum;i+=numThreadsPerBlock){
                GatomInfo[bidStartAtom + i].position[0] = GTempFrag[bid].atoms[i].position[0];
                GatomInfo[bidStartAtom + i].position[1] = GTempFrag[bid].atoms[i].position[1];
                GatomInfo[bidStartAtom + i].position[2] = GTempFrag[bid].atoms[i].position[2];
            }


        }

        bool move_add(const InfoStruct *info, InfoStruct *Ginfo, AtomArray *fragmentInfo, AtomArray *GfragmentInfo, residue *GresidueInfo, Atom *GatomInfo, const float *Ggrid, const float *Gff,
                    const int moveFragType, AtomArray *GTempFrag, Atom *TempInfo, Atom *GTempInfo, curandState *d_rng_states){

            if (fragmentInfo[moveFragType].totalNum == fragmentInfo[moveFragType].maxNum)
                return false;
            
            const int nBlock = fragmentInfo[moveFragType].confBias;



            

            // printf("The nbar value is %f\n", nbar);
            // printf("The B value is %f\n", B);

            // printf("The size of a AtomArray is %d\n", sizeof(AtomArray));

            Gmove_add<<<nBlock, numThreadsPerBlock>>>(Ginfo, GfragmentInfo, GresidueInfo, GatomInfo, Ggrid, Gff, moveFragType, GTempFrag, GTempInfo, d_rng_states);
            
            cudaMemcpy(TempInfo, GTempInfo, sizeof(Atom)*nBlock, cudaMemcpyDeviceToHost);


            const float *period = info->cryst;

            const int waterNum = fragmentInfo[info->fragTypeNum - 1].totalNum;


            // use number of waters to determin nbar value; this allows to take into account excluded volume by solutes

            const float nbar = waterNum / 55.0 * fragmentInfo[moveFragType].conc;


            const float B = beta * fragmentInfo[moveFragType].muex + log(nbar);
            std::unordered_set<unsigned int> conf_index_unused;
            std::unordered_set<unsigned int> conf_index_used;

            std::vector<double> conf_p;

            conf_p.resize(nBlock);

            int conf_index;

            bool needUpdate = false;

            


            
            for (int i = 0;i < nBlock; i++)
                conf_index_unused.insert(i);


            // int addNum = 0;


            while (conf_index_unused.size()>0){
                // printf("conf_index_unused.size() = %d\n",conf_index_unused.size());

                auto it = *conf_index_unused.begin();
                // printf("it = %d \n",it);

                conf_index_used.clear();

                conf_index_used.insert(it);
                conf_index_unused.erase(it);
                

                double sum_p = 0;
                float energy_min = TempInfo[it].charge;
                for (auto iit = conf_index_unused.begin(); iit != conf_index_unused.end();){
                    if ( distanceP(TempInfo[it].position, TempInfo[*iit].position, period) 
                        <= info->cutoff ){
                        // 	printf("sqrd%d:%f\n",iit,(conf_com[*it][0] - conf_com[iit][0] ) * (conf_com[*it][0] - conf_com[iit][0] )  +
                        // (conf_com[*it][1] - conf_com[iit][1] ) * (conf_com[*it][1] - conf_com[iit][1] )  +
                        // (conf_com[*it][2] - conf_com[iit][2] ) * (conf_com[*it][2] - conf_com[iit][2] ));

                            if (TempInfo[*iit].charge < energy_min)
                                energy_min = TempInfo[*iit].charge;

                            
                            conf_index_used.insert(*iit);
                            iit = conf_index_unused.erase(iit);

                            
                    }
                    else
                        iit++;
                }

                // printf("conf_index_unused.size() = %d\n",conf_index_unused.size());
                // printf("conf_index_used.size() = %d\n",conf_index_used.size());


                for (auto iit: conf_index_used) {
                    conf_p[iit] = exp(- beta * (TempInfo[iit].charge - energy_min ));
                    sum_p += conf_p[iit];
                }

                
                if (sum_p == 0) {
                    conf_index = it;
                } else {
                    
                    // conf_p[it] = conf_p[it] / sum_p;
                    double conf_p_sum = 0;
                    for (auto iit : conf_index_used)
                    {
                        conf_p_sum += conf_p[iit] / sum_p;
                        conf_p[iit] = conf_p_sum;
                    }
                    float ran = (float)rand() / (float)RAND_MAX;

                    for (auto iit : conf_index_used){
                        conf_index = iit;

                        if (ran < conf_p[iit] ){				
                            break;
                        }

                    }
                }

                // Fragment frag_tmp = frag;

                // q = conf_q[conf_index];
                // new_com = conf_com[conf_index];
                float energy_new = TempInfo[conf_index].charge;

                // frag_tmp.translate_to(new_com);
                // frag_tmp.rotate(q);

                conf_p[conf_index] =  exp(-beta * (TempInfo[conf_index].charge - energy_min )) / sum_p;



                float fn_tmp = info->cavityFactor / ( conf_index_used.size() * conf_p[conf_index] );


                // printf("fn = %f\n",fn_tmp);


                float diff = energy_new;


                float n =  fragmentInfo[moveFragType].totalNum;
                float p = Min(1, fn_tmp / (n + 1) * exp(B - beta * diff));
                float ran = (float) rand() / (float)RAND_MAX;


                if (ran < p)
                {
                    // success		
                    
                    for (auto iit = conf_index_unused.begin();iit != conf_index_unused.end();){

                    if (distanceP(TempInfo[conf_index].position, TempInfo[*iit].position, period) 
                        <= info->cutoff){
                    // if ( (conf_com[conf_index][0] - conf_com[*iit][0] ) * (conf_com[conf_index][0] - conf_com[*iit][0] )  +
                    // 	(conf_com[conf_index][1] - conf_com[*iit][1] ) * (conf_com[conf_index][1] - conf_com[*iit][1] )  +
                    // 	(conf_com[conf_index][2] - conf_com[*iit][2] ) * (conf_com[conf_index][2] - conf_com[*iit][2] )  
                    // 	<= gcmc->options->energy_cutoff_frag_sqrd){
                        // 	printf("sqrd%d:%f\n",iit,(conf_com[*it][0] - conf_com[iit][0] ) * (conf_com[*it][0] - conf_com[iit][0] )  +
                        // (conf_com[*it][1] - conf_com[iit][1] ) * (conf_com[*it][1] - conf_com[iit][1] )  +
                        // (conf_com[*it][2] - conf_com[iit][2] ) * (conf_com[*it][2] - conf_com[iit][2] ));
                            // printf("Delete %d\n",*iit);

                            iit = conf_index_unused.erase(iit);
                            
                    }
                    else
                        iit ++ ;
                    }

                    // printf("Accept %d energy = %f \n",conf_index, TempInfo[conf_index].charge);


                    // Accept:
                    if (fragmentInfo[moveFragType].totalNum < fragmentInfo[moveFragType].maxNum){
                        
                        TempInfo[conf_index].type = fragmentInfo[moveFragType].totalNum;
                        fragmentInfo[moveFragType].totalNum += 1;

                        // addNum++;
                        
                        needUpdate = true;

                    }
                        



                    // printf("Accept\n");

                    // printf("\nAccept\n");
                    // return_value = true;
                }
                else
                {
                    // printf("Reject %d energy = %f \n",conf_index, TempInfo[conf_index].charge);
                    //printf("\nReject\n");
                    // reject
                    //return false;
                }

                        // printf("conf_index_unused.size() = %d\n",conf_index_unused.size());

                }


 

            if (needUpdate){
                
                cudaMemcpy(GTempInfo, TempInfo, sizeof(Atom)*nBlock, cudaMemcpyHostToDevice);
                GupdateAdd<<<nBlock, numThreadsPerBlock>>>(GfragmentInfo,GresidueInfo, GatomInfo ,GTempFrag ,GTempInfo, moveFragType,fragmentInfo[moveFragType].totalNum);

                


            }

            

            // cudaMemcpy(fragmentInfo, GfragmentInfo, sizeof(AtomArray)*info->fragTypeNum, cudaMemcpyDeviceToHost);
            // int newwater = fragmentInfo[info->fragTypeNum - 1].totalNum;

            // printf("waterNum = %d\n",waterNum);
            // // printf("addNum = %d\n",addNum);
            // printf("newwater = %d\n",newwater);
            // for (int i = 0;i< nBlock;i++){
            //     printf("The energy of the %d fragment is %f\t", i, TempInfo[i].charge);
            //     printf("The type is %d\t", TempInfo[i].type);
            //     printf("The position is %f, %f, %f\n", TempInfo[i].position[0], TempInfo[i].position[1], TempInfo[i].position[2]);
            // }

            //cudaMemcpy(GTempInfo, TempInfo, sizeof(Atom)*maxConf, cudaMemcpyHostToDevice);


            // Gmove_add_check<<<1, numThreadsPerBlock>>>(Ginfo, GfragmentInfo, GresidueInfo, GatomInfo, moveFragType, GTempFrag, GTempInfo, d_rng_states);

            // computeEnergy<<<nBlock, numThreadsPerBlock>>>(Ginfo, GfragmentInfo, GresidueInfo, GatomInfo, Ggrid, Gff, moveFragType);

            
            return needUpdate;
        
        
        }
        

        __global__ void Gmove_del(const InfoStruct *Ginfo, AtomArray *GfragmentInfo, residue *GresidueInfo, 
                Atom *GatomInfo, const float *Ggrid, const float *Gff, const int moveFragType,
                AtomArray *GTempFrag, Atom *GTempInfo, curandState *d_rng_states) {
                
            __shared__ InfoStruct SharedInfo;
            __shared__ AtomArray SharedFragmentInfo;

            int threadId = numThreadsPerBlock * blockIdx.x + threadIdx.x;
            curandState *rng_states = &d_rng_states[threadId];


            int tid = threadIdx.x;

            if (threadIdx.x == 0) {
                SharedInfo = Ginfo[0];
            }
            
            if (threadIdx.x == 1) {
                SharedFragmentInfo = GfragmentInfo[moveFragType];
                SharedFragmentInfo.startRes = GTempInfo[blockIdx.x].type + GfragmentInfo[moveFragType].startRes;
                GTempInfo[blockIdx.x].position[0] = GresidueInfo[SharedFragmentInfo.startRes].position[0];
                GTempInfo[blockIdx.x].position[1] = GresidueInfo[SharedFragmentInfo.startRes].position[1];
                GTempInfo[blockIdx.x].position[2] = GresidueInfo[SharedFragmentInfo.startRes].position[2];

                

                // if (blockIdx.x == 636){
                //     printf("The center of the %d fragment is %8.3f%8.3f%8.3f\n", blockIdx.x, GTempInfo[blockIdx.x].position[0], GTempInfo[blockIdx.x].position[1], GTempInfo[blockIdx.x].position[2]);
                //     printf("The energy of the fragment is %f\n", GTempInfo[blockIdx.x].charge);
                //     printf("The type of the fragment is %d\n", GTempInfo[blockIdx.x].type);
                //     for (int i = 0; i < SharedFragmentInfo.num_atoms; i++){
                //         float total = SharedFragmentInfo.atoms[i].position[0] + SharedFragmentInfo.atoms[i].position[1] + SharedFragmentInfo.atoms[i].position[2] + 1;
                //         printf("%8.3f%8.3f%8.3f %8.3f\n",  SharedFragmentInfo.atoms[i].position[0], SharedFragmentInfo.atoms[i].position[1], SharedFragmentInfo.atoms[i].position[2], total);
                //     }
                // }
            }

            __syncthreads();

            for (int i=tid;i<SharedFragmentInfo.num_atoms;i+=numThreadsPerBlock){
                SharedFragmentInfo.atoms[i].position[0] = GatomInfo[GresidueInfo[SharedFragmentInfo.startRes].atomStart + i].position[0];
                SharedFragmentInfo.atoms[i].position[1] = GatomInfo[GresidueInfo[SharedFragmentInfo.startRes].atomStart + i].position[1];
                SharedFragmentInfo.atoms[i].position[2] = GatomInfo[GresidueInfo[SharedFragmentInfo.startRes].atomStart + i].position[2];
                // SharedFragmentInfo.atoms[i].charge = GatomInfo[GresidueInfo[SharedFragmentInfo.startRes].atomStart + i].charge;
                // SharedFragmentInfo.atoms[i].type = GatomInfo[GresidueInfo[SharedFragmentInfo.startRes].atomStart + i].type;
            }
            
            // randomFragment(SharedInfo, SharedFragmentInfo, &GTempInfo[blockIdx.x], Ggrid, rng_states);

            __syncthreads();

            // if (tid == 0 && blockIdx.x == 636){
            //     printf("The center of the %d fragment is %8.3f%8.3f%8.3f\n", blockIdx.x, GTempInfo[blockIdx.x].position[0], GTempInfo[blockIdx.x].position[1], GTempInfo[blockIdx.x].position[2]);
            //     printf("The energy of the fragment is %f\n", GTempInfo[blockIdx.x].charge);
            //     printf("The type of the fragment is %d\n", GTempInfo[blockIdx.x].type);
            //     for (int i = 0; i < SharedFragmentInfo.num_atoms; i++){
            //         printf("%8.3f%8.3f%8.3f\n", SharedFragmentInfo.atoms[i].position[0], SharedFragmentInfo.atoms[i].position[1], SharedFragmentInfo.atoms[i].position[2]);
            //     }
            // }

            if (tid == 0)
                GTempFrag[blockIdx.x] = SharedFragmentInfo;

            // if (SharedInfo.PME == 0)
                calcEnergy(SharedInfo, SharedFragmentInfo, GfragmentInfo, GresidueInfo, GatomInfo, Gff, &GTempInfo[blockIdx.x]);



                

        }

        

        
        __global__ void GupdateDel(AtomArray *GfragmentInfo, residue *GresidueInfo, 
                        Atom *GatomInfo,
                        AtomArray *GTempFrag, Atom *GTempInfo, const int moveFragType, const int totalNum, const int conf_index) {
            
            int tid = threadIdx.x;
            // int bid = blockIdx.x;

            __shared__ int bidStartRes;
            __shared__ int bidStartAtom;
            __shared__ int bidAtomNum;

            __shared__ int bidEndRes;
            __shared__ int bidEndAtom;

            if (tid == 0){
                GfragmentInfo[moveFragType].totalNum = totalNum;

            }

            if (totalNum == 0)
                return;

            if (tid == 1){

                bidStartRes = GfragmentInfo[moveFragType].startRes + GTempInfo[conf_index].type;
                bidAtomNum = GresidueInfo[bidStartRes].atomNum;
                bidStartAtom = GresidueInfo[bidStartRes].atomStart;

                bidEndRes = GfragmentInfo[moveFragType].startRes + totalNum;
                bidEndAtom = GresidueInfo[bidEndRes].atomStart;

            }
            __syncthreads();



            if (tid == 0){
                
                GresidueInfo[bidStartRes].position[0] = GresidueInfo[bidEndRes].position[0];
                GresidueInfo[bidStartRes].position[1] = GresidueInfo[bidEndRes].position[1];
                GresidueInfo[bidStartRes].position[2] = GresidueInfo[bidEndRes].position[2];
            }
            // __syncthreads();
            for (int i = tid;i<bidAtomNum;i+=numThreadsPerBlock){
                GatomInfo[bidStartAtom + i].position[0] = GatomInfo[bidEndAtom + i].position[0];
                GatomInfo[bidStartAtom + i].position[1] = GatomInfo[bidEndAtom + i].position[1];
                GatomInfo[bidStartAtom + i].position[2] = GatomInfo[bidEndAtom + i].position[2];
            }

            // printf("The energy of the fragment %d is %f\n", GTempInfo[conf_index].type, GTempInfo[conf_index].charge);

        }



        bool move_del(const InfoStruct *info, InfoStruct *Ginfo, AtomArray *fragmentInfo, AtomArray *GfragmentInfo, residue *GresidueInfo, Atom *GatomInfo, const float *Ggrid, const float *Gff,
                    const int moveFragType, AtomArray *GTempFrag, Atom *TempInfo, Atom *GTempInfo, curandState *d_rng_states){

            
            const int nBlock = min(fragmentInfo[moveFragType].confBias, fragmentInfo[moveFragType].totalNum);

            if (nBlock == 0){
                return false;
            }

            
            std::vector<int> nums(fragmentInfo[moveFragType].totalNum);
            for (int i = 0; i < fragmentInfo[moveFragType].totalNum; ++i) {
                nums[i] = i;
            }

            for (int i = 0; i < nBlock; ++i) {
                int j = i + rand() % (fragmentInfo[moveFragType].totalNum - i);

                std::swap(nums[i], nums[j]);

            }

            for (int i=0;i < nBlock;i++){
                TempInfo[i].type = nums[i];
            }

            cudaMemcpy(GTempInfo, TempInfo, sizeof(Atom)*nBlock, cudaMemcpyHostToDevice);
            

            // printf("The nbar value is %f\n", nbar);
            // printf("The B value is %f\n", B);

            // printf("The size of a AtomArray is %d\n", sizeof(AtomArray));

            Gmove_del<<<nBlock, numThreadsPerBlock>>>(Ginfo, GfragmentInfo, GresidueInfo, GatomInfo, Ggrid, Gff, moveFragType, GTempFrag, GTempInfo, d_rng_states);

            cudaMemcpy(TempInfo, GTempInfo, sizeof(Atom)*nBlock, cudaMemcpyDeviceToHost);

            std::vector<double> conf_p(nBlock);
            int conf_index = 0;
            float energy_max = TempInfo[0].charge; 
            double sum_p = 0;

            for (int i=1;i< nBlock;i++){
                if (TempInfo[i].charge > energy_max){
                    energy_max = TempInfo[i].charge;
                }
            }

            for (int i=0;i< nBlock;i++){
                conf_p[i] = exp(beta * (TempInfo[i].charge - energy_max )); // reverse the sign of the energy
                sum_p += conf_p[i];
            }

            if (sum_p == 0) {
                conf_index = 0;
            } else {
                
                // conf_p[it] = conf_p[it] / sum_p;
                double conf_p_sum = 0;
                for (int i=0;i< nBlock;i++)
                {
                    conf_p_sum += conf_p[i] / sum_p;
                    conf_p[i] = conf_p_sum;
                }
                float ran = (float)rand() / (float)RAND_MAX;

                for (int i=0;i< nBlock;i++){
                    conf_index = i;

                    if (ran < conf_p[i] ){				
                        break;
                    }

                }
            }

            // printf("The energy of the fragment %d, %d is %f\n", conf_index, TempInfo[conf_index].type ,TempInfo[conf_index].charge);
            // printf("The max energy is %f\n", energy_max);


            // for (int i=0;i < nBlock;i++){
            //     printf("The energy of the fragment %d, %d is %f\t", i, TempInfo[i].type ,TempInfo[i].charge);
            //     printf("The position of the fragment is %f %f %f\n", TempInfo[i].position[0], TempInfo[i].position[1], TempInfo[i].position[2]);
            // }


            const int waterNum = fragmentInfo[info->fragTypeNum - 1].totalNum;


            // use number of waters to determin nbar value; this allows to take into account excluded volume by solutes

            const float nbar = waterNum / 55.0 * fragmentInfo[moveFragType].conc;


            const float B = beta * fragmentInfo[moveFragType].muex + log(nbar);
                    
            float diff = -TempInfo[conf_index].charge;
            
            conf_p[conf_index] =  exp(beta * (TempInfo[conf_index].charge - energy_max )) / sum_p;

            float fn = 1;
 

            fn /= nBlock * conf_p[conf_index];

            float n =  fragmentInfo[moveFragType].totalNum;


            float p = Min(1, n / fn * exp(-B - beta * diff));
            // cout << p << endl;
            float ran = rand() / (float)RAND_MAX;

            // printf("B = %f, nbar = %f, diff = %f, p = %f, ran = %f, fn * exp(-B - beta * diff) = %f\n", B, nbar, diff, p, ran, fn * exp(-B - beta * diff));

            if (ran < p)
            {
                // success
                // printf("Delete!\n");
                fragmentInfo[moveFragType].totalNum -= 1;

                GupdateDel<<<1, numThreadsPerBlock>>>(GfragmentInfo,GresidueInfo, GatomInfo ,GTempFrag ,GTempInfo, moveFragType,fragmentInfo[moveFragType].totalNum, conf_index);
                
                return true;

            }
            else
            {
                // reject
                return false;
            }

            
            
        }






        






    } 