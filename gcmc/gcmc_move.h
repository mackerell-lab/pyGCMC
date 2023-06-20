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
            
            // for (int i=0;i<2000;i++)
            
            // rotate_atoms(SharedFragmentInfo.atoms, SharedFragmentInfo.num_atoms, randomThi, randomPhi);
            rotate_atoms_shared(SharedFragmentInfo.atoms, SharedFragmentInfo.num_atoms, randomThi, randomPhi);

            __syncthreads();

            if (tid < SharedFragmentInfo.num_atoms){
                SharedFragmentInfo.atoms[tid].position[0] += randomR[0];
                SharedFragmentInfo.atoms[tid].position[1] += randomR[1];
                SharedFragmentInfo.atoms[tid].position[2] += randomR[2];
            }

            if (tid < 3){
                GTempInfo->position[tid] = randomR[tid];
            }





        }

        __device__ inline float fast_round(const float a) {
            return a >= 0 ? (int)(a + 0.5f) : (int)(a - 0.5f);
        }

        __device__ inline float distanceP(const float x[3], const float y[3], const float period[3]){
            float dx = x[0] - y[0];
            float dy = x[1] - y[1];
            float dz = x[2] - y[2];

            dx -= fast_round(dx / period[0]) * period[0];
            dy -= fast_round(dy / period[1]) * period[1];
            dz -= fast_round(dz / period[2]) * period[2];

            return sqrtf(dx * dx + dy * dy + dz * dz);
        }

        __device__ inline void calcProtEnergy(const InfoStruct SharedInfo, AtomArray &SharedFragmentInfo , AtomArray *GfragmentInfo,
                                    residue *GresidueInfo, Atom *GatomInfo, const float *Gff, Atom *GTempInfo, float *sh_energy){
            
            int tid = threadIdx.x;
            __shared__ int maxResidueNum;
            if (tid == 0)
                maxResidueNum = GfragmentInfo->startRes;
            
            
            __syncthreads();

            for (int resi = tid;resi < maxResidueNum; resi+= blockDim.x){

                // if (distanceP(GTempInfo->position, GresidueInfo[resi].position, SharedInfo.cryst) > SharedInfo.cutoff) {
                if (distanceP(GTempInfo->position, GresidueInfo[resi].position, SharedInfo.cryst) > SharedInfo.cutoff) {
                    continue;
                }

                // int resiStart = GresidueInfo[resi].atomStart;
                // int resiEnd = GresidueInfo[resi].atomStart + GresidueInfo[resi].atomNum;
                // float resiEnergy = 0;
                // for (int atomi = resiStart; atomi < resiEnd; atomi++){
                //     int atomType = GatomInfo[atomi].type;
                //     float atomCharge = GatomInfo[atomi].charge;
                //     float atomEnergy = 0;
                //     for (int atomj = 0; atomj < SharedFragmentInfo.num_atoms; atomj++){
                //         float distance = sqrt((GatomInfo[atomi].position[0] - SharedFragmentInfo.atoms[atomj].position[0]) * (GatomInfo[atomi].position[0] - SharedFragmentInfo.atoms[atomj].position[0]) +
                //                             (GatomInfo[atomi].position[1] - SharedFragmentInfo.atoms[atomj].position[1]) * (GatomInfo[atomi].position[1] - SharedFragmentInfo.atoms[atomj].position[1]) +
                //                             (GatomInfo[atomi].position[2] - SharedFragmentInfo.atoms[atomj].position[2]) * (GatomInfo[atomi].position[2] - SharedFragmentInfo.atoms[atomj].position[2]));
                //         float energy = Gff[atomType * SharedInfo.atomTypeNum + SharedFragmentInfo.atoms[atomj].type] * atomCharge * SharedFragmentInfo.atoms[atomj].charge / distance;
                //         atomEnergy += energy;
                //     }
                //     resiEnergy += atomEnergy;
                // }
                // sh_energy[tid] += resiEnergy;
            }


        }
        

        __device__ void calcEnergy(const InfoStruct SharedInfo, AtomArray &SharedFragmentInfo , AtomArray *GfragmentInfo,
                                    residue *GresidueInfo, Atom *GatomInfo, const float *Gff, Atom *GTempInfo){

            __shared__ float sh_energy[128];

            calcProtEnergy(SharedInfo, SharedFragmentInfo, GfragmentInfo, GresidueInfo, GatomInfo, Gff, GTempInfo, sh_energy);
            


        }

        __global__ void Gmove_add(const InfoStruct *Ginfo, AtomArray *GfragmentInfo, residue *GresidueInfo, 
                        Atom *GatomInfo, const float *Ggrid, const float *Gff, const int moveFragType,
                        AtomArray *GTempFrag, Atom *GTempInfo, curandState *d_rng_states) {
                    
            __shared__ InfoStruct SharedInfo;
            __shared__ AtomArray SharedFragmentInfo;

            int threadId = blockDim.x * blockIdx.x + threadIdx.x;
            curandState *rng_states = &d_rng_states[threadId];


            

            if (threadIdx.x == 0) {
                SharedInfo = Ginfo[0];
                SharedFragmentInfo = GfragmentInfo[moveFragType];
            }

            __syncthreads();

            randomFragment(SharedInfo, SharedFragmentInfo, &GTempInfo[blockIdx.x], Ggrid, rng_states);

            __syncthreads();

            calcEnergy(SharedInfo, SharedFragmentInfo, GfragmentInfo, GresidueInfo, GatomInfo, Gff, &GTempInfo[blockIdx.x]);

            __syncthreads();

            // if (threadIdx.x == 0 && blockIdx.x == 0)

            //     printf("The position of the first atom is %f, %f, %f\n", SharedFragmentInfo.atoms[0].position[0], SharedFragmentInfo.atoms[0].position[1], SharedFragmentInfo.atoms[0].position[2]);







        }


        bool move_add(const InfoStruct *Ginfo, AtomArray *fragmentInfo, AtomArray *GfragmentInfo, residue *GresidueInfo, Atom *GatomInfo, const float *Ggrid, const float *Gff,
                    const int moveFragType, AtomArray *GTempFrag, Atom *TempInfo, Atom *GTempInfo, curandState *d_rng_states){

            
            const int nBlock = fragmentInfo[moveFragType].confBias;

            // printf("The size of a AtomArray is %d\n", sizeof(AtomArray));

            Gmove_add<<<nBlock, numThreadsPerBlock>>>(Ginfo, GfragmentInfo, GresidueInfo, GatomInfo, Ggrid, Gff, moveFragType, GTempFrag, GTempInfo, d_rng_states);

            // computeEnergy<<<nBlock, numThreadsPerBlock>>>(Ginfo, GfragmentInfo, GresidueInfo, GatomInfo, Ggrid, Gff, moveFragType);

            
        
        
        
        }
        
        
    } 