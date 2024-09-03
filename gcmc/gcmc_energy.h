/*
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved
        Mingtian Zhao, Alexander D. MacKerell Jr.
    E-mail:
        zhaomt@outerbanks.umaryland.edu
        alex@outerbanks.umaryland.edu
*/



#include "gcmc.h"


#ifndef GCMC_ENERGY_H_
#define GCMC_ENERGY_H_




extern "C"{

        // Returns the minimum of two float values
        __device__ __host__ inline float Min(float a,float b)
        {
            return !(b<a)?a:b;	
        }

        // Rotates atoms in shared memory
        __device__ inline void rotate_atoms_shared(Atom *atoms, int num_atoms, float axis[3], float angle) {
            __shared__ float sh_axis[3];
            __shared__ float sh_R[3][3];

            int idx = threadIdx.x;

            if (idx < 3) {
                sh_axis[idx] = axis[idx];
            }

            if (idx == 0) {
                float norm = sqrt(sh_axis[0] * sh_axis[0] + sh_axis[1] * sh_axis[1] + sh_axis[2] * sh_axis[2]);
                for (int i = 0; i < 3; i++) sh_axis[i] /= norm;

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

            __syncthreads();

            for (int i = idx; i < num_atoms; i += blockDim.x) {
                float atom_position[3] = {atoms[i].position[0], atoms[i].position[1], atoms[i].position[2]};
                
                for (int j = 0; j < 3; j++) {
                    atoms[i].position[j] = 0;
                    for (int k = 0; k < 3; k++) {
                        atoms[i].position[j] += atom_position[k] * sh_R[j][k];
                    }
                }
            }
        }

        // Generates a random fragment
        __device__ void randomFragment(const InfoStruct &SharedInfo, AtomArray &SharedFragmentInfo, Atom *GTempInfo, const float *Ggrid, curandState *rng_states) {
            int tid = threadIdx.x;

            __shared__ float randomR[3], randomThi[3];
            __shared__ float randomPhi;
            __shared__ int gridN;

            if (tid < 3) {
                randomR[tid] = curand_uniform(rng_states) * SharedInfo.grid_dx;
                randomThi[tid] = curand_uniform(rng_states);
            }
            if (tid == 6) randomPhi = curand_uniform(rng_states) * 2 * PI;
            if (tid == 7) gridN = curand(rng_states) % SharedInfo.totalGridNum;

            __syncthreads();

            if (tid < 3) randomR[tid] += Ggrid[gridN * 3 + tid];

            __syncthreads();
            
            rotate_atoms_shared(SharedFragmentInfo.atoms, SharedFragmentInfo.num_atoms, randomThi, randomPhi);

            __syncthreads();

            for (int i = tid; i < SharedFragmentInfo.num_atoms; i += blockDim.x) {
                for (int j = 0; j < 3; j++) {
                    SharedFragmentInfo.atoms[i].position[j] += randomR[j];
                }
            }

            if (tid < 3) GTempInfo->position[tid] = randomR[tid];
            if (tid == 4) GTempInfo->type = -1;
        }

        // Fast rounding of float to integer
        __device__ __host__ inline float fast_round(const float a) {
            return a >= 0 ? (int)(a + 0.5f) : (int)(a - 0.5f);
        }

        // Calculates distance with periodic boundary conditions
        __device__ __host__ inline float distanceP(const float x[3], const float y[3], const float period[3]){
            float d[3];
            for (int i = 0; i < 3; i++) {
                d[i] = x[i] - y[i];
                d[i] -= fast_round(d[i] / period[i]) * period[i];
            }
            return sqrtf(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        }

        
        // Calculates van der Waals energy
        __device__ inline float calc_vdw_nbfix (float sigma, float epsilon, float dist_sqrd)
        {
            float sigma_sqrd = sigma * sigma * 100;
            float sigma_dist_sqrd = sigma_sqrd / dist_sqrd;
            return 4*epsilon * ( pow(sigma_dist_sqrd, 6) - pow(sigma_dist_sqrd, 3) );
        }

        
        // Calculates electrostatic energy
        __device__ inline float calc_elec (float charge1, float charge2, float dist)
        {
            return 1388.431112 * (charge1 * charge2) / dist;
        }


        // Calculates total energy
        __device__ inline void calcEnergy(const InfoStruct &SharedInfo, AtomArray &SharedFragmentInfo, 
                                          AtomArray *GfragmentInfo, residue *GresidueInfo, Atom *GatomInfo, 
                                          const float *Gff, Atom *GTempInfo) {
            __shared__ float sh_energy[numThreadsPerBlock];
            int tid = threadIdx.x;
            sh_energy[tid] = 0;

            __syncthreads();

            // Calculate protein energy
            int maxResidueNum = GfragmentInfo->startRes;
            for (int resi = tid; resi < maxResidueNum; resi += numThreadsPerBlock) {
                if (distanceP(GTempInfo->position, GresidueInfo[resi].position, SharedInfo.cryst) <= SharedInfo.cutoff) {
                    int resiStart = GresidueInfo[resi].atomStart;
                    int resiEnd = resiStart + GresidueInfo[resi].atomNum;
                    for (int atomi = resiStart; atomi < resiEnd; atomi++) {
                        for (int atomj = 0; atomj < SharedFragmentInfo.num_atoms; atomj++) {
                            float distance = distanceP(GatomInfo[atomi].position, SharedFragmentInfo.atoms[atomj].position, SharedInfo.cryst);
                            int typeij = SharedFragmentInfo.atoms[atomj].type * SharedInfo.ffYNum + GatomInfo[atomi].type;
                            sh_energy[tid] += calc_vdw_nbfix(Gff[typeij * 2], Gff[typeij * 2 + 1], distance * distance);
                            sh_energy[tid] += calc_elec(SharedFragmentInfo.atoms[atomj].charge, GatomInfo[atomi].charge, distance);
                        }
                    }
                }
            }

            __syncthreads();

            // Calculate fragment energy
            for (int fragType = 0; fragType < SharedInfo.fragTypeNum; fragType++) {
                int startResidueNum = GfragmentInfo[fragType].startRes;
                int endResidueNum = startResidueNum + GfragmentInfo[fragType].totalNum;

                for (int resi = tid + startResidueNum; resi < endResidueNum; resi += numThreadsPerBlock) {
                    if (resi != SharedFragmentInfo.startRes && 
                        distanceP(GTempInfo->position, GresidueInfo[resi].position, SharedInfo.cryst) <= SharedInfo.cutoff) {
                        int resiStart = GresidueInfo[resi].atomStart;
                        int resiEnd = resiStart + GresidueInfo[resi].atomNum;
                        for (int atomi = resiStart; atomi < resiEnd; atomi++) {
                            for (int atomj = 0; atomj < SharedFragmentInfo.num_atoms; atomj++) {
                                float distance = distanceP(GatomInfo[atomi].position, SharedFragmentInfo.atoms[atomj].position, SharedInfo.cryst);
                                int typeij = SharedFragmentInfo.atoms[atomj].type * SharedInfo.ffYNum + GatomInfo[atomi].type;
                                sh_energy[tid] += calc_vdw_nbfix(Gff[typeij * 2], Gff[typeij * 2 + 1], distance * distance);
                                sh_energy[tid] += calc_elec(SharedFragmentInfo.atoms[atomj].charge, GatomInfo[atomi].charge, distance);
                            }
                        }
                    }
                }
                __syncthreads();
            }

            // Reduction sum
            for (int s = numThreadsPerBlock / 2; s > 0; s >>= 1) {
                if (tid < s) {
                    sh_energy[tid] += sh_energy[tid + s];
                }
                __syncthreads();
            }

            if (tid == 0) {
                GTempInfo->charge = sh_energy[0];
            }
        }

    }

#endif
