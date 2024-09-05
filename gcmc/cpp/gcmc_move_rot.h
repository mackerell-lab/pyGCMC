/*
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved
        Mingtian Zhao, Alexander D. MacKerell Jr.
    E-mail:
        zhaomt@outerbanks.umaryland.edu
        alex@outerbanks.umaryland.edu
*/

#include "gcmc.h"
#include "gcmc_energy.h"

extern "C" {

// Kernel function to perform rotation move on GPU
__global__ void Gmove_rot(const InfoStruct *Ginfo, AtomArray *GfragmentInfo, Residue *GresidueInfo, 
                          Atom *GatomInfo, const float *Ggrid, const float *Gff, const int moveFragType,
                          AtomArray *GTempFrag, Atom *GTempInfo, curandState *d_rng_states) {
    __shared__ InfoStruct SharedInfo;
    __shared__ AtomArray SharedFragmentInfo;
    __shared__ float energyTemp;
    __shared__ float randomR[3];
    __shared__ float randomThi[3];
    __shared__ float randomPhi;
    __shared__ float center[3];

    int threadId = numThreadsPerBlock * blockIdx.x + threadIdx.x;
    curandState *rng_states = &d_rng_states[threadId];
    int tid = threadIdx.x;

    // Initialize shared memory
    if (threadIdx.x == 0) {
        SharedInfo = Ginfo[0];
    }
    
    if (threadIdx.x == 1) {
        SharedFragmentInfo = GfragmentInfo[moveFragType];
        SharedFragmentInfo.startRes = GTempInfo[blockIdx.x].type + GfragmentInfo[moveFragType].startRes;
        GTempInfo[blockIdx.x].position[0] = GresidueInfo[SharedFragmentInfo.startRes].position[0];
        GTempInfo[blockIdx.x].position[1] = GresidueInfo[SharedFragmentInfo.startRes].position[1];
        GTempInfo[blockIdx.x].position[2] = GresidueInfo[SharedFragmentInfo.startRes].position[2];
    }

    __syncthreads();

    // Copy atom positions to shared memory
    for (int i = tid; i < SharedFragmentInfo.num_atoms; i += numThreadsPerBlock) {
        SharedFragmentInfo.atoms[i].position[0] = GatomInfo[GresidueInfo[SharedFragmentInfo.startRes].atomStart + i].position[0];
        SharedFragmentInfo.atoms[i].position[1] = GatomInfo[GresidueInfo[SharedFragmentInfo.startRes].atomStart + i].position[1];
        SharedFragmentInfo.atoms[i].position[2] = GatomInfo[GresidueInfo[SharedFragmentInfo.startRes].atomStart + i].position[2];
    }
    
    __syncthreads();

    calcEnergy(SharedInfo, SharedFragmentInfo, GfragmentInfo, GresidueInfo, GatomInfo, Gff, &GTempInfo[blockIdx.x]);

    __syncthreads();

    if (tid == 3) {
        energyTemp = GTempInfo[blockIdx.x].charge;
    }

    // Generate random rotation parameters
    if (tid < 3) {
        randomThi[tid] = curand_uniform(rng_states);
        center[tid] = 0;
    }
    if (tid == 3) {
        randomPhi = curand_uniform(rng_states) * 2 * PI;
    }

    __syncthreads();

    // Calculate center of mass
    if (tid < 3) {
        for (int i = 0; i < SharedFragmentInfo.num_atoms; i++) {
            center[tid] += SharedFragmentInfo.atoms[i].position[tid] / SharedFragmentInfo.num_atoms;
        }
    }

    __syncthreads();

    // Move atoms to origin
    for (int i = tid; i < SharedFragmentInfo.num_atoms; i += numThreadsPerBlock) {
        SharedFragmentInfo.atoms[i].position[0] -= center[0];
        SharedFragmentInfo.atoms[i].position[1] -= center[1];
        SharedFragmentInfo.atoms[i].position[2] -= center[2];
    }
    __syncthreads();

    // Perform rotation
    rotate_atoms_shared(SharedFragmentInfo.atoms, SharedFragmentInfo.num_atoms, randomThi, randomPhi);

    __syncthreads();

    // Move atoms back to original position
    for (int i = tid; i < SharedFragmentInfo.num_atoms; i += numThreadsPerBlock) {
        SharedFragmentInfo.atoms[i].position[0] += center[0];
        SharedFragmentInfo.atoms[i].position[1] += center[1];
        SharedFragmentInfo.atoms[i].position[2] += center[2];
    }

    __syncthreads();

    // Calculate new energy
    calcEnergy(SharedInfo, SharedFragmentInfo, GfragmentInfo, GresidueInfo, GatomInfo, Gff, &GTempInfo[blockIdx.x]);

    __syncthreads();

    if (tid == 3) {
        GTempInfo[blockIdx.x].charge = GTempInfo[blockIdx.x].charge - energyTemp;
    }
    
    if (tid == 0)
        GTempFrag[blockIdx.x] = SharedFragmentInfo;
}

// Kernel function to update rotated fragments
__global__ void GupdateRot(AtomArray *GfragmentInfo, Residue *GresidueInfo, 
                           Atom *GatomInfo,
                           AtomArray *GTempFrag, Atom *GTempInfo, const int moveFragType) {
    int tid = threadIdx.x;
    int bid = blockIdx.x;

    __shared__ int bidStartRes;
    __shared__ int bidStartAtom;
    __shared__ int bidAtomNum;

    if (GTempInfo[bid].type == -1)
        return;

    if (tid == 1) {
        bidStartRes = GfragmentInfo[moveFragType].startRes + GTempInfo[bid].type;
    }
    __syncthreads();
    if (tid == 0) {
        bidAtomNum = GresidueInfo[bidStartRes].atomNum;
        bidStartAtom = GresidueInfo[bidStartRes].atomStart;
        GresidueInfo[bidStartRes].position[0] = GTempInfo[bid].position[0];
        GresidueInfo[bidStartRes].position[1] = GTempInfo[bid].position[1];
        GresidueInfo[bidStartRes].position[2] = GTempInfo[bid].position[2];
    }
    __syncthreads();
    for (int i = tid; i < bidAtomNum; i += numThreadsPerBlock) {
        GatomInfo[bidStartAtom + i].position[0] = GTempFrag[bid].atoms[i].position[0];
        GatomInfo[bidStartAtom + i].position[1] = GTempFrag[bid].atoms[i].position[1];
        GatomInfo[bidStartAtom + i].position[2] = GTempFrag[bid].atoms[i].position[2];
    }
}

// Function to perform rotation move
bool move_rot(const InfoStruct *infoHost, InfoStruct *infoDevice, AtomArray *fragmentInfoHost, AtomArray *fragmentInfoDevice, 
              Residue *residueInfoDevice, Atom *atomInfoDevice, const float *gridDevice, const float *ffDevice,
              const int moveFragType, AtomArray *tempFragDevice, Atom *tempInfoHost, Atom *tempInfoDevice, curandState *rngStatesDevice) {

    const int numBlocks = min(fragmentInfoHost[moveFragType].confBias, fragmentInfoHost[moveFragType].totalNum);

    if (numBlocks == 0) {
        return false;
    }

    // Initialize and shuffle fragment indices
    std::vector<int> nums(fragmentInfoHost[moveFragType].totalNum);
    for (int i = 0; i < fragmentInfoHost[moveFragType].totalNum; ++i) {
        nums[i] = i;
    }
    for (int i = 0; i < numBlocks; ++i) {
        int j = i + rand() % (fragmentInfoHost[moveFragType].totalNum - i);
        std::swap(nums[i], nums[j]);
    }

    // Assign shuffled indices to TempInfo
    for (int i = 0; i < numBlocks; i++) {
        tempInfoHost[i].type = nums[i];
    }

    // Copy TempInfo to device memory and perform rotation move
    cudaMemcpy(tempInfoDevice, tempInfoHost, sizeof(Atom)*numBlocks, cudaMemcpyHostToDevice);
    Gmove_rot<<<numBlocks, numThreadsPerBlock>>>(infoDevice, fragmentInfoDevice, residueInfoDevice, atomInfoDevice, 
                                                 gridDevice, ffDevice, moveFragType, tempFragDevice, tempInfoDevice, rngStatesDevice);
    cudaMemcpy(tempInfoHost, tempInfoDevice, sizeof(Atom)*numBlocks, cudaMemcpyDeviceToHost);

    const float *period = infoHost->cryst;
    std::unordered_set<unsigned int> confIndexUnused;
    std::unordered_set<unsigned int> confIndexUsed;
    std::vector<double> confProbabilities(numBlocks);
    int confIndex;
    bool needUpdate = false;

    // Main loop for configuration selection
    for (int i = 0; i < numBlocks; i++)
        confIndexUnused.insert(i);
    while (!confIndexUnused.empty()) {
        auto it = *confIndexUnused.begin();
        confIndexUsed.clear();
        confIndexUsed.insert(it);
        confIndexUnused.erase(it);

        double sumP = 0;
        float energyMin = tempInfoHost[it].charge;

        // Group configurations within cutoff distance
        for (auto iit = confIndexUnused.begin(); iit != confIndexUnused.end();) {
            if (distanceP(tempInfoHost[it].position, tempInfoHost[*iit].position, period) <= infoHost->cutoff) {
                if (tempInfoHost[*iit].charge < energyMin)
                    energyMin = tempInfoHost[*iit].charge;
                confIndexUsed.insert(*iit);
                iit = confIndexUnused.erase(iit);
            } else {
                ++iit;
            }
        }

        // Calculate probabilities
        for (auto iit : confIndexUsed) {
            confProbabilities[iit] = exp(-beta * (tempInfoHost[iit].charge - energyMin));
            sumP += confProbabilities[iit];
        }

        // Select configuration
        if (sumP == 0) {
            confIndex = it;
        } else {
            double confPSum = 0;
            for (auto iit : confIndexUsed) {
                confPSum += confProbabilities[iit] / sumP;
                confProbabilities[iit] = confPSum;
            }
            float ran = (float)rand() / (float)RAND_MAX;
            for (auto iit : confIndexUsed) {
                confIndex = iit;
                if (ran < confProbabilities[iit]) {
                    break;
                }
            }
        }

        // Calculate acceptance probability
        float energyNew = tempInfoHost[confIndex].charge;
        confProbabilities[confIndex] = exp(-beta * (tempInfoHost[confIndex].charge - energyMin)) / sumP;
        float fnTmp = 1 / (confIndexUsed.size() * confProbabilities[confIndex]);
        float diff = energyNew;
        float p = Min(1, fnTmp * exp(-beta * diff));

        float ran = (float)rand() / (float)RAND_MAX;
        int tempInfoHostType = tempInfoHost[confIndex].type;

        for (auto iit : confIndexUsed) {
            tempInfoHost[iit].type = -1;
        }

        // Accept or reject move
        if (ran < p) {
            for (auto iit = confIndexUnused.begin(); iit != confIndexUnused.end();) {
                if (distanceP(tempInfoHost[confIndex].position, tempInfoHost[*iit].position, period) <= infoHost->cutoff) {
                    tempInfoHost[*iit].type = -1;
                    iit = confIndexUnused.erase(iit);
                } else {
                    ++iit;
                }
            }

            tempInfoHost[confIndex].type = tempInfoHostType;
            needUpdate = true;
        }
    }

    // Update if move is accepted
    if (needUpdate) {
        cudaMemcpy(tempInfoDevice, tempInfoHost, sizeof(Atom)*numBlocks, cudaMemcpyHostToDevice);
        GupdateRot<<<numBlocks, numThreadsPerBlock>>>(fragmentInfoDevice, residueInfoDevice, atomInfoDevice, 
                                                      tempFragDevice, tempInfoDevice, moveFragType);
    }

    return needUpdate;
}

} // extern "C"
