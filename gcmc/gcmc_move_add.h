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

// CUDA kernel for adding a fragment
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
    }

    __syncthreads();

    randomFragment(SharedInfo, SharedFragmentInfo, &GTempInfo[blockIdx.x], Ggrid, rng_states);

    __syncthreads();

    if (tid == 0)
        GTempFrag[blockIdx.x] = SharedFragmentInfo;

    calcEnergy(SharedInfo, SharedFragmentInfo, GfragmentInfo, GresidueInfo, GatomInfo, Gff, &GTempInfo[blockIdx.x]);
}

// CUDA kernel for updating after adding a fragment
__global__ void GupdateAdd(AtomArray *GfragmentInfo, residue *GresidueInfo, 
                           Atom *GatomInfo, AtomArray *GTempFrag, Atom *GTempInfo, 
                           const int moveFragType, const int totalNum) {
    int tid = threadIdx.x;
    int bid = blockIdx.x;

    __shared__ int bidStartRes;
    __shared__ int bidStartAtom;
    __shared__ int bidAtomNum;

    if (tid == 0 && bid == 0) {
        GfragmentInfo[moveFragType].totalNum = totalNum;
    }

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
        GatomInfo[bidStartAtom + i].type = GTempFrag[bid].atoms[i].type;
        GatomInfo[bidStartAtom + i].charge = GTempFrag[bid].atoms[i].charge;
    }
}

// Function to perform the move_add operation
bool move_add(const InfoStruct *infoHost, InfoStruct *infoDevice, AtomArray *fragmentInfoHost, 
              AtomArray *fragmentInfoDevice, residue *residueInfoDevice, Atom *atomInfoDevice, 
              const float *gridDevice, const float *ffDevice, const int moveFragType, 
              AtomArray *tempFragDevice, Atom *tempInfoHost, Atom *tempInfoDevice, 
              curandState *rngStatesDevice) {

    if (fragmentInfoHost[moveFragType].totalNum == fragmentInfoHost[moveFragType].maxNum)
        return false;

    const int numBlocks = fragmentInfoHost[moveFragType].confBias;

    Gmove_add<<<numBlocks, numThreadsPerBlock>>>(infoDevice, fragmentInfoDevice, residueInfoDevice, 
                                                 atomInfoDevice, gridDevice, ffDevice, moveFragType, 
                                                 tempFragDevice, tempInfoDevice, rngStatesDevice);

    cudaMemcpy(tempInfoHost, tempInfoDevice, sizeof(Atom)*numBlocks, cudaMemcpyDeviceToHost);

    const float *period = infoHost->cryst;
    const float numBars = period[0] * period[1] * period[2] * fragmentInfoHost[moveFragType].conc * MOLES_TO_MOLECULES;
    const float B = beta * fragmentInfoHost[moveFragType].muex + log(numBars);

    std::unordered_set<unsigned int> confIndexUnused;
    std::unordered_set<unsigned int> confIndexUsed;
    std::vector<double> confProbabilities(numBlocks);

    int confIndex;
    bool needUpdate = false;

    for (int i = 0; i < numBlocks; i++)
        confIndexUnused.insert(i);

    // Main loop for processing configurations
    while (confIndexUnused.size() > 0) {
        auto it = *confIndexUnused.begin();
        confIndexUsed.clear();
        confIndexUsed.insert(it);
        confIndexUnused.erase(it);

        double sumP = 0;
        float energyMin = tempInfoHost[it].charge;

        // Process unused configurations
        for (auto iit = confIndexUnused.begin(); iit != confIndexUnused.end();) {
            if (distanceP(tempInfoHost[it].position, tempInfoHost[*iit].position, period) <= infoHost->cutoff) {
                if (tempInfoHost[*iit].charge < energyMin)
                    energyMin = tempInfoHost[*iit].charge;
                confIndexUsed.insert(*iit);
                iit = confIndexUnused.erase(iit);
            } else {
                iit++;
            }
        }

        // Calculate probabilities
        for (auto iit: confIndexUsed) {
            confProbabilities[iit] = exp(- beta * (tempInfoHost[iit].charge - energyMin));
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

        float energyNew = tempInfoHost[confIndex].charge;
        confProbabilities[confIndex] = exp(-beta * (tempInfoHost[confIndex].charge - energyMin)) / sumP;

        float fnTmp = infoHost->cavityFactor / (confIndexUsed.size() * confProbabilities[confIndex]);
        float diff = energyNew;
        float n = fragmentInfoHost[moveFragType].totalNum;
        float p = Min(1, fnTmp / (n + 1) * exp(B - beta * diff));

        float ran = (float) rand() / (float)RAND_MAX;

        if (ran < p) {
            // Remove unused configurations within cutoff
            for (auto iit = confIndexUnused.begin(); iit != confIndexUnused.end();) {
                if (distanceP(tempInfoHost[confIndex].position, tempInfoHost[*iit].position, period) <= infoHost->cutoff) {
                    iit = confIndexUnused.erase(iit);
                } else {
                    iit++;
                }
            }

            // Update fragment information if possible
            if (fragmentInfoHost[moveFragType].totalNum < fragmentInfoHost[moveFragType].maxNum) {
                tempInfoHost[confIndex].type = fragmentInfoHost[moveFragType].totalNum;
                fragmentInfoHost[moveFragType].totalNum += 1;
                needUpdate = true;
            }
        }
    }

    if (needUpdate) {
        cudaMemcpy(tempInfoDevice, tempInfoHost, sizeof(Atom)*numBlocks, cudaMemcpyHostToDevice);
        GupdateAdd<<<numBlocks, numThreadsPerBlock>>>(fragmentInfoDevice, residueInfoDevice, 
                                                      atomInfoDevice, tempFragDevice, tempInfoDevice, 
                                                      moveFragType, fragmentInfoHost[moveFragType].totalNum);
    }

    return needUpdate;
}

}
