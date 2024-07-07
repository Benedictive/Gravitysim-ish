
#include <iostream>
#include <math.h>
#include <vector>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "gravitysim/kernel.hpp"

// function to add the elements of two arrays
__global__ void add(int n, float* x, float* y)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = index; i < n; i += stride)
        y[i] = x[i] + y[i];
}

int cudaAdd(void)
{
    int N = 1 << 20; // 1M elements

    float* x, * y;
    cudaMallocManaged(&x, N * sizeof(float));
    cudaMallocManaged(&y, N * sizeof(float));

    // initialize x and y arrays on the host
    for (int i = 0; i < N; i++) {
        x[i] = 1.0f;
        y[i] = 2.0f;
    }

    int blockSize = 256;
    int numBlocks = (N + blockSize - 1) / blockSize;
    // Run kernel on 1M elements on the CPU
    add <<<numBlocks, blockSize>>> (N, x, y);

    // Wait for GPU to finish before accessing on host
    cudaDeviceSynchronize();

    // Check for errors (all values should be 3.0f)
    float maxError = 0.0f;
    for (int i = 0; i < N; i++)
        maxError = fmax(maxError, fabs(y[i] - 3.0f));
    std::cout << "Max error: " << maxError << std::endl;

    // Free memory
    cudaFree(x);
    cudaFree(y);

    return 0;
}

__global__ void runGravity(size_t objCount, Spaceobject* spaceobjects, Spaceobject* spaceobjectsStaticCopy, int* collisionList, double simDeltaTime) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < objCount; i+=stride) {
        Spaceobject& spobject = spaceobjects[i];
        for (size_t j = 0; j < objCount; j++) {
            if (i == j) continue;
            Spaceobject& otherSpobject = spaceobjectsStaticCopy[j];

            double xDist = otherSpobject.x - spobject.x;
            double yDist = otherSpobject.y - spobject.y;
            double distance = sqrt((xDist*xDist) + (yDist*yDist));
            double minDistance = spobject.radius + otherSpobject.radius;

            if (distance < minDistance) {
                // if self is not already designated for removal, designate other object for removal
                collisionList[i] = j;
                spobject.speedX = (1)*(((spobject.mass*spobject.speedX)+(otherSpobject.mass*otherSpobject.speedX))/(spobject.mass+otherSpobject.mass));
                spobject.speedY = (1)*(((spobject.mass*spobject.speedY)+(otherSpobject.mass*otherSpobject.speedY))/(spobject.mass+otherSpobject.mass));
                spobject.mass += otherSpobject.mass;
            } else {
                double unitDirX = (otherSpobject.x - spobject.x) / distance;
                double unitDirY = (otherSpobject.y - spobject.y) / distance;
                const double bigG = 6.674e-11;
                double force = bigG * ((spobject.mass*otherSpobject.mass)/(distance*distance));

                spobject.speedX += unitDirX * ((force / spobject.mass) * simDeltaTime);
                spobject.speedY += unitDirY * ((force / spobject.mass) * simDeltaTime);
            }
        }
        spobject.x += spobject.speedX * simDeltaTime;
        spobject.y += spobject.speedY * simDeltaTime;            
    }
}

void setupGravityMemory(size_t objectCount, Spaceobject** spaceobjects, Spaceobject** spaceobjectsStaticCopy, int** objectcollisions) {
    cudaMallocManaged(spaceobjects, objectCount * sizeof(Spaceobject));
    cudaMallocManaged(spaceobjectsStaticCopy, objectCount * sizeof(Spaceobject));
    cudaMallocManaged(objectcollisions, objectCount * sizeof(int));
}

void freeGravityMemory(Spaceobject* spaceobjects, Spaceobject* spaceobjectsStaticCopy, int* objectcollisions) {
    cudaFree(spaceobjects);
    cudaFree(spaceobjectsStaticCopy);
    cudaFree(objectcollisions);
}

void cudaGravity(std::vector<Spaceobject>& spaceobjectsRef, std::vector<int>& collisions, double simDeltaTime, Spaceobject* spaceobjects, Spaceobject* spaceobjectsStaticCopy, int* objectcollisions) {
    for (size_t i = 0; i < spaceobjectsRef.size(); i++) {
        spaceobjects[i] = spaceobjectsRef[i];
        spaceobjectsStaticCopy[i] = spaceobjectsRef[i];
        objectcollisions[i] = -1;
    }

    unsigned int blockSize = 256;
    unsigned int numBlocks = (static_cast<unsigned int>(spaceobjectsRef.size()) + blockSize - 1) / blockSize;

    runGravity <<<numBlocks, blockSize>>> (spaceobjectsRef.size(), spaceobjects, spaceobjectsStaticCopy, objectcollisions, simDeltaTime);

    // Wait for GPU to finish before accessing on host
    cudaDeviceSynchronize();

    for (size_t i = 0; i < spaceobjectsRef.size(); i++) {
        spaceobjectsRef[i] = spaceobjects[i];
        collisions[i] = objectcollisions[i];
    }
}