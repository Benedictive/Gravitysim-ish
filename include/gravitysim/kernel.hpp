#pragma once

#include "gravitysim.hpp"

int cudaAdd(void);
void setupGravityMemory(size_t objectCount, Spaceobject** spaceobjects, Spaceobject** spaceobjectsStaticCopy, int** objectcollisions);
void freeGravityMemory(Spaceobject* spaceobjects, Spaceobject* spaceobjectsStaticCopy, int* objectcollisions);
void cudaGravity(std::vector<Spaceobject>& spaceobjectsRef, std::vector<int>& collisions, double simDeltaTime, Spaceobject* spaceobjects, Spaceobject* spaceobjectsStaticCopy, int* objectcollisions);