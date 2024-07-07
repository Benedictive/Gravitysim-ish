//#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <deque>
#include <thread>
#include <semaphore>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>

#include "simplecanvas/simple_canvas.hpp"

#include "gravitysim/gravitysim.hpp"
#include "gravitysim/kernel.hpp"

const unsigned int WIDTH = 1920;
const unsigned int HEIGHT = 1080;
const double screenRatio = static_cast<double>(HEIGHT)/static_cast<double>(WIDTH);

double xOffset = 0;
double yOffset = 0;
const double offsetFactor = 1000;
double zoom = 1;

std::vector<unsigned int> imageData;
std::binary_semaphore imageDataLock{1};

bool running = true;

void handleInput(SDL_Event& event) {
    if (event.type == SDL_KEYDOWN) {
        if (event.key.state == SDL_PRESSED) {
            switch (event.key.keysym.sym) {
                case SDLK_w:
                    //move up
                    yOffset += ((offsetFactor*screenRatio)/zoom)/2;
                    break;
                case SDLK_s:
                    //move down
                    yOffset -= ((offsetFactor*screenRatio)/zoom)/2;
                    break;
                case SDLK_a:
                    //move left
                    xOffset += (offsetFactor/zoom)/2;
                    break;
                case SDLK_d:
                    //move right
                    xOffset -= (offsetFactor/zoom)/2;
                    break;
                case SDLK_PLUS:
                    //zoom in
                    zoom *= 1.5;
                    break;
                case SDLK_MINUS:
                    //zoom out
                    zoom /= 1.5;
                    break;
            }
        }
    }
}

double vecMangnitude(double x, double y) {
    return std::sqrt(std::pow(x, 2) + std::pow(y, 2));
}

#define CPU 1

#if CPU
void Gravity(unsigned int height, unsigned int width) {
    const double simSecsPerRealSec = 50000000;
    // how much time must pass per iteration in the simulation to achieve the desired simulationspeed
    double simDeltaTime = 0.1;

    // determines how efficient the collision conserves kinetic energy
    // 0 atm since collisions absorb eachother
    double objEllasticity = 0;
    std::chrono::duration<double> onesecond(1);

    const unsigned int loopsPerSimIter = 1;

    const double bigG = 6.674e-11;

    const unsigned int objectcount = 1000;
    std::vector<Spaceobject> spaceobjects;

    const int worldScale = 1000;
    const int worldWidth = width * worldScale;
    const int worldHeight = height * worldScale;

    std::srand(567890);
    int randomInt = std::rand();

    // prepopulate space with objects
    for (int i = 0; i < objectcount; i++) {
        spaceobjects.push_back(Spaceobject((std::rand() % width) * worldScale * 100, (std::rand() % height) * worldScale * 100));
    }
    Spaceobject newObj(spaceobjects[objectcount/2].x, spaceobjects[objectcount/2].y, 500*200);
    spaceobjects[objectcount/2] = newObj;

    Spaceobject newObj2(spaceobjects[objectcount/3].x, spaceobjects[objectcount/3].y, 500*100);
    spaceobjects[objectcount/3] = newObj2;

    Spaceobject newObj3(spaceobjects[objectcount/4].x, spaceobjects[objectcount/4].y, 500*100);
    spaceobjects[objectcount/4] = newObj3;

    Spaceobject bfp(2e8, 2e8, 2e5);
    bfp.speedX += 0.1;
    bfp.speedY -= 0.03;
    spaceobjects.push_back(bfp);

    std::vector<Spaceobject> spaceobjectsWorkCopy(spaceobjects);

    double prevDistance = 0;
    double prevDDistance = 0;
    while (running) {
        std::copy(spaceobjects.begin(), spaceobjects.end(), spaceobjectsWorkCopy.begin());
        std::vector<size_t> indToRemove;
        auto startTime = std::chrono::high_resolution_clock::now();
        double summedDist = 0;
        for (unsigned int i = 0; i < loopsPerSimIter; i++) {
            for (size_t i = 0; i < spaceobjectsWorkCopy.size(); i++) {
                Spaceobject& spobject = spaceobjectsWorkCopy[i];
                for (size_t j = 0; j < spaceobjects.size(); j++) {
                    if (i == j) continue;
                    Spaceobject& otherSpobject = spaceobjects[j];

                    double distance = vecMangnitude(otherSpobject.x - spobject.x, otherSpobject.y - spobject.y);
                    summedDist += distance;
                    double minDistance = spobject.radius + otherSpobject.radius;

                    if (distance < minDistance) {
                        // if self is not already designated for removal, designate other object for removal
                        if (std::find(indToRemove.begin(), indToRemove.end(), i) == indToRemove.end()) {
                            indToRemove.push_back(j);

                            spobject.speedX = (1+objEllasticity)*(((spobject.mass*spobject.speedX)+(otherSpobject.mass*otherSpobject.speedX))/(spobject.mass+otherSpobject.mass))-(objEllasticity*spobject.speedX);
                            spobject.speedY = (1+objEllasticity)*(((spobject.mass*spobject.speedY)+(otherSpobject.mass*otherSpobject.speedY))/(spobject.mass+otherSpobject.mass))-(objEllasticity*spobject.speedY);
                            spobject.mass += otherSpobject.mass;
                        }
                    } else {
                        double unitDirX = (otherSpobject.x - spobject.x) / distance;
                        double unitDirY = (otherSpobject.y - spobject.y) / distance;
                        double force = bigG * ((spobject.mass*otherSpobject.mass)/(distance*distance));

                        spobject.speedX += unitDirX * ((force / spobject.mass) * simDeltaTime);
                        spobject.speedY += unitDirY * ((force / spobject.mass) * simDeltaTime);
                    }
                }
                spobject.x += spobject.speedX * simDeltaTime;
                spobject.y += spobject.speedY * simDeltaTime;            
            }
        }

        spaceobjects.clear();
        for (size_t idx = 0; idx < spaceobjectsWorkCopy.size(); idx++) {
            if (std::find(indToRemove.begin(), indToRemove.end(), idx) == indToRemove.end()) {
                spaceobjects.push_back(spaceobjectsWorkCopy[idx]);
            }
        }
        spaceobjectsWorkCopy = spaceobjects;

        imageDataLock.acquire();
        std::fill(imageData.begin(), imageData.end(), 0);
        
        for (Spaceobject& spobject : spaceobjects) {
            //draw
            int finY = static_cast<int>(((spobject.y/worldScale)+yOffset)*zoom);
            int finX = static_cast<int>(((spobject.x/worldScale)+xOffset)*zoom);
            if ((finY < 0 || finY > static_cast<int>(height-1)) || (finX < 0 || finX > static_cast<int>(width-1))) continue;
            imageData[finY * width + finX] = 0x00FFFFFF;
        }
        imageDataLock.release();
        
        auto endTime = std::chrono::high_resolution_clock::now();
	    std::chrono::duration<double> executionTime{endTime - startTime};
        simDeltaTime = executionTime.count() * simSecsPerRealSec;

        auto ups = onesecond / executionTime;
        std::cout << std::format("UPS = {0:10.5}", ups) << std::endl;
    }
}

#else

void Gravity(unsigned int height, unsigned int width) {
    const double simSecsPerRealSec = 50000000;
    // how much time must pass per iteration in the simulation to achieve the desired simulationspeed
    double simDeltaTime = 0.1;

    // determines how efficient the collision conserves kinetic energy
    // 0 atm since collisions absorb eachother
    std::chrono::duration<double> onesecond(1);

    const unsigned int objectcount = 2000;
    std::vector<Spaceobject> spaceobjects;

    const int worldScale = 1000;
    const int worldWidth = width * worldScale;
    const int worldHeight = height * worldScale;

    // prepopulate space with objects
    for (int i = 0; i < objectcount; i++) {
        spaceobjects.push_back(Spaceobject((std::rand() % width) * worldScale * 100, (std::rand() % height) * worldScale * 100));
    }
    Spaceobject newObj(spaceobjects[objectcount/2].x, spaceobjects[objectcount/2].y, 500*200);
    spaceobjects[objectcount/2] = newObj;

    Spaceobject newObj2(spaceobjects[objectcount/3].x, spaceobjects[objectcount/3].y, 500*100);
    spaceobjects[objectcount/3] = newObj2;

    Spaceobject newObj3(spaceobjects[objectcount/4].x, spaceobjects[objectcount/4].y, 500*100);
    spaceobjects[objectcount/4] = newObj3;

    Spaceobject bigMass(2e8, 2e8, 2e5);
    bigMass.speedX += 0.1;
    bigMass.speedY -= 0.03;
    spaceobjects.push_back(bigMass);

    Spaceobject* spaceobjectsCuda = nullptr;
    Spaceobject* spaceobjectsStaticCopy = nullptr;
    int* objectcollisions = nullptr;

    setupGravityMemory(spaceobjects.size(), &spaceobjectsCuda, &spaceobjectsStaticCopy, &objectcollisions);

    while (running) {
        std::vector<int> collisionMap(spaceobjects.size());
        std::vector<size_t> indToRemove;
        auto startTime = std::chrono::high_resolution_clock::now();
        
        cudaGravity(spaceobjects, collisionMap, simDeltaTime, spaceobjectsCuda, spaceobjectsStaticCopy, objectcollisions);

        for (size_t i = 0; i < collisionMap.size(); i++) {
            if (collisionMap[i] >= 0 && std::find(indToRemove.begin(), indToRemove.end(), i) == indToRemove.end()) {
                indToRemove.push_back(collisionMap[i]);
            }
        }
        
        std::vector<Spaceobject> spaceobjectsTemp(spaceobjects);
        spaceobjects.clear();
        for (size_t idx = 0; idx < spaceobjectsTemp.size(); idx++) {
            if (std::find(indToRemove.begin(), indToRemove.end(), idx) == indToRemove.end()) {
                spaceobjects.push_back(spaceobjectsTemp[idx]);
            }
        }

        imageDataLock.acquire();
        std::fill(imageData.begin(), imageData.end(), 0);
        
        for (Spaceobject& spobject : spaceobjects) {
            //draw
            int finY = static_cast<int>(((spobject.y/worldScale)+yOffset)*zoom);
            int finX = static_cast<int>(((spobject.x/worldScale)+xOffset)*zoom);
            if ((finY < 0 || finY > static_cast<int>(height-1)) || (finX < 0 || finX > static_cast<int>(width-1))) continue;
            imageData[finY * width + finX] = 0x00FFFFFF;
        }
        imageDataLock.release();
        
        auto endTime = std::chrono::high_resolution_clock::now();
	    std::chrono::duration<double> executionTime{endTime - startTime};
        simDeltaTime = executionTime.count() * simSecsPerRealSec;

        auto ups = onesecond / executionTime;
        std::cout << std::format("UPS = {0:10.5}", ups) << std::endl;
    }

    freeGravityMemory(spaceobjectsCuda, spaceobjectsStaticCopy, objectcollisions);
}
#endif

void workerThreadFunc(unsigned int width, unsigned int height, unsigned int maxElevation, int continentCount) {
	auto startTime = std::chrono::high_resolution_clock::now();
	
    Gravity(height, width);

	auto endTime = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> executionTime{endTime - startTime};
	std::cout << std::format("GenerateMap took {0:} seconds", executionTime) << std::endl;
}

void displayData(SimpleCanvas* drawer) {
	unsigned int winHeight = drawer->WINDOW_HEIGHT;
	unsigned int winWidth = drawer->WINDOW_WIDTH;

	std::vector<int>& screenPixels = drawer->GetPixelStore();

	std::vector<unsigned int> localScreenData(winHeight * winWidth);

	imageDataLock.acquire();
	std::copy(imageData.begin(), imageData.end(), localScreenData.begin());
	imageDataLock.release();
	for (unsigned int y=0; y<winHeight; y++) {
		for (unsigned int x=0; x<winWidth; x++) {
			unsigned int pxlVal = localScreenData[y*winWidth + x];
			screenPixels[y * winWidth + x] = (pxlVal<<8);
		}
	}
}

int main(int argc, char* argv[]) {

	SimpleCanvas drawer(HEIGHT, WIDTH);
    drawer.RegisterEventHandler(handleInput);

	imageData = std::vector<unsigned int>(drawer.WINDOW_WIDTH * drawer.WINDOW_HEIGHT);
	
	std::thread workerThread(workerThreadFunc, drawer.WINDOW_WIDTH, drawer.WINDOW_HEIGHT, 0x00FFFFFF, 3);

	drawer.Run(displayData);

    running = false;
	workerThread.join();

	return 0;
}