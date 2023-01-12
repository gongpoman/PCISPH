#include<iostream>
#include<vector>

#include "myutil.h"
#include "particle.h"

// TODO 잘 돌아갈 때 변수 쳐내는 식으로 가자
// TODO 일단 조금 비효율적인 부분이 전체 particle을 가지고 논다는점?

PCISPH::PCISPH(glm::ivec3 numOfFluidParticles, glm::vec3 boundarySideLen,bool drawAll) {

    drawWall = drawAll;

	restDensity = 997.0f;
	radius = 0.025f;
	particleMass = restDensity * 4.0f / 3.0f * 3.141592f * radius * radius * radius;
	coreRad = radius * 5.0f;

    numFluidParticlesX = numOfFluidParticles.x;
    numFluidParticlesY = numOfFluidParticles.y;
    numFluidParticlesZ = numOfFluidParticles.z;
    numFluidParticles = numFluidParticlesX * numFluidParticlesY * numFluidParticlesZ;

	boundaryX = boundarySideLen.x;
	boundaryY = boundarySideLen.y;
	boundaryZ = boundarySideLen.z;

	nGridDivX = (unsigned int)(boundaryX / coreRad);
	nGridDivY = (unsigned int)(boundaryY / coreRad);
	nGridDivZ = (unsigned int)(boundaryZ / coreRad);

	upperMaxGridDiv = (nGridDivX > nGridDivY) ?
		((nGridDivX > nGridDivZ) ? nGridDivX : nGridDivZ) :
		((nGridDivY > nGridDivZ) ? nGridDivY : nGridDivZ);

	upperMaxGridDiv = __nthDigit(upperMaxGridDiv);

	eta = 0.01f;
	minIter = 0;
	maxIter = 100;

    particlesSorter = Z_Sort();
    fluidParticlesSorter = Z_Sort();
    wallParticlesSorter = Z_Sort();

	__sceneInitialize();
}
void PCISPH::__sceneInitialize() {
	__brickInit();
	__gridInit();
    __particleInit(1);
    __sorterInit();

}

void PCISPH::__brickInit() {
	const float wallDensity = restDensity; // wall의 전체 density가 물의 density와 같도록 설정.

	Particle wallParticle = Particle();
	wallParticle.isWall = true;
	wallParticle.density = wallDensity;
	wallParticle.pos = glm::vec3(0);

	wallParticle.pressure = 0.0f; // TODO calc pressure 여기 박아넣을 값 찾아야 됨. 만약 안되면 음.. estimation 계속 해줘야 됨.

	int numRow = (int)(coreRad / radius); // 5.0

	for (unsigned int i = 0; i < numRow; i++) {
		for (unsigned int j = 0; j < numRow; j++) {
			for (unsigned int k = 0; k < numRow; k++) {
				wallParticle.pos = glm::vec3(0) + glm::vec3(i * radius, j * radius, k * radius);
				wallParticle.px = wallParticle.pos;
				brick.push_back(wallParticle);
			}
		}
	}
}
void PCISPH::__gridInit() {
//    neighborIdices = new std::vector<glm::ivec3>[nGridDivX * nGridDivY * nGridDivZ]; 

    unsigned int numBoundary = 0;
    std::vector<glm::vec3> boundaryOrigins;

    for (unsigned int i = 0; i < nGridDivX; i++) {
        for (unsigned int j = 0; j < nGridDivY; j++) {
            for (unsigned int k = 0; k < nGridDivZ; k++) {
                //unsigned int cidx = mortonEncode(i, j, k);

                bool isBoundaryCell = false;

                for (int ii = -1; ii < 2; ii++) {
                    for (int jj = -1; jj < 2; jj++) {
                        for (int kk = -1; kk < 2; kk++) {

                            bool x1 = i == 0 && ii == -1;
                            bool x2 = i == nGridDivX - 1 && ii == 1;
                            bool y1 = j == 0 && jj == -1;
                            bool y2 = j == nGridDivY - 1 && jj == 1;
                            bool z1 = k == 0 && kk == -1;
                            bool z2 = k == nGridDivZ - 1 && kk == 1;
                            /*
                            if (!(x1 || x2 || y1 || y2 || z1 || z2)) {
                                glm::ivec3 coord = glm::ivec3(i + ii, j + jj, k + kk);
                                unsigned int neighborCell = __cubeMorton(coord.x, coord.y, coord.z);

                                if (neighborCell != (unsigned int)(-1))    
                                    neighborIdices[cidx].push_back(neighborCell);  
                            }
                            else {
                                isBoundaryCell = true;
                            }
                            */
                            if ((x1 || x2 || y1 || y2 || z1 || z2)) 
                                isBoundaryCell = true;
                        }
                    }
                }

                if (isBoundaryCell) {
                    numBoundary++;
                    boundaryOrigins.push_back(glm::vec3(i * coreRad, j * coreRad, k * coreRad));
                }
            }
        }
    }

    numWallParticles = brick.size() * numBoundary; // determined number of boundary particle 

    numParticles = numWallParticles + numFluidParticles;
    
    numDrawParticles = (drawWall) ? numParticles : numFluidParticles;

    __particlesArrInit(&wallParticles,numWallParticles,WALL);
//  predParticles = new Particle[numParticle];

    for (auto iter = boundaryOrigins.begin(); iter != boundaryOrigins.end(); iter++) {
        __addWallParticles(*iter);
    }

}
void PCISPH::__particleInit(int mode) {

    __particlesArrInit(&particles, numParticles,ALL);
    __particlesArrInit(&fluidParticles, numFluidParticles,FLUID);

    float totMass = numFluidParticles* particleMass / restDensity;
    float spacing = std::pow(totMass/ (float)((numFluidParticlesX - 1)* (numFluidParticlesY - 1)* (numFluidParticlesZ - 1)),1.0f/3.0f);

    sideLenX = (float)numFluidParticlesX * spacing;
    sideLenY = (float)numFluidParticlesY * spacing;
    sideLenZ = (float)numFluidParticlesZ * spacing;

    if (mode == 0) {
        for (unsigned int i = 0; i < numFluidParticlesX; i++) {
            for (unsigned int j = 0; j < numFluidParticlesY; j++) {
                for (unsigned int k = 0; k < numFluidParticlesZ; k++) {

                    Particle tempParticle = Particle();

                    float x = boundaryX / 2.0f + sideLenX * (((float)i / (float)numFluidParticlesX) - 0.5f);
                    float y = sideLenY/(float)numFluidParticlesY * j;
                    float z = boundaryZ / 2.0f + sideLenZ * (((float)k / (float)numFluidParticlesZ) - 0.5f);

                    tempParticle.pos = glm::vec3(x,y,z);
                    glm::ivec3 cellCoordXYZ = __getCellCoord(tempParticle.pos);
                    tempParticle.cellIdx = mortonEncode(cellCoordXYZ.x, cellCoordXYZ.y, cellCoordXYZ.z );

                    __appendParticle(&fluidParticles, tempParticle,FLUID);

                }
            }
        }
    }

    if (mode == 1) {

        for (unsigned int i = 0; i < numFluidParticlesX; i++) {
            for (unsigned int j = 0; j < numFluidParticlesY; j++) {
                for (unsigned int k = 0; k < numFluidParticlesZ; k++) {

                    Particle tempParticle = Particle();


                    float x = sideLenX * ((float)i / (float)numFluidParticlesX);
                    float y = sideLenY / (float)numFluidParticlesY * j;
                    float z = sideLenZ * ((float)k / (float)numFluidParticlesZ);

                    tempParticle.pos = glm::vec3(x, y, z);
                    glm::ivec3 cellCoordXYZ = __getCellCoord(tempParticle.pos);
                    tempParticle.cellIdx = mortonEncode(cellCoordXYZ.x, cellCoordXYZ.y, cellCoordXYZ.z);

                    __appendParticle(&fluidParticles, tempParticle,FLUID);

                }
            }
        }

    }
}


void PCISPH::__addWallParticles(glm::vec3 cellOrigin) {

    std::vector<Particle> tempBrick = brick;

    for (unsigned int i = 0; i < brick.size(); i++) {
        tempBrick[i].pos += cellOrigin;
        tempBrick[i].px = tempBrick[i].pos;

        glm::ivec3 cellIdxCoord = __getCellCoord(tempBrick[i].pos);
        tempBrick[i].cellIdx = mortonEncode(cellIdxCoord.x, cellIdxCoord.y, cellIdxCoord.z);
        tempBrick[i].cellIdxPred = tempBrick[i].cellIdx;

        __appendParticle(&wallParticles,tempBrick[i],WALL);
    }

}
glm::ivec3 PCISPH::__getCellCoord(glm::vec3 position) {
    float x = position.x / boundaryX * (float) nGridDivX;
    float y = position.y / boundaryY * (float)nGridDivY;
    float z = position.z / boundaryZ * (float)nGridDivZ;

    if (x <= -FLT_EPSILON || y <= -FLT_EPSILON || z <= -FLT_EPSILON)
        return glm::ivec3(-1, -1, -1);

    return glm::ivec3((int)x, (int)y, (int)z);
}
glm::vec3 PCISPH::__gridLocalPos(glm::vec3 pos) {
    float x = (pos.x / coreRad < FLT_EPSILON) ? ((int)(pos.x / coreRad) - 1) * coreRad : (int)(pos.x / coreRad) * coreRad;
    float y = (pos.y / coreRad < FLT_EPSILON) ? ((int)(pos.y / coreRad) - 1) * coreRad : (int)(pos.y / coreRad) * coreRad;
    float z = (pos.z / coreRad < FLT_EPSILON) ? ((int)(pos.z / coreRad) - 1) * coreRad : (int)(pos.z / coreRad) * coreRad;

    return pos - glm::vec3(x, y, z);
}

void PCISPH::__particlesArrInit(ParticlesArray* particleArr,unsigned int numbers,ARRTYPE arrtype) {

    particleArr->count = 0;
    particleArr->numMaxParticles = numbers;

    particleArr->isWall = new bool[numbers];

    particleArr->cellIdx = new int[numbers];
    particleArr->cellIdxPred = new int[numbers];

    particleArr->pressure = new float[numbers];
    particleArr->density = new float[numbers];
    particleArr->densityVar = new float[numbers];

    particleArr->pos= new glm::vec3[numbers];
    particleArr->px= new glm::vec3[numbers];

    particleArr->force_g = glm::vec3(0, -9.8f, 0) * particleMass;

    if (arrtype != WALL) {
        particleArr->vel = new glm::vec3[numbers];
        particleArr->pv = new glm::vec3[numbers];

        particleArr->force_p = new glm::vec3[numbers];
        particleArr->force_ext = new glm::vec3[numbers];
        particleArr->force_vis = new glm::vec3[numbers];
    }
}

void PCISPH::__appendParticle(ParticlesArray* particleArr, Particle particle, ARRTYPE arrtype) {

    if (particleArr->count >= particleArr->numMaxParticles) {
        std::cout << "ERROR :: PARTICLE ARRAY IS FULL" << std::endl;
        system("PAUSE");
        return;
    }

    unsigned int idx = particleArr->count;

    particleArr->isWall[idx] = particle.isWall;

    particleArr->cellIdx[idx] = particle.cellIdx;
    particleArr->cellIdxPred[idx] = particle.cellIdxPred;

    particleArr->pressure[idx] = particle.pressure;
    particleArr->density[idx] = particle.density;
    particleArr->densityVar[idx] = particle.densityVar;

    particleArr->pos[idx] = particle.pos;
    particleArr->px[idx] = particle.px;

    if (arrtype != WALL) {
        particleArr->vel[idx] = particle.vel;
        particleArr->pv[idx] = particle.pv;
        particleArr->force_p[idx] = particle.force_p;
        particleArr->force_ext[idx] = particle.force_ext;
        particleArr->force_vis[idx] = particle.force_vis;
    }


    particleArr->count++;

}
void PCISPH::__deleteParticleArr(ParticlesArray* partArr, ARRTYPE arrtype) {
    delete[] partArr->isWall;

    delete[] partArr->cellIdx;
    delete[] partArr->cellIdxPred;

    delete[] partArr->pressure;
    delete[] partArr->density;
    delete[] partArr->densityVar;

    delete[] partArr->pos;
    delete[] partArr->px;

    if (arrtype != WALL) {
        delete[] partArr->vel;
        delete[] partArr->pv;
        delete[] partArr->force_p;
        delete[] partArr->force_ext;
        delete[] partArr->force_vis;
    }
}

void PCISPH::__sorterInit() {
    unsigned int maxCIdx = mortonEncode(nGridDivX - 1, nGridDivY - 1, nGridDivZ - 1);

    particlesSorter = Z_Sort(numParticles, &particles, ALL, maxCIdx);
    fluidParticlesSorter = Z_Sort(numFluidParticles, &fluidParticles, FLUID, maxCIdx);
    wallParticlesSorter = Z_Sort(numWallParticles, &wallParticles, WALL, maxCIdx);


    wallParticlesSorter.sortby(ACTUAL_POS);             //wall particle 정렬

    fluidParticlesSorter.sortby(ACTUAL_POS);            //fluid particle 정렬

    /*
    //TEST
    std::cout << "파티클 갯수 : " << numFluidParticles << std::endl;
    for (int i = 0; i < numFluidParticles; i++) {
        std::cout << i << " 번쨰 particle : " << fluidParticles.cellIdx[i] << " : " << mortonDecode(fluidParticles.cellIdx[i]).x << " : " << mortonDecode(fluidParticles.cellIdx[i]).y << " : " << mortonDecode(fluidParticles.cellIdx[i]).z << " : " << fluidParticles.pos[i].x << " : " << fluidParticles.pos[i].y << " : " << fluidParticles.pos[i].z << std::endl;
    }
    system("PAUSE");

    std::cout << "이건 WALL 몇개만 출력해본거" << std::endl;

    for (int i = 0; i < 1001; i++) {
        std::cout << i << " 번쨰 particle : " << wallParticles.cellIdx[i] << " : " << mortonDecode(wallParticles.cellIdx[i]).x << " : " << mortonDecode(wallParticles.cellIdx[i]).y << " : " << mortonDecode(wallParticles.cellIdx[i]).z << " : " << wallParticles.pos[i].x << " : " << wallParticles.pos[i].y << " : " << wallParticles.pos[i].z << std::endl;
    }
    system("PAUSE");

    */

    // 정렬된 두개로 전체에 대한 Arr을 만든다.
    particlesSorter.integrateSortedArr(wallParticlesSorter, fluidParticlesSorter, ACTUAL_POS);

}

PCISPH::~PCISPH() {
    __deleteParticleArr(&particles,ALL);
    __deleteParticleArr(&fluidParticles,FLUID);
    __deleteParticleArr(&wallParticles,WALL);
}





void PCISPH::update() {

    ////////////////////////////////////////////////////////////
    /// initial state calculation
    ////////////////////////////////////////////////////////////       
    for (int i = 0; i < numFluidParticles; i++) {
        fluidParticles.pressure[i] = 0.0f;
        fluidParticles.force_p[i] = glm::vec3(0);
        fluidParticles.force_ext[i] = glm::vec3(0);
        fluidParticles.force_vis[i] = glm::vec3(0);
    }
    //////////////////////////////////////////////////////////// 
    fluidParticlesSorter.sortby(ACTUAL_POS);
    particlesSorter.integrateSortedArr(wallParticlesSorter, fluidParticlesSorter, ACTUAL_POS);

    // accumulation counter
    // neighbor search

}

//==================================================================================
PCISPH::Z_Sort::Z_Sort() {
    count = 0;
    particleArr = NULL;
}
PCISPH::Z_Sort::Z_Sort(int count, ParticlesArray* pParticleArr,ARRTYPE type,unsigned int maxCellIndex) {
    arrT = type;
    this->count = count;
    this->particleArr = pParticleArr;
    maxCellIdx = maxCellIndex;
}

void PCISPH::Z_Sort::sortby(SORTMODE mode) {
    
    currentMode = mode;
    __buildCounter(); // offset을 만들어 놨음. cellIdx에 있는 particle이 몇개있는지 앎.

    unsigned int* sortedIdx = __buildSortedIdx();

    /*  TO TEST
    if(arrT == FLUID)
    for(int i=0; i<count; i++)
        std::cout << sortedIdx[i] << " : " << mortonDecode(particleArr->cellIdx[sortedIdx[i]]).x * 0.125f << " : " << mortonDecode(particleArr->cellIdx[sortedIdx[i]]).y * 0.125f << " : " << mortonDecode(particleArr->cellIdx[sortedIdx[i]]).z * 0.125f << std::endl;
    system("PAUSE ");
    */

    //sortedIdx를 이용해서 정렬 해야 됨.

    ParticlesArray tempArr;

    __particlesArrInit(&tempArr, count, arrT);
    __swapAndSort(&tempArr,sortedIdx);

    delete[] sortedIdx;
}

// fluidSorter를 이용해서 정렬을 한 다음 전체 particle을 만들기 위해 이 함수를 호출함.
void PCISPH::Z_Sort::integrateSortedArr(Z_Sort wallSorter, Z_Sort fluidSorter,SORTMODE mode) {

    currentMode = mode;

    __integrateCounter(wallSorter, fluidSorter);

    int* sortedIdx = __integrateSortedIdx(wallSorter,fluidSorter);

    __integrate(wallSorter, fluidSorter, sortedIdx);


    delete[] sortedIdx;
}
void PCISPH::Z_Sort::__integrateCounter(Z_Sort wallSorter, Z_Sort fluidSorter) {
    //wallParticle과 fluidParticle의 hashMap합친다
    particleCounter = wallSorter.particleCounter;

    // mode에 따라서 cellIdx또는 cellIdxPred에 있는 particle의 갯수를 샘.
    for (auto ele : fluidSorter.particleCounter) {
        if (particleCounter.count(ele.first)) {
            particleCounter[ele.first] += ele.second;
        }
        else
            particleCounter[ele.first] = ele.second;
    }

    //occupiedCellIdx 만든다.
    occupiedCellIdx.clear();
    occupiedCellIdxInv.clear();

    unsigned int occIdx = 0;
    for (unsigned int i = 0; i <= maxCellIdx; i++) {
        if (particleCounter.count(i)) {
            occupiedCellIdx.push_back(i);
            occupiedCellIdxInv[i] = occIdx++;
        }
    }

}
int* PCISPH::Z_Sort::__integrateSortedIdx(Z_Sort wallSorter, Z_Sort fluidSorter) {

    int* sortIdxArr = new int[count];

    std::vector<std::vector<int>> sortParticleIdx;
       
    for (unsigned int i = 0; i < occupiedCellIdx.size(); i++)
        sortParticleIdx.push_back(std::vector<int>());

    if (currentMode == ACTUAL_POS) {
        for (unsigned int i = 0; i < wallSorter.count; i++) {
            unsigned int idx = occupiedCellIdxInv[wallSorter.particleArr->cellIdx[i]];
            sortParticleIdx[idx].push_back(i);
        }
        for (unsigned int i = 1; i < fluidSorter.count+1; i++) {
            unsigned int idx = occupiedCellIdxInv[fluidSorter.particleArr->cellIdx[i]];
            sortParticleIdx[idx].push_back(-1* i);
        }
    }
    else if (currentMode == PRED_POS) {
        for (unsigned int i = 0; i < wallSorter.count; i++) {
            unsigned int idx = occupiedCellIdxInv[wallSorter.particleArr->cellIdxPred[i]];
            sortParticleIdx[idx].push_back(i);
        }
        for (unsigned int i = 1; i < fluidSorter.count+1; i++) {
            unsigned int idx = occupiedCellIdxInv[fluidSorter.particleArr->cellIdxPred[i]];
            sortParticleIdx[idx].push_back(-1 * i);
        }
    }

    unsigned int iii = 0;
    for (auto iter = sortParticleIdx.begin(); iter != sortParticleIdx.end(); iter++) {
        for (auto iteriter = iter->begin(); iteriter != iter->end(); iteriter++) {
            sortIdxArr[iii] = *iteriter;
            iii++;
        }
    }
    return sortIdxArr;

}
void PCISPH::Z_Sort::__integrate(Z_Sort wallSorter, Z_Sort fluidSorter, int* sortedIdx) {
    
    for (unsigned int i = 0; i < count; i++) {
        int curIdx = sortedIdx[i];
        bool isFluid = (curIdx < 0);
        if (isFluid) {
            curIdx = -curIdx -1;

            particleArr->isWall[i] = false;

            particleArr->cellIdx[i] = fluidSorter.particleArr->cellIdx[curIdx];
            particleArr->cellIdxPred[i] = fluidSorter.particleArr->cellIdxPred[curIdx];

            particleArr->pressure[i] = fluidSorter.particleArr->pressure[curIdx];
            particleArr->density[i] = fluidSorter.particleArr->density[curIdx];
            particleArr->densityVar[i] = fluidSorter.particleArr->densityVar[curIdx];

            particleArr->vel[i] = fluidSorter.particleArr->vel[curIdx];
            particleArr->pos[i] = fluidSorter.particleArr->pos[curIdx];
            particleArr->px[i] = fluidSorter.particleArr->px[curIdx];
            particleArr->pv[i] = fluidSorter.particleArr->pv[curIdx];

            particleArr->force_p[i] = fluidSorter.particleArr->force_p[curIdx];
            particleArr->force_ext[i] = fluidSorter.particleArr->force_ext[curIdx];
            particleArr->force_vis[i] = fluidSorter.particleArr->force_vis[curIdx];
        }
    }
    for (unsigned int i = 0; i < count; i++) {
        int curIdx = sortedIdx[i];

        bool isFluid = (curIdx < 0);
        if (!isFluid) {
            particleArr->isWall[i] = true;

            particleArr->cellIdx[i] = wallSorter.particleArr->cellIdx[curIdx];
            particleArr->cellIdxPred[i] = wallSorter.particleArr->cellIdxPred[curIdx];

            particleArr->pressure[i] = wallSorter.particleArr->pressure[curIdx];
            particleArr->density[i] = wallSorter.particleArr->density[curIdx];
            particleArr->densityVar[i] = wallSorter.particleArr->densityVar[curIdx];

            particleArr->pos[i] = wallSorter.particleArr->pos[curIdx];
            particleArr->px[i] = wallSorter.particleArr->px[curIdx];

        }
    }
}



void PCISPH::Z_Sort::__buildCounter() {

    particleCounter.clear();
    
    if (currentMode == ACTUAL_POS) {

        for (unsigned int i = 0; i < particleArr->count; i++) {

            if (particleCounter.count(particleArr->cellIdx[i]))
                particleCounter[particleArr->cellIdx[i]]++;
            else
                particleCounter[particleArr->cellIdx[i]] = 1;
        }
    }
    else if (currentMode == PRED_POS) {

        for (unsigned int i = 0; i < particleArr->count; i++) {

            if (particleCounter.count(particleArr->cellIdxPred[i]))
                particleCounter[particleArr->cellIdxPred[i]]++;
            else
                particleCounter[particleArr->cellIdxPred[i]] = 1;
        }
    }



    occupiedCellIdx.clear();
    occupiedCellIdxInv.clear();

    unsigned int occIdx = 0;
    for (unsigned int i = 0; i <= maxCellIdx; i++) {
        if (particleCounter.count(i)) {
            occupiedCellIdx.push_back(i);
            occupiedCellIdxInv[i] = occIdx++;
        }
    }


    /* TODO 이거 나중에 neighbor search 할 때 필요함. 두개 합치기 전에는 이거를 하면 안됨.
    for (auto iter = occupiedCellIdx.begin() + 1; iter != occupiedCellIdx.end(); iter++) {
        particleCounter[*iter] += particleCounter[(*iter) - 1];
    }
    */




    //test
        /*
    if (arrT == WALL) {
        std::cout << "WALL" << std::endl;
    

        std::cout << occupiedCellIdx.size() << std::endl;
        for (auto ele : occupiedCellIdx) {
            std::cout << ele << std::endl;
        }
        std::cout << maxCellIdx << std::endl;

        int sumsum = 0;
        for (auto element : particleCounter) {
            sumsum++;
            std::cout << "CELLIDX : " << element.first << " COUNT :  " << element.second << std::endl;
            if (element.second  != 125) {
                std::cout <<"왜 cell에 125개가 아닌거지."<< element.first<< " : " << element.second << std::endl;
                std::cout << mortonDecode(element.first).x * 0.125f<<" : " << mortonDecode(element.first).y * 0.125f << " : " << mortonDecode(element.first).z * 0.125f << std::endl;
                // TODOTEST demorton 해서 origin 어딘지. 거기 속한 particle들 좌표 뭔지. 몇번째 particle들인지.
                system("PAUSE");
            }

        }
        system("PAUSE");
    }

        */
}
unsigned int* PCISPH::Z_Sort::__buildSortedIdx() {

    unsigned int *sortIdxArr = new unsigned int[count];

    std::vector<std::vector<unsigned int>> sortParticleIdx;

    for (unsigned int i = 0; i < occupiedCellIdx.size(); i++)
        sortParticleIdx.push_back(std::vector<unsigned int>());


    if (currentMode == ACTUAL_POS) {
        for (unsigned int i = 0; i < count; i++) {
            unsigned int idx = occupiedCellIdxInv[particleArr->cellIdx[i]];
            sortParticleIdx[idx].push_back(i);
        }
    }
    else if (currentMode == PRED_POS) {
        for (unsigned int i = 0; i < count; i++) {
            unsigned int idx = occupiedCellIdxInv[particleArr->cellIdxPred[i]];
            sortParticleIdx[idx].push_back(i);
        }
    }

    unsigned int iii = 0;
    for (auto iter = sortParticleIdx.begin(); iter != sortParticleIdx.end(); iter++) {
        for (auto iteriter = iter->begin(); iteriter != iter->end(); iteriter++) {
            sortIdxArr[iii] = *iteriter;
            iii++;
        }
    }
    return sortIdxArr;
}


void PCISPH::Z_Sort::__particlesArrInit(ParticlesArray* particleArr,unsigned int numbers,ARRTYPE arrtype) {

    particleArr->isWall = new bool[numbers];

    particleArr->cellIdx = new int[numbers];
    particleArr->cellIdxPred = new int[numbers];

    particleArr->pressure = new float[numbers];
    particleArr->density = new float[numbers];
    particleArr->densityVar = new float[numbers];

    particleArr->pos= new glm::vec3[numbers];
    particleArr->px= new glm::vec3[numbers];

    if (arrtype != WALL) {
        particleArr->pv = new glm::vec3[numbers];
        particleArr->vel = new glm::vec3[numbers];
        particleArr->force_p = new glm::vec3[numbers];
        particleArr->force_ext = new glm::vec3[numbers];
        particleArr->force_vis = new glm::vec3[numbers];
    }
}
void PCISPH::Z_Sort::__swapAndSort(ParticlesArray* temp, unsigned int* sortedIdx) {
    // temp에 정렬해서 만들기.

    if (arrT == WALL) {
        for (unsigned int i = 0; i < count; i++) {

            unsigned int curIdx = sortedIdx[i];

            temp->isWall[i] = particleArr->isWall[curIdx];

            temp->cellIdx[i] = particleArr->cellIdx[curIdx];
            temp->cellIdxPred[i] = particleArr->cellIdxPred[curIdx];

            temp->pressure[i] = particleArr->pressure[curIdx];
            temp->density[i] = particleArr->density[curIdx];
            temp->densityVar[i] = particleArr->densityVar[curIdx];

            temp->pos[i] = particleArr->pos[curIdx];
            temp->px[i] = particleArr->px[curIdx];
        }


        delete[] particleArr->pos;
        delete[] particleArr->px;
        
        delete[] particleArr->pressure;
        delete[] particleArr->density;
        delete[] particleArr->densityVar;
         
        delete[] particleArr->isWall ;
         
        delete[] particleArr->cellIdx;
        delete[] particleArr->cellIdxPred;


        particleArr->pos = temp->pos;
        particleArr->px = temp->px;

        particleArr->pressure = temp->pressure;
        particleArr->density = temp->density;
        particleArr->densityVar = temp->densityVar;

        particleArr->isWall = temp->isWall;

        particleArr->cellIdx = temp->cellIdx;
        particleArr->cellIdxPred = temp->cellIdxPred;

    }
    else {
        for (unsigned int i = 0; i < count; i++) {

            unsigned int curIdx = sortedIdx[i];

            temp->isWall[i] = particleArr->isWall[curIdx];
                                   
            temp->cellIdx[i] = particleArr->cellIdx[curIdx];
            temp->cellIdxPred[i] = particleArr->cellIdxPred[curIdx];
                                   
            temp->pressure[i] = particleArr->pressure[curIdx];
            temp->density[i] = particleArr->density[curIdx];
            temp->densityVar[i] = particleArr->densityVar[curIdx];
                                   
            temp->vel[i] = particleArr->vel[curIdx];
            temp->pos[i] = particleArr->pos[curIdx];
            temp->pv[i] = particleArr->pv[curIdx];
            temp->px[i] = particleArr->px[curIdx];
                                   
            temp->force_p[i] = particleArr->force_p[curIdx];
            temp->force_ext[i] = particleArr->force_ext[curIdx];
            temp->force_vis[i] = particleArr->force_vis[curIdx];
        }

        delete[] particleArr->force_p;
        delete[] particleArr->force_vis;
        delete[] particleArr->force_ext;

        delete[] particleArr->pos;
        delete[] particleArr->px ;
        delete[] particleArr->vel;
        delete[] particleArr->pv ;

        delete[] particleArr->pressure;
        delete[] particleArr->density;
        delete[] particleArr->densityVar;

        delete[] particleArr->isWall;
        
        delete[] particleArr->cellIdx;
        delete[] particleArr->cellIdxPred;


        particleArr->force_p = temp->force_p;
        particleArr->force_vis = temp->force_vis;
        particleArr->force_ext = temp->force_ext;

        particleArr->pos = temp->pos;
        particleArr->px = temp->px;
        particleArr->vel = temp->vel;
        particleArr->pv = temp->pv;

        particleArr->pressure = temp->pressure;
        particleArr->density = temp->density;
        particleArr->densityVar = temp->densityVar;

        particleArr->isWall = temp->isWall;

        particleArr->cellIdx = temp->cellIdx;
        particleArr->cellIdxPred = temp->cellIdxPred;
    }


    // 기존에 있었던 것 삭제하기.
}