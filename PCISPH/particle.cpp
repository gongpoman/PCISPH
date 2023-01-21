#include<iostream>
#include<vector>
#include<unordered_map>

#include "myutil.h"
#include "particle.h"


void printVector(glm::vec3 x) {
    std::cout << x.x << " : " << x.y << " : " << x.z << std::endl;
}


// TODO �� ���ư� �� ���� �ĳ��� ������ ����
// TODO �ϴ� ���� ��ȿ������ �κ��� ��ü particle�� ������ ��ٴ���?

PCISPH::PCISPH(glm::ivec3 numOfFluidParticles, glm::vec3 boundarySideLen,bool drawAll) {

    drawWall = drawAll;

	restDensity = 997.0f;
	radius = 0.025f;
	particleMass = restDensity * 4.0f * 3.141592f * radius * radius * radius / 3.0f ;
	coreRad = radius * 3.0f;
    wallMassRatio = 0.0f;
    numBrickRow = 2;

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
	minIter = 3;
	maxIter = 30;


//    particlesSorter = Z_Sort();
    fluidParticlesSorter = Z_Sort();
    wallParticlesSorter = Z_Sort();


	__sceneInitialize();
}
void PCISPH::__sceneInitialize() {
	__brickInit();
	__gridInit();
    __particleInit(0);
    __sorterInit();
    __calcDelta();
}
void PCISPH::__brickInit() {
	const float wallDensity = restDensity;//4.0f/3.0f*3.141592f*restDensity; // wall�� ��ü density�� ���� density�� ������ ����.

	Particle wallParticle = Particle();
	wallParticle.isWall = true;
	wallParticle.density = wallDensity;
	wallParticle.pos = glm::vec3(0);

	wallParticle.pressure = 0.0f; // TODO calc pressure ���� �ھƳ��� �� ã�ƾ� ��. ���� �ȵǸ� ��.. estimation ��� ����� ��.

	int numBrickRow = 2;
    float spacingRatio = coreRad/ radius /((float)numBrickRow);

    wallMassRatio = spacingRatio * spacingRatio * spacingRatio * 3.0f / 4.0f / 3.141592f ;
    
	for (unsigned int i = 0; i < numBrickRow; i++) {
		for (unsigned int j = 0; j < numBrickRow; j++) {
			for (unsigned int k = 0; k < numBrickRow; k++) {
				wallParticle.pos = glm::vec3(0) + spacingRatio*radius*glm::vec3(i , j , k );
				wallParticle.px = wallParticle.pos;
				brick.push_back(wallParticle);
			}
		}
	}
}
void PCISPH::__gridInit()  {

    unsigned int numBoundaryFull = 0;
    unsigned int numBoundaryHH = 0;
    unsigned int numBoundaryHHH = 0;

    std::vector<glm::vec3> boundaryOrigins;

    for (unsigned int i = 0; i < nGridDivX; i++) {
        for (unsigned int j = 0; j < nGridDivY; j++) {
            for (unsigned int k = 0; k < nGridDivZ; k++) {

                int isBoundaryCell = 0;

                for (int ii = -1; ii < 2; ii++) {
                    for (int jj = -1; jj < 2; jj++) {
                        for (int kk = -1; kk < 2; kk++) {

                            bool x1 = i == 0 && ii == -1;
                            bool x2 = i == nGridDivX - 1 && ii == 1;
                            bool y1 = j == 0 && jj == -1;
                            bool y2 = j == nGridDivY - 1 && jj == 1;
                            bool z1 = k == 0 && kk == -1;
                            bool z2 = k == nGridDivZ - 1 && kk == 1;

                            if (x1 || x2 || y1 || y2 || z1 || z2) {
                                    isBoundaryCell++;
                            }
                        }
                    }
                }
                if (isBoundaryCell > 0) {

                    if (isBoundaryCell == 9) {

                        boundaryOrigins.push_back(glm::vec3(i * coreRad, j * coreRad, k * coreRad));
                        numBoundaryFull++;
                    }
                    else if (isBoundaryCell == 15) {
                        boundaryOrigins.push_back(glm::vec3(i * coreRad, j * coreRad, k * coreRad));
                        numBoundaryHH++;
                    }
                    else if (isBoundaryCell == 19){
                        boundaryOrigins.push_back(glm::vec3(i * coreRad, j * coreRad, k * coreRad));
                        numBoundaryHHH++;
                        }
                     
                }
            }
        }
    }

    numWallParticles = brick.size() * numBoundaryFull + numBoundaryHH * brick.size() / (numBrickRow * numBrickRow) + numBoundaryHHH * brick.size() / (numBrickRow * numBrickRow * numBrickRow); // determined number of boundary particle 
    numWallParticles = brick.size() * (numBoundaryFull + numBoundaryHH  + numBoundaryHHH); // determined number of boundary particle 

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

    float spacing;
    if (numFluidParticlesX == 1 && numFluidParticlesY == 1 && numFluidParticlesZ == 1)
        spacing = 0;
    else
    {
        spacing = std::pow(particleMass / restDensity, 1.0f / 3.0f);
    }

    sideLenX = (float)(numFluidParticlesX+1) * spacing;
    sideLenY = (float)(numFluidParticlesY+1) * spacing;
    sideLenZ = (float)(numFluidParticlesZ+1) * spacing;


    unsigned int pID = 0;

    if (mode == 0) {
        for (unsigned int i = 0; i < numFluidParticlesX; i++) {
            for (unsigned int j = 0; j < numFluidParticlesY; j++) {
                for (unsigned int k = 0; k < numFluidParticlesZ; k++) {

                    Particle tempParticle = Particle();

                    float x = boundaryX / 2.0f + sideLenX * (((float)i / (float)numFluidParticlesX) - 0.5f);
                    float y = sideLenY/(float)numFluidParticlesY * j  + coreRad;
                    float z = boundaryZ / 2.0f + sideLenZ * (((float)k / (float)numFluidParticlesZ) - 0.5f);
                    /*
                    x = 0.0f;
                    z = 0.0f;
                    y = coreRad * 5;
                    */
                    tempParticle.pos = glm::vec3(x,y,z);
                    glm::ivec3 cellCoordXYZ = __getCellCoord(tempParticle.pos);
                    tempParticle.cellIdx = mortonEncode(cellCoordXYZ.x, cellCoordXYZ.y, cellCoordXYZ.z );
                    tempParticle.pID = pID++;
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


                    float x = sideLenX * ((float)i / (float)numFluidParticlesX) + coreRad;
                    float y = sideLenY / (float)numFluidParticlesY * j +coreRad;
                    float z = sideLenZ * ((float)k / (float)numFluidParticlesZ) + coreRad;

                    tempParticle.pos = glm::vec3(x, y, z);
                    glm::ivec3 cellCoordXYZ = __getCellCoord(tempParticle.pos);
                    tempParticle.cellIdx = mortonEncode(cellCoordXYZ.x, cellCoordXYZ.y, cellCoordXYZ.z);
                    tempParticle.pID = pID++;
                    __appendParticle(&fluidParticles, tempParticle,FLUID);

                }
            }
        }
    }


    for (unsigned int i = 0; i < numFluidParticles; i++) {
        glm::ivec3 cellCoord = __getCellCoord(fluidParticles.pos[i]);
        fluidParticles.cellIdx[i] = mortonEncode(cellCoord.x, cellCoord.y, cellCoord.z);
    }

    for (unsigned int i = 0; i < numFluidParticles; i++) {
        fluidParticles.px[i] = fluidParticles.pos[i];
        glm::ivec3 cellCoord = __getCellCoord(fluidParticles.px[i]);
        fluidParticles.cellIdxPred[i] = mortonEncode(cellCoord.x, cellCoord.y, cellCoord.z);
    }
}

void PCISPH::__addWallParticles(glm::vec3 cellOrigin) {

    std::vector<Particle> tempBrick = brick;

    const glm::ivec3 cellCoord = __getCellCoord(cellOrigin);


    int BXFLAG = (cellCoord.x == 0)?-1 : ((cellCoord.x == nGridDivX - 1)? 1 : 0);
    int BYFLAG = (cellCoord.y == 0)?-1 :((cellCoord.y == nGridDivY - 1)? 1 : 0);
    int BZFLAG = (cellCoord.z == 0)?-1 : ((cellCoord.z == nGridDivZ - 1)? 1 : 0);


    if (std::abs(BXFLAG) + std::abs(BYFLAG) + std::abs(BZFLAG) == 1) {

        for (unsigned int i = 0; i < brick.size(); i++) {

            tempBrick[i].pos += cellOrigin;
            tempBrick[i].px = tempBrick[i].pos;

            glm::ivec3 cellIdxCoord = __getCellCoord(tempBrick[i].pos);
            tempBrick[i].cellIdx = mortonEncode(cellIdxCoord.x, cellIdxCoord.y, cellIdxCoord.z);
            tempBrick[i].cellIdxPred = tempBrick[i].cellIdx;


            //////////////////////////////////////////

            tempBrick[i].densityVar = wallMassRatio * 0.95f ;
            
            //////////////////////////////////////////

            __appendParticle(&wallParticles, tempBrick[i], WALL);
        }
    }
    else if ((std::abs(BXFLAG) + std::abs(BYFLAG) + std::abs(BZFLAG)) == 2) {

        for (unsigned int i = 0; i < brick.size(); i++) {

            /*
            if (BXFLAG == 1 && tempBrick[i].pos.x > 10e-5)
                continue;
            if (BXFLAG == -1 && tempBrick[i].pos.x < coreRad * ((float)(numBrickRow - 1)) / ((float)numBrickRow) - 10e-5)
                continue;

            if (BYFLAG == 1 && tempBrick[i].pos.y > 10e-5)
                continue;
            if (BYFLAG == -1 && tempBrick[i].pos.y < coreRad * ((float)(numBrickRow - 1)) / ((float)numBrickRow) - 10e-5)
                continue;

            if (BZFLAG == 1 && tempBrick[i].pos.z > 10e-5)
                continue;
            if (BZFLAG == -1 && tempBrick[i].pos.z < coreRad * ((float)(numBrickRow - 1)) / ((float)numBrickRow) - 10e-5)
                continue;
            */

            tempBrick[i].pos += cellOrigin;
            tempBrick[i].px = tempBrick[i].pos;

            glm::ivec3 cellIdxCoord = __getCellCoord(tempBrick[i].pos);
            tempBrick[i].cellIdx = mortonEncode(cellIdxCoord.x, cellIdxCoord.y, cellIdxCoord.z);
            tempBrick[i].cellIdxPred = tempBrick[i].cellIdx;



            //////////////////////////////
            tempBrick[i].densityVar = wallMassRatio * 2.0f - 1.0f;
            
            // TODO �̷� �������� ���Ǻη� density ���� �ټ��� �ְ�. �̺κ��� �߿���. �̰��� �ֺ�������. �е��� ���������� ���� ���� �ʿ䰡 �־��.
            if (tempBrick[i].pos.y < coreRad * 2.1f)
                tempBrick[i].densityVar = 0.6f;
            else
                tempBrick[i].densityVar = 0.5;
            //////////////////////////////



            __appendParticle(&wallParticles, tempBrick[i], WALL);
        }


    }
    else if ((std::abs(BXFLAG) + std::abs(BYFLAG) + std::abs(BZFLAG)) == 3) {

        for (unsigned int i = 0; i < brick.size(); i++) {
            /*
            if (BXFLAG == 1 && tempBrick[i].pos.x > 10e-5)
                continue;
            if (BXFLAG == -1 && tempBrick[i].pos.x < coreRad * ((float)(numBrickRow - 1)) / ((float)numBrickRow) - 10e-5)
                continue;

            if (BYFLAG == 1 && tempBrick[i].pos.y > 10e-5)
                continue;
            if (BYFLAG == -1 && tempBrick[i].pos.y < coreRad * ((float)(numBrickRow - 1)) / ((float)numBrickRow) - 10e-5)
                continue;

            if (BZFLAG == 1 && tempBrick[i].pos.z > 10e-5)
                continue;
            if (BZFLAG == -1 && tempBrick[i].pos.z < coreRad * ((float)(numBrickRow - 1)) / ((float)numBrickRow) - 10e-5)
                continue;
            */    

            tempBrick[i].pos += cellOrigin;
            tempBrick[i].px = tempBrick[i].pos;

            glm::ivec3 cellIdxCoord = __getCellCoord(tempBrick[i].pos);
            tempBrick[i].cellIdx = mortonEncode(cellIdxCoord.x, cellIdxCoord.y, cellIdxCoord.z);
            tempBrick[i].cellIdxPred = tempBrick[i].cellIdx;


            //////////////////////////////////////////s
            tempBrick[i].densityVar = 2.0f - wallMassRatio;
            tempBrick[i].densityVar = 0.7f; // ������ �������̴� Ű��ϱ� �� �������µ�? ( �ƴ� )
                                            //������ ū �ǹ� ���µ�... �ִ��� ���������. ���� 2root2 r �̶�.

            //////////////////////////////////////////



            __appendParticle(&wallParticles, tempBrick[i], WALL);
        }
    }

}
glm::ivec3 PCISPH::__getCellCoord(glm::vec3 position) {

    float xc = position.x * (float)nGridDivX / boundaryX + 1e-7;
    float yc = position.y * (float)nGridDivY / boundaryY + 1e-7;
    float zc = position.z * (float)nGridDivZ / boundaryZ + 1e-7;

    int xx = (xc <= -FLT_EPSILON) ? -1 : xc;
    int yy = (yc <= -FLT_EPSILON) ? -1 : yc;
    int zz = (zc <= -FLT_EPSILON) ? -1 : zc;

    // �̰� xc �׳� ����ϸ� 13 �����µ� static cast�� �ϰ� �Ǹ� �� 12������. ���ؾ��ϳ�.

    // ����ȯ �� �� ����� �׷� �����ΰŰ����� ��� ����?

    return glm::ivec3(xx,yy,zz);
//    return glm::ivec3((int)xc, (int)yc, (int)zc);
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

        particleArr->particleID = new unsigned int[numbers];
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

        particleArr->particleID[idx] = particle.pID;
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
        delete[] partArr->particleID;
    }
}

void PCISPH::__sorterInit() {
    unsigned int maxCIdx = mortonEncode(nGridDivX - 1, nGridDivY - 1, nGridDivZ - 1);

//    particlesSorter = Z_Sort(numParticles, &particles, ALL, maxCIdx);
    fluidParticlesSorter = Z_Sort(numFluidParticles, &fluidParticles, FLUID, maxCIdx);
    wallParticlesSorter = Z_Sort(numWallParticles, &wallParticles, WALL, maxCIdx);

    wallParticlesSorter.sortby(ACTUAL_POS);             //wall particle ����


    fluidParticlesSorter.sortby(ACTUAL_POS);            //fluid particle ����

    /*
    //TEST
    std::cout << "��ƼŬ ���� : " << numFluidParticles << std::endl;
    for (int i = 0; i < numFluidParticles; i++) {
        std::cout << i << " ���� particle : " << fluidParticles.cellIdx[i] << " : " << mortonDecode(fluidParticles.cellIdx[i]).x << " : " << mortonDecode(fluidParticles.cellIdx[i]).y << " : " << mortonDecode(fluidParticles.cellIdx[i]).z << " : " << fluidParticles.pos[i].x << " : " << fluidParticles.pos[i].y << " : " << fluidParticles.pos[i].z << std::endl;
    }
    system("PAUSE");

    std::cout << "�̰� WALL ��� ����غ���" << std::endl;

    for (int i = 0; i < 1001; i++) {
        std::cout << i << " ���� particle : " << wallParticles.cellIdx[i] << " : " << mortonDecode(wallParticles.cellIdx[i]).x << " : " << mortonDecode(wallParticles.cellIdx[i]).y << " : " << mortonDecode(wallParticles.cellIdx[i]).z << " : " << wallParticles.pos[i].x << " : " << wallParticles.pos[i].y << " : " << wallParticles.pos[i].z << std::endl;
    }
    system("PAUSE");

    */

    // ���ĵ� �ΰ��� ��ü�� ���� Arr�� �����.

    //particlesSorter.integrateSortedArr(wallParticlesSorter, fluidParticlesSorter, ACTUAL_POS);

}

PCISPH::~PCISPH() {
    __deleteParticleArr(&particles,ALL);
    __deleteParticleArr(&fluidParticles,FLUID);
    __deleteParticleArr(&wallParticles,WALL);
}
//===================================================================================

void PCISPH::update() {
    /// initial state calculation
    __sortFluidParticles(ACTUAL_POS);
    __makeNCellMap();


    __initializeFrameStates();

    /* TEST
    for (auto occCellIdxIter = fluidParticlesSorter.occupiedCellIdx.begin();
        occCellIdxIter != fluidParticlesSorter.occupiedCellIdx.end(); occCellIdxIter++) {
        std::cout << *occCellIdxIter << " cellIdx �� neighbor" << std::endl;

        for (auto ele : fluidNCellIdx[*occCellIdxIter]) {
            std::cout << ele << " COORD : " << mortonDecode(ele).x << " : " << mortonDecode(ele).y << " : " << mortonDecode(ele).z << std::endl;
        }
        for (auto ele : wallNCellIdx[*occCellIdxIter]) {
            std::cout << ele << " COORD : " << mortonDecode(ele).x << " : " << mortonDecode(ele).y << " : " << mortonDecode(ele).z << std::endl;
        }
        std::cout << "===================================" << std::endl;
    }
    system("PAUSE");
    */

    bool densityOverFlag = true;
    unsigned int iteration = 0;

    while ((iteration < minIter || densityOverFlag) && iteration < maxIter) {

        /////////////////////////////////////////////
        //predictive
        __pos_vel_Predictive();
        /////////////////////////////////////////////


        __sortFluidParticles(PRED_POS);
        __makePredNCellMap();

        /////////////////////////////////////////////
        // pred Density calc 
        densityOverFlag = false;
        unsigned int overCount = __density_Predictive();
        if(overCount > 0)
            densityOverFlag = true;
        /////////////////////////////////////////////


//ACTUAL        __sortFluidParticles(ACTUAL_POS);

        /////////////////////////////////////////////
        //corrective
        __corrective();

        /////////////////////////////////////////////



        std::cout << iteration+1<<"th iteration " << "larger part : " << overCount << std::endl;
        iteration++;

        /*
        for (int i = 0; i < numFluidParticles; i++) {
            if (fluidParticles.particleID[i]==0) {
                std::cout << "ID " << fluidParticles.particleID[i] << " particles profile" << std::endl;
                std::cout << "density : " << fluidParticles.density[i] << std::endl;
                std::cout << "pressure : " << fluidParticles.pressure[i] << std::endl;
                std::cout << "pos : ";
                printVector(fluidParticles.pos[i]);
                std::cout << "vel : ";
                printVector(fluidParticles.vel[i]);
                std::cout << "px : ";
                printVector(fluidParticles.px[i]);
                std::cout << "pv : ";
                printVector(fluidParticles.pv[i]);
                std::cout << "FORCE_P : ";
                printVector(fluidParticles.force_p[i]);
                std::cout << "FORCE_NET : ";
                printVector(fluidParticles.force_p[i] + fluidParticles.force_g);
            }
        }
        */
    }

    for (unsigned int i = 0; i < numFluidParticles; i++) {
        fluidParticles.vel[i] = fluidParticles.pv[i];
        fluidParticles.pos[i] = fluidParticles.px[i];
    }

//    system("PAUSE");

//    std::cout << fluidParticles.force_p[0].x<< " : " << fluidParticles.force_p[0].y << " : " << fluidParticles.force_p[0].z << std::endl;
}

float PCISPH::calcPredDensity(unsigned int idx) {

    const glm::vec3 px = fluidParticles.px[idx];
    const unsigned int cIdx = fluidParticles.cellIdxPred[idx];
    float sum = 0;


    // fluid particles NEIGHBOR
    for ( auto nCellIdx : fluidPredNCellIdx[cIdx]) { 

        unsigned int endIdx = fluidParticlesSorter.particleCounter[nCellIdx];

        int occIdx = fluidParticlesSorter.occupiedCellIdxInv[nCellIdx] - 1;
        unsigned int startIdx = (occIdx<0)? 0 : fluidParticlesSorter.particleCounter[fluidParticlesSorter.occupiedCellIdx[occIdx]];

        // �� cell�ȿ� �ִ� ��� particle�鿡 ���ؼ�.
        for (unsigned int i = startIdx; i < endIdx; i++) {
            float len = glm::length(px - fluidParticles.px[i]);
            sum += W(len);
        }
    }

    // wall particles NEIGHBOR
    for (auto nCellIdx : wallPredNCellIdx[cIdx]) { // �� �̿� cell�� ���ؼ�

        unsigned int endIdx = wallParticlesSorter.particleCounter[nCellIdx];

        int occIdx = wallParticlesSorter.occupiedCellIdxInv[nCellIdx] - 1;
        unsigned int startIdx = (occIdx < 0) ? 0 : wallParticlesSorter.particleCounter[wallParticlesSorter.occupiedCellIdx[occIdx]];

        // �� cell�ȿ� �ִ� ��� particle�鿡 ���ؼ�.
        for (unsigned int i = startIdx; i < endIdx; i++) {
            float len = glm::length(px - wallParticles.pos[i]);
            sum += W(len) * wallParticles.densityVar[i];
        }

    }
    sum *= particleMass;


    return sum;
}
float PCISPH::calcDensity(unsigned int idx) {

    const glm::vec3 pos = fluidParticles.pos[idx];
    const unsigned int cIdx = fluidParticles.cellIdx[idx];
    float sum = 0;
    
    // fluid particles NEIGHBOR
    for (auto nCellIdx : fluidNCellIdx[cIdx]) { // �� �̿� cell�� ���ؼ�

        unsigned int endIdx = fluidParticlesSorter.particleCounter[nCellIdx];

        int occIdx = fluidParticlesSorter.occupiedCellIdxInv[nCellIdx] - 1;
        unsigned int startIdx = (occIdx < 0) ? 0 : fluidParticlesSorter.particleCounter[fluidParticlesSorter.occupiedCellIdx[occIdx]];

        // �� cell�ȿ� �ִ� ��� particle�鿡 ���ؼ�.
        for (unsigned int i = startIdx; i < endIdx; i++) {

            const float len = glm::length(pos - fluidParticles.pos[i]);
            sum += W(len);
        }
    }
        // wall particles NEIGHBOR
    for (auto nCellIdx : wallNCellIdx[cIdx]) { // �� �̿� cell�� ���ؼ�

        unsigned int endIdx = wallParticlesSorter.particleCounter[nCellIdx];
        int occIdx = wallParticlesSorter.occupiedCellIdxInv[nCellIdx] - 1;
        unsigned int startIdx = (occIdx < 0) ? 0 : wallParticlesSorter.particleCounter[wallParticlesSorter.occupiedCellIdx[occIdx]];

        // �� cell�ȿ� �ִ� ��� particle�鿡 ���ؼ�.
        for (unsigned int i = startIdx; i < endIdx; i++) {
            float len = glm::length(pos - wallParticles.pos[i]);
            sum += W(len) * wallParticles.densityVar[i];
        }
    }
    sum *= particleMass;
    return sum;
}
glm::vec3 PCISPH::forceP(unsigned int idx) {
    /*
    const glm::vec3 pos = fluidParticles.pos[idx];
    const unsigned int cIdx = fluidParticles.cellIdx[idx];
    glm::vec3 netForceP = glm::vec3(0);


    const float dens = fluidParticles.density[idx];
    const float pres = fluidParticles.density[idx];
    

    for (auto nCellIdx : fluidNCellIdx[cIdx]) {

        const unsigned int endIdx = fluidParticlesSorter.particleCounter[nCellIdx];

        const int occIdx = fluidParticlesSorter.occupiedCellIdxInv[nCellIdx] - 1;
        const unsigned int startIdx = (occIdx < 0) ? 0
            : fluidParticlesSorter.particleCounter[fluidParticlesSorter.occupiedCellIdx[occIdx]];

        for (unsigned int i = startIdx; i < endIdx; i++) {

            const glm::vec3 dist = fluidParticles.pos[i] - pos;
            const glm::vec3 grad_Wp = gradW(dist);

            const float nD = fluidParticles.density[i];
            const float nP = fluidParticles.pressure[i];

            const float coef = particleMass * particleMass * (nP / nD / nD + pres / dens / dens);
            //const float coef = particleMass * particleMass * 2.0f * pres / dens / dens;


            netForceP += coef * grad_Wp;
        }
    }
    for (auto nCellIdx : wallNCellIdx[cIdx]) {

        const unsigned int endIdx = wallParticlesSorter.particleCounter[nCellIdx];

        const int occIdx = wallParticlesSorter.occupiedCellIdxInv[nCellIdx] - 1;
        const unsigned int startIdx = (occIdx < 0) ? 0
            : wallParticlesSorter.particleCounter[wallParticlesSorter.occupiedCellIdx[occIdx]];

        for (unsigned int i = startIdx; i < endIdx; i++) {

            const glm::vec3 dist = wallParticles.pos[i] - pos;
            const glm::vec3 grad_Wp = gradW(dist);

            const float nD = wallParticles.density[i];
            const float nP = wallParticles.pressure[i];

            const float coef = particleMass * particleMass * (nP / nD / nD + pres / dens / dens);
            //const float coef = particleMass * particleMass * 2.0f * pres / dens / dens;

            netForceP += coef * grad_Wp;
            }
    }

    /*
    if (glm::length(netForceP) > 100000.f) {

        for (auto nCellIdx : fluidNCellIdx[cIdx]) {

            const unsigned int endIdx = fluidParticlesSorter.particleCounter[nCellIdx];

            const int occIdx = fluidParticlesSorter.occupiedCellIdxInv[nCellIdx] - 1;
            const unsigned int startIdx = (occIdx < 0) ? 0
                : fluidParticlesSorter.particleCounter[fluidParticlesSorter.occupiedCellIdx[occIdx]];

            for (unsigned int i = startIdx; i < endIdx; i++) {

                const glm::vec3 dist = fluidParticles.pos[i] - pos;
                const glm::vec3 grad_Wp = gradW(dist);
                if (glm::length(grad_Wp) > 100000.f)
                std::cout << glm::length(grad_Wp) << std::endl;
            }
        }

    }
        */

    const glm::vec3 pos = fluidParticles.px[idx];
    const unsigned int cIdx = fluidParticles.cellIdxPred[idx];
    glm::vec3 netForceP = glm::vec3(0);


    const float dens = fluidParticles.density[idx];
    const float pres = fluidParticles.pressure[idx];


    for (auto nCellIdx : fluidPredNCellIdx[cIdx]) {

        const unsigned int endIdx = fluidParticlesSorter.particleCounter[nCellIdx];

        const int occIdx = fluidParticlesSorter.occupiedCellIdxInv[nCellIdx] - 1;
        const unsigned int startIdx = (occIdx < 0) ? 0
            : fluidParticlesSorter.particleCounter[fluidParticlesSorter.occupiedCellIdx[occIdx]];

        for (unsigned int i = startIdx; i < endIdx; i++) {

            const glm::vec3 dist = fluidParticles.px[i] - pos;
            const glm::vec3 grad_Wp = gradW(dist);

            const float nD = fluidParticles.density[i];
            const float nP = fluidParticles.pressure[i];

            const float coef = particleMass * particleMass * (nP / nD / nD + pres / dens / dens);
            //const float coef = particleMass * particleMass * 2.0f * pres / dens / d ens;


            netForceP += coef * grad_Wp;
        }
    }

    //////////////////////

    glm::vec3 tempWallForce = glm::vec3(0);

    float smoothing = 1.0f;

    for (auto nCellIdx : wallPredNCellIdx[cIdx]) {

        const unsigned int endIdx = wallParticlesSorter.particleCounter[nCellIdx];

        const int occIdx = wallParticlesSorter.occupiedCellIdxInv[nCellIdx] - 1;
        const unsigned int startIdx = (occIdx < 0) ? 0
            : wallParticlesSorter.particleCounter[wallParticlesSorter.occupiedCellIdx[occIdx]];

        for (unsigned int i = startIdx; i < endIdx; i++) {

            const glm::vec3 dist = wallParticles.px[i] - pos;
            const glm::vec3 grad_Wp = gradW(dist);
            /*
            const float nD = wallParticles.density[i];
            const float nP = wallParticles.pressure[i];

            const float coef = particleMass * particleMass * (nP / nD / nD + pres / dens / dens);
            */
            const float coef = wallParticles.densityVar[i] * wallParticles.densityVar[i] * particleMass * particleMass * pres / dens / dens;

            tempWallForce += coef * grad_Wp;

            //netForceP += coef * grad_Wp;
        }
    }

    //////////////////////////////

    return netForceP + tempWallForce;

    /*

    float smoothing=1.0f;
    unsigned int smoothingCount = 0;
    for (auto nCellIdx : wallPredNCellIdx[cIdx]) {

        const unsigned int endIdx = wallParticlesSorter.particleCounter[nCellIdx];

        const int occIdx = wallParticlesSorter.occupiedCellIdxInv[nCellIdx] - 1;
        const unsigned int startIdx = (occIdx < 0) ? 0
            : wallParticlesSorter.particleCounter[wallParticlesSorter.occupiedCellIdx[occIdx]];

        for (unsigned int i = startIdx; i < endIdx; i++) {

            const glm::vec3 dist = wallParticles.px[i] - pos;
            const glm::vec3 grad_Wp = gradW(dist);
            if (grad_Wp != glm::vec3(0))
                smoothingCount++;
            /*
            const float nD = wallParticles.density[i];
            const float nP = wallParticles.pressure[i];

            const float coef = particleMass * particleMass * (nP / nD / nD + pres / dens / dens);
            *

            const float coef = particleMass * particleMass * pres / dens / dens;

            netForceP += coef * grad_Wp;
        }
    }


    return netForceP;
    */

}
void PCISPH::__initializeFrameStates() {

    for (auto occCellIdxIter = fluidParticlesSorter.occupiedCellIdx.begin();
        occCellIdxIter != fluidParticlesSorter.occupiedCellIdx.end(); occCellIdxIter++) {

        unsigned int endIdx = fluidParticlesSorter.particleCounter[*occCellIdxIter];    // ���� Cell�� fluid particle���� fluidParticleArr�� ���ĵǾ�����. �� ���� ������ ��ġ.
        unsigned int startIdx = (occCellIdxIter == fluidParticlesSorter.occupiedCellIdx.begin()) ? 0 : fluidParticlesSorter.particleCounter[*(occCellIdxIter - 1)]; // ó�� ��ġ

        for (unsigned int idx = startIdx; idx < endIdx; idx++) {
            fluidParticles.pressure[idx] = 0.0f;
            fluidParticles.force_p[idx] = glm::vec3(0);
            fluidParticles.force_ext[idx] = forceExt(idx);
            fluidParticles.force_vis[idx] = forceVis(idx); // TODO �ϴ� SKIP
        }
    }
}


void PCISPH::__pos_vel_Predictive() {

    for (unsigned int i = 0; i < numFluidParticles; i++) {

        /*
        fluidParticles.pv[i] = fluidParticles.vel[i] + deltaTime * ( (fluidParticles.force_p[i]+ fluidParticles.force_g) / particleMass);
        fluidParticles.px[i] = fluidParticles.pos[i] + deltaTime * fluidParticles.pv[i];
        */
//        fluidParticles.px[i] += deltaTime * deltaTime * ((fluidParticles.force_p[i] + fluidParticles.force_g) / particleMass);

        fluidParticles.pv[i] = fluidParticles.vel[i] + deltaTime * ((fluidParticles.force_p[i] + fluidParticles.force_g) / particleMass);
        fluidParticles.px[i] = fluidParticles.pos[i] + fluidParticles.pv[i] * deltaTime;

        __predOutGridResol(i);
    }
//    system("PAUSE");
}
void PCISPH::__predOutGridResol(unsigned int i) {

    glm::ivec3 cellIdx = __getCellCoord(fluidParticles.px[i]);
    glm::ivec3 flag;

    flag.x = (cellIdx.x == -1) ? -1 : (cellIdx.x >= nGridDivX) ? 1 : 0;
    flag.y = (cellIdx.y == -1) ? -1 : (cellIdx.y >= nGridDivY) ? 1 : 0;
    flag.z = (cellIdx.z == -1) ? -1 : (cellIdx.z >= nGridDivZ) ? 1 : 0;


    if (flag.x == -1) {
        fluidParticles.pv[i].x *= -0.999f;
        fluidParticles.px[i].x += 2*(-fluidParticles.px[i].x);
    }
    if (flag.x == 1) {
        fluidParticles.pv[i].x *= -0.999f;
        fluidParticles.px[i].x += 2 * (boundaryX - fluidParticles.px[i].x);
    }

    if (flag.y == -1) {
        fluidParticles.pv[i].y *= -0.999f;
        fluidParticles.px[i].y += 2 * (-fluidParticles.px[i].y);
    }
    if (flag.y == 1) {
        fluidParticles.pv[i].y *= -0.999f;
        fluidParticles.px[i].y += 2 * (boundaryY - fluidParticles.px[i].y);
    }

    if (flag.z == -1) {
        fluidParticles.pv[i].z *= -0.999f;
        fluidParticles.px[i].z += 2 * (-fluidParticles.px[i].z);
    }
    if (flag.z == 1) {
        fluidParticles.pv[i].z *= -0.999f;
        fluidParticles.px[i].z += 2 * (boundaryZ - fluidParticles.px[i].z);
    }

    // after push and reflection if there is outBoundary particle, that means simulation explode;
    cellIdx = __getCellCoord(fluidParticles.px[i]);
    flag.x = (cellIdx.x == -1) ? -1 : (cellIdx.x >= nGridDivX) ? 1 : 0;
    flag.y = (cellIdx.y == -1) ? -1 : (cellIdx.y >= nGridDivY) ? 1 : 0;
    flag.z = (cellIdx.z == -1) ? -1 : (cellIdx.z >= nGridDivZ) ? 1 : 0;

    if (flag.x != 0 || flag.y != 0 || flag.z != 0) {

        std::cout << fluidParticles.particleID[i]<<"this partciel EXPLODE" << std::endl;
        std::cout << fluidParticles.px[i].x << " : " << fluidParticles.px[i].y << " : " << fluidParticles.px[i].z << std::endl;
        std::cout << "PCISPH::__PREDOUTGRIDRESOL() SIMULATION EXPLODE!!! : "<<i <<" particle"<< std::endl;
        system("PAUSE");
    }
}

unsigned int PCISPH::__density_Predictive() {
    // ���� px, predIdx ����. fluidParticleArr ���� �� ����.
    unsigned int overThresholdCount = 0;

    for (auto occCellIdxIter = fluidParticlesSorter.occupiedCellIdx.begin();
        occCellIdxIter != fluidParticlesSorter.occupiedCellIdx.end(); occCellIdxIter++) {

        unsigned int endIdx = fluidParticlesSorter.particleCounter[*occCellIdxIter];    // ���� Cell�� fluid particle���� fluidParticleArr�� ���ĵǾ�����. �� ���� ������ ��ġ.
        unsigned int startIdx = (occCellIdxIter == fluidParticlesSorter.occupiedCellIdx.begin()) ? 0 : fluidParticlesSorter.particleCounter[*(occCellIdxIter - 1)]; // ó�� ��ġ

        for (unsigned int idx = startIdx; idx < endIdx; idx++) {

            // 1��
            fluidParticles.density[idx] = calcPredDensity(idx);
            fluidParticles.densityVar[idx] = (fluidParticles.density[idx] - restDensity);
            
//            fluidParticles.densityVar[idx] = (calcPredDensity(idx) - restDensity);//2��

            if (fluidParticles.densityVar[idx] > eta * restDensity) {
                overThresholdCount++;
            }
            
        }
    }

    return overThresholdCount;
}
void PCISPH::__corrective() {

    /// actual density   
    /// delta for fluidParticles only
    for (auto occCellIdxIter = fluidParticlesSorter.occupiedCellIdx.begin();
        occCellIdxIter != fluidParticlesSorter.occupiedCellIdx.end(); occCellIdxIter++) {

        unsigned int endIdx = fluidParticlesSorter.particleCounter[*occCellIdxIter];    // ���� Cell�� fluid particle���� fluidParticleArr�� ���ĵǾ�����. �� ���� ������ ��ġ.
        unsigned int startIdx = (occCellIdxIter == fluidParticlesSorter.occupiedCellIdx.begin()) ? 0 : fluidParticlesSorter.particleCounter[*(occCellIdxIter - 1)]; // ó�� ��ġ

        for (unsigned int idx = startIdx; idx < endIdx; idx++) {


            if (fluidParticles.densityVar[idx] < FLT_EPSILON )
                fluidParticles.pressure[idx] = 0.0f;
            else {


                float smoothingConst = (1.0f - std::abs(fluidParticles.densityVar[idx]) / restDensity);


                //smoothingConst = std::pow(smoothingConst, 10.0f);

                smoothingConst = 0.7f * std::pow(smoothingConst, 10.0f) + 0.3f;

//                smoothingConst = (1.0f - fluidParticles.densityVar[idx] / restDensity);


                fluidParticles.pressure[idx] += smoothingConst * delta * fluidParticles.densityVar[idx];//1��
                fluidParticles.pressure[idx] = std::max(fluidParticles.pressure[idx], 0.0f);
                
            }

            
        }
    }
     
    /// Pressure Force
    for (auto occCellIdxIter = fluidParticlesSorter.occupiedCellIdx.begin();
        occCellIdxIter != fluidParticlesSorter.occupiedCellIdx.end(); occCellIdxIter++) {

        unsigned int endIdx = fluidParticlesSorter.particleCounter[*occCellIdxIter];    // ���� Cell�� fluid particle���� fluidParticleArr�� ���ĵǾ�����. �� ���� ������ ��ġ.
        unsigned int startIdx = (occCellIdxIter == fluidParticlesSorter.occupiedCellIdx.begin()) ? 0 : fluidParticlesSorter.particleCounter[*(occCellIdxIter - 1)]; // ó�� ��ġ

        for (unsigned int idx = startIdx; idx < endIdx; idx++) {


            float smoothingConst = (1.0f - std::abs(fluidParticles.densityVar[idx]) / restDensity);

            smoothingConst = 0.7f * std::pow(smoothingConst, 10.0f) + 0.3f;

            fluidParticles.force_p[idx] += smoothingConst*forceP(idx);
            
            //fluidParticles.force_p[idx] = forceP(idx);


        }
        
    }

}

float PCISPH::calcDelta(unsigned int idx) {

//    const glm::vec3 pos = fluidParticles.pos[idx];
    const glm::vec3 pos = fluidParticles.pos[idx];
    glm::vec3 sumOfGrad = glm::vec3(0);
    float sum = 0.0f;

    float beta = deltaTime * deltaTime * particleMass * particleMass * 2.0f / restDensity / restDensity;

    //const unsigned int cIdx = fluidParticles.cellIdx[idx];
    const unsigned int cIdx = fluidParticles.cellIdx[idx];

    // fluid particles NEIGHBOR
    for (auto nCellIdx : fluidNCellIdx[cIdx]) { // �� �̿� cell�� ���ؼ�

        const unsigned int endIdx = fluidParticlesSorter.particleCounter[nCellIdx];

        const int occIdx = fluidParticlesSorter.occupiedCellIdxInv[nCellIdx] - 1;
        const unsigned int startIdx = (occIdx < 0) ? 0 : fluidParticlesSorter.particleCounter[fluidParticlesSorter.occupiedCellIdx[occIdx]];

        for (unsigned int i = startIdx; i < endIdx; i++) {

            const glm::vec3 grad = fluidParticles.pos[i] - pos;
            sumOfGrad += gradW(grad);
            sum += glm::dot(grad, grad);
        }
    }

    for (auto nCellIdx : wallNCellIdx[cIdx]) {

        unsigned int endIdx = wallParticlesSorter.particleCounter[nCellIdx];

        int occIdx = wallParticlesSorter.occupiedCellIdxInv[nCellIdx] - 1;
        unsigned int startIdx = (occIdx < 0) ? 0 : wallParticlesSorter.particleCounter[wallParticlesSorter.occupiedCellIdx[occIdx]];

        // �� cell�ȿ� �ִ� ��� particle�鿡 ���ؼ�.
        for (unsigned int i = startIdx; i < endIdx; i++) {

            const glm::vec3 grad = wallParticles.pos[i] - pos;

            sumOfGrad += gradW(grad);
            sum += glm::dot(grad, grad);

        }
    }

    if (sum <= FLT_EPSILON && sum >= -FLT_EPSILON)
        return 0;

    sum += glm::dot(sumOfGrad, sumOfGrad);
    beta *= sum;
    beta = 1.0f / beta;


    return beta;
}



void PCISPH::__sortFluidParticles(SORTMODE mode) {

    if (mode == ACTUAL_POS) {
        for (unsigned int i = 0; i < numFluidParticles; i++) {
            glm::ivec3 cellCoord = __getCellCoord(fluidParticles.pos[i]);
            fluidParticles.cellIdx[i] = mortonEncode(cellCoord.x, cellCoord.y, cellCoord.z);
        }
        fluidParticlesSorter.sortby(ACTUAL_POS);
    }
    else if (mode == PRED_POS) {
        for (unsigned int i = 0; i < numFluidParticles; i++) {
            glm::ivec3 cellCoord = __getCellCoord(fluidParticles.px[i]);
            fluidParticles.cellIdxPred[i] = mortonEncode(cellCoord.x, cellCoord.y, cellCoord.z);
        }
        fluidParticlesSorter.sortby(PRED_POS);
    }

}
void PCISPH::__makeNCellMap() {
    //�� cell�� ���ؼ� neighbor cell�� particle�� �ִ� cell idx�� vector�� �����Ѵ�. �׸��� �װ��� ���� cellIdx�� mapping ��.

    fluidNCellIdx.clear();  //std::unordered_map<unsigned int, std::vector<unsigned int>> fluidNCellIdx;
    wallNCellIdx.clear();   //std::unordered_map<unsigned int, std::vector<unsigned int>> wallNCellIdx;

    for (auto occCellIdxIter = fluidParticlesSorter.occupiedCellIdx.begin();
        occCellIdxIter != fluidParticlesSorter.occupiedCellIdx.end(); occCellIdxIter++) {

        std::vector<unsigned int> fluidNCell;
        std::vector<unsigned int> wallNCell;
        __nCellSearch(&fluidNCell, &wallNCell, mortonDecode(*occCellIdxIter));

        fluidNCellIdx[*occCellIdxIter] = fluidNCell;
        wallNCellIdx[*occCellIdxIter] = wallNCell;
        /* TEST
        std::cout << "PCISPH::update() test neighborSearch" << std::endl;
        std::cout << "current cellIdx : " << occCellIdx << " coord" << mortonDecode(occCellIdx).x << " : " << mortonDecode(occCellIdx).y << " : " << mortonDecode(occCellIdx).z << std::endl;
        std::cout << std::endl << "fluidNCell" << std::endl;
        for (auto ele : fluidNCell)
            std::cout << ele << " : " << mortonDecode(ele).x << " : " << mortonDecode(ele).y << " : " << mortonDecode(ele).z << std::endl;
        std::cout << std::endl << "wallNCell" << std::endl;
        for (auto ele : wallNCell)
            std::cout << ele << " : " << mortonDecode(ele).x << " : " << mortonDecode(ele).y << " : " << mortonDecode(ele).z << std::endl;
        system("PAUSE");
        */
    }
}
void PCISPH::__nCellSearch(std::vector<unsigned int>* fluidNIdx, std::vector<unsigned int>* wallNIdx, glm::ivec3 center) {
    // ���� cell���� �ֺ� cell�� particle�� �ִ��� ����.
    // ������ vector�� �־��ٰǵ�.
    // �ֺ� cell�� �ִ� particle�� wall ����, fluid���� ���� �־��ٰ���.
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            for (int k = -1; k <= 1; k++) {
                glm::ivec3 steper = center + glm::ivec3(i, j, k);
                int nCell = mortonEncode(steper.x, steper.y, steper.z);

                if (nCell >= 0 && fluidParticlesSorter.particleCounter.count(nCell))
                    fluidNIdx->push_back(nCell);

                if (nCell >= 0 && wallParticlesSorter.particleCounter.count(nCell))
                    wallNIdx->push_back(nCell);


            }
        }
    }
}
void PCISPH::__makePredNCellMap() {
    //�� cell�� ���ؼ� neighbor cell�� particle�� �ִ� cell idx�� vector�� �����Ѵ�. �׸��� �װ��� ���� cellIdx�� mapping ��.

    fluidPredNCellIdx.clear();  //std::unordered_map<unsigned int, std::vector<unsigned int>> fluidNCellIdx;
    wallPredNCellIdx.clear();   //std::unordered_map<unsigned int, std::vector<unsigned int>> wallNCellIdx;

    for (auto occCellIdxIter = fluidParticlesSorter.occupiedCellIdx.begin();
        occCellIdxIter != fluidParticlesSorter.occupiedCellIdx.end(); occCellIdxIter++) {

        std::vector<unsigned int> fluidNCell;
        std::vector<unsigned int> wallNCell;
        __nCellSearch(&fluidNCell, &wallNCell, mortonDecode(*occCellIdxIter));

        fluidPredNCellIdx[*occCellIdxIter] = fluidNCell;
        wallPredNCellIdx[*occCellIdxIter] = wallNCell;
    }
}

void PCISPH::__calcDelta() {

    glm::vec3 sumGradW = glm::vec3(0);
    float sumGradW2 = 0.0f;

    glm::vec3 xj = glm::vec3(-coreRad, -coreRad, -coreRad);
    const glm::vec3 xi(0, 0, 0);
    while (xj[0] <= coreRad + 10e-7)
    {
        while (xj[1] <= coreRad + 10e-7)
        {
            while (xj[2] <= coreRad + 10e-7)
            {
                // check if xj is in the support of xi
                if (glm::length(xi - xj) < coreRad)
                {
                    const glm::vec3 grad = gradW(xi - xj);
                    sumGradW += grad;
                    sumGradW2 += glm::dot(grad, grad);
                }
                xj[2] += radius *2.0f;
            }
            xj[1] += radius * 2.0f;
            xj[2] = -coreRad;
        }
        xj[0] += radius * 2.0f;
        xj[1] = -coreRad;
        xj[2] = -coreRad;
    }

    const float beta = deltaTime * deltaTime * particleMass * particleMass * 2.0f / restDensity / restDensity;

    delta = 1.0f / (beta * sumGradW2); // �����ϰ� ������ sumGradW = 0��.
    std::cout << ((glm::length(sumGradW) * glm::length(sumGradW))) << " : " << sumGradW2  << " : " << delta << std::endl;

    system("PAUSE");


    return;

    /*
    // delta���� �̿� particle�� 2�� �̻� �������� delta ���� ����ϰ� ����.
    // �ٵ� ������ �и� 0�Ǹ鼭 ����. 1���� �������� 4������ ũ�� ����.

    
    glm::vec3 asdf[8] = {
        glm::vec3(-radius, -radius, -radius),glm::vec3(-radius, -radius, radius),
        glm::vec3(radius, -radius, radius),glm::vec3(radius, -radius, -radius),
        glm::vec3(-radius, radius, -radius),glm::vec3(-radius, radius, radius),
        glm::vec3(radius, radius, radius),glm::vec3(radius, radius, -radius),
    };

    sumGradW = glm::vec3(0);
    sumGradW2 = 0.0f;

    for (int i = 0; i < 4; i++) {     

        if (i == 7 || i == 6 || i == 4 || i == 3)
            continue;

        xj = asdf[i];
        const glm::vec3 grad = gradW(glm::vec3(0) - xj);
        sumGradW += grad;
        sumGradW2 += glm::dot(grad, grad);

        std::cout << glm::dot(grad, grad) << std::endl;
    }
    
    delta = 1.0f / (beta*(glm::length(sumGradW) * glm::length(sumGradW) + sumGradW2)); // �����ϰ� ������ sumGradW = 0��.

    std::cout<< (glm::length(sumGradW) * glm::length(sumGradW))<< " : " << sumGradW2 <<  " : " << delta << std::endl;

    system("PAUSE");
    */
}


//==================================================================================
//==================================================================================
//==================================================================================
//==================================================================================
//==================================================================================
//==================================================================================
//==================================================================================
//==================================================================================
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

    __buildCounter(); // offset�� ����� ����. cellIdx�� �ִ� particle�� ��ִ��� ��.

    unsigned int* sortedIdx = __buildSortedIdx();

    /*  TO TEST
    if(arrT == FLUID)
    for(int i=0; i<count; i++)
        std::cout << sortedIdx[i] << " : " << mortonDecode(particleArr->cellIdx[sortedIdx[i]]).x * 0.125f << " : " << mortonDecode(particleArr->cellIdx[sortedIdx[i]]).y * 0.125f << " : " << mortonDecode(particleArr->cellIdx[sortedIdx[i]]).z * 0.125f << std::endl;
    system("PAUSE ");
    */

    //sortedIdx�� �̿��ؼ� ���� �ؾ� ��.

    ParticlesArray tempArr;

    __particlesArrInit(&tempArr, count, arrT);

    __swapAndSort(&tempArr,sortedIdx);


    delete[] sortedIdx;

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

    for (auto iter = occupiedCellIdx.begin() + 1; iter != occupiedCellIdx.end(); iter++) {
        particleCounter[*iter] += particleCounter[*(iter - 1)];
    }

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
                std::cout <<"�� cell�� 125���� �ƴѰ���."<< element.first<< " : " << element.second << std::endl;
                std::cout << mortonDecode(element.first).x * 0.125f<<" : " << mortonDecode(element.first).y * 0.125f << " : " << mortonDecode(element.first).z * 0.125f << std::endl;
                // TODOTEST demorton �ؼ� origin �����. �ű� ���� particle�� ��ǥ ����. ���° particle������.
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
        particleArr->particleID = new unsigned int[numbers];
    }
}
void PCISPH::Z_Sort::__swapAndSort(ParticlesArray* temp, unsigned int* sortedIdx) {
    // temp�� �����ؼ� �����.

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


            temp->particleID[i] = particleArr->particleID[curIdx];
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
        delete[] particleArr->particleID;


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

        particleArr->particleID = temp->particleID;
    }


    // ������ �־��� �� �����ϱ�.
}



//=========================================== integrate �ʿ� ����.
/*

// fluidSorter�� �̿��ؼ� ������ �� ���� ��ü particle�� ����� ���� �� �Լ��� ȣ����.
void PCISPH::Z_Sort::integrateSortedArr(Z_Sort wallSorter, Z_Sort fluidSorter, SORTMODE mode) {

    currentMode = mode;

    __integrateCounter(wallSorter, fluidSorter);

    int* sortedIdx = __integrateSortedIdx(wallSorter, fluidSorter);

    __integrate(wallSorter, fluidSorter, sortedIdx);


    delete[] sortedIdx;
}
void PCISPH::Z_Sort::__integrateCounter(Z_Sort wallSorter, Z_Sort fluidSorter) {
    //wallParticle�� fluidParticle�� hashMap��ģ��
    particleCounter = wallSorter.particleCounter;

    // mode�� ���� cellIdx�Ǵ� cellIdxPred�� �ִ� particle�� ������ ��.
    for (auto ele : fluidSorter.particleCounter) {
        if (particleCounter.count(ele.first)) {
            particleCounter[ele.first] += ele.second;
        }
        else
            particleCounter[ele.first] = ele.second;
    }

    //occupiedCellIdx �����.
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
        for (unsigned int i = 1; i < fluidSorter.count + 1; i++) {
            unsigned int idx = occupiedCellIdxInv[fluidSorter.particleArr->cellIdx[i]];
            sortParticleIdx[idx].push_back(-1 * i);
        }
    }
    else if (currentMode == PRED_POS) {
        for (unsigned int i = 0; i < wallSorter.count; i++) {
            unsigned int idx = occupiedCellIdxInv[wallSorter.particleArr->cellIdxPred[i]];
            sortParticleIdx[idx].push_back(i);
        }
        for (unsigned int i = 1; i < fluidSorter.count + 1; i++) {
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
            curIdx = -curIdx - 1;

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

*/


void PCISPH::nSearchTest() {

    for (int i = 0; i < 8; i++) {
        fluidParticles.pos[i] = glm::vec3(1.0f, 1.0f, 1.0f) + float(i) * glm::vec3(coreRad, 0.0f, 0.0f);
        fluidParticles.pressure[i] = i * 100.0f;
    }
    for (int i = 0; i < 8; i++)
        fluidParticles.px[i] = glm::vec3(1.0f, 1.0f, 1.0f) - float(i) * glm::vec3(coreRad,0.0f, 0.0f);

    __sortFluidParticles(ACTUAL_POS);
    __makeNCellMap();
    
    // i �� neighbor�� pressure�� �� ����غ�.
    
    for (int i = 0; i < 8; i++) {
        const glm::vec3 pos = fluidParticles.pos[i];
        const unsigned int cIdx = fluidParticles.cellIdx[i];


        std::cout << i << " th particles ( P ) : " << fluidParticles.pressure[i] << " (cellIdx): " << fluidParticles.cellIdx[i] << " : " << mortonDecode(cIdx).x << " : " << mortonDecode(cIdx).y << " : " << mortonDecode(cIdx).z << std::endl;

        for (auto nCellIdx : fluidNCellIdx[cIdx]) {

            unsigned int endIdx = fluidParticlesSorter.particleCounter[nCellIdx];

            int occIdx = fluidParticlesSorter.occupiedCellIdxInv[nCellIdx] - 1;
            unsigned int startIdx = (occIdx < 0) ? 0
                : fluidParticlesSorter.particleCounter[fluidParticlesSorter.occupiedCellIdx[occIdx]];

            for (unsigned int i = startIdx; i < endIdx; i++) {
                unsigned int c = fluidParticles.cellIdx[i];
                std::cout<<"(P) : "<<fluidParticles.pressure[i]<< " ||| cellCoord : " << mortonDecode(c).x << " : " << mortonDecode(c).y << " : " << mortonDecode(c).z << std::endl;
            }
        }
        std::cout << std::endl;

    }
    std::cout << "================================================================================ PREDPREDPREDPREDPREDPREDPREDPRED" << std::endl;
    
    __sortFluidParticles(PRED_POS);
    __makePredNCellMap();

    for (int i = 0; i < 8; i++) {
        const glm::vec3 pos = fluidParticles.px[i];
        const unsigned int cIdx = fluidParticles.cellIdxPred[i];

        std::cout << i << " th particles ( P ) : " << fluidParticles.pressure[i] << " (cellIdx): " << fluidParticles.cellIdxPred[i] << " : " << mortonDecode(cIdx).x << " : " << mortonDecode(cIdx).y << " : " << mortonDecode(cIdx).z << std::endl;
        printVector(fluidParticles.px[i]);

        for (auto nCellIdx : fluidPredNCellIdx[cIdx]) {

            unsigned int endIdx = fluidParticlesSorter.particleCounter[nCellIdx];

            int occIdx = fluidParticlesSorter.occupiedCellIdxInv[nCellIdx] - 1;
            unsigned int startIdx = (occIdx < 0) ? 0
                : fluidParticlesSorter.particleCounter[fluidParticlesSorter.occupiedCellIdx[occIdx]];

            for (unsigned int i = startIdx; i < endIdx; i++) {
                unsigned int c = fluidParticles.cellIdxPred[i];
                std::cout << "(P) : " << fluidParticles.pressure[i] << " ||| cellCoord : " << mortonDecode(c).x << " : " << mortonDecode(c).y << " : " << mortonDecode(c).z << std::endl;
            }
        }
        std::cout << std::endl;
    }
    std::cout << "================================================================================" << std::endl;


    __sortFluidParticles(ACTUAL_POS);
    
    for (int i = 0; i < 8; i++) {
        const glm::vec3 pos = fluidParticles.pos[i];
        const unsigned int cIdx = fluidParticles.cellIdx[i];

        std::cout << i << " th particles ( P ) : " << fluidParticles.pressure[i] << " (cellIdx): " << fluidParticles.cellIdx[i] << " : " << mortonDecode(cIdx).x << " : " << mortonDecode(cIdx).y << " : " << mortonDecode(cIdx).z << std::endl;

        for (auto nCellIdx : fluidNCellIdx[cIdx]) {

            unsigned int endIdx = fluidParticlesSorter.particleCounter[nCellIdx];

            int occIdx = fluidParticlesSorter.occupiedCellIdxInv[nCellIdx] - 1;
            unsigned int startIdx = (occIdx < 0) ? 0
                : fluidParticlesSorter.particleCounter[fluidParticlesSorter.occupiedCellIdx[occIdx]];

            for (unsigned int i = startIdx; i < endIdx; i++) {
                unsigned int c = fluidParticles.cellIdx[i];
                std::cout << "(P) : " << fluidParticles.pressure[i] << " ||| cellCoord : " << mortonDecode(c).x << " : " << mortonDecode(c).y << " : " << mortonDecode(c).z << std::endl;
            }
        }
        std::cout << std::endl;

    }
    
    system("PAUSE");
}
