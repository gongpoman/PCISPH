#include<iostream>
#include<vector>

#include "myutil.h"
#include "particle.h"

PCISPH::PCISPH(unsigned int numX, unsigned int numY, unsigned int numZ,bool drawAll) {

    drawWall = drawAll;

	restDensity = 997.0f;
	radius = 0.025f;
	particleMass = restDensity * 4.0f / 3.0f * 3.141592f * radius * radius * radius;
	coreRad = radius * 5.0f;


    numWaterParticlesX = numX;
    numWaterParticlesY = numY;
    numWaterParticlesZ = numZ;
    numFluidParticles = numWaterParticlesX * numWaterParticlesY * numWaterParticlesZ;

	boundaryX = 3.0f;
	boundaryY = 1.5f;
	boundaryZ = 3.0f;

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

	__sceneInitialize();
}
void PCISPH::__sceneInitialize() {
	__brickInit();
	__gridInit();
    __particleInit(0);
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

    __particlesArrInit(&wallParticles,numWallParticles);
//  predParticles = new Particle[numParticle];

    for (auto iter = boundaryOrigins.begin(); iter != boundaryOrigins.end(); iter++) {
        __addWallParticles(*iter);
    }

}

void PCISPH::__particleInit(int mode) {

    __particlesArrInit(&particles, numParticles);
    __particlesArrInit(&fluidParticles, numFluidParticles);

    //TODO================================================================================================================================================================================================================================================
    if (mode == 0) {
        for (unsigned int i = 0; i < numWaterParticlesX; i++) {
            for (unsigned int j = 0; j < numWaterParticlesY; j++) {
                for (unsigned int k = 0; k < numWaterParticlesZ; k++) {

                    fluidParticles[i * numRow * numRow + numRow * j + k].pos = glm::vec3(sideLen / 2.0f - sideLen * (i / (float)numRow)
                        , sideLen - sideLen * (j / (float)numRow) + coreRad
                        , sideLen / 2.0f - sideLen * (k / (float)numRow));
                    particles[i * numRow * numRow + numRow * j + k].density = restDensity;

                }
            }
        }
    }

    if (mode == 1) {

        for (unsigned int i = 0; i < numRow; i++) {
            for (unsigned int j = 0; j < numRow; j++) {
                for (unsigned int k = 0; k < numRow; k++) {

                    particles[i * numRow * numRow + numRow * j + k].pos = glm::vec3(sideLen * (i / (float)numRow) - cubicBoundaryLen / 2.0f + coreRad
                        , sideLen - sideLen * (j / (float)numRow) + coreRad
                        , sideLen * (k / (float)numRow) - cubicBoundaryLen / 2.0f + coreRad);
                    particles[i * numRow * numRow + numRow * j + k].density = restDensity;
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

        __appendParticle(&wallParticles,tempBrick[i]);
    }

}
glm::ivec3 PCISPH::__getCellCoord(glm::vec3 position) {
    float x = position.x / boundaryX * nGridDivX;
    float y = position.y / boundaryY * nGridDivY;
    float z = position.z / boundaryZ * nGridDivZ;

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

void PCISPH::__particlesArrInit(ParticlesArray* particleArr,unsigned int numbers) {

    particleArr->count = 0;
    particleArr->numMaxParticles = numbers;

    particleArr->isWall = new bool[numbers];

    particleArr->cellIdx = new unsigned int[numbers];
    particleArr->cellIdxPred = new unsigned int[numbers];

    particleArr->pressure = new float[numbers];
    particleArr->density = new float[numbers];
    particleArr->densityVar = new float[numbers];

    particleArr->vel= new glm::vec3[numbers];
    particleArr->pos= new glm::vec3[numbers];
    particleArr->pv= new glm::vec3[numbers];
    particleArr->px= new glm::vec3[numbers];
    
    particleArr->force_p= new glm::vec3[numbers];
    particleArr->force_ext= new glm::vec3[numbers];
    particleArr->force_vis= new glm::vec3[numbers];

    particleArr->force_g = glm::vec3(0, -9.8f, 0) * particleMass;

}
void PCISPH::__appendParticle(ParticlesArray* particleArr, Particle particle) {

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

    particleArr->vel[idx] = particle.vel;
    particleArr->pos[idx] = particle.pos;
    particleArr->pv[idx] = particle.pv;
    particleArr->px[idx] = particle.px;

    particleArr->force_p[idx] = particle.force_p;
    particleArr->force_ext[idx] = particle.force_ext;
    particleArr->force_vis[idx] = particle.force_vis;

    particleArr->count++;
}