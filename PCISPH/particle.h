#pragma once

#include<vector>

#include<glm/glm.hpp>


class Particle {
public:
	bool isWall;

	unsigned int cellIdx;
	glm::vec3 vel;
	glm::vec3 pos;

	unsigned int cellIdxPred;
	glm::vec3 pv;
	glm::vec3 px;

	glm::vec3 force_p;
	glm::vec3 force_vis;
	glm::vec3 force_ext;

	float pressure;
	float density;
	float densityVar;

	Particle() {
		cellIdx = -1;

		isWall = false;
		vel = glm::vec3(0);
		pos = glm::vec3(0);

		pv = glm::vec3(0);
		px = glm::vec3(0);

		force_p = glm::vec3(0);
		force_vis = glm::vec3(0);
		force_ext = glm::vec3(0);

		pressure = 0.0f;
		density = 0.0f;
		densityVar = 0.0f;
	}
};


class PCISPH {
public:
	unsigned int numParticles;
	unsigned int numFluidParticles;
	unsigned int numWallParticles;

	unsigned int numDrawParticles;
	bool drawWall;
	

	Particle *particles;

	Particle *fluidParticles;
	Particle *wallParticles;

	float restDensity;
	float radius;
	float particleMass;
	float coreRad;

	unsigned int numWaterParticlesX;
	unsigned int numWaterParticlesY;
	unsigned int numWaterParticlesZ;

	float sideLenX;
	float sideLenY;
	float sideLenZ;
	
	float boundaryX;
	float boundaryY;
	float boundaryZ;

	unsigned int nGridDivX;
	unsigned int nGridDivY;
	unsigned int nGridDivZ;
	unsigned int upperMaxGridDiv;

	//std::vector<glm::ivec3>* neighborIdices;

	float eta;
	unsigned int minIter;
	unsigned int maxIter;

	PCISPH(unsigned int numX = 20, unsigned int numY =20, unsigned int numZ=20, bool drawAll = true);
	~PCISPH();

private :

	std::vector<Particle> brick;


	void __sceneInitialize();
	void __brickInit();
	void __gridInit();
	void __addWallParticles(glm::vec3, unsigned int);
	void __particleInit(int mode);

	glm::ivec3 __getCellCoord(glm::vec3);

	glm::vec3 __gridLocalPos(glm::vec3 pos);

};


