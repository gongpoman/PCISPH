#pragma once

#include<vector>
#include<unordered_map>

#include<glm/glm.hpp>


extern const float deltaTime;

enum ARRTYPE {
	WALL,
	FLUID,
	ALL
};

enum SORTMODE {
	ACTUAL_POS,
	PRED_POS
};

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

struct ParticlesArray {
	unsigned int count;	// ���� � particle �ִ���.
	unsigned int numMaxParticles; //�ִ� particle ����.

	bool *isWall;	// �̰� wall�� particle����.

	int *cellIdx; // actual position ���� ��� cell�� �����ִ��� morton code�� ��Ÿ�� ��.
	glm::vec3 *vel;
	glm::vec3 *pos;

	int *cellIdxPred; //  pred position ���� ��� cell�� �����ִ��� morton code�� ��Ÿ�� ��.
	glm::vec3 *pv;
	glm::vec3 *px;

	glm::vec3 *force_p;
	glm::vec3 *force_vis;
	glm::vec3 *force_ext;

	glm::vec3 force_g;

	float *pressure;
	float *density;
	float *densityVar;

	/* dependent...
	void __particlesArrInit(ParticlesArray*, unsigned int numbers);
	void __appendParticle(ParticlesArray*, Particle);
	*/
};


class PCISPH {
public:

	class Z_Sort {
	public :

		ARRTYPE arrT;
		SORTMODE currentMode;

		unsigned int maxCellIdx;

		int count;
		
		std::vector<unsigned int> occupiedCellIdx; // ordered
		std::unordered_map<unsigned int, unsigned int> occupiedCellIdxInv;
		std::unordered_map<unsigned int, unsigned int> particleCounter;

		ParticlesArray* particleArr; 

		Z_Sort();
		Z_Sort(int, ParticlesArray*,ARRTYPE,unsigned int);

		void sortby(SORTMODE mode = ACTUAL_POS); // particleArr�� /particle Counter�� offset���� �Ͽ�/ mode�� �������� ����.

		void integrateSortedArr(Z_Sort wallSorter, Z_Sort fluidSorter, SORTMODE mode); // ���ĵ� particlesArr ��ġ�°� 
	private:
		void __integrateCounter(Z_Sort wallSorter, Z_Sort fluidSorter);
		int* __integrateSortedIdx(Z_Sort wallSorter, Z_Sort fluidSorter);
		void __integrate(Z_Sort wallSorter, Z_Sort fluidSorter, int*);

		void __buildCounter();
		unsigned int* __buildSortedIdx();
		void __particlesArrInit(ParticlesArray* particleArr, unsigned int numbers, ARRTYPE arrtype);
		void __swapAndSort(ParticlesArray* temp,unsigned int*);

	};

public:


	unsigned int numParticles;
	unsigned int numFluidParticles;
	unsigned int numWallParticles;

	unsigned int numDrawParticles;
	bool drawWall;
	

	ParticlesArray particles;	// �̰� �ణ temporary������ ����. fulidParticle ������ ���Ҷ��� �ʿ���?

	ParticlesArray fluidParticles;
	ParticlesArray wallParticles;

	Z_Sort particlesSorter;
	Z_Sort fluidParticlesSorter;
	Z_Sort wallParticlesSorter;

	//particle constant
	float restDensity;
	float radius;
	float particleMass;

	//grid constant
	float boundaryX;
	float boundaryY;
	float boundaryZ;

	float coreRad;

	//std::vector<glm::ivec3>* neighborIdices;

	float eta;
	unsigned int minIter;
	unsigned int maxIter;

	PCISPH(glm::ivec3 numOfFluidParticles,glm::vec3 boundarySideLen, bool drawAll = true);
	~PCISPH();

	void update();


private :

	std::vector<unsigned int> wallCellIdx;

	std::vector<Particle> brick;


	//particle constant
	unsigned int numFluidParticlesX; // for particle initial position
	unsigned int numFluidParticlesY;
	unsigned int numFluidParticlesZ;
	float sideLenX;					// for particle initial position
	float sideLenY;
	float sideLenZ;


	//grid constant
	unsigned int nGridDivX;
	unsigned int nGridDivY;
	unsigned int nGridDivZ;
	unsigned int upperMaxGridDiv;


	void __sceneInitialize();
	void __brickInit();
	void __gridInit();
	void __sorterInit(); // Z_sort ����. particles�̰� ����ִµ� �̰� build�ϱ�.

	void __addWallParticles(glm::vec3);
	void __particleInit(int mode);

	glm::ivec3 __getCellCoord(glm::vec3);
	glm::vec3 __gridLocalPos(glm::vec3 pos);




	void __particlesArrInit(ParticlesArray*,unsigned int numbers,ARRTYPE arrtype);
	void __appendParticle(ParticlesArray*, Particle, ARRTYPE arrtype);
	void __deleteParticleArr(ParticlesArray*, ARRTYPE arrtype);





};


