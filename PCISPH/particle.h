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

	unsigned int pID;

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

		pID = 0;
	}
};

struct ParticlesArray {
	unsigned int count;	// 현재 몇개 particle 있는지.
	unsigned int numMaxParticles; //최대 particle 갯수.

	unsigned int *particleID;

	bool *isWall;	// 이게 wall의 particle인지.

	int *cellIdx; // actual position 기준 어느 cell에 속해있는지 morton code로 나타낸 것.
	glm::vec3 *vel;
	glm::vec3 *pos;

	int *cellIdxPred; //  pred position 기준 어느 cell에 속해있는지 morton code로 나타낸 것.
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

		void sortby(SORTMODE mode = ACTUAL_POS); // particleArr을 /particle Counter를 offset으로 하여/ mode를 기준으로 정렬.

//		void integrateSortedArr(Z_Sort wallSorter, Z_Sort fluidSorter, SORTMODE mode);
	private:
		void __buildCounter();
		unsigned int* __buildSortedIdx();
		void __particlesArrInit(ParticlesArray* particleArr, unsigned int numbers, ARRTYPE arrtype);
		void __swapAndSort(ParticlesArray* temp,unsigned int*);

		//		void __integrateCounter(Z_Sort wallSorter, Z_Sort fluidSorter);
//		int* __integrateSortedIdx(Z_Sort wallSorter, Z_Sort fluidSorter);
//		void __integrate(Z_Sort wallSorter, Z_Sort fluidSorter, int*);

	};

public:

	unsigned int numParticles;
	unsigned int numFluidParticles;
	unsigned int numWallParticles;

	unsigned int numDrawParticles;
	bool drawWall;
	

	ParticlesArray particles;	// 이거 약간 temporary느낌이 있음. fulidParticle 물리량 구할때만 필요한?

	ParticlesArray fluidParticles;
	ParticlesArray wallParticles;

//	Z_Sort particlesSorter;
	Z_Sort fluidParticlesSorter;
	Z_Sort wallParticlesSorter;

	//particle constant
	float restDensity;
	float radius;
	float particleMass;
	float wallMassRatio;
	int numBrickRow;

	//grid constant
	float boundaryX;
	float boundaryY;
	float boundaryZ;


	//grid constant
	unsigned int nGridDivX;
	unsigned int nGridDivY;
	unsigned int nGridDivZ;
	unsigned int upperMaxGridDiv;

	float coreRad;
	float delta;

	//std::vector<glm::ivec3>* neighborIdices;
	float eta;
	unsigned int minIter;
	unsigned int maxIter;

	PCISPH(glm::ivec3 numOfFluidParticles,glm::vec3 boundarySideLen, bool drawAll = true);
	~PCISPH();

	void update();
	void nSearchTest();

private :

	std::unordered_map<unsigned int, std::vector<unsigned int>> fluidNCellIdx, fluidPredNCellIdx;
	std::unordered_map<unsigned int, std::vector<unsigned int>> wallNCellIdx, wallPredNCellIdx;
	
	std::vector<Particle> brick;

	//particle constant
	unsigned int numFluidParticlesX; // for particle initial position
	unsigned int numFluidParticlesY;
	unsigned int numFluidParticlesZ;
	float sideLenX;					// for particle initial position
	float sideLenY;
	float sideLenZ;



	//vectors

	glm::vec3 forceVis(unsigned int idx) {
		return glm::vec3(0);
	}
	glm::vec3 forceExt(unsigned int idx) {
		return glm::vec3(0);
	}

	glm::vec3 forceP(unsigned int particleIdx);
	//scalers

	float calcDelta(unsigned int particleIdx);
	float calcDensity(unsigned int particleIdx);
	float calcPredDensity(unsigned int particleIdx);


	//update sequence

	void __sortFluidParticles(SORTMODE);

	void __makeNCellMap();
	void __nCellSearch(std::vector<unsigned int>* fluidNIdx, std::vector<unsigned int>* wallNIdx, glm::ivec3);
	void __makePredNCellMap();

	void __initializeFrameStates();

	void __pos_vel_Predictive();
	void __predOutGridResol(unsigned int i);

	unsigned int __density_Predictive();

	void __corrective();



	void __calcDelta();

	void __sceneInitialize();
	void __brickInit();
	void __gridInit();
	void __sorterInit(); // Z_sort 선언. particles이거 비어있는데 이거 build하기.

	void __addWallParticles(glm::vec3);
	void __brickDensityModificationInnerCorner(Particle*);
	void __particleInit(int mode);

	glm::ivec3 __getCellCoord(glm::vec3);

	void __particlesArrInit(ParticlesArray*,unsigned int numbers,ARRTYPE arrtype);
	void __appendParticle(ParticlesArray*, Particle, ARRTYPE arrtype);
	void __deleteParticleArr(ParticlesArray*, ARRTYPE arrtype);



	float W(float dist) {
		float res = 0.0f;
		const float q = dist / coreRad;
		

		const float m_k = 8.0f / (3.141592f * coreRad * coreRad * coreRad);

		if (q <= 1.0f) {
			
			if (q <= 0.5f) {
				const float q2 = q * q;
				const float q3 = q * q * q;
				res = m_k * (6.0f * q3 - 6.0f * q2 + 1.0f);
			}
			else {
				res = m_k * (2.0f * std::pow(1.0f - q, 3));
			}

		}
		return res;
	}

	glm::vec3 gradW(glm::vec3 dist) {


		glm::vec3 res;
		const float rl = glm::length(dist);
		const float  q = rl / coreRad;

		const float m_l = 48.0f / (3.141592f * coreRad * coreRad * coreRad);

		if ((rl > 1.0e-7) && (q <= 1.0))
		{
			glm::vec3 gradq = dist / rl;
			gradq /= coreRad;
			if (q <= 0.5)
			{
				res = m_l * q * (3.0f * q - 2.0f) * gradq;
			}
			else
			{
				const float factor = 1.0f - q;
				res = m_l * (-factor * factor) * gradq;
			}
		}
		else
			res = glm::vec3(0);

		return res;
	}



};




/*
   int a[1000] = { 0 };
   int temp = 0;


   a[idx] += 1;
   temp++;


   bool chkFlag = false;
   for (int tempidx = 0; tempidx < 1000; tempidx++) {
	   if (a[tempidx] != 1) {
		   std::cout << "tempIdxtempIdxtempIdxtempIdxtempIdxtempIdxtempIdxtempIdxtempIdxtempIdxtempIdxtempIdxtempIdxtempIdxtempIdxtempIdxtempIdxtempIdxtempIdxtempIdx : " << tempidx << std::endl;
		   chkFlag = true;
	   }
   }

   std::cout << temp << std::endl;
   system("PAUSE");
   */