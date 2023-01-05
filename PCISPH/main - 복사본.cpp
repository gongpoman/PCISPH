#include<iostream>
#include<vector>
#include<tuple>

#include<map>
#include<cmath>

#include<glad/glad.h>
#include<GLFW/glfw3.h>

#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"


#include"camera.h"
#include"shader.h"
#include"model.h"
#include "myutil.h"
#include "particle.h"


void particleInit(int mode);
void addWallParticles(glm::vec3 cellOrigin,unsigned int);
void gridInit();
void brickInit();


void instanceMat();

void neighborSearch(unsigned int Idx);
void neighborSearchPred(unsigned int Idx);


void __cubeBoundaryCellIdx();
void __cubeBoundaryCellIdxPred();
glm::ivec3 __cubeCellCoord(glm::vec3 pos);
void predoutBoundarySolution(Particle*);

void PCIupdate();

float calcPredDensity(unsigned int particleIdx);
float calcDensity(unsigned int particleIdx);
float calcDelta(unsigned int particleIdx);

glm::vec3 forcePressureSpiky(unsigned int particleIdx);
glm::vec3 forceVis(unsigned int particleIdx);
glm::vec3 forceSurfaceTension(unsigned int particleIdx);




extern const unsigned int SCR_WIDTH = 1280;
extern const unsigned int SCR_HEIGHT = 720;
extern Camera cam(glm::vec3(2.0f, 1.3f, 2.0f));

extern float lastX, lastY;
extern bool isFirstMove = true;

constexpr float linearVisc = 0.25f;
constexpr float quadVisc = .5f;

extern constexpr float deltaTime = 1 / 180.0f;
glm::mat4* instWorlds;


// particle setting.
constexpr unsigned int numWaterParticleX = 15;
constexpr unsigned int numWaterParticleY = 15;
constexpr unsigned int numWaterParticleZ = 15;

constexpr unsigned int numWaterParticle = numWaterParticleX * numWaterParticleY * numWaterParticleZ;
unsigned int numWallParticle=0;
unsigned int numParticle=0;
unsigned int numDrawParticle;


constexpr float restDensity = 997.0f;
constexpr float radius = 0.025f;
constexpr float particleMass = restDensity * 4.0f / 3.0f * 3.141592f * radius * radius * radius;
constexpr float coreRad = radius * 4.0f;

const float sideLenX =  numWaterParticleX / numWaterParticle * std::pow(numWaterParticle * particleMass / restDensity, 1.0f / 3.0f);
const float sideLenY = numWaterParticleY / numWaterParticle * std::pow(numWaterParticle * particleMass / restDensity, 1.0f / 3.0f);
const float sideLenZ = numWaterParticleZ / numWaterParticle * std::pow(numWaterParticle * particleMass / restDensity, 1.0f / 3.0f);

constexpr glm::vec3 force_g = glm::vec3(0.0f, -9.8f, 0.0f) * particleMass;

// boundary setting
constexpr float boundaryX = 3.0f;
constexpr float boundaryY = 1.5f;
constexpr float boundaryZ = 3.0f;

constexpr unsigned int nGridDivX = (unsigned int)(boundaryX / coreRad);
constexpr unsigned int nGridDivY = (unsigned int)(boundaryY / coreRad);
constexpr unsigned int nGridDivZ = (unsigned int)(boundaryZ / coreRad);

unsigned int upperMaxGridDiv;

std::vector<unsigned int>* neighborIdices;
vector<Particle> brick; // cell boundary �� ĭ.



// predict - correct setting
constexpr float eta = 0.01f; 
constexpr unsigned int minIter = 0;
constexpr unsigned int maxIter = 100;






Particle* particles = NULL;
vector<unsigned int> mortonCount;
std::map<unsigned int, unsigned int> idxMap;
std::vector<unsigned int> particleIndices;

Particle* predParticles;
vector<unsigned int> mortonCountPred;
std::map<unsigned int, unsigned int> predIdxMap;
std::vector<unsigned int> predParticleIndices;




constexpr bool drawWallParticle = true;
GLFWwindow* glInitialize();
void initialStateLog();
void sceanInitialize();
std::tuple<Shader*,Model*,unsigned int> renderInitialize();



int main(){
    
    static unsigned int frame = 0;

    //glInit
    GLFWwindow* window = glInitialize();
    if (!window)
        return -1;


    // particle, grid init
    sceanInitialize();



    //Sphere render setting
    Shader* particleShader;
    Model* sphere;
    unsigned int instVBO;
    std::tie(particleShader, sphere,instVBO) = renderInitialize();

    // print scean Log
    initialStateLog();

    while (!glfwWindowShouldClose(window))
    {
        frame++;
        std::cout<<frame <<" frame start=====================================================================================================================" << std::endl;

        processInput(window);
        glClearColor(0.9f, 0.9f, 0.9f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        ////////////Update
        PCIupdate();
        instanceMat();
        __mapBuffer(instVBO, instWorlds, numDrawParticle* sizeof(glm::mat4));


        ////////////render 
        particleShader->use();
        for (unsigned int i = 0; i < sphere->meshes.size(); i++) {
            glBindVertexArray(sphere->meshes[i].VAO);
            glDrawElementsInstanced(GL_TRIANGLES, sphere->meshes[i].indices.size(), GL_UNSIGNED_INT, 0, numDrawParticle);
            glBindVertexArray(0);
        }
        glUseProgram(0);

        glfwSwapBuffers(window);
        glfwPollEvents();

        std::cout << "end=====================================================================================================================" << std::endl;
        std::cout << " =====================================================================================================================" << std::endl;
        std::cout << " =====================================================================================================================" << std::endl;
        std::cout << " =====================================================================================================================" << std::endl;
        std::cout << " =====================================================================================================================" << std::endl;

        /* To Test
        for (unsigned int i = 0; i < numParticle; i++) {
            if (particles[i].pos.y > 3.0f)
                std::cout << i<<" : "<<particles[i].pos.y << std::endl;
        }
        */
    }

    std::cout << frame << "frames rendered" << std::endl;

    //global variables 
    delete[] predParticles;
    delete[] particles; 
    delete[] instWorlds;
    // local variables
    delete sphere;
    delete particleShader;

    glfwTerminate();
    return 0;
}








float calcDelta(unsigned int particleIdx) {

    const float kernelConst = 15.0f / 3.141592f / (coreRad * coreRad * coreRad * coreRad * coreRad * coreRad);

    float beta = deltaTime* deltaTime * particleMass*particleMass * 2.0f / restDensity/ restDensity;

    glm::vec3 sumOfGrad = glm::vec3(0);

    Particle target = particles[particleIdx];

    float sum = 0.0f;

    unsigned int temptemptemp = 0;

    for (auto iter = neighbors->begin(); iter != neighbors->end(); iter++) {
        glm::vec3 dist = (*iter).pos - target.pos;
        float distLen = glm::length(dist);

        if (distLen - coreRad < FLT_EPSILON && distLen>FLT_EPSILON) {
            glm::vec3 grad_Wp = glm::normalize(dist);
            grad_Wp *= (coreRad - distLen) * (coreRad - distLen) * 3.0f * kernelConst;

            sumOfGrad += grad_Wp;
            sum += glm::dot(grad_Wp, grad_Wp);
            temptemptemp++;
        }
    }
    sum += glm::dot(sumOfGrad, sumOfGrad);
    beta *= sum;
    beta = 1.0f / beta;
    
    if (sum <= FLT_EPSILON) {

        std::cout<< particleIdx << "th particle : " << "DELTA EXPLODE" << std::endl;
        
        std::cout << " ��� ���Ҵµ�?" << temptemptemp << std::endl;

        std::cout << " neighbor ���� : " << neighbors->size() <<"density Var : " << particles[particleIdx].densityVar << std::endl;

        if (neighbors->size() == 1)
            std::cout<<"�̿��� idx" << particles[particleIdx].actualIdx << std::endl;

        std::cout << particles[particleIdx].pos.x << ": " << particles[particleIdx].pos.y << " : " << particles[particleIdx].pos.z << std::endl;

        system("PAUSE");
    }
    // �̰� ���� ������ predicted �� �༮�� �������� ���� �ִٰ�. �ٵ�. ���� dist ���� position���θ� ���ݾ�. px�� �ƴ϶�. �׷��ϱ� ������.
    // ���� actual cellIdx�� px�� ���� �� cellIdx�� �ٸ��ž�.

    return beta;
}

float calcPredDensity(unsigned int particleIdx) {

    const glm::vec3 pos = predParticles[particleIdx].px;
    const float coreRad2 = coreRad * coreRad;

    float sum = 0;


    for (auto iter = neighbors->begin(); iter != neighbors->end(); iter++) {

        float len = glm::length(pos - (*iter).px);
        if (len - coreRad < FLT_EPSILON) {
            len = std::pow((coreRad2 - len * len), 3.0f);
            sum += len;
        }
    }

    sum *= 1.566682f / std::pow(coreRad, 9);// 315.0f / 64.0f / 3.141592f = 1.5666817...;
    sum *= particleMass;

    return sum;
}

float calcDensity(unsigned int particleIdx) {

    const glm::vec3 pos = particles[particleIdx].pos;
    const float coreRad2 = coreRad * coreRad;
    float sum = 0;

    for (auto iter = neighbors->begin(); iter != neighbors->end(); iter++) {

        float len = glm::length(pos - (*iter).pos);
        if (len - coreRad < FLT_EPSILON) {
            len = std::pow((coreRad2 - len * len), 3.0f);
            sum += len;
        }
    }

    sum *= 1.566682f / std::pow(coreRad, 9);// 315.0f / 64.0f / 3.141592f = 1.5666817...;
    sum *= particleMass;

    return sum;
}

glm::vec3 forceVis(unsigned int particleIdx) {

    const float kernelConst = 2.0f / (coreRad * coreRad * coreRad * 1.772454f); // root pi = 1.7724538...

    Particle target = particles[particleIdx];
    glm::vec3 netForceVis = glm::vec3(0);

    for (auto iter = neighbors->begin(); iter != neighbors->end(); iter++) {
        glm::vec3 dist = ((*iter).pos - target.pos);
        float distLen = glm::length(dist);

        float distRatio = distLen / coreRad;

        float mu = glm::dot((*iter).vel - target.vel, dist);

        if (mu < FLT_EPSILON && distLen - coreRad < FLT_EPSILON && distLen > FLT_EPSILON) {
            float pipi = 0;
            mu /= (0.01f * coreRad + distLen * distRatio);
            pipi = -linearVisc * mu + quadVisc * mu * mu;
            //pipi /= restDensity;
            pipi /= (target.density + (*iter).density) / 2.0f;
            glm::vec3 grad_Wvis = pipi * dist * kernelConst / std::pow(2.718282f, distRatio);
            netForceVis -= grad_Wvis;
        }
    }
    return netForceVis;
}

glm::vec3 forceSurfaceTension(unsigned int particleIdx) {
    return glm::vec3(0);
}

// ALTER
glm::vec3 forcePressureSpiky(unsigned int particleIdx) {
    const float kernelConst = 15.0f / 3.141592f / (coreRad * coreRad * coreRad * coreRad * coreRad * coreRad);
    Particle target = particles[particleIdx];
    glm::vec3 netForceP = glm::vec3(0);

    for (auto iter = neighbors->begin(); iter != neighbors->end(); iter++) {
        glm::vec3 dist = (*iter).pos - target.pos;
        float distLen = glm::length(dist);

        if (distLen - coreRad < FLT_EPSILON && distLen>FLT_EPSILON) {
            glm::vec3 grad_Wp = glm::normalize(dist);


            grad_Wp *= (coreRad - distLen) * (coreRad - distLen) * 3.0f * kernelConst;

            grad_Wp *= particleMass * particleMass * ((*iter).pressure / (*iter).density / (*iter).density + particles[particleIdx].pressure / particles[particleIdx].density / particles[particleIdx].density);
            //grad_Wp *= particleMass * particleMass * ((*iter).pressure + particles[particleIdx].pressure) / restDensity / restDensity;

            if (glm::length(grad_Wp) > 10000000000000000000.0f) {
                std::cout << "3. " << (*iter).pressure << " : " << particles[particleIdx].pressure << std::endl;
                system("PAUSE");
            }

            netForceP -= grad_Wp;
        }
    }

    return netForceP;
}



glm::vec3 gridLocalPos(glm::vec3 pos) {
    float x = (pos.x / coreRad < FLT_EPSILON) ? ((int)(pos.x / coreRad) - 1) * coreRad : (int)(pos.x / coreRad) * coreRad;
    float y = (pos.y / coreRad < FLT_EPSILON) ? ((int)(pos.y / coreRad) - 1) * coreRad : (int)(pos.y / coreRad) * coreRad;
    float z = (pos.z / coreRad < FLT_EPSILON) ? ((int)(pos.z / coreRad) - 1) * coreRad : (int)(pos.z / coreRad) * coreRad;

    return pos - glm::vec3(x, y, z);
}



// TODO water particle�� ���ؼ��� �ؾߵǴ� ������, wall particle ���� ����ؾ� �ϴ� ������ cut�� �� �ؾߵ�.
void PCIupdate() {


    ////////////////////////////////////////////////////////////
    /// initial state calculation
    ////////////////////////////////////////////////////////////       

    __cubeBoundaryCellIdx();
    unsigned int currentCellIdx = -1;
    unsigned int prevCellIdx = -1;

    for (auto particleIdx = particleIndices.begin(); particleIdx != particleIndices.end(); particleIdx++) {

        currentCellIdx = particles[*particleIdx].cellIdx;
        if (prevCellIdx != currentCellIdx) {
            neighbors->clear();
            neighborSearch(*particleIdx);
        }
        particles[*particleIdx].pressure = 0.0f;
        particles[*particleIdx].force_p = glm::vec3(0);
        particles[*particleIdx].force_vis = glm::vec3(0);
        particles[*particleIdx].force_ext = glm::vec3(0);
        prevCellIdx = currentCellIdx;
    }
    ////////////////////////////////////////////////////////////       



    bool densityOverFlag = true;
    unsigned int iteration = 0;

    while ((iteration < minIter || densityOverFlag) && iteration < maxIter) {

        int largerParticleNum = 0;
        densityOverFlag = false;


        ////////////////////////////////////////////////////////////
        /// pos, vel prediction. 
        ////////////////////////////////////////////////////////////       

        for (auto particleIdx = particleIndices.begin(); particleIdx != particleIndices.end(); particleIdx++) {
            particles[*particleIdx].pv = particles[*particleIdx].vel + (( particles[*particleIdx].force_ext + particles[*particleIdx].force_p + particles[*particleIdx].force_vis) / particleMass + acc_g) * deltaTime;
            particles[*particleIdx].px = particles[*particleIdx].pos + particles[*particleIdx].pv * deltaTime;
            predoutBoundarySolution(&particles[*particleIdx]);
        }
        ////////////////////////////////////////////////////////////       


        __cubeBoundaryCellIdxPred();

        ////////////////////////////////////////////////////////////
        /// Pred density variation.
        ////////////////////////////////////////////////////////////


        prevCellIdx = -1;
        for (auto predParticleIdx = predParticleIndices.begin(); predParticleIdx != predParticleIndices.end(); predParticleIdx++) {

            currentCellIdx = predParticles[*predParticleIdx].cellIdxPred;
            if (prevCellIdx != currentCellIdx) {
                neighbors->clear();
                neighborSearchPred(*predParticleIdx);
            }

            predParticles[*predParticleIdx].density = calcPredDensity(*predParticleIdx);

            predParticles[*predParticleIdx].densityVar = (predParticles[*predParticleIdx].density - restDensity >= 0) ? predParticles[*predParticleIdx].density - restDensity : 0.0f;

            if ((predParticles[*predParticleIdx].densityVar / restDensity) > eta) {
                largerParticleNum++;
                densityOverFlag = true;


            }

            prevCellIdx = currentCellIdx;
        }

        for (auto predParticleIdx = predParticleIndices.begin(); predParticleIdx != predParticleIndices.end(); predParticleIdx++) {
            particles[predParticles[*predParticleIdx].actualIdx].densityVar = predParticles[*predParticleIdx].densityVar;
        }
        ////////////////////////////////////////////////////////////       



        ////////////////////////////////////////////////////////////
        /// actual density
        ////////////////////////////////////////////////////////////       

        prevCellIdx = -1;
        for (auto particleIdx = particleIndices.begin(); particleIdx != particleIndices.end(); particleIdx++) {

            currentCellIdx = particles[*particleIdx].cellIdx;
            if (prevCellIdx != currentCellIdx) {
                neighbors->clear();
                neighborSearch(*particleIdx);
            }
            particles[*particleIdx].density = calcDensity(*particleIdx);


            if ((particles[*particleIdx].densityVar / restDensity) > eta)                                       
                particles[*particleIdx].pressure += particles[*particleIdx].densityVar * calcDelta(*particleIdx);


            prevCellIdx = currentCellIdx;
        }
        ////////////////////////////////////////////////////////////       



        ////////////////////////////////////////////////////////////
        /// Pressure Force
        ////////////////////////////////////////////////////////////        �ڸ��ϰ� waterParticle�� �ٲ�� ��������.

        prevCellIdx = -1;
        for (auto particleIdx = particleIndices.begin(); particleIdx != particleIndices.end(); particleIdx++) {

            currentCellIdx = particles[*particleIdx].cellIdx;
            if (prevCellIdx != currentCellIdx) {
                neighbors->clear();
                neighborSearch(*particleIdx);
            }
            particles[*particleIdx].force_p = forcePressureSpiky(*particleIdx);

            prevCellIdx = currentCellIdx;
        }
        ////////////////////////////////////////////////////////////       

        std::cout << "============================================================larger part : " << largerParticleNum << std::endl;
        iteration++;
    }

    for (auto particleIdx = particleIndices.begin(); particleIdx != particleIndices.end(); particleIdx++) {
        particles[*particleIdx].vel = particles[*particleIdx].pv;
        particles[*particleIdx].pos = particles[*particleIdx].px;
    }

    std::cout << iteration << std::endl;

}


//=============================================================== �Ƹ� �� ���� ���� �ʾ�����...


void instanceMat() {
    int idx = 0;
    for (unsigned int i = 0; i < numParticle; i++) {
        if (!particles[i].isWall || drawWallParticle) {
            instWorlds[idx] = glm::translate(glm::mat4(1.0f), particles[i].pos);
            instWorlds[idx] = glm::scale(instWorlds[idx], glm::vec3(0.025f));
            idx++;
        }
    }
}


void predoutBoundarySolution(Particle* particle) {
    glm::ivec3 cellIdx = __cubeCellCoord(particle->px);
    glm::ivec3 flag;

    flag.x = (cellIdx.x == -1) ? -1 : (cellIdx.x >= nGridDivisoin) ? 1 : 0;
    flag.z = (cellIdx.z == -1) ? -1 : (cellIdx.z >= nGridDivisoin) ? 1 : 0;
    flag.y = (cellIdx.y == -1) ? -1 : (cellIdx.y >= nGridDivisoin) ? 1 : 0;


    if (flag.x == -1) {
        particle->pv.x *= -1;
        particle->px.x += 2 * (-cubicBoundaryLen / 2.0f - particle->px.x);
    }
    if (flag.x == 1) {
        particle->pv.x *= -1;
        particle->px.x += 2 * (cubicBoundaryLen / 2.0f - particle->px.x);
    }

    if (flag.z == -1) {
        particle->pv.z *= -1;
        particle->px.z += 2 * (-cubicBoundaryLen / 2.0f - particle->px.z);
    }
    if (flag.z == 1) {
        particle->pv.z *= -1;
        particle->px.z += 2 * (cubicBoundaryLen / 2.0f - particle->px.z);
    }

    if (flag.y == -1) {
        particle->pv.y *= -0.999f;
        particle->px.y = -(particle->px.y);
    }
    if (flag.y == 1) {
        particle->pv.y *= -1;
        particle->px.y += 2 * (+cubicBoundaryLen - particle->px.y);
    }
}


void neighborSearch(unsigned int idx) {

    for (unsigned int i = 0; i < neighborIdices[particles[idx].cellIdx].size(); i++) {
        unsigned int cidx = neighborIdices[particles[idx].cellIdx][i];
        if (idxMap.count(cidx) != 0) {
            cidx = idxMap[cidx];
            unsigned int startPoint = (cidx == 0) ? 0 : mortonCount[cidx - 1];
            unsigned int endPoint = mortonCount[cidx];

            for (unsigned int j = startPoint; j < endPoint; j++) {
                neighbors->push_back(particles[j]);
            }
        }
    }
}

void neighborSearchPred(unsigned int idx) {

    for (unsigned int i = 0; i < neighborIdices[predParticles[idx].cellIdxPred].size(); i++) {
        unsigned int cidx = neighborIdices[predParticles[idx].cellIdxPred][i];
        if (predIdxMap.count(cidx) != 0) {
            cidx = predIdxMap[cidx];
            unsigned int startPoint = (cidx == 0) ? 0 : mortonCountPred[cidx - 1];
            unsigned int endPoint = mortonCountPred[cidx];

            for (unsigned int j = startPoint; j < endPoint; j++) {
                neighbors->push_back(predParticles[j]);
            }
        }
    }
}

void __cubeBoundaryCellIdx() {          // z indexing.�� �µ��� ������ �ϰ�, ���ĵ� particles array���� ��� idx���� �������� � cellIdx�� ������ particle���� �� ������ ���� �ڵ�.

    particleIndices.clear();
    mortonCount.clear(); // ��� idx���� d��� idx����.
    idxMap.clear();

    unsigned int* mortonCounter = new unsigned[upperNGridDiv * upperNGridDiv * upperNGridDiv]; // temp
    unsigned int outBoundardyCount = 0;

    for (unsigned int i = 0; i < upperNGridDiv * upperNGridDiv * upperNGridDiv; i++)
        mortonCounter[i] = 0;


    for (unsigned int i = 0; i < numParticle; i++) {
        if (!particles[i].isWall) {
            glm::ivec3 cellCoord = __cubeCellCoord(particles[i].pos);

            particles[i].cellIdx = __cubeMorton(cellCoord.x, cellCoord.y, cellCoord.z);

            if (particles[i].cellIdx == (unsigned int)(-1)) {
                outBoundardyCount++;
                std::cout << "ACTUAL CASE : " << i << "th particle" << std::endl;
                std::cout << particles[i].pos.x << " : " << particles[i].pos.y << " : " << particles[i].pos.z << std::endl;
                std::cout << particles[i].vel.x << " : " << particles[i].vel.y << " : " << particles[i].vel.z << std::endl;
                
            }
            else
                mortonCounter[particles[i].cellIdx]++;
        }
        else {
            mortonCounter[particles[i].cellIdx]++;
        }
    }

    if (outBoundardyCount > 0) {
        std::cout << "out boundary particle : " << outBoundardyCount << std::endl;
        system("PAUSE");
    }

    unsigned int invIdx = 0;


    for (unsigned int i = 0; i < upperNGridDiv * upperNGridDiv * upperNGridDiv; i++)
        if (mortonCounter[i] > 0) {
            mortonCount.push_back(mortonCounter[i]);
            idxMap.insert(std::pair<unsigned int, unsigned int>(i, invIdx++));              // morton code cell index�� �־��ָ�. �� cell idx�� morton counter�� ���° index������ �����ϴ� map  
        }

    for (auto iter = mortonCount.begin() + 1; iter != mortonCount.end(); iter++)
        *(iter) += *(iter - 1);

    std::vector<std::vector<Particle>> tempParticles;
    for (unsigned int i = 0; i < mortonCount.size(); i++)
        tempParticles.push_back(std::vector<Particle>());

    for (unsigned int i = 0; i < numParticle; i++) {
        unsigned int idx = idxMap[particles[i].cellIdx];
        tempParticles[idx].push_back(particles[i]);
    }

    unsigned int i = 0;
    for (auto iter = tempParticles.begin(); iter != tempParticles.end(); iter++) {
        for (auto iteriter = iter->begin(); iteriter != iter->end(); iteriter++) {
            particles[i] = *iteriter;

            if (!particles[i].isWall) {
                particles[i].actualIdx = i;
                particleIndices.push_back(i);
            }

            i++;
        }
    }

    delete[] mortonCounter;
}
void __cubeBoundaryCellIdxPred() {          // z indexing.�� �µ��� ������ �ϰ�, ���ĵ� particles array���� ��� idx���� �������� � cellIdx�� ������ particle���� �� ������ ���� �ڵ�.



    predParticleIndices.clear();
    mortonCountPred.clear(); // ��� idx���� d��� idx����.
    predIdxMap.clear();

    unsigned int* mortonCounter = new unsigned[upperNGridDiv * upperNGridDiv * upperNGridDiv]; // temp
    unsigned int outBoundardyCount = 0;

    for (unsigned int i = 0; i < upperNGridDiv * upperNGridDiv * upperNGridDiv; i++)
        mortonCounter[i] = 0;


    for (unsigned int i = 0; i < numParticle; i++) {
        if (!particles[i].isWall) {
            glm::ivec3 cellCoord = __cubeCellCoord(particles[i].px);

            particles[i].cellIdxPred = __cubeMorton(cellCoord.x, cellCoord.y, cellCoord.z);

            if (particles[i].cellIdxPred == (unsigned int)(-1)) {
                outBoundardyCount++;
                std::cout << "Pred CASE : " << i << "th particle" << std::endl;
                std::cout << particles[i].px.x << " : " << particles[i].px.y << " : " << particles[i].px.z << std::endl;
                std::cout << particles[i].pv.x << " : " << particles[i].pv.y << " : " << particles[i].pv.z << std::endl;

            }
            else
                mortonCounter[particles[i].cellIdxPred]++;
        }
        else {
            mortonCounter[particles[i].cellIdxPred]++;
        }
    }

    if (outBoundardyCount > 0) {
        std::cout << "out boundary particle : " << outBoundardyCount << std::endl;
        system("PAUSE");
    }

    unsigned int invIdx = 0;


    for (unsigned int i = 0; i < upperNGridDiv * upperNGridDiv * upperNGridDiv; i++)
        if (mortonCounter[i] > 0) {
            mortonCountPred.push_back(mortonCounter[i]);
            predIdxMap.insert(std::pair<unsigned int, unsigned int>(i, invIdx++));              // morton code cell index�� �־��ָ�. �� cell idx�� morton counter�� ���° index������ �����ϴ� map  
        }

    for (auto iter = mortonCountPred.begin() + 1; iter != mortonCountPred.end(); iter++)
        *(iter) += *(iter - 1);

    std::vector<std::vector<Particle>> tempParticles;
    for (unsigned int i = 0; i < mortonCountPred.size(); i++)
        tempParticles.push_back(std::vector<Particle>());

    for (unsigned int i = 0; i < numParticle; i++) {
        unsigned int idx = predIdxMap[particles[i].cellIdxPred];
        tempParticles[idx].push_back(particles[i]);
    }

    unsigned int i = 0;

    for (auto iter = tempParticles.begin(); iter != tempParticles.end(); iter++) {
        for (auto iteriter = iter->begin(); iteriter != iter->end(); iteriter++) {
            predParticles[i] = *iteriter;

            if (!predParticles[i].isWall) {
                predParticleIndices.push_back(i);
            }
            i++;
        }
    }

    delete[] mortonCounter;
}


glm::ivec3 __cubeCellCoord(glm::vec3 pos) {

    float x = ((pos.x + cubicBoundaryLen / 2.0f) / cubicBoundaryLen * nGridDivisoin);
    float y = (pos.y / cubicBoundaryLen * nGridDivisoin);
    float z = ((pos.z + cubicBoundaryLen / 2.0f) / cubicBoundaryLen * nGridDivisoin);

    x = (x + FLT_EPSILON <= 0) ? -1.0f : x;
    z = (z + FLT_EPSILON <= 0) ? -1.0f : z;
    y = (pos.y + FLT_EPSILON <= 0) ? -1.0f : y;

    return glm::ivec3((int)x, (int)y, (int)z);
}





GLFWwindow* glInitialize() {
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "PCISPH", NULL, NULL);

    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return NULL;
    }
    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return NULL;
    }

    glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
    glEnable(GL_DEPTH_TEST);

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);


    //    glfwSetCursorPosCallback(window, mouse_callback);
    //    glfwSetScrollCallback(window, scroll_callback);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    return window;
}
void initialStateLog() {
    std::cout << "========================================" << std::endl;

    std::cout << "totalParticles : " << numParticle << std::endl <<
        "Water Particles : " << numWaterParticle << std::endl <<
        "Wall Particles : " << numWallParticle << std::endl;

    std::cout << "========================================" << std::endl;
    system("PAUSE");
}
std::tuple < Shader*, Model*,unsigned int > renderInitialize() {

    Shader* particleShader = new Shader("resources/shader/paricleL_vs.txt", "resources/shader/paricleL_fs.txt");
    Model* sphere = new Model("resources/objects/sphere.obj");

    particleShader->use();
    glm::mat4 viewMat = glm::lookAt(cam.position, glm::vec3(0.0f, 0.5f, 0.0f), cam.up);
    glm::mat4 projMat = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
    glUniformMatrix4fv(glGetUniformLocation(particleShader->ID, "view"), 1, GL_FALSE, &viewMat[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(particleShader->ID, "proj"), 1, GL_FALSE, &projMat[0][0]);
    glUniform3f(glGetUniformLocation(particleShader->ID, "lightPos"), cam.position.x, cam.position.y, cam.position.z);
    glUseProgram(0);

    unsigned int instVBO;

    for (unsigned int i = 0; i < sphere->meshes.size(); i++) {  // Definitely, sphere->meshes.size() = 1

        glBindVertexArray(sphere->meshes[i].VAO);

        glGenBuffers(1, &instVBO);
        glBindBuffer(GL_ARRAY_BUFFER, instVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(glm::mat4) * numDrawParticle, &instWorlds[0][0], GL_DYNAMIC_DRAW);

        glEnableVertexAttribArray(3);
        glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), (void*)0);
        glEnableVertexAttribArray(4);
        glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), (void*)(1 * sizeof(glm::vec4)));
        glEnableVertexAttribArray(5);
        glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), (void*)(2 * sizeof(glm::vec4)));
        glEnableVertexAttribArray(6);
        glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), (void*)(3 * sizeof(glm::vec4)));

        glVertexAttribDivisor(3, 1);
        glVertexAttribDivisor(4, 1);
        glVertexAttribDivisor(5, 1);
        glVertexAttribDivisor(6, 1);

        glBindVertexArray(0);
    }

    return { particleShader,sphere,instVBO};
}
void sceanInitialize() {

    brickInit();
    upperNGridDiv = __nthDigit(nGridDivisoin);
    neighborIdices = new std::vector<unsigned int>[upperNGridDiv * upperNGridDiv * upperNGridDiv]; // nGridDiv^3 ��ŭ�� �ƴ϶�. nGridDiv�� 2���������� �ø���ŭ ��ƾ� ��. ngriddiv= 27 -> 32^3�� ��ƾ� ��.
    gridInit();

    particleInit(0);

    numDrawParticle = (drawWallParticle) ? numParticle : numWaterParticle;

    instWorlds = new glm::mat4[numDrawParticle];
    instanceMat();

}


void brickInit() {
    const float wallDensity = restDensity; // wall�� ��ü density�� ���� density�� ������ ����.

    Particle wallParticle = Particle();
    wallParticle.isWall = true;
    wallParticle.density = wallDensity;
    wallParticle.pos = glm::vec3(0);

    wallParticle.pressure = 0.0f; // TODO calc pressure ���� �ھƳ��� �� ã�ƾ� ��. ���� �ȵǸ� ��.. estimation ��� ����� ��.

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
// TODO ������ü�� boundary�ƴϸ� ���� ���ľ� ��.
void gridInit() {

    ////////////////////////////////////////////////
    // �������� ����. morton code�� �� ��. ������� 3���� ���� cell�� ������ 27���� �ִٰ� ����. 27*27*27. �׷��� ��� ���� �� �ִ� index�� ����Ǵ� 27^3 �� ����. �ٵ� morton code�� encoding���� �� ������ index���� 0~32*32*32�� ������ ������.
    // �׷��� cell idx�� ���� �� ���� ����ؾ� �� �͵��� ����.
    // �̰͋����� ���� r = 0.025
    // coreRad = 5*r ���� �� �� ������ 32���� cell�� ���ͼ� radius, coreRad �̷��� ������.
    ////////////////////////////////////////////////
    unsigned int numBoundary = 0;

    std::vector<glm::vec3> boundaryOrigins;

    for (unsigned int i = 0; i < nGridDivX; i++) {
        for (unsigned int j = 0; j < nGridDivY; j++) {
            for (unsigned int k = 0; k < nGridDivZ; k++) {

                unsigned int cidx = __cubeMorton(i, j, k);

                bool isBoundaryCell = false;

                //////////////////////////////////////////////// 
                /// mortonIdx neighbor
                //////////////////////////////////////////////// 
                for (int ii = -1; ii < 2; ii++) {
                    for (int jj = -1; jj < 2; jj++) {
                        for (int kk = -1; kk < 2; kk++) {

                            bool x1 = i == 0 && ii == -1;
                            bool x2 = i == nGridDivX - 1 && ii == 1;
                            bool y1 = j == 0 && jj == -1;
                            bool y2 = j == nGridDivY - 1 && jj == 1;
                            bool z1 = k == 0 && kk == -1;
                            bool z2 = k == nGridDivZ - 1 && kk == 1;

                            if (!(x1 || x2 || y1 || y2 || z1 || z2)) {
                                glm::ivec3 coord = glm::ivec3(i + ii, j + jj, k + kk);
                                unsigned int neighborCell = __cubeMorton(coord.x, coord.y, coord.z);

                                if (neighborCell != (unsigned int)(-1))     // �ٸ� ������� ���� ��. ������ 32*32�� ������ �� �������ϱ�. �̷��� ���� �ص� �Ǵµ�. coreRad�ٲٰų� �׷��� �ٷ� �� ������ �������� ��.
                                    neighborIdices[cidx].push_back(neighborCell);  // �ٵ� �ϴ� ������ ���� 4*4 ���̷� �ϰ� �Ǹ� boundary cell�� particle���� �ʹ� ���Ƽ� ������ �ȵ�µ�  grid�� ���̷��� ��Ե� �ؾߵ�.
                            }
                            else {
                                isBoundaryCell = true;
                            }
                        }
                    }
                }
                ////////////////////////////////////////////////


                if (isBoundaryCell && (j + 1) * coreRad <= 1.0f) {
                    numBoundary++;
                    boundaryOrigins.push_back(glm::vec3(i * coreRad - boundaryX / 2.0f, j * coreRad, k * coreRad - boundaryZ / 2.0f));
                }
            }
        }
    }

    numWallParticle = brick.size() * numBoundary; // determined number of boundary particle 

    numParticle = numWallParticle + numWaterParticle;

    particles = new Particle[numParticle];
    predParticles = new Particle[numParticle];

    //boundary origins ���鼭. numParticle ���� numWallParticle-1���� ä���ֱ�.
    unsigned int offsetIdx = 0;

    for (auto iter = boundaryOrigins.begin(); iter != boundaryOrigins.end(); iter++) {
        addWallParticles(*iter, offsetIdx++);
    }

}
void addWallParticles(glm::vec3 cellOrigin, unsigned int offsetIdx) {

    vector<Particle> tempBrick = brick;

    for (unsigned int i = 0; i < brick.size(); i++) {
        tempBrick[i].pos += cellOrigin;
        tempBrick[i].px = tempBrick[i].pos;

        glm::ivec3 cellIdxCoord = __cubeCellCoord(tempBrick[i].pos);
        tempBrick[i].cellIdx = __cubeMorton(cellIdxCoord.x, cellIdxCoord.y, cellIdxCoord.z);
        tempBrick[i].cellIdxPred = tempBrick[i].cellIdx;

        particles[numWaterParticle + i + offsetIdx * brick.size()] = tempBrick[i];
    }

}

void particleInit(int mode) {

    if (mode == 0) {
        int numRow = pow(numWaterParticle, 1 / 3.0f);


        for (unsigned int i = 0; i < numRow; i++) {
            for (unsigned int j = 0; j < numRow; j++) {
                for (unsigned int k = 0; k < numRow; k++) {

                    particles[i * numRow * numRow + numRow * j + k].pos = glm::vec3(sideLen / 2.0f - sideLen * (i / (float)numRow)
                        , sideLen - sideLen * (j / (float)numRow) + coreRad
                        , sideLen / 2.0f - sideLen * (k / (float)numRow));
                    particles[i * numRow * numRow + numRow * j + k].density = restDensity;

                }
            }
        }
    }
    if (mode == 1) {
        int numRow = pow(numWaterParticle, 1 / 3.0f);

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
