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


void instanceMat();

extern const unsigned int SCR_WIDTH = 1280;
extern const unsigned int SCR_HEIGHT = 720;
extern Camera cam(glm::vec3(0.0f, 1.5f, 3.0f));

PCISPH pcisph(glm::ivec3(20,20,20),glm::vec3(3.0f,1.5f,3.0f), false);

extern float lastX, lastY;
extern bool isFirstMove = true;

extern const float deltaTime = 1 / 180.0f;s

glm::mat4* instWorlds;


GLFWwindow* glInitialize();
void initialStateLog();

std::tuple<Shader*,Model*,unsigned int> renderInitialize();

int main(){
    
    static unsigned int frame = 0;

    //glInit
    GLFWwindow* window = glInitialize();
    if (!window)
        return -1;

    instWorlds = new glm::mat4[pcisph.numDrawParticles];
    instanceMat();

    //particle render setting
    Shader* particleShader;
    Model* sphere;
    unsigned int instVBO;
    std::tie(particleShader, sphere,instVBO) = renderInitialize();

    // print scean Log
    initialStateLog();

    while (!glfwWindowShouldClose(window))
    {
        frame++;
        std::cout<<frame <<" frame start=======================================================================================" << std::endl;

        processInput(window);
        glClearColor(0.9f, 0.9f, 0.9f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        ////////////Update
        pcisph.update();
        instanceMat();
        __mapBuffer(instVBO, instWorlds, pcisph.numDrawParticles* sizeof(glm::mat4));


        ////////////render 
        particleShader->use();
        for (unsigned int i = 0; i < sphere->meshes.size(); i++) {
            glBindVertexArray(sphere->meshes[i].VAO);
            glDrawElementsInstanced(GL_TRIANGLES, sphere->meshes[i].indices.size(), GL_UNSIGNED_INT, 0, pcisph.numDrawParticles);
            glBindVertexArray(0);
        }
        glUseProgram(0);

        glfwSwapBuffers(window);
        glfwPollEvents();

        std::cout << frame << " frame end=======================================================================================" << std::endl;
        std::cout << " =====================================================================================================================" << std::endl;
        std::cout << " =====================================================================================================================" << std::endl;
 
    }

    std::cout << frame << "frames rendered" << std::endl;

    //global variables 
    delete[] instWorlds;
    // local variables
    delete sphere;
    delete particleShader;

    glfwTerminate();
    return 0;
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

    std::cout << "totalParticles : " << pcisph.numParticles << std::endl <<
        "Water Particles : " << pcisph.numFluidParticles << std::endl <<
        "Wall Particles : " << pcisph.numWallParticles << std::endl << std::endl << std::endl;

    std::cout << "DRAWALL ? T/F : " << ((pcisph.drawWall) ? "T" : "F") << std::endl;

    std::cout << "========================================" << std::endl;
    system("PAUSE");
}
std::tuple < Shader*, Model*, unsigned int > renderInitialize() {

    Shader* particleShader = new Shader("resources/shader/paricleL_vs.txt", "resources/shader/paricleL_fs.txt");
    Model* sphere = new Model("resources/objects/sphere.obj");

    cam.position += glm::vec3(pcisph.boundaryX / 2.0f, 0.0f, pcisph.boundaryZ / 2.0f);

    particleShader->use();
    glm::mat4 viewMat = glm::lookAt(cam.position, glm::vec3(0.0f, 0.5f, 0.0f) + glm::vec3(pcisph.boundaryX / 2.0f, 0.0f, pcisph.boundaryZ / 2.0f), cam.up);
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
        glBufferData(GL_ARRAY_BUFFER, sizeof(glm::mat4) * pcisph.numDrawParticles, &instWorlds[0][0], GL_DYNAMIC_DRAW);

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

    return { particleShader,sphere,instVBO };
}

void instanceMat() {

    int idx = 0;
    const float rad = pcisph.radius;

    if (pcisph.drawWall) {
        for (unsigned int i = 0; i < pcisph.numFluidParticles; i++) {
            instWorlds[idx] = glm::translate(glm::mat4(1.0f), pcisph.fluidParticles.pos[i]);
            instWorlds[idx] = glm::scale(instWorlds[idx], glm::vec3(rad));
            idx++;
        }
        for (unsigned int i = 0; i < pcisph.numWallParticles; i++) {
            instWorlds[idx] = glm::translate(glm::mat4(1.0f), pcisph.wallParticles.pos[i]);
            instWorlds[idx] = glm::scale(instWorlds[idx], glm::vec3(rad));
            idx++;
        }
    }
    else {
        for (unsigned int i = 0; i < pcisph.numFluidParticles; i++) {
            instWorlds[idx] = glm::translate(glm::mat4(1.0f), pcisph.fluidParticles.pos[i]);
            instWorlds[idx] = glm::scale(instWorlds[idx], glm::vec3(rad));
            idx++;

        }
    }


    //LOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOG
    //LOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOG
    //LOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOG
    //LOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOG
    //LOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOGLOG

        // TODO particle.h 참조. 여기 써놓은 것을 해야 log를 찍을 수 있고, 잘 parpticle이 박혔는지 확인할 수 있음. 그리고 그 부분이 아마 핵심부분임.
        /*
        if (pcisph.drawWall) {
            for (unsigned int i = 0; i < pcisph.numParticles; i++)
                std::cout << pcisph.particles.pos[i].x<<", " << pcisph.particles.pos[i].y << ", " << pcisph.particles.pos[i].z << std::endl;
        }
        else{
            for (unsigned int i = 0; i < pcisph.numFluidParticles; i++)
                std::cout << pcisph.fluidParticles.pos[i].x << ", " << pcisph.fluidParticles.pos[i].y << ", " << pcisph.fluidParticles.pos[i].z << std::endl;

        }
        */
}

//================================================================================
//================================================================================
//================================================================================
//================================================================================
//================================================================================

/*
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
        
        std::cout << " 몇번 돌았는데?" << temptemptemp << std::endl;

        std::cout << " neighbor 갯수 : " << neighbors->size() <<"density Var : " << particles[particleIdx].densityVar << std::endl;

        if (neighbors->size() == 1)
            std::cout<<"이웃의 idx" << particles[particleIdx].actualIdx << std::endl;

        std::cout << particles[particleIdx].pos.x << ": " << particles[particleIdx].pos.y << " : " << particles[particleIdx].pos.z << std::endl;

        system("PAUSE");
    }
    // 이게 도는 시점이 predicted 된 녀석들 기준으로 돌고 있다고. 근데. 지금 dist 보면 position으로만 하잖아. px가 아니라. 그러니까 터지지.
    // 지금 actual cellIdx와 px에 있을 때 cellIdx가 다른거야.

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

*/

/*
// TODO water particle에 대해서만 해야되는 놈인지, wall particle 까지 고려해야 하는 놈인지 cut을 잘 해야돼.
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
        ////////////////////////////////////////////////////////////        자명하게 waterParticle만 바뀌는 물리량임.

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


//=============================================================== 아마 다 구현 되지 않았을까...






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


void __cubeBoundaryCellIdx() {          // z indexing.에 맞도록 정렬을 하고, 정렬된 particles array에서 어느 idx부터 어디까지가 어떤 cellIdx를 가지는 particle인지 다 저장해 놓는 코드.

    particleIndices.clear();
    mortonCount.clear(); // 어느 idx부터 d어느 idx까지.
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
            idxMap.insert(std::pair<unsigned int, unsigned int>(i, invIdx++));              // morton code cell index를 넣어주면. 그 cell idx가 morton counter의 몇번째 index인지를 저장하는 map  
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

*/