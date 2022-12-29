#include<iostream>
#include<vector>
#include<map>
#include<cmath>

#include<glad/glad.h>
#include<GLFW/glfw3.h>
#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"


#include"camera.h"
#include"shader.h"
#include"model.h"
#include "particle.h"



void processInput(GLFWwindow* window);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void __mapBuffer(unsigned int& vbo, void* data, unsigned int size);

void updateParticles();
void particleInit(int mode);
void gridInit(std::vector<unsigned int>*);
void instanceMat(glm::mat4*);
void neighborSearch(std::vector<Particle>* neighbors,unsigned int Idx);


void __cubeBoundaryCellIdx();
glm::ivec3 __cubeCellCoord(glm::vec3 pos);
unsigned int __cubeMorton(int x, int y, int z);
inline uint64_t mortonEncode_LUT(unsigned int x, unsigned int y, unsigned int z);
void outBoundarySolution(Particle*);


// TODO 안쓰이는거 지우기.
float calcDensity(unsigned int particleIdx);
glm::vec3 forceVis(unsigned int particleIdx);
glm::vec3 forceSurfaceTension(unsigned int particleIdx);

unsigned int __nthDigit(unsigned int nGridDivisoin);



extern const unsigned int SCR_WIDTH = 1280;
extern const unsigned int SCR_HEIGHT = 720;
Camera cam(glm::vec3(0.0f, 7.0f, -7.0f));

float lastX, lastY;
bool isFirstMove = true;


extern const float deltaTime = 1/90.0f;
const unsigned int numParticle = 30*30*30 ;
const glm::vec3 acc_g = glm::vec3(0.0f,-9.8f,0.0f);

const float restDensity = 900.0f;
const float radius = 0.03f;
const float particleMass = restDensity * 4.0f / 3.0f * 3.141592f * radius * radius * radius;
const float coreRad = radius * 6.0f; //  (= grid side length, h )

const float linearVisc = 0.25f;
const float quadVisc = .5f;



Particle* particles = NULL;
vector<unsigned int> mortonCount;
std::map<unsigned int, unsigned int> idxMap;

vector<unsigned int>* neighborIdices;
std::vector<Particle>* neighbors;


// 정육면체 grid 가정함.
const float cubicBoundaryLen = 5.0f;
unsigned int nGridDivisoin = (unsigned int)(cubicBoundaryLen / coreRad);
unsigned int upperNGridDiv;

//viscosity

//surface tension

// predict - correct
const float eta = 0.01f;
const unsigned int minIter = 10;


int main()
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "PCISPH", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
    glEnable(GL_DEPTH_TEST);

//    glfwSetCursorPosCallback(window, mouse_callback);
//    glfwSetScrollCallback(window, scroll_callback);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);


    //Sphere render setting

    Shader particleShader("resources/shader/paricle_vs.txt", "resources/shader/paricle_fs.txt");
    Model sphere("resources/objects/sphere.obj");

    particleShader.use();
    glm::mat4 viewMat = glm::lookAt(cam.position, glm::vec3(0.0f), cam.up);
    glm::mat4 projMat = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
    glUniformMatrix4fv(glGetUniformLocation(particleShader.ID, "view"), 1, GL_FALSE, &viewMat[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(particleShader.ID, "proj"), 1, GL_FALSE, &projMat[0][0]);
    glUseProgram(0);

    // particle init
    particleInit(0);
    glm::mat4* instWorlds = new glm::mat4[numParticle];
    instanceMat(instWorlds);

    //Sphere render setting
    unsigned int instVBO;
    for (unsigned int i = 0; i < sphere.meshes.size(); i++) {

        glBindVertexArray(sphere.meshes[i].VAO);

        glGenBuffers(1, &instVBO);
        glBindBuffer(GL_ARRAY_BUFFER, instVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(glm::mat4) * numParticle, &instWorlds[0][0], GL_DYNAMIC_DRAW);

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

    //grid Init
    upperNGridDiv = std::pow(2,__nthDigit(nGridDivisoin));
    neighborIdices = new std::vector<unsigned int>[upperNGridDiv*upperNGridDiv*upperNGridDiv]; // nGridDiv^3 만큼이 아니라. nGridDiv의 2진수형태의 올림만큼 잡아야 됨. ngriddiv= 27 -> 32^3개 잡아야 됨.
    gridInit(neighborIdices);

    neighbors = new std::vector<Particle>();


    /*/grid test==============================================================================================================================================================================
    std::cout << nGridDivisoin * nGridDivisoin * nGridDivisoin << std::endl;
        for (unsigned int i = 0; i < nGridDivisoin; i++) {
            for (unsigned int j = 0; j < nGridDivisoin; j++) {
                for (unsigned int k = 0; k < nGridDivisoin; k++) {

                    unsigned int cellIDX = __cubeMorton(i, j, k);
                    std::cout << "cellIdx : " << i<< " , " << j << " , " << k << " : " << cellIDX;
                    std::cout << std::endl;
                }
            }
        }

    system("PAUSE"); 

    //==============================================================================================================================================================================
    */

    while (!glfwWindowShouldClose(window))
    {
        std::cout << "loop start =====================================================================================================================" << std::endl;
        processInput(window);

        glClearColor(0.9f, 0.9f, 0.9f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // state update
        updateParticles();
        instanceMat(instWorlds);
        __mapBuffer(instVBO, instWorlds, numParticle* sizeof(glm::mat4));

        //rendering
        particleShader.use();
        for (unsigned int i = 0; i < sphere.meshes.size(); i++) {
            glBindVertexArray(sphere.meshes[i].VAO);
            glDrawElementsInstanced(GL_TRIANGLES, sphere.meshes[i].indices.size(), GL_UNSIGNED_INT, 0, numParticle);
            glBindVertexArray(0);
        }
        glUseProgram(0);

        glfwSwapBuffers(window);
        glfwPollEvents();

        std::cout << "loop end =====================================================================================================================" << std::endl;
    }

    delete[] particles;
    delete[] instWorlds;
    delete neighbors;

    glfwTerminate();
    return 0;
}


void updateParticles() {

    //1 set grid. particles sorting. offset is saved in mortonCount, mortonIdx
    std::cout << " asdf " << std::endl;
    __cubeBoundaryCellIdx(); 

    std::cout << " asdf1 " << std::endl;
    unsigned int currentCellIdx = -1;
    unsigned int prevCellIdx = -1;
    for (unsigned int particleIdx = 0; particleIdx < numParticle; particleIdx++) {

        // neighbor search할건데. particle이 정렬해놨음. 그러니까. 이전 particle과 같은 cellIdx에 속해있는 녀석에 대해서 이웃들을 구할 때 다시 neighbor 찾을 필요 없음.
        currentCellIdx = particles[particleIdx].cellIdx;
        if (prevCellIdx != currentCellIdx) {
            neighbors->clear();
            neighborSearch(neighbors, particleIdx);
        }

        // gravity는 predict correct 할때 그냥 알아서 계산해주는 방식으로 ㄱㄱ
        particles[particleIdx].force_vis = forceVis(particleIdx);
        particles[particleIdx].force_ext = forceSurfaceTension(particleIdx);


        particles[particleIdx].pressure = 0.0f;
        particles[particleIdx].force_p = glm::vec3(0);

        prevCellIdx = currentCellIdx;
    }

    //predictive, corrective
    for (unsigned int particleIdx = 0; particleIdx < numParticle; particleIdx++) {

    }

    for (unsigned int particleIdx = 0; particleIdx < numParticle; particleIdx++) {
        glm::vec3 acc = particles[particleIdx].force_ext + particles[particleIdx].force_p + particles[particleIdx].force_vis;
        acc /= particleMass;
        particles[particleIdx].vel += (acc + acc_g) * deltaTime;
        particles[particleIdx].pos+= particles[particleIdx].vel* deltaTime;
        outBoundarySolution(&particles[particleIdx]);
    }


}

const float kernelConst = 2.0f / (coreRad * coreRad * coreRad * 1.772454f); // root pi = 1.7724538...

glm::vec3 forceVis(unsigned int particleIdx) { //TODO linear Visc, quadVisc, kernel 조절하기.

    Particle target = particles[particleIdx];
    glm::vec3 netForceVis = glm::vec3(0);  

    for (auto iter = neighbors->begin(); iter != neighbors->end(); iter++) {
        glm::vec3 dist = ((*iter).pos - target.pos);
        float distLen = glm::length(dist);

        float distRatio = distLen / coreRad;

        float mu = glm::dot((*iter).vel - target.vel, dist);

        if (mu < FLT_EPSILON && distLen - coreRad < FLT_EPSILON) {
            float pipi = 0;
            mu /= (0.01f * coreRad + distLen * distRatio);
            pipi = -linearVisc * mu + quadVisc * mu * mu;
            pipi /= (target.density + (*iter).density);
            glm::vec3 grad_Wvis = 2.0f * pipi * dist * kernelConst / std::pow(2.718282f, distRatio);
            netForceVis -= grad_Wvis;
        }
    }

    return netForceVis;
}


glm::vec3 forceSurfaceTension(unsigned int particleIdx) {
    return glm::vec3(0);
}

float calcDensity(unsigned int particleIdx) { // 주변의 밀도 구할때 쓰임 ㅇㅇ
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


void outBoundarySolution(Particle* particle) {
    glm::ivec3 cellIdx = __cubeCellCoord(particle->pos);
    glm::ivec3 flag;

    flag.x = (cellIdx.x == -1) ? -1 : (cellIdx.x >= nGridDivisoin) ? 1 : 0;
    flag.z = (cellIdx.z == -1) ? -1 : (cellIdx.z >= nGridDivisoin) ? 1 : 0;
    flag.y = (cellIdx.y == -1) ? -1 : (cellIdx.y >= nGridDivisoin) ? 1 : 0;


    if (flag.x == -1) {
        particle->vel.x *= -1;
        particle->pos.x += 2*(-cubicBoundaryLen / 2.0f - particle->pos.x);
    }
    if (flag.x == 1) {
        particle->vel.x *= -1;
        particle->pos.x += 2 * (cubicBoundaryLen / 2.0f - particle->pos.x);
    }

    if (flag.z == -1) {
        particle->vel.z *= -1;
        particle->pos.z += 2 * (-cubicBoundaryLen / 2.0f - particle->pos.z);
    }
    if (flag.z == 1) {
        particle->vel.z *= -1;
        particle->pos.z += 2 * (cubicBoundaryLen / 2.0f - particle->pos.z);
    }

    if (flag.y == -1) {
        particle->vel.y *= -1;
        particle->pos.y += 2 * (-1.0f - particle->pos.y);
    }
    if (flag.y == 1) {
        particle->vel.y *= -1;
        particle->pos.y += 2 * (-1.0f + cubicBoundaryLen - particle->pos.y);
    }

}

void __cubeBoundaryCellIdx() {          // z indexing.에 맞도록 정렬을 하고, 정렬된 particles array에서 어느 idx부터 어디까지가 어떤 cellIdx를 가지는 particle인지 다 저장해 놓는 코드.

    mortonCount.clear(); // 어느 idx부터 d어느 idx까지.

    unsigned int* mortonCounter = new unsigned[upperNGridDiv* upperNGridDiv* upperNGridDiv]; // temp
    unsigned int outBoundardyCount = 0;

    for (unsigned int i = 0; i < upperNGridDiv * upperNGridDiv * upperNGridDiv; i++)
        mortonCounter[i] = 0;


    for (unsigned int i = 0; i < numParticle; i++) {
        glm::ivec3 cellCoord = __cubeCellCoord(particles[i].pos);

        particles[i].cellIdx = __cubeMorton(cellCoord.x, cellCoord.y, cellCoord.z);

        if (particles[i].cellIdx == (unsigned int)( - 1))
            outBoundardyCount++;
        else
            mortonCounter[particles[i].cellIdx]++;
        
    }

    if (outBoundardyCount > 0) {
        std::cout << "out boundary particle exist" << std::endl;
        system("PAUSE");
    }

    unsigned int invIdx = 0;
    
    idxMap.clear();

    for (unsigned int i = 0; i < upperNGridDiv * upperNGridDiv * upperNGridDiv; i++)
        if (mortonCounter[i] > 0) {
            mortonCount.push_back(mortonCounter[i]);
            idxMap.insert(std::pair<unsigned int, unsigned int>(i, invIdx++));
        }
    delete[] mortonCounter;


    if(mortonCount.size()>0)
        for (auto iter = mortonCount.begin()+1; iter != mortonCount.end(); iter++) 
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
            particles[i++] = *iteriter;
        }
    }

}


void gridInit(std::vector<unsigned int>* indices) {
    // 문제점이 있음. morton code를 쓸 때. 예를들어 3개의 축의 cell의 갯수가 27개씩 있다고 하자. 27*27*27. 그러면 사실 나올 수 있는 index의 경우의는 27^3 이 맞음. 근데 morton code로 encoding했을 때 나오는 index값은 0~32*32*32의 범위를 가진다.
    // 그래서 cell idx를 만들 때 조금 고려해야 할 것들이 있음.
    for (unsigned int i = 0; i < nGridDivisoin; i++) {
        for (unsigned int j = 0; j < nGridDivisoin; j++) {
            for (unsigned int k = 0; k < nGridDivisoin; k++) {

                unsigned int cidx = __cubeMorton(i, j, k);

                for (int ii = -1; ii < 2; ii++) {
                    for (int jj = -1; jj < 2; jj++) {
                        for (int kk = -1; kk < 2; kk++) {

                            glm::ivec3 coord = glm::ivec3(i + ii, j + jj, k + kk);
                            unsigned int neighborCell = __cubeMorton(coord.x, coord.y, coord.z);

                            if (neighborCell != (unsigned int)(-1))
                                indices[cidx].push_back(neighborCell);
                        }
                    }
                }
            }
        }
    }
    
}

void neighborSearch(std::vector<Particle>* neighbors,unsigned int idx) {
    for (unsigned int i = 0; i < neighborIdices[particles[idx].cellIdx].size(); i++) {

        unsigned int cidx = neighborIdices[particles[idx].cellIdx][i];


        if (idxMap.count(cidx) != 0) {
            cidx = idxMap[cidx];
            unsigned int startPoint = mortonCount[cidx];
            unsigned int endPoint = (cidx >= mortonCount.size() - 1) ? numParticle : mortonCount[cidx + 1];

            for (unsigned int j = startPoint; j < endPoint; j++) {
                neighbors->push_back(particles[j]);
            }
        }
    }
}


glm::ivec3 __cubeCellCoord(glm::vec3 pos) {

    float x = ((pos.x + cubicBoundaryLen / 2.0f) / cubicBoundaryLen * nGridDivisoin);
    float y = ((pos.y + 1.0f) / cubicBoundaryLen * nGridDivisoin);
    float z = ((pos.z + cubicBoundaryLen / 2.0f) / cubicBoundaryLen * nGridDivisoin);

    x = (x + FLT_EPSILON <= 0) ? -1.0f : x;
    z = (z + FLT_EPSILON <= 0) ? -1.0f : z;
    y = (y + FLT_EPSILON <= 0) ? -1.0f : y;

    return glm::ivec3((int)x, (int)y, (int)z);
}

unsigned int __cubeMorton(int x, int y, int z) {
    x = (x >= 0) ? x : -1;
    y = (y >= 0) ? y : -1;
    z = (z >= 0) ? z : -1;

    if ((x < 0) || (y < 0) || (z < 0)) // out of bounary
        return -1;


    return  mortonEncode_LUT(x, y, z);
}



void __mapBuffer(unsigned int& vbo, void* data, unsigned int size) {
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    void* ptr = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

    memcpy(ptr, data, size);

    glUnmapBuffer(GL_ARRAY_BUFFER);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        cam.processKeyboard(FORWARD);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        cam.processKeyboard(LEFT);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        cam.processKeyboard(BACKWARD);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        cam.processKeyboard(RIGHT);
}
void mouse_callback(GLFWwindow* window, double xpos, double ypos) {

    if (isFirstMove) {
        isFirstMove = false;
        lastX = xpos;
        lastY = ypos;
    }

    float xoffset = static_cast<float>(xpos) - lastX;
    float yoffset = static_cast<float>(ypos) - lastY;
    lastX = xpos;
    lastY = ypos;

    cam.processMovement(xoffset, yoffset);
}
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    cam.processMouseScroll(yoffset);
}

void particleInit(int mode) {

    particles = new Particle[numParticle];

    if (mode == 0) {
        int numRow = pow(numParticle, 1 / 3.0f);

        const float sideLen = 3.0f;

        for (unsigned int i = 0; i < numRow; i++) {
            for (unsigned int j = 0; j < numRow; j++) {
                for (unsigned int k = 0; k < numRow; k++) {

                    particles[i * numRow * numRow + numRow * j + k].pos = glm::vec3(sideLen / 2.0f - sideLen * (i / (float)numRow)
                        , sideLen - sideLen * (j / (float)numRow)
                        , sideLen / 2.0f - sideLen * (k / (float)numRow));
                    particles[i * numRow * numRow + numRow * j + k].density = restDensity;

                }
            }
        }
    }
}
void instanceMat(glm::mat4* instWorlds) {
    for (unsigned int i = 0; i < numParticle; i++) {
        instWorlds[i] = glm::translate(glm::mat4(1.0f), particles[i].pos);
        instWorlds[i] = glm::scale(instWorlds[i], glm::vec3(0.03f));
    }
}


static const uint32_t morton256_x[256] =
{
0x00000000,
0x00000001, 0x00000008, 0x00000009, 0x00000040, 0x00000041, 0x00000048, 0x00000049, 0x00000200,
0x00000201, 0x00000208, 0x00000209, 0x00000240, 0x00000241, 0x00000248, 0x00000249, 0x00001000,
0x00001001, 0x00001008, 0x00001009, 0x00001040, 0x00001041, 0x00001048, 0x00001049, 0x00001200,
0x00001201, 0x00001208, 0x00001209, 0x00001240, 0x00001241, 0x00001248, 0x00001249, 0x00008000,
0x00008001, 0x00008008, 0x00008009, 0x00008040, 0x00008041, 0x00008048, 0x00008049, 0x00008200,
0x00008201, 0x00008208, 0x00008209, 0x00008240, 0x00008241, 0x00008248, 0x00008249, 0x00009000,
0x00009001, 0x00009008, 0x00009009, 0x00009040, 0x00009041, 0x00009048, 0x00009049, 0x00009200,
0x00009201, 0x00009208, 0x00009209, 0x00009240, 0x00009241, 0x00009248, 0x00009249, 0x00040000,
0x00040001, 0x00040008, 0x00040009, 0x00040040, 0x00040041, 0x00040048, 0x00040049, 0x00040200,
0x00040201, 0x00040208, 0x00040209, 0x00040240, 0x00040241, 0x00040248, 0x00040249, 0x00041000,
0x00041001, 0x00041008, 0x00041009, 0x00041040, 0x00041041, 0x00041048, 0x00041049, 0x00041200,
0x00041201, 0x00041208, 0x00041209, 0x00041240, 0x00041241, 0x00041248, 0x00041249, 0x00048000,
0x00048001, 0x00048008, 0x00048009, 0x00048040, 0x00048041, 0x00048048, 0x00048049, 0x00048200,
0x00048201, 0x00048208, 0x00048209, 0x00048240, 0x00048241, 0x00048248, 0x00048249, 0x00049000,
0x00049001, 0x00049008, 0x00049009, 0x00049040, 0x00049041, 0x00049048, 0x00049049, 0x00049200,
0x00049201, 0x00049208, 0x00049209, 0x00049240, 0x00049241, 0x00049248, 0x00049249, 0x00200000,
0x00200001, 0x00200008, 0x00200009, 0x00200040, 0x00200041, 0x00200048, 0x00200049, 0x00200200,
0x00200201, 0x00200208, 0x00200209, 0x00200240, 0x00200241, 0x00200248, 0x00200249, 0x00201000,
0x00201001, 0x00201008, 0x00201009, 0x00201040, 0x00201041, 0x00201048, 0x00201049, 0x00201200,
0x00201201, 0x00201208, 0x00201209, 0x00201240, 0x00201241, 0x00201248, 0x00201249, 0x00208000,
0x00208001, 0x00208008, 0x00208009, 0x00208040, 0x00208041, 0x00208048, 0x00208049, 0x00208200,
0x00208201, 0x00208208, 0x00208209, 0x00208240, 0x00208241, 0x00208248, 0x00208249, 0x00209000,
0x00209001, 0x00209008, 0x00209009, 0x00209040, 0x00209041, 0x00209048, 0x00209049, 0x00209200,
0x00209201, 0x00209208, 0x00209209, 0x00209240, 0x00209241, 0x00209248, 0x00209249, 0x00240000,
0x00240001, 0x00240008, 0x00240009, 0x00240040, 0x00240041, 0x00240048, 0x00240049, 0x00240200,
0x00240201, 0x00240208, 0x00240209, 0x00240240, 0x00240241, 0x00240248, 0x00240249, 0x00241000,
0x00241001, 0x00241008, 0x00241009, 0x00241040, 0x00241041, 0x00241048, 0x00241049, 0x00241200,
0x00241201, 0x00241208, 0x00241209, 0x00241240, 0x00241241, 0x00241248, 0x00241249, 0x00248000,
0x00248001, 0x00248008, 0x00248009, 0x00248040, 0x00248041, 0x00248048, 0x00248049, 0x00248200,
0x00248201, 0x00248208, 0x00248209, 0x00248240, 0x00248241, 0x00248248, 0x00248249, 0x00249000,
0x00249001, 0x00249008, 0x00249009, 0x00249040, 0x00249041, 0x00249048, 0x00249049, 0x00249200,
0x00249201, 0x00249208, 0x00249209, 0x00249240, 0x00249241, 0x00249248, 0x00249249
};

// pre-shifted table for Y coordinates (1 bit to the left)
static const uint32_t morton256_y[256] = {
0x00000000,
0x00000002, 0x00000010, 0x00000012, 0x00000080, 0x00000082, 0x00000090, 0x00000092, 0x00000400,
0x00000402, 0x00000410, 0x00000412, 0x00000480, 0x00000482, 0x00000490, 0x00000492, 0x00002000,
0x00002002, 0x00002010, 0x00002012, 0x00002080, 0x00002082, 0x00002090, 0x00002092, 0x00002400,
0x00002402, 0x00002410, 0x00002412, 0x00002480, 0x00002482, 0x00002490, 0x00002492, 0x00010000,
0x00010002, 0x00010010, 0x00010012, 0x00010080, 0x00010082, 0x00010090, 0x00010092, 0x00010400,
0x00010402, 0x00010410, 0x00010412, 0x00010480, 0x00010482, 0x00010490, 0x00010492, 0x00012000,
0x00012002, 0x00012010, 0x00012012, 0x00012080, 0x00012082, 0x00012090, 0x00012092, 0x00012400,
0x00012402, 0x00012410, 0x00012412, 0x00012480, 0x00012482, 0x00012490, 0x00012492, 0x00080000,
0x00080002, 0x00080010, 0x00080012, 0x00080080, 0x00080082, 0x00080090, 0x00080092, 0x00080400,
0x00080402, 0x00080410, 0x00080412, 0x00080480, 0x00080482, 0x00080490, 0x00080492, 0x00082000,
0x00082002, 0x00082010, 0x00082012, 0x00082080, 0x00082082, 0x00082090, 0x00082092, 0x00082400,
0x00082402, 0x00082410, 0x00082412, 0x00082480, 0x00082482, 0x00082490, 0x00082492, 0x00090000,
0x00090002, 0x00090010, 0x00090012, 0x00090080, 0x00090082, 0x00090090, 0x00090092, 0x00090400,
0x00090402, 0x00090410, 0x00090412, 0x00090480, 0x00090482, 0x00090490, 0x00090492, 0x00092000,
0x00092002, 0x00092010, 0x00092012, 0x00092080, 0x00092082, 0x00092090, 0x00092092, 0x00092400,
0x00092402, 0x00092410, 0x00092412, 0x00092480, 0x00092482, 0x00092490, 0x00092492, 0x00400000,
0x00400002, 0x00400010, 0x00400012, 0x00400080, 0x00400082, 0x00400090, 0x00400092, 0x00400400,
0x00400402, 0x00400410, 0x00400412, 0x00400480, 0x00400482, 0x00400490, 0x00400492, 0x00402000,
0x00402002, 0x00402010, 0x00402012, 0x00402080, 0x00402082, 0x00402090, 0x00402092, 0x00402400,
0x00402402, 0x00402410, 0x00402412, 0x00402480, 0x00402482, 0x00402490, 0x00402492, 0x00410000,
0x00410002, 0x00410010, 0x00410012, 0x00410080, 0x00410082, 0x00410090, 0x00410092, 0x00410400,
0x00410402, 0x00410410, 0x00410412, 0x00410480, 0x00410482, 0x00410490, 0x00410492, 0x00412000,
0x00412002, 0x00412010, 0x00412012, 0x00412080, 0x00412082, 0x00412090, 0x00412092, 0x00412400,
0x00412402, 0x00412410, 0x00412412, 0x00412480, 0x00412482, 0x00412490, 0x00412492, 0x00480000,
0x00480002, 0x00480010, 0x00480012, 0x00480080, 0x00480082, 0x00480090, 0x00480092, 0x00480400,
0x00480402, 0x00480410, 0x00480412, 0x00480480, 0x00480482, 0x00480490, 0x00480492, 0x00482000,
0x00482002, 0x00482010, 0x00482012, 0x00482080, 0x00482082, 0x00482090, 0x00482092, 0x00482400,
0x00482402, 0x00482410, 0x00482412, 0x00482480, 0x00482482, 0x00482490, 0x00482492, 0x00490000,
0x00490002, 0x00490010, 0x00490012, 0x00490080, 0x00490082, 0x00490090, 0x00490092, 0x00490400,
0x00490402, 0x00490410, 0x00490412, 0x00490480, 0x00490482, 0x00490490, 0x00490492, 0x00492000,
0x00492002, 0x00492010, 0x00492012, 0x00492080, 0x00492082, 0x00492090, 0x00492092, 0x00492400,
0x00492402, 0x00492410, 0x00492412, 0x00492480, 0x00492482, 0x00492490, 0x00492492
};

// Pre-shifted table for z (2 bits to the left)
static const uint32_t morton256_z[256] = {
0x00000000,
0x00000004, 0x00000020, 0x00000024, 0x00000100, 0x00000104, 0x00000120, 0x00000124, 0x00000800,
0x00000804, 0x00000820, 0x00000824, 0x00000900, 0x00000904, 0x00000920, 0x00000924, 0x00004000,
0x00004004, 0x00004020, 0x00004024, 0x00004100, 0x00004104, 0x00004120, 0x00004124, 0x00004800,
0x00004804, 0x00004820, 0x00004824, 0x00004900, 0x00004904, 0x00004920, 0x00004924, 0x00020000,
0x00020004, 0x00020020, 0x00020024, 0x00020100, 0x00020104, 0x00020120, 0x00020124, 0x00020800,
0x00020804, 0x00020820, 0x00020824, 0x00020900, 0x00020904, 0x00020920, 0x00020924, 0x00024000,
0x00024004, 0x00024020, 0x00024024, 0x00024100, 0x00024104, 0x00024120, 0x00024124, 0x00024800,
0x00024804, 0x00024820, 0x00024824, 0x00024900, 0x00024904, 0x00024920, 0x00024924, 0x00100000,
0x00100004, 0x00100020, 0x00100024, 0x00100100, 0x00100104, 0x00100120, 0x00100124, 0x00100800,
0x00100804, 0x00100820, 0x00100824, 0x00100900, 0x00100904, 0x00100920, 0x00100924, 0x00104000,
0x00104004, 0x00104020, 0x00104024, 0x00104100, 0x00104104, 0x00104120, 0x00104124, 0x00104800,
0x00104804, 0x00104820, 0x00104824, 0x00104900, 0x00104904, 0x00104920, 0x00104924, 0x00120000,
0x00120004, 0x00120020, 0x00120024, 0x00120100, 0x00120104, 0x00120120, 0x00120124, 0x00120800,
0x00120804, 0x00120820, 0x00120824, 0x00120900, 0x00120904, 0x00120920, 0x00120924, 0x00124000,
0x00124004, 0x00124020, 0x00124024, 0x00124100, 0x00124104, 0x00124120, 0x00124124, 0x00124800,
0x00124804, 0x00124820, 0x00124824, 0x00124900, 0x00124904, 0x00124920, 0x00124924, 0x00800000,
0x00800004, 0x00800020, 0x00800024, 0x00800100, 0x00800104, 0x00800120, 0x00800124, 0x00800800,
0x00800804, 0x00800820, 0x00800824, 0x00800900, 0x00800904, 0x00800920, 0x00800924, 0x00804000,
0x00804004, 0x00804020, 0x00804024, 0x00804100, 0x00804104, 0x00804120, 0x00804124, 0x00804800,
0x00804804, 0x00804820, 0x00804824, 0x00804900, 0x00804904, 0x00804920, 0x00804924, 0x00820000,
0x00820004, 0x00820020, 0x00820024, 0x00820100, 0x00820104, 0x00820120, 0x00820124, 0x00820800,
0x00820804, 0x00820820, 0x00820824, 0x00820900, 0x00820904, 0x00820920, 0x00820924, 0x00824000,
0x00824004, 0x00824020, 0x00824024, 0x00824100, 0x00824104, 0x00824120, 0x00824124, 0x00824800,
0x00824804, 0x00824820, 0x00824824, 0x00824900, 0x00824904, 0x00824920, 0x00824924, 0x00900000,
0x00900004, 0x00900020, 0x00900024, 0x00900100, 0x00900104, 0x00900120, 0x00900124, 0x00900800,
0x00900804, 0x00900820, 0x00900824, 0x00900900, 0x00900904, 0x00900920, 0x00900924, 0x00904000,
0x00904004, 0x00904020, 0x00904024, 0x00904100, 0x00904104, 0x00904120, 0x00904124, 0x00904800,
0x00904804, 0x00904820, 0x00904824, 0x00904900, 0x00904904, 0x00904920, 0x00904924, 0x00920000,
0x00920004, 0x00920020, 0x00920024, 0x00920100, 0x00920104, 0x00920120, 0x00920124, 0x00920800,
0x00920804, 0x00920820, 0x00920824, 0x00920900, 0x00920904, 0x00920920, 0x00920924, 0x00924000,
0x00924004, 0x00924020, 0x00924024, 0x00924100, 0x00924104, 0x00924120, 0x00924124, 0x00924800,
0x00924804, 0x00924820, 0x00924824, 0x00924900, 0x00924904, 0x00924920, 0x00924924
};

inline uint64_t mortonEncode_LUT(unsigned int x, unsigned int y, unsigned int z) {
    uint64_t answer = 0;
    answer = morton256_z[(z >> 16) & 0xFF] | // we start by shifting the third byte, since we only look at the first 21 bits
        morton256_y[(y >> 16) & 0xFF] |
        morton256_x[(x >> 16) & 0xFF];
    answer = answer << 48 | morton256_z[(z >> 8) & 0xFF] | // shifting second byte
        morton256_y[(y >> 8) & 0xFF] |
        morton256_x[(x >> 8) & 0xFF];
    answer = answer << 24 |
        morton256_z[(z) & 0xFF] | // first byte
        morton256_y[(y) & 0xFF] |
        morton256_x[(x) & 0xFF];
    return answer;
}


unsigned int __nthDigit(unsigned int nGridDivisoin) {
    unsigned int digit = 0;
    while (nGridDivisoin != 0) {
        nGridDivisoin /= 2;
        digit++;
    }
    return digit;
}