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

extern const unsigned int SCR_WIDTH = 1280;
extern const unsigned int SCR_HEIGHT = 720;
Camera cam(glm::vec3(0.0f, 7.0f, 7.0f));

float lastX, lastY;
bool isFirstMove = true;


extern const float deltaTime = 1/60.0f;
const unsigned int numParticle = 50*50*50 ;
const glm::vec3 acc_g = glm::vec3(0.0f,-9.8f,0.0f);

const float restDensity = 997.0f;
const float radius = 0.03f;
const float particleMass = restDensity*4.0f / 3.0f * 3.141592f * radius * radius * radius;
const float coreRad = radius * 6.0f; //  (= grid side length, h )


Particle* particles = NULL;
vector<unsigned int> mortonCount;
std::map<unsigned int, unsigned int> idxMap;

vector<unsigned int>* neighborIdices;


// 정육면체 grid 가정함.
const float cubicBoundaryLen = 5.0f;
unsigned int nGridDivisoin = (unsigned int)(cubicBoundaryLen / coreRad);

void outBoundarySolution(Particle*);


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
    neighborIdices = new std::vector<unsigned int>[nGridDivisoin * nGridDivisoin * nGridDivisoin];
    gridInit(neighborIdices);


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

    glfwTerminate();
    return 0;
}


void particleInit(int mode) {

    particles = new Particle[numParticle];

    if (mode == 0) {
        int numRow = pow(numParticle, 1/3.0f);

        const float sideLen = 3.0f;

        for (unsigned int i = 0; i < numRow; i++) {
            for (unsigned int j = 0; j < numRow; j++) {
                for (unsigned int k = 0; k < numRow; k++) {
                    
                    particles[i * numRow * numRow + numRow * j + k].pos = glm::vec3(sideLen / 2.0f - sideLen * (i / (float)numRow)
                        ,sideLen - sideLen * (j / (float)numRow) + 3.9f
                    ,sideLen / 2.0f - sideLen * (k / (float)numRow));
                    
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

void updateParticles() {
    // 2 TODOexternal force, boundary...

    //1 set grid. particles sorting. offset is saved in mortonCount, mortonIdx
    __cubeBoundaryCellIdx(); 

    std::vector<Particle>* neighbors = new std::vector<Particle>();


    unsigned int currentCellIdx = -1;
    unsigned int prevCellIdx = -1;
    for (unsigned int particleIdx = 0; particleIdx < numParticle; particleIdx++) {

        // neighbor search할건데. particle이 정렬해놨음. 그러니까. 이전 particle과 같은 cellIdx에 속해있는 녀석에 대해서 이웃들을 구할 때 다시 neighbor 찾을 필요 없음.
        currentCellIdx = particles[particleIdx].cellIdx;
        if (prevCellIdx != currentCellIdx) {
            neighbors->clear();
            neighborSearch(neighbors, particleIdx);
        }

        particles[particleIdx].force_vis = glm::vec3(0);

        particles[particleIdx].force_ext = glm::vec3(0);// 2 TODOexternal force, boundary...

        particles[particleIdx].pressure = 0.0f;
        particles[particleIdx].force_p = glm::vec3(0);

        prevCellIdx = currentCellIdx;
    }

    for (unsigned int particleIdx = 0; particleIdx < numParticle; particleIdx++) {
        glm::vec3 acc = particles[particleIdx].force_ext + particles[particleIdx].force_p + particles[particleIdx].force_vis;
        acc /= particleMass;
        particles[particleIdx].vel += (acc + acc_g) * deltaTime;
        particles[particleIdx].pos+= particles[particleIdx].vel* deltaTime;
        outBoundarySolution(&particles[particleIdx]);
    }


    delete neighbors;
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

    unsigned int* mortonCounter = new unsigned[nGridDivisoin * nGridDivisoin * nGridDivisoin]; // temp
    unsigned int outBoundardyCount = 0;

    for (unsigned int i = 0; i < nGridDivisoin * nGridDivisoin * nGridDivisoin; i++)
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

    for (unsigned int i = 0; i < nGridDivisoin * nGridDivisoin * nGridDivisoin; i++)
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

    unsigned int idx = 0;
    for (unsigned int i = 0; i < nGridDivisoin; i++) {
        for (unsigned int j = 0; j < nGridDivisoin; j++) {
            for (unsigned int k = 0; k < nGridDivisoin; k++) {

                glm::ivec3 coord;


                for (int ii = -1; ii < 2; ii++) {
                    for (int jj = -1; jj < 2; jj++) {
                        for (int kk = -1; kk < 2; kk++) {

                            coord = glm::ivec3(i + ii, j + jj, k + kk);
                            unsigned int cellIdx = __cubeMorton(coord.x, coord.y, coord.z);

                            if (cellIdx != (unsigned int)(-1))
                                indices[idx].push_back(cellIdx);
                        }
                    }
                }

                idx++;
            }
        }
    }
    
}


// TODO 이거 한번 어떤 cellIdx에 대해서 구한 neighbor들은 다시 쓸 수 있을거니까. 어떻게든 여기서 시간 조금 더 줄일 수도 있을 것 같음.
void neighborSearch(std::vector<Particle>* neighbors,unsigned int idx) {
    for (unsigned int i = 0; i < neighborIdices[particles[idx].cellIdx].size(); i++) {
        unsigned int cidx = neighborIdices[particles[idx].cellIdx][i];
        cidx = idxMap[cidx];
        unsigned int startPoint = mortonCount[cidx];
        unsigned int endPoint = (cidx>=mortonCount.size()-1)? numParticle :mortonCount[cidx + 1];

        for (unsigned int j = startPoint; j < endPoint; j++) {
            neighbors->push_back(particles[j]);
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

unsigned int __cubeMorton(int x, int y, int z) { // TODO 이러면 된다는데... 비트연산자, 저거 16진수도 뭔지 잘 모르겠고. https://stackoverflow.com/questions/1024754/how-to-compute-a-3d-morton-number-interleave-the-bits-of-3-ints

    x = (x >= 0) ? x : -1;
    y = (y >= 0) ? y : -1;
    z = (z >= 0) ? z : -1;

    if ((x < 0) || (y < 0) || (z < 0)) // out of bounary
        return -1;

    static const unsigned int B[] = { 0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF };
    static const unsigned int S[] = { 1, 2, 4, 8 };

    x = (x | (x << S[3])) & B[3];
    x = (x | (x << S[2])) & B[2];
    x = (x | (x << S[1])) & B[1];
    x = (x | (x << S[0])) & B[0];

    y = (y | (y << S[3])) & B[3];
    y = (y | (y << S[2])) & B[2];
    y = (y | (y << S[1])) & B[1];
    y = (y | (y << S[0])) & B[0];

    z = (z | (z << S[3])) & B[3];
    z = (z | (z << S[2])) & B[2];
    z = (z | (z << S[1])) & B[1];
    z = (z | (z << S[0])) & B[0];

    return x | (y << 1) | (z<<2);

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
