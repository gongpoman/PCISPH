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
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void __mapBuffer(unsigned int& vbo, void* data, unsigned int size);

void updateParticles();
void particleInit(int mode);
void instanceMat(glm::mat4*);


void __cubeBoundaryCellIdx();
glm::ivec3 __cubeCellCoord(glm::vec3 pos);
unsigned int __cubeMorton(int x, int y, int z);

extern const unsigned int SCR_WIDTH = 1280;
extern const unsigned int SCR_HEIGHT = 720;
Camera cam(glm::vec3(3.0f, 7.0f, 7.0f));

float lastX, lastY;
bool isFirstMove = true;




extern const float deltaTime = 1/60.0f;
const unsigned int numParticle = 40*40*40;
const glm::vec3 acc_g = glm::vec3(0.0f,-9.8f,0.0f);

const float restDensity = 997.0f;
const float radius = 0.03f;
const float particleMass = restDensity*4.0f / 3.0f * 3.141592f * radius * radius * radius;
const float coreRad = radius * 6.0f; //  (= grid side length, h )

Particle* particles = NULL;
vector<unsigned int> mortonIdx;
vector<unsigned int> mortonCount;

// 정육면체 grid 가정함.
const float cubicBoundaryLen = 5.0f;
unsigned int nGridDivisoin = (unsigned int)(cubicBoundaryLen / coreRad);



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

    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);
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
                        ,sideLen - sideLen * (j / (float)numRow)
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

    //1 set grid. particles are sorted. offseted with mortonCounter
    __cubeBoundaryCellIdx(); 


    for (unsigned int particleIdx = 0; particleIdx < numParticle; particleIdx++) {

        particles[particleIdx].force_vis = glm::vec3(0);
        particles[particleIdx].force_ext = glm::vec3(0);// 2 TODOexternal force, boundary...
        particles[particleIdx].pressure = 0.0f;
        particles[particleIdx].force_p = glm::vec3(0);
    }

    for (unsigned int particleIdx = 0; particleIdx < numParticle; particleIdx++) {
        glm::vec3 acc = particles[particleIdx].force_ext + particles[particleIdx].force_p + particles[particleIdx].force_vis;
        acc /= particleMass;
        particles[particleIdx].vel += (acc + acc_g) * deltaTime;
        particles[particleIdx].pos+= particles[particleIdx].vel* deltaTime;
    }



}


void __cubeBoundaryCellIdx() {          // z indexing.에 맞도록 정렬을 하고, 정렬된 particles array에서 어느 idx부터 어디까지가 어떤 cellIdx를 가지는 particle인지 다 저장해 놓는 코드.

    mortonCount.clear(); // 어느 idx부터 d어느 idx까지.
    mortonIdx.clear(); // 위의 idx 범위가 나타내는게 어느 cellIdx인지.

    unsigned int* mortonCounter = new unsigned[nGridDivisoin * nGridDivisoin * nGridDivisoin]; // temp
    unsigned int outBounardyCount = 0;

    for (unsigned int i = 0; i < nGridDivisoin * nGridDivisoin * nGridDivisoin; i++)
        mortonCounter[i] = 0;


    for (unsigned int i = 0; i < numParticle; i++) {
        glm::ivec3 cellCoord = __cubeCellCoord(particles[i].pos);

        particles[i].cellIdx = __cubeMorton(cellCoord.x, cellCoord.y, cellCoord.z);

        if (particles[i].cellIdx == (unsigned int)( - 1))
            outBounardyCount++;
        else
            mortonCounter[particles[i].cellIdx]++;
        
    }

    unsigned int invIdx = 0;
    std::map<unsigned int, unsigned int> idxMap;

    for (unsigned int i = 0; i < nGridDivisoin * nGridDivisoin * nGridDivisoin; i++)
        if (mortonCounter[i] > 0) {
            mortonIdx.push_back(i);
            mortonCount.push_back(mortonCounter[i]);
            idxMap.insert(std::pair<unsigned int, unsigned int>(i, invIdx));
        }
    if (outBounardyCount > 0) {
        mortonCount.push_back(outBounardyCount);
        mortonIdx.push_back((unsigned int)(- 1));
        idxMap.insert(std::pair<unsigned int, unsigned int>( ((unsigned int)(-1)) , invIdx )) ;
    }
    delete[] mortonCounter;



    if(mortonCount.size()>0)
        for (auto iter = mortonCount.begin()+1; iter != mortonCount.end(); iter++) 
            *(iter) += *(iter - 1);

    std::vector<std::vector<Particle>> tempParticles;
    for (unsigned int i = 0; i < mortonIdx.size(); i++) 
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


glm::ivec3 __cubeCellCoord(glm::vec3 pos) {
    unsigned int x = (unsigned int)((pos.x + cubicBoundaryLen / 2.0f) / cubicBoundaryLen * nGridDivisoin);
    unsigned int y = (unsigned int)(pos.y / cubicBoundaryLen * nGridDivisoin);
    unsigned int z = (unsigned int)((pos.z + cubicBoundaryLen / 2.0f) / cubicBoundaryLen * nGridDivisoin);
    return glm::ivec3(x, y, z);
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

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
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
