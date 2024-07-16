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

extern Camera cam(glm::vec3(0.0f,1.8f, 3.f));


//PCISPH pcisph(glm::ivec3(25, 16, 25), glm::vec3(2.5f, 1.5f, 2.5f), false);
PCISPH pcisph(glm::ivec3(10, 15, 20), glm::vec3(1.5f, 1.8f, 1.f), false);



extern float lastX, lastY;
extern bool isFirstMove = true;

bool mouseEnabel = false;
bool keyboardEnable = false;

extern const float deltaTime = 1/180.0f;

glm::mat4* instWorlds;
glm::vec3* instVel;
float* instDen;

GLFWwindow* glInitialize();
void initialStateLog();
std::tuple<Shader*,Model*, unsigned int, unsigned int, unsigned int> renderInitialize();

int main(){
    
    static unsigned int frame = 0;

    //glInit
    GLFWwindow* window = glInitialize();
    if (!window)
        return -1;

    instWorlds = new glm::mat4[pcisph.numDrawParticles];
    instDen = new float[pcisph.numDrawParticles];
    instVel = new glm::vec3[pcisph.numDrawParticles];

    instanceMat();

    //particle render setting
    Shader* particleShader;
    Model* sphere;
    unsigned int instVBO;
    unsigned int instDenVBO;
    unsigned int instVelVBO;
    std::tie(particleShader, sphere,instVBO,instDenVBO,instVelVBO) = renderInitialize();

    // print scean Log
    initialStateLog();

    while (!glfwWindowShouldClose(window))
    {
        frame++;
        std::cout<<frame <<" frame start=======================================================================================" << std::endl;

        if(keyboardEnable)
            processInput(window);

        glClearColor(0.9f, 0.9f, 0.9f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        ////////////Update
        pcisph.update();
//        pcisph.nSearchTest();
        instanceMat();
        __mapBuffer(instVBO, instWorlds, pcisph.numDrawParticles* sizeof(glm::mat4));
        __mapBuffer(instVelVBO, instVel, pcisph.numDrawParticles* sizeof(glm::vec3));

        ////////////render 
        particleShader->use();
        
        if (keyboardEnable) {
            glm::mat4 viewMat = cam.getViewMatrix();
            glUniformMatrix4fv(glGetUniformLocation(particleShader->ID, "view"), 1, GL_FALSE, &viewMat[0][0]);
        }

        glUniform3f(glGetUniformLocation(particleShader->ID, "lightPos"), cam.position.x, cam.position.y, cam.position.z);
        // DRAW CALL
        #pragma omp parallel for
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
 
    }

    std::cout << frame << "frames rendered" << std::endl;

    //global
    delete[] instWorlds;
    delete[] instVel;
    delete[] instDen;

    // local 
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

    if (mouseEnabel) {
        glfwSetCursorPosCallback(window, mouse_callback);
        glfwSetScrollCallback(window, scroll_callback);
    }
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

std::tuple < Shader*, Model*, unsigned int, unsigned int, unsigned int > renderInitialize() {

    Shader* particleShader = new Shader("resources/shader/paricleL_vs.txt", "resources/shader/paricleL_fs.txt");
    Model* sphere = new Model("resources/objects/sphere.obj");

    cam.position += glm::vec3(pcisph.boundaryX / 2.0f, 0.0f, pcisph.boundaryZ / 2.0f);

    particleShader->use();
    glm::mat4 viewMat = glm::lookAt(cam.position, glm::vec3(0.0f, 0.7f, 0.0f) + glm::vec3(pcisph.boundaryX / 2.0f, 0.0f, pcisph.boundaryZ / 2.0f), cam.up);
    glm::mat4 projMat = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
    glUniformMatrix4fv(glGetUniformLocation(particleShader->ID, "view"), 1, GL_FALSE, &viewMat[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(particleShader->ID, "proj"), 1, GL_FALSE, &projMat[0][0]);
    glUniform3f(glGetUniformLocation(particleShader->ID, "lightPos"), cam.position.x, cam.position.y, cam.position.z);
    glUseProgram(0);

    unsigned int instVBO;

    unsigned int instDenVBO; //#Density CHECK

    unsigned int instVelVBO;


    for (unsigned int i = 0; i < sphere->meshes.size(); i++) {  // Definitely, sphere->meshes.size() = 1

        glBindVertexArray(sphere->meshes[i].VAO);

        glGenBuffers(1, &instVBO);
        glBindBuffer(GL_ARRAY_BUFFER, instVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(glm::mat4) * pcisph.numDrawParticles, &instWorlds[0][0], GL_DYNAMIC_DRAW);

        // 4*4 matrix∏¶ vertex shader 3~6ø° ≥—∞‹¡‹
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


        glGenBuffers(1, &instDenVBO);
        glBindBuffer(GL_ARRAY_BUFFER, instDenVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * pcisph.numDrawParticles, &instDen[0], GL_DYNAMIC_DRAW);

        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(float), (void*)0);
        glVertexAttribDivisor(2, 1);


        glGenBuffers(1, &instVelVBO);
        glBindBuffer(GL_ARRAY_BUFFER, instVelVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * pcisph.numDrawParticles, &instVel[0][0], GL_DYNAMIC_DRAW);

        glEnableVertexAttribArray(7);
        glVertexAttribPointer(7, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
        glVertexAttribDivisor(7, 1);



        glBindVertexArray(0);
    }

    return { particleShader,sphere,instVBO,instDenVBO,instVelVBO };
}
void instanceMat() {

    int idx = 0;
    const float rad = pcisph.radius ;

    if (pcisph.drawWall) {
#pragma omp parallel for
        for (unsigned int i = 0; i < pcisph.numFluidParticles; i++) {
            instWorlds[idx] = glm::translate(glm::mat4(1.0f), pcisph.fluidParticles.pos[i]);
            instWorlds[idx] = glm::scale(instWorlds[idx], glm::vec3(rad));

            instDen[idx] = pcisph.fluidParticles.density[i];

            instVel[idx] = pcisph.fluidParticles.vel[i];

            idx++;
        }
#pragma omp parallel for
        for (unsigned int i = 0; i < pcisph.numWallParticles; i++) {
            instWorlds[idx] = glm::translate(glm::mat4(1.0f), pcisph.wallParticles.pos[i]);
            instWorlds[idx] = glm::scale(instWorlds[idx], glm::vec3(rad));

            instDen[idx] = pcisph.wallParticles.density[i];
            
            instVel[idx] = glm::vec3(0.0f);

            idx++;
        }
    }
    else {
#pragma omp parallel for
        for (unsigned int i = 0; i < pcisph.numFluidParticles; i++) {
            instWorlds[idx] = glm::translate(glm::mat4(1.0f), pcisph.fluidParticles.pos[i]);
            instWorlds[idx] = glm::scale(instWorlds[idx], glm::vec3(rad));

            instDen[idx] = pcisph.fluidParticles.density[i];
            
            instVel[idx] = pcisph.fluidParticles.vel[i];

            
            //printf("%f\n", instVel[idx][0] * instVel[idx][0] + instVel[idx][1] * instVel[idx][1] + instVel[idx][2] * instVel[idx][2]);

            idx++;

        }
    }

}

