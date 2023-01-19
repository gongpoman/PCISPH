
#include<glad/glad.h>
#include<GLFW/glfw3.h>


#include"myutil.h"
#include"camera.h"

extern Camera cam;

extern float lastX=0, lastY=0;
extern bool isFirstMove;

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

int mortonEncode(int x, int y, int z) {
    x = (x <0) ?  -1:x;
    y = (y <0) ?  -1:y;
    z = (z <0) ? -1:z;

    if ((x < 0) || (y < 0) || (z < 0)) // out of bounary
        return -1;


    return  mortonEncode_LUT(x, y, z);
}

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
    if (nGridDivisoin <= 256)
        digit = 256;
    if (nGridDivisoin <= 128)
        digit = 128;
    if (nGridDivisoin <= 64)
        digit = 64;
    if (nGridDivisoin <= 32)
        digit = 32;
    if (nGridDivisoin <= 16)
        digit = 16;
    if (nGridDivisoin <= 8)
        digit = 8;
    if (nGridDivisoin <= 4)
        digit = 4;
    if (nGridDivisoin <= 2)
        digit = 2;

    return digit;
}


uint32_t Compact1By2(uint32_t x)
{
    x &= 0x09249249;                  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
    x = (x ^ (x >> 2)) & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x ^ (x >> 4)) & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x ^ (x >> 8)) & 0xff0000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x ^ (x >> 16)) & 0x000003ff; // x = ---- ---- ---- ---- ---- --98 7654 3210
    return x;
}

uint32_t DecodeMorton3X(uint32_t code)
{
    return Compact1By2(code >> 0);
}

uint32_t DecodeMorton3Y(uint32_t code)
{
    return Compact1By2(code >> 1);
}

uint32_t DecodeMorton3Z(uint32_t code)
{
    return Compact1By2(code >> 2);
}

glm::ivec3 mortonDecode(unsigned int morton) {
    return glm::ivec3(DecodeMorton3X(morton), DecodeMorton3Y(morton), DecodeMorton3Z(morton));
}


