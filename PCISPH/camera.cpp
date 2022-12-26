
#include<glad/glad.h>
#include<GLFW/glfw3.h>

#include<glm/glm.hpp>
#include<glm/gtc/matrix_transform.hpp>

#include"camera.h"


extern const float deltaTime;

glm::mat4 Camera::getViewMatrix() {
	return glm::lookAt(position, position + front, up);
}
void Camera::processKeyboard(Camera_Movement dir) {
	float velocity = movementSpeed * deltaTime;
	if (dir == FORWARD) {
		position += velocity * front;
	}
	if (dir == BACKWARD) {
		position -= velocity * front;
	}
	if (dir == LEFT) {
		position -= velocity * right;
	}
	if (dir == RIGHT) {
		position += velocity * right;
	}
}
void Camera::processMovement(float xoffset, float yoffset, GLboolean constraintPitch) {
	xoffset *= mouseSensitivity;
	yoffset *= mouseSensitivity;

	yaw += xoffset;
	pitch -= yoffset;

	if (constraintPitch) {
		if (pitch> 89.0f)
			pitch = 89.0f;
		if (pitch < -89.0f)
			pitch = -89.0f;
	}
	updateCameraVectors();
}
void Camera::processMouseScroll(float yoffset) {
	zoom -= yoffset;
	if (zoom < 1.0f) {
		zoom = 1.0f;
	}
	if (zoom > 45.0f)
		zoom = 45.0f;
}

void Camera::updateCameraVectors() {
	using namespace glm;
	vec3 front;
	front.x = cos(radians(pitch)) * cos(radians(yaw));
	front.y = sin(radians(pitch));
	front.z = cos(radians(pitch)) * sin(radians(yaw));
	this->front=front;
	this->right = normalize(cross(worldUp,-front));
	this->up = normalize(cross(right, front));
}