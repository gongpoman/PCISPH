#pragma once

#include<glm/glm.hpp>


class Particle {
public:
	unsigned int cellIdx;

	glm::vec3 vel;
	glm::vec3 pos;
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

		vel = glm::vec3(0);
		pos = glm::vec3(0);

		pv = glm::vec3(0);
		px = glm::vec3(0);

		force_p = glm::vec3(0);
		force_vis = glm::vec3(0);
		force_ext = glm::vec3(0);

		pressure = 0.0f;
		density = 0.0f;
		densityVar=0.0f;
	}

};