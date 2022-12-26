#include"shader.h"

#include<glad/glad.h>
#include<GLFW/glfw3.h>

#include<string>
#include<fstream>
#include<sstream>
#include <iostream>

Shader::Shader(const char* vertexPath, const char* fragPath) {
	std::string vertexCode;
	std::string fragCode;

	std::ifstream vShaderFile;
	std::ifstream fShaderFile;

	vShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	fShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	try {
		vShaderFile.open(vertexPath);
		fShaderFile.open(fragPath);

		std::stringstream vShaderStream,fShaderStream;

		vShaderStream << vShaderFile.rdbuf();
		fShaderStream << fShaderFile.rdbuf();

		vShaderFile.close();
		fShaderFile.close();

		vertexCode = vShaderStream.str();
		fragCode = fShaderStream.str();
	}
	catch(std::ifstream::failure e) {
		std::cout << "ERROR::SHADER:FILsE_NOT_SUCCESFULLY_READ" << std::endl;
	}

	const char* vShaderCode = vertexCode.c_str();
	const char* fShaderCode = fragCode.c_str();

	unsigned int vertex, fragment;

	int success;
	char infoLog[512];

	vertex = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertex, 1, &vShaderCode, NULL);
	glCompileShader(vertex);

	glGetShaderiv(vertex, GL_COMPILE_STATUS, &success);
	if (!success) {
		glGetShaderInfoLog(vertex, 512, NULL, infoLog);
		std::cout << "ERROR::SHADER::VERETEX::COMPILATION_FAILED \n" << infoLog << std::endl;
	}
	fragment = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragment, 1, &fShaderCode, NULL);
	glCompileShader(fragment);

	glGetShaderiv(fragment, GL_COMPILE_STATUS, &success);
	if (!success) {
		glGetShaderInfoLog(fragment, 512, NULL, infoLog);
		std::cout<< "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED \n" << infoLog << std::endl;
	}
	
	ID = glCreateProgram();
	glAttachShader(ID, vertex);
	glAttachShader(ID, fragment);
	glLinkProgram(ID);

	glGetProgramiv(ID, GL_LINK_STATUS, &success);
	if (!success) {
		glGetProgramInfoLog(ID, 512, NULL, infoLog);
		std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED \n" << infoLog << std::endl;
	}

	glDeleteShader(fragment);
	glDeleteShader(vertex);
}

void Shader::use() {
	glUseProgram(ID);
}