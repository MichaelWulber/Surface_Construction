#include <imebra/imebra.h>
#include <iostream>
#include <vector>
#include <string>

#define GLEW_STATIC
#include <Gl/glew.h>

#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Shader.h"
#include "LookUpTable.h"
#include "Mesh.h"

#include "SphereFunc.h"
#include "TorusFunc.h"
#include "GenusFunc.h"
#include "TrefoilFunc.h"
#include "HyperbolicParaboloidFunc.h"
#include "DoubleConeFunc.h"

std::vector<GLfloat> genMesh(ImplicitFunc* function, GLfloat cubeSize);

std::vector<GLfloat> genCTMesh(const std::vector<std::vector<std::vector<int>>> &vLuminance, const std::vector<std::vector<std::vector<int>>> &vCondition, double aspect_ratio);
GLfloat CTinterpolate(GLfloat a, GLfloat b);
std::vector<GLfloat> findVerticesCT(int i, int j, int k, int index, GLfloat* vertex[3], const std::vector<std::vector<std::vector<int>>> &vLuminance);

int edgeListIndex(const bool arr[8]);
std::vector<GLfloat> findVertices(int i, int j, int k, int index, GLfloat* vertex[3], GLfloat*** vals);
GLfloat interpolate(GLfloat a, GLfloat aVal, GLfloat b, GLfloat bVal);

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);

std::string getFilePath(int index);

const GLint WIDTH = 1200, HEIGHT = 1200;
int screenWidth, screenHeight;

// Mesh Objects
Mesh sphere;
Mesh torus;
Mesh genus;
Mesh trefoil;
Mesh hyperbolicParaboloid;
Mesh doubleCones;
Mesh ct;
Mesh current;

ImplicitFunc* sphereFunc;
ImplicitFunc* torusFunc;
ImplicitFunc* genusFunc;
ImplicitFunc* trefoilFunc;
ImplicitFunc* hyperbolicParaboloidFunc;
ImplicitFunc* doubleConesFunc;

int main() {
	glfwInit();

	// set the version of openGL to 3.3
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

	// use the modern stuff
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

	GLFWwindow *window = glfwCreateWindow(WIDTH, HEIGHT, "Marching Cubes", nullptr, nullptr);

	// ensures pixels coordinates are mapped to the screen correctly (accounts for pixel density)
	int screenWidth, screenHeight;
	glfwGetFramebufferSize(window, &screenWidth, &screenHeight);

	if (nullptr == window) {
		std::cout << "Failed to create GLFW window" << std::endl;
		return EXIT_FAILURE;
	}

	// set the current context to the window we just created
	glfwMakeContextCurrent(window);
	glfwSetKeyCallback(window, key_callback);

	// use modern approach to obtain function pointers
	glewExperimental = GL_TRUE;

	if (GLEW_OK != glewInit()) {
		std::cout << "Failed to initialize GLEW" << std::endl;
		return EXIT_FAILURE;
	}

	glViewport(0, 0, screenWidth, screenHeight);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	Shader ourShader("core.vert", "core.frag");

	// create torus mesh
	torusFunc = new TorusFunc(0.5f, 0.3f);
	std::vector<GLfloat> torusVertices = genMesh(torusFunc, 1.0f);
	torus = Mesh(0.95f, 0.52f, 0.16f);
	torus.setVPositions(torusVertices);
	torus.genVNormals();
	torus.genBuffer();

	// create sphere mesh
	sphereFunc = new SphereFunc(1.0f);
	std::vector<GLfloat> sphereVertices = genMesh(sphereFunc, 1.0f);
	sphere = Mesh(0.95f, 0.52f, 0.16f);
	sphere.setVPositions(sphereVertices);
	sphere.genVNormals();
	sphere.genBuffer();

	// create Genus 2 mesh
	genusFunc = new GenusFunc();
	std::vector<GLfloat> genusVertices = genMesh(genusFunc, 2.0f);
	genus = Mesh(0.95f, 0.52f, 0.16f);
	genus.setVPositions(genusVertices);
	genus.genVNormals();
	genus.genBuffer();

	// create trefoil mesh
	trefoilFunc = new TrefoilFunc(0.15f, 0.10f);
	std::vector<GLfloat> trefoilVertices = genMesh(trefoilFunc, 2.0f);
	trefoil = Mesh(0.95f, 0.52f, 0.16f);
	trefoil.setVPositions(trefoilVertices);
	trefoil.genVNormals();
	trefoil.genBuffer();

	// create double cones mesh
	doubleConesFunc = new DoubleConeFunc(1.0f, 1.0f);
	std::vector<GLfloat> doubleConesVertices = genMesh(doubleConesFunc, 1.0f);
	doubleCones = Mesh(0.95f, 0.52f, 0.16f);
	doubleCones.setVPositions(doubleConesVertices);
	doubleCones.genVNormals();
	doubleCones.genBuffer();

	// create hyperbolic paraboloid  mesh
	hyperbolicParaboloidFunc = new HyperbolicParaboloidFunc(1.0f, -0.5f);
	std::vector<GLfloat> hyperbolicParaboloidVertices = genMesh(hyperbolicParaboloidFunc, 1.0f);
	hyperbolicParaboloid = Mesh(0.95f, 0.52f, 0.16f);
	hyperbolicParaboloid.setVPositions(hyperbolicParaboloidVertices);
	hyperbolicParaboloid.genVNormals();
	hyperbolicParaboloid.genBuffer();

	// --- create mesh from DICOM ---
	std::cout << "reading DICOM...  ";

	const int dimX = 64;
	const int dimY = 64;
	const int numFrames = 64;
	float frame_scale = 729 / numFrames;
	int start = 0;

	std::vector<std::vector<std::vector<int>>> vLuminance(numFrames, std::vector<std::vector<int>>(dimX, std::vector<int>(dimY, 0)));
	std::vector<std::vector<std::vector<int>>> vCondition(numFrames, std::vector<std::vector<int>>(dimX, std::vector<int>(dimY, 0)));
	for (int index = 1; index <= numFrames; ++index) {
		// Load DICOM file
		std::unique_ptr<imebra::DataSet> loadedDataSet(imebra::CodecFactory::load(getFilePath(index * frame_scale + start)));

		// Get Image
		std::unique_ptr<imebra::Image> image(loadedDataSet->getImageApplyModalityTransform(0));

		// Get color space
		std::string colorSpace = image->getColorSpace();
		int width = image->getWidth();
		int height = image->getHeight();

		// Retrieve data handler
		std::unique_ptr<imebra::ReadingDataHandlerNumeric> dataHandler(image->getReadingDataHandler());

		// Obtain Threshold Frame
		float x_factor = (float)width / (float)dimX;
		float y_factor = (float)height / (float)dimY;
		for (int i = 0; i < width; ++i) {
			for (int j = 0; j < height; ++j) {
				std::int32_t luminance = dataHandler->getSignedLong(j * width + i);
				vLuminance[index - 1][i / x_factor][j / y_factor] = luminance;

				if (luminance > 0) {
					vCondition[index - 1][i / x_factor][j / y_factor] = 1;
				}
				else {
					vCondition[index - 1][i / x_factor][j / y_factor] = 0;
				}
			}
		}

		std::cout << index << std::endl;
	}
	std::cout << "finished" << std::endl;
	// --- end ---

	std::vector<GLfloat> CTVertices = genCTMesh(vLuminance, vCondition, 1.25);
	ct = Mesh(0.95f, 0.52f, 0.16f);
	ct.setVPositions(CTVertices);
	ct.genVNormals();
	ct.genBuffer();

	// set current mesh
	current = ct;

	// create openGL buffer and attribute objects
	GLuint VBO, VAO;
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	current.bindBuffer();

	// create projection transformation
	glm::mat4 projection;
	projection = glm::perspective(glm::radians(45.0f), (GLfloat)(screenWidth) / (GLfloat)(screenHeight), 0.1f, 1000.0f);

	// GAME LOOP
	while (!glfwWindowShouldClose(window)) {
		glfwPollEvents();

		glClearColor(0.0f, 0.3f, 0.8f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		ourShader.use();

		// set up MVP matrix
		glm::mat4 model(1.0f);
		model = glm::rotate(model, (GLfloat)glfwGetTime() * -0.75f, glm::vec3(0.0f, 0.3f, 1.0f));
		glm::mat4 view = glm::lookAt(
			glm::vec3(3, 3, 3), // Camera is at (3,3,3), in World Space
			glm::vec3(0, 0, 0), // and looks at the origin
			glm::vec3(0, 1, 0)  // Head is up (set to 0,-1,0 to look upside-down)
			);
		glm::mat4 MVP = projection * view * model;

		GLint MVPLoc = glGetUniformLocation(ourShader.Program, "MVP");
		glUniformMatrix4fv(MVPLoc, 1, GL_FALSE, glm::value_ptr(MVP));

		GLint modelLoc = glGetUniformLocation(ourShader.Program, "model");
		glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

		glBindVertexArray(VAO);
		glDrawArrays(GL_TRIANGLES, 0, current.vPositions.size() / 3);
		glBindVertexArray(0);

		glfwSwapBuffers(window);
	}
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);

	glfwTerminate();

	return EXIT_SUCCESS;
}

std::vector<GLfloat> genMesh(ImplicitFunc* function, GLfloat cubeSize) {
	std::cout << "generating mesh..." << std::endl;
	GLfloat minX = -cubeSize;
	GLfloat minY = -cubeSize;
	GLfloat minZ = -cubeSize;
	GLfloat maxX = cubeSize;
	GLfloat maxY = cubeSize;
	GLfloat maxZ = cubeSize;
	GLfloat x, y, z, a;
	bool byteArray[8];

	const GLint dim = 64;

	bool vertices[dim][dim][dim];


	// array of values for x, y, z
	// vertexCoord[0][] = x's, vertexCoord[1][] = y's, vertex Coord[2][] = z's
	GLfloat* vertexCoord[3] = { new GLfloat[dim], new GLfloat[dim], new GLfloat[dim] };
	for (GLint i = 0; i < dim; ++i) {
		a = ((GLfloat)i / ((GLfloat)dim - 1));
		x = maxX * a + minX * (1.0f - a);
		y = maxY * a + minY * (1.0f - a);
		z = maxZ * a + minZ * (1.0f - a);
		vertexCoord[0][i] = x;
		vertexCoord[1][i] = y;
		vertexCoord[2][i] = z;
	}

	// vertices stores 0 or 1 depending on whether vertex is inside sphere or not
	// vertexVals stores the actual value from the implicit function
	GLfloat*** vertexVals = new GLfloat**[dim];
	for (GLint i = 0; i < dim; ++i) {
		vertexVals[i] = new GLfloat*[dim];

		for (GLint j = 0; j < dim; ++j) {
			vertexVals[i][j] = new GLfloat[dim];

			for (GLint k = 0; k < dim; ++k) {
				x = vertexCoord[0][i];
				y = vertexCoord[1][j];
				z = vertexCoord[2][k];
				vertices[i][j][k] = function->isInside(x, y, z);
				vertexVals[i][j][k] = function->function(x, y, z);
			}
		}
	}


	// Go through every cube and check vertices;
	// triangleVertices stores vertices of facets as {x0, y0, z0, x1, y1, z1, ..., xn, yn, zn}
	std::vector<GLfloat> triangleVertices;
	std::vector<GLfloat> temp;
	for (GLint i = 0; i < dim - 1; ++i) {
		for (GLint j = 0; j < dim - 1; ++j) {
			for (GLint k = 0; k < dim - 1; ++k) {
				byteArray[0] = vertices[i][j][k];
				byteArray[1] = vertices[i + 1][j][k];
				byteArray[2] = vertices[i + 1][j][k + 1];
				byteArray[3] = vertices[i][j][k + 1];
				byteArray[4] = vertices[i][j + 1][k];
				byteArray[5] = vertices[i + 1][j + 1][k];
				byteArray[6] = vertices[i + 1][j + 1][k + 1];
				byteArray[7] = vertices[i][j + 1][k + 1];
				int index = edgeListIndex(byteArray);

				temp = findVertices(i, j, k, index, vertexCoord, vertexVals);
				triangleVertices.insert(triangleVertices.end(), temp.begin(), temp.end());
			}
		}
	}
	std::cout << "mesh complete" << std::endl;

	return triangleVertices;
}

int edgeListIndex(const bool arr[8]) {
	int index = 0;
	int factor = 1;
	for (int i = 0; i < 8; ++i) {
		index += arr[i] * factor;
		factor = 2 * factor;
	}
	return index;
}

std::vector<GLfloat> findVertices(int i, int j, int k, int index,
	GLfloat* vertex[3], GLfloat*** vals) {
	std::vector<GLfloat> triangleVertices;
	int edgeNum;
	GLfloat intersection;
	GLfloat aVal, bVal;
	GLfloat a, b;
	GLfloat x, y, z;



	for (int e = 0; e < 13; ++e) {
		edgeNum = aCases[index][e];
		switch (edgeNum) {
		case -1:
			return triangleVertices;
		case 0:
			y = vertex[1][j];
			z = vertex[2][k];

			a = vertex[0][i];
			aVal = vals[i][j][k];
			b = vertex[0][i + 1];
			bVal = vals[i + 1][j][k];
			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(intersection);
			triangleVertices.push_back(y);
			triangleVertices.push_back(z);
			break;
		case 1:
			x = vertex[0][i + 1];
			y = vertex[1][j];

			a = vertex[2][k];
			aVal = vals[i + 1][j][k];
			b = vertex[2][k + 1];
			bVal = vals[i + 1][j][k + 1];
			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(x);
			triangleVertices.push_back(y);
			triangleVertices.push_back(intersection);
			break;
		case 2:
			y = vertex[1][j];
			z = vertex[2][k + 1];

			a = vertex[0][i];
			aVal = vals[i][j][k + 1];
			b = vertex[0][i + 1];
			bVal = vals[i + 1][j][k + 1];
			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(intersection);
			triangleVertices.push_back(y);
			triangleVertices.push_back(z);
			break;
		case 3:
			x = vertex[0][i];
			y = vertex[1][j];

			a = vertex[2][k];
			aVal = vals[i][j][k];
			b = vertex[2][k + 1];
			bVal = vals[i][j][k + 1];
			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(x);
			triangleVertices.push_back(y);
			triangleVertices.push_back(intersection);
			break;
		case 4:
			y = vertex[1][j + 1];
			z = vertex[2][k];

			a = vertex[0][i];
			aVal = vals[i][j + 1][k];
			b = vertex[0][i + 1];
			bVal = vals[i + 1][j + 1][k];
			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(intersection);
			triangleVertices.push_back(y);
			triangleVertices.push_back(z);
			break;
		case 5:
			x = vertex[0][i + 1];
			y = vertex[1][j + 1];

			a = vertex[2][k];
			aVal = vals[i + 1][j + 1][k];
			b = vertex[2][k + 1];
			bVal = vals[i + 1][j + 1][k + 1];
			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(x);
			triangleVertices.push_back(y);
			triangleVertices.push_back(intersection);
			break;
		case 6:
			y = vertex[1][j + 1];
			z = vertex[2][k + 1];

			a = vertex[0][i];
			aVal = vals[i][j + 1][k + 1];
			b = vertex[0][i + 1];
			bVal = vals[i + 1][j + 1][k + 1];
			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(intersection);
			triangleVertices.push_back(y);
			triangleVertices.push_back(z);
			break;
		case 7:
			x = vertex[0][i];
			y = vertex[1][j + 1];

			a = vertex[2][k];
			aVal = vals[i][j + 1][k];
			b = vertex[2][k + 1];
			bVal = vals[i][j + 1][k + 1];
			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(x);
			triangleVertices.push_back(y);
			triangleVertices.push_back(intersection);
			break;
		case 8:
			x = vertex[0][i];
			z = vertex[2][k];

			a = vertex[1][j];
			aVal = vals[i][j][k];
			b = vertex[1][j + 1];
			bVal = vals[i][j + 1][k];
			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(x);
			triangleVertices.push_back(intersection);
			triangleVertices.push_back(z);
			break;
		case 9:
			x = vertex[0][i + 1];
			z = vertex[2][k];

			a = vertex[1][j];
			aVal = vals[i + 1][j][k];
			b = vertex[1][j + 1];
			bVal = vals[i + 1][j + 1][k];
			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(x);
			triangleVertices.push_back(intersection);
			triangleVertices.push_back(z);
			break;
		case 10:
			x = vertex[0][i + 1];
			z = vertex[2][k + 1];

			a = vertex[1][j];
			aVal = vals[i + 1][j][k + 1];
			b = vertex[1][j + 1];
			bVal = vals[i + 1][j + 1][k + 1];
			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(x);
			triangleVertices.push_back(intersection);
			triangleVertices.push_back(z);
			break;
		case 11:
			x = vertex[0][i];
			z = vertex[2][k + 1];

			a = vertex[1][j];
			aVal = vals[i][j][k + 1];
			b = vertex[1][j + 1];
			bVal = vals[i][j + 1][k + 1];
			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(x);
			triangleVertices.push_back(intersection);
			triangleVertices.push_back(z);
			break;
		}
	}

	return triangleVertices;
}

GLfloat interpolate(GLfloat a, GLfloat aVal, GLfloat b, GLfloat bVal) {
	return a + ((-aVal) * (b - a) / (bVal - aVal));
}

// handle input
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_1 && action == GLFW_PRESS) {
		current = sphere;
		current.bindBuffer();
	}

	else if (key == GLFW_KEY_2 && action == GLFW_PRESS) {
		current = torus;
		current.bindBuffer();
	}

	else if (key == GLFW_KEY_3 && action == GLFW_PRESS) {
		current = genus;
		current.bindBuffer();
	}

	else if (key == GLFW_KEY_4 && action == GLFW_PRESS) {
		current = trefoil;
		current.bindBuffer();
	}

	else if (key == GLFW_KEY_5 && action == GLFW_PRESS) {
		current = hyperbolicParaboloid;
		current.bindBuffer();
	}

	else if (key == GLFW_KEY_6 && action == GLFW_PRESS) {
		current = doubleCones;
		current.bindBuffer();
	}

	else if (key == GLFW_KEY_7 && action == GLFW_PRESS) {
		current = ct;
		current.bindBuffer();
	}
}

std::string getFilePath(int index) {
	std::string path = ".\\alligator_skull\\WitmerLab_VisInt-Alligator_OUVC10606_head_";
	std::string type = ".dcm";

	if (index < 10) {
		return path + "00" + std::to_string(index) + type;
	}

	else if (index < 100) {
		return path + "0" + std::to_string(index) + type;
	}

	else {
		return path + std::to_string(index) + type;
	}
}

std::vector<GLfloat> genCTMesh(const std::vector<std::vector<std::vector<int>>> &vLuminance, const std::vector<std::vector<std::vector<int>>> &vCondition, double aspect_ratio) {
	std::cout << "generating mesh..." << std::endl;
	GLfloat minX = -1.0f * aspect_ratio;
	GLfloat minY = -1.0f;
	GLfloat minZ = -1.0f;
	GLfloat maxX = 1.0f * aspect_ratio;
	GLfloat maxY = 1.0f;
	GLfloat maxZ = 1.0f;
	GLfloat x, y, z, a;
	bool byteArray[8];

	const GLint dimX = vCondition.size();
	const GLint dimY = vCondition[0].size();
	const GLint dimZ = vCondition[0][0].size();

	GLfloat* vertexCoord[3] = { new GLfloat[dimX], new GLfloat[dimY], new GLfloat[dimZ] };
	for (GLint i = 0; i < dimX; ++i) {
		a = ((GLfloat)i / ((GLfloat)dimX - 1));
		x = maxX * a + minX * (1.0f - a);
		vertexCoord[0][i] = x;
	}

	for (GLint i = 0; i < dimY; ++i) {
		a = ((GLfloat)i / ((GLfloat)dimX - 1));
		y = maxY * a + minY * (1.0f - a);
		vertexCoord[1][i] = y;
	}

	for (GLint i = 0; i < dimZ; ++i) {
		a = ((GLfloat)i / ((GLfloat)dimX - 1));
		z = maxZ * a + minZ * (1.0f - a);
		vertexCoord[2][i] = z;
	}

	//for (GLint i = 0; i < dimX - 1; ++i) {
	//	for (GLint j = 0; j < dimY - 1; ++j) {
	//		if (vCondition[i][j][dimZ / 2] == 1 && vLuminance[i][j][dimZ / 2] < 100) {
	//			std::cout << "( " << vCondition[i][j][dimZ / 2] << ", " << vLuminance[i][j][dimZ / 2] << " )"<< std::endl;
	//		}

	//		else if (vCondition[i][j][dimZ / 2] == 0 && vLuminance[i][j][dimZ / 2] > 100) {
	//			std::cout << "( " << vCondition[i][j][dimZ / 2] << ", " << vLuminance[i][j][dimZ / 2] << " )" << std::endl;
	//		}
	//	}
	//}

	// Go through every cube and check vertices;
	// triangleVertices stores vertices of facets as {x0, y0, z0, x1, y1, z1, ..., xn, yn, zn}
	std::vector<GLfloat> triangleVertices;
	std::vector<GLfloat> temp;
	for (GLint i = 0; i < dimX - 1; ++i) {
		for (GLint j = 0; j < dimY - 1; ++j) {
			for (GLint k = 0; k < dimZ - 1; ++k) {
				byteArray[0] = vCondition[i][j][k];
				byteArray[1] = vCondition[i + 1][j][k];
				byteArray[2] = vCondition[i + 1][j][k + 1];
				byteArray[3] = vCondition[i][j][k + 1];
				byteArray[4] = vCondition[i][j + 1][k];
				byteArray[5] = vCondition[i + 1][j + 1][k];
				byteArray[6] = vCondition[i + 1][j + 1][k + 1];
				byteArray[7] = vCondition[i][j + 1][k + 1];
				int index = edgeListIndex(byteArray);

				temp = findVerticesCT(i, j, k, index, vertexCoord, vLuminance);
				triangleVertices.insert(triangleVertices.end(), temp.begin(), temp.end());
			}
		}
	}
	std::cout << "mesh complete" << std::endl;

	return triangleVertices;
}

GLfloat CTinterpolateA(GLfloat a, GLfloat aVal, GLfloat b, GLfloat bVal) {
	return a + ((0 - aVal) * (b - a) / (bVal - aVal));
}

std::vector<GLfloat> findVerticesCT(int i, int j, int k, int index, GLfloat* vertex[3], const std::vector<std::vector<std::vector<int>>> &vLuminance) {
	std::vector<GLfloat> triangleVertices;
	int edgeNum;
	GLfloat intersection;
	GLfloat aVal, bVal;
	GLfloat a, b;
	GLfloat x, y, z;



	for (int e = 0; e < 13; ++e) {
		edgeNum = aCases[index][e];
		switch (edgeNum) {
		case -1:
			return triangleVertices;
		case 0:
			y = vertex[1][j];
			z = vertex[2][k];

			a = vertex[0][i];
			aVal = vLuminance[i][j][k];
			b = vertex[0][i + 1];
			bVal = vLuminance[i + 1][j][k];

			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(intersection);
			triangleVertices.push_back(y);
			triangleVertices.push_back(z);
			break;
		case 1:
			x = vertex[0][i + 1];
			y = vertex[1][j];

			a = vertex[2][k];
			aVal = vLuminance[i + 1][j][k];
			b = vertex[2][k + 1];
			bVal = vLuminance[i + 1][j][k + 1];

			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(x);
			triangleVertices.push_back(y);
			triangleVertices.push_back(intersection);
			break;
		case 2:
			y = vertex[1][j];
			z = vertex[2][k + 1];

			a = vertex[0][i];
			aVal = vLuminance[i][j][k + 1];
			b = vertex[0][i + 1];
			bVal = vLuminance[i + 1][j][k + 1];

			intersection = CTinterpolateA(a, aVal, b, bVal);

			triangleVertices.push_back(intersection);
			triangleVertices.push_back(y);
			triangleVertices.push_back(z);
			break;
		case 3:
			x = vertex[0][i];
			y = vertex[1][j];

			a = vertex[2][k];
			aVal = vLuminance[i][j][k];
			b = vertex[2][k + 1];
			bVal = vLuminance[i][j][k + 1];

			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(x);
			triangleVertices.push_back(y);
			triangleVertices.push_back(intersection);
			break;
		case 4:
			y = vertex[1][j + 1];
			z = vertex[2][k];

			a = vertex[0][i];
			aVal = vLuminance[i][j + 1][k];
			b = vertex[0][i + 1];
			bVal = vLuminance[i + 1][j + 1][k];

			intersection = CTinterpolateA(a, aVal, b, bVal);

			triangleVertices.push_back(intersection);
			triangleVertices.push_back(y);
			triangleVertices.push_back(z);
			break;
		case 5:
			x = vertex[0][i + 1];
			y = vertex[1][j + 1];

			a = vertex[2][k];
			aVal = vLuminance[i + 1][j + 1][k];
			b = vertex[2][k + 1];
			bVal = vLuminance[i + 1][j + 1][k + 1];

			intersection = CTinterpolateA(a, aVal, b, bVal);

			triangleVertices.push_back(x);
			triangleVertices.push_back(y);
			triangleVertices.push_back(intersection);
			break;
		case 6:
			y = vertex[1][j + 1];
			z = vertex[2][k + 1];

			a = vertex[0][i];
			aVal = vLuminance[i][j + 1][k + 1];
			b = vertex[0][i + 1];
			bVal = vLuminance[i + 1][j + 1][k + 1];

			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(intersection);
			triangleVertices.push_back(y);
			triangleVertices.push_back(z);
			break;
		case 7:
			x = vertex[0][i];
			y = vertex[1][j + 1];

			a = vertex[2][k];
			aVal = vLuminance[i][j + 1][k];
			b = vertex[2][k + 1];
			bVal = vLuminance[i][j + 1][k + 1];

			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(x);
			triangleVertices.push_back(y);
			triangleVertices.push_back(intersection);
			break;
		case 8:
			x = vertex[0][i];
			z = vertex[2][k];

			a = vertex[1][j];
			aVal = vLuminance[i][j][k];
			b = vertex[1][j + 1];
			bVal = vLuminance[i][j + 1][k];

			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(x);
			triangleVertices.push_back(intersection);
			triangleVertices.push_back(z);
			break;
		case 9:
			x = vertex[0][i + 1];
			z = vertex[2][k];

			a = vertex[1][j];
			aVal = vLuminance[i + 1][j][k];
			b = vertex[1][j + 1];
			bVal = vLuminance[i + 1][j + 1][k];

			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(x);
			triangleVertices.push_back(intersection);
			triangleVertices.push_back(z);
			break;
		case 10:
			x = vertex[0][i + 1];
			z = vertex[2][k + 1];

			a = vertex[1][j];
			aVal = vLuminance[i + 1][j][k + 1];
			b = vertex[1][j + 1];
			bVal = vLuminance[i + 1][j + 1][k + 1];

			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(x);
			triangleVertices.push_back(intersection);
			triangleVertices.push_back(z);
			break;
		case 11:
			x = vertex[0][i];
			z = vertex[2][k + 1];

			a = vertex[1][j];
			aVal = vLuminance[i][j][k + 1];
			b = vertex[1][j + 1];
			bVal = vLuminance[i][j + 1][k + 1];

			intersection = interpolate(a, aVal, b, bVal);

			triangleVertices.push_back(x);
			triangleVertices.push_back(intersection);
			triangleVertices.push_back(z);
			break;
		}
	}

	return triangleVertices;
}