#include <imebra/imebra.h>
#include <iostream>
#include <vector>

#define GLEW_STATIC
#include <Gl/glew.h>

#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Shader.h"
#include "LookUpTable.h"
#include "Mesh.h"

bool isInsideSphere(GLfloat x, GLfloat y, GLfloat z, GLfloat radius);
GLfloat sphereImplicitFunction(GLfloat x, GLfloat y, GLfloat z);
std::vector<GLfloat> genSphereMesh();

bool isInsideTorus(GLfloat x, GLfloat y, GLfloat z, GLfloat R, GLfloat a);
GLfloat torusImplicitFunction(GLfloat x, GLfloat y, GLfloat z, GLfloat R, GLfloat a);
std::vector<GLfloat> genTorusMesh();

int edgeListIndex(const bool arr[8]);
std::vector<GLfloat> findVertices(int i, int j, int k, int index, GLfloat* vertex[3], GLfloat*** vals);
GLfloat interpolate(GLfloat a, GLfloat aVal, GLfloat b, GLfloat bVal);

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);

const GLint WIDTH = 1200, HEIGHT = 1200;
int screenWidth, screenHeight;

// Mesh Objects
Mesh sphere;
Mesh torus;
Mesh current;

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
	std::vector<GLfloat> torusVertices = genTorusMesh();
	torus = Mesh(0.1f, 0.3f, 0.8f, 0.9f, 0.9f, 1.0f);
	torus.setVPositions(torusVertices);
	torus.genBuffer();

	// create sphere mesh
	std::vector<GLfloat> sphereVertices = genSphereMesh();
	sphere = Mesh(0.1f, 0.3f, 0.8f, 0.9f, 0.9f, 1.0f);
	sphere.setVPositions(sphereVertices);
	sphere.genBuffer();

	// --- create mesh from DICOM ---

	// FOR EACH FILE:
	// Load DICOM file
	std::unique_ptr<imebra::DataSet> loadedDataSet(imebra::CodecFactory::load("C:\\Users\\wulbe_000\\Desktop\\alligator_head\\test.dcm"));

	// Get Image
	std::unique_ptr<imebra::Image> image(loadedDataSet->getImageApplyModalityTransform(0));

	// Get dicom data
	double slice_thickness = loadedDataSet->getDouble(imebra::TagId(imebra::tagId_t::SliceThickness_0018_0050), 0, -1);

	// Get the color space
	std::string colorSpace = image->getColorSpace();
	int width = image->getWidth();
	int height = image->getHeight();

	// Retrieve the data handler
	std::unique_ptr<imebra::ReadingDataHandlerNumeric> dataHandler(image->getReadingDataHandler());

	std::cout << "width: " << width << std::endl;
	std::cout << "height: " << width << std::endl;

	//// Obtain Threshold Frame
	std::vector<std::vector<int> > frame(width, std::vector<int>(height, 0));
	for (int i = 0; i < width; ++i) {
		for (int j = 0; j < height; ++j) {
			// For monochrome images
			std::int32_t luminance = dataHandler->getSignedLong(j * width + i);

			if (luminance > 100 && luminance <= 300) {
				frame[i][j] = 1;
			}
		}
	}

	// print frame values
	for (int i = width/2; i < width/2 + width/10; ++i) {
		for (int j = height/2; j < height/2 + height/10; ++j) {
			std::cout << frame[i][j] << " ";
		}
		std::cout << std::endl;
	}

	// --- end ---

	// set current mesh
	current = torus;

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

		glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		ourShader.use();

		// set up MVP matrix
		glm::mat4 model(1.0f);
		model = glm::rotate(model, (GLfloat)glfwGetTime() * -1.0f, glm::vec3(1.0f, 1.0f, 0.0f));
		glm::mat4 view = glm::lookAt(
			glm::vec3(3, 3, 3), // Camera is at (3,3,3), in World Space
			glm::vec3(0, 0, 0), // and looks at the origin
			glm::vec3(0, 1, 0)  // Head is up (set to 0,-1,0 to look upside-down)
			);
		glm::mat4 MVP = projection * view * model;

		GLint MVPLoc = glGetUniformLocation(ourShader.Program, "MVP");
		glUniformMatrix4fv(MVPLoc, 1, GL_FALSE, glm::value_ptr(MVP));

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

GLfloat isovalue = 0.5f;		// In this case, isovalue = radius;


bool isInsideSphere(GLfloat x, GLfloat y, GLfloat z, GLfloat radius) {
	return (x * x + y * y + z * z <= radius);
}

GLfloat sphereImplicitFunction(GLfloat x, GLfloat y, GLfloat z) {
	return (x * x + y * y + z * z);
}

std::vector<GLfloat> genSphereMesh() {
	std::cout << "generating mesh..." << std::endl;
	GLfloat minX = -1.0f;
	GLfloat minY = -1.0f;
	GLfloat minZ = -1.0f;
	GLfloat maxX = 1.0f;
	GLfloat maxY = 1.0f;
	GLfloat maxZ = 1.0f;
	GLfloat x, y, z, a;
	bool byteArray[8];

	isovalue = 0.5f;

	const GLint dim = 5; // number of vertices on bounding box edge
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
				vertices[i][j][k] = isInsideSphere(x, y, z, isovalue);
				vertexVals[i][j][k] = sphereImplicitFunction(x, y, z);
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
	return a + ((isovalue - aVal) * (b - a) / (bVal - aVal));
}

bool isInsideTorus(GLfloat x, GLfloat y, GLfloat z, GLfloat R, GLfloat a) {
	return ((x*x + y*y + z*z + R*R - a*a) * (x*x + y*y + z*z + R*R - a*a) - 4 * R*R*(x*x + y*y) < 0);
}

GLfloat torusImplicitFunction(GLfloat x, GLfloat y, GLfloat z, GLfloat R, GLfloat a) {
	return (x*x + y*y + z*z + R*R - a*a) * (x*x + y*y + z*z + R*R - a*a) - 4 * R*R*(x*x + y*y);
}

std::vector<GLfloat> genTorusMesh() {
	std::cout << "generating mesh..." << std::endl;
	GLfloat minX = -1.0f;
	GLfloat minY = -1.0f;
	GLfloat minZ = -1.0f;
	GLfloat maxX = 1.0f;
	GLfloat maxY = 1.0f;
	GLfloat maxZ = 1.0f;
	GLfloat x, y, z, a;
	bool byteArray[8];

	const GLfloat r1 = 0.5f;
	const GLfloat r2 = 0.1f;

	const GLint dim = 50; // number of vertices on bounding box edge
	bool vertices[dim][dim][dim];

	isovalue = 0.0f;

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
				vertices[i][j][k] = isInsideTorus(x, y, z, r1, r2);
				vertexVals[i][j][k] = torusImplicitFunction(x, y, z, r1, r2);
			}
		}
	}

	//for (int i = 0; i < dim; ++i) {
	//	for (int j = 0; j < dim; ++j) {
	//		std::cout << vertices[i][j][dim / 2] << " ";
	//	}
	//	std::cout << std::endl;
	//}

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

// handle input
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_S && action == GLFW_PRESS) {
		current = sphere;
		current.bindBuffer();
	}

	else if (key == GLFW_KEY_T && action == GLFW_PRESS) {
		current = torus;
		current.bindBuffer();
	}
}