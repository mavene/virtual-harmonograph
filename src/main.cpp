#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <list>
#include <map>
#include <numeric>

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>

#include "shaderSource.h"
#include <imgui_impl_opengl3.h>
#include <imgui_impl_glfw.h>

using namespace std;

#define _ROTATE_FACTOR 0.005f
#define _SCALE_FACTOR 0.005f
#define _TRANS_FACTOR 0.003f
#define _Z_NEAR 0.001f
#define _Z_FAR 100.0f
#define M_PI 3.14159265358979323846

#define NUMBER_OF_VERTICES 2000 // Define the number of vertices

/***********************************************************************/
/**************************   global variables   ***********************/
/***********************************************************************/

// declaration
void processInput(GLFWwindow *window);
void framebuffer_size_callback(GLFWwindow *window, int width, int height);
void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods);
void mouse_button_callback(GLFWwindow *window, int button, int action, int mods);
void scroll_callback(GLFWwindow *window, double xoffset, double yoffset);
void cursor_pos_callback(GLFWwindow *window, double xpos, double ypos);

// Window size
unsigned int winWidth = 1600;
unsigned int winHeight = 2400;

// Camera
glm::vec3 camera_position = glm::vec3(0.0f, 0.0f, 2.5f);
glm::vec3 camera_target = glm::vec3(0.0f, 0.0f, 0.0f);
glm::vec3 camera_up = glm::vec3(0.0f, 1.0f, 0.0f);
float camera_fovy = 45.0f;
glm::mat4 projection;

// Mouse interaction
bool leftMouseButtonHold = false;
bool isFirstMouse = true;
float prevMouseX;
float prevMouseY;
glm::mat4 modelMatrix = glm::mat4(1.0f);

// Colour
int lightColorID = 0;
int meshColorID1 = 0;
int meshColorID2 = 0;
// Note: Colors have to be readjusted to achieve the gradients in the palette (not exact RGB)

// Light color table
glm::vec3 lightColorTable[3] =
    {
        glm::vec3(0.8745, 0.749, 0.549),
        glm::vec3(0.749, 0.721, 0.596),
        glm::vec3(0.796, 0.678, 0.581)};

// Mesh color table
glm::vec3 meshColorTable1[3] =
    {
        glm::vec3(0.87, 0.7, 0.49),
        glm::vec3(0.368, 0.357, 0.349),
        glm::vec3(0.76, 0.678, 0.581)};

glm::vec3 meshColorTable2[3] =
    {
        glm::vec3(0.511, 0.364, 0.43),
        glm::vec3(0.267, 0.556, 0.552),
        glm::vec3(0.532, 0.38, 0.398)};

// Animation Control
bool isAnimating = true;
bool isExported = false;
float animationTime = 0.0f;
int reconstructMethod = -1; // old surface reconstruction is 1

// Parameters
float amplitude = 0.5f;
float animationIncrement = 0.01f;
float *freqPtr1 = new float[3];
float *freqPtr2 = new float[3];
float *dampPtr1 = new float[3];
float *dampPtr2 = new float[3];
float *phasePtr1 = new float[3];
float *phasePtr2 = new float[3];

// Data buffers
std::vector<float> vertices;

///=========================================================================================///
///                            Functions for Manipulating 3D Model
///=========================================================================================///

void RotateModel(float angle, glm::vec3 axis)
{
    glm::vec3 rotateCenter = glm::vec3(modelMatrix[3][0], modelMatrix[3][1], modelMatrix[3][2]);

    glm::mat4 rotateMatrix = glm::mat4(1.0f);
    rotateMatrix = glm::translate(rotateMatrix, rotateCenter);
    rotateMatrix = glm::rotate(rotateMatrix, angle, axis);
    rotateMatrix = glm::translate(rotateMatrix, -rotateCenter);

    modelMatrix = rotateMatrix * modelMatrix;
}

void TranslateModel(glm::vec3 transVec)
{
    glm::mat4 translateMatrix = glm::mat4(1.0f);
    translateMatrix = glm::translate(translateMatrix, transVec);

    modelMatrix = translateMatrix * modelMatrix;
}

void ScaleModel(float scale)
{
    glm::vec3 scaleCenter = glm::vec3(modelMatrix[3][0], modelMatrix[3][1], modelMatrix[3][2]);

    glm::mat4 scaleMatrix = glm::mat4(1.0f);
    scaleMatrix = glm::translate(scaleMatrix, scaleCenter);
    scaleMatrix = glm::scale(scaleMatrix, glm::vec3(scale, scale, scale));
    scaleMatrix = glm::translate(scaleMatrix, -scaleCenter);

    modelMatrix = scaleMatrix * modelMatrix;
}

void SetMeshColor()
{
    // Increment color IDs (cycle between 0-2)
    lightColorID = (lightColorID + 1) % 3;
    meshColorID1 = (meshColorID1 + 1) % 3;
    meshColorID2 = (meshColorID2 + 1) % 3;
}

void SetPresets(int id)
{
    switch (id)
    {
    // Preset 1 - butterfly shaped
    case 0:
        freqPtr1[0] = 2.0f; // X component of 1st pendulum
        freqPtr1[1] = 2.0f; // Y component
        freqPtr1[2] = 2.1f; // Z component

        freqPtr2[0] = 3.0f; // X component of 2nd pendulum
        freqPtr2[1] = 3.0f; // Y component
        freqPtr2[2] = 3.1f; // Z component

        dampPtr1[0] = 0.002f;
        dampPtr1[1] = 0.002f;
        dampPtr1[2] = 0.002f;

        dampPtr2[0] = 0.002f;
        dampPtr2[1] = 0.002f;
        dampPtr2[2] = 0.002f;

        phasePtr1[0] = 0;
        phasePtr1[1] = M_PI / 4;
        phasePtr1[2] = 0;

        phasePtr2[0] = M_PI / 2;
        phasePtr2[1] = M_PI / 4;
        phasePtr2[2] = M_PI / 2;
        break;
    // Preset 2 - croissant shaped
    case 1:
        freqPtr1[0] = 1.618f;
        freqPtr1[1] = 2.618f;
        freqPtr1[2] = 4.236f;

        freqPtr2[0] = 1.618f;
        freqPtr2[1] = 2.618f;
        freqPtr2[2] = 4.236f;

        dampPtr1[0] = 0.008f;
        dampPtr1[1] = 0.008f;
        dampPtr1[2] = 0.008f;

        dampPtr2[0] = 0.008f;
        dampPtr2[1] = 0.008f;
        dampPtr2[2] = 0.008f;

        phasePtr1[0] = 0;
        phasePtr1[1] = M_PI / 2;
        phasePtr1[2] = M_PI;

        phasePtr2[0] = M_PI / 4;
        phasePtr2[1] = (3 * M_PI) / 4;
        phasePtr2[2] = (5 * M_PI) / 4;
        break;
    // Preset 3 - jet shaped
    case 2:
        freqPtr1[0] = 1.3502f;
        freqPtr1[1] = 2.6835f;
        freqPtr1[2] = 2.7173f;

        freqPtr2[0] = 1.3502f;
        freqPtr2[1] = 2.6835f;
        freqPtr2[2] = 2.7173f;

        dampPtr1[0] = 0.0f;
        dampPtr1[1] = 0.0295f;
        dampPtr1[2] = 0.0169f;

        dampPtr2[0] = 0.0f;
        dampPtr2[1] = 0.0295f;
        dampPtr2[2] = 0.00169f;

        phasePtr1[0] = 1.3324f;
        phasePtr1[1] = 1.6826f;
        phasePtr1[2] = 1.6031f;

        phasePtr2[0] = 1.3324f;
        phasePtr2[1] = 1.6826f;
        phasePtr2[2] = 1.6031f;
        break;
    default:
        break;
    }
}

///=========================================================================================///
///                                    Callback Functions
///=========================================================================================///

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow *window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and
    // height will be significantly larger than specified on retina displays.

    glViewport(0, 0, width, height);

    winWidth = width;
    winHeight = height;
}
// glfw: whenever the mouse button is clicked, this callback is called
// ---------------------------------------------------------
void mouse_button_callback(GLFWwindow *window, int button, int action, int mods)
{
    auto &io = ImGui::GetIO();
    if (io.WantCaptureMouse || io.WantCaptureKeyboard)
    {
        return;
    }
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        leftMouseButtonHold = true;
    }
    else
    {
        leftMouseButtonHold = false;
    }
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow *window, double xoffset, double yOffset)
{
    float scale = 1.0f + _SCALE_FACTOR * yOffset;

    ScaleModel(scale);
}

// glfw: whenever the cursor moves, this callback is called
// ---------------------------------------------------------
void cursor_pos_callback(GLFWwindow *window, double mouseX, double mouseY)
{
    float dx, dy;
    float nx, ny, scale, angle;

    if (leftMouseButtonHold)
    {
        if (isFirstMouse)
        {
            prevMouseX = mouseX;
            prevMouseY = mouseY;
            isFirstMouse = false;
        }

        else
        {
            if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
            {
                float dx = _TRANS_FACTOR * (mouseX - prevMouseX);
                float dy = -1.0f * _TRANS_FACTOR * (mouseY - prevMouseY); // reversed since y-coordinates go from bottom to top

                prevMouseX = mouseX;
                prevMouseY = mouseY;

                TranslateModel(glm::vec3(dx, dy, 0));
            }

            else
            {
                float dx = mouseX - prevMouseX;
                float dy = -(mouseY - prevMouseY); // reversed since y-coordinates go from bottom to top

                prevMouseX = mouseX;
                prevMouseY = mouseY;

                // Rotation
                nx = -dy;
                ny = dx;
                scale = sqrt(nx * nx + ny * ny);

                // We use "ArcBall Rotation" to compute the rotation axis and angle based on the mouse motion
                nx = nx / scale;
                ny = ny / scale;
                angle = scale * _ROTATE_FACTOR;

                RotateModel(angle, glm::vec3(nx, ny, 0.0f));
            }
        }
    }

    else
    {
        isFirstMouse = true;
    }
}

// glfw: whenever a key is pressed, this callback is called
// ----------------------------------------------------------------------
void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_C && action == GLFW_PRESS)
    {
        SetMeshColor();
    }
}

///=========================================================================================///
///                          Vertex Normals + Surfaces + Extrusions
///=========================================================================================///

// Function to calculate proper normals (from tangents) for each vertex using neighbours
std::vector<glm::vec3> calculateNormals(const std::vector<float> &vertices)
{
    std::vector<glm::vec3> normals;
    // Calculate tangent vectors
    for (size_t i = 0; i < vertices.size(); i += 3)
    {
        // Calculate the direction from the current point to the next point
        glm::vec3 currentPoint(vertices[i], vertices[i + 1], vertices[i + 2]);
        glm::vec3 nextPoint;
        if (i + 3 < vertices.size())
        {
            nextPoint = glm::vec3(vertices[i + 3], vertices[i + 4], vertices[i + 5]);
        }
        else
        {
            // If it's the last point, use the previous point as the next point
            nextPoint = glm::vec3(vertices[i - 3], vertices[i - 2], vertices[i - 1]);
        }
        glm::vec3 tangent = glm::normalize(nextPoint - currentPoint);
        normals.push_back(tangent); // For now, it looks nicer with the directional so its not strictly normals

        // Calculate the normal (90deg rotation)
        // glm::vec3 normal = glm::vec3(-tangent.y, tangent.x, 0.0f);
        // normals.push_back(normal);
    }
    return normals;
}

void calculateTriangleStripNormals(const std::vector<glm::vec3> &vertices, const std::vector<unsigned int> &indices, std::vector<glm::vec3> &normals, bool invertNormals = false)
{
    normals.clear();
    normals.resize(vertices.size(), glm::vec3(0.0f));

    // Calculate normals for each triangle in the triangle strip
    for (size_t i = 0; i < indices.size() - 2; ++i)
    {
        const glm::vec3 &v0 = vertices[indices[i]];
        const glm::vec3 &v1 = vertices[indices[i + 1]];
        const glm::vec3 &v2 = vertices[indices[i + 2]];

        glm::vec3 edge1 = v1 - v0;
        glm::vec3 edge2 = v2 - v0;
        glm::vec3 triangleNormal = glm::cross(edge1, edge2);

        // Ensure correct winding order by flipping normals if necessary
        if (glm::dot(triangleNormal, normals[indices[i]]) < 0.0f)
        {
            triangleNormal = -triangleNormal;
        }

        normals[indices[i]] += triangleNormal;
        normals[indices[i + 1]] += triangleNormal;
        normals[indices[i + 2]] += triangleNormal;
    }

    // Normalize the accumulated normals
    for (auto &normal : normals)
    {
        if (glm::length(normal) > 0.0f)
        {
            normal = glm::normalize(normal);
        }

        // Invert when necessary
        if (invertNormals)
        {
            normal = -normal;
        }
    }
}

void extrudeSurface(const std::vector<glm::vec3> &surfaceVertices, std::vector<glm::vec3> &surfaceNormals, float extrusionDistance, std::vector<glm::vec3> &extrudedVertices, std::vector<unsigned int> &topSurfaceIndices, std::vector<unsigned int> &bottomSurfaceIndices, std::vector<unsigned int> &frontSurfaceIndices, std::vector<unsigned int> &endSurfaceIndices, std::vector<unsigned int> &sideSurface1Indices, std::vector<unsigned int> &sideSurface2Indices)
{
    size_t numVertices = surfaceVertices.size();
    size_t numExtrudedVertices = 2 * numVertices;

    extrudedVertices.reserve(numExtrudedVertices);

    // Create vertices for the top surface of the extrusion
    for (size_t i = 0; i < numVertices; ++i)
    {
        // Extrude each vertex along its normal direction
        glm::vec3 extrudedVertex = surfaceVertices[i] + extrusionDistance * surfaceNormals[i];
        extrudedVertices.push_back(extrudedVertex);
    }

    // Create vertices for the bottom surface of the extrusion
    for (size_t i = 0; i < numVertices; ++i)
    {
        // Bottom surface vertices are the same as the original surface vertices but at a lower position
        glm::vec3 bottomVertex = surfaceVertices[i];
        extrudedVertices.push_back(bottomVertex);
    }

    for (size_t i = 0; i < numVertices; ++i)
    {
        topSurfaceIndices.push_back(i);
    }

    for (size_t i = 0; i < numVertices; ++i)
    {
        bottomSurfaceIndices.push_back(numVertices + i);
    }

    // create indices for front + end surface
    frontSurfaceIndices.push_back(0);
    frontSurfaceIndices.push_back(numVertices);
    frontSurfaceIndices.push_back(1);
    frontSurfaceIndices.push_back(numVertices + 1);

    endSurfaceIndices.push_back(numVertices - 1);
    endSurfaceIndices.push_back(2 * (numVertices)-1);
    endSurfaceIndices.push_back(numVertices - 2);
    endSurfaceIndices.push_back(2 * (numVertices)-2);

    // create indices for side faces
    for (size_t i = 0; i < numVertices; i += 2)
    {
        sideSurface1Indices.push_back(i);
        sideSurface1Indices.push_back(i + numVertices);
    }
    for (size_t i = 1; i < numVertices; i += 2)
    {
        sideSurface2Indices.push_back(i);
        sideSurface2Indices.push_back(i + numVertices);
    }
}

///=========================================================================================///
///                                       Helper Functions for Distance
///=========================================================================================///

float calculateAverageDistance(const std::vector<float> &vertices)
{
    if (vertices.size() < 6)
        return 1.0f;

    float totalDist = 0.0f;
    int numDistances = 0;

    // for now just leave it at 6
    for (size_t i = 0; i < vertices.size() - 3; i += 3)
    {
        glm::vec3 p1(vertices[i], vertices[i + 1], vertices[i + 2]);
        glm::vec3 p2(vertices[i + 3], vertices[i + 4], vertices[i + 5]);

        float dist = glm::distance(p1, p2);
        totalDist += dist;
        numDistances++;
    }

    return numDistances > 0 ? totalDist / numDistances : 1.0f;
}

///=========================================================================================///
///                                       Export OBJ
///=========================================================================================///

void exportToObj(const std::vector<glm::vec3> &extrudedVertices, const std::string &filename, std::vector<glm::vec3> allNormals, int startingIndices[],
                 std::vector<unsigned int> &topSurfaceIndices, std::vector<unsigned int> &bottomSurfaceIndices, std::vector<unsigned int> &frontSurfaceIndices, std::vector<unsigned int> &endSurfaceIndices, std::vector<unsigned int> &sideSurface1Indices, std::vector<unsigned int> &sideSurface2Indices)
{

    size_t numVertices = extrudedVertices.size() / 2;

    std::ofstream outputFile(filename);
    if (!outputFile.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << " for writing\n";
        return;
    }

    // iterate over vertices
    for (const auto &vertex : extrudedVertices)
    {
        outputFile << "v " << vertex.x << " " << vertex.y << " " << vertex.z << "\n";
    }

    // iterate over normals
    for (const auto &normal : allNormals)
    {
        outputFile << "vn " << normal.x << " " << normal.y << " " << normal.z << "\n";
    }

    auto exportFaces = [&outputFile](const std::vector<unsigned int> &surfaceIndices, int start)
    {
        // Iterate through each triangle strip
        for (int i = 0; i < surfaceIndices.size() - 2; ++i)
        {
            int v1, v2, v3;
            int vn1, vn2, vn3;
            if (i % 2 == 0) // this accounts for winding order
            {
                v1 = surfaceIndices[i] + 1;
                v2 = surfaceIndices[i + 1] + 1;
                v3 = surfaceIndices[i + 2] + 1;

                // determine normal indices for current face
                vn1 = start + i + 1;
                vn2 = start + i + 2;
                vn3 = start + i + 3;
            }
            else
            {
                v1 = surfaceIndices[i + 2] + 1;
                v2 = surfaceIndices[i + 1] + 1;
                v3 = surfaceIndices[i] + 1;

                // determine normal indices for current face
                vn1 = start + i + 3;
                vn2 = start + i + 2;
                vn3 = start + i + 1;
            }

            // Export face with vertex indices
            outputFile << "f " << v1 << "//" << vn1 << " " << v2 << "//" << vn2 << " " << v3 << "//" << vn3 << "\n";
        }
    };

    exportFaces(topSurfaceIndices, startingIndices[0]);
    exportFaces(bottomSurfaceIndices, startingIndices[1]);
    exportFaces(frontSurfaceIndices, startingIndices[2]);
    exportFaces(endSurfaceIndices, startingIndices[3]);
    exportFaces(sideSurface1Indices, startingIndices[4]);
    exportFaces(sideSurface2Indices, startingIndices[5]);

    outputFile.close();
}

// Had to add own export OBJ
void exportTubeToObj(const std::vector<glm::vec3> &vertices,
                     const std::vector<unsigned int> &indices,
                     const std::vector<glm::vec3> &normals,
                     const std::string &filename)
{
    std::ofstream outputFile(filename);
    if (!outputFile.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << " for writing\n";
        return;
    }

    // Write header/comment
    outputFile << "# Harmonograph tube export\n";
    outputFile << "# Vertices: " << vertices.size() << "\n";
    outputFile << "# Faces: " << indices.size() / 3 << "\n\n";

    // Write vertices
    for (const auto &vertex : vertices)
    {
        outputFile << "v " << vertex.x << " " << vertex.y << " " << vertex.z << "\n";
    }

    // Write normals
    for (const auto &normal : normals)
    {
        outputFile << "vn " << normal.x << " " << normal.y << " " << normal.z << "\n";
    }

    // Write faces (triangles)
    // OBJ is 1-indexed, so add 1 to each index
    for (size_t i = 0; i < indices.size(); i += 3)
    {
        outputFile << "f ";
        // For each vertex in the triangle
        for (size_t j = 0; j < 3; ++j)
        {
            // Format: vertex_index/texture_index/normal_index
            // We're using the same index for vertices and normals
            outputFile << (indices[i + j] + 1) << "//" << (indices[i + j] + 1) << " ";
        }
        outputFile << "\n";
    }

    outputFile.close();
    std::cout << "Tube successfully exported to " << filename << std::endl;
}

///=========================================================================================///
///                             Helper Functions for VBO
///=========================================================================================///

void generateAndBindIndexVBO(unsigned int &indexVBO, const std::vector<unsigned int> &indices)
{
    glGenBuffers(1, &indexVBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);
}

// Function to generate and bind a VBO for normals
void generateAndBindNormalVBO(unsigned int &normalVBO, const std::vector<glm::vec3> &normals)
{
    glGenBuffers(1, &normalVBO);
    glBindBuffer(GL_ARRAY_BUFFER, normalVBO);
    glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(glm::vec3), normals.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(1); // Assuming the attribute location for normals is 1
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void *)0);
}

// Draw surface using indices and normals VBOs
void drawSurface(unsigned int indexVBO, unsigned int normalVBO, const std::vector<unsigned int> &indices, GLsizei count)
{
    glBindBuffer(GL_ARRAY_BUFFER, normalVBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO);
    glDrawElements(GL_TRIANGLE_STRIP, count, GL_UNSIGNED_INT, 0);
}

///=========================================================================================///
///                                      Harmonograph Function
///=========================================================================================///

void drawHarmonograph(float animationTime, bool reconstructSurface, unsigned int shaderProgram)
{
    // Buffers for harmonograph line segments
    unsigned int VBO, VAO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);

    // Draw the harmonograph
    std::vector<float> vertices;
    float x, y, z, time;

    for (time = 0; time < animationTime; time += 0.01)
    {
        x = amplitude * sin(time * freqPtr1[0] + phasePtr1[0]) * exp(-dampPtr1[0] * time) + sin(time * freqPtr2[0] + phasePtr2[0]) * exp(-dampPtr2[0] * time);
        y = amplitude * sin(time * freqPtr1[1] + phasePtr1[1]) * exp(-dampPtr1[1] * time) + sin(time * freqPtr2[1] + phasePtr2[1]) * exp(-dampPtr2[1] * time);
        z = amplitude * sin(time * freqPtr1[2] + phasePtr1[2]) * exp(-dampPtr1[2] * time) + sin(time * freqPtr2[2] + phasePtr2[2]) * exp(-dampPtr2[2] * time);

        vertices.push_back(x);
        vertices.push_back(y);
        vertices.push_back(z);
    }

    // Upload vertices data to GPU
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), &vertices[0], GL_STATIC_DRAW);

    // Set attribute pointers
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);

    // Draw line segments -> originally uncommented
    // glDrawArrays(GL_LINE_STRIP, 0, vertices.size() / 3);

    // Draw point cloud
    glPointSize(3); // 5
    glDrawArrays(GL_POINTS, 0, vertices.size() / 3);

    // Extrude surface
    if (reconstructSurface)
    {
        switch (reconstructMethod)
        {
        // Extrusion along normals (Previous implementation from module)
        case 0:
        {
            std::vector<glm::vec3> lineVertices;
            for (size_t i = 0; i < vertices.size(); i += 3)
            {
                lineVertices.push_back(glm::vec3(vertices[i], vertices[i + 1], vertices[i + 2]));
            }

            // Calculate normals
            std::vector<glm::vec3> normals = calculateNormals(vertices);

            // buffers for normals
            unsigned int normalVBO, normalVAO;
            glGenVertexArrays(1, &normalVAO);
            glGenBuffers(1, &normalVBO);
            glBindVertexArray(normalVAO);
            glBindBuffer(GL_ARRAY_BUFFER, normalVBO);

            // Store normal vectors' data in VBO
            glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(glm::vec3), &normals[0], GL_STATIC_DRAW);

            // Set vertex attribute pointer for position
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void *)0);
            glEnableVertexAttribArray(0);

            // Draw the normal vectors
            // Calculate vertices for line segments (control points to normals)
            std::vector<glm::vec3> surfaceVertices; // to store vertices of point x on line, normal of said point, point x+1 on line (cont)
            for (size_t i = 0; i < lineVertices.size() - 4; ++i)
            {
                surfaceVertices.push_back(lineVertices[i]);
                surfaceVertices.push_back(lineVertices[i] + normals[i]); // 0.5f * normals[i]);
            }

            // Store surface vertices' data in VBO
            glBufferData(GL_ARRAY_BUFFER, surfaceVertices.size() * sizeof(glm::vec3), &surfaceVertices[0], GL_STATIC_DRAW);

            // Set vertex attribute pointer for line segment vertices
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void *)0);
            glEnableVertexAttribArray(0);

            // glDrawArrays(GL_LINES, 0, surfaceVertices.size());
            // glPointSize(5.0f);
            // glDrawArrays(GL_POINTS, 0, normals.size() * 2); // this draws out the normal end points

            std::vector<glm::vec3> surfaceNormals(surfaceVertices.size());
            std::vector<unsigned int> indices(surfaceVertices.size());
            std::iota(indices.begin(), indices.end(), 0);
            calculateTriangleStripNormals(surfaceVertices, indices, surfaceNormals);

            // Extrude surface
            std::vector<glm::vec3> extrudedVertices;
            float extrusionDistance = 0.1; // adjust extrusion here
            std::vector<unsigned int> topSurfaceIndices, bottomSurfaceIndices, frontSurfaceIndices, endSurfaceIndices, sideSurface1Indices, sideSurface2Indices;

            extrudeSurface(surfaceVertices, surfaceNormals, extrusionDistance, extrudedVertices, topSurfaceIndices, bottomSurfaceIndices, frontSurfaceIndices, endSurfaceIndices, sideSurface1Indices, sideSurface2Indices);

            std::vector<glm::vec3> topSurfaceNormals, bottomSurfaceNormals, frontSurfaceNormals, endSurfaceNormals, sideSurface1Normals, sideSurface2Normals;

            calculateTriangleStripNormals(extrudedVertices, topSurfaceIndices, topSurfaceNormals, true);
            calculateTriangleStripNormals(extrudedVertices, bottomSurfaceIndices, bottomSurfaceNormals);
            calculateTriangleStripNormals(extrudedVertices, frontSurfaceIndices, frontSurfaceNormals);
            calculateTriangleStripNormals(extrudedVertices, endSurfaceIndices, endSurfaceNormals);
            calculateTriangleStripNormals(extrudedVertices, sideSurface1Indices, sideSurface1Normals);
            calculateTriangleStripNormals(extrudedVertices, sideSurface2Indices, sideSurface2Normals);

            // Create and bind VAO and VBO for extruded surface
            unsigned int extrudedVAO, extrudedVBO;
            glGenVertexArrays(1, &extrudedVAO);
            glGenBuffers(1, &extrudedVBO);
            glBindVertexArray(extrudedVAO);
            glBindBuffer(GL_ARRAY_BUFFER, extrudedVBO);
            glBufferData(GL_ARRAY_BUFFER, extrudedVertices.size() * sizeof(glm::vec3), extrudedVertices.data(), GL_STATIC_DRAW);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void *)0);
            glEnableVertexAttribArray(0);

            // Create and bind VBOs for surface indices
            unsigned int topSurfaceIndexVBO, bottomSurfaceIndexVBO, frontSurfaceIndexVBO, endSurfaceIndexVBO, sideSurface1IndexVBO, sideSurface2IndexVBO;
            unsigned int topSurfaceNormalsVBO, bottomSurfaceNormalsVBO, frontSurfaceNormalsVBO, endSurfaceNormalsVBO, sideSurface1NormalsVBO, sideSurface2NormalsVBO;

            // Generate and bind VBOs for surface indices
            generateAndBindIndexVBO(topSurfaceIndexVBO, topSurfaceIndices);
            generateAndBindIndexVBO(bottomSurfaceIndexVBO, bottomSurfaceIndices);
            generateAndBindIndexVBO(frontSurfaceIndexVBO, frontSurfaceIndices);
            generateAndBindIndexVBO(endSurfaceIndexVBO, endSurfaceIndices);
            generateAndBindIndexVBO(sideSurface1IndexVBO, sideSurface1Indices);
            generateAndBindIndexVBO(sideSurface2IndexVBO, sideSurface2Indices);

            // Generate and bind VBOs for surface normals
            generateAndBindNormalVBO(topSurfaceNormalsVBO, topSurfaceNormals);
            generateAndBindNormalVBO(bottomSurfaceNormalsVBO, bottomSurfaceNormals);
            generateAndBindNormalVBO(frontSurfaceNormalsVBO, frontSurfaceNormals);
            generateAndBindNormalVBO(endSurfaceNormalsVBO, endSurfaceNormals);
            generateAndBindNormalVBO(sideSurface1NormalsVBO, sideSurface1Normals);
            generateAndBindNormalVBO(sideSurface2NormalsVBO, sideSurface2Normals);

            // Draw each surface
            drawSurface(topSurfaceIndexVBO, topSurfaceNormalsVBO, topSurfaceIndices, topSurfaceIndices.size());
            drawSurface(bottomSurfaceIndexVBO, bottomSurfaceNormalsVBO, bottomSurfaceIndices, bottomSurfaceIndices.size());
            drawSurface(frontSurfaceIndexVBO, frontSurfaceNormalsVBO, frontSurfaceIndices, frontSurfaceIndices.size());
            drawSurface(endSurfaceIndexVBO, endSurfaceNormalsVBO, endSurfaceIndices, endSurfaceIndices.size());
            drawSurface(sideSurface1IndexVBO, sideSurface1NormalsVBO, sideSurface1Indices, sideSurface1Indices.size());
            drawSurface(sideSurface2IndexVBO, sideSurface2NormalsVBO, sideSurface2Indices, sideSurface2Indices.size());

            // Cleanup after rendering the extruded surface
            glDeleteVertexArrays(1, &extrudedVAO);
            glDeleteBuffers(1, &extrudedVBO);
            glDeleteBuffers(1, &topSurfaceIndexVBO);
            glDeleteBuffers(1, &bottomSurfaceIndexVBO);

            // combining all normals into allNormals - this is for obj export
            std::vector<std::vector<glm::vec3>> allNormalsList = {
                topSurfaceNormals, bottomSurfaceNormals,
                frontSurfaceNormals, endSurfaceNormals,
                sideSurface1Normals, sideSurface2Normals};

            std::vector<glm::vec3> allNormals;

            int currentIndex = 0;

            constexpr int numGroups = 6;
            int startingIndices[numGroups] = {0};
            allNormals.clear();

            for (int i = 0; i < numGroups; ++i)
            {
                startingIndices[i] = currentIndex;
                allNormals.insert(allNormals.end(), allNormalsList[i].begin(), allNormalsList[i].end());
                currentIndex += allNormalsList[i].size();
            }

            if (isExported)
            {
                // std::cout << "exporting";
                exportToObj(extrudedVertices, "harmonograph_object.obj", allNormals, startingIndices, topSurfaceIndices, bottomSurfaceIndices, frontSurfaceIndices, endSurfaceIndices, sideSurface1Indices, sideSurface2Indices);
                isExported = false;
            }
            break;
        }

        // Reference: https://legacy.cs.indiana.edu/ftp/techreports/TR425.pdf
        // Extrusion along path
        case 1:
        {
            // First approach.... has some crinkles...?
            // std::vector<glm::vec3> tubeVertices;
            // std::vector<unsigned int> tubeIndices;
            // std::vector<glm::vec3> tubeNormals;

            // std::vector<glm::vec3> pathPoints;
            // for (size_t i = 0; i < vertices.size(); i += 3)
            // {
            //     pathPoints.push_back(glm::vec3(vertices[i], vertices[i + 1], vertices[i + 2]));
            // }

            // size_t numPathPoints = pathPoints.size();
            // float tubeRadius = 0.02f;
            // int segments = 24;

            // for (size_t i = 0; i < numPathPoints; ++i)
            // {
            //     glm::vec3 currentPoint = pathPoints[i];
            //     glm::vec3 tangent;
            //     if (i < numPathPoints - 1)
            //     {
            //         tangent = glm::normalize(pathPoints[i + 1] - currentPoint);
            //     }
            //     else
            //     {
            //         tangent = glm::normalize(currentPoint - pathPoints[i - 1]);
            //     }

            //     glm::vec3 perpendicular;

            //     if (std::abs(tangent.x) < std::abs(tangent.y) && std::abs(tangent.x) < std::abs(tangent.z))
            //     {
            //         perpendicular = glm::normalize(glm::cross(tangent, glm::vec3(1, 0, 0)));
            //     }
            //     else if (std::abs(tangent.y) < std::abs(tangent.z))
            //     {
            //         perpendicular = glm::normalize(glm::cross(tangent, glm::vec3(0, 1, 0)));
            //     }
            //     else
            //     {
            //         perpendicular = glm::normalize(glm::cross(tangent, glm::vec3(0, 0, 1)));
            //     }

            //     glm::vec3 binormal = glm::cross(tangent, perpendicular);
            //     binormal = glm::normalize(binormal);

            //     for (int j = 0; j < segments; ++j)
            //     {
            //         float angle = 2.0f * M_PI * j / segments;

            //         glm::vec3 circlePos = currentPoint +
            //                               tubeRadius * cos(angle) * perpendicular +
            //                               tubeRadius * sin(angle) * binormal;

            //         tubeVertices.push_back(circlePos);

            //         glm::vec3 normal = glm::normalize(circlePos - currentPoint);
            //         tubeNormals.push_back(normal);
            //     }
            // }

            // for (size_t i = 0; i < numPathPoints - 1; ++i)
            // {
            //     int baseIndex = i * segments;

            //     for (int j = 0; j < segments; ++j)
            //     {
            //         int nextJ = (j + 1) % segments;

            //         tubeIndices.push_back(baseIndex + j);
            //         tubeIndices.push_back(baseIndex + segments + j);
            //         tubeIndices.push_back(baseIndex + nextJ);

            //         tubeIndices.push_back(baseIndex + nextJ);
            //         tubeIndices.push_back(baseIndex + segments + j);
            //         tubeIndices.push_back(baseIndex + segments + nextJ);
            //     }
            // }

            // New approach (resample + somewhat parallel transport?)
            std::vector<glm::vec3> tubeVertices;
            std::vector<unsigned int> tubeIndices;
            std::vector<glm::vec3> tubeNormals;

            // This is just for visualisation!
            // std::vector<glm::vec3> frameTangents;
            // std::vector<glm::vec3> frameNormals;
            // std::vector<glm::vec3> frameBinormals;
            // std::vector<glm::vec3> framePositions;

            static float tubeRadius = 0.05f; // prev 0.005f for visualizing the TNB
            static int segments = 12;

            std::vector<glm::vec3> pathPoints;
            for (size_t i = 0; i < vertices.size(); i += 3)
            {
                pathPoints.push_back(glm::vec3(vertices[i], vertices[i + 1], vertices[i + 2]));
            }

            std::vector<glm::vec3> resampledPath;
            float desiredSpacing = calculateAverageDistance(vertices) * 0.5f;

            resampledPath.push_back(pathPoints[0]);

            for (size_t i = 1; i < pathPoints.size(); ++i)
            {
                glm::vec3 prev = resampledPath.back();
                glm::vec3 current = pathPoints[i];
                float dist = glm::distance(prev, current);

                if (dist > desiredSpacing)
                {
                    int numSegments = std::ceil(dist / desiredSpacing);
                    for (int j = 1; j < numSegments; ++j)
                    {
                        float t = static_cast<float>(j) / numSegments;
                        glm::vec3 interpolated = prev * (1.0f - t) + current * t;
                        resampledPath.push_back(interpolated);
                    }
                }
                resampledPath.push_back(current);
            }

            pathPoints = resampledPath;
            size_t numPathPoints = pathPoints.size();

            if (numPathPoints < 2)
            {
                std::cout << "Not enough points to create tube" << std::endl;
                break;
            }

            glm::vec3 tangent = glm::normalize(pathPoints[1] - pathPoints[0]);

            glm::vec3 normal;
            if (std::abs(tangent.y) > 0.9f)
            {
                normal = glm::normalize(glm::cross(tangent, glm::vec3(1, 0, 0)));
            }
            else
            {
                normal = glm::normalize(glm::cross(tangent, glm::vec3(0, 1, 0)));
            }

            glm::vec3 binormal = glm::normalize(glm::cross(tangent, normal));

            // This is just for visualization!
            // framePositions.push_back(pathPoints[0]);
            // frameTangents.push_back(tangent);
            // frameNormals.push_back(normal);
            // frameBinormals.push_back(binormal);

            for (size_t i = 0; i < numPathPoints; ++i)
            {
                glm::vec3 currentPoint = pathPoints[i];

                if (i < numPathPoints - 1)
                {
                    glm::vec3 nextTangent = glm::normalize(pathPoints[i + 1] - currentPoint);

                    if (i > 0)
                    {
                        float dot = glm::clamp(glm::dot(tangent, nextTangent), -1.0f, 1.0f);
                        float angle = std::acos(dot);

                        if (angle > 0.001f)
                        {
                            glm::vec3 rotationAxis = glm::normalize(glm::cross(tangent, nextTangent));

                            if (glm::length(rotationAxis) > 0.001f)
                            {
                                glm::mat4 rotation = glm::rotate(glm::mat4(1.0f), angle, rotationAxis);

                                normal = glm::vec3(rotation * glm::vec4(normal, 0.0f));
                                binormal = glm::vec3(rotation * glm::vec4(binormal, 0.0f));

                                // normal = glm::normalize(normal);
                                // binormal = glm::normalize(glm::cross(nextTangent, normal));
                                // normal = glm::normalize(glm::cross(binormal, nextTangent));
                            }
                        }
                    }

                    tangent = nextTangent;

                    // This is just for visualization!
                    // Store this frame (add this line)
                    // if (i > 0)
                    // { // We already stored the first frame
                    //     framePositions.push_back(currentPoint);
                    //     frameTangents.push_back(tangent);
                    //     frameNormals.push_back(normal);
                    //     frameBinormals.push_back(binormal);
                    // }
                }

                for (int j = 0; j < segments; ++j)
                {
                    float angle = 2.0f * M_PI * j / segments;

                    glm::vec3 circlePos = currentPoint + tubeRadius * (cos(angle) * normal + sin(angle) * binormal);

                    tubeVertices.push_back(circlePos);

                    glm::vec3 vertexNormal = glm::normalize(circlePos - currentPoint);
                    tubeNormals.push_back(vertexNormal);
                }
            }

            for (size_t i = 0; i < numPathPoints - 1; ++i)
            {
                int baseIndex = i * segments;

                for (int j = 0; j < segments; ++j)
                {
                    int nextJ = (j + 1) % segments;

                    tubeIndices.push_back(baseIndex + j);            // 0
                    tubeIndices.push_back(baseIndex + segments + j); // 2
                    tubeIndices.push_back(baseIndex + nextJ);        // 1

                    tubeIndices.push_back(baseIndex + nextJ);            // 1
                    tubeIndices.push_back(baseIndex + segments + j);     // 2
                    tubeIndices.push_back(baseIndex + segments + nextJ); // 3
                }
            }

            // Turned off for the frames viz
            unsigned int tubeVAO, tubeVBO, tubeEBO, tubeNormalVBO;
            glGenVertexArrays(1, &tubeVAO);
            glGenBuffers(1, &tubeVBO);
            glGenBuffers(1, &tubeEBO);
            glGenBuffers(1, &tubeNormalVBO);

            glBindVertexArray(tubeVAO);

            glBindBuffer(GL_ARRAY_BUFFER, tubeVBO);
            glBufferData(GL_ARRAY_BUFFER, tubeVertices.size() * sizeof(glm::vec3), tubeVertices.data(), GL_STATIC_DRAW);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void *)0);
            glEnableVertexAttribArray(0);

            glBindBuffer(GL_ARRAY_BUFFER, tubeNormalVBO);
            glBufferData(GL_ARRAY_BUFFER, tubeNormals.size() * sizeof(glm::vec3), tubeNormals.data(), GL_STATIC_DRAW);
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void *)0);
            glEnableVertexAttribArray(1);

            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, tubeEBO);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, tubeIndices.size() * sizeof(unsigned int), tubeIndices.data(), GL_STATIC_DRAW);

            glDrawElements(GL_TRIANGLES, tubeIndices.size(), GL_UNSIGNED_INT, 0);

            glBindBuffer(GL_ARRAY_BUFFER, 0);
            glBindVertexArray(0);
            glDeleteVertexArrays(1, &tubeVAO);
            glDeleteBuffers(1, &tubeVBO);
            glDeleteBuffers(1, &tubeEBO);
            glDeleteBuffers(1, &tubeNormalVBO);

            // This is just for visualization!
            // if (!framePositions.empty())
            // {
            //     std::vector<glm::vec3> frameLines;
            //     float lineLength = tubeRadius * 2.5f;

            //     // Avoid clutter
            //     size_t step = 5; // 1; (shows everything, be careful)

            //     for (size_t i = 0; i < frameTangents.size(); i += step)
            //     {
            //         glm::vec3 point = framePositions[i];

            //         // Add tangent line (red)
            //         frameLines.push_back(point);
            //         frameLines.push_back(point + lineLength * frameTangents[i]);

            //         // Add normal line (green)
            //         frameLines.push_back(point);
            //         frameLines.push_back(point + lineLength * frameNormals[i]);

            //         // Add binormal line (blue)
            //         frameLines.push_back(point);
            //         frameLines.push_back(point + lineLength * frameBinormals[i]);
            //     }

            //     // Draw frame lines
            //     unsigned int frameVAO, frameVBO;
            //     glGenVertexArrays(1, &frameVAO);
            //     glGenBuffers(1, &frameVBO);
            //     glBindVertexArray(frameVAO);
            //     glBindBuffer(GL_ARRAY_BUFFER, frameVBO);
            //     glBufferData(GL_ARRAY_BUFFER, frameLines.size() * sizeof(glm::vec3), frameLines.data(), GL_STATIC_DRAW);
            //     glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void *)0);
            //     glEnableVertexAttribArray(0);
            //     glLineWidth(5.0f);
            //     for (size_t i = 0; i < frameLines.size() / 6; i++)
            //     {
            //         // Tangent (red)
            //         glUniform3f(glGetUniformLocation(shaderProgram, "meshColor1"), 1.0f, 0.0f, 0.0f);
            //         glDrawArrays(GL_LINES, i * 6, 2);

            //         // Normal (green)
            //         glUniform3f(glGetUniformLocation(shaderProgram, "meshColor1"), 0.0f, 1.0f, 0.0f);
            //         glDrawArrays(GL_LINES, i * 6 + 2, 2);

            //         // Binormal (blue)
            //         glUniform3f(glGetUniformLocation(shaderProgram, "meshColor1"), 0.0f, 0.0f, 1.0f);
            //         glDrawArrays(GL_LINES, i * 6 + 4, 2);
            //     }

            //     glBindBuffer(GL_ARRAY_BUFFER, 0);
            //     glBindVertexArray(0);
            //     glDeleteVertexArrays(1, &frameVAO);
            //     glDeleteBuffers(1, &frameVBO);

            if (isExported)
            {
                exportTubeToObj(tubeVertices, tubeIndices, tubeNormals, "harmonograph_smooth_tube.obj");
                isExported = false;
            }

            break;
        }
            //}
        default:
            break;
        }
    }

    // Clean up
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    glDeleteBuffers(1, &VBO);
    glDeleteVertexArrays(1, &VAO);
}

///=========================================================================================///
///                                      Main Function
///=========================================================================================///

int main(void)
{
    GLFWwindow *window;

    // Initialize the library
    if (!glfwInit())
    {
        return -1;
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); // Request OpenGL 3.x
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // Create a windowed mode window and its OpenGL context
    window = glfwCreateWindow(winWidth, winHeight, "Harmonograph", NULL, NULL);

    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);

    ///// setting up the shaders
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback); // correct resize
    glfwSetScrollCallback(window, scroll_callback);                    // scale
    glfwSetCursorPosCallback(window, cursor_pos_callback);             // translate OR rotate
    glfwSetKeyCallback(window, key_callback);                          // change color
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    // tell GLFW to capture the mouse
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);

    // configure multi-sampling
    glEnable(GL_MULTISAMPLE);

    //// build and compile our shader program
    //// ------------------------------------
    // vertex shader
    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);

    // check for shader compile errors
    int success;
    char infoLog[512];
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n"
                  << infoLog << std::endl;
        return -1;
    }

    // fragment shader
    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);

    // check for shader compile errors
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n"
                  << infoLog << std::endl;
        return -1;
    }

    // link shaders
    unsigned int shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    // check for linking errors
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success)
    {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n"
                  << infoLog << std::endl;
        return -1;
    }
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    // Set uniform locations
    int projectionLoc = glGetUniformLocation(shaderProgram, "projection");
    int viewLoc = glGetUniformLocation(shaderProgram, "view");
    int modelLoc = glGetUniformLocation(shaderProgram, "model");
    int lightColorLoc = glGetUniformLocation(shaderProgram, "lightColor");
    int meshColorLoc1 = glGetUniformLocation(shaderProgram, "meshColor1");
    int meshColorLoc2 = glGetUniformLocation(shaderProgram, "meshColor2");
    int viewPosLoc = glGetUniformLocation(shaderProgram, "viewPos");

    unsigned int VBO, VAO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);

    // Enable shader
    glUseProgram(shaderProgram);

    // Set attribute pointers
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);

    // Initialize animation time
    float animationTime = 0.0f;
    // printf("%s\n", glGetString(GL_VERSION));

    const char *glsl_version = "#version 330";

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO &io = ImGui::GetIO();
    //(void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard; // Enable Keyboard Controls
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;  // Enable Gamepad Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    // ImGui::StyleColorsLight();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
#ifdef __EMSCRIPTEN__
    ImGui_ImplGlfw_InstallEmscriptenCanvasResizeCallback("#canvas");
#endif
    ImGui_ImplOpenGL3_Init(glsl_version);

    // Initialize harmonograph values
    freqPtr1[0] = 2.0f; // X component of 1st pendulum
    freqPtr1[1] = 2.0f; // Y component
    freqPtr1[2] = 2.1f; // Z component

    freqPtr2[0] = 3.0f; // X component of 2nd pendulum
    freqPtr2[1] = 3.0f; // Y component
    freqPtr2[2] = 3.1f; // Z component

    dampPtr1[0] = 0.002f;
    dampPtr1[1] = 0.002f;
    dampPtr1[2] = 0.002f;

    dampPtr2[0] = 0.002f;
    dampPtr2[1] = 0.002f;
    dampPtr2[2] = 0.002f;

    phasePtr1[0] = 0;
    phasePtr1[1] = M_PI / 4;
    phasePtr1[2] = 0;

    phasePtr2[0] = M_PI / 2;
    phasePtr2[1] = M_PI / 4;
    phasePtr2[2] = M_PI / 2;

    int freeze = 1;

    // Loop until the user closes the window
    while (!glfwWindowShouldClose(window))
    {
        // Process inputs
        processInput(window);
        glfwPollEvents();

        // Start Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        static int clicked = 0;
        ImGui::Begin("Harmonograph Controls");

        ImGui::SliderFloat("Speed", &animationIncrement, 0.01f, 0.5f, "Speed: %.4f", ImGuiSliderFlags_AlwaysClamp);
        ImGui::SliderFloat("Amplitude", &amplitude, 0.01f, 5.0f, "Amplitude: %.4f", ImGuiSliderFlags_AlwaysClamp);
        ImGui::Text("Pendulum 1");
        ImGui::SliderFloat("Frequency - x", &freqPtr1[0], 0.0f, 4.0f, "%.4f");
        ImGui::SliderFloat("Frequency - y", &freqPtr1[1], 0.0f, 4.0f, "%.4f");
        ImGui::SliderFloat("Frequency - z", &freqPtr1[2], 0.0f, 4.0f, "%.4f");
        ImGui::SliderFloat("Damping - x", &dampPtr1[0], 0.0f, 1.0f, "%.4f");
        ImGui::SliderFloat("Damping - y", &dampPtr1[1], 0.0f, 1.0f, "%.4f");
        ImGui::SliderFloat("Damping - z", &dampPtr1[2], 0.0f, 1.0f, "%.4f");
        ImGui::SliderFloat("Phase - x", &phasePtr1[0], 0.0f, 3.14f, "%.4f");
        ImGui::SliderFloat("Phase - y", &phasePtr1[1], 0.0f, 3.14f, "%.4f");
        ImGui::SliderFloat("Phase - z", &phasePtr1[2], 0.0f, 3.14f, "%.4f");

        ImGui::Text("Pendulum 2");
        ImGui::SliderFloat("Frequency - x", &freqPtr2[0], 0.0f, 4.0f, "%.4f");
        ImGui::SliderFloat("Frequency - y", &freqPtr2[1], 0.0f, 4.0f, "%.4f");
        ImGui::SliderFloat("Frequency - z", &freqPtr2[2], 0.0f, 4.0f, "%.4f");
        ImGui::SliderFloat("Damping - x", &dampPtr2[0], 0.0f, 1.0f, "%.4f");
        ImGui::SliderFloat("Damping - y", &dampPtr2[1], 0.0f, 1.0f, "%.4f");
        ImGui::SliderFloat("Damping - z", &dampPtr2[2], 0.0f, 1.0f, "%.4f");
        ImGui::SliderFloat("Phase - x", &phasePtr2[0], 0.0f, 3.14f, "%.4f");
        ImGui::SliderFloat("Phase - y", &phasePtr2[1], 0.0f, 3.14f, "%.4f");
        ImGui::SliderFloat("Phase - z", &phasePtr2[2], 0.0f, 3.14f, "%.4f");

        ImGui::Text("Presets");
        if (ImGui::Button("Option 1"))
        {
            SetPresets(0);
        }
        ImGui::SameLine();
        if (ImGui::Button("Option 2"))
        {
            SetPresets(1);
        }
        ImGui::SameLine();
        if (ImGui::Button("Option 3"))
        {
            SetPresets(2);
        }

        ImGui::Text("Controls");
        if (ImGui::Button("Freeze"))
        {
            freeze++;
        }
        ImGui::SameLine();
        if (ImGui::Button("Reset"))
        {
            animationTime = 0;
            freeze = 1;
            isAnimating = true;
            reconstructMethod = -1;
        }
        ImGui::SameLine();
        if (ImGui::Button("Change Colour"))
        {
            SetMeshColor();
        }
        ImGui::SameLine();
        if (ImGui::Button("2D <-> 3D"))
        {
            isAnimating = !isAnimating;
        }
        ImGui::SameLine();
        if (ImGui::Button("Export"))
        {
            isExported = true;
        }

        ImGui::Text("Mesh Generation Methods");
        ImGui::RadioButton("Normals Extrusion", &reconstructMethod, 0);
        ImGui::SameLine();
        ImGui::RadioButton("Tube Extrusion", &reconstructMethod, 1);
        ImGui::SameLine();
        ImGui::End();

        // Render OpenGL here
        glClearColor(0.95f, 0.95f, 0.95f, 1.0f); // change background colour
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Set up the viewing transformation
        projection = glm::perspective(glm::radians(camera_fovy), (float)winWidth / (float)winHeight, _Z_NEAR, _Z_FAR);
        glm::mat4 view = glm::lookAt(camera_position, camera_target, camera_up);

        glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, &projection[0][0]);
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, &view[0][0]);
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, &modelMatrix[0][0]);
        glUniform3fv(lightColorLoc, 1, &lightColorTable[lightColorID][0]);
        glUniform3fv(meshColorLoc1, 1, &meshColorTable1[meshColorID1][0]);
        glUniform3fv(meshColorLoc2, 1, &meshColorTable2[meshColorID2][0]);
        glUniform3fv(viewPosLoc, 1, &camera_position[0]);

        drawHarmonograph(animationTime, !isAnimating, shaderProgram);

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        if (freeze & 1)
        {
            animationTime += animationIncrement;
        }

        // Swap front and back buffers
        glfwSwapBuffers(window);
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    // Deallocate all resources
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteProgram(shaderProgram);

    glfwTerminate();

    return 0;
}