

#ifndef SHADERSOURCE_H
#define SHADERSOURCE_H

/******************************************************************************/
/******************************   Default Shader ******************************/
/******************************************************************************/

const char* vertexShaderSource =
"#version 330\n"
"uniform mat4 model;\n"
"uniform mat4 view;\n"
"uniform mat4 projection;\n"
"uniform vec3 lightColor;\n"
"uniform vec3 meshColor1;\n"
"uniform vec3 meshColor2;\n"
"layout (location = 0) in vec3 aPos;\n"
"layout (location = 1) in vec3 normal;\n"
"out vec3 lcolor;\n"
"out vec3 mcolor1;\n"
"out vec3 mcolor2;\n"
"out vec3 facenormal;\n"
"out vec3 FragPos;\n"
"void main()\n"
"{\n"
"   gl_Position = projection * view * model * vec4(aPos, 1.0);\n"
"	lcolor = lightColor;\n"
"	mcolor1 = meshColor1;\n"
"	mcolor2 = meshColor2;\n"
"FragPos = vec3(model * vec4(aPos, 1.0));\n"
"   facenormal = mat3(transpose(inverse(model))) * normal;\n"
"}\n";

const char* fragmentShaderSource =
"#version 330\n"
"uniform vec3 viewPos;\n"
"in vec3 lcolor;\n"
"in vec3 mcolor1;\n"
"in vec3 mcolor2;\n"
"in vec3 facenormal;\n"
"in vec3 FragPos;\n"
"out vec4 Fragment;\n"
"void main()\n"
"{\n"

// Calculate lighting
"float ambientStrength = 0.6;"
"vec3 ambient = ambientStrength * lcolor;"

"float diffuseStrength = 0.85;"
"vec3 lightPos = vec3(1.0,1.0,10.0);"
"vec3 norm = normalize(facenormal);"
"vec3 lightDir = normalize(lightPos - FragPos);"
"float diff = max(dot(norm, lightDir), 0.0);"
"vec3 diffuse = diffuseStrength * diff * lcolor;"

"float specularStrength = 0.45;"
"vec3 viewDir = normalize(viewPos - FragPos);"
"vec3 reflectDir = reflect(-lightDir, norm);"
"float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);"
"vec3 specular = specularStrength * spec * lcolor;"

// Calculate gradient interpolation from given mesh colors (using only x and y, no depth cos looks sus)
"vec2 FragCoord = vec2(FragPos[0], FragPos[1]);\n"
"vec2 st = FragCoord/0.8;\n"
"float mixValue = distance(st, vec2(0,1));\n"
"vec3 gradientColor = mix(mcolor1, mcolor2, mixValue);\n"

// Apply noise (using smoothed coordinates and random values -> shader cannot do rand() so hardcode)
"float noise = fract(sin(dot(FragCoord, vec2(12.9898, 78.233)) * 43758.5453)) * 0.1;" // note 0.1 adjusts the magnitude/degree of noise
"vec3 noisedGradientColor = gradientColor + vec3(noise); \n"
"vec3 result = (ambient + diffuse + specular) * noisedGradientColor;\n"

"Fragment = vec4(result, 1.0);\n"
"}\n";


#endif //SHADERSOURCE_H
