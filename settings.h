#pragma once
#include "Vec3D.h"

#ifndef SETTINGS_H_
#define SETTINGS_H_

const float backgroundlighting = 0.34f;
const float diffusePower = 1.5f;
//const float specularHardness = 67.f; // this is stored in the .mtl file
const float ShadowFactor = 0.5f;

const float shadowRadius = 0.1f; // This is the "radius" of the light source
const int shadowSamples = 4; // This results in shadowSamples ^ 3 shadowrays being shot, set to 1 for hard shadows

const int MSAA = 16; //valid: 1, 4, 16; default: 1

#define backgroundColor Vec3Df(0.2f, 0.5f, 1.0f)
#define endOfReflection Vec3Df(0.0f, 0.0f, 0.0f)

#endif //SETTINGS_H_
