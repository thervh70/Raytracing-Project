#pragma once
#include "Vec3D.h"

#ifndef SETTINGS_H_
#define SETTINGS_H_

const float backgroundlighting = 0.34f;
const float diffusePower = 0.5f;
//const float specularHardness = 67.f; // this is stored in the .mtl file
const float ShadowFactor = 0.5f;

const float shadowRadius = 0.9f; // This is the "radius" of the light source
const int shadowSamples = 1; // This results in shadowSamples ^ 3 shadowrays being shot, set to 1 for hard shadows

const int MSAA = 1; //valid: 1, 4, 16; default: 1
const float depthPower = 2.5f; // default: 4.0f

#define backgroundColor Vec3Df(0.2f, 0.5f, 1.0f)
#define endOfReflection Vec3Df(0.0f, 0.0f, 0.0f)

#endif //SETTINGS_H_
