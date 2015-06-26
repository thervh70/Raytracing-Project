#ifndef RAYTRACING_Hjdslkjfadjfasljf
#define RAYTRACING_Hjdslkjfadjfasljf
#include <vector>
#include "mesh.h"
#include "KDtree.h"

//Welcome to your MAIN PROJECT...
//THIS IS THE MOST relevant code for you!
//this is an important file, raytracing.cpp is what you need to fill out
//In principle, you can do the entire project ONLY by working in these two files

extern Mesh MyMesh; //Main mesh
extern std::vector<Vec3Df> MyLightPositions;
extern Vec3Df MyCameraPosition; //currCamera
extern const unsigned int WindowSize_X;//window resolution width
extern const unsigned int WindowSize_Y;//window resolution height
extern unsigned int RayTracingResolutionX;  // largeur fenetre
extern unsigned int RayTracingResolutionY;  // largeur fenetre

// built KD tree
extern AccelTreeNode treeRoot;

//use this function for any preprocessing of the mesh.
void init();

//you can use this function to transform a click to an origin and destination
//the last two values will be changed. There is no need to define this function.
//it is defined elsewhere
void produceRay(int x_I, int y_I, Vec3Df & origin, Vec3Df & dest);


//your main function to rewrite
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest, int i, float prev_Ni, float &minT);

//a function to debug --- you can draw in OpenGL here
void yourDebugDraw();

//want keyboard interaction? Here it is...
void yourKeyboardFunc(char t, int x, int y, const Vec3Df & rayOrigin, const Vec3Df & rayDestination);

struct Hitpair {
	bool bHit;
	//res = [a, b, t]
	Vec3Df res;
	Vec3Df hitPoint;
};

struct HitTriangle {
	bool bHit;
	float a, b;
	Vec3Df hitPoint;
	Triangle triangle;
	Material material;
};

struct TestRay {
	Vec3Df origin;
	Vec3Df destination;
	Vec3Df color;
	int type;
	int k;

	TestRay()
		: origin(Vec3Df()), destination(Vec3Df()), color(Vec3Df(1, 1, 1)), type(0), k(0) {};

	TestRay(Vec3Df o, Vec3Df d)
		: origin(o), destination(d), color(Vec3Df(1, 1, 1)), type(0), k(0) {};

	TestRay(Vec3Df o, Vec3Df d, Vec3Df c)
		: origin(o), destination(d), color(c), type(0), k(0) {};

	TestRay(Vec3Df o, Vec3Df d, Vec3Df c, int t, int k)
		: origin(o), destination(d), color(c), type(t), k(k) {};
};

HitTriangle checkHit(const Vec3Df & origin, const Vec3Df & dest);
HitTriangle checkHit(const Vec3Df & origin, const Vec3Df & dest, float &minT);
inline Hitpair checkHit(const Triangle & triangle, const Vec3Df & origin, const Vec3Df & dest, float minT);
inline Vec3Df calculateIntersectionPoint(const Vec3Df & rayOrigin, const Vec3Df & rayDest);

extern inline AccelTreeNode findChildNode(const AccelTreeNode &parent, int axis, const Vec3Df &position);
extern inline AccelTreeNode findNextNode(const AccelTreeNode &curN, const Vec3Df &position, const Vec3Df &destination);
extern inline bool contains(const std::vector<AccelTreeNode*> &vec, const AccelTreeNode &element);

#endif