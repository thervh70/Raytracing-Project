#ifndef RAYTRACING_Hjdslkjfadjfasljf
#define RAYTRACING_Hjdslkjfadjfasljf
#include <vector>
#include "mesh.h"

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

//use this function for any preprocessing of the mesh.
void init();

//you can use this function to transform a click to an origin and destination
//the last two values will be changed. There is no need to define this function.
//it is defined elsewhere
void produceRay(int x_I, int y_I, Vec3Df & origin, Vec3Df & dest);


//your main function to rewrite
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest, int i, float prev_Ni);

//a function to debug --- you can draw in OpenGL here
void yourDebugDraw();

//want keyboard interaction? Here it is...
void yourKeyboardFunc(char t, int x, int y, const Vec3Df & rayOrigin, const Vec3Df & rayDestination);

struct Hitpair {
	bool bHit;
	Vec3Df res;
	Vec3Df hitPoint;
};

struct AccelTreeNode {
	// A list of all parents of this node and this node (first: root, last: this node)
	// Can be used to check which nodes you have traversed already while checking an adjacent node
	// Tradeof: linear instead of logarithmic space complexity (around 3 * triangles bytes are needed)
	// For a logarithmic speedup (around 0.32 ln(triangles) + 0.51 times faster than with logarithmic space complexity)
	// As long as we don't use over 1 billion triangles, we should be fine.
	std::vector<AccelTreeNode*> parentList;

	// The left and right child of this node
	// nullplt if there is no right child
	AccelTreeNode *leftChild, *rightChild;

	// The bounds defining the dimensions of the space of this node
	// nullplt if there is no left child
	float xStart, xEnd, yStart, yEnd, zStart, zEnd;

	// The triangles contained in this node
	std::vector<Triangle> triangles;

	
	// Constructor with triangles
	AccelTreeNode(std::vector<Triangle> tri, float xS, float xE, float yS, float yE, float zS, float zE)
	{
		triangles = tri;
		xStart = xS;
		xEnd = xE;
		yStart = yS;
		yEnd = yE;
		zStart = zS;
		zEnd = zE;
	}

	// Constructor without triangles
	AccelTreeNode(float xS, float xE, float yS, float yE, float zS, float zE)
	{
		xStart = xS;
		xEnd = xE;
		yStart = yS;
		yEnd = yE;
		zStart = zS;
		zEnd = zE;
	}

	// basic constructor
	AccelTreeNode() {};
};

inline Hitpair checkHit(const Triangle & triangle, const Vec3Df & origin, const Vec3Df & dest, float minT);
inline Vec3Df calculateIntersectionPoint(const Vec3Df & rayOrigin, const Vec3Df & rayDest);

// KD tree functions
void buildKDtree();
void splitSpaces(AccelTreeNode& tree, int axis);
inline AccelTreeNode findChildNode(const AccelTreeNode &parent, int axis, const Vec3Df &position);
inline AccelTreeNode findNextNode(const AccelTreeNode &curN, const Vec3Df &position, const Vec3Df &destination);
//inline Vec3Df findNextNode(const AccelTreeNode &curN, const Vec3Df &position, const Vec3Df &destination);
inline bool contains(const std::vector<AccelTreeNode*> &vec, const AccelTreeNode &element);
inline void projectOriginOnRoot(Vec3Df &origin, Vec3Df &dest);
inline Vec3Df calculateProjectionOnRoot(Vec3Df &origin, Vec3Df &dest);

#endif