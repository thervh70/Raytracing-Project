#pragma once

#ifndef KDTREE_H_
#define KDTREE_H_

#include <vector>
#include "Vec3D.h"

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
	std::vector<int> triangles;


	// Constructor with triangles
	AccelTreeNode(std::vector<int> tri, float xS, float xE, float yS, float yE, float zS, float zE)
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

void buildKDtree();
void splitSpaces(AccelTreeNode& tree, const int axis);
float calcBestSplit(AccelTreeNode &tree, int axis);
inline AccelTreeNode findChildNode(const AccelTreeNode &parent, int axis, const Vec3Df &position);
inline AccelTreeNode findNextNode(const AccelTreeNode &curN, const Vec3Df &position, const Vec3Df &destination);
inline bool contains(const std::vector<AccelTreeNode*> &vec, const AccelTreeNode &element);
inline void projectOriginOnRoot(Vec3Df &origin, Vec3Df &dest);
inline Vec3Df calculateProjectionOnRoot(Vec3Df &origin, Vec3Df &dest);

#endif //KDTREE_H_