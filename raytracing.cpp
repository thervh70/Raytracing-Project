#include <stdio.h>
#ifdef WIN32
#include <windows.h>
#endif
#include "GL/glut.h"
#include "raytracing.h"
#include "config.h"
#include "Vec3D.h"
#include "Matrix33.h"


struct TestRay {
	Vec3Df origin;
	Vec3Df destination;
	Vec3Df color;

	TestRay() 
		: origin(Vec3Df()), destination(Vec3Df()), color(Vec3Df(1, 1, 1)) {};

	TestRay(Vec3Df o, Vec3Df d) 
		: origin(o), destination(d), color(Vec3Df(1, 1, 1)) {};

	TestRay(Vec3Df o, Vec3Df d, Vec3Df c) 
		: origin(o), destination(d), color(c) {};
};
// All the rays in testRay will be drawn in yourDebugDraw().
std::vector<TestRay> testRay;

// built KD tree
AccelTreeNode treeRoot;

// remember if the acceleration tree has been built
bool builtAccelTree = false;
int treeDepth = 0;
int treeNodes = 1;

//use this function for any preprocessing of the mesh.
void init()
{
	//load the mesh file
	//please realize that not all OBJ files will successfully load.
	//Nonetheless, if they come from Blender, they should, if they 
	//are exported as WavefrontOBJ.
	//PLEASE ADAPT THE LINE BELOW TO THE FULL PATH OF THE dodgeColorTest.obj
	//model, e.g., "C:/temp/myData/GraphicsIsFun/dodgeColorTest.obj", 
	//otherwise the application will not load properly
	MyMesh.loadMesh(OBJ_PATH, true);
	MyMesh.computeVertexNormals();

	//one first move: initialize the first light source
	//at least ONE light source has to be in the scene!!!
	//here, we set it to the current location of the camera
	MyLightPositions.push_back(*(new Vec3Df(10, 0, 0)));

	testRay.push_back(TestRay());

	/* FOR TESTING ONLY ~ Maarten
	Matrix33f m(Vec3Df(3, 4, 9), Vec3Df(5, 12, 8), Vec3Df(9, 3, 1));
	std::cout << m << " det:" << m.det() << std::endl;
	std::cout << "X = " << m.solve(Vec3Df(1, 2, 3)) << std::endl;
	system("pause");*/

	// FOR TESTING ONLY ~ Mathias
	/*testRayOrigin = Vec3Df(2.0f, 5.0f, 1.0f);
	testRayDestination = Vec3Df(8.0f, 7.0f, 4.0f);
	performRayTracing(testRayOrigin, testRayDestination);*/

	// Build the acceleration tree
	if (!builtAccelTree)
	{
		buildKDtree();
		builtAccelTree = true;
	}
}

//return the color of your pixel.
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest, int k)
{
	Vec3Df resCol, point, n, v0, v1, v2, vec1, vec2, dir = dest - origin;
	Hitpair hitpair;
	Triangle triangle;
	float t, D, minT = std::numeric_limits<float>::max();
	std::vector<Triangle> triangles;// = MyMesh.triangles;
	Material material;

	// get starting node
	AccelTreeNode currNode = findChildNode(treeRoot, 0, origin);
	AccelTreeNode hitNode = treeRoot;
	Vec3Df currentPoint = origin, currentDestination = dest;
	std::vector<AccelTreeNode*> oldParentList;
	bool gotHit = false;

	while (true)
	{

		// Check all triangles of the parent's which you havem't checked yet and
		// of the current node

		// Create an iterator for the current and old parent lists
		std::vector<AccelTreeNode*> curPar = currNode.parentList;
		int n = 0;

		// skip all parents which you've checked already
		while (n < curPar.size() && n < oldParentList.size() && curPar[n] == oldParentList[n])
		{
			++n;
		}

		// iterate over all triangles which you haven't checked yet
		while (n < curPar.size())
		{
			triangles = (*curPar[n]).triangles;
			for (int i = 0; i < triangles.size(); ++i)
			{

				hitpair = checkHit(triangles[i], origin, dest, minT);

				if (!hitpair.bHit)
					continue;
				hitNode = (*curPar[n]);
				gotHit = true;
				minT = hitpair.res[2];
				triangle = triangles[i];
			}
			++n;
		}

		if (gotHit && (/*currNode.parentList.back() == &hitNode ||*/ !contains(currNode.parentList, hitNode)))
			//found hit
			break;

		oldParentList = currNode.parentList;
		currNode = findNextNode(currNode, currentPoint, currentDestination);

		
/*		if(currentPoint.p[0] < treeRoot.xStart ||
			currentPoint.p[0] > treeRoot.xEnd ||
			currentPoint.p[1] < treeRoot.yStart ||
			currentPoint.p[1] > treeRoot.yEnd ||
			currentPoint.p[2] < treeRoot.zStart ||
			currentPoint.p[2] > treeRoot.zEnd)*/
		if (currNode.parentList.size() == 1)
			//outside of root
			break;
	}

	// get the right material
	for (int i = 0; i < MyMesh.triangles.size();++i)
	{
		if (MyMesh.triangles[i] == triangle)
		{
			material = MyMesh.materials[MyMesh.triangleMaterials[i]];
			break;
		}
	}
	

	/**
	Work In progress, Youri Arkesteijn, trying to get this to work on the same way as it does in scene previeuw with spacebar.
	**/
	Vec3Df intersectionPoint = origin + minT * (dest - origin);
	resCol = material.Kd()*0.04f;
	float angle;

	for (Vec3Df v : MyLightPositions) {
		angle = Vec3Df::cosAngle(intersectionPoint - v, intersectionPoint - origin);
		if (angle < 0)
			angle *= -1;

		//if (angle > 0) {
			resCol += material.Kd()*(1-angle)/MyLightPositions.size();
		//}
	}

//WIP Mathias
	if (material.has_illum()) {
		int illum = material.illum();
		if (illum == 3) {
			// Vertices of the triangle
			v0 = MyMesh.vertices[triangle.v[0]].p;
			v1 = MyMesh.vertices[triangle.v[1]].p;
			v2 = MyMesh.vertices[triangle.v[2]].p;

			// Vector from v2 to v0
			vec1 = v0 - v2;
			// Vector from v2 to v1
			vec2 = v1 - v2;

			n = Vec3Df::crossProduct(vec1, vec2);
			point = dir - 2 * Vec3Df::dotProduct(n, dir) * n;

			performRayTracing(hitpair.hitPoint, point, k++);
		}
	}


	// f(x,y) = (1 - x)*v1 + (x - y)*v2 + y*v3

	return resCol;
}



void yourDebugDraw()
{
	//draw open gl debug stuff
	//this function is called every frame

	//let's draw the mesh
	MyMesh.draw();
	
	//let's draw the lights in the scene as points
	glPushAttrib(GL_ALL_ATTRIB_BITS); //store all GL attributes
	glDisable(GL_LIGHTING);
	glColor3f(1,1,1);
	glPointSize(10);
	glBegin(GL_POINTS);
	for (int i=0;i<MyLightPositions.size();++i)
		glVertex3fv(MyLightPositions[i].pointer());
	glEnd();
	glPopAttrib();//restore all GL attributes
	//The Attrib commands maintain the state. 
	//e.g., even though inside the two calls, we set
	//the color to white, it will be reset to the previous 
	//state after the pop.


	//as an example: we draw the test ray, which is set by the keyboard function
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glDisable(GL_LIGHTING);

	glBegin(GL_LINES);
	for (int i = 0; i < testRay.size(); ++i) {
		glColor3fv(testRay[i].color.pointer());
		glVertex3fv(testRay[i].origin.pointer());
		glVertex3fv(testRay[i].destination.pointer());
	}
	glEnd();

	//we also draw the light positions
	glPointSize(10);
	glBegin(GL_POINTS);
		for (int i = 0; i < MyLightPositions.size(); ++i)
			glVertex3fv(MyLightPositions[i].pointer());
	glEnd();

	glPopAttrib();
	
	//draw whatever else you want...
	////glutSolidSphere(1,10,10);
	////allows you to draw a sphere at the origin.
	////using a glTranslate, it can be shifted to whereever you want
	////if you produce a sphere renderer, this 
	////triangulated sphere is nice for the preview
}


//yourKeyboardFunc is used to deal with keyboard input.
//t is the character that was pressed
//x,y is the mouse position in pixels
//rayOrigin, rayDestination is the ray that is going in the view direction UNDERNEATH your mouse position.
//
//A few keys are already reserved: 
//'L' adds a light positioned at the camera location to the MyLightPositions vector
//'l' modifies the last added light to the current 
//    camera position (by default, there is only one light, so move it with l)
//    ATTENTION These lights do NOT affect the real-time rendering. 
//    You should use them for the raytracing.
//'r' calls the function performRaytracing on EVERY pixel, using the correct associated ray. 
//    It then stores the result in an image "result.ppm".
//    Initially, this function is fast (performRaytracing simply returns 
//    the target of the ray - see the code above), but once you replaced 
//    this function and raytracing is in place, it might take a 
//    while to complete...
void yourKeyboardFunc(char t, int x, int y, const Vec3Df & rayOrigin, const Vec3Df & rayDestination)
{
	AccelTreeNode test;
	AccelTreeNode testNode;
	Vec3Df position = Vec3Df(1.2, 1.2, 1.2), destination = Vec3Df(2,2,2);
	switch (t)
	{
	case ' ':
		//here, as an example, I use the ray to fill in the values for my upper global ray variable
		//I use these variables in the debugDraw function to draw the corresponding ray.
		//try it: Press a key, move the camera, see the ray that was launched as a line.
		testRay[0].origin = rayOrigin;

		// calculate the coordinates where the ray will hit an object
		// and use those coordinates as ray destination
		testRay[0].destination = calculateIntersectionPoint(rayOrigin, rayDestination);

		// make the ray the color of the intersection point (and slightly brighter)
		testRay[0].color = performRayTracing(testRay[0].origin, testRay[0].destination, 0) + Vec3Df(0.25, 0.25, 0.25);
/*
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;		

		std::cout << "Tree: " << std::endl;
		std::cout << "kdTree location: " << &treeRoot << std::endl;
		std::cout << "Parent list size: " << treeRoot.parentList.size() << std::endl;
		std::cout << "Parent list [0]: " << treeRoot.parentList[0] << std::endl;
		std::cout << "Left child location: " << treeRoot.leftChild << std::endl;
		std::cout << "Right child location: " << treeRoot.rightChild << std::endl;
		std::cout << "Triangle list size: " << treeRoot.triangles.size() << std::endl;
		std::cout << "xStart: " << treeRoot.xStart << std::endl;
		std::cout << "yStart: " << treeRoot.yStart << std::endl;
		std::cout << "zStart: " << treeRoot.zStart << std::endl;
		std::cout << "xEnd: " << treeRoot.xEnd << std::endl;
		std::cout << "yEnd: " << treeRoot.yEnd << std::endl;
		std::cout << "zEnd: " << treeRoot.zEnd << std::endl;

		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;

		std::cout << "Tree -> Left: " << std::endl;
		std::cout << "Parent list size: " << (*treeRoot.leftChild).parentList.size() << std::endl;
		std::cout << "Parent list [0]: " << (*treeRoot.leftChild).parentList[0] << std::endl;
		std::cout << "Parent list [1]: " << (*treeRoot.leftChild).parentList[1] << std::endl;
		std::cout << "Triangle list size: " << (*treeRoot.leftChild).triangles.size() << std::endl;
		std::cout << "xStart: " << (*treeRoot.leftChild).xStart << std::endl;
		std::cout << "yStart: " << (*treeRoot.leftChild).yStart << std::endl;
		std::cout << "zStart: " << (*treeRoot.leftChild).zStart << std::endl;
		std::cout << "xEnd: " << (*treeRoot.leftChild).xEnd << std::endl;
		std::cout << "yEnd: " << (*treeRoot.leftChild).yEnd << std::endl;
		std::cout << "zEnd: " << (*treeRoot.leftChild).zEnd << std::endl;
		std::cout << "left child: " << (*treeRoot.leftChild).leftChild << std::endl;
		std::cout << "right child: " << (*treeRoot.leftChild).rightChild << std::endl;

		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;

		std::cout << "Tree -> Right: " << std::endl;
		std::cout << "Parent list size: " << (*treeRoot.rightChild).parentList.size() << std::endl;
		std::cout << "Parent list [0]: " << (*treeRoot.rightChild).parentList[0] << std::endl;
		std::cout << "Parent list [1]: " << (*treeRoot.rightChild).parentList[1] << std::endl;
		std::cout << "Triangle list size: " << (*treeRoot.rightChild).triangles.size() << std::endl;
		std::cout << "xStart: " << (*treeRoot.rightChild).xStart << std::endl;
		std::cout << "yStart: " << (*treeRoot.rightChild).yStart << std::endl;
		std::cout << "zStart: " << (*treeRoot.rightChild).zStart << std::endl;
		std::cout << "xEnd: " << (*treeRoot.rightChild).xEnd << std::endl;
		std::cout << "yEnd: " << (*treeRoot.rightChild).yEnd << std::endl;
		std::cout << "zEnd: " << (*treeRoot.rightChild).zEnd << std::endl;
		std::cout << "left child: " << (*treeRoot.rightChild).leftChild << std::endl;
		std::cout << "right child: " << (*treeRoot.rightChild).rightChild << std::endl;
		*/
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;

		std::cout << "Tree depth: " << treeDepth << std::endl;
		std::cout << "Amount of nodes: " << treeNodes << std::endl;


		testNode = findChildNode(treeRoot, 0, Vec3Df(1.2, 1.2, 1.2));
		test = findNextNode(testNode, position, destination); //findChildNode(treeRoot, 0, position);

		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;

		for (Vec3Df v : MyLightPositions) {
			std::cout << "Light position: " << v << std::endl;
		}

		std::cout << "Origin      " << testRay[0].origin << std::endl;
		std::cout << "Destination " << testRay[0].destination << std::endl;
		std::cout << "Color       " << testRay[0].color << std::endl;
		
		for (Vec3Df v : MyLightPositions) {
			std::cout << "  Light " << Vec3Df::cosAngle(testRay[0].destination - v, testRay[0].destination - rayOrigin) << " from " << v << std::endl;
		}
		break;

	case 'c':
		//This will generate a SINGLE light source on the position of the camera.
		MyLightPositions.clear();
		MyLightPositions.push_back(MyCameraPosition);
		std::cout << "Cleared light sources " << std::endl;
		break;
	}
}

// Return the intersection point of the ray with the object in the scene
// return the rayDest if no object was hit
inline Vec3Df calculateIntersectionPoint(const Vec3Df & rayOrigin, const Vec3Df & rayDest)
{
	// Get the triangles
	std::vector<Triangle> triangles = MyMesh.triangles;
	Triangle *closestTriangle = nullptr;
	float minT = std::numeric_limits<float>::max();
	Hitpair hitpair;

	// Check for each triangle if the ray hits it,
	// after this loop, the closest triangle will be stored in closestTriangle
	for (int i = 0; i < triangles.size(); ++i)
	{
		hitpair = checkHit(triangles[i], rayOrigin, rayDest, minT);
		if (!hitpair.bHit)
			continue;
		minT = hitpair.res[2];
		closestTriangle = &triangles[i];
	}

	// If no intersection, return the rayDest
	if (closestTriangle == nullptr)
	{
		return rayDest;
	}

	// Return the intersecton point with the triangle
	return rayOrigin + minT * (rayDest - rayOrigin);
}

// check if the triangle is hit by the ray from origin to dest and closer to the origin than minT
inline Hitpair checkHit(const Triangle & triangle, const Vec3Df & origin, const Vec3Df & dest, float minT)
{
	// local variables
	Vec3Df v0, v1, v2, vec1, vec2, vec3, dir = dest - origin;
	Hitpair result;

	// Vertices of the triangle
	v0 = MyMesh.vertices[triangle.v[0]].p;
	v1 = MyMesh.vertices[triangle.v[1]].p;
	v2 = MyMesh.vertices[triangle.v[2]].p;

	// Vector from v2 to v0
	vec1 = v0 - v2;
	// Vector from v2 to v1
	vec2 = v1 - v2;
	// Vector from origin to v0
	vec3 = v0 - origin;

	Matrix33f matrix(vec1, vec2, dir);
	Vec3Df res2 = matrix.solve(origin - v2);

	/*		if (&res2 == Matrix33f::getBadVec()) {
	continue;
	}*/

	Vec3Df res(res2[0], res2[1], -res2[2]); // res[2] = -res[2].
											// This is done in this way because you can't modify vector entries in a thread apparently

	//Check if hit
	//res = [a, b, t]
	if (res[0] < 0 || res[1] < 0 || res[0] + res[1] > 1 || res[2] < 0 || res[2] > minT)
		result.bHit = false;
	else
		result.bHit = true;
	
	result.res = res;
	result.hitPoint = origin + res[2]*dir;
	return result;
}


// Build a KD tree and store the planes with the triangles in a vector of vectors
void buildKDtree()
{
	std::cout << "Building tree ..." << std::endl;

	float xMin, xMax, yMin, yMax, zMin, zMax;

	// Calculate the bounding box of the scene
	std::vector<Triangle>::const_iterator it = MyMesh.triangles.begin();
	xMin = MyMesh.vertices[(*it).v[0]].p[0];
	xMax = MyMesh.vertices[(*it).v[0]].p[0];
	yMin = MyMesh.vertices[(*it).v[0]].p[1];
	yMax = MyMesh.vertices[(*it).v[0]].p[1];
	zMin = MyMesh.vertices[(*it).v[0]].p[2];
	zMax = MyMesh.vertices[(*it).v[0]].p[2];

	// Get the smallest and largest x, y and z values of the triangles of the mesh
	while (it != MyMesh.triangles.end()) {
		for (int i = 0; i < 3; ++i)
		{
			// check x
			if (MyMesh.vertices[(*it).v[i]].p[0] < xMin)
				xMin = MyMesh.vertices[(*it).v[i]].p[0];
			else if (MyMesh.vertices[(*it).v[i]].p[0]> xMax)
				xMax = MyMesh.vertices[(*it).v[i]].p[0];

			// check y
			if (MyMesh.vertices[(*it).v[i]].p[1] < yMin)
				yMin = MyMesh.vertices[(*it).v[i]].p[1];
			else if (MyMesh.vertices[(*it).v[i]].p[1] > yMax)
				yMax = MyMesh.vertices[(*it).v[i]].p[1];

			// check z
			if (MyMesh.vertices[(*it).v[i]].p[2] < zMin)
				zMin = MyMesh.vertices[(*it).v[i]].p[2];
			else if (MyMesh.vertices[(*it).v[i]].p[2] > zMax)
				zMax = MyMesh.vertices[(*it).v[i]].p[2];
		}
		++it;
	}

	// Assign the found values to the kd tree
	treeRoot = AccelTreeNode(MyMesh.triangles, xMin, xMax, yMin, yMax, zMin, zMax);

	// Add the root node to it's own parent list.
	treeRoot.parentList.push_back(&treeRoot);

	// Split the main node into smaller nodes
	splitSpaces(treeRoot, 0);
	std::cout << "... done building tree" << std::endl;
}

void splitSpaces(AccelTreeNode& tree, int axis) {

	if (axis > treeDepth)
		++treeDepth;

//	std::cout << "Start split" << std::endl;
	std::cout << "Node triangle size: " << tree.triangles.size() << std::endl;
	// split the KDtreeCube into subspaces and recursevily recall on those subspaces

	// stop if there are 40 or less triangles in the current level (!!! 40 is randomly chosen !!!)
	if (tree.triangles.size() < 41)
	{
//		std::cout << "return: too few triangles left: " << &tree << std::endl;
		tree.leftChild = nullptr;
		tree.rightChild = nullptr;
		return;
	}


	// stop if subspace axis which is being divided is smaller than 0.05 (!!! 0.05 is randomly chosen !!!)
	if ((axis % 3 == 0 && tree.xEnd - tree.xStart < 0.05) ||
		(axis % 3 == 1 && tree.yEnd - tree.yStart < 0.05) ||
		(tree.zEnd - tree.zStart < 0.05))
	{
//		std::cout << "return: too small subspace: " << &tree << std::endl;
		tree.leftChild = nullptr;
		tree.rightChild = nullptr;
		return;
	}

	AccelTreeNode *left = new AccelTreeNode(), *right = new AccelTreeNode();
	std::vector<Triangle> *leftTri = new std::vector<Triangle>(), *rightTri = new std::vector<Triangle>(), *bothTri = new std::vector<Triangle>();
	
	float mid;

	// split the spacce in 2 smaller subspaces and recursively call this function on those subspaces
	if (axis % 3 == 0) {
		//split x-axis
		mid = (tree.xStart + tree.xEnd) / 2.0f;

		// split the space in half
		(*left) = AccelTreeNode(tree.xStart, mid, tree.yStart, tree.yEnd, tree.zStart, tree.zEnd);
		(*right) = AccelTreeNode(mid, tree.xEnd, tree.yStart, tree.yEnd, tree.zStart, tree.zEnd);
	}
	else if (axis % 3 == 1) {
		//split y-axis
		mid = (tree.yStart + tree.yEnd) / 2.0f;

		// split the space in half
		(*left) = AccelTreeNode(tree.xStart, tree.xEnd, tree.yStart, mid, tree.zStart, tree.zEnd);
		(*right) = AccelTreeNode(tree.xStart, tree.xEnd, mid, tree.yEnd, tree.zStart, tree.zEnd);
	}
	else if (axis % 3 == 2) {
		//split z-axis
		mid = (tree.zStart + tree.zEnd) / 2.0f;

		// split the space in half
		(*left) = AccelTreeNode(tree.xStart, tree.xEnd, tree.yStart, tree.yEnd, tree.zStart, mid);
		(*right) = AccelTreeNode(tree.xStart, tree.xEnd, tree.yStart, tree.yEnd, mid, tree.zEnd);
	}

	// calculate which triangles should be in the left, right, or both parts.
	for (std::vector<Triangle>::const_iterator it = tree.triangles.begin(); it < tree.triangles.end(); ++it)
	{
		if (MyMesh.vertices[(*it).v[0]].p[axis % 3] < mid && MyMesh.vertices[(*it).v[1]].p[axis % 3] < mid &&
			MyMesh.vertices[(*it).v[2]].p[axis % 3] < mid)
			(*leftTri).push_back(*it);
		else if (MyMesh.vertices[(*it).v[0]].p[axis % 3] > mid && MyMesh.vertices[(*it).v[1]].p[axis % 3] > mid &&
			MyMesh.vertices[(*it).v[2]].p[axis % 3] > mid)
			(*rightTri).push_back(*it);
		else
			(*bothTri).push_back(*it);
	}

	// Assign the triangles to the right node
	tree.triangles = *bothTri;
	(*left).triangles = *leftTri;
	(*right).triangles = *rightTri;

	// Set the parent lists of the left/right child nodes
	(*left).parentList = tree.parentList;
	(*left).parentList.push_back(left);
	(*right).parentList = tree.parentList;
	(*right).parentList.push_back(right);

	// Set left/right nodes as child nodes
	tree.leftChild = left;
	tree.rightChild = right;

	// recursively calculate the KDtree subspaces of the two parts
	splitSpaces((*left), axis + 1);
	++treeNodes;
	splitSpaces((*right), axis + 1);
	++treeNodes;

//	std::cout << "End split" << std::endl;
}


// Recursively find in which child node of the given parent the given position is
// If you don't know a parent, use the treeRoot variable (it's the root of the tree) and use axis = 0
inline AccelTreeNode findChildNode(const AccelTreeNode &parent, int axis, const Vec3Df &position)
{
	// if this node doesn't have any children, we found the right node
	if (parent.leftChild == nullptr)
		return parent;

	// else continue looking
	if (axis % 3 == 0)
	{
		// check x-axis
		if ((*parent.leftChild).xEnd > position.p[0])
			return findChildNode((*parent.leftChild), axis + 1, position);
		else //if ((*parent.rightChild).xEnd < position.p[0])
			return findChildNode((*parent.rightChild), axis + 1, position);
	}
	else if (axis % 3 == 1)
	{
		// check y-axis
		if ((*parent.leftChild).yEnd > position.p[1])
			return findChildNode((*parent.leftChild), axis + 1, position);
		else //if((*parent.rightChild).yEnd < position.p[1])
			return findChildNode((*parent.rightChild), axis + 1, position);
	}
	else //if (axis % 3 == 2)
	{
		// check z-axis
		if ((*parent.leftChild).zEnd > position.p[2])
			return findChildNode((*parent.leftChild), axis + 1, position);
		else //if ((*parent.rightChild).zEnd < position.p[2])
			return findChildNode((*parent.rightChild), axis + 1, position);
	}
	
	//return parent;
}


// Find the next node which the ray will hit
inline AccelTreeNode findNextNode(const AccelTreeNode &curN, Vec3Df &position, Vec3Df &destination)
{
	Vec3Df hitCube = findNodeBoxHitPoint(curN, position, destination);

	// avoid floating point error
	// !!! NEEDS BETTER SOLUTION !!!
	hitCube += (destination / 1000.0f);

	destination += hitCube - position;
	position = hitCube;

	// check if we are outside of the root
	if (position.p[0] > treeRoot.xEnd || position.p[0] < treeRoot.xStart ||
		position.p[1] > treeRoot.yEnd || position.p[1] < treeRoot.yStart ||
		position.p[2] > treeRoot.zEnd || position.p[2] < treeRoot.zStart)
		// return root if we are outside of the root.
		return treeRoot;

		// if not outside root, return the right node
	return findChildNode(treeRoot, 0, hitCube);
}


inline Vec3Df findNodeBoxHitPoint(const AccelTreeNode &curN, const Vec3Df &position, const Vec3Df &destination)
{
	float i, tempx, tempy, tempz;

	if (position.p[0] != destination.p[0])
	{
		// Check the y-z faces of the space
		tempx = curN.xStart - position.p[0];
		i = tempx / (destination.p[0] - position.p[0]);
		if (i > 0) 
		{
			tempy = i * (destination.p[1] - position.p[1]) + position.p[1];
			tempz = i * (destination.p[2] - position.p[2]) + position.p[2];

			// check if succesful
			if (tempy > curN.yStart && tempy < curN.yEnd && tempz > curN.zStart && tempz < curN.zEnd)
				return Vec3Df(curN.xStart, tempy, tempz);
		}
		else
		{
			// check other y-z face
			tempx = curN.xEnd - position.p[0];
			i = tempx / (destination.p[0] - position.p[0]);
			tempy = i * (destination.p[1] - position.p[1]) + position.p[1];
			tempz = i * (destination.p[2] - position.p[2]) + position.p[2];

			// check if succesful
			if (tempy > curN.yStart && tempy < curN.yEnd && tempz > curN.zStart && tempz < curN.zEnd)
				return Vec3Df(curN.xEnd, tempy, tempz);
		}
	}

	if (position.p[1] != destination.p[1])
	{
		// Check the x-z faces of the space
		tempy = curN.yStart - position.p[1];
		i = tempy / (destination.p[1] - position.p[1]);
		if (i > 0)
		{
			tempx = i * (destination.p[0] - position.p[0]) + position.p[0];
			tempz = i * (destination.p[2] - position.p[2]) + position.p[2];

			// check if succesful
			if (tempx > curN.xStart && tempx < curN.xEnd && tempz > curN.zStart && tempz < curN.zEnd)
				return Vec3Df(tempx, curN.yStart, tempz);
		}
		else
		{
			// check other x-z face
			tempy = curN.yEnd - position.p[1];
			i = tempy / (destination.p[1] - position.p[1]);
			tempx = i * (destination.p[0] - position.p[0]) + position.p[0];
			tempz = i * (destination.p[2] - position.p[2]) + position.p[2];

			// check if succesful
			if (tempx > curN.xStart && tempx < curN.xEnd && tempz > curN.zStart && tempz < curN.zEnd)
				return Vec3Df(tempx, curN.yEnd, tempz);
		}
	}

	// Check the x-y faces of the space
	tempz = curN.yStart - position.p[2];
	i = tempz / (destination.p[2] - position.p[2]);
	if (i > 0)
	{
		tempx = i * (destination.p[0] - position.p[0]) + position.p[0];
		tempy = i * (destination.p[1] - position.p[1]) + position.p[1];

		// we must be succesfull (or something went wrong...)
//		if (tempx > curN.xStart && tempx < curN.xEnd && tempy > curN.yStart && tempy < curN.yEnd)
			return Vec3Df(tempx, tempy, curN.zStart);
	}
	else
	{
		// check other x-y face
		tempz = curN.yEnd - position.p[2];
		i = tempz / (destination.p[2] - position.p[2]);
		tempx = i * (destination.p[0] - position.p[0]) + position.p[0];
		tempy = i * (destination.p[1] - position.p[1]) + position.p[1];

		// we must be succesfull (or something went wrong...)
		return Vec3Df(tempx, tempy, curN.zEnd);
	}
}

inline bool contains(const std::vector<AccelTreeNode*> &vec, const AccelTreeNode &element)
{
	for (std::vector<AccelTreeNode*>::const_iterator it = vec.begin(); it != vec.end(); ++it)
		if ((**it).xStart == element.xStart &&
			(**it).yStart == element.yStart &&
			(**it).zStart == element.zStart &&
			(**it).xEnd == element.xEnd &&
			(**it).yEnd == element.yEnd &&
			(**it).zEnd == element.zEnd)
		{
			return true;
		}
	return false;
}
