#include <stdio.h>
#ifdef WIN32
#include <windows.h>
#endif
#include "GL/glut.h"
#include "raytracing.h"
#include "config.h"
#include "Vec3D.h"
#include "Matrix33.h"
#include "settings.h"


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
bool debug = false;

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
	MyLightPositions.push_back(*(new Vec3Df(0, 10, 0)));

	testRay.push_back(TestRay());

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
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest, int k, float prev_Ni)
{

	if (k > 5) return Vec3Df();

	Vec3Df resCol = Vec3Df(), hitPoint, dir = dest - origin;
	Hitpair hitpair;
	Triangle triangle;
	float t, minT = std::numeric_limits<float>::max();
	int triangleIndex;
	std::vector<Triangle> triangles = MyMesh.triangles;

	Material material;
	bool hit = false;
	float a, b;

	// get starting node
	AccelTreeNode currNode = findChildNode(treeRoot, 0, origin);
	AccelTreeNode hitNode = treeRoot;
	std::vector<AccelTreeNode*> oldParentList;
	bool gotHit = false;

	// If current point is outside of the root (camera is far away), then put it inside the root.
/*	if (currentPoint.p[0] < treeRoot.xStart || currentPoint.p[0] > treeRoot.xEnd ||
		currentPoint.p[1] < treeRoot.yStart || currentPoint.p[1] > treeRoot.yEnd ||
		currentPoint.p[2] < treeRoot.zStart || currentPoint.p[2] > treeRoot.zEnd)
	{
		projectOriginOnRoot(currentPoint, currentDestination);
	}*/

	while (true)
	{
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
				hitPoint = hitpair.hitPoint;
				
				a = hitpair.res[0];
				b = hitpair.res[1];
			}
			++n;
		}
		if (gotHit && (currNode.parentList.back() == &hitNode || !contains(currNode.parentList, hitNode)))
			//found hit
			break;

		oldParentList = currNode.parentList;
		currNode = findNextNode(currNode, origin, dest);
		
		if (currNode.parentList.size() == 1)
			//outside of root
			break;
	}

	// get the right material
	for (int i = 0; i < MyMesh.triangles.size();++i)
	{
		if (MyMesh.triangles[i] == triangle)
		{
			hit = true;
			triangleIndex = i;
			material = MyMesh.materials[MyMesh.triangleMaterials[i]];
			break;
		}
	}
	
	
	if (debug)
		testRay[k].destination = hit ? hitPoint : dest;

	if (!hit) return backgroundColor;

	// Normals of three vectors of triangle
	Vec3Df
		normalA = MyMesh.vertices[triangle.v[0]].n,
		normalB = MyMesh.vertices[triangle.v[1]].n,
		normalC = MyMesh.vertices[triangle.v[2]].n;

	// Barycentric vertex normal biliniear interpolation shading mode
	float u = a, v = b, w = 1.f - u - v;
	Vec3Df interpolatedNormal = normalA*u + normalB*v + normalC*w;
	interpolatedNormal.normalize();

	// default lighting in all parts that even are in shadow everywhere.
	// Maarten - This should be Ka, ambient light, but our .mtl files have Ka = (0,0,0).
	resCol = material.Kd()*backgroundlighting;

	float angle, distanceToLight;
	Vec3Df lightDir, viewDir, halfwayVector;
	float specularHighlight = material.Ns();

	// Full shadow reduces light by (0.9,0.9,0.9) when ShadowFactor = 0.9
	Vec3Df shadowRGB = Vec3Df(1.f, 1.f, 1.f) * ShadowFactor / MyLightPositions.size();
	Hitpair shadowHit;

	viewDir = origin - hitPoint;
	viewDir.normalize();

	for (Vec3Df v : MyLightPositions) {
		lightDir = v - hitPoint;
		distanceToLight = lightDir.normalize();

		// Diffuse lighting
		angle = Vec3Df::dotProduct(interpolatedNormal, lightDir); // cos(phi) = a . b / 1 / 1 (vectors are normalized)
		if (angle > 0)
			resCol += material.Kd()*angle * diffusePower; // / distanceToLight;

		// Self Shadows (Dark side of object
		if (Vec3Df::dotProduct(lightDir, interpolatedNormal) <= 0) {
			resCol -= shadowRGB;
		} else {
			// Specular lighting
			halfwayVector = (lightDir + viewDir);
			halfwayVector.normalize();
			angle = Vec3Df::dotProduct(interpolatedNormal, halfwayVector);
			if (angle > 0)
				resCol += material.Ks()*std::pow(angle, specularHighlight);

			for (Triangle t : MyMesh.triangles) {
				//Shadow casted by an object
				shadowHit = checkHit(t, (hitPoint + 0.001f * lightDir), v, std::numeric_limits<float>::max());
				if (shadowHit.bHit) {
					resCol -= shadowRGB;
				}
			}
		}
	}

	// These are the indices of refraction.
	// n1 is the material from where the ray comes from, which is the Ni value of the previous hitPoint,
	// if it has a Ni value, else float 1 (standard for air).
	float n1 = (prev_Ni != 0.0f) ? prev_Ni : 1.0f;
	// n2 is the material of the hitPoint, which is the Ni value if it has a Ni value, else float 1.
	float n2 = (material.has_Ni()) ? material.Ni() : 1.0f;

	if (material.illum() == 3) {

		// reflectdir = dir_of_ray - 2 * ray_projected_on_normal
		const Vec3Df reflectdir = dir - 2 * Vec3Df::dotProduct(interpolatedNormal, dir) * interpolatedNormal,
			newOrigin = hitPoint + 0.001f * reflectdir,
			newDest = hitPoint + reflectdir;

		if (debug)
			testRay.push_back(TestRay(hitPoint, Vec3Df(), Vec3Df()));

		Vec3Df newCol = performRayTracing(newOrigin, newDest, ++k, n2);
		resCol = 0.5*resCol + 0.5*newCol;

		if (debug)
			testRay[k].color = newCol;
	}

	// Normalized viewing vector, 
	// which will be considered as the vector of incidence.
	Vec3Df vIncidence = dir;
	vIncidence.normalize();

	// If the refraction index of a material is greater than 1, 
	// then refraction has to be taken into account.
	if (material.Ni() > 1.0f) {
		// Transmission Ray calculation using a simplified version of the gigantic formula given in the slides:
		// (1); tRay = n1/n2 * vIncidence + ( n1/n2 * cos(thetaIncidence) - sqrt(1 - sin^2(thetaTransmitted) ) ) * normal.
		// with (2); sin^2(thetaTransmitted) = (n1/n2)^2 * (1 - cos^2(thetaIncidence))

		// Refractionindex is based on the division of the two refractive indices 
		// of the involved materials.
		float refractIndex = ( n1 != n2 ) ? n2 : (1 / n2);

		// Cos(Theta) with theta as the angle of incidence calculation, by
		// calculating the dotProduct of the vector of indence and the normal of the hitPoint.
		float cosThetaIncidence = Vec3Df::dotProduct(vIncidence, interpolatedNormal);

		// Sin^2(Theta) calculation, with theta as the angle of transmittance, by using formula 2.
		float sin2ThetaTransmitted = pow(refractIndex, 2) * (1 - (pow(cosThetaIncidence, 2)));

		Vec3Df newOriginR, newDestR;

		// Important: only shoot the tRay if the angle is smaller than the critical angle.
		if (sin2ThetaTransmitted <= 1) {
			// The result is a transmitted ray, by using formula 1.
			const Vec3Df tRay = refractIndex * vIncidence + (refractIndex * cosThetaIncidence - sqrt(1 - sin2ThetaTransmitted)) * interpolatedNormal,
				newOriginR = hitPoint + 0.001f * tRay,
				newDestR = hitPoint + tRay;

			// Debug and Tracing		
			if (debug)
				testRay.push_back(TestRay(hitPoint, Vec3Df(), Vec3Df()));

			Vec3Df newCol = performRayTracing(newOriginR, newDestR, ++k, n2);
			resCol = 0.3*resCol + 0.7*newCol;

			if (debug)
				testRay[k].color = newCol;

		}
		else {
			// The result is just a reflected ray.
			const Vec3Df refRay = vIncidence - 2 * Vec3Df::dotProduct(vIncidence, interpolatedNormal) * interpolatedNormal,
				newOriginR = hitPoint + 0.001f * refRay,
				newDestR = hitPoint + refRay;

			// Debug and Tracing		
			if (debug)
				testRay.push_back(TestRay(hitPoint, Vec3Df(), Vec3Df()));

			Vec3Df newCol = performRayTracing(newOriginR, newDestR, ++k, n2);
			resCol = 0.5*resCol + 0.5*newCol;

			if (debug)
				testRay[k].color = newCol;

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
	Vec3Df position = Vec3Df(0.274594f, 1.89004f, 3.5032f), destination = Vec3Df(0.526746f,0.308055f,0.47311f);
	switch (t)
	{
	case ' ':
		debug = true;
		testRay.clear();
		testRay.push_back(TestRay());

		//here, as an example, I use the ray to fill in the values for my upper global ray variable
		//I use these variables in the debugDraw function to draw the corresponding ray.
		//try it: Press a key, move the camera, see the ray that was launched as a line.
		testRay[0].origin = rayOrigin;

		// calculate the coordinates where the ray will hit an object
		// and use those coordinates as ray destination
		testRay[0].destination = rayDestination;

		// make the ray the color of the intersection point (and slightly brighter)
		testRay[0].color = performRayTracing(testRay[0].origin, testRay[0].destination, 0, 0.0f);

		std::cout << "Test ray trace:" << std::endl;

		for (Vec3Df v : MyLightPositions) {
			std::cout << "Light position: " << v << std::endl;
		}

		std::cout << std::endl;
		for (TestRay r : testRay) {
			std::cout << "Origin      " << r.origin << std::endl;
			std::cout << "Destination " << r.destination << std::endl;
			std::cout << "Color       " << r.color << std::endl << std::endl;
		}
		
		debug = false;
		break;

	case 'c':
		//This will generate a SINGLE light source on the position of the camera.
		MyLightPositions.clear();
		MyLightPositions.push_back(MyCameraPosition);
		std::cout << "Cleared light sources " << std::endl;
		break;

	case 'o':
		std::cout << "Opening result.bmp";
		int res = system("result.bmp");
		std::cout << " exited with code " << res << std::endl;
		break;
	}
}

// check if the triangle is hit by the ray from origin to dest and closer to the origin than minT
inline Hitpair checkHit(const Triangle & triangle, const Vec3Df & origin, const Vec3Df & dest, float minT)
{
	// local variables
	Vec3Df v0, v1, v2, vec1, vec2, dir = dest - origin;
	Hitpair result;

	// Vertices of the triangle
	v0 = MyMesh.vertices[triangle.v[0]].p;
	v1 = MyMesh.vertices[triangle.v[1]].p;
	v2 = MyMesh.vertices[triangle.v[2]].p;

	// Vector from v2 to v0
	vec1 = v0 - v2;
	// Vector from v2 to v1
	vec2 = v1 - v2;

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
	std::cout << std::endl;
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

	// Assign the found values to the kd tree (with a small margin)
	// !!! IMPORTANT: NEED BETTER SOLUTION THAN INCREASING MARGIN TO ALLOW CAMERA TO BE OUTSIDE OF THE ROOT !!!
	treeRoot = AccelTreeNode(MyMesh.triangles, xMin - 5.1f, xMax + 5.1f, yMin - 5.1f, yMax + 5.1f, zMin - 5.1f, zMax + 5.1f);

	// Add the root node to it's own parent list.
	treeRoot.parentList.push_back(&treeRoot);

	// Split the main node into smaller nodes
	splitSpaces(treeRoot, 0);
	std::cout << "... done building tree" << std::endl;
	std::cout << "Tree depth: " << treeDepth << std::endl;
	std::cout << "Amount of nodes: " << treeNodes << std::endl;
	std::cout << std::endl;
}

void splitSpaces(AccelTreeNode& tree, int axis) {

	if (axis > treeDepth)
		++treeDepth;


	// split the KDtreeCube into subspaces and recursevily recall on those subspaces

	// stop if there are 50 or less triangles in the current level (!!! 50 is randomly chosen !!!)
	if (tree.triangles.size() < 51)
	{
		tree.leftChild = nullptr;
		tree.rightChild = nullptr;
		return;
	}


	// stop if subspace axis which is being divided is smaller than 0.01 (!!! 0.01 is randomly chosen !!!)
	if ((axis % 3 == 0 && tree.xEnd - tree.xStart < 0.01) ||
		(axis % 3 == 1 && tree.yEnd - tree.yStart < 0.01) ||
		(tree.zEnd - tree.zStart < 0.01))
	{
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
		else
			return findChildNode((*parent.rightChild), axis + 1, position);
	}
	else if (axis % 3 == 1)
	{
		// check y-axis
		if ((*parent.leftChild).yEnd > position.p[1])
			return findChildNode((*parent.leftChild), axis + 1, position);
		else
			return findChildNode((*parent.rightChild), axis + 1, position);
	}
	else
	{
		// check z-axis
		if ((*parent.leftChild).zEnd > position.p[2])
			return findChildNode((*parent.leftChild), axis + 1, position);
		else
			return findChildNode((*parent.rightChild), axis + 1, position);
	}
}

// Find the point where the ray hits the outside of the current node
inline AccelTreeNode findNextNode(const AccelTreeNode &curN, const Vec3Df &position, const Vec3Df &destination)
{

	Vec3Df dir = destination - position;
	float f, tempx, tempy, tempz;
	AccelTreeNode node = treeRoot;
	
	// Check the y-z faces of the space
	if (dir.p[0] != 0)
	{
		if (dir.p[0] < 0)
		{
			// negative direction:
			// check xStart
			f = (curN.xStart - position.p[0]) / dir.p[0];
			tempy = f * dir.p[1] + position.p[1];
			tempz = f * dir.p[2] + position.p[2];

			// check if succesful
			if (tempy > curN.yStart && tempy < curN.yEnd && tempz > curN.zStart && tempz < curN.zEnd)
			{
				// return the root if we would otherwise go outside of the root.
				if (curN.xStart == treeRoot.xStart)
					return treeRoot;
				
				// return the right node
				for (int axis = 0;; ++axis)
				{
					if (axis % 3 == 0)
					{
						// check x-axis
						if((*node.leftChild).xEnd >= curN.xStart)
							node = *node.leftChild;
						else
							node = *node.rightChild;
					}
					else if (axis % 3 == 1)
					{
						// check y-axis
						if(tempy < (*node.leftChild).yEnd)
							node = *node.leftChild;
						else
							node = *node.rightChild;
					}
					else
					{
						// check z-axis
						if (tempz < (*node.leftChild).zEnd)
							node = *node.leftChild;
						else
							node = *node.rightChild;
					}
					if (node.leftChild == nullptr)
						return node;
				}
			}
		}
		else
		{
			// positive direction:
			// check xEnd
			f = (curN.xEnd - position.p[0]) / dir.p[0];
			tempy = f * dir.p[1] + position.p[1];
			tempz = f * dir.p[2] + position.p[2];

			// check if succesful
			if (tempy > curN.yStart && tempy < curN.yEnd && tempz > curN.zStart && tempz < curN.zEnd)
			{
				// return the root if we would otherwise go outside of the root.
				if (curN.xEnd == treeRoot.xEnd)
					return treeRoot;

				// return the right node
				for (int axis = 0;; ++axis)
				{
					if (axis % 3 == 0)
					{
						// check x-axis
						if ((*node.leftChild).xEnd <= curN.xEnd)
							node = *node.rightChild;
						else
							node = *node.leftChild;
					}
					else if (axis % 3 == 1)
					{
						// check y-axis
						if (tempy < (*node.leftChild).yEnd)
							node = *node.leftChild;
						else
							node = *node.rightChild;
					}
					else
					{
						// check z-axis
						if (tempz < (*node.leftChild).zEnd)
							node = *node.leftChild;
						else
							node = *node.rightChild;
					}
					if (node.leftChild == nullptr)
						return node;
				}
			}
		}
	}

	// Check the x-z faces of the space
	if (dir.p[1] != 0)
	{
		if (dir.p[1] < 0)
		{
			// negative direction:
			// check yStart
			f = (curN.yStart - position.p[1]) / dir.p[1];
			tempx = f * dir.p[0] + position.p[0];
			tempz = f * dir.p[2] + position.p[2];

			// check if succesful
			if (tempx > curN.xStart && tempx < curN.xEnd && tempz > curN.zStart && tempz < curN.zEnd)
			{
				// return the root if we would otherwise go outside of the root.
				if (curN.yStart == treeRoot.yStart)
					return treeRoot;

				// return the right node
				for (int axis = 0;; ++axis)
				{
					if (axis % 3 == 0)
					{
						// check x-axis
						if (tempx < (*node.leftChild).xEnd)
							node = *node.leftChild;
						else
							node = *node.rightChild;
					}
					else if (axis % 3 == 1)
					{
						// check y-axis
						if ((*node.leftChild).yEnd >= curN.yStart)
							node = *node.leftChild;
						else
							node = *node.rightChild;
					}
					else
					{
						// check z-axis
						if (tempz < (*node.leftChild).zEnd)
							node = *node.leftChild;
						else
							node = *node.rightChild;
					}
					if (node.leftChild == nullptr)
						return node;
				}
			}
		}
		else
		{
			// positive direction:
			// check yEnd
			f = (curN.yEnd - position.p[1]) / dir.p[1];
			tempx = f * dir.p[0] + position.p[0];
			tempz = f * dir.p[2] + position.p[2];

			// check if succesful
			if (tempx > curN.xStart && tempx < curN.xEnd && tempz > curN.zStart && tempz < curN.zEnd)
			{
				// return the root if we would otherwise go outside of the root.
				if (curN.yEnd == treeRoot.yEnd)
					return treeRoot;

				// return the right node
				for (int axis = 0;; ++axis)
				{
					if (axis % 3 == 0)
					{
						// check x-axis
						if (tempx < (*node.leftChild).xEnd)
							node = *node.leftChild;
						else
							node = *node.rightChild;
					}
					else if (axis % 3 == 1)
					{
						// check y-axis
						if ((*node.leftChild).yEnd <= curN.yEnd)
							node = *node.rightChild;
						else
							node = *node.leftChild;
					}
					else
					{
						// check z-axis
						if (tempz < (*node.leftChild).zEnd)
							node = *node.leftChild;
						else
							node = *node.rightChild;
					}
					if (node.leftChild == nullptr)
						return node;
				}
			}
		}
	}

	// Check the x-y faces of the space
	if (dir.p[2] != 0)
	{
		if (dir.p[2] < 0)
		{
			// negative direction:
			// check yStart
			f = (curN.zStart - position.p[2]) / dir.p[2];
			tempx = f * dir.p[0] + position.p[0];
			tempy = f * dir.p[1] + position.p[1];

			// return the root if we would otherwise go outside of the root.
			if (curN.zStart == treeRoot.zStart)
				return treeRoot;

			// return the right node
			for (int axis = 0;; ++axis)
			{
				if (axis % 3 == 0)
				{
					// check x-axis
					if (tempx < (*node.leftChild).xEnd)
						node = *node.leftChild;
					else
						node = *node.rightChild;
				}
				else if (axis % 3 == 1)
				{
					// check y-axis
					if (tempy < (*node.leftChild).yEnd)
						node = *node.leftChild;
					else
						node = *node.rightChild;
				}
				else
				{
					// check z-axis
					if ((*node.leftChild).zEnd >= curN.zStart)
						node = *node.leftChild;
					else
						node = *node.rightChild;
				}
				if (node.leftChild == nullptr)
					return node;
			}
		}
		else
		{
			// positive direction:
			// check yEnd
			f = (curN.zEnd - position.p[2]) / dir.p[2];
			tempx = f * dir.p[0] + position.p[0];
			tempy = f * dir.p[1] + position.p[1];

			// return the root if we would otherwise go outside of the root.
			if (curN.zEnd == treeRoot.zEnd)
				return treeRoot;

			// return the right node
			for (int axis = 0;; ++axis)
			{
				if (axis % 3 == 0)
				{
					// check x-axis
					if (tempx < (*node.leftChild).xEnd)
						node = *node.leftChild;
					else
						node = *node.rightChild;
				}
				else if (axis % 3 == 1)
				{
					// check y-axis
					if (tempy < (*node.leftChild).yEnd)
						node = *node.leftChild;
					else
						node = *node.rightChild;
				}
				else
				{
					// check z-axis
					if ((*node.leftChild).zEnd <= curN.zEnd)
						node = *node.rightChild;
					else
						node = *node.leftChild;
				}
				if (node.leftChild == nullptr)
					return node;
			}
		}
	}

	return treeRoot;
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

inline void projectOriginOnRoot(Vec3Df &origin, Vec3Df &dest)
{
	Vec3Df hitRoot = calculateProjectionOnRoot(origin, dest);

	// avoid floating point error
	// !!! NEEDS BETTER SOLUTION !!!
	hitRoot += (dest / 1000.f);

	dest += hitRoot - origin;
	origin = hitRoot;
}

// If the camera is outside of the root, use this to move the origin inside of the root.
// this function looks a lot like the findNodeBoxHitPoint function
inline Vec3Df calculateProjectionOnRoot(Vec3Df &origin, Vec3Df &dest)
{
	float i, tempx, tempy, tempz;

	// check y-z faces of the root
	if (origin.p[0] != dest.p[0] && !(origin.p[0] > treeRoot.xStart && origin.p[0] < treeRoot.xEnd))
	{
		if (std::abs(treeRoot.xStart - origin.p[0]) < std::abs(treeRoot.xEnd - origin.p[0]))
		{
			// closer to xStart
			tempx = treeRoot.xStart - origin.p[0];
			i = tempx / (dest.p[0] - origin.p[0]);

			tempy = i * (dest.p[1] - origin.p[1]) + origin.p[1];
			tempz = i * (dest.p[2] - origin.p[2]) + origin.p[2];

			// check if succesful
			if (tempy > treeRoot.yStart && tempy < treeRoot.yEnd && tempz > treeRoot.zStart && tempz < treeRoot.zEnd)
				return Vec3Df(treeRoot.xStart, tempy, tempz);
		}
		else
		{
			// closer to xEnd
			tempx = treeRoot.xEnd - origin.p[0];
			i = tempx / (dest.p[0] - origin.p[0]);

			tempy = i * (dest.p[1] - origin.p[1]) + origin.p[1];
			tempz = i * (dest.p[2] - origin.p[2]) + origin.p[2];

			// check if succesful
			if (tempy > treeRoot.yStart && tempy < treeRoot.yEnd && tempz > treeRoot.zStart && tempz < treeRoot.zEnd)
				return Vec3Df(treeRoot.xEnd, tempy, tempz);
		}
	}

	// check x-z faces of the root
	if (origin.p[1] != dest.p[1] && !(origin.p[1] > treeRoot.yStart && origin.p[1] < treeRoot.yEnd))
	{
		if (std::abs(treeRoot.yStart - origin.p[1]) < std::abs(treeRoot.yEnd - origin.p[1]))
		{
			// closer to yStart
			tempy = treeRoot.yStart - origin.p[1];
			i = tempy / (dest.p[1] - origin.p[1]);

			tempx = i * (dest.p[0] - origin.p[0]) + origin.p[0];
			tempz = i * (dest.p[2] - origin.p[2]) + origin.p[2];

			// check if succesful
			if (tempx > treeRoot.xStart && tempx < treeRoot.xEnd && tempz > treeRoot.zStart && tempz < treeRoot.zEnd)
				return Vec3Df(tempx, treeRoot.yStart, tempz);
		}
		else
		{
			// closer to yEnd
			tempy = treeRoot.yEnd - origin.p[1];
			i = tempy / (dest.p[1] - origin.p[1]);

			tempx = i * (dest.p[0] - origin.p[0]) + origin.p[0];
			tempz = i * (dest.p[2] - origin.p[2]) + origin.p[2];

			// check if succesful
			if (tempx > treeRoot.xStart && tempx < treeRoot.xEnd && tempz > treeRoot.zStart && tempz < treeRoot.zEnd)
				return Vec3Df(tempx, treeRoot.yEnd, tempz);
		}
	}

	// check x-y faces of the root
	if (origin.p[2] != dest.p[2] && !(origin.p[2] > treeRoot.zStart && origin.p[2] < treeRoot.zEnd))
	{
		if (std::abs(treeRoot.zStart - origin.p[2]) < std::abs(treeRoot.zEnd - origin.p[2]))
		{
			// closer to zStart
			tempz = treeRoot.zStart - origin.p[2];
			i = tempz / (dest.p[2] - origin.p[2]);

			if (tempx > treeRoot.xStart && tempx < treeRoot.xEnd && tempy > treeRoot.yStart && tempy < treeRoot.yEnd)
				return Vec3Df(tempx, tempy, treeRoot.zStart);
		}
		else
		{
			// closer to zEnd
			tempz = treeRoot.zEnd - origin.p[2];
			i = tempz / (dest.p[2] - origin.p[2]);

			tempx = i * (dest.p[0] - origin.p[0]) + origin.p[0];
			tempy = i * (dest.p[1] - origin.p[1]) + origin.p[1];

			if (tempx > treeRoot.xStart && tempx < treeRoot.xEnd && tempy > treeRoot.yStart && tempy < treeRoot.yEnd)
				return Vec3Df(tempx, tempy, treeRoot.zEnd);
		}
	}

	return origin;
}
