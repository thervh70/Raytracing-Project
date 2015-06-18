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
std::vector<KDtreeCube> kdTree;

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

	// FOR TESTING ONLY ~ Mathias
	/*testRayOrigin = Vec3Df(2.0f, 5.0f, 1.0f);
	testRayDestination = Vec3Df(8.0f, 7.0f, 4.0f);
	performRayTracing(testRayOrigin, testRayDestination);*/
}

//return the color of your pixel.
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest, int k)
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

	for (int i = 0; i < triangles.size(); ++i)
	{
		hitpair = checkHit(triangles[i], origin, dest, minT);

		if (!hitpair.bHit)
			continue;

		hit = true;

		minT = hitpair.res[2];
		triangleIndex = i;
		material = MyMesh.materials[MyMesh.triangleMaterials[i]];
		triangle = triangles[i];
		hitPoint = hitpair.hitPoint;

		/**
		Shadows
		**/

		//Find vectors from current triangle to all lights
		//If dotproduct of surface normal with surface-light vector is < 0
		//Object is shadowing itself so therefore no need to trace.

		//Check whether this vector does hit any other objects (blocking light)
		//If (hitObject) shadow from that light

	}
	if (debug)
		testRay[k].destination = hit ? hitPoint : dest;

	if (!hit) return backgroundColor;

	// Calculate intersection point with triangle from origin
	Vec3Df intersectionPoint = origin + minT * (dest - origin);

	// Normals of three vectors of triangle
	Vec3Df
		a = MyMesh.vertices[MyMesh.triangles[triangleIndex].v[0]].p,
		b = MyMesh.vertices[MyMesh.triangles[triangleIndex].v[1]].p,
		c = MyMesh.vertices[MyMesh.triangles[triangleIndex].v[2]].p,
		NormalA = MyMesh.vertices[MyMesh.triangles[triangleIndex].v[0]].n,
		NormalB = MyMesh.vertices[MyMesh.triangles[triangleIndex].v[1]].n,
		NormalC = MyMesh.vertices[MyMesh.triangles[triangleIndex].v[2]].n;

	/**
	Barycentric vertex normal biliniear interpolation shading mode
	**/
	Vec3Df 
		vp0 = b - a,
		vp1 = c - a,
		vp2 = intersectionPoint - a;
	float d00 = Vec3Df::dotProduct(vp0, vp0);
	float d01 = Vec3Df::dotProduct(vp0, vp1);
	float d11 = Vec3Df::dotProduct(vp1, vp1);
	float d20 = Vec3Df::dotProduct(vp2, vp0);
	float d21 = Vec3Df::dotProduct(vp2, vp1);
	float denom = d00 * d11 - d01 * d01;
	float v = (d11 * d20 - d01 * d21) / denom;
	float w = (d00 * d21 - d01 * d20) / denom;
	float u = 1.0f - v - w;
	Vec3Df TriangleNormal = NormalA*u + NormalB*v + NormalC*w;

	// default lighting in all parts that even are in shadow everywhere.
	// Maarten - This should be Ka, ambient light, but our .mtl files have Ka = (0,0,0).
	resCol = material.Kd()*backgroundlighting;

	float angle, distanceToLight;
	Vec3Df lightToIntersect, viewToIntersect, halfwayVector;

	for (Vec3Df v : MyLightPositions) {
		lightToIntersect = hitPoint - v;
		distanceToLight = lightToIntersect.normalize();

		viewToIntersect = hitPoint - origin;

		halfwayVector = (lightToIntersect + viewToIntersect);
		halfwayVector.normalize();

		/*
		std::cout << "Normal: " << TriangleNormal << std::endl 
			<< " diffuse angle: " << Vec3Df::cosAngle(TriangleNormal, lightToIntersect) << std::endl 
			<< " specular angle:  " << Vec3Df::cosAngle(TriangleNormal, halfwayVector) << std::endl;
		*/
		// Diffuse lighting
		angle = -Vec3Df::cosAngle(TriangleNormal, lightToIntersect);

		if (angle > 0)
			resCol += material.Kd()*angle*diffusePower / distanceToLight;
		
		angle = -Vec3Df::cosAngle(TriangleNormal, halfwayVector);
		if (angle > 0)
			resCol += material.Ks()*std::pow(angle, specularHardness);
			//resCol += Vec3Df(1.0f, 1.0f, 1.0f)*std::pow(angle, specularHardness);
	}

	if (material.illum() == 3) {

		// reflectdir = dir_of_ray - 2 * ray_projected_on_normal
		const Vec3Df reflectdir = dir - 2 * Vec3Df::dotProduct(TriangleNormal, dir) * TriangleNormal,
			newOrigin = hitPoint + 0.001 * reflectdir,
			newDest = hitPoint + reflectdir;

		if (debug)
			testRay.push_back(TestRay(newOrigin, Vec3Df(), Vec3Df()));

		resCol = 0.5*resCol + 0.5*performRayTracing(newOrigin, newDest, ++k);

		if (debug)
			testRay[k].color = resCol;
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
		testRay[0].color = performRayTracing(testRay[0].origin, testRay[0].destination, 0);


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

// Return the intersection point of the ray with the object in the scene
// return the rayDest if no object was hit
// Maarten - I think we don't need this function anymore, as this was used for debug, 
// but the debug stuff has been moved into performRayTracing via bool debug
/*inline Vec3Df calculateIntersectionPoint(const Vec3Df & rayOrigin, const Vec3Df & rayDest)
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
}*/

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
	float xMin, xMax, yMin, yMax, zMin, zMax;

	// Calculate the bounding box of the scene
	std::vector<Triangle>::const_iterator it = MyMesh.triangles.begin();
	xMin = (*it).v[0];
	xMax = (*it).v[0];
	yMin = (*it).v[1];
	yMax = (*it).v[1];
	zMin = (*it).v[2];
	zMax = (*it).v[2];
	++it;

	// Get the smallest and largest x, y and z values of the triangles of the mesh
	while (it != MyMesh.triangles.end()) {
		++it;
		// check x
		if ((*it).v[0] < xMin)
			xMin = (*it).v[0];
		else if ((*it).v[0] > xMax)
			xMax = (*it).v[0];

		// check y
		if ((*it).v[1] < yMin)
			yMin = (*it).v[1];
		else if ((*it).v[1] > yMax)
			yMax = (*it).v[1];

		// check z
		if ((*it).v[2] < zMin)
			zMin = (*it).v[2];
		else if ((*it).v[2] > zMax)
			zMax = (*it).v[2];
	}

	//TODO: split the scene into spaces which contain triangles.
	KDtreeCube scene = KDtreeCube(MyMesh.triangles, xMin, xMax, yMin, yMax, zMin, zMax);
	unsigned int minTriangles = (int) std::sqrtf( MyMesh.triangles.size());
	kdTree = splitSpace(scene, 0, minTriangles);
}

// recursively split spaces in half, till the spaces are small enough
std::vector<KDtreeCube> splitSpace(KDtreeCube cube, unsigned int axis, unsigned int minTriangles)
{
	std::vector<KDtreeCube> result;

	// don't divide subspaces with less than minTriangles triangles <--- !!! IMPORTANT: minTriangles is randomly chosen, needs testing !!!
	if (cube.triangles.size() < minTriangles)
	{
		result.push_back(cube);
		return result;
	}
	KDtreeCube left, right;
	// split the spacce in 2 smaller subspaces and recursively call this function on those subspaces
	if (axis % 3 == 0) {
		//split x-axis
		float xMid = (cube.xStart + cube.xEnd) / 2;

		// split the space in half
		left = KDtreeCube(cube.triangles, cube.xStart, xMid, cube.yStart, cube.yEnd, cube.zStart, cube.zEnd);
		right = KDtreeCube(cube.triangles, xMid, cube.xEnd, cube.yStart, cube.yEnd, cube.zStart, cube.zEnd);

	}
	else if (axis % 3 == 1) {
		//split y-axis
		float yMid = (cube.yStart + cube.yEnd) / 2;

		// split the space in half
		left = KDtreeCube(cube.triangles, cube.xStart, cube.xEnd, cube.yStart, yMid, cube.zStart, cube.zEnd);
		right = KDtreeCube(cube.triangles, cube.xStart, cube.xEnd, yMid, cube.yEnd, cube.zStart, cube.zEnd);

	}
	else if (axis % 3 == 2) {
		//split z-axis
		float zMid = (cube.yStart + cube.yEnd) / 2;

		// split the space in half
		left = KDtreeCube(cube.triangles, cube.xStart, cube.xEnd, cube.yStart, cube.yEnd, cube.zStart, zMid);
		right = KDtreeCube(cube.triangles, cube.xStart, cube.xEnd, cube.yStart, cube.yEnd, zMid, cube.zEnd);
	}

	// remove the triangles which are not inside of the cube
	removeTrianglesNotInSubSpace(left);
	removeTrianglesNotInSubSpace(right);

	// recursively calculate the KDtree subspaces of the two parts
	std::vector<KDtreeCube> leftres = splitSpace(left, axis + 1, minTriangles);
	std::vector<KDtreeCube> rightres = splitSpace(right, axis + 1, minTriangles);

	// concatenate the resulting vectors
	result.reserve(leftres.size() + rightres.size());
	result.insert(result.end(), leftres.begin(), leftres.end());
	result.insert(result.end(), rightres.begin(), rightres.end());

	
	// return the resulting vector
	return result;
}

// Check which of the triangles in space really lie within the space and delete the other traingles
void removeTrianglesNotInSubSpace(KDtreeCube cube)
{
	//TODO: find which triangles (partly) lie within the defined cube
}
