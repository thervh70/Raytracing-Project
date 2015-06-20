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
	MyLightPositions.push_back(*(new Vec3Df(0, 10, 0)));

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
	float a, b;

	for (int i = 0; i < triangles.size(); ++i)
	{
		hitpair = checkHit(triangles[i], origin, dest, minT);

		if (!hitpair.bHit)
			continue;

		hit = true;

		a = hitpair.res[0];
		b = hitpair.res[1];
		minT = hitpair.res[2];
		triangleIndex = i;
		material = MyMesh.materials[MyMesh.triangleMaterials[i]];
		triangle = triangles[i];
		hitPoint = hitpair.hitPoint;
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

		halfwayVector = (lightDir + viewDir);
		halfwayVector.normalize();

		// Diffuse lighting
		angle = Vec3Df::dotProduct(interpolatedNormal, lightDir); // cos(phi) = a . b / 1 / 1 (vectors are normalized)
		if (angle > 0)
			resCol += material.Kd()*angle * diffusePower; // / distanceToLight;

		// Shadows
		if (Vec3Df::dotProduct(lightDir, interpolatedNormal) <= 0.2) {
			resCol -= shadowRGB;
		} else {
			// Specular lighting
			angle = Vec3Df::dotProduct(interpolatedNormal, halfwayVector);
			if (angle > 0)
				resCol += material.Ks()*std::pow(angle, specularHighlight);

			for (Triangle t : triangles) {
				shadowHit = checkHit(t, (hitPoint + 0.001 * lightDir), v, std::numeric_limits<float>::max());
				if (shadowHit.bHit) {
					resCol -= shadowRGB;
				}
			}
		}
	}

	if (material.illum() == 3) {

		// reflectdir = dir_of_ray - 2 * ray_projected_on_normal
		const Vec3Df reflectdir = dir - 2 * Vec3Df::dotProduct(interpolatedNormal, dir) * interpolatedNormal,
			newOrigin = hitPoint + 0.001f * reflectdir,
			newDest = hitPoint + reflectdir;

		if (debug)
			testRay.push_back(TestRay(hitPoint, Vec3Df(), Vec3Df()));

		Vec3Df newCol = performRayTracing(newOrigin, newDest, ++k);
		resCol = 0.5*resCol + 0.5*newCol;

		if (debug)
			testRay[k].color = newCol;
	}

	// Normalized viewing vector, 
	// which will be considered as the vector of incidence.
	Vec3Df vIncidence = dir;
	vIncidence.normalize();

	// These are the indices of refraction.
	// n1 first material, n2 second material.
	float n1, n2;

	// If the refraction index of a material is greater than 1, 
	// then refraction has to be taken into account.
	if (material.Ni() > 1.0f) {
		// Transmission Ray calculation using a simplified version of the gigantic formula given in the slides:
		// (1); tRay = n1/n2 * vIncidence + ( n1/n2 * cos(thetaIncidence) - sqrt(1 - sin^2(thetaTransmitted) ) ) * normal.
		// with (2); sin^2(thetaTransmitted) = (n1/n2)^2 * (1 - cos^2(thetaIncidence))

		// Refractionindex is based on the division of the two refractive indices 
		// of the involved materials.
		float refractIndex = n1 / n2;

		// Cos(Theta) with theta as the angle of incidence calculation, by
		// calculating the dotProduct of the vector of indence and the normal of the hitPoint.
		float cosThetaIncidence = Vec3Df::dotProduct(vIncidence, interpolatedNormal);

		// Sin^2(Theta) calculation, with theta as the angle of incidence, by using formula 2.
		float sin2ThetaTransmitted = pow(refractIndex, 2) * (1 - (pow(cosThetaIncidence, 2)));

		// The result is a transmitted ray, by using formula 1.
		Vec3Df tRay = refractIndex * vIncidence + (refractIndex * cosThetaIncidence - sqrt(1 - sin2ThetaTransmitted)) * interpolatedNormal;

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
