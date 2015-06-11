#include <stdio.h>
#ifdef WIN32
#include <windows.h>
#endif
#include "GL/glut.h"
#include "raytracing.h"
#include "config.h"
#include "Vec3D.h"
#include "Matrix33.h"


//temporary variables
//these are only used to illustrate 
//a simple debug drawing. A ray 
Vec3Df testRayOrigin;
Vec3Df testRayDestination;


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
	MyLightPositions.push_back(*(new Vec3Df(100, 42, 17)));

	/* FOR TESTING ONLY ~ Maarten
	Matrix33f m(Vec3Df(3, 4, 9), Vec3Df(5, 12, 8), Vec3Df(9, 3, 1));
	std::cout << m << " det:" << m.det() << std::endl;
	std::cout << "X = " << m.solve(Vec3Df(1, 2, 3)) << std::endl;
	system("pause");*/

	// FOR TESTING ONLY ~ Mathias
	/*testRayOrigin = Vec3Df(2.0f, 5.0f, 1.0f);
	testRayDestination = Vec3Df(8.0f, 7.0f, 4.0f);
	performRayTracing(testRayOrigin, testRayDestination);*/
}

//return the color of your pixel.
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest, int k)
{
	Vec3Df resCol, point, n, dir = dest - origin;
	Hitpair hitpair;
	float t, D, minT = std::numeric_limits<float>::max();
	std::vector<Triangle> triangles = MyMesh.triangles;
	Material material;

	for (int i = 0; i < triangles.size(); ++i)
	{
		hitpair = checkHit(triangles[i], origin, dest, minT);

		if (!hitpair.bHit)
			continue;

		minT = hitpair.res[2];

		material = MyMesh.materials[MyMesh.triangleMaterials[i]];
	}

	Vec3Df intersectionPoint = origin + minT * (dest - origin);
	resCol = Vec3Df(0.0f, 0.0f, 0.0f);
	int angle;
	for (Vec3Df v : MyLightPositions) {
		angle = Vec3Df::cosAngle(intersectionPoint-v, dest - origin);

		if(angle != -2147483648 && angle != 0)
			std::cout << angle << std::endl;

		if (angle > 0) {
			resCol = material.Kd()*angle/ MyLightPositions.size();
		}
	}
//WIP
/*	n = 0; // Interpolate over vertices
	point = dir - 2 * Vec3Df::dotProduct(n, dir) * n;
	if (material.has_illum()) {
		int illum = material.illum();
		if (illum == 3) {

			performRayTracing(hitpair.hitPoint, point, k++);
		}
	}*/


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
	glColor3f(0,1,1);
	glVertex3f(testRayOrigin[0], testRayOrigin[1], testRayOrigin[2]);
	glColor3f(0,0,1);
	glVertex3f(testRayDestination[0], testRayDestination[1], testRayDestination[2]);
	glEnd();
	glPointSize(10);
	glBegin(GL_POINTS);
	glVertex3fv(MyLightPositions[0].pointer());
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
		//here, as an example, I use the ray to fill in the values for my upper global ray variable
		//I use these variables in the debugDraw function to draw the corresponding ray.
		//try it: Press a key, move the camera, see the ray that was launched as a line.
		testRayOrigin = rayOrigin;

		// calculate the coordinates where the ray will hit an object
		// and use those coordinates as ray destination
		testRayDestination = calculateIntersectionPoint(rayOrigin, rayDestination);

		std::cout << "Origin      " << testRayOrigin << std::endl;
		std::cout << "Destination " << testRayDestination << std::endl;
		std::cout << "Color       " << performRayTracing(testRayOrigin, testRayDestination, 0) << std::endl;
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
