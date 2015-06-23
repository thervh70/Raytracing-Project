#include <stdio.h>
#ifdef WIN32
#include <windows.h>
#endif
#include "GL/glut.h"
#include "raytracing.h"
#include "KDtree.h"
#include "config.h"
#include "Vec3D.h"
#include "Matrix33.h"
#include "settings.h"

// All the rays in testRay will be drawn in yourDebugDraw().
std::vector<TestRay> testRay;
bool debug = false;
float focusDepth;

// remember if the acceleration tree has been built
bool builtAccelTree = false;

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

Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest, int k, float prev_Ni) {
	float depth;
	return performRayTracing(origin, dest, k, prev_Ni, depth);
}

//return the color of your pixel.
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest, int k, float prev_Ni, float &minT)
{
	if (k > 5) return endOfReflection;

	Vec3Df resCol;
	Vec3Df dir = dest - origin;

	// The struct "hit" contains all information about the hitpoint
	HitTriangle hit = checkHit(origin, dest, minT);

	if (debug)
		testRay[k].destination = hit.bHit ? hit.hitPoint : dest;

	if (!hit.bHit) return backgroundColor;

	// Normals of three vectors of triangle
	Vec3Df
		normalA = MyMesh.vertices[hit.triangle.v[0]].n,
		normalB = MyMesh.vertices[hit.triangle.v[1]].n,
		normalC = MyMesh.vertices[hit.triangle.v[2]].n;

	// Barycentric vertex normal biliniear interpolation shading mode
	float u = hit.a, v = hit.b, w = 1.f - u - v;
	Vec3Df interpolatedNormal = normalA*u + normalB*v + normalC*w;
	interpolatedNormal.normalize();

	// default lighting in all parts that even are in shadow everywhere.
	// Maarten - This should be Ka, ambient light, but our .mtl files have Ka = (0,0,0).
	resCol = hit.material.Kd()*backgroundlighting;

	float angle, distanceToLight;
	Vec3Df lightDir, viewDir, halfwayVector;
	float specularHighlight = hit.material.Ns();

	// Full shadow reduces light by (0.9,0.9,0.9) when ShadowFactor = 0.9
	float shadowPart = ShadowFactor / MyLightPositions.size();
	Vec3Df shadowRGB = Vec3Df(1.f, 1.f, 1.f) * shadowPart;
	float shadowMax = 0.f;
	HitTriangle shadowHit;

	viewDir = origin - hit.hitPoint;
	viewDir.normalize();

	for (Vec3Df v : MyLightPositions) {
		lightDir = v - hit.hitPoint;
		distanceToLight = lightDir.normalize();

		// Diffuse lighting
		angle = Vec3Df::dotProduct(interpolatedNormal, lightDir); // cos(phi) = a . b / 1 / 1 (vectors are normalized)
		if (angle > 0)
			resCol += hit.material.Kd()*angle * diffusePower; // / distanceToLight;

		angle = Vec3Df::dotProduct(lightDir, interpolatedNormal);
		// Self Shadows (Dark side of object)
		if (angle <= 0.f) {
			shadowMax = 1.0f;
		} 
		else if (angle > 0.f & angle < 0.3f) {
			shadowMax = ((1.0f/0.3f) * (0.3f - angle));
		} 
		if (angle > 0.f) {
			// Specular lighting
			halfwayVector = (lightDir + viewDir);
			halfwayVector.normalize();
			angle = Vec3Df::dotProduct(interpolatedNormal, halfwayVector);
			if (angle > 0)
				resCol += hit.material.Ks()*std::pow(angle, specularHighlight);
			
			//Shadow casted by an object
			if (shadowSamples <= 1) {
				if (checkHit((hit.hitPoint + 0.001f * lightDir), v).bHit)
					shadowMax = 1.0f;
			}
			if (shadowSamples > 1) {
				float softShadow;
				// The + 0.000...0001f prevents the compiler saying "divide by zero"
				float delta = 2.f / (float)(shadowSamples - 1 + 0.000000000000000000001f);
				for (float x = -1.f; x <= 1.f; x += delta)
				for (float y = -1.f; y <= 1.f; y += delta)
				for (float z = -1.f; z <= 1.f; z += delta) {
					float rx = 0.f;//delta / 2.f * (double)rand() / (double)RAND_MAX - delta / 4.f;
					float ry = 0.f;//delta / 2.f * (double)rand() / (double)RAND_MAX - delta / 4.f;
					float rz = 0.f;//delta / 2.f * (double)rand() / (double)RAND_MAX - delta / 4.f;
					shadowHit = checkHit((hit.hitPoint + 0.5f * lightDir), v + Vec3Df((x + rx) * shadowRadius, (y + ry) * shadowRadius, (z + rz) * shadowRadius));
					if (shadowHit.bHit) {
						softShadow += 1.0f / shadowSamples / shadowSamples / shadowSamples;
					}
				}
				shadowMax = softShadow > shadowMax ? softShadow : shadowMax;
			}
		}
		resCol -= shadowMax * shadowRGB;
	}

	// These are the indices of refraction.
	// n1 is the material from where the ray comes from, which is the Ni value of the previous hitPoint,
	// if it has a Ni value, else float 1 (standard for air).
	float n1 = (prev_Ni != 0.0f) ? prev_Ni : 1.0f;
	// n2 is the material of the hitPoint, which is the Ni value if it has a Ni value, else float 1.
	float n2 = (hit.material.has_Ni()) ? hit.material.Ni() : 1.0f;

	if (hit.material.illum() == 3) {

		// reflectdir = dir_of_ray - 2 * ray_projected_on_normal
		const Vec3Df reflectdir = dir - 2 * Vec3Df::dotProduct(interpolatedNormal, dir) * interpolatedNormal,
			newOrigin = hit.hitPoint + 0.001f * reflectdir,
			newDest = hit.hitPoint + reflectdir;

		if (debug)
			testRay.push_back(TestRay(hit.hitPoint, Vec3Df(), Vec3Df()));

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
	if (hit.material.Ni() > 1.0f) {
		// Transmission Ray calculation using a simplified version of the gigantic formula given in the slides:
		// (1); tRay = n1/n2 * vIncidence + ( n1/n2 * cos(thetaIncidence) - sqrt(1 - sin^2(thetaTransmitted) ) ) * normal.
		// with (2); sin^2(thetaTransmitted) = (n1/n2)^2 * (1 - cos^2(thetaIncidence))

		// Refractionindex calculation is based on the current refractionindex and the previous one.
		// If they are not equal, then use the current one (a ray will now be shot into the object).
		// If they are equal, then calculate n2^-1 (a ray will now be shot from the object, in order to get out of it).
		float refractIndex = (n1 != n2) ? n2 : (1 / n2);

		// Cos(Theta) with theta as the angle of incidence calculation, by
		// calculating the dotProduct of the vector of indence and the normal of the hitPoint.
		float cosThetaIncidence = Vec3Df::dotProduct(vIncidence, interpolatedNormal);

		// Sin^2(Theta) calculation, with theta as the angle of transmittance, by using formula 2.
		float sin2ThetaTransmitted = pow(refractIndex, 2) * (1 - (pow(cosThetaIncidence, 2)));

		// Important: only shoot the tRay if the angle is smaller than the critical angle.
		if (sin2ThetaTransmitted <= 1) {
			// The result is a transmitted ray, by using formula 1.
			const Vec3Df tRay = refractIndex * vIncidence + (refractIndex * cosThetaIncidence - sqrt(1 - sin2ThetaTransmitted)) * interpolatedNormal,
				newOriginR = hit.hitPoint + 0.001f * tRay,
				newDestR = hit.hitPoint + tRay;

			// Debug and Tracing		
			if (debug)
				testRay.push_back(TestRay(hit.hitPoint, Vec3Df(), Vec3Df()));

			Vec3Df newCol = performRayTracing(newOriginR, newDestR, ++k, n2);
			resCol = 0.3*resCol + 0.7*newCol;

			if (debug)
				testRay[k].color = newCol;

		}
		else {
			// The result is just a reflected ray.
			const Vec3Df refRay = vIncidence - 2 * Vec3Df::dotProduct(vIncidence, interpolatedNormal) * interpolatedNormal,
				newOriginR = hit.hitPoint + 0.001f * refRay,
				newDestR = hit.hitPoint + refRay;

			// Debug and Tracing		
			if (debug)
				testRay.push_back(TestRay(hit.hitPoint, Vec3Df(), Vec3Df()));

			Vec3Df newCol = performRayTracing(newOriginR, newDestR, ++k, n2);
			resCol = 0.5*resCol + 0.5*newCol;

			if (debug)
				testRay[k].color = newCol;

		}
	}

	return resCol;
}

HitTriangle checkHit(const Vec3Df & origin, const Vec3Df & dest) {
	float d;
	return checkHit(origin, dest, d);
}

HitTriangle checkHit(const Vec3Df & origin, const Vec3Df & dest, float &minT) {
	HitTriangle res;
	Hitpair hitpair;
	Vec3Df dir = dest - origin;
	minT = std::numeric_limits<float>::max();
	std::vector<int> triangles;

	// get starting node
	AccelTreeNode currNode = findChildNode(treeRoot, 0, origin);
	AccelTreeNode hitNode = treeRoot;
	std::vector<AccelTreeNode*> oldParentList;
	bool gotHit = false;
	res.bHit = false;

	/* If current point is outside of the root (camera is far away), then put it inside the root.
	if (currentPoint.p[0] < treeRoot.xStart || currentPoint.p[0] > treeRoot.xEnd ||
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

				hitpair = checkHit(MyMesh.triangles[triangles[i]], origin, dest, minT);

				if (!hitpair.bHit)
					continue;
				hitNode = (*curPar[n]);
				res.bHit = gotHit = true;
				res.triangle = MyMesh.triangles[triangles[i]];
				res.material = MyMesh.materials[MyMesh.triangleMaterials[triangles[i]]];
				res.hitPoint = hitpair.hitPoint;
				res.a = hitpair.res[0];
				res.b = hitpair.res[1];

				minT = hitpair.res[2];
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

	return res;
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

	Vec3Df res(res2[0], res2[1], -res2[2]); // res[2] = -res[2].
											// This is done in this way because you can't modify vector entries in a thread apparently

											//Check if hit
											//res = [a, b, t]
	if (res[0] < 0 || res[1] < 0 || res[0] + res[1] > 1 || res[2] < 0 || res[2] > minT)
		result.bHit = false;
	else
		result.bHit = true;

	result.res = res;
	result.hitPoint = origin + res[2] * dir;
	return result;
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
	int res;
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
		testRay[0].color = performRayTracing(testRay[0].origin, testRay[0].destination, 0, 0.0f, focusDepth);

		std::cout << "DEBUG RAY TRACE" << std::endl;

		for (Vec3Df v : MyLightPositions) {
			std::cout << "Light position: " << v << std::endl;
		}
		
		std::cout << std::endl;
		for (TestRay r : testRay) {
			std::cout << "Origin        " << r.origin << std::endl;
			std::cout << "Destination   " << r.destination << std::endl;
			std::cout << "Color         " << r.color << std::endl;
			std::cout << "(Focus)-Depth " << focusDepth << std::endl << std::endl;
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
		res = system("result.bmp");
		std::cout << " exited with code " << res << std::endl;
		break;
	}
}
