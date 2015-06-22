#ifdef WIN32
#include <windows.h>
#endif
#include "GL/glut.h"
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <thread>
#include <algorithm> //random_shuffle
#include <assert.h>
#include "config.h"
#include "raytracing.h"
#include "mesh.h"
#include "traqueboule.h"
#include "imageWriter.h"
#include "settings.h"



//This is the main application
//Most of the code in here, does not need to be modified.
//It is enough to take a look at the function "drawFrame",
//in case you want to provide your own different drawing functions



Vec3Df MyCameraPosition;

//MyLightPositions stores all the light positions to use
//for the ray tracing. Please notice, the light that is 
//used for the real-time rendering is NOT one of these, 
//but following the camera instead.
std::vector<Vec3Df> MyLightPositions;

//Main mesh 
Mesh MyMesh; 




class RayTracer {
public:
	RayTracer(void) {
		Image result(WindowSize_X, WindowSize_Y);

		for (int i = 0; i < 16; ++i)
			pixelsdone[i] = 0;

		for (int i = 0; i < WindowSize_Y; ++i)
			linedone[i] = false;

		produceRay(0, 0, &origin00, &dest00);
		produceRay(0, WindowSize_Y - 1, &origin01, &dest01);
		produceRay(WindowSize_X - 1, 0, &origin10, &dest10);
		produceRay(WindowSize_X - 1, WindowSize_Y - 1, &origin11, &dest11);
	};
	Vec3Df raytrace(double x, double y, float &depth);
	static void produceRay(int x_I, int y_I, Vec3Df * origin, Vec3Df * dest);
	void threadmethod(int threadID);
	double doDaRayTracingShizz();
private:
	clock_t begin;
	int pixelsdone[16];
	std::thread t[16];
	unsigned int tstart[16];
	unsigned int tend[16];
	unsigned int tcurrent[16];
	bool linedone[3000];

	//Setup an image with the size of the current image.
	Image result;

	//produce the rays for each pixel, by first computing
	//the rays for the corners of the frustum.
	Vec3Df origin00, dest00;
	Vec3Df origin01, dest01;
	Vec3Df origin10, dest10;
	Vec3Df origin11, dest11;
};


/**
 * Main function, which is drawing an image (frame) on the screen
*/
void drawFrame( )
{
	yourDebugDraw();
}

//animation is called for every image on the screen once
void animate()
{
	MyCameraPosition=getCameraPosition();
	glutPostRedisplay();
}



void display(void);
void reshape(int w, int h);
void keyboard(unsigned char key, int x, int y);

/**
 * Main Programme
 */
int main(int argc, char** argv)
{	
    glutInit(&argc, argv);

    //framebuffer setup
    glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );

    // positioning and size of window
    glutInitWindowPosition(WindowPos_X, WindowPos_Y);
    glutInitWindowSize(WindowSize_X,WindowSize_Y);
    glutCreateWindow(argv[0]);	

    //initialize viewpoint
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0,0,-4);
    tbInitTransform();     // This is for the trackball, please ignore
    tbHelp();             // idem
	MyCameraPosition=getCameraPosition();

	//activate the light following the camera
    glEnable( GL_LIGHTING );
    glEnable( GL_LIGHT0 );
    glEnable(GL_COLOR_MATERIAL);
    int LightPos[4] = {0,0,2,0};
    int MatSpec [4] = {1,1,1,1};
    glLightiv(GL_LIGHT0,GL_POSITION,LightPos);

	//normals will be normalized in the graphics pipeline
	glEnable(GL_NORMALIZE);
    //clear color of the background is black.
	glClearColor (backgroundColor[0], backgroundColor[1], backgroundColor[2], 0.0);

	
	// Activate rendering modes
    //activate depth test
	glEnable( GL_DEPTH_TEST ); 
    //draw front-facing triangles filled
	//and back-facing triangles as wires
    glPolygonMode(GL_FRONT,GL_FILL);
    glPolygonMode(GL_BACK,GL_LINE);
    //interpolate vertex colors over the triangles
	glShadeModel(GL_SMOOTH);


	// glut setup... to ignore
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutDisplayFunc(display);
    glutMouseFunc(tbMouseFunc);    // trackball
    glutMotionFunc(tbMotionFunc);  // uses mouse
    glutIdleFunc( animate);
	

	init();

	

    
	//main loop for glut... this just runs your application
    glutMainLoop();
        
    return 0;  // execution never reaches this point
}











/**
 * OpenGL setup - functions do not need to be changed! 
 * you can SKIP AHEAD TO THE KEYBOARD FUNCTION
 */
//what to do before drawing an image
 void display(void)
{
	glPushAttrib(GL_ALL_ATTRIB_BITS);//store GL state
    // Effacer tout
    glClear( GL_COLOR_BUFFER_BIT  | GL_DEPTH_BUFFER_BIT); // clear image
    
    glLoadIdentity();  

    tbVisuTransform(); // init trackball

    drawFrame( );    //actually draw

    glutSwapBuffers();//glut internal switch
	glPopAttrib();//return to old GL state
}
//Window changes size
void reshape(int w, int h)
{
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //glOrtho (-1.1, 1.1, -1.1,1.1, -1000.0, 1000.0);
    gluPerspective (50, (float)w/h, 0.01, 10);
    glMatrixMode(GL_MODELVIEW);
}













// react to keyboard input
void keyboard(unsigned char key, int x, int y)
{
	double times = 0;
    printf("key %d '%c' pressed at %d,%d\n", key, key, x, y);
    fflush(stdout);
    switch (key)
    {
	//add/update a light based on the camera position.
	case 'L':
		MyLightPositions.push_back(getCameraPosition());
		break;
	case 'l':
		MyLightPositions[MyLightPositions.size()-1]=getCameraPosition();
		break;
	case 'r':
	{
		//Pressing r will launch the raytracing.
		RayTracer r;
		r.doDaRayTracingShizz();
		break;
	}
	case 'R':
		for (int i = 0; i < 5; ++i) {
			RayTracer r;
			times += r.doDaRayTracingShizz() / 5;
		}
		printf("\nAverage RayTracing time: %f", times);
		break;
	case 27:     // touche ESC
        exit(0);
    }

	
	//produce the ray for the current mouse position
	Vec3Df testRayOrigin, testRayDestination;
	RayTracer::produceRay(x, y, &testRayOrigin, &testRayDestination);

	yourKeyboardFunc(key,x,y, testRayOrigin, testRayDestination);
}

/* @return amount of seconds the raytracing took */
double RayTracer::doDaRayTracingShizz() {
	printf("\nRaytracing a complete %i x %i image \n", WindowSize_X, WindowSize_Y);
	// Without the following line, changing sizes will mess up the preview.
	result.writeImageBMP("result.bmp");

	begin = clock();

	unsigned n = Thread_Amount;
	n = (n < 1 ? 1 : n); //not less than  1, duh. xD
	n = (n > 16 ? 16 : n); //not more than 16, as thread array is hardcoded length 16
	printf("There are %d threads", n);
	float linesperthread = float(WindowSize_Y) / float(n);

	for (int j = 0; j < n; ++j) {
		tcurrent[j] = tstart[j] = j * linesperthread;
		tend[j] = (j + 1) * linesperthread;
		t[j] = std::thread(&RayTracer::threadmethod, this, j);
	}
	// make sure the last thread won't get out of bounds
	tend[n - 1] = WindowSize_Y;

	int currpx;
	int totalsize = WindowSize_X * WindowSize_Y;
	do {
		//print progress once per second
		Sleep(1000);
		currpx = 0;
		for (int j = 0; j < n; ++j)
			currpx += pixelsdone[j];

		clock_t now = clock();
		double elapsed   = double(now - begin);
		double totaltime = double(totalsize) / double(currpx) * elapsed;
		double remaining = (totaltime - elapsed) / CLOCKS_PER_SEC;
		int r = (int)round(remaining);

		printf("\n%3d%%\tPixel %8d/%8d  %5d s remaining", 100 * currpx / totalsize, currpx, totalsize, r);

	} while (currpx < result._width * result._height);
	
	for (int j = 0; j < n; ++j)
		t[j].join();

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	printf("\nRaytracing comlete in %g seconds, storing result...\n", elapsed_secs);
	result.writeImagePPM("result.ppm");
	printf("Stored result in result.ppm");
	result.writeImageBMP("result.bmp");
	printf(" and result.bmp\n");
	return elapsed_secs;
}

void depthIncrease(float &depth, float &depthHolder, int &hitNumber) {
	if (depthHolder != std::numeric_limits<float>::max()) {
		depth += depthHolder;
		++hitNumber;
	}
}


void RayTracer::threadmethod(int threadID)
{
	bool done = false;
	unsigned y = tcurrent[threadID];
	while (!done)
	{
		// Perform raytracing on the line
		for (unsigned int x = 0; x < WindowSize_X; ++x)
		{
			float depth, depthHolder;
			int hitNumber = 0;
			Vec3Df rgb = raytrace(x, y);
			if (rgb != backgroundColor) {
				// MSAA is done using a rotated (square) grid,
				// and can therefore only be 4x or 16x.
				switch (MSAA) {
				case 4:
					rgb  = raytrace(x - 1. / 8., y - 3. / 8., depthHolder);	// +#++
					depthIncrease(depth, depthHolder, hitNumber);
					rgb += raytrace(x + 3. / 8., y - 1. / 8., depthHolder);	// +++#
					depthIncrease(depth, depthHolder, hitNumber);
					rgb += raytrace(x - 3. / 8., y + 1. / 8., depthHolder);	// #+++
					depthIncrease(depth, depthHolder, hitNumber);
					rgb += raytrace(x + 1. / 8., y + 3. / 8., depthHolder);	// ++#+
					depthIncrease(depth, depthHolder, hitNumber);
					rgb /= 4;
					if (hitNumber > 0)
						depth /= hitNumber;
					else
						depth = std::numeric_limits<float>::max();
					break;
				case 16:
					rgb  = raytrace(x -  9. / 32., y - 15. / 32., depthHolder);	// +++#++++ ++++++++
					depthIncrease(depth, depthHolder, hitNumber);
					rgb += raytrace(x -  1. / 32., y - 13. / 32., depthHolder);	// +++++++# ++++++++
					depthIncrease(depth, depthHolder, hitNumber);
					rgb += raytrace(x +  7. / 32., y - 11. / 32., depthHolder);	// ++++++++ +++#++++
					depthIncrease(depth, depthHolder, hitNumber);
					rgb += raytrace(x + 15. / 32., y -  9. / 32., depthHolder);	// ++++++++ +++++++#
					depthIncrease(depth, depthHolder, hitNumber);
					rgb += raytrace(x - 11. / 32., y -  7. / 32., depthHolder);	// ++#+++++ ++++++++
					depthIncrease(depth, depthHolder, hitNumber);
					rgb += raytrace(x -  3. / 32., y -  5. / 32., depthHolder);	// ++++++#+ ++++++++
					depthIncrease(depth, depthHolder, hitNumber);
					rgb += raytrace(x +  5. / 32., y -  3. / 32., depthHolder);	// ++++++++ ++#+++++
					depthIncrease(depth, depthHolder, hitNumber);
					rgb += raytrace(x + 13. / 32., y -  1. / 32., depthHolder);	// ++++++++ ++++++#+
					depthIncrease(depth, depthHolder, hitNumber);

					rgb += raytrace(x - 13. / 32., y +  1. / 32., depthHolder);	// +#++++++ ++++++++
					depthIncrease(depth, depthHolder, hitNumber);
					rgb += raytrace(x -  5. / 32., y +  3. / 32., depthHolder);	// +++++#++ ++++++++
					depthIncrease(depth, depthHolder, hitNumber);
					rgb += raytrace(x +  3. / 32., y +  5. / 32., depthHolder);	// ++++++++ +#++++++
					depthIncrease(depth, depthHolder, hitNumber);
					rgb += raytrace(x + 11. / 32., y +  7. / 32., depthHolder);	// ++++++++ +++++#++
					depthIncrease(depth, depthHolder, hitNumber);
					rgb += raytrace(x - 15. / 32., y +  9. / 32., depthHolder);	// #+++++++ ++++++++
					depthIncrease(depth, depthHolder, hitNumber);
					rgb += raytrace(x -  7. / 32., y + 11. / 32., depthHolder);	// ++++#+++ ++++++++
					depthIncrease(depth, depthHolder, hitNumber);
					rgb += raytrace(x +  1. / 32., y + 13. / 32., depthHolder);	// ++++++++ #+++++++
					depthIncrease(depth, depthHolder, hitNumber);
					rgb += raytrace(x +  9. / 32., y + 15. / 32., depthHolder);	// ++++++++ ++++#+++
					depthIncrease(depth, depthHolder, hitNumber);
					rgb /= 16.;
					if (hitNumber > 0)
						depth /= hitNumber;
					else
						depth = std::numeric_limits<float>::max();
					break;
				}
			}

			//store the result in an image 
			result.setPixel(x, y, RGBValue(rgb[0], rgb[1], rgb[2]), depth);
			++pixelsdone[threadID];
		}
		result.writeImageBMP("result.bmp", 0, y, result._width, 1);

		// search for next line to do
		if (!linedone[tcurrent[threadID] + 1]) {
			linedone[++tcurrent[threadID]] = true;
			++y;
		} else {
			done = true;

			std::vector<int> randomlist;
			for (int i = 0; i < Thread_Amount; ++i)
				randomlist.push_back(i); // fill list with threadIDs
			std::random_shuffle(randomlist.begin(), randomlist.end());

			for (int i = 0; i < randomlist.size(); ++i) {
				int tid = randomlist[i];
				if (done && tcurrent[tid]+1 < tend[tid] && !linedone[tend[tid]-1]) {
					linedone[tend[tid]-1] = true;
					y = tend[tid]-1;
					--tend[tid];
					done = false;
				}
			}
		}
	}
	printf("   %d done", threadID + 1);
}
//perform the raytracing for one x/y position
Vec3Df RayTracer::raytrace(double x, double y, float &depth)
{
	Vec3Df origin, dest;
	//produce the rays for each pixel, by interpolating 
	//the four rays of the frustum corners.
	double xscale = 1.0 - x / (double)(WindowSize_X - 1);
	double yscale = 1.0 - y / (double)(WindowSize_Y - 1);

	origin = yscale*(xscale*origin00 + (1 - xscale)*origin10) +
		(1 - yscale)*(xscale*origin01 + (1 - xscale)*origin11);
	dest = yscale*(xscale*dest00 + (1 - xscale)*dest10) +
		(1 - yscale)*(xscale*dest01 + (1 - xscale)*dest11);

	//launch raytracing for the given ray.
	return performRayTracing(origin, dest, 0, depth);
}

//transform the x, y position on the screen into the corresponding 3D world position
void RayTracer::produceRay(int x_I, int y_I, Vec3Df * origin, Vec3Df * dest)
{
	int viewport[4];
	double modelview[16];
	double projection[16];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview); //recuperer matrices
	glGetDoublev(GL_PROJECTION_MATRIX, projection); //recuperer matrices
	glGetIntegerv(GL_VIEWPORT, viewport);//viewport
	int y_new = viewport[3] - y_I;

	double x, y, z;

	gluUnProject(x_I, y_new, 0, modelview, projection, viewport, &x, &y, &z);
	origin->p[0] = float(x);
	origin->p[1] = float(y);
	origin->p[2] = float(z);
	gluUnProject(x_I, y_new, 1, modelview, projection, viewport, &x, &y, &z);
	dest->p[0] = float(x);
	dest->p[1] = float(y);
	dest->p[2] = float(z);
}
