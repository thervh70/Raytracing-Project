#ifdef WIN32
#include <windows.h>
#endif
#include "GL/glut.h"
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <thread>
#include <assert.h>
#include "raytracing.h"
#include "mesh.h"
#include "traqueboule.h"
#include "imageWriter.h"


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

unsigned int WindowSize_X = 200;  // resolution X
unsigned int WindowSize_Y = 200;  // resolution Y




class RayTracer {
public:
	RayTracer(void) {
		Image result(WindowSize_X, WindowSize_Y);

		begin = clock();
		prevsec = begin;

		Vec3Df a, b, c, d, e, f, g, h;
		produceRay(0, 0, &a, &b);
		produceRay(0, WindowSize_Y - 1, &c, &d);
		produceRay(WindowSize_X - 1, 0, &e, &f);
		produceRay(WindowSize_X - 1, WindowSize_Y - 1, &g, &h);
		origin00 = a; dest00 = b;
		origin01 = c; dest01 = d;
		origin10 = e; dest10 = f;
		origin11 = g; dest11 = h;
	};
	static void produceRay(int x_I, int y_I, Vec3Df * origin, Vec3Df * dest);
	void threadmethod(unsigned int y);
	void doDaRayTracingShizz();
private:
	clock_t begin, prevsec;

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
    glutInitWindowPosition(200, 100);
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
	glClearColor (0.0, 0.0, 0.0, 0.0);

	
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
	case 27:     // touche ESC
        exit(0);
    }

	
	//produce the ray for the current mouse position
	Vec3Df testRayOrigin, testRayDestination;
	RayTracer::produceRay(x, y, &testRayOrigin, &testRayDestination);

	yourKeyboardFunc(key,x,y, testRayOrigin, testRayDestination);
}

void RayTracer::doDaRayTracingShizz() {
	printf("\nRaytracing a complete %i x %i image \n", WindowSize_X, WindowSize_Y);
	//std::cout << "origin00 " << origin00 << "    dest00 " << dest00 << std::endl;

	unsigned n = std::thread::hardware_concurrency();
	n = (n == 0 ? 2 : n);
	n = (n > 8 ? 8 : n);
	n *= 2; // More efficient than *1, apparently, or at least for me
	std::thread t[16]; // = 8 * 2
	printf("There are %d threads, %d are fired at a time \n", n/2, n);

	for (unsigned int y = 0; y < WindowSize_Y; y += n) {
		for (int j = 0; j < n; ++j)
			t[j] = std::thread(&RayTracer::threadmethod, this, y + j);
		
		for (int j = 0; j < n; ++j)
			t[j].join();
		
		//print progress once per second
		if (clock() - prevsec > CLOCKS_PER_SEC) {
			int currpx = y*WindowSize_Y;
			int s = WindowSize_X*WindowSize_Y;
			printf("%3d%%\tPixel %8d/%8d\n", 100 * currpx / s, currpx, s);
			prevsec += CLOCKS_PER_SEC;
		}
	}

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	printf("Raytracing comlete in %g seconds, storing result...\n", elapsed_secs);
	result.writeImagePPM("result.ppm");
	printf("Stored result in result.ppm");
	result.writeImageBMP("result.bmp");
	printf(" and result.bmp\n");
}

void RayTracer::threadmethod(unsigned int y) {
	//std::cout << "origin00 " << origin00 << "    dest00 " << dest00 << std::endl;
	for (unsigned int x = 0; x<WindowSize_X;++x)
	{
		Vec3Df origin, dest;
		//produce the rays for each pixel, by interpolating 
		//the four rays of the frustum corners.
		float xscale = 1.0f - float(x) / (WindowSize_X - 1);
		float yscale = 1.0f - float(y) / (WindowSize_Y - 1);

		origin = yscale*(xscale*origin00 + (1 - xscale)*origin10) +
			(1 - yscale)*(xscale*origin01 + (1 - xscale)*origin11);
		dest = yscale*(xscale*dest00 + (1 - xscale)*dest10) +
			(1 - yscale)*(xscale*dest01 + (1 - xscale)*dest11);

		//launch raytracing for the given ray.
		Vec3Df rgb = performRayTracing(origin, dest);
		//store the result in an image 
		result.setPixel(x, y, RGBValue(rgb[0], rgb[1], rgb[2]));
	}
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
