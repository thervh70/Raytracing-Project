# Bullet list of extra implemented features for the Raytracing-Project

## Group 6

![Final image for the competition](http://i.imgur.com/eQAUTMn.jpg)

### Features (chronologically ordered)

 - Diffuse light
 - Multi-threading
 - Preview in BMP-format
 - Reflection
 - Blinn-Phong Shading
 - Specular highlights
 - Shadows
 - Multiple light sources
 - KD-Tree
 - Refraction
 - MSAA
 - Textures
 - Depth of Field
 - Smooth transition between reflection and refraction

### Debug features
 - Raytracing
    - 'r' to launch the multi-threaded RayTracer, rendering an image. 
    - 'R' to launch the RayTracer 5 times and calculate the average time of these runs.
 - Debug Ray
    - ' ' to shoot a debug ray through the scene, console output is plenty
     Will also draw all lines reflected and refracted in the scene.
     Known bug: first ray shot will always return a white color (1.0, 1.0, 1.0)
 - Lighting
    - 'L' to add a light source at the camera position
    - 'l' to move the last added light source to the camera position
    - 'c' to clear all light sources and add a light source at camera
 - Focal point
    - 'f' to set the focal point to the depth of the object the mouse is pointed at
(our final render is made without setting the focal point, as this was still a bit bugged)
 - Opening result
    - 'o' to open the result of the raytracing, 'result.bmp'.
     If opened before rendering, Windows Photo Viewer will refresh it when progress is written to the file.
    - 'O' to open the result of the Depth of Field post-processing, 'blurred.bmp'.
 - Save lights and camera position
    - 's' to save the current state of the camera and the light positions.
    - 'S' to load these settings back in the preview (warning: crashes if you haven't saved before)
 - Exit
    - 'Esc' to exit the program.


Special thanks goes out to hjmediastudios for the permission to use and adapt the Lego Soldier Object.
Source: http://www.blendswap.com/blends/view/6437

We did not use any code from the internet, except for an image reading library called Corona (library included with source-code). We did use formulas from Wikipedia, the course slides, and other Computer Graphics related websites.

Copyright © 2015 TI1806 - Group 6 (see LICENSE.md)
