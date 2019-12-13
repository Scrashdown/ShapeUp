# Projective Dynamics and ShapeUp Framework

## What it is

A minimal framework that implements ShapeUp style mesh editing.
The GUI is implemented using **nanogui** and **OpenGL** and ShapeUp is implemented using **Eigen** and **OpenMP**.
The mesh is intermediately loaded as a **Surface_mesh** object, but converted to simple vertex positions / triangle indices matrices for ShapeUp internally.

## Build
### Unix
This project uses Wenzel Jakob's [nanogui](https://github.com/wjakob/nanogui), which is included under externals,
but has additional dependencies, which can be resolved on Ubuntu/Debian, via:

    $ apt-get install cmake xorg-dev libglu1-mesa-dev

and on Redhead/Fedora, via:

	$ sudo dnf install cmake mesa-libGLU-devel libXi-devel libXcursor-devel libXinerama-devel libXrandr-devel xorg-x11-server-devel

This project also relies on the libigl library to work, imported as a submodule. It may be done with the following command if needs be:

    $ git submodule add https://github.com/libigl/libigl.git external/libigl

After that, simply build the project as usual:

	mkdir build
	cd build
	cmake ..
	make

The ShapeUp editor can then be run using

    ./shapeup/shapeup

### Windows/MVSC 19
	Run CMake
	Set the source code dir to the base dir
	Set the build dir to the base dir + "/build"
	Press Configure
	Press Generate
	Open the generated solution in the build directory
	Build the solution and run the shapeup project
    
## ShapeUp Controls
    Open a mesh.
    Vertices can be selected using ctrl + left-click + drag, which shows a selection box.
    Selected vertices (blue) can be made into position constraint groups using the button “Fix Selection”. Create two such groups on different places on the mesh.
    Vertices marked in green can be moved (as a group) by alt/option + left-click + drag, which changes the target positions of the constraints.
    The shape is maintained by using either edge springs or triangle strain (or, once implemented, Iso 1-Ring constraints), which can be selected using the corresponding radio button. Note how much nicer the shape is maintained using the triangle strain + bending constraints, see the next slides for details.
    Additional constraints can be added to selected vertices by the buttons to the lower left. For now, only smoothing constraints are available.
    The currently used constraints are shown as groups in the window to the top right. There, the constraint weights can be changed using the sliders. Hovering over a slider, shows the vertices on which these constraints are defined.

