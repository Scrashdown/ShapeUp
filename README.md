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

    $ git submodule update --init --recursive

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

## Sample Scene
In this section we show how one can recreate a melting candle scene.
Start by running the application using the following command in the *simulation* folder:

	./simulation

Then one can load the candle object by clicking on *Open from disk* and selecting the candle.obj from ShapeUp/data/ folder.
At any moment it is possible to view the temperature of the surface vertices by toggling *Temperature*. We want our mesh to be tetrahedralized so it is necessary to toggle the *Tetrahedralize* button.
Start by adding diffusion to the mesh. Clicking on *Temp. Diffuse* and choosing the appropriate diffusion and number of iterations. 
The next step is to add constraints to our simulation. Start by adding a *Floor* to the scene. Then adding *Temperature Elasticity* and *Triangle Bending* constraints, these do not require vertices to be selected as they will be applied to every vertex by default. By default these constraints will have sliders set to 1, these sliders can be modified at the user's discretion in order to achieve the optimal behaviour. The default values will do fine for testing purposes. At this point 3 constraints should be showing up on the top right panel. Additional constraints can be added such as positional constraints if some vertices are required to be stationnary.

By running the simulation the candle will fall to the floor and stabilize. At this point user can modify vertex temperatures using the panel on the left, under *Set Temperature*. First select the vertices to be modified and then change the temperature value. When everything is ready for simulation, clicking on *Simulate* will start the simulation. If temperatures are visible it is possible to notice the diffusion and consequent modification of the mesh's geometry - melting -. Try to select the bottom vertices and making them very hot (0K is cold and 2000K is extremely hot) to see what happens!

Plasticity is achieved after a certain threshold and thus spring constraints resting lengths will not reset to original values. This is expected. Resetting positions after this event may cause undefined behaviour therefore we advise users to restart the program if they wish to reset constraints and positions after melting.

Enjoy our work as much as we have enjoyed working on it! 



