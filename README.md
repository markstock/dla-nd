# dla-nd
An arbitrary-dimensional diffusion-limited aggregation simulator.

## Features
`dla-nd` supports the following features

* Automatic, persistent, adaptive 2^D-tree space subdivision
* Statistically-correct particle motion bias (bulk motion) numerically solves quartic equation each step
* Stubbornness (Nittman and Stanley's artificial surface tension model)
* Tunable stickiness and grippiness
* 2-d or 3-d density field construction
* Raw, Radiance, and generic OBJ output for particle positions and connectivities
* PNG and ASCII PGM output for particle locations and density field
* Restart from an old simulation, or merge multiple simulations using raw .part files 

See some images and more at [the old homepage](http://markjstock.org/dla-nd/).

## Compilation
This software uses the CMake build system, which is available on Windows, OSX, and Linux. The only external dependency is libpng, which is usually installed on most systems. If not, you can get both dependencies with (Red Hat/Fedora, Debian/Ubuntu, OSX/brew):

    sudo dnf install cmake libpng-devel
    sudo apt-get install cmake libpng-dev
    brew install cmake libpng

Once that's settled, from the directory with this README, run:

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release
    make

and several binary executables will be made. If you want to build a specific executable for 1D through 5D, just ask for it directly:

    make dla-3d

## Usage
You no longer need to use an input file to run dla-nd, *but* I haven't documented the options well. See the `setup.c` file for those. But to get started with some fun stuff, try:

    dla-2d -p 0 0 -every 1000 -end 100000 -zone -2 2 -2 2 -dens -res 1000

This will run a 2D sim which starts with a single particle at 0,0, writes a 1000x1000 pixel PNG image out every 1000 steps, in which the particles are rendered onto a density field that spans -2..2 in x and y.

    dla-3d -p 0 0 0 -every 10000 -end 100000 -zone -1 1 -1 1 -1 1 -slices -res 200 -seg

The above command will perform a similar job, but will write 200 PNG files every 10000 steps (slices) and also write an OBJ-like "seg" file containing the vector geometry. This `parts0001.seg` file can be manipulated with my [stickkit](https://github.com/markstock/stickkit) program.

Windows users can also run "cmd" to enter the DOS-like command prompt. From there, you can run any of the commands above.

## Options
To do: fill this out with descriptions of `-grip` and `-stick` and bulk motion, etc. Maybe some images, too.

## Change Log

v1.7	2018-12-09	2D and 3D density output supported, more command-line options supported

v1.6	2015-03-21	Enabling multiple types of bulk velocity

v1.5	2015-01-01	Smooth junctions each step

v1.4	2014-03-05	Multiple methods for segment rejection

v1.3	2008-10-18	Support file-input options on command-line

v1.2	2008-03-30	Support for "irradiation planes" and other enhancements to allow Lichtenberg figures

v1.1	2006-04-06	Chiral growth

v1.0	2005-09-14	Windows and Linux executables included

v0.9    2005-09-02	includes timer, line-at-a-time statistics file, constantly-
                overwritten pnd output file, random seed
        much faster due to elimination of slow and unnecessary routine

v0.8	2005-02-26	fixed bulk motion algorithm; now uses numeric root-finder to
		determine actual dt for diffusing particles
	added experimental rotex growth using source/sink/vortex influences

v0.7	2005-02-24	added libpng support, removed GIF support and related system calls
	added bit depth parameter to output file
	merged image_size and density pixel sizes---all output images use
		image_size now
	fixed bounds on sim-> construction variables
	fixed moving_particle problem in replace_all_particles()
	allowed automatic resizing of tree to accommodate all particles

v0.6	2004-08-31	upgraded dla3d to dla-nd
	other minor improvements

v0.5	2004-08-26	added grip, reorganized parameters, added obj output

v0.4	2004-02-27	corrected directionality, added bias, added stickiness

v0.3	2004-01-20	added directionality, reorganized

v0.2	2003-12-04	new output, correct aggregation

v0.1	2003-12-04	improved order over O(N^2)

v0.0	2003-12-03	copied part3d source, begin work

