# dla-nd

An arbitrary-dimensional diffusion-limited aggregation simulator.
See some images and more at (http://markjstock.org/dla-nd/)


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


## Building the software
In the directory that this file is in, just type

    make

and several binary executables will be made. If you want to build a 
specific executable for 1D through 5D, just type

    make dla-5d

replacing the 5 with 1, 2, 3, or 4, when appropriate. If you want to
build an executable for higher dimensions, you'll need to edit the 
`structs.h` file and add your own code.


## Usage

Edit the `sample.dnd` file, making sure to match all coordinate input with
the number of dimensions in the problem.

Unix users can then run the desired program as shown:

    dla-2d sample.dnd
    dla-2d < sample.dnd

Windows users have several options. You can drag the .dnd file onto the
appropriate .exe file, whereupon the output files will either appear in
the current working directory, or your home directory (which in my
case was C:\Documents and Settings\mstock).

You can also run "cmd" to enter the DOS-like command prompt. From there,
you can run 

    dla-2d.exe sample.dnd

Lastly, you can download and install MSYS, a Unix-like environment,
and run the program from within MSYS using

    dla-2d.exe sample.dnd
    dla-2d.exe < sample.dnd

If your .dnd file is set up to read in another file, you need to make
sure that the other file is in either (or both) your current working 
directory or your home directory. As a windows novice, this confused me.


## Advanced use

Feel free to monkey with the code. If you look in the `Makefile`, and the 
first part of `structs.h`, you'll quickly see how to make executables that
can track ANY number of dimensions. So, if you want to run a 27-dimensional
DLA, all you need is a large-enough computer.


## Conclusion

Thanks for taking the time to check this out. If you are having fun, I'd 
love to hear from you.

2005-09-14


## Change Log

v1.1	Chiral growth

v1.0	Windows and Linux executables included

v0.9    includes timer, line-at-a-time statistics file, constantly-
                overwritten pnd output file, random seed
        much faster due to elimination of slow and unnecessary routine

v0.8	fixed bulk motion algorithm; now uses numeric root-finder to
		determine actual dt for diffusing particles
	added experimental rotex growth using source/sink/vortex influences

v0.7	added libpng support, removed GIF support and related system calls
	added bit depth parameter to output file
	merged image_size and density pixel sizes---all output images use
		image_size now
	fixed bounds on sim-> construction variables
	fixed moving_particle problem in replace_all_particles()
	allowed automatic resizing of tree to accommodate all particles

v0.6	upgraded dla3d to dla-nd
	other minor improvements
