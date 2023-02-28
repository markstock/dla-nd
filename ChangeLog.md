# dla-nd
An arbitrary-dimensional diffusion-limited aggregation simulator.

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

