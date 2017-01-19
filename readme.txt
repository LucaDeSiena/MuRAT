% PROGRAM MuRAT2.0 for Multi-resolution attenuation tomography analysis
%
% PROGRAM AIM: Multi-step attenuation tomography using the
% coda-normalization method on high frequency datasets.
%
% Author:  L. De Siena.
%           
%  March 2016
% 
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
Installation procedure:

Execute the main program from the command line 

>>MuRAT

The program works using the Matlab user interface and processes
files in ascii format.

List of Programs:

1. MuRAT: Main program with GUI.

2. function segments.m, which creates the data used to build the inversion
   matrix

3. commands.m: change to set your own parameters

%--------------------------------------------------------------------------
% Structure of the data (recording geometry and travel-times)
%--------------------------------------------------------------------------

The program requires the following files/instructions:

1. Update the file 'commands.m', with the names of the velocity models you
   are using and the values of the quantity (number of steps, central
   frequency, lapse time, etc.) you want to set. See also Point 7.

2. One or two velocity models parametrized in a regular constant Cartesian
   grids. The first three columns are the Cartesian coordinates, with
   negative depths. The fourth column is the value of the velocity.

%--------------------------------------------------------------------------   

   Example:
   The first coordinate to change must be Z, then Y and finally X:

             X     Y     Z     V    N
             x1    y1    z1    v1   1
             x1    y1    z2    v2   2
             .........................
             x1    y1    zm    vm   m
             x1    y2    z1    vm1  m+1
             x1    y2    z2    vm2  m+2
             .........................
             x1    yk    zm    vkm  km
             x2    y1    z1    vkm1 km+1
             .........................

   Two example files of velocity models
   ('modv1.txt' and 'modv2.txt') are inside the folder.

%--------------------------------------------------------------------------

3. The ray-files in meters or km. Each file must be propagated in the
   reference system of the velocity models. The first column counts
   the steps. The second and third columns are X and Y. The fourth column is the
   coordinate Z of the ray, and must be positive. It can include topography,
   if the ray is produced in a reference system with topography (see
   example). The first point must be the source, the last the station.

   The length of each ray-step (the difference in length between a
   point and its following) must be much less then the step of the velocity
   grid. The fourth column is the travel time between the two points.
   The program is set to work with positive ray-depth, while the reference
   system of the velocity has negative depths.
   The fuction segments automatically changes positive depths into negative. 

%--------------------------------------------------------------------------

   Example: if the step of the velocity model with topography
            is 500 m the one of the ray should be below 100:

Topography = 3350 = originz
Source = -600 = -0.6 km in the velocity model = 600 m depth
Station = 1000 = 1 km in the velocity model = -1000 m depth


            modv                       ray1
         
	   0   0   -1000           1 2050   1525   3950 <= source depth+originz 
         0   0   -500            2 2120   1515   3974
         0   0   -0              3 2122   1545   3988
         0   0   500             4 1954   1540   3999
         0   0   1000            5 1945   1515   3987
         ............            .....................
				      1980 1203   3450   2350 <= station depth+originz
   



   See also the example files for the rays in the folder example/rays.

%--------------------------------------------------------------------------

4. The program works on the vertical and/or horizontal
   registrations of three-component stations. Each ray
   corresponds to three seismic tracest.
   The order of the ray files must be the same of the traces,
   starting with the WE, then SN, and vertical.

%--------------------------------------------------------------------------

   Example: The code provides the lists, which must be checked.
		2 rays recorded at 1 three component station.
            list file for the rays: 6 rows
            list file for the traces: 6 rows (1 for each component).
          
        Rays-list file     Trace-list file
 Row1   ./rays/ray1        ./traces/ray1-WE
 Row2   ./rays/ray1        ./traces/ray1-SN
 Row3   ./rays/ray1        ./traces/ray1-V
 Row4   ./rays/ray2        ./traces/ray2-WE
 Row5   ./rays/ray2        ./traces/ray2-SN
 Row6   ./rays/ray2        ./traces/ray2-V

   
%--------------------------------------------------------------------------

5. The trace-files. The program reads two columns ascii files.
   The first column is the time. The second column is the seismic trace
   (displacement, velocity, or acceleration).
   The program can read up to 3 registrations per ray.
   The files must be included in the folder /traces.

6. The time-picking file. The program needs the picking of P- and/or
   S-wave direct arrival and the origin time of the earthquake. The order
   of the time-picking file must be the same of the list of the traces.
   The first column is the origin time of the earthquake.
   The second column is the P wave arrival (if not used, can be set to zero).
   The third column is the S-wave arrival (if not used, can be set to zero).
   
%--------------------------------------------------------------------------

      EXAMPLE: 2 rays recorded at 1 three component station.
               list file for the rays: 6 rows
               list file for the traces: 6 rows (1 for each component).
               time picking file: 6 rows (1 for each component even if
               the pickings are identical)

           Rays-list file  Trace-list file           Time file
 Row1      ./rays/ray1     ./traces/ray1-WE    t0ray1  tPray1  tSray1
 Row2      ./rays/ray1     ./traces/ray1-SN    t0ray1  tPray1  tSray1
 Row3      ./rays/ray1     ./traces/ray1-V     t0ray1  tPray1  tSray1
 Row4      ./rays/ray2     ./traces/ray2-WE    t0ray2  tPray2  tSray2
 Row5      ./rays/ray2     ./traces/ray2-SN    t0ray2  tPray2  tSray2
 Row6      ./rays/ray2     ./traces/ray2-V     t0ray2  tPray2  tSray2
    
   See also the example file 'tempi.txt' in the folder /example 

%--------------------------------------------------------------------------

7. Set the variables in the file 'commands.m':

   a. PorS = 2,3;
 
	(2) P-wave tomography

	(3) S-wave tomography 

   b. compon = 1,2,3;
      
	(1)
	
	1 component (usually vertical) recording - fit to work with P-waves

	(2)

	2 component (horizontals) recording- to work with S-waves

	(3)

	3 components recording - the horizontal energy are averaged first

   c. nstep = 1,2;

      (1)

      Allows the operator to obtain the average inverse quality factor
      of the medium, the geometrical spreading, and the constant of the
      inversion problem, as well as the inverse quality factors in a grid
	  of step given by the first velocity model,

      (2)

      the inverse quality factors in a grid of step given by the second
      velocity model, and constrained by the previous inversion, are
      included in the results

   d. modv = load('modv.txt','%f');
      v0= mean(modv(:,4));
	  resol2 = abs(modv(2,3)-modv(1,3))/2;

	  The velocity model with corresponding average velocity and resolution.
   
   e. tempi = load('tempi.txt','%f');

	  The file containing the enucleation time (column 1), the P-wave
	  arrival (second) and S-wave arrival (third). If not used, second or
      third columns can be set to zero.

   f. fin = From 1 to M seconds - Real number

      Length of the window used to measure S- and coda- waves spectra in
      seconds.

   g. cf = From 1 to 18 Hz - Real number

      Central frequency of the frequency band [cf-cf/3 cf+cf/3] where
      spectral amplitude is measured.

   h. tc >> 2*tS - Real number

      The lapse time we use to measure the coda energy. Must be more than
      twice the S-wave travel time. If tC is too large, coda wave energy is
      equal to noise energy and the S-wave spectra are not normalized.
      tC must be set after a signal-to-noise analysis.

   i. tW - Integer number

      Total duration of the coda window used to measure coda energy


   j. originWE, originSN, originz;

	If not included in the velocity coordinates, this is the origin in UTM
	(x = WE and y = SN). The positive topography of the region, if not
	included in the ray. Set them to zero if you are working with  your
	own coordinates.

   k. smoot;

	In order to check the images, the operator can directly choose his
    smoothing parameter. In this case, set it to 1, in the command window
    the program will ask for input, a real number between 0 and Inf. 

   k. modv2 = load('modv2.txt','%f');
	resol22 = abs(modv2(2,3)-modv2(1,3))/2;
    smoot

      The same as before, for the second velocity model.

   The operator can change any of these values and see the differences in
      the attenuation features. The reliability of the images obtained
      must be checked through resolution analysis.

%--------------------------------------------------------------------------
% Output of MuRAT
%--------------------------------------------------------------------------

The outputs of the program depend on the variable nstep:

    nstep	output
      1		Qmean.txt, Q3D.txt, Q_input.txt, Q_check.txt
      2		Qmean.txt, Q3D.txt, Q3D_final.txt, Q_input2.txt, Q_check2.txt

Qmean.txt : The constant, the geometrical spreading, and the inverse average Q
		(first row) and the corresponding errors given by the covariance
		matrix (second).

Q3D.txt : The first Q structure (WE, SN, Z, 1/Q).

Q_input.txt: a simple input structure to check the reliability of the results

Q_check.txt: the result of the inversion on the input structure and he diagonals
	     of the resolution matrix

Q3D_final.txt :  This file is produced if nstep = 2 and shows the combined
		     results of the double inversion.

Q_check2.txt: Output of the tests for the second inversion.


%--------------------------------------------------------------------------
% Plots produced by MuRAT
%--------------------------------------------------------------------------

Figure 1: Up) The logarithm of the direct-to-coda energy ratios on travel time.
	      The fit of the inversion (red continuous line) and the uncertaintiies
	        (from the  covariance matrix of the least square inversion -
		  black dotted lines) are also shown.

	    Down) The log of the coda-to-noise energy ratios.
           The eventual geometrical spreading is also indicated,
           as well as the percentage of measurements with ratio below
           the user-defined soil.

Figure 2: The L-curveof the first inverse problem with corresponding
	    smoothing parameter - see Aster et al. 2005, l_curve function.

Figure 3: The behaviors of the singular elements of the SVD problem
	    (blue line) is compared with the data space (red dots and green
	    crosses in order to evaluate the correct number of blocks
	    efficiently solved by the inversion - see Aster et al. 2005,
	    picard function).

If nstep = 2:

Figure 4 and 5: For the second inversion, corresponding to Figures 2 and 3. 
