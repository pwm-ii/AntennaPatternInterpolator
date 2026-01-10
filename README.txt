
                     ANTENNA PATTERN INTERPOLATOR README
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


1.) PURPOSE
---------------------------------------------------------------------------
This tool is designed to reconstruct a full 3D antenna radiation pattern 
using only two 2D principal cuts:
  1. The Azimuth Cut
  2. The Elevation Cut


2.) REQUIREMENTS
---------------------------------------------------------------------------
Language: Python 3.x
Required Libraries: numpy, scipy, matplotlib, tkinter


3.) INPUT FILE FORMAT
---------------------------------------------------------------------------
The tool accepts text files (.ant, .txt, .csv) containing a SINGLE column 
of gain values (in dB).

Structure Requirements:
  - The total number of lines MUST be an EVEN number.
  - The lines are split in half between Azimuth and Elevation Data.
        - Top Half of file = Azimuth Data.
        - Bottom Half of file = Elevation Data.
  - Minimum Length: The file must have at least 10 lines total (5 for Azimuth, 5 for Elevation).
  - Data should ideally be NORMALIZED (Peak Gain = 0 dB) before loading.

Example Layout (file.ant):
  -0.5    <-- Azimuth 0 deg
  -0.2    <-- Azimuth 1 deg

          [...]

  -0.5    <-- Azimuth 359 deg
  -10.0   <-- Elevation 0 deg
  -10.1   <-- Elevation 1 deg
   
          [...]

  -10.0   <-- Elevation 359 deg


4.) INTERPOLATION METHODS
---------------------------------------------------------------------------
A. Summing
   - Methodology: Mathematically adds the logarithmic/decibel patterns 
     (equivalent to multiplying the linear patterns).
   - Requirement: Inputs MUST be normalized to 0 dB, or result will be 
     double the actual gain.
   - Perfectly reconstructs omni-directional patterns with no error due to 
     axial symmetry.
   - Relatively good at approximating the main lobe of directive antennas 
     but can fail with side lobes.
   - Systematically underestimates gain in complex directional antennas.
   - Tends to create "creases" where the two cuts intersect.

B. Approximation (Cross-Weighted)
   - Methodology: A type of summing algorithm with geometric cross-weighting 
     (Azimuth pattern weights the Elevation samples and Elevation pattern 
     weights the Azimuth samples).
   - Best for generally directional antennas where the Summing algorithm 
     fails to capture side-lobe details.
   - This method generally produces more optimistic (higher gain) estimates 
     than summing.
   - Normalization Parameter (k): Controls the mathematical "locus" 
     of the weights for the Approximation algorithm.
       - k=1: Linear weighting.
       - k=2 (DEFAULT): Minimizes approximation error for most antennas.
       - Increasing k produces a more conservative (lower gain) estimate.
       - k=infinty: Approximation algorithm = Summing algorithm

C. Hybrid (Optimized)
   - Methodology: A weighted mix that uses the Summing algorithm for the 
     Main Lobe (where it is most accurate) and the Approximation algorithm 
     for the side lobes.
   - Best for high-gain directional arrays.
   - Transition Slope Parameter (n): Controls the "slope" of the 
     bridging function between the two algorithms.
       - Higher n: Summing Algorithm (which is less accurate in side lobes) 
         dominates a larger portion of the pattern. In effect, this improves 
         the Main Lobe accuracy but can degrade overall mean error.
       - Recommended Range: 3.5 to 6.0.
           - Default is set to 6.0
           - Yagi-Uda Antenna: Optimal n = 3.5
           - Panel Antenna: Optimal n = 5.0


5.) SETTINGS & OPTIONS
---------------------------------------------------------------------------
[ ] Auto-Center Peaks
    - If checked, the tool shifts the data so the main lobe (peak gain) is 
      aligned to 0 degrees.
    - Keep on unless your measurements are already perfectly aligned 
      to the coordinate system.

[ ] Enforce Loop Closure
    - Ensures the value at 0 degrees matches the value at 360 degrees.
    - Keep on to prevent gaps in the 3D mesh.

[ ] Enable 3D Surface Smoothing (Gaussian Filter)
    - Applies a low-pass filter to the generated 3D data.

6.) EXPORT FORMAT
---------------------------------------------------------------------------
The tool exports a CSV file with three columns: Phi, Theta, and Gain.
  - Phi range:   0 to 359
  - Theta range: 0 to 179
