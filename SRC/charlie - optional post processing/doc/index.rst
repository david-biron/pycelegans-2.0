.. PyCelegans PostProcess documentation master file, created by
   sphinx-quickstart on Thu Apr 11 17:06:08 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyCelegans PostProcess Documentation
====================================

.. toctree::
   :maxdepth: 2

Overview and Usage
------------------
The purpose of the functions in this module are to check the output of 
PyCelegans and to format it for downstream data analysis. The post-processing 
code is run in two independent modes: **check** and **spline**.

``postprocess.py check ...``

The new integrity checks introduced are:

- Find frames where the midline length changes drastically
- Find frames where the head and tail have been misidentified

In addition, integrity checks from PyCelegans are included:

- Frames of looped worms

To identify frames that do not fit the length requirements, the ``checklength`` module calculates the average length of the worm midline over a sliding window, and locates frames whose length deviates significantly from this value. This comparison can be performed in two modes: *relative* (thresholded by number of standard deviations away from the mean) or *absolute* (thresholded by actual deviation in length away from the mean).

To locate frames where the head/tail orientation is flipped, the ``flipheadtail`` module compares the current frame either (A) to the surrounding frames in a sliding window or (B) to just the previous frame. If (A), it fits each midline to a polynomial, resulting in *m* x 2 (where *m* is the degree of the polynomial) coefficients for each frame. These coefficients are averaged over a sliding window, and the values for each individual frame compared to this mean to find frames where the posture deviates significantly from the local value. This process is performed iteratively to account for outliers, but still relies upon the assumption that the majority of the frames within the window are correct. If the score calculate by comparing the coefficients for one frame to the local average is greater than a maximum threshold, the midline is flipped. If it is less than a minimum threshold, nothing happens. If the score falls between the minimum and maximum thresholds, it is (optionally) marked as "indeterminate" and rejected from further analysis. If (B), it compares the total distance between each point along the current midline and the midline of the previous frame, oriented either head-to-tail or tail-to-head. If the latter number is smaller, a flip is identified.

After running the post-process module in ``check`` mode, the number of frames that passed each check will be printed. The outputs of each integrity check are saved as a .txt file in the original ``properties`` directory.

``postprocess.py spline ...``

The midline is fit to a smoothing B-spline, with the option to include the results of the integrity checks. The integrity checks are included in two steps: bad frames may be rejected, and flipped worms may optionally be reoriented. If any integrity checks are missing, they are simply ignored; i.e., if ``postprocess.py spline`` if run before ``postprocess.py check``, only the ``is_loop`` criterion will be used. The length of the midline is also calculated by evaluating the integral of (*dx*:sup:`2` + *dy*:sup:`2`) :sup:`1/2` along the midline. The sidepaths may be optionally splined, resulting in a closed B-spline that defines the perimeter of the worm. In this case the worm area is also calculated, by evaluating the integral of (*xdy* -- *ydx*) / 2 along the worm perimeter. These values, along with their corresponding frame number, are finally exported to a .mat file for analysis in Matlab. *Note:* the format of the Matlab ouput will depend upon how the splines are generated; if a constant number of points are used then the midline and sidepaths will be matrices of size *m* x *n* x 2 where *m* is the number of frames (that passed all integrity checks, if used) and *n* is the number of points in the spline. If  constant spacing is used, *n* may vary so the results will be cells.

After running the post-process module in ``spline`` mode, the number of frames that threw no error will be printed. If analyzing data from the ``./properties`` directory, the splined .txt files will be saved to a separate directory as well as to a .mat file (by default, ``./processed`` and ``./processed.mat``, respectively).

Run integrity checks
--------------------

To view the help section for the integrity checks, type

``postprocess.py check --help``

which will display::
    
    usage: postprocess.py check [-h] [-p PROPERTIES] [-f FRAME_RANGE FRAME_RANGE]
                                [-v {0,1}] [-ls LEN_SPAN]
                                [-lv LEN_MAX_VAL | -ld LEN_MAX_DEV]
                                [-fd FLIP_DEGREE] [-fi FLIP_ITERS] [-fs FLIP_SPAN]
                                [-fl FLIP_MIN_DEV] [-fu FLIP_MAX_DEV]

    optional arguments:
      -h, --help            show this help message and exit
      -p PROPERTIES, --properties PROPERTIES
                            properties directory; default = ./properties
      -f FRAME_RANGE FRAME_RANGE, --frame_range FRAME_RANGE FRAME_RANGE
                            frame range, given as [first frame, last frame);
                            default = '0 inf'
      -v {0,1}, --version {0,1}
                            version of process code used to find midline (0 =
                            Nick's, 1 = Marc's); default = 0
      -ls LEN_SPAN, --len_span LEN_SPAN
                            size of window used to check lengths, over which the
                            average worm length should be constant; default = 6000
      -lv LEN_MAX_VAL, --len_max_val LEN_MAX_VAL
                            maximum allowed value of length outside of average
                            value in the moving window set by the span, set as a
                            percentage; set this value but not len_max_dev to use
                            an absolute value of the threshold instead of
                            comparing to a deviation about the average; default =
                            0.05
      -ld LEN_MAX_DEV, --len_max_dev LEN_MAX_DEV
                            maximum allowed deviation of length from average value
                            of the standard score in the moving window set by the
                            span; default = 0.00
      -fd FLIP_DEGREE, --flip_degree FLIP_DEGREE
                            degree used in polynomial fit to determine correct
                            orientation; default = 10
      -fi FLIP_ITERS, --flip_iters FLIP_ITERS
                            number of times to run flip-check procedure; default =
                            3
      -fs FLIP_SPAN, --flip_span FLIP_SPAN
                            size of window used to check flips, over which the
                            posture of the worm should be approximately constant;
                            default = 20
      -fl FLIP_MIN_DEV, --flip_min_dev FLIP_MIN_DEV
                            minimum deviation of coefficients used to fit worm
                            midline in moving window, below which the orientation
                            is considered correct; default = 2.00
      -fu FLIP_MAX_DEV, --flip_max_dev FLIP_MAX_DEV
                            maximum deviation of coefficients used to fit worm
                            midline in moving window, above which the orientation
                            is considered flipped; default = 3.00

To check the first 100 frames of an experiment in ``./properties`` 
using Nick's midline and default parameters:

``postprocess.py check -p ./properties -f 0 100``

To check the frames 3000 through 4000 (inclusive), set the ``frame_range``:

``postprocess.py check -f 3000 4001``

To check every frame of an experiment using Marc's midline, set the ``version`` flag = 1:

``postprocess.py check -v 1``

To use a more restrictive threshold for the length-check, decrease the maximum allowed deviations by setting ``len_max_val``:

``postprocess.py check -lv 0.02``

To use a less restrictive threshold for the length-check, increase the maximum allowed deviation by setting ``len_max_val``:

``postprocess.py check -lv 0.1``

To use the standard score instead of the fractional change in length, set ``len_max_dev``:

``postprocess.py check -ld 3.0``

To change the window size used to check for bad midline lengths (e.g., to 10 min at a frame rate of 10 fps), set ``length_span``:

``postprocess.py check -ls 6000``

To correct for head/tail flips without rejecting "indeterminate" frames, set ``flip_min_dev`` and ``flip_max_dev`` (the lower and upper thresholds, respectively) equal:

``postprocess.py check -fl 3.0 -fu 3.0``

To change the window size used to check for flips (e.g., to 3 sec at a frame rate of 10 fps), set ``flip_span`` to any integer > 1:

``postprocess.py check -fs 30``

To check for flips by comparing just to the previous frame (calculating distances between points along a splined midline, instead of comparing the coefficients of the polynomial fit), set ``flip_span`` = 1:

``postprocess.py check -fs 1``

Fit to smoothing spline
-----------------------

To view the help section for splining the data, type

``postprocess.py spline --help``

which will display::
    
    usage: postprocess.py spline [-h] [-p PROPERTIES] [-f FRAME_RANGE FRAME_RANGE]
                                 [-v {0,1}] [-c] [-r] [-k] [-o OUTPUT]
                                 [-s SMOOTHING] [-n NUM_POINTS | -sp SPACING]
                                 [-l LENGTH]

    optional arguments:
      -h, --help            show this help message and exit
      -p PROPERTIES, --properties PROPERTIES
                            properties directory; default = ./properties
      -f FRAME_RANGE FRAME_RANGE, --frame_range FRAME_RANGE FRAME_RANGE
                            frame range, given as [first frame, last frame);
                            default = '0 inf'
      -v {0,1}, --version {0,1}
                            version of process code used to find midline (0 =
                            Nick's, 1 = Marc's); default = 0
      -c, --use_checks      include this flag to use info from various post-
                            processing checks when calculating splines
      -r, --reorient        include this flag to reorient tail-to-head any worms
                            marked as incorrectly flipped; this flag will be
                            ignored if the '--use_checks' flag is not used
      -k, --keep_sides      include this flag to read the sidepaths and save a
                            single splined worm outline to file
      -o OUTPUT, --output OUTPUT
                            name of directory and .mat file to save splined data;
                            default = ./processed
      -s SMOOTHING, --smoothing SMOOTHING
                            level of smoothing in creating splined version of worm
                            midline; should be of the order of the number of
                            points in the midline; default = 25
      -n NUM_POINTS, --num_points NUM_POINTS
                            number of points in the each splined worm; default =
                            100
      -sp SPACING, --spacing SPACING
                            distance between consecutive (x, y) points in each
                            splined worm; this value is exclusive of the length
                            and specifying it may result in splines with unequal
                            numbers of points; default = 0.00
      -l LENGTH, --length LENGTH
                            desired length of each splined midline; if specified
                            all splined worms will be approximately the same
                            length (depending on whether a uniform number of
                            points or spacing is given); default = 0.00

To spline the first 100 frames of an experiment in ``./properties`` 
using Nick's midline and default parameters, and without using any automatic
corrections from the integrity checks:

``postprocess.py spline -f 0 100``

To spline every frame using Marc's midline without input from integrity checks:

``postprocess.py spline -v 1``

To include integrity checks (but ignore flips), just set the ``use_checks`` flag:

``postprocess.py spline -c``

To include integrity checks (and automatically reorient tail-to-head flips), set both the ``use_checks`` and ``reorient`` flags:

``postprocess.py spline -cr``

To also include sidepaths and worm areas, set the ``keep_sides`` flag:

``postprocess.py spline -crk``

To make the spline smoother, increase the value of ``smoothing``:

``postprocess.py spline -s 50``

To use half as many points in the midline, set the value of ``num_points``:

``postprocess.py spline -n 50``

To use a spline in which the points are spaced an equal number of pixels apart (this will produce results with a non-uniform number of points across splines), set the ``spacing``:

``postprocess.py spline -sp 5.0``

*Note:* if the sidepaths are also splined, they are automatically generated with twice the number of points in the midline, but with twice the smoothing factor (which should result in approximately the same amount of local smoothing).

:mod:`postprocess` Module
-------------------------

.. automodule:: postprocess
    :members:
    :undoc-members:
    :show-inheritance:

:mod:`importdata` Module
------------------------

.. automodule:: importdata
    :members:
    :undoc-members:
    :show-inheritance:

:mod:`checklength` Module
-------------------------

.. automodule:: checklength
    :members:
    :undoc-members:
    :show-inheritance:

:mod:`flipheadtail` Module
--------------------------

.. automodule:: flipheadtail
    :members:
    :undoc-members:
    :show-inheritance:

:mod:`splinefit` Module
-----------------------

.. automodule:: splinefit
    :members:
    :undoc-members:
    :show-inheritance:

:mod:`calc` Module
------------------

.. automodule:: calc
    :members:
    :undoc-members:
    :show-inheritance:

:mod:`parseargs` Module
-----------------------

.. automodule:: parseargs
    :members:
    :undoc-members:
    :show-inheritance:

:mod:`readwrite` Module
-----------------------

.. automodule:: readwrite
    :members:
    :undoc-members:
    :show-inheritance:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`