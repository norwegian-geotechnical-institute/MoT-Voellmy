com9MoTVoellmy: Voellmy-Type Model for Dense Flow Avalanches
============================================================
<br/>


General remarks
---------------

`com9MoT-Voellmy` is a module for dense flow (snow) avalanche computations (DFA) using the model MoT-Voellmy developed at the Norwegian Geotechnical Institute (NGI). The computational kernel is an ISO-C code, with compiled executables available for Linux, MS Windows and macOS. These executables can also be run from the command line, with a run configuration file (RCF) as input. The `com9MoT-Voellmy` module in `AvaFrame` provides a graphical user interface that helps the user with pre-preprocessing the input data, sets up the RCF, runs the executable once or multiple times, and visualizes the simulation results on the screen and in a comprehensive report.

The configuration can be modified in order to change any of the default settings and also allows to perform simulations for varying parameters all at once.

The present document uses the terms *thickness averaged/integrated* for consistency with the `AvaFrame` documentation, whereas in the literature on gravity mass flow models and in the [NGI Technical Note on MoT-Voellmy](x-pdf:./20230100-06-TN.MoT-Voellmy_eqs_numerics.pdf) the terms *depth averaged/integrated* are usually used.

Dense flow avalanche simulations with MoT-Voellmy can be performed for different release area scenarios, with or without
entrainment. The friction parameters of the Voellmy model, *μ* and *k* (equivalent to *g*/*ξ* as used in many models and publications), can be specified as fixed values in the entire avalanche path or read from a raster file. It is possible to account for the braking effect of forest in parts of the avalanche path due to the extra drag proportional to the stand density and local average trunk diameter. Furthermore, where and when the bending moment exerted by the avalanche exceeds the bending strength of the trees, tree breaking and the concomitant reduction of the forest braking effect can be modeled.

**Notes**:

- As with all numerical avalanche models, MoT-Voellmy must be used with caution and the results must be carefully and critically analyzed. If the initial conditions (release area, fracture depth, snow conditions if entrainment is included) or friction parameters are not chosen according to the specific situation, the simulated run-out, velocity and pressure can be misleading, with potentially grave consequences.
- The simulation procedures and recommended parameter values elaborated for Austria in the `com1DFA` module need not be directly applicable in countries with different climate, such as Norway.
- An extensive model comparison on a simplified "hockeystick" slope and a real-world avalanche path with identical initial conditions and friction parameter values has shown tha `com1DFA` and `MoT-Voellmy` give broadly similar results. However, there can be non-negligible discrepancies due to differences in the numerical schemes and, in particular, the stopping criteria, that are applied.
- The configuration provided with com1DFA is well-tested and applied for hazard mapping in Austria. It is therefore expected to give a good starting point for simulations with MoT-Voellmy under similar conditions. As for `com1DFA`, unwanted/unexpected/spurious side-effects may appear if the model is used in a context it is not developed and tested for, like slushflows, debris flows or rock avalanches.
<br/><br/>
   

Theory and modeling options
---------------------------

For a more comprehensive description of the geometrical setting, the fundamental balance equations and modeling options of `MoT-Voellmy`, see the [Technical Note on MoT-Voellmy](x-pdf:./20230100-06-TN.MoT-Voellmy_eqs_numerics.pdf). The geometrical formulation and the balance equations (but not the numerical scheme) are almost the same as for RAMMS::AVALANCHE (not the *Extended* version) and are well described by [Christen et al. (2010)](https://doi.org/10.1016/j.coldregions.2010.04.005). While `com1DFA` and MoT-Voellmy share many features, there are some differences that users should be aware of. The largest differences concern the numerics, but these are less visible to the users unless numerical instabilities arise. Both models solve similar equations, but `com1DFA` offers some extensions to the Voellmy expression for the shear stress, *τ<sub>b</sub>*, at the interface between the flowing avalanche and the snow cover,

.. math::
     \mathbf{\tau}_b = -\frac{\|}{\|\mathbf{u}}
                        \left[ \mu \rho_f h g_\perp + k \rho_f \mathbf{u}^2 \right],

where *ρ<sub>f</sub>* is the flow density (assumed constant in both models), *h* the flow thickness, *g*<sub>⟂</sub> the component of the gravitational acceleration perpendicular to the slope, **u** the (thickness-averaged) flow velocity parallel to the slope, and *μ* and *k* are the dimensionless dry-friction and drag coefficients, respectively. For details on the extensions offered by `com1DFA`, see its documentation[documentation](https://docs.avaframe.org/en/latest/moduleCom1DFA.html) .

In simulations with **entrainment**, MoT-Voellmy ensures that in each cell at most the snow quantity declared as erodible is entrained. The model offers five different choices for the formula determining the entrainment rate:

* `none` sets the entrainment rate to 0, no erodible snow-cover thickness or other snow parameters need to be specified.
* `RAMMS` implements a formula available in some early releases of RAMMS::AVALANCHE, with the entrainment rate (= entrained mass per unit time and unit area measured parallel to the slope, in kg/(m² s)) given by .. math:: q_e = c \rho_b \|\mathbf{u}\|. *ρ<sub>b</sub>* is the density of the snow cover, *c* a dimensionless entrainment coefficient, to be specified in the run configutration file. In RAMMS::AVALANCHE, suggested values were 0.1 < *c* < 1, which essentially lead to frontal entrainment. With values *c* < 0.01, scouring entrainment along the avalanche body is obtained. A raster file containing the erodible snow thicknesses in each cell (in m) must be supplied.<br/>
The shear strength of the snow cover is not considered in this model. A very thin but fast flow layer may therefore entrain unrealistically large quantities of snow.
* `AvaFrame` implements the basal entrainment formula of `com1DFA` (the plowing formula would require front tracking, which is not implemented in MoT-Voellmy):<br/>
.. math::
    q_e = \frac{\|\mathbf{\tau}(\mathbf{u}, h)\|}{e_b} \|\mathbf{u}\| .
<br/>
In addition to the erodible snow thickness, the specific erosion energy (units m²/s²) must be specified either as a constant or in an AAIGrid raster file. This option has not been tested extensively, and users are directed to the documentation of `com1DFA` and related publications for guidance about choosing *e<sub>b</sub>*.
* `TJEM` (Tangential Jump Entrainment Model) determines the basal entrainment rate from a balance between the shear stress that the flow would exert on a non-erodible bed and the shear strength *τ<sub>c</sub>* of the snow cover. The difference max(*τ<sub>b</sub>*'−*τ<sub>c</sub>*, 0) how much eroded snow (per unit area) can be accelerated to the thickness-averaged avalanche velocity, resulting in<br/>
.. math::
    q_e = \frac{max(*τ<sub>b</sub>*'−*τ<sub>c</sub>*, 0)}{\rho_f \|\mathbf{u}\|} .
<br/>
The true shear stress at the interface between snow cover and avalanche is limited to *τ<sbu>c</sub>* if there is entrainment.<br/>
As in the other models, the erodible snow depth must be given in an AAIGrid raster file. In addition, *τ<sbu>c</sub> (in Pa) is to be specified, either as a constant value or through an AAIGrid raster file.<br/>
Optionally, the snow shear strength may be set to increase linearly with depth into the snow cover. The user can either set a constant value (in Pa/m) or provide an AAIGrid file with spatially variable values.
* `GOEM` (Grigorian–Ostroumov erosion model) provides the basal scour model proposed by Grigorian and Ostroumov in the mid-1970s, here extended to a 2D setting. See the [Technical Note on MoT-Voellmy](x-pdf:./20230100-06-TN.MoT-Voellmy_eqs_numerics.pdf) and the references therein for an explanation of the theory behind the model.<br/>
Similar to the `TJEM`, the `GOEM` requires a constant value or an AAIGrid file for the snow-cover strength (more precisely, a combination of compressive and shear strength). Additionally, a dimensionless and constant empirical coefficient *c* must be specified; it is used to mimick the shear stresses at the erosion front, which are not explicitly included in this model. The entrainment rate is very sensitive to *c* if *c* is not much smaller than 1. As in `RAMMS`, values above about 0.1 lead to quasi-frontal entrainment; for basal scour, *c* ~ *k* is recommended.

**Forest effects** are modeled quite differently from `com1DFA`. Intact trees are treated as a source of extra friction. The cell-averaged extra shear stress is equal to the drag force on a single tree of (local) average diameter times the number of trees per unit slope area. This implies that interactions between the flow disturbances created by neighboring trees are neglected, which is expected to lead to an under-estimation of the effect. Detrainment, which is a central mechaniscm in `com1DFA`, is disregarded because the model by [Feistl et al. (2014)](https://doi.org/10.3189/2014JoG13J055) is only applicable to small avalanches. The drag coefficient of cylindrical obstacles like trees depends strongly on the Reynolds number of the flow, which closely correlates with the Froude number for friction-dominated granular flows. Using results from laboratory experiments and numerical simulations, the braking effect of a forest is modeled by modifying *μ* and *k* at each time step by Δ*μ* = 1.25 *n D h* and Δ*k* = 0.5 *n D h*, where *h* is the instantaneous local flow thickness, *D* the average local trunk diameter at breast height, and *n* is the number the number of trees per unit slope area (i.e., measured obliquely).

The bending moment exerted on a tree by the avalanche is approximated as *M<sub>a</sub>* ~ *ρ<sub>f</sub>* **u**<sup>2</sup> *D h*<sup>2</sup>. The resistance of the trees is determined by the tree diameter and the modulus of rupture (MoR) as *M<sub>c</sub>* = (π/32) MoR *D*<sup>3</sup>. Once *M<sub>a</sub>* > *M<sub>c</sub>* in a cell, all trees in the cell are assumed to break (or be uprooted). This is captured by reducing *n* in the cell exponentially with time, *n*(*t*−*t*<sub>1</sub>) = *n*(*t*<sub>1</sub>)·e<sup>−*λ*·(*t*−*t*<sub>1</sub>)</sup>. The decay constant *λ* is assumed to be inversely proportional to the average trunk diameter, *λ* = 0.15 s/m / *D* (users can change 0.15 s/m to another value). The modulus of rupture is of the order of magnitude of 50 MPa, but can vary by at least a factor 2 in either direction depending on the tree species and the water content of the trees. This value can also be modified by the user.
<br/><br/>


Numerics
--------

Those equations are solved numerically using a **SPH** method (:cite:`LiLi2010,Sa2007`).
**SPH**  is a mesh free method where the basic idea is to divide the avalanche into
small mass particles. The particles interact with each other according to the
equation of motion described in :ref:`moduleCom1DFA:Theory` and the chosen kernel function.
This kernel function describes the domain of influence of a particle (through the smoothing length parameter).
See theory :ref:`theoryCom1DFA:com1DFA DFA-Kernel theory` for further details.
<br/><br/>


Input
-----

Simulations are performed within an avalanche directory, organized with the folder structure described below.

.. Note::  An avalanche directory can be created by running: :py:mod:`runInitializeProject.py`, which creates the required folder structure:

  ::

    NameOfAvalanche/
      Inputs/
        REL/      - release area scenario
        RES/      - resistance areas
        ENT/      - entrainment areas
        POINTS/   - split points
        LINES/    - avalanche paths
        POLYGONS/ - crop shapes
        SECREL/   - secondary release areas (not used in com9MoTVoellmy)
        RASTERS/  – friction parameter fields
      Outputs/
      Work/


In the directory ``Inputs``, the files listed below are required. Be aware that ALL inputs have to be provided in the same
projection.

* **digital elevation model** as raster file in `ESRI ASCII grid format` (AAIGrid) <https://desktop.arcgis.com/en/arcmap/10.3/manage-data/raster-and-images/esri-ascii-raster-format.htm>,
* **release area scenario** in Inputs/REL, either as (multi-)polygon shapefile (multiple features are possible) or as AAIGrid raster file with the same extent and resolution as the DEM

If a shapefile is used for the release area scenario, the following restrictions apply:

* The release area polygon must not contain any "holes" or inner rings.
* The release area name should not contain an underscore, if so '_AF' is added.
* The recommended attributes of the features are *name*, *thickness* (see :ref:`moduleCom1DFA:Release-, entrainment thickness settings`) and *ci95* (see :ref:`moduleAna4Stats:probAna - Probability maps`).
* ALL features within one shapefile are released at the same time (and interact), this is what we refer to as *scenario*.
* If you want to simulate different scenarios with the same features, you have to copy them to separate shapefiles.

With AAIGrid files, there are fewer restrictions:

- The release area name should not contain an underscore, if so '_AF' is added.
- ALL cells with non-zero release thickness are released at the same time (and interact), this is what we refer to as *scenario*.
- It is possible to have multiple, disconnected release areas, and they may contain holes.
- The release thickness may be spatially varying.


The following files are optional. Please note: in the standard configuration (i.e., ``simTypeList = available``),
the *null* variant is always run! This means that, if a resistance and/or an entrainment file is given (as described below),
at least two results are generated: the *null* variant and the variant with entrainment and/or resistance.

Polygons (shapefiles):

* one entrainment area (multi-) polygon shapefile (in Inputs/ENT)
       - marks the (multiple) areas where entrainment can occur.
       - attribute *thickness* (see :ref:`moduleCom1DFA:Release-, entrainment thickness settings`)
       - must not contain any "holes" or inner rings
       - Alternatively, the user can provide a AAIGrid raster file of erodible snow depth and— depending on the specified entrainment model—additional raster files for snow shear strength or erosion energy and for the bed-normal gradient of snow shear strength.
* one resistance area (multi-)polygon shapefile (in Inputs/RES)
      - The features mark the (multiple) areas where resistance due to forest is considered.
      - The resistance areas must not contain any "holes" or inner rings.
      - Please consider the information about resistance below :ref:`moduleCom1DFA:Resistance setup`
* AAIGrid raster files for the Voellmy friction parameters *μ* and *k* (in Inputs/RASTERS)
    - The rasters must have the same extent as the DEM.
    - The file names must end with ``_mu.asc`` and ``_k.asc``.
    - If ``meshCellSize`` is different from simulation ``meshCellSize``, the fields will be remeshed.
    - Only used if ``Parameters`` is set to ``variable``.
* one ``_cropshape.shp`` shape file (in Inputs/POLYGONS)
    - Provides a polygon located inside the DEM that defines the area  (bounds of polygon) for which the peak fields are plotted in the simulation report.
    - If no polygon is provided, the peak fields are shown for the extent where peak field values are nonzero.



Release-, entrainment thickness settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note::
    Thickness is unambiguous: it is measured normal to the slope.

Release, entrainment and secondary release thickness can be specified in two different ways:

1. Via **shape file**:

    - Add an attribute called `thickness` for each feature.
    - Important: ALL features have to have a single thickness value, which can differ between features.
    - For entrainment area only: If the thickness value is missing, the thickness value is taken from `entThIfMissingInShp` (default 0.3 m) in the configuration file. If multiple features are in the entrainment file, the thickness attribute has to be set either for ALL or NONE of the features.
    - For backwards compatibility, the attribute 'd0' also works, but we suggest to use `thickness` in new projects.
    - Set the flag `THICKNESSFromShp` (i.e. relThFromShp, entThFromShp, secondaryRelthFromShp) to True in the configuration file (default is True).
    - A parameter variation can be added with the `THICKNESSPercentVariation` parameter in the configuration file in the form of ``+-percentage$numberOfSteps``. Provided a `+`, a positive variation will be performed, if `-` is given, only a negative variation is performed. If no sign is given, both directions will be used. Additionally, a variation can be added with the `THICKNESSRangeVariation` parameter in the configuration file in the form of ``+-range$numberOfSteps``. Provided a `+`, a positive variation will be performed, if `-` is given, only a negative variation is performed. If no sign is given, both directions will be used. Furthermore, there is the option to vary the thickness in a range of ± the 95% confidence interval value, which is also read from the shape file (requires an attribute called ci95). In order to use this variation, set the 'THICKESSRangeFromCiVariation' to ``ci95$numberOfSteps``.

2. Via **configuration file (ini)**:

    - set the flag 'THICKNESSFromShp' to False
    - provide your desired thickness value in the respective THICKNESS parameter (i.e. relTh, entTh or secondaryRelth)
    - in addition to the `THICKNESSPercentVariation` and `THICKNESSRangeVariation`
    options (see option 1) and the standard variation options in
    :ref:`configuration:Configuration`, you can also directly set e.g. `relTh =
    1.$50$2`, ``referenceValue$+-percentage$numberOfSteps``, resulting in a
    variation of relTh from 0.5 to 1.5m in two steps.

Only available for release thickness:

3. Via **release thickness file**:

    - set the flag 'relThFromShp' to False
    - set the flag 'relThFromFile' to True
    - save a raster file with info on release thickness as raster file in
    ``Inputs/RELTH`` the number of rows and columns must match the DEM raster
    with desired meshCellSize (recommended)
  - if the cellsize does not match the requested meshCellSize, the file is
    remeshed.


Friction parameters
^^^^^^^^^^^^^^^^^^^

By default the friction parameter set *samosATAuto* is active. This uses the calculated release volume (including
secondary release areas) to determine the parameters used for the samosAT friction model.
See :ref:`samosatfrict` for the limits regarding release volumes.


DEM input data
^^^^^^^^^^^^^^^^
Regarding the DEM data: if the DEM in ``Inputs`` is not of cell size 5 meters, it is remeshed to a
cell size of 5 meters. However, it is also possible to specify a desired cell size in the
configuration file (parameter `meshCellSize`). In this case, also consider reading :ref:`FAQ:Can the spatial resolution of simulations performed with com1DFA (dense flow) be changed?`.
If the cell size of the DEM in ``Inputs`` is equal to the desired mesh cell size, the DEM is used without modification. If the cell sizes do not match, several options are available:

    - cleanremeshedRasters = True, directory ``Inputs/remeshedRasters`` is cleaned, and the DEM in Inputs/
      is remeshed to the desired cell size - this is the default setting

    - cleanremeshedRasters = False and a DEM including the name of the DEM in Inputs/ and the desired cell size is found
      in Inputs/remeshedRasters - this DEM is used without modification

    - cleanremeshedRasters = False and no matching DEM is found in Inputs/remeshedRasters - the DEM in Inputs/ is remeshed
      to the desired cell size

If the DEM in Inputs/ is remeshed, it is then saved to ``Inputs/remeshedRasters`` and available for subsequent
simulations.


Dam input
^^^^^^^^^

The com1DFA module provides the option to take the effect of dams into account.
This is done using a ad-hoc method based on particles being reflected/deflected by a dam wall.

The dam is described by the crown line, the slope and the restitution coefficient:

  - crown line as shape file (use the line type and enable the "additional dimensions" option in order
    to specify the z coordinate).
    The z coordinate corresponds to the absolute height (terrain elevation plus dam height).
    The dam is then located on the left side of the dam (when one travels from the first point to the last
    point of the shapefile line).
    The dam shape files live in the ``avaDir/Inputs/DAM/`` directory (only one file is allowed).

  - the ``slope`` of the dam (in degrees °) between the horizontal plane and the wall to be provided in the shape file
    as an attribute (default value is 60° in the provided examples: avaSlide, avaKot and avaBowl)

  - the restitution coefficient (:math:`\alpha_\text{rest}`), a float between 0 (no reflection
    in the normal direction) and 1 (full reflection) to be specified in the ini file (default value is 0)
<br/></br/>



Model configuration
-------------------

The model configuration is read from a configuration file: ``com1DFA/com1DFACfg.ini``. In this file,
all model parameters are listed and can be modified. We recommend to create a local copy
and keep the default configuration in ``com1DFA/com1DFACfg.ini`` untouched.
For this purpose, in ``AvaFrame/avaframe/`` run:

  ::

    cp com1DFA/com1DFACfg.ini com1DFA/local_com1DFACfg.ini

and modify the parameter values in there. For more information see :ref:`configuration:Configuration`.

It is also possible to perform multiple simulations at once, with varying input parameters.
<br/></br/>


Output
------

Using the default configuration, the simulation results are saved to: *Outputs/com9MoTVoellmy* and include:

* raster files of the peak values for pressure, flow thickness and flow velocity (*Outputs/com9MoTVoellmy/peakFiles*)
* raster files of the peak values for pressure, flow thickness and flow velocity for the initial time step (*Outputs/com9MoTVoellmy/peakFiles/timeSteps*)
* markdown report including figures for all simulations (*Outputs/com9MoTVoellmy/reports*)
    - If a ``_cropshape.shp`` file is provided in Inputs/POLYGONS, plots are cropped to the rectangular bounds of the polygon.
    - If ``showOnlineBackground = True`` in avaFrameCfg.ini and a suitable ``mapProvider`` is set, peak fields are plotted onto the corresponding map.
* log files of all simulations (*Outputs/com* configuration files for all simulations (*Outputs/com9MoTVoellmy1DFA/configurationFiles*)
    - All configuration files that were created for a simulation to be run are stored in (*Outputs/com9MoTVoellmy/configurationFiles*)
    - One file for each simulation that has actually been performed is saved in (*Outputs/com9MoTVoellmy/configurationFiles/configurationFilesDone*)
    - One file for each simulation that has actually been performed by the latest call of ``runCom9MoTVoellmy.py`` is saved in (*Outputs/com9MoTVoellmy/configurationFiles/configurationFilesLatest*)

    .. Note::
        This kind of storage of configurations from actually performed simulations allows a run that has been terminated
        to be resumed without re-running simulations that have already been performed. For this, just restart the run.

The naming of the output files has the following structure, shown with the example of
*relAlr_ff5f9b78c6_C_L_null_dfa_ppr*:

* *relAlr* – release area name, usually the name of the shapefile
* *ff5f9b78c6* – individual hash of the configuration file used for the simulation. All files related to this simulation have the same hash in their name. This allows to identify which files belong to which simulation.
* *C* – indicator of the setup used: D for default setup, C for custom setup, i.e., something was changed in the configuration file.
* *L* – indicator of the size category used for the friction model: L for large, M for medium, S for small
* *null* – indicator of the run type: null for null variant, ent for entrainment variant, res for resistance variant, etc
* *dfa* – indicator of the simulation type: dfa for dense flow avalanche
* *ppr* – indicator of the result type: ppr for peak pressure, pfv for peak flow velocity, pft for peak flow thickness, etc


Optional outputs

In the configuration file, it is possible to select the frequency of storing flow configurations and the flow variables that shall be exported.
The result types that can be chosen to be exported are (all correspond to fields):

* ppr – peak pressure
  (:math:`pressure = \mathbf{\rho}  \mathbf{u}²` with :math:`\rho` snow density and :math:`\mathbf{u}` flow velocity)
* pfv – peak flow velocity
* pft – peak flow thickness
* pta – peak travel angle
* FV – flow velocity
* FT – flow thickness
* P – pressure
* FM - flow mass
* Vx, Vy – *x* and *y* components of velocity (parallel to the slope)
* TA – travel angle
* demAdapted – modified DEM if the option `Evolving geometry` is set to `yes`

Have a look at the designated subsection Output in ``com9MoTVoellmy/com9MoTVoellmyCfg.ini``.


Parallel computation
--------------------

If multiple runs of com1DFA are to be executed, these will be calulated in parallel via
multiprocessing. So each task itself is calculated on only one core, but different tasks
are run at the same time.

This happens if you have one of the following (or a combination of them):

* multiple scenarios (multiple input release shapefiles)
* multiple runtypes, i.e null variant and entrainment/resistance variant (e.g.: simTypeList = null|ent)
* some kind of parameter variation (e.g.: relTh = 1.0|1.5|1.7)

The number of CPU cores is controlled in the main ``avaframeCfg.ini`` file. By default a
maximimum of 50 percent of your available cores is being utilized. However you can set
a different number if needed. For sequential execution set nCPU to 1.
<br/><br/>


To run
--------

* First go to ``AvaFrame/avaframe``
* Copy ``avaframeCfg.ini`` to ``local_avaframeCfg.ini`` and set your desired avalanche directory name
* Create an avalanche directory with required input files - for this task you can use :ref:`moduleIn3Utils:Initialize Project`
* Copy ``com1DFA/com1DFACfg.ini`` to ``com1DFA/local_com1DFACfg.ini`` and if desired change configuration settings
* If you are on a develop installation, make sure you have an updated compilation, see :ref:`developinstall:Setup AvaFrame`
* Run:
  ::

    python3 runCom9MoTVoellmy.py

<br/>


