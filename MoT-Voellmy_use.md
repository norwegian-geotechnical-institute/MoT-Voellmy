# Notes on Using MoT-Voellmy

</br>

### General features

MoT-Voellmy is NGI's model for calculating the run-out of snow avalanches. It is a _continuum model_, i.e., the avalanche is treated as a (fluid) body rather than as a point mass as in the block models. Thus, the length, width and thickness of the flow are variable both in space and time. The dynamical equations governing the flow are derived from the requirement that mass and momentum must be conserved (if erosion and deposition along the path and the external forces of gravity and friction are properly accounted for).

The "trick" behind numerical continuum models is to discretize space and time into sufficiently small finite time steps, line segments (in 1D), subareas (in 2D) or elementary volumes (in 3D). This procedure reduces the number of degrees of freedom from infinite to a finite number, so that the fundamental equations no longer are partial differential equations but (usually highly non-linear) algebraic equations. These can be solved numerically on a computer.

In many respects, MoT-Voellmy is similar to the widely used (commercial) code [RAMMS::AVALANCHE](https://doi.org/10.1016/j.coldregions.2010.04.005), but its numerics are simpler, which makes it very fast. Both programs do not solve the basic equations for the 3D velocity vector $(u, v, w)$ in 3D, but instead assume the flow to always follow the terrain surface. This makes it possible to reduce the equations to 2D equations for the slope-parallel _depth-averaged_ velocity $(\bar{u}, \bar{v})$ and the flow depth (or flow thickness) $h$, which is measured perpendicular to the terrain. Thus, one loses all information about the variability of $u$ and $v$ in the direction normal to the terrain. On the upside, the computational effort is reduced by some two to three orders-of-magnitude. Knowledge of the depth-averaged velocity is sufficient in most practical applications, like hazard mapping or dimensioning of mitigation measures.  
</br>
  

### Required input data

MoT-Voellmy is started from the command line with one argument, the name of the _run configuration file_ (RCF) as

    <path to the code>/MoT-Voellmy.2025-02-10 <RCF>

on Linux and macOS, or as

    <path to the code>\MoT-Voellmy.2025-02-10.exe <RCF>
	
on MS Windows.

In the RCF, the user specifies the spatial input data by listing the corresponding raster files and fixes the model parameters. The meaning and reasonable value ranges of these parameters are briefly explained below. For more detailed information and the mathematical formulation of the model components, see [Basic Equations and Numerical Methods in MoT-Voellmy](file:./MoT-Voellmy_use.md).

In the folder for MoT-Voellmy, there is a template RCF. All the data can be changed in a text editor, but all lines must be present and their sequence preserved (including comment lines starting with #). The keys on the left side must not be changed. The RCF should be renamed so that the user will still remember the purpose of the simulation some years later. It is recommended to save the RCFs in the folder or subfolder of the current project.  


#### Auxiliary information:

General information on the present run, not relevant for the program but for the user:

    MoT-Voellmy input file version          2024-09-10   (do not change!)
    Area of Interest                        Grasdalen
    UTM zone                                33N
    EPSG geodetic datum code                25833
    Run name                                Ryggfonn_2021-04-11_02

MoT-Voellmy version 2025-02-10 is designed for the input file version 2024-09-10 but will also accept the input file formats 2021-10-25 and 2020-06-23.
  
#### Flags:

These parameters take the values `yes` or `no`, except for the flag `Entrainment`, which has the options `none`, `RAMMS`, `TJEM`, `AvaFrame` or `GOEM`, the flag `Rheology`, which is presently limited to `Voellmy`, and the flag `Parameters` with the alternatives `constant` and `variable`. With `constant`, spatially constant values of $\mu$ and $k$ are assumed (see below under _Physical parameters_ and _Raster input data_).

    Centrifugal effects             {yes|no}
If "yes", the terrain-normal component of gravity is modified from $g \cos\theta$ to $\max(0, g \cos\theta + \kappa_{\mathbf u} \mathbf{u}^2)$, where $\kappa_{\mathbf{u}}$ is the terrain curvature in the flow direction $\mathbf{u}$. This accounts for the modification of dry friction due to centrifugal forces in curved terrain.

    Forest effects                  {yes|no}
Enables or disables braking of the flow by the forest and breaking of the forest by the flow. See below for the additional data files that must be supplied.

    Entrainment                     {none|RAMMS|TJEM|AvaFrame|GOEM}
Choose one of the five options. With `none`, entrainment is switched off. `RAMMS` enables the simple entrainment model available in early versions of RAMMS::AVALANCHE. `TJEM` is the shorthand for *Tangential-Jump Entrainment Model*. `AvaFrame` implements the basal erosion model of the AvaFrame computational module com1DFA, while `GOEM` stands for *Grigorian–Ostroumov Erosion Model*. Each of the four models has its specific requirements for setting (empirical) parameters and for additional raster data for snow depth and snow properties, to be discussed below.

    Bed strength profile            {constant|global|local}
In the case of the erosion models TJEM, GOEM and AvaFrame, where the entrainment rate depends on the snow-cover properties (shear strength, compressive strength and comminution energy, respectively), this property can be set to be constant or to vary linearly with depth in the snow cover. In the latter case, the gradient can either be specified globally by a single number, or locally through a raster file.

    Deposition                      {no}
In the present version of MoT-Voellmy, no deposition model is implemented.

    Evolving geometry               {yes|no}
If entrainment is enabled, the change of the terrain slope and curvature due to erosion and deposition is tracked and may induce self-channeling of the flow.

    Write velocity vectors          {yes|no}
If enabled, the code writes not only the modulus of the velocity, $\sqrt{{\bar{\mathbf{u}}}^2}$ but also its components $\bar{u}$, $\bar{v}$ to separate files.

    Write maximum pressure          {yes|no}
If enabled, the maximum impact pressure $p_{\mathrm{imp}} = \rho {\bar{\mathbf{u}}}^2$ in each cell attained over the entire simulation is written to file.

    Write instant. pressure         {yes|no}
If enabled, the instantaneous impact pressure $p_{\mathrm{imp}} = \rho \bar{\mathbf{u}}^2$ in each cell is written to file at each intermediate output time.


#### Physical parameters:

    Gravitational acceleration   (m/s^2)    9.81
This value can be changed to simulate extra-terrestrial or sub-marine landslides.

    Flow density                (kg/m^3)    250.0
The flow density, $\rho_f$, matters only for the impact pressure because the friction force is assumed proportional to $\rho_f$.

    Bed density                 (kg/m^3)    140.0
At present, the bed (snow-cover) density $\rho_b$ is used for converting the erosion rate ($\mathrm{kg\,m}^{-2}\,\mathrm{s}^{-1}$) to an erosion speed ($\mathrm{m\,s}^{-1}$).

    Deposit density             (kg/m^3)    450.0
This parameter has been introduced in anticipation of a deposition model; it is not used in MoT-Voellmy version 2025-02-10.

    Dry-friction coefficient         (-)    0.35
    Turbulent drag coefficient       (-)    0.0023
If the flag `Parameters` is set to `constant`, the global values for $\mu$ and $k$ must be set here. Else, raster file names containing the spatially variable values are to be specified instead.

    Effective drag height            (m)    0.0
This parameter modifies the flow-depth dependence of the drag term. If set to a very small positive value $h_{\mathrm{drag}} \ll 1\,\mathrm{m}$, the model behaves like the PCM model with $M/D = k/h_{\mathrm{drag}}$. The original Voellmy formula is recovered in the limit $h_{\mathrm{drag}} \rightarrow \infty$, but can also be enforced by setting $h_{\mathrm{drag}} \leq 0$.

    Passive earth-pressure coeff.    (-)    1.0
Values &lt; 1.0 can be used to stabilize the simulation in very rugged terrain by reducing the hydrostatic pressure between cells. Keep it at 1.0 whenever possible.

    Tree drag coefficient           (-)     1.0
This constant multiplies the velocity-dependent part of the force exerted by the trees. It is recommended to leave it at its default value, 1.

    Modulus of rupture            (MPa)    50.0
This parameter measures the average strength of trees against breaking. Measured values vary significantly between tree species and with the tree age and the season. For beech, 50 MPa is a reasonable average value.

    Forest decay coefficient     (m/s)      0.15
Together with the tree diameter, this coefficient determines how quickly trees fall down (and the forest loses its braking effect) once the avalanche has broken them.

    Erosion coefficient            (-)      0.01
In the RAMMS erosion model, this is the proportionality coefficient between avalanche flow speed and erosion speed. SLF suggested values in the range 0.1--1, but basal erosion is obtained only with values in the range 0.001--0.01.</br>
In the Grigorian–Ostroumov model (GOEM), the dynamic component of the erosive pressure is proportional to this coefficient and the stagnation pressure. Reasonable values are expected to be in the range $10^{-3}$--$10^{-2}$, but Grigorian and Ostroumov used much larger values.</br>
This parameter is not used in the TJEM and AvaFrame erosion model.

    Bed friction coefficient      (-)      0.25
If the bed strength profile is declared constant, MoT-Voellmy sets this parameter to 0. If the strength profile is declared `global`, the strength is assumed to increase with depth $b_0 - z$ according to $\tau(z) = \tau(b_0) + \mu_s \rho_b (b_0 - z) g \cos\theta$. In the case `local`, the values are read from a raster file, the name of which must be given on this line.


#### Numerical parameters

    Output filename root                    02
All generated output files will be given this prefix so different runs can be distinguished. It may also be of the form `<subfolder>/<some prefix>` to place the files in a separate subfolder.

    Output format                           ESRI_ASCII_Grid
The output raster files can be written either in ESRI ASCII Grid format (bulky but human readable) or in <code>BinaryTerrain</code> (version 1.3) format (faster and more compact). Both formats are understood by QGIS and ArcGIS.

    Simulation time                  (s)    180.0
This value specifies after which simulated time the computation is stopped, even though some parts of the flow (typically in the tail) may still be moving.

    Minimum time step                (s)    0.0001
    Maximum time step                (s)    0.5000
    Initial CFL number               (-)    0.8
In the time discretization, MoT-Voellmy uses variable timesteps to balance numerical stability (short timesteps) against efficiency (long timesteps). The timestep must be shorter than the time it takes a wave to traverse a computational cell. The CFL number specifies how close to this theoretical limit the time step is allowed to be; this number must be &lt; 1.0. The code keeps it, however, between the minimum and maximum time step chosen by the user. If the time step chosen by the code falls below the minimum value, the simulation is considered divergent and is stopped.

    Output interval                  (s)    1.0
Intermediate results are written to files at fixed intervals. At the end of the simulation, the maximum values, over the entire simulated interval, of $\bar{u}$, $\bar{v}$, $h$ and optionally of some other variables are written to files as well.

    Minimum flow depth               (m)    0.05
    Minimum speed                  (m/s)    0.1
Cells in which the flow depth or the speed fall below the specified values will be considered stopped. They may, however, be remobilized if mass from neighboring cells pushes them strongly enough to overcome friction.

    Momentum threshold          (kg m/s)    100.0
A stopping criterion for the entire flow simulation rather than for individual cells can be activated by choosing a suitable value for the threshold of momentum integrated over the entire flow mass, below which the flow is considered stopped in practice. The flow momentum is calculated by summing the _absolute values_ of all cell momenta. Reasonable values can be inferred from the release volume, the assumed flow density and a mean velocity threshold of the order 0.1–1 m/s. This would give 250,000 kg m/s for a small avalanche of 10,000 m^3^ with a density of 250 kg/m^3^ and a threshold velocity of 0.1 m/s.


#### Raster input data:

All the raster input data must be in ESRI ASCII Grid format and cover the same rectangular area at the same resolution. This means that some clipping and interpolation may be required when preparing the data for a simulation.

- <em>Digital terrain model</em> (DTM) over a rectangular area with a spatial resolution of typically 5–10 m (for snow avalanches). If the resolution is poorer, important terrain features are not properly resolved. If the resolution is much finer, the flow moves over terrain that is much more hummocky than the snow cover in wintertime. The grid must be regular, i.e., consist of square cells when projected onto a horizontal plane.
    <pre>Grid filename                           RA6157_dem.asc</pre>
    
- <em>Release area and fracture depth</em> are specified by the value 0 in non-release-area cells and the slab depth (in m, measured normal to the terrain) in the release area. There may be several disconnected release areas if one assumes they are released simultaneously. The fracture depth may vary inside a release area.
    <pre>Release depth filename                  RA6157_wf_h.asc</pre>

Depending on the desired simulation, additional raster files may need to be specified. All these files must be congruent with the digital terrain model, i.e., have the same number of rows and columns, the same cell size and lower-left corner coordinates and the same number indicating missing data:

- <em>Friction parameters $\mu$ and $k$</em>: If the friction parameters are assumed to vary spatially variable (reflecting spatial variation of snow properties), they must be specified in two raster files. The dimensionless parameter $k$ relates to the traditional Voellmy parameter $\xi$ with units m/s² as $k = g/\xi$.
    <pre>Dry-friction coefficient         (-)    RA6157_mu.asc
Turbulent drag coefficient       (-)    RA6157_k.asc</pre>
In addition, set `Parameters variable`
  
- <em>Start velocity</em>: To start a simulation with the flow mass already in motion, one can specify the velocity components $\bar{u}(t=0)$ and $\bar{v}(t=0)$ (in m/s) in two raster files. If starting from rest, use hyphens instead of file names.
    <pre>Start velocity u filename               -
Start velocity v filename               -</pre>

- If <em>entrainment of bed material</em> shall be simulated, two raster files specifying the local initial depth $b_0$ (in m) and shear strength $\tau_c$ (in Pa) of the erodible material must be supplied.
    <pre>Bed depth filename                      RA6157_b0.asc
Bed shear strength filename             RA6157_tauc.asc</pre>

- <em>Interaction of the flow with forest</em> requires two raster files giving the local values of $n\cdot D$ and $D$, where $n$ in 1/m² is the local tree density and $D$ (in m) is the local average tree diameter. In addition a global value for the <em>modulus of rupture</em> (MoR, in MPa, default 50 MPa) characterizing the average rupture strength of the wood must be set in the RCF. Set this parameter to a very high value if you wish to disable breaking or uprooting of the trees.
    <pre>Forest density filename                 RA6157_nD.asc
Tree diameter filename                  RA6157_tD.asc</pre></br>

Here is an example of a complete RCF:
<pre># Run information
#
MoT-Voellmy input file version          2024-09-10
Area of Interest                        Grasdalen
UTM zone                                33N
EPSG geodetic datum code                25833
Run name                                Ryggfonn_2021-04-11_02
#
# File names
#
Grid filename                           /data/Tests/Ryggfonn/Input/dhm.asc
Release depth filename                  /data/Tests/Ryggfonn/Input/h0.asc
Bed depth filename                      -
Bed shear strength filename             -
Forest density filename                 -
Tree diameter filename                  -
Start velocity u filename               -
Start velocity v filename               -
Output filename root                    /data/Tests/Ryggfonn/Run_02/02
Output format                           ESRI_ASCII_Grid
#
# Physical parameters
#
Gravitational acceleration   (m/s^2)    9.81
Flow density                (kg/m^3)    250.0
Bed density                 (kg/m^3)    140.0
Deposit density             (kg/m^3)    450.0
Rheology                                Voellmy
Parameters                              constant
Dry-friction coefficient         (-)    0.40
Turbulent drag coefficient       (-)    0.001
Effective drag height            (m)    3.0
Centrifugal effects                     yes
Passive earth-pressure coeff.    (-)    1.0
#
Forest effects                          no
Tree drag coefficient            (-)    1.0
Modulus of rupture             (MPa)    50.0
Forest decay coefficient       (m/s)    0.15
#
Entrainment                             none
Erosion coefficient              (-)    0.0
Bed strength profile                    global
Bed friction coefficient         (-)    0.25
Deposition                              no
Evolving geometry                       no
#
# Numerical parameters
#
Simulation time                  (s)    100.0
Minimum time step                (s)    0.001
Maximum time step                (s)    0.2
Output interval                  (s)    1.0
Write velocity vectors                  no
Write maximum pressure                  yes
Write instant. pressure                 no
Minimum flow depth               (m)    0.01
Minimum speed                  (m/s)    0.01
Momentum threshold          (kg m/s)    100.0
Initial CFL number               (-)    0.8

</pre></br>


## Running a simulation

There are different ways a simulation with MoT-Voellmy can be launched:

1. The standard way is manually from a terminal (in MS Windows: `cmd.exe` or `PowerShell.exe`). This method is simple and fast, with minimal overhead.</br>
If the directory containing the MoT-Voellmy executable is in the user's `PATH`, the simulation can be started by the command
    <pre>MoT-Voellmy-linux.&lt;version&gt; &lt;path to RCF&gt;&lt;RCF&gt;</pre>
on Linux,
    <pre>MoT-Voellmy-macos.&lt;version&gt; &lt;path to RCF&gt;&lt;RCF&gt;</pre>
on macOS systems, and with
    <pre>MoT-Voellmy-win64.&lt;version&gt;.exe; &lt;path to RCF&gt;&lt;RCF&gt;</pre>
on MS Windows.

2. This repository also features a Jupyter Notebook, which allows to run MoT-Voellmy in specific cells within the Notebook and to document all work in the Notebook as well. In principle, complete project reports including interactive elements, maps and animations can be produced as Notebooks. The user interacts with the Notebook through a web browser. This method requires Jupyter Notebook or JupyterLab to be installed on the local machine or on a server, together with the Jupyter Server application and so-called kernels, which interactively execute the code entered by the user. The notebook in this repository is set up to run MoT-Voellmy on the test case Ryggfonn, but it is straightforward to  change this to another example provided in the repository or a case specified by the user.

3. A third approach is to use an application written for the user's preferred GIS. Such an application was written at NGI for ArcMap several years ago and is being updated to the newest versions of MoT-Voellmy and ArcGIS Pro. The project [AvaFrame](https://docs.avaframe.org/) develops a similar graphical user interface for the open-source program [QGIS](https://qgis.org/) and the computational modules offered in the AvaFrame framework. Work on integrating MoT-Voellmy in AvaFrame is ongoing and will make a QGIS graphical user interface (GUI) available to users of MoT-Voellmy.</br>
Such  GUIs typically include functions to
    - select the simulation area,
    - specify the release area,
    - potentially create further raster input files like the erodible snow depth, the snow-cover shear strength, as well as the number per unit area and the mean diameter of trees,
    - put together the RCF in a form,
    - run MoT-Voellmy with the RCF,
    - visualize the simulation results in maps, and
    - support the generation of suitably symbolized maps for reporting.  
</br>


## Visualization of results

MoT-Voellmy prints information on the progress of the simulation to the screen (unless it is running within NAKSIN 4). This allows the user to judge whether there are any irregularities that require closer scrutiny.

The simulation results are output as raster files:

- The file format can be chosen in the run configuration file as `ESRI_ASCII_Grid` or `BinaryTerrain 1.3`. Both formats can be read by QGIS or ArcGIS. BinaryTerrain 1.3 is read and written somewhat faster and the files are more compact, but ESRI ASCII Grid files can usually be compressed better than BinaryTerrain 1.3 files.

- Two types of files can be written to disk:
    - <em>Time slices</em> showing the values of one field variable at a given time, e.g., $h(x_i, y_j, t_*)$, the flow depth at time $t_*$ in every grid cell $(i,j)$
    - <em>Maximum values</em> of one field variable over the entire simulation, i.e., $\max_{t_n} h(x_i,y_j, t_n)$ at each cell $(i,j)$

- The initial conditions $h(x_i, y_j, t_0)$ and $s(x_i, y_j, t_0)$ and the maximum fields $h_{\mathrm{max}}(x_i, y_j)$ and $s_{\mathrm{max}}(x_i,y_j)$ are always written. $s = \sqrt{\mathbf{u}^2}$ is the speed.

- In addition, the user can specify which other variables among the velocity components $u$, $v$, the pressure $p$, the erodible snow-cover depth $b$ and the forest stand density $n$ should be output.

- Time slices are written at times $k \Delta t_{\mathrm{out}}$, where $\Delta t_{\mathrm{out}}$ is the output interval in seconds. If this interval is smaller than the minimum timestep, time slices are written after each timestep. If $\Delta_{\mathrm{out}}$ is larger than the simulation time $t_{\mathrm{max}}$, only the initial, final and maximum values are written to file.

How the results are visualized, will depend on the purpose of the simulation and personal preferences. Besides the full-featured geographic information systems (GIS) like QGIS or ArcGIS Pro, simpler programs like Python scripts using `matplotlib` can also produce basic plots and even animations quickly and easily.

Often it is helpful to load at least a topographic background map (if available), the DTM, the initial condition $h(x_i, y_j, t_0)$ and the maxima of $h$ and $s$. If no background map is available, contour lines can be created in QGIS from the DTM. (Raster – Extract – Contour, typically thin 10 m contour lines and thicker 100 m contour lines help the user grasp the terrain properties). To facilitate analysis of the results, suitable choices for the color palettes and the values of layer transparency should be found and saved as styles.

In the folder `Test_cases`, there is a file `Ryggfonn.qgz`. Open it in QGIS to visualize your simulation results (there may be a warning that some files were not found – press the button `Browse` to find them. You can add additional layers to view intermediate snapshots of the avalanche descent.

If one wishes to explore the model further, one will modify the RCF—*saving it under a different name*—and run MoT-Voellmy with it. For example, to find out how the coefficient <em>k</em> affects the speed and run-out, one chooses another constant value for <em>k</em> or replaces the values for <em>µ</em> and <em>k</em> by the names of the files contained in the same folder.

MoT-Voellmy also allows modeling snow entrainment along the path. To see this in action, one can enter the file names `RGF_1987-01-28_b0.asc`, `RGF_1987-01-28_tauc.asc` in lines 13 and 14 of the RCF, and set the value of `Entrainment` to `TJEM`. To study the effect of the fracture depth, one may modify the raster file `RGF_1987-01-28_h0.asc`, e.g., by loading the file into an editor and replacing all instances of "1.0" by some other fracture depth.