
* Create & run LaMEM models from julia
* Available functions

* Available functions

[GitHub](https://github.com/JuliaGeodynamics/LaMEM.jl "View the repository on GitHub")[](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/main/docs/src/LaMEM_ModelFunctions.md "Edit source on GitHub")[](# "Settings")

[List functions](#List-functions)[](#List-functions "Permalink")
================================================================

These are all the functions that are provided for the LaMEM Julia Setup interface

[`LaMEM.LaMEM_Model.BCBlock`](#LaMEM.LaMEM_Model.BCBlock) — Type

    LaMEM boundary condition `BCBlock` object

* `npath::Int64`: Number of path points of Bezier curve (path-points only!)
    
* `theta::Vector{Float64}`: # Orientation angles at path points (counter-clockwise positive)
    
* `time::Vector{Float64}`: Times at path points
    
* `path::Vector{Float64}`: Path points x-y coordinates
    
* `npoly::Int64`: Number of polygon vertices
    
* `poly::Vector{Float64}`: Polygon x-y coordinates at initial time
    
* `bot::Float64`: Polygon bottom coordinate
    
* `top::Float64`: Polygon top coordinate
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/BoundaryConditions.jl#L97-L101)

[`LaMEM.LaMEM_Model.BoundaryConditions`](#LaMEM.LaMEM_Model.BoundaryConditions) — Type

    Structure that contains the LaMEM boundary conditions information.

* `noslip::Vector{Int64}`: No-slip boundary flag mask (left right front back bottom top)
    
* `open_top_bound::Int64`: Stress-free (free surface/infinitely fast erosion) top boundary flag
    
* `temp_top::Float64`: Constant temperature on the top boundary
    
* `temp_bot::Float64`: Constant temperature on the bottom boundary
    
* `exx_num_periods::Int64`: number intervals of constant background strain rate (x-axis)
    
* `exx_time_delims::Vector{Float64}`: time delimiters (one less than number of intervals, not required for one interval)
    
* `exx_strain_rates::Vector{Float64}`: strain rates for each interval
    
* `eyy_num_periods::Int64`: eyy_num_periods
    
* `eyy_time_delims::Vector{Float64}`: eyy_time_delims
    
* `eyy_strain_rates::Vector{Float64}`: eyy_strain_rates
    
* `exy_num_periods::Int64`: exy_num_periods
    
* `exy_time_delims::Vector{Float64}`: exy_time_delims
    
* `exy_strain_rates::Vector{Float64}`: exy_strain_rates
    
* `exz_num_periods::Int64`: exz_num_periods
    
* `exz_time_delims::Vector{Float64}`: exz_time_delims
    
* `exz_strain_rates::Vector{Float64}`: exz_strain_rates
    
* `eyz_num_periods::Int64`: eyz_num_periods
    
* `eyz_time_delims::Vector{Float64}`: eyz_time_delims
    
* `eyz_strain_rates::Vector{Float64}`: eyz_strain_rates
    
* `bg_ref_point::Vector{Float64}`: background strain rate reference point (fixed)
    
* `VelocityBoxes::Vector{VelocityBox}`: List of added velocity boxes
    
* `BCBlocks::Vector{BCBlock}`: List of added Bezier blocks
    
* `VelCylinders::Vector{VelCylinder}`: List of added velocity cylinders
    
* `bvel_face::Union{Nothing, String}`: Face identifier (Left; Right; Front; Back; CompensatingInflow)
    
* `bvel_face_out::Union{Nothing, Int64}`: Velocity on opposite side: -1 for inverted velocity; 0 for no velocity; 1 for the same direction of velocity
    
* `bvel_bot::Union{Nothing, Float64}`: Bottom coordinate of inflow window
    
* `bvel_top::Union{Nothing, Float64}`: Top coordinate of inflow window
    
* `velin_num_periods::Union{Nothing, Int64}`: Number of periods when velocity changes (Optional)
    
* `velin_time_delims::Union{Nothing, Vector}`: Change velocity at 2 and 5 Myrs (one less than number of intervals, not required for one interval) (Optional)
    
* `bvel_velin::Union{Nothing, Vector}`: inflow velocity for each time interval(Multiple values required if velin_num_periods>1)
    
* `bvel_velout::Union{Nothing, Float64}`: outflow velocity (if not specified, computed from mass balance)
    
* `bvel_relax_d::Union{Nothing, Float64}`: vert.distance from bvel_bot and bvel_top over which velocity is reduced linearly
    
* `bvel_velbot::Union{Nothing, Int64}`: bottom inflow velocity for use with bvel_face=CompensatingInflow
    
* `bvel_veltop::Union{Nothing, Int64}`: top inflow velocity for use with bvel_face=CompensatingInflow
    
* `bvel_temperature_inflow::Union{Nothing, String}`: bvel_temperature_inflow: Thermal age of the plate, which can be constant if set to Fixed_thermal_age or Constant_T_inflow (Temperature of the inflow material is constant everywhere)
    
* `bvel_thermal_age::Union{Nothing, Float64}`: In dimensional unit. If the user specify this value, he needs to specify the temperature of the mantle and top as well
    
* `bvel_temperature_mantle::Union{Nothing, Float64}`: In dimensional unit. Temperature of the mantle
    
* `bvel_temperature_top::Union{Nothing, Float64}`: In dimensional unit. temperature of the top
    
* `bvel_temperature_constant::Union{Nothing, Float64}`: Constant temperature inflow.
    
* `bvel_num_phase::Union{Nothing, Int64}`: Imposes a stratigraphy of phase injected in the inflow boundary \[if undefined, it uses the phase close to the boundary\]
    
* `bvel_phase::Union{Nothing, Vector{Int64}}`: phase number of inflow material \[if undefined, it uses the phase close to the boundary\] from bottom to top
    
* `bvel_phase_interval::Union{Nothing, Vector{Float64}}`: Depth interval of injection of the phase (the interval is defined by num_phase+1 coordinates). e.g. \[-120 -100 -10 0 \]
    
* `open_bot_bound::Union{Nothing, Int64}`: # Permeable lower boundary flag
    
* `permeable_phase_inflow::Union{Nothing, Int64}`: Phase of the inflow material from the bottom (The temperature of the inflow phase it is the same of the bottom boundary) in case of open_bot_bound=1
    
* `fix_phase::Union{Nothing, Int64}`: fixed phase (no-flow condition)
    
* `fix_cell::Union{Nothing, Int64}`: fixed cells (no-flow condition)
    
* `fix_cell_file::Union{Nothing, String}`: fixed cells input file (extension is .xxxxxxxx.dat)
    
* `temp_bot_num_periods::Union{Nothing, Int64}`: How many periods with different temp_bot do we have?
    
* `temp_bot_time_delim::Union{Nothing, Vector{Float64}}`: At which time do we switch from one to the next period?
    
* `Plume_InflowBoundary::Union{Nothing, Int64}`: # have a plume-like inflow boundary @ bottom
    
* `Plume_Type::Union{Nothing, String}`: Type of plume inflow boundary.
    
    * `"Inflow_type"` or
    * `"Pressure_type"` (circular) or
    * `"Permeable_Type"` which combines the open bot boundary with the plume boundary condition (the option herein listed overwrites open_bot, so do not activate that)

* `Plume_Dimension::Union{Nothing, String}`: 2D or 3D (circular)
    
* `Plume_areaFrac::Union{Nothing, Float64}`: how much of the plume is actually in the model. This usually 1 (default) but lower if the plume is in a corner of a symmetric setup and matters for the outflow
    
* `Plume_Phase::Union{Nothing, Int64}`: phase of plume material
    
* `Plume_Depth::Union{Nothing, Float64}`: # depth of provenience of the plume (i.e. how far from the bottom of the model the plume source is)
    
* `Plume_Mantle_Phase::Union{Nothing, Int64}`: # Astenosphere phase (if the inflow occurs outside the plume radius)
    
* `Plume_Temperature::Union{Nothing, Float64}`: # temperature of inflow plume
    
* `Plume_Inflow_Velocity::Union{Nothing, Float64}`: # Inflow velocity (not required if Pressure_Type) in cm/year if using GEOunits
    
* `Plume_VelocityType::Union{Nothing, String}`: `"Gaussian"` or `"Poiseuille"`
    
* `Plume_Center::Union{Nothing, Vector{Float64}}`: # \[X,Y\] of center (2nd only in case of 3D plume)
    
* `Plume_Radius::Union{Nothing, Float64}`: # Width/Radius of plume
    
* `Plume_Phase_Mantle::Union{Nothing, Int64}`: # Inflow phase. If the velocity happens to be positive in the domain, the inflow material has a constant phase and the temperature of the bottom
    
* `pres_top::Union{Nothing, Float64}`: Pressure on the top boundary
    
* `pres_bot::Union{Nothing, Float64}`: Pressure on the bottom boundary
    
* `init_pres::Union{Nothing, Int64}`: pressure initial guess flag; linear profile between pres_top and pres_bot in the unconstrained cells
    
* `init_temp::Union{Nothing, Int64}`: temperature initial guess flag; linear profile between temp_top and temp_bot
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/BoundaryConditions.jl#L251-L256)

[`LaMEM.LaMEM_Model.Dike`](#LaMEM.LaMEM_Model.Dike) — Type

    Defines the properties related to inserting dikes

* `ID::Int64`: Material phase ID
    
* `Mf::Float64`: value for dike/magma- accommodated extension, between 0 and 1, in the front of the box, for phase dike
    
* `Mc::Float64`: \[optional\] value for dike/magma- accommodate extension, between 0 and 1, for dike phase; M is linearly interpolated between Mf & Mc and Mc & Mb, if not set, Mc default is set to -1 so it is not used
    
* `y_Mc::Union{Nothing, Float64}`: \[optional\], location for Mc, must be between front and back boundaries of dike box, if not set, default value to 0.0, but not used
    
* `Mb::Union{Nothing, Float64}`: value for dike/magma-accommodated extension, between 0 and 1, in the back of the box, for phase dike
    
* `PhaseID::Union{Nothing, Int64}`: Phase ID
    
* `PhaseTransID::Union{Nothing, Int64}`: Phase transition ID
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Materials.jl#L604-L608)

[`LaMEM.LaMEM_Model.FreeSurface`](#LaMEM.LaMEM_Model.FreeSurface) — Type

    Structure that contains the LaMEM free surface information.

* `surf_use::Int64`: Free surface activation flag
    
* `surf_corr_phase::Int64`: air phase ratio correction flag (phases in an element that contains are modified based on the surface position)
    
* `surf_level::Union{Nothing, Float64}`: initial level of the free surface
    
* `surf_air_phase::Union{Nothing, Int64}`: phase ID of sticky air layer
    
* `surf_max_angle::Float64`: maximum angle with horizon (smoothed if larger)
    
* `surf_topo_file::String`: initial topography file (redundant)
    
* `erosion_model::Int64`: erosion model \[0-none (default), 1-infinitely fast, 2-prescribed rate with given level\]
    
* `er_num_phases::Int64`: number of erosion phases
    
* `er_time_delims::Vector{Float64}`: erosion time delimiters (one less than number)
    
* `er_rates::Vector{Float64}`: constant erosion rates in different time periods
    
* `er_levels::Vector{Int64}`: levels above which we apply constant erosion rates in different time periods
    
* `sediment_model::Int64`: sedimentation model \[0-none (dafault), 1-prescribed rate with given level, 2-cont. margin\]
    
* `sed_num_layers::Int64`: number of sediment layers
    
* `sed_time_delims::Vector{Float64}`: sediment layers time delimiters (one less than number)
    
* `sed_rates::Vector{Float64}`: sediment rates in different time periods
    
* `sed_levels::Vector{Float64}`: levels below which we apply constant sediment rates in different time periods
    
* `sed_phases::Vector{Int64}`: sediment layers phase numbers in different time periods
    
* `marginO::Vector{Float64}`: lateral coordinates of continental margin - origin
    
* `marginE::Vector{Float64}`: lateral coordinates of continental margin - 2nd point
    
* `hUp::Float64`: up dip thickness of sediment cover (onshore)
    
* `hDown::Float64`: down dip thickness of sediment cover (off shore)
    
* `dTrans::Float64`: half of transition zone
    
* `Topography::Union{Nothing, GeophysicalModelGenerator.CartData}`: Topography grid
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/FreeSurface.jl#L7-L11)

[`LaMEM.LaMEM_Model.GeomBox`](#LaMEM.LaMEM_Model.GeomBox) — Type

    LaMEM geometric primitive `Box` object

* `phase::Int64`: phase
    
* `bounds::Vector{Float64}`: box bound coordinates: `left`, `right`, `front`, `back`, `bottom`, `top`
    
* `Temperature::Union{Nothing, String}`: optional: Temperature structure. possibilities: \[constant, linear, halfspace\]
    
* `cstTemp::Union{Nothing, Float64}`: required in case of \[`constant`\]: temperature value \[in Celcius in case of GEO units\]
    
* `topTemp::Union{Nothing, Float64}`: required in case of \[`linear,halfspace`\]: temperature @ top \[in Celcius in case of GEO units\]
    
* `botTemp::Union{Nothing, Float64}`: required in case of \[`linear,halfspace`\]: temperature @ top \[in Celcius in case of GEO units\]
    
* `thermalAge::Union{Nothing, Float64}`: required in case of \[`halfspace`\]: thermal age of lithosphere \[in Myrs if GEO units are used\]
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ModelSetup.jl#L207-L211)

[`LaMEM.LaMEM_Model.GeomCylinder`](#LaMEM.LaMEM_Model.GeomCylinder) — Type

    LaMEM geometric primitive `Cylinder` object

* `phase::Int64`: phase
    
* `radius::Float64`: radius of cylinder
    
* `base::Vector{Float64}`: center of base of cylinder
    
* `cap::Vector{Float64}`: center of cap of cylinder
    
* `Temperature::Union{Nothing, String}`: optional: Temperature structure. possibilities: \[constant\]
    
* `cstTemp::Union{Nothing, Float64}`: required in case of \[`constant`\]: temperature value \[in Celcius in case of GEO units\]
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ModelSetup.jl#L428-L432)

[`LaMEM.LaMEM_Model.GeomEllipsoid`](#LaMEM.LaMEM_Model.GeomEllipsoid) — Type

    LaMEM geometric primitive `Ellipsoid` object

* `phase::Int64`: phase
    
* `axes::Vector{Float64}`: semi-axes of ellipsoid in `x`, `y` and `z`
    
* `center::Vector{Float64}`: center of sphere
    
* `Temperature::Union{Nothing, String}`: optional: Temperature of the sphere. possibilities: \[constant, or nothing\]
    
* `cstTemp::Union{Nothing, Float64}`: required in case of \[constant\]: temperature value \[in Celcius in case of GEO units\]
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ModelSetup.jl#L159-L163)

[`LaMEM.LaMEM_Model.GeomHex`](#LaMEM.LaMEM_Model.GeomHex) — Type

    LaMEM geometric primitive `Hex` object to define hexahedral elements

* `phase::Int64`: phase
    
* `coord::Vector{Float64}`: `x`-`y`-`z` coordinates for each of 8 nodes (24 parameters) (counter)-clockwise for an arbitrary face, followed by the opposite face
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ModelSetup.jl#L322-L326)

[`LaMEM.LaMEM_Model.GeomLayer`](#LaMEM.LaMEM_Model.GeomLayer) — Type

    LaMEM geometric primitive `Layer` object

* `phase::Int64`: phase
    
* `top::Float64`: top of layer
    
* `bottom::Float64`: bottom of layer
    
* `cosine::Union{Nothing, Int64}`: optional: add a cosine perturbation on top of the interface (if 1)
    
* `wavelength::Union{Nothing, Float64}`: required if cosine: wavelength in x-direction
    
* `amplitude::Union{Nothing, Float64}`: required if cosine: amplitude of perturbation
    
* `Temperature::Union{Nothing, String}`: optional: Temperature structure. possibilities: \[constant, linear, halfspace\]
    
* `cstTemp::Union{Nothing, Float64}`: required in case of \[`constant`\]: temperature value \[in Celcius in case of GEO units\]
    
* `topTemp::Union{Nothing, Float64}`: required in case of \[`linear,halfspace`\]: temperature @ top \[in Celcius in case of GEO units\]
    
* `botTemp::Union{Nothing, Float64}`: required in case of \[`linear,halfspace`\]: temperature @ top \[in Celcius in case of GEO units\]
    
* `thermalAge::Union{Nothing, Float64}`: required in case of \[`halfspace`\]: thermal age of lithosphere \[in Myrs if GEO units are used\]
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ModelSetup.jl#L363-L367)

[`LaMEM.LaMEM_Model.GeomRidgeSeg`](#LaMEM.LaMEM_Model.GeomRidgeSeg) — Type

    LaMEM geometric primitive `RidgeSeg` object

* `phase::Int64`: phase
    
* `bounds::Vector{Float64}`: box bound coordinates: `left`, `right`, `front`, `back`, `bottom`, `top`
    
* `ridgeseg_x::Vector{Float64}`: coordinate order: left, right \[can be different for oblique ridge\]
    
* `ridgeseg_y::Vector{Float64}`: coordinate order: front, back \[can be different for oblique ridge\]
    
* `Temperature::String`: initial temperature structure \[ridge must be set to `halfspace_age` –\> setTemp=4\]
    
* `topTemp::Float64`: required in case of \[`linear,halfspace`\]: temperature @ top \[in Celcius in case of GEO units\]
    
* `botTemp::Float64`: required in case of \[`linear,halfspace`\]: temperature @ top \[in Celcius in case of GEO units\]
    
* `age0::Float64`: minimum age of seafloor at ridge \[in `Myr` in case of GEO units\]
    
* `maxAge::Union{Nothing, Float64}`: \[optional\] parameter that indicates the maximum thermal age of a plate
    
* `v_spread::Union{Nothing, Float64}`: \[optional\] parameter that indicates the spreading velocity of the plate; if not defined it uses bvel_velin specified elsewhere
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ModelSetup.jl#L260-L264)

[`LaMEM.LaMEM_Model.GeomSphere`](#LaMEM.LaMEM_Model.GeomSphere) — Type

    LaMEM geometric primitive `sphere` object

* `phase::Int64`: phase
    
* `radius::Float64`: radius of sphere
    
* `center::Vector{Float64}`: center of sphere
    
* `Temperature::Union{Nothing, String}`: optional: Temperature of the sphere. possibilities: \[constant, or nothing\]
    
* `cstTemp::Union{Nothing, Float64}`: required in case of \[constant\]: temperature value \[in Celcius in case of GEO units\]
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ModelSetup.jl#L110-L114)

[`LaMEM.LaMEM_Model.Grid`](#LaMEM.LaMEM_Model.Grid) — Type

    Structure that contains the LaMEM grid information

* `nmark_x::Int64`: number of markers/element in x-direction
    
* `nmark_y::Int64`: number of markers/element in y-direction
    
* `nmark_z::Int64`: number of markers/element in x-direction
    
* `nel_x::Vector{Int64}`: number of elements in x-direction
    
* `nel_y::Vector{Int64}`: number of elements in y-direction
    
* `nel_z::Vector{Int64}`: number of elements in z-direction
    
* `coord_x::Vector{Float64}`: coordinates in x-direction
    
* `coord_y::Vector{Float64}`: coordinates in y-direction
    
* `coord_z::Vector{Float64}`: coordinates in z-direction
    
* `nseg_x::Int64`: number of segments in x-direction (if we employ variable grid spacing in x-direction)
    
* `nseg_y::Int64`: number of segments in y-direction (if we employ variable grid spacing in y-direction)
    
* `nseg_z::Int64`: number of segments in z-direction (if we employ variable grid spacing in z-direction)
    
* `bias_x::Vector{Float64}`: bias in x-direction (if we employ variable grid spacing in x-direction)
    
* `bias_y::Vector{Float64}`: bias in y-direction (if we employ variable grid spacing in y-direction)
    
* `bias_z::Vector{Float64}`: bias in z-direction (if we employ variable grid spacing in z-direction)
    
* `Grid::GeophysicalModelGenerator.LaMEM_grid`: Contains the LaMEM Grid object
    
* `Phases::Array{Int32}`: Phases; 3D phase information
    
* `Temp::Array{Float64}`: Temp; 3D phase information
    

**Example 1**

    julia> d=LaMEM.Grid(coord_x=[0.0, 0.7, 0.8, 1.0], bias_x=[0.3,1.0,3.0], nel_x=[10,4,2])
    LaMEM grid with 1D refinement: 
      nel         : ([10, 4, 2], [16], [16])
      marker/cell : (3, 3, 3)
      x           ϵ [0.0, 0.7, 0.8, 1.0], bias=[0.3, 1.0, 3.0], nseg=3, Δmin=0.025000000000000022, Δmax=0.1499999999999999
      y           ϵ [-10.0 : 0.0]
      z           ϵ [-10.0 : 0.0]

**Example 2**

    julia> d=LaMEM.Grid(nel=(10,20))
    LaMEM grid with constant Δ: 
      nel         : ([10], [1], [20])
      marker/cell : (3, 3, 3)
      x           ϵ [-10.0 : 10.0]
      y           ϵ [-10.0 : 0.0]
      z           ϵ [-10.0 : 0.0]

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Grid.jl#L4-L32)

[`LaMEM.LaMEM_Model.Materials`](#LaMEM.LaMEM_Model.Materials) — Type

    Structure that contains the material properties in the current simulation

* `Phases::Vector{Phase}`: Different Materials implemented
    
* `SofteningLaws::Vector{Softening}`: Softening laws implemented
    
* `PhaseTransitions::Vector{PhaseTransition}`: Internal Phase Transitions (that change the ID of markers) implemented
    
* `Dikes::Vector{Dike}`: Dikes implemented (mostly for MOR simulations)
    
* `PhaseAggregates::Vector{PhaseAggregate}`: Phase aggregates (combines different phases such as upper_lower crust into one for visualization purposes)
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Materials.jl#L662-L666)

[`LaMEM.LaMEM_Model.Model`](#LaMEM.LaMEM_Model.Model) — Type

    Model

Structure that holds all the information to create a LaMEM input file

* `Scaling::Scaling`: Scaling parameters
    
* `Grid::Grid`: LaMEM Grid
    
* `Time::Any`: Time options
    
* `FreeSurface::Any`: Free surface options
    
* `BoundaryConditions::Any`: Boundary conditions
    
* `SolutionParams::Any`: Global solution parameters
    
* `Solver::Any`: Solver options and optional PETSc options
    
* `ModelSetup::Any`: Model setup
    
* `Output::Any`: Output options
    
* `PassiveTracers::Any`: Passive tracers
    
* `Materials::Any`: Material parameters for each of the phases
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Model.jl#L9-L15)

[`LaMEM.LaMEM_Model.Model`](#LaMEM.LaMEM_Model.Model-Tuple) — Method

    Model(args...)

Allow to define a model setup by specifying some of the basic objects

**Example**

    julia> d = Model(Grid(nel=(10,1,20)), Scaling(NO_units()))
    LaMEM Model setup
    |
    |-- Scaling             :  GeoParams.Units.GeoUnits{GeoParams.Units.NONE}
    |-- Grid                :  nel=(10, 1, 20); xϵ(-10.0, 10.0), yϵ(-10.0, 0.0), zϵ(-10.0, 0.0) 
    |-- Time                :  nstep_max=50; nstep_out=1; time_end=1.0; dt=0.05
    |-- Boundary conditions :  noslip=[0, 0, 0, 0, 0, 0]
    |-- Solution parameters :  
    |-- Solver options      :  direct solver; superlu_dist; penalty term=10000.0
    |-- Model setup options :  Type=geom; 
    |-- Output options      :  filename=output; pvd=1; avd=0; surf=0
    |-- Materials           :  1 phases;  
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Model.jl#L93-L116)

[`LaMEM.LaMEM_Model.Model`](#LaMEM.LaMEM_Model.Model-Tuple{}) — Method

    Model(;
        Scaling=Scaling(GEO_units()),
        Grid=Grid(), 
        Time=Time(),
        FreeSurface=FreeSurface(),
        BoundaryConditions=BoundaryConditions(),
        SolutionParams=SolutionParams(),
        Solver=Solver(),
        ModelSetup=ModelSetup(),
        Output=Output(),
        PassiveTracers=PassiveTracers(),
        Materials=Materials()
        )

Creates a LaMEM Model setup.

* `Scaling::Scaling`
    
* `Grid::Grid`
    
* `Time::Any`
    
* `FreeSurface::Any`
    
* `BoundaryConditions::Any`
    
* `SolutionParams::Any`
    
* `Solver::Any`
    
* `ModelSetup::Any`
    
* `Output::Any`
    
* `PassiveTracers::Any`
    
* `Materials::Any`
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Model.jl#L53-L72)

[`LaMEM.LaMEM_Model.ModelSetup`](#LaMEM.LaMEM_Model.ModelSetup) — Type

    Structure that contains the LaMEM Model Setup and Advection options

* `msetup::String`: Setup type - can be `geom` (phases are assigned from geometric primitives, using `add_geom!(model, ...)`), `files` (from julia input), `polygons` (from geomIO input, which requires `poly_file` to be specified)
    
* `rand_noise::Int64`: add random noise to the particle location
    
* `rand_noiseGP::Int64`: random noise flag, subsequently applied to geometric primitives
    
* `bg_phase::Int64`: background phase ID
    
* `save_mark::Int64`: save marker to disk flag
    
* `mark_load_file::String`: marker input file (extension is .xxxxxxxx.dat), if using `msetup`=`files`
    
* `mark_save_file::String`: marker output file (extension is .xxxxxxxx.dat)
    
* `poly_file::String`: polygon geometry file (redundant), if using `msetup`=`polygons`
    
* `temp_file::String`: initial temperature file (redundant), if not set on markers
    
* `advect::String`: advection scheme; options=`none` (no advection); `basic` (Euler classical implementation \[default\]); `Euler` (Euler explicit in time); `rk2` (Runge-Kutta 2nd order in space)
    
* `interp::String`: velocity interpolation scheme; options = `stag` (trilinear interpolation from FDSTAG points), `minmod` ( MINMOD interpolation to nodes, trilinear interpolation to markers + correction), `stagp` ( STAG_P empirical approach by T. Gerya)
    
* `stagp_a::Float64`: STAG_P velocity interpolation parameter
    
* `mark_ctrl::String`: marker control type; options are `subgrid` (default; marker control enforced over fine scale grid), `none` (none), `basic` (AVD for cells + corner insertion), and `avd` (pure AVD for all control volumes)
    
* `nmark_lim::Vector{Int64}`: min/max number per cell (marker control)
    
* `nmark_avd::Vector{Int64}`: x-y-z AVD refinement factors (avd marker control)
    
* `nmark_sub::Int64`: max number of same phase markers per subcell (subgrid marker control)
    
* `geom_primitives::Vector`: Different geometric primitives that can be selected if we `msetup``=`geom`; see`GeomSphere`
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ModelSetup.jl#L7-L11)

[`LaMEM.LaMEM_Model.Multigrid`](#LaMEM.LaMEM_Model.Multigrid) — Type

    Structure that has info about setting up multigrid for LaMEM

* `nel::Tuple{Int64, Int64, Int64}`: Number of elements at the fine level

* `levels::Int64`: Number of levels

* `smooth::Int64`: number of smoothening steps per level

* `smooth_jacobi_factor::Float64`: factor for jacbi smoothener oer level

* `smoother::String`: smoother used at every level

* `coarse_ksp::String`: coarse grid ksp type preonly or fgmres

* `coarse_pc::String`: coarse grid pc type \["superlu_dist", "mumps", "gamg", "telescope","redundant"\]

* `coarse_coarse_pc::String`: coarse coarse grid solver in case we use redundant or telescope coarse grid solves

* `coarse_coarse_ksp::String`: coarse coarse grid solver in case we use redundant or telescope coarse grid solves

* `cores::Int64`: number of cores used in the simulation

* `cores_coarse::Int64`: number of cores used for coarse grid solver (in case we use pctelescope)

* `gamg_threshold::Float64`: GAMG threshold

* `gamg_coarse_eq_limit::Int64`: GAMG coarse grid equation limit

* `gamg_repartition::Bool`: GAMG repartition coarse grids? (default=false)

* `gamg_parallel_coarse::Bool`: GAMG parallel coarse grid solver? (default=false)

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Multigrid.jl#L5-L10)

[`LaMEM.LaMEM_Model.Output`](#LaMEM.LaMEM_Model.Output) — Type

    Structure that contains the LaMEM output options

* `out_file_name::Any`: output file name
    
* `out_dir::Any`: output directory
    
* `param_file_name::Any`: parameter filename
    
* `write_VTK_setup::Any`: write VTK initial model setup
    
* `out_pvd::Any`: activate writing .pvd file
    
* `out_phase::Any`: dominant phase
    
* `out_density::Any`: density
    
* `out_visc_total::Any`: total (viscoelastoplastic) viscosity
    
* `out_visc_creep::Any`: creep viscosity
    
* `out_velocity::Any`: velocity
    
* `out_pressure::Any`: (dynamic) pressure
    
* `out_tot_press::Any`: total pressure
    
* `out_eff_press::Any`: effective pressure
    
* `out_over_press::Any`: out_over_press
    
* `out_litho_press::Any`: lithospheric pressure
    
* `out_pore_press::Any`: pore pressure
    
* `out_temperature::Any`: temperature
    
* `out_dev_stress::Any`: deviatoric strain rate tensor
    
* `out_j2_dev_stress::Any`: second invariant of deviatoric stress tensor
    
* `out_strain_rate::Any`: deviatoric strain rate tensor
    
* `out_j2_strain_rate::Any`: second invariant of strain rate tensor
    
* `out_shmax::Any`: sh max
    
* `out_ehmax::Any`: eh max
    
* `out_yield::Any`: yield stress
    
* `out_rel_dif_rate::Any`: relative proportion of diffusion creep strainrate
    
* `out_rel_dis_rate::Any`: relative proportion of dislocation creep strainrate
    
* `out_rel_prl_rate::Any`: relative proportion of peierls creep strainrate
    
* `out_rel_pl_rate::Any`: relative proportion of plastic strainrate
    
* `out_plast_strain::Any`: accumulated plastic strain
    
* `out_plast_dissip::Any`: plastic dissipation
    
* `out_tot_displ::Any`: total displacement
    
* `out_moment_res::Any`: momentum residual
    
* `out_cont_res::Any`: continuity residual
    
* `out_energ_res::Any`: energy residual
    
* `out_melt_fraction::Any`: Melt fraction
    
* `out_fluid_density::Any`: fluid density
    
* `out_conductivity::Any`: conductivity
    
* `out_vel_gr_tensor::Any`: velocity gradient tensor
    
* `out_surf::Any`: activate surface output
    
* `out_surf_pvd::Any`: activate writing .pvd file
    
* `out_surf_velocity::Any`: surface velocity
    
* `out_surf_topography::Any`: surface topography
    
* `out_surf_amplitude::Any`: amplitude of topography (=topo-average(topo))
    
* `out_mark::Any`: activate marker output
    
* `out_mark_pvd::Any`: activate writing .pvd file
    
* `out_avd::Any`: activate AVD phase output
    
* `out_avd_pvd::Any`: activate writing .pvd file
    
* `out_avd_ref::Any`: AVD grid refinement factor
    
* `out_ptr::Any`: activate
    
* `out_ptr_ID::Any`: ID of the passive tracers
    
* `out_ptr_phase::Any`: phase of the passive tracers
    
* `out_ptr_Pressure::Any`: interpolated pressure
    
* `out_ptr_Temperature::Any`: temperature
    
* `out_ptr_MeltFraction::Any`: melt fraction computed using P-T of the marker
    
* `out_ptr_Active::Any`: option that highlight the marker that are currently active
    
* `out_ptr_Grid_Mf::Any`: option that allow to store the melt fraction seen within the cell
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Output.jl#L6-L10)

[`LaMEM.LaMEM_Model.PassiveTracers`](#LaMEM.LaMEM_Model.PassiveTracers) — Type

    Structure that contains the LaMEM passive tracers parameters.

* `Passive_Tracer::Int64`: activate passive tracers?"

* `PassiveTracer_Box::Union{Nothing, Vector{Float64}}`: Dimensions of box in which we distribute passive tracers \[Left, Right, Front, Back, Bottom, Top\]

* `PassiveTracer_Resolution::Vector{Int64}`: The number of passive tracers in every direction

* `PassiveTracer_ActiveType::Union{Nothing, String}`: Under which condition are they activated? \["Always"\], "Melt_Fraction", "Temperature", "Pressure", "Time"

* `PassiveTracer_ActiveValue::Union{Nothing, Float64}`: The value to activate them

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/PassiveTracers.jl#L8-L13)

[`LaMEM.LaMEM_Model.Phase`](#LaMEM.LaMEM_Model.Phase) — Type

    Defines the material properties for each of the phases

* `ID::Union{Nothing, Int64}`: Material phase ID
    
* `Name::Union{Nothing, String}`: Description of the phase
    
* `rho::Union{Nothing, Float64}`: Density \[kg/m^3\]
    
* `eta::Union{Nothing, Float64}`: Linear viscosity \[Pas\]
    
* `visID::Union{Nothing, Int64}`: material ID for phase visualization (default is ID)
    
* `diff_prof::Union{Nothing, String}`: Build-in DIFFUSION creep profiles:
    
    Example: `"Dry__Olivine_diff_creep-Hirth_Kohlstedt_2003"`
    
    Available build-in diffusion creep rheologies are:
    
    1.  From \[Hirth, G. and Kohlstedt D. (2003), Rheology of the upper mantle and the mantle wedge: A view from the experimentalists\]:
    
    * `"Dry_Olivine_diff_creep-Hirth_Kohlstedt_2003"`
    * `"Wet_Olivine_diff_creep-Hirth_Kohlstedt_2003_constant_C_OH"`
    * `"Wet_Olivine_diff_creep-Hirth_Kohlstedt_2003"`
    
    1.  From \[Rybacki and Dresen, 2000, JGR\]:
    
    * `"Dry_Plagioclase_RybackiDresen_2000"`
    * `"Wet_Plagioclase_RybackiDresen_2000"`
    
    Note that you can always specify your own, by setting `Bd`, `Ed`, `Vd` accordingly.
    

* `disl_prof::Union{Nothing, String}`: Build-in DISLOCATION creep profiles:
    
    Example: `"Granite-Tirel_et_al_2008"`
    
    Available build-in dislocation creep rheologies are:
    
    1.  From \[Ranalli 1995\]:
    
    * `"Dry_Olivine-Ranalli_1995"`
    * `"Wet_Olivine-Ranalli_1995"`
    * `"Wet_Quarzite-Ranalli_1995"`
    * `"Quarzite-Ranalli_1995"`
    * `"Mafic_Granulite-Ranalli_1995"`
    * `"Plagioclase_An75-Ranalli_1995"`
    
    1.  From \[Carter and Tsenn (1986). Flow properties of continental lithosphere - page 18\]:
    
    * `"Quartz_Diorite-Hansen_Carter_1982"`
    
    1.  From \[J. de Bremond d'Ars et al. Tectonophysics (1999). Hydrothermalism and Diapirism in the Archaean: gravitational instability constrains. - page 5\]
    
    * `"Diabase-Caristan_1982"`
    * `"Tumut_Pond_Serpentinite-Raleigh_Paterson_1965"`
    
    1.  From \[Mackwell, Zimmerman & Kohlstedt (1998). High-temperature deformation\]:
    
    * `"Maryland_strong_diabase-Mackwell_et_al_1998"`
    
    1.  From \[Ueda et al (PEPI 2008)\]:
    
    * `"Wet_Quarzite-Ueda_et_al_2008"`
    
    1.  From \[Huismans et al 2001\]:
    
    * `"Diabase-Huismans_et_al_2001"`
    * `"Granite-Huismans_et_al_2001"`
    
    1.  From \[Burg And Podladchikov (1999)\]:
    
    * `"Dry_Upper_Crust-Schmalholz_Kaus_Burg_2009"`
    * `"Weak_Lower_Crust-Schmalholz_Kaus_Burg_2009"`
    * `"Olivine-Burg_Podladchikov_1999"`
    
    1.  From \[Rybacki and Dresen, 2000, JGR\]:
    
    * `"Dry_Plagioclase_RybackiDresen_2000"`
    * `"Wet_Plagioclase_RybackiDresen_2000"`
    
    1.  From \[Hirth, G. & Kohlstedt (2003), D. Rheology of the upper mantle and the mantle wedge: A view from the experimentalists\]:
    
    * `"Wet_Olivine_disl_creep-Hirth_Kohlstedt_2003"`
    * `"Wet_Olivine_disl_creep-Hirth_Kohlstedt_2003_constant_C_OH"`
    * `"Dry_Olivine_disl_creep-Hirth_Kohlstedt_2003"`
    
    1.  From \[SchmalholzKausBurg(2009), Geology (wet olivine)\]:
    
    * `"Wet_Upper_Mantle-Burg_Schmalholz_2008"`
    * `"Granite-Tirel_et_al_2008"`
    
    1.  From \[Urai et al.(2008)\]:
    
    * `"Ara_rocksalt-Urai_et_al.(2008)"`
    
    1.  From \[Bräuer et al. (2011) Description of the Gorleben site (PART 4): Geotechnical exploration of the Gorleben salt dome - page 126\]:
    
    * `"RockSaltReference_BGRa_class3-Braeumer_et_al_2011"`
    
    1.  From \[Mueller and Briegel (1978)\]:
    
    * `"Polycrystalline_Anhydrite-Mueller_and_Briegel(1978)"`
    
    Note that you can always specify your own, by setting `Bn`, `En`, `Vn`, and `n` accordingly.
    

* `peir_prof::Union{Nothing, String}`: Build-in PEIERLS creep profiles:
    
    example: `"Olivine_Peierls-Kameyama_1999"`
    
    Available profiles:
    
    * `"Olivine_Peierls-Kameyama_1999"`

* `rho_n::Union{Nothing, Float64}`: depth-dependent density model parameter
    
* `rho_c::Union{Nothing, Float64}`: depth-dependent density model parameter
    
* `beta::Union{Nothing, Float64}`: pressure-dependent density model parameter
    
* `G::Union{Nothing, Float64}`: shear modulus
    
* `Kb::Union{Nothing, Float64}`: bulk modulus
    
* `E::Union{Nothing, Float64}`: Young's modulus
    
* `nu::Union{Nothing, Float64}`: Poisson's ratio
    
* `Kp::Union{Nothing, Float64}`: pressure dependence parameter
    
* `Bd::Union{Nothing, Float64}`: DIFFUSION creep pre-exponential constant
    
* `Ed::Union{Nothing, Float64}`: activation energy
    
* `Vd::Union{Nothing, Float64}`: activation volume
    
* `eta0::Union{Nothing, Float64}`: POWER LAW reference viscosity
    
* `e0::Union{Nothing, Float64}`: reference strain rate
    
* `Bn::Union{Nothing, Float64}`: DISLOCATION creep pre-exponential constant
    
* `En::Union{Nothing, Float64}`: activation energy
    
* `Vn::Union{Nothing, Float64}`: activation volume
    
* `n::Union{Nothing, Float64}`: power law exponent
    
* `Bp::Union{Nothing, Float64}`: PEIERLS creep pre-exponential constant
    
* `Ep::Union{Nothing, Float64}`: activation energy
    
* `Vp::Union{Nothing, Float64}`: activation volume
    
* `taup::Union{Nothing, Float64}`: scaling stress
    
* `gamma::Union{Nothing, Float64}`: approximation parameter
    
* `q::Union{Nothing, Float64}`: stress-dependence parameter
    
* `eta_fk::Union{Nothing, Float64}`: reference viscosity for Frank-Kamenetzky viscosity
    
* `gamma_fk::Union{Nothing, Float64}`: gamma parameter for Frank-Kamenetzky viscosity
    
* `TRef_fk::Union{Nothing, Float64}`: reference Temperature for Frank-Kamenetzky viscosity (if not set it is 0°C)
    
* `ch::Union{Nothing, Float64}`: cohesion
    
* `fr::Union{Nothing, Float64}`: friction angle
    
* `eta_st::Union{Nothing, Float64}`: stabilization viscosity (default is eta_min)
    
* `eta_vp::Union{Nothing, Float64}`: viscoplastic plasticity regularisation viscosity
    
* `rp::Union{Nothing, Float64}`: pore-pressure ratio
    
* `chSoftID::Union{Nothing, Int64}`: friction softening law ID
    
* `frSoftID::Union{Nothing, Int64}`: cohesion softening law ID
    
* `healID::Union{Nothing, Int64}`: healing ID, points to healTau in Softening
    
* `alpha::Union{Nothing, Float64}`: thermal expansivity
    
* `Cp::Union{Nothing, Float64}`: specific heat (capacity), J⋅K−1⋅kg−1
    
* `k::Union{Nothing, Float64}`: thermal conductivity
    
* `A::Union{Nothing, Float64}`: radiogenic heat production
    
* `T::Union{Nothing, Float64}`: optional temperature to set within the phase
    
* `Latent_hx::Union{Nothing, Float64}`: optional, used for dike heating, J/kg
    
* `T_liq::Union{Nothing, Float64}`: optional, used for dike heating, liquidus temperature of material, celsius
    
* `T_sol::Union{Nothing, Float64}`: optional, used for dike heating, solidus temperature of material, celsius
    
* `T_Nu::Union{Nothing, Float64}`: default value for thermal conductivity boundary
    
* `nu_k::Union{Nothing, Float64}`: optional parameter, Nusselt number for use with conductivity
    
* `rho_ph::Union{Nothing, String}`: name of the phase diagram you want to use (still needs rho to be defined for the initial guess of pressure)
    
* `rho_ph_dir::Union{Nothing, String}`: in case the phase diagram has a different path provide the path (without the name of the actual PD) here
    
* `mfc::Union{Nothing, Float64}`: melt fraction viscosity correction factor (positive scalar)
    
* `GeoParams::Union{Nothing, Vector{GeoParams.MaterialParameters.ConstitutiveRelationships.AbstractCreepLaw}}`: GeoParams creeplaws
    
    Set diffusion or dislocation creeplaws as provided by the GeoParams package:
    
        julia> using GeoParams
        julia> a = SetDiffusionCreep(GeoParams.Diffusion.dry_anorthite_Rybacki_2006);
        julia> p = Phase(ID=1,Name="test", GeoParams=[a]);
    
    Note that GeoParams should be a vector, as you could, for example, have diffusion and dislocation creep parameters
    
    Note also that this will overwrite any other creeplaws provided in the Phase struct.
    

* `grainsize::Union{Nothing, Float64}`: grainsize [m](not%20used%20in%20LaMEM) This is not actually used in LaMEM, but is required when setting diffusion creep parameters by using GeoParams

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Materials.jl#L7-L11)

[`LaMEM.LaMEM_Model.PhaseAggregate`](#LaMEM.LaMEM_Model.PhaseAggregate) — Type

    Defines phase aggregates, which can be useful for visualization purposes

* `name::String`: Name of the phase aggregate
    
* `phaseID::Union{Nothing, Vector{Int64}}`: Phases to be combined
    
* `numPhase::Union{Nothing, Int64}`: number of aggregated phases
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Materials.jl#L442-L446)

[`LaMEM.LaMEM_Model.PhaseTransition`](#LaMEM.LaMEM_Model.PhaseTransition) — Type

    Defines phase transitions on markers (that change the Phase ID of a marker depending on some conditions)

* `ID::Int64`: Phase_transition law ID
    
* `Type::String`: \[Constant, Clapeyron, Box\]: Constant - the phase transition occurs only at a fixed value of the parameter; Clapeyron - clapeyron slope
    
* `Name_Clapeyron::Union{Nothing, String}`: Type of predefined Clapeyron slope, such as Mantle_Transition_660km
    
* `PTBox_Bounds::Union{Nothing, Vector{Float64}}`: box bound coordinates: \[left, right, front, back, bottom, top\]
    
* `BoxVicinity::Union{Nothing, Int64}`: 1: only check particles in the vicinity of the box boundaries (2: in all directions)
    
* `Parameter_transition::Union{Nothing, String}`: \[T = Temperature, P = Pressure, Depth = z-coord, X=x-coord, Y=y-coord, APS = accumulated plastic strain, MeltFraction, t = time\] parameter that triggers the phase transition
    
* `ConstantValue::Union{Nothing, Float64}`: Value of the parameter \[unit of T,P,z, APS\]
    
* `number_phases::Union{Nothing, Int64}`: The number of involved phases \[default=1\]
    
* `PhaseAbove::Union{Nothing, Vector{Int64}}`: Above the chosen value the phase is 1, below it, the value is PhaseBelow
    
* `PhaseBelow::Union{Nothing, Vector{Int64}}`: Below the chosen value the phase is PhaseBelow, above it, the value is 1
    
* `PhaseInside::Union{Nothing, Vector{Int64}}`: Phase within the box \[use -1 if you don't want to change the phase inside the box\]
    
* `PhaseOutside::Union{Nothing, Vector{Int64}}`: Phase outside the box \[use -1 if you don't want to change the phase outside the box. If combined with OutsideToInside, all phases that come in are set to PhaseInside\]
    
* `PhaseDirection::Union{Nothing, String}`: \[BothWays=default; BelowToAbove; AboveToBelow\] Direction in which transition works
    
* `ResetParam::Union{Nothing, String}`: \[APS\] Parameter to reset on particles below PT or within box
    
* `PTBox_TempType::Union{Nothing, String}`: # Temperature condition witin the box \[none, constant, linear, halfspace\]
    
* `PTBox_topTemp::Union{Nothing, Float64}`: Temp @ top of box \[for linear & halfspace\]
    
* `PTBox_botTemp::Union{Nothing, Float64}`: Temp @ bottom of box \[for linear & halfspace\]
    
* `PTBox_thermalAge::Union{Nothing, Float64}`: Thermal age, usually in geo-units \[Myrs\] \[only in case of halfspace\]
    
* `PTBox_cstTemp::Union{Nothing, Float64}`: Temp within box \[only for constant T\]
    
* `v_box::Union{Nothing, Float64}`: \[optional\] only for NotInAirBox, velocity with which box moves in cm/yr
    
* `t0_box::Union{Nothing, Float64}`: \[optional\] beginning time of movemen in Myr
    
* `t1_box::Union{Nothing, Float64}`: \[optional\] end time of movement in Myr
    
* `clapeyron_slope::Union{Nothing, Float64}`: \[optional\] clapeyron slope of phase transition \[in K/MPa\]; `P = ( T - T0_clapeyron ) * clapeyron_slope + P0_clapeyron`
    
* `P0_clapeyron::Union{Nothing, Float64}`: \[optional\] P0_clapeyron \[Pa\]
    
* `T0_clapeyron::Union{Nothing, Float64}`: \[optional\] T0_clapeyron \[C\]
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Materials.jl#L490-L494)

[`LaMEM.LaMEM_Model.Scaling`](#LaMEM.LaMEM_Model.Scaling) — Type

    Scaling{T} is a structure that contains the scaling info, employed in the current simulation

* `Scaling::Any`: Scaling object (as in GeoParams), which can be `GEO_units()`, `NO_units()`, or `SI_units()`

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Scaling.jl#L3-L8)

[`LaMEM.LaMEM_Model.Softening`](#LaMEM.LaMEM_Model.Softening) — Type

    Defines strain softening parameters

* `ID::Int64`: softening law ID
    
* `APS1::Float64`: Begin of softening, in units of accumulated plastic strain (APS)
    
* `APS2::Float64`: End of softening, in units of accumulated plastic strain (APS)
    
* `A::Float64`: Reduction ratio
    
* `Lm::Union{Nothing, Float64}`: Material length scale (in selected units, e.g. km in geo)
    
* `APSheal2::Union{Nothing, Float64}`: APS when healTau2 activates
    
* `healTau::Union{Nothing, Float64}`: healing timescale parameter \[Myr\]
    
* `healTau2::Union{Nothing, Float64}`: healing timescale parameter \[Myr\] starting at APS=APSheal2
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Materials.jl#L380-L384)

[`LaMEM.LaMEM_Model.SolutionParams`](#LaMEM.LaMEM_Model.SolutionParams) — Type

    Structure that contains the LaMEM global solution parameters.

* `gravity::Vector{Float64}`: gravitational acceleration vector
    
* `FSSA::Float64`: free surface stabilization parameter \[0 - 1\]; The value has to be between 0 and 1
    

* `FSSA_allVel::Int64`: free surface stabilization parameter applied to all velocity components? Default is yes; if not it is only applied to the z-component

* `shear_heat_eff::Float64`: shear heating efficiency parameter \[0 - 1\]
    
* `Adiabatic_Heat::Float64`: Adiabatic Heating activation flag and efficiency. [0.0 - 1.0](e.g.,%200.5%20means%20that%20only%2050%20percent%20of%20the%20potential%20adiabatic%20heating%20affects%20the%20energy%20equation)
    
* `act_temp_diff::Int64`: temperature diffusion activation flag
    
* `act_therm_exp::Int64`: thermal expansion activation flag
    
* `act_steady_temp::Int64`: steady-state temperature initial guess activation flag
    
* `steady_temp_t::Float64`: time for (quasi-)steady-state temperature initial guess
    
* `nstep_steady::Int64`: number of steps for (quasi-)steady-state temperature initial guess (default = 1)
    
* `act_heat_rech::Int64`: recharge heat in anomalous bodies after (quasi-)steady-state temperature initial guess (=2: recharge after every diffusion step of initial guess)
    
* `init_lith_pres::Int64`: sets initial pressure to be the lithostatic pressure (stabilizes compressible setups in the first steps)
    
* `init_guess::Int64`: create an initial guess step (using constant viscosity `eta_ref` before starting the simulation
    
* `p_litho_visc::Int64`: use lithostatic instead of dynamic pressure for creep laws
    
* `p_litho_plast::Int64`: use lithostatic pressure for plasticity
    
* `p_lim_plast::Int64`: limit pressure at first iteration for plasticity
    
* `p_shift::Int64`: add a constant value \[MPa\] to the total pressure field, before evaluating plasticity (e.g., when the domain is located @ some depth within the crust)
    
* `act_p_shift::Int64`: pressure shift activation flag (enforce zero pressure on average in the top cell layer); note: this overwrites p_shift above!
    
* `eta_min::Float64`: viscosity lower bound \[Pas\]
    
* `eta_max::Float64`: viscosity upper limit \[Pas\]
    
* `eta_ref::Float64`: Reference viscosity (used for the initial guess) \[Pas\]
    
* `T_ref::Float64`: Reference temperature \[C\]
    
* `RUGC::Float64`: universal gas constant (you need to change this only for non-dimensional setups)
    
* `min_cohes::Float64`: cohesion lower bound \[Pa\]
    
* `min_fric::Float64`: friction lower bound \[degree\]
    
* `tau_ult::Float64`: ultimate yield stress \[Pa\]
    
* `rho_fluid::Float64`: fluid density for depth-dependent density model
    
* `gw_level_type::String`: ground water level type for pore pressure computation (see below)
    
* `gw_level::Float64`: ground water level at the free surface (if defined)
    
* `biot::Float64`: Biot pressure parameter
    
* `get_permea::Float64`: effective permeability computation activation flag
    
* `rescal::Float64`: stencil rescaling flag (for internal constraints, for example while computing permeability)
    
* `mfmax::Float64`: maximum melt fraction affecting viscosity reduction
    
* `lmaxit::Int64`: maximum number of local rheology iterations
    
* `lrtol::Float64`: local rheology iterations relative tolerance
    
* `act_dike::Int64`: dike activation flag (additonal term in divergence)
    
* `useTk::Int64`: switch to use T-dependent conductivity, 0: not active
    
* `dikeHeat::Int64`: switch to use Behn & Ito heat source in the dike
    
* `adiabatic_gradient::Float64`: Adiabatic gradient in combination with Behn & Ito dike
    
* `Compute_velocity_gradient::Int64`: compute the velocity gradient tensor 1: active, 0: not active. If active, it automatically activates the output in the .pvd file
    
* `Phasetrans::Int64`: Activate Phase Transitions on Particles or not, 0: not.
    
* `Passive_Tracer::Int64`: Activate Passive Tracers or not?
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/SolutionParams.jl#L8-L13)

[`LaMEM.LaMEM_Model.Solver`](#LaMEM.LaMEM_Model.Solver) — Type

    Structure that contains the LaMEM solver options

* `SolverType::String`: solver employed \[`"direct"` or `"multigrid"`\]
    
* `DirectSolver::String`: mumps/superlu_dist/pastix/umfpack (requires these external PETSc packages to be installed!)
    
* `DirectPenalty::Float64`: penalty parameter \[employed if we use a direct solver\]
    
* `MGLevels::Int64`: number of MG levels \[default=3\]
    
* `MGSweeps::Int64`: number of MG smoothening steps per level \[default=10\]
    
* `MGSmoother::String`: type of smoothener used \[chebyshev or jacobi\]
    
* `MGJacobiDamp::Float64`: Dampening parameter \[only employed for Jacobi smoothener; default=0.6\]
    
* `MGCoarseSolver::String`: coarse grid solver if using multigrid \[`"direct"` / `"mumps"` / `"superlu_dist"` or `"redundant"` \- more options specifiable through the command-line options `-crs_ksp_type` & `-crs_pc_type`\]
    
* `MGRedundantNum::Int64`: How many times do we copy the coarse grid? \[only employed for redundant solver; default is 4\]
    
* `MGRedundantSolver::String`: The coarse grid solver for each of the redundant solves \[only employed for redundant; options are `"mumps"`/`"superlu_dist"` with default `"superlu_dist"`\]
    
* `PETSc_options::Vector{String}`: List with (optional) PETSc options
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Solver.jl#L5-L9)

[`LaMEM.LaMEM_Model.Time`](#LaMEM.LaMEM_Model.Time) — Type

    Structure that contains the LaMEM timestepping information. An explanation of the paramneters is given in the struct `Time_info`

* `time_end::Float64`: simulation end time
    
* `dt::Float64`: initial time step
    
* `dt_min::Float64`: minimum time step (declare divergence if lower value is attempted)
    
* `dt_max::Float64`: maximum time step
    
* `dt_out::Float64`: output step (output at least at fixed time intervals)
    
* `inc_dt::Float64`: time step increment per time step (fraction of unit)
    
* `CFL::Float64`: CFL (Courant-Friedrichs-Lewy) criterion
    
* `CFLMAX::Float64`: CFL criterion for elasticity
    
* `nstep_max::Int64`: maximum allowed number of steps (lower bound: time_end/dt_max)
    
* `nstep_out::Int64`: save output every n steps; Set this to -1 to deactivate saving output
    
* `nstep_rdb::Int64`: save restart database every n steps
    
* `num_dt_periods::Int64`: number of time stepping periods
    
* `time_dt_periods::Vector{Int64}`: timestamps where timestep should be fixed (first entry has to 0)
    
* `step_dt_periods::Vector{Float64}`: target timesteps ar timestamps above
    
* `nstep_ini::Int64`: save output for n initial steps
    
* `time_tol::Float64`: relative tolerance for time comparisons
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Time.jl#L5-L9)

[`LaMEM.LaMEM_Model.VelCylinder`](#LaMEM.LaMEM_Model.VelCylinder) — Type

    LaMEM boundary condition internal velocty cylinder `VelCylinder` object

* `baseX::Float64`: X-coordinate of base of cylinder
    
* `baseY::Float64`: Y-coordinate of base of cylinder
    
* `baseZ::Float64`: Z-coordinate of base of cylinder
    
* `capX::Float64`: X-coordinate of cap of cylinder
    
* `capY::Float64`: Y-coordinate of cap of cylinder
    
* `capZ::Float64`: Z-coordinate of cap of cylinder
    
* `radius::Float64`: radius of cylinder
    
* `vx::Union{Nothing, Float64}`: Vx velocity of cylinder (default is unconstrained)
    
* `vy::Union{Nothing, Float64}`: Vy velocity of cylinder (default is unconstrained)
    
* `vz::Union{Nothing, Float64}`: Vz velocity of cylinder (default is unconstrained)
    
* `advect::Int64`: cylinder advection flag
    
* `vmag::Float64`: magnitude of velocity applied along the cylinder's axis of orientation
    
* `type::String`: velocity profile \[uniform or parabolic\]
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/BoundaryConditions.jl#L167-L171)

[`LaMEM.LaMEM_Model.VelocityBox`](#LaMEM.LaMEM_Model.VelocityBox) — Type

    Define velocity regions within the modelling region, by specifying its center point and width along the three axis.

* `cenX::Float64`: X-coordinate of center of box
    
* `cenY::Float64`: Y-coordinate of center of box
    
* `cenZ::Float64`: Z-coordinate of center of box
    
* `widthX::Float64`: Width of box in x-direction
    
* `widthY::Float64`: Width of box in y-direction
    
* `widthZ::Float64`: Width of box in Z-direction
    
* `vx::Union{Nothing, Float64}`: Vx velocity of box (default is unconstrained)
    
* `vy::Union{Nothing, Float64}`: Vx velocity of box (default is unconstrained)
    
* `vz::Union{Nothing, Float64}`: Vx velocity of box (default is unconstrained)
    
* `advect::Int64`: box advection flag
    

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/BoundaryConditions.jl#L9-L14)

[`GeophysicalModelGenerator.above_surface`](#GeophysicalModelGenerator.above_surface-Tuple{Model,%20GeophysicalModelGenerator.CartData}) — Method

    above_surface(model::Model, DataSurface_Cart::CartData)

Returns a boolean grid that is `true` if the `Phases/Temp` grid are above the surface

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/GMG_interface.jl#L110-L114)

[`GeophysicalModelGenerator.add_box!`](#GeophysicalModelGenerator.add_box!-Tuple{Model}) — Method

    add_box!(model::Model; xlim=Tuple{2}, [ylim=Tuple{2}], zlim=Tuple{2},
            Origin=nothing, StrikeAngle=0, DipAngle=0,
            phase = ConstantPhase(1),
            T=nothing )

Adds a box with phase & temperature structure to a 3D model setup. This simplifies creating model geometries in geodynamic models See the documentation of the GMG routine for the full options.

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/GMG_interface.jl#L9-L18)

[`GeophysicalModelGenerator.add_cylinder!`](#GeophysicalModelGenerator.add_cylinder!-Tuple{Model}) — Method

    add_cylinder!(model::Model;                                      # required input
                    base=Tuple{3}, cap=Tuple{3}, radius=Tuple{1},   # center and radius of the sphere
                    phase = ConstantPhase(1),                       # Sets the phase number(s) in the sphere
                    T=nothing )                                     # Sets the thermal structure (various fucntions are available)

See the documentation of the GMG routine

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/GMG_interface.jl#L41-L50)

[`GeophysicalModelGenerator.add_ellipsoid!`](#GeophysicalModelGenerator.add_ellipsoid!-Tuple{Model}) — Method

    add_ellipsoid!(model::Model;                                 # required input
                    cen=Tuple{3}, axes=Tuple{3},                # center and semi-axes of the ellpsoid
                    Origin=nothing, StrikeAngle=0, DipAngle=0,  # origin & dip/strike
                    phase = ConstantPhase(1),                   # Sets the phase number(s) in the box
                    T=nothing )

See the documentation of the GMG routine

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/GMG_interface.jl#L54-L63)

[`GeophysicalModelGenerator.add_layer!`](#GeophysicalModelGenerator.add_layer!-Tuple{Model}) — Method

    add_layer!(model::Model; xlim, ylim, zlim=Tuple{2},
            phase = ConstantPhase(1),
            T=nothing )

Adds a layer with phase & temperature structure to a 3D model setup. This simplifies creating model geometries in geodynamic models See the documentation of the GMG routine for the full options.

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/GMG_interface.jl#L21-L29)

[`GeophysicalModelGenerator.add_polygon!`](#GeophysicalModelGenerator.add_polygon!-Tuple{Model}) — Method

    add_polygon!(model::Model;                                 # required input
                    xlim::Vector, 
                    ylim=Vector,
                    zlim=Vector(), 
                    phase = ConstantPhase(1),                 # Sets the phase number(s) in the box
                    T=nothing)

See the documentation of the GMG routine

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/GMG_interface.jl#L67-L77)

[`GeophysicalModelGenerator.add_slab!`](#GeophysicalModelGenerator.add_slab!-Tuple{Model,%20GeophysicalModelGenerator.Trench}) — Method

    add_slab!(model::Model;                                 # required input
                    trench::Trench; 
                    phase = ConstantPhase(1),                 # Sets the phase number(s) in the box
                    T=nothing)

See the documentation of the GMG routine

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/GMG_interface.jl#L81-L89)

[`GeophysicalModelGenerator.add_sphere!`](#GeophysicalModelGenerator.add_sphere!-Tuple{Model}) — Method

    add_sphere!(model::Model; cen=Tuple{3}, radius=Tuple{1}, phase = ConstantPhase(1), T=nothing)

See the documentation of the GMG routine

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/GMG_interface.jl#L33-L38)

[`GeophysicalModelGenerator.add_stripes!`](#GeophysicalModelGenerator.add_stripes!-Tuple{Model}) — Method

    add_stripes!(Phase, Grid::AbstractGeneralGrid;
                stripAxes       = (1,1,0),
                stripeWidth     =  0.2,
                stripeSpacing   =  1,
                Origin          =  nothing,
                StrikeAngle     =  0,
                DipAngle        =  10,
                phase           =  ConstantPhase(3),
                stripePhase     =  ConstantPhase(4))

See the documentation of the GMG routine

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/GMG_interface.jl#L92-L105)

[`GeophysicalModelGenerator.below_surface`](#GeophysicalModelGenerator.below_surface-Tuple{Model,%20GeophysicalModelGenerator.CartData}) — Method

    below_surface(model::Model, DataSurface_Cart::CartData)

Returns a boolean grid that is `true` if the `Phases/Temp` grid are below the surface

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/GMG_interface.jl#L139-L143)

[`LaMEM.IO_functions.passivetracer_time`](#LaMEM.IO_functions.passivetracer_time) — Function

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Model.jl#L244)

[`LaMEM.IO_functions.passivetracer_time`](#LaMEM.IO_functions.passivetracer_time-Tuple{Union{Int64,%20Vector{Int64}},%20Model}) — Method

    PT = passivetracer_time(ID::Union{Vector{Int64},Int64}, model::Model)

This reads passive tracers with `ID` from a LaMEM simulation specified by `model`, and returns a named tuple with the temporal evolution of these passive tracers. We return `x`,`y`,`z` coordinates and all fields specified in `FileName` for particles number `ID`.

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Model.jl#L250-L256)

[`LaMEM.IO_functions.read_LaMEM_simulation`](#LaMEM.IO_functions.read_LaMEM_simulation-Tuple{Model}) — Method

    Timestep, FileNames, Time = read_LaMEM_simulation(model::Model; phase=false, surf=false, passive_tracers=false)

Reads a LaMEM simulation as specified in `model` and returns the timesteps, times and filenames of that simulation once it is finished.

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L256-L260)

[`LaMEM.IO_functions.read_LaMEM_timestep`](#LaMEM.IO_functions.read_LaMEM_timestep) — Function

    data, time = read_LaMEM_timestep(model::Model, TimeStep::Int64=0; fields=nothing, phase=false, surf=false, last=true)

Reads a specific `Timestep` from a simulation specified in `model`

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L263-L267)

[`LaMEM.LaMEM_Model.Check_LaMEM_Model`](#LaMEM.LaMEM_Model.Check_LaMEM_Model-Tuple{Model}) — Method

    Check_LaMEM_Model(m::Model)

Checks the LaMEM Setup Model `m` for errors

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ErrorChecking.jl#L5-L9)

[`LaMEM.LaMEM_Model.Create_Grid`](#LaMEM.LaMEM_Model.Create_Grid-NTuple{15,%20Any}) — Method

This creates a LaMEM grid

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Grid.jl#L147-L149)

[`LaMEM.LaMEM_Model.UpdateDefaultParameters`](#LaMEM.LaMEM_Model.UpdateDefaultParameters-Tuple{Model}) — Method

    model = UpdateDefaultParameters(model::Model)

This updates the default parameters depending on some of the input parameters. If you activate passive tracers, for example, it will also activate output for that

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/DefaultParams.jl#L3-L8)

[`LaMEM.LaMEM_Model.above_surface!`](#LaMEM.LaMEM_Model.above_surface!-Tuple{Model,%20GeophysicalModelGenerator.CartData}) — Method

    above_surface!(model::Model, DataSurface_Cart::CartData; phase::Int64=nothing, T::Number=nothing)

Sets the `Temp` or `Phases` above the surface `DataSurface_Cart` to a constant value.

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/GMG_interface.jl#L118-L122)

[`LaMEM.LaMEM_Model.add_dike!`](#LaMEM.LaMEM_Model.add_dike!-Tuple{Model,%20Dike}) — Method

    add_dike!(model::Model, dike::Dike)

This adds a phase transition `phase_trans` to `model`

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L182-L185)

[`LaMEM.LaMEM_Model.add_geom!`](#LaMEM.LaMEM_Model.add_geom!-Tuple{Model,%20Any}) — Method

    add_geom!(model::Model, geom_object)

This adds an internal geometric primitive object `geom_object` to the LaMEM Model Setup `model`.

Currently available primitive geom objects are:

* `GeomSphere`
* `GeomEllipsoid`
* `GeomBox`
* `GeomLayer`
* `GeomCylinder`
* `GeomRidgeSeg`
* `GeomHex`

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L192-L205)

[`LaMEM.LaMEM_Model.add_geom!`](#LaMEM.LaMEM_Model.add_geom!-Tuple{Model,%20Vararg{Any}}) — Method

    add_geom!(model::Model, geom_object)

Add several geometric objects @ once.

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L216-L219)

[`LaMEM.LaMEM_Model.add_petsc!`](#LaMEM.LaMEM_Model.add_petsc!-Tuple{Model,%20Vararg{Any}}) — Method

    add_petsc!(model::Model, option::String)

Adds one or more PETSc options to the model

**Example**

    julia> d = Model()
    julia> add_petsc!(d,"-snes_npicard 3")

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L127-L139)

[`LaMEM.LaMEM_Model.add_phase!`](#LaMEM.LaMEM_Model.add_phase!-Tuple{Model,%20Phase}) — Method

    add_phase!(model::Model, phase::Phase)

This adds a `phase` (with material properties) to `model`

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L54-L57)

[`LaMEM.LaMEM_Model.add_phase!`](#LaMEM.LaMEM_Model.add_phase!-Tuple{Model,%20Vararg{Any}}) — Method

    add_phase!(model::Model, phases...)

Add several phases @ once.

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L64-L67)

[`LaMEM.LaMEM_Model.add_phaseaggregate!`](#LaMEM.LaMEM_Model.add_phaseaggregate!-Tuple{Model,%20PhaseAggregate}) — Method

    add_phaseaggregate!(model::Model, phaseagg::PhaseAggregate)

This adds a phase aggregate law `phaseagg` to `model`

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L164-L167)

[`LaMEM.LaMEM_Model.add_phasetransition!`](#LaMEM.LaMEM_Model.add_phasetransition!-Tuple{Model,%20PhaseTransition}) — Method

    add_phasetransition!(model::Model, phase_trans::PhaseTransition)

This adds a phase transition `phase_trans` to `model`

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L173-L176)

[`LaMEM.LaMEM_Model.add_softening!`](#LaMEM.LaMEM_Model.add_softening!-Tuple{Model,%20Softening}) — Method

    add_softening!(model::Model, soft::Softening)

This adds a plastic softening law `soft` to `model`

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L155-L158)

[`LaMEM.LaMEM_Model.add_topography!`](#LaMEM.LaMEM_Model.add_topography!-Tuple{Model,%20GeophysicalModelGenerator.CartData}) — Method

    add_topography!(model::Model, topography::CartData; surf_air_phase=0, surf_topo_file="topography.txt", open_top_bound=1,  surf_level=0.0)

Adds the topography surface to the model

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L381-L385)

[`LaMEM.LaMEM_Model.add_vbox!`](#LaMEM.LaMEM_Model.add_vbox!-Tuple{Model,%20Vararg{Any}}) — Method

    add_vbox!(model::Model, vboxes...)

Add several phases @ once.

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L23-L26)

[`LaMEM.LaMEM_Model.add_vbox!`](#LaMEM.LaMEM_Model.add_vbox!-Tuple{Model,%20VelocityBox}) — Method

add_vbox!(model::Model, vbox::VelocityBox) This adds a `vbox` (with its properties) to `model`

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L14-L17)

[`LaMEM.LaMEM_Model.below_surface!`](#LaMEM.LaMEM_Model.below_surface!-Tuple{Model,%20GeophysicalModelGenerator.CartData}) — Method

    below_surface!(model::Model, DataSurface_Cart::CartData; phase::Union{Int64,Nothing}=nothing, T::Union{Number,Nothing}=nothing)

Sets the `Temp` or `Phases` below the surface `DataSurface_Cart` to a constant value.

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/GMG_interface.jl#L147-L151)

[`LaMEM.LaMEM_Model.compute_dof`](#LaMEM.LaMEM_Model.compute_dof-Tuple{Tuple{Int64,%20Int64,%20Int64}}) — Method

Returns the total degrees of freedom for a LaMEM simulation

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Multigrid.jl#L103-L105)

[`LaMEM.LaMEM_Model.copy_phase`](#LaMEM.LaMEM_Model.copy_phase-Tuple{Phase}) — Method

    copy_phase(phase::Phase; kwargs...)

This copies a phase with material properties, while allowing to change some parameters

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L365-L369)

[`LaMEM.LaMEM_Model.create_initialsetup`](#LaMEM.LaMEM_Model.create_initialsetup) — Function

    create_initialsetup(model::Model, cores::Int64=1, args::String=""; verbose=verbose)

Creates the initial model setup of LaMEM from `model`, which includes:

* Writing the LaMEM (*.dat) input file

and in case we do not employt geometric primitives to create the setup:

* Write the VTK file (if requested when `model.Output.write_VTK_setup=true`)
* Write the marker files to disk (if `model.ModelSetup.msetup="files"`)

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Model.jl#L261-L272)

[`LaMEM.LaMEM_Model.cross_section`](#LaMEM.LaMEM_Model.cross_section) — Function

    Cross = cross_section(cart::CartData, field::Symbol =:phase; x=nothing, y=nothing, z=nothing)

Creates a cross-section through the data and returns `x,z` coordinates

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L298-L302)

[`LaMEM.LaMEM_Model.cross_section`](#LaMEM.LaMEM_Model.cross_section) — Function

    data_tuple, axes_str = cross_section(model::LaMEM.Model, field=:phases; x=nothing, y=nothing, z=nothing)

This creates a cross-section through the initial model setup & returns a 2D array

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L285-L289)

[`LaMEM.LaMEM_Model.digitsep`](#LaMEM.LaMEM_Model.digitsep-Tuple{Integer}) — Method

    digitsep(value::Integer; separator=",", per_separator=3)

Convert an integer to a string, separating each `per_separator` digits by `separator`.

    digitsep(12345678)  # "12,345,678"
    digitsep(12345678, seperator= "'")  # "12'345'678"
    digitsep(12345678, seperator= "-", per_separator=4)  # "1234-5678"

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Multigrid.jl#L120-L129)

[`LaMEM.LaMEM_Model.flatten`](#LaMEM.LaMEM_Model.flatten-Tuple{GeophysicalModelGenerator.CartData,%20Symbol,%20Any,%20Any,%20Any}) — Method

Creates a 2D array out of a cross-section and a specified data field

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L315-L317)

[`LaMEM.LaMEM_Model.hasplasticity`](#LaMEM.LaMEM_Model.hasplasticity-Tuple{Phase}) — Method

    hasplasticity(p::Phase)

`true` if `p` contains plastic parameters (cohesion or friction angle)

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L427-L431)

[`LaMEM.LaMEM_Model.is_rectilinear`](#LaMEM.LaMEM_Model.is_rectilinear-Tuple{GeophysicalModelGenerator.CartData}) — Method

    is_rectilinear(topography::CartData)

Checks whether `topography` is rectilinear

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ErrorChecking.jl#L33-L37)

[`LaMEM.LaMEM_Model.isdefault`](#LaMEM.LaMEM_Model.isdefault-Tuple{Any,%20Any}) — Method

    isdefault(s1::S, s_default::S)

Checks whether a struct `s1` has default parameters `s_default`

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L410-L414)

[`LaMEM.LaMEM_Model.prepare_lamem`](#LaMEM.LaMEM_Model.prepare_lamem) — Function

    prepare_lamem(model::Model, cores::Int64=1, args:String=""; verbose=false)

Prepares a LaMEM run for the parameters that are specified in `model`, without running the simulation 1) Create the `*.dat` file 2) Write markers to disk in case we use a "files" setup

This is useful if you want to prepare a model on one machine but run it on another one (e.g. a cluster)

Set `model.Output.write_VTK_setup` to `true` if you want to write a `VTK` file of the model setup

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Model.jl#L213-L223)

[`LaMEM.LaMEM_Model.print_short`](#LaMEM.LaMEM_Model.print_short-Tuple{Multigrid}) — Method

This creates a single string, so we can use it in the command line

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Multigrid.jl#L263-L265)

[`LaMEM.LaMEM_Model.replace_phase!`](#LaMEM.LaMEM_Model.replace_phase!-Tuple{Model,%20Phase}) — Method

    replace_phase!(model::Model, phase_new::Phase; ID::Int64=nothing, Name::String=nothing)

This replaces a `phase` within a LaMEM Model Setup `model` with `phase_new` either based on its `Name` or `ID`. Note that it is expected that only one such phase is present in the current setup.

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L107-L112)

[`LaMEM.LaMEM_Model.rm_geom!`](#LaMEM.LaMEM_Model.rm_geom!-Tuple{Model}) — Method

    rm_geom!(model::Model)

This removes all existing geometric objects from `model`

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L226-L229)

[`LaMEM.LaMEM_Model.rm_last_phase!`](#LaMEM.LaMEM_Model.rm_last_phase!-Tuple{Model}) — Method

    rm_last_phase!(model::Model, phase::Phase)

This removes the last added `phase` from `model`

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L76-L79)

[`LaMEM.LaMEM_Model.rm_last_vbox!`](#LaMEM.LaMEM_Model.rm_last_vbox!-Tuple{Model}) — Method

rm_last_vbox!(model::Model) This removes the last added `vbox` from `model`

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L34-L37)

[`LaMEM.LaMEM_Model.rm_phase!`](#LaMEM.LaMEM_Model.rm_phase!-Tuple{Model,%20Int64}) — Method

    rm_phase!(model::Model, ID::Int64)

This removes a phase with `ID` from `model`

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L87-L90)

[`LaMEM.LaMEM_Model.rm_phase!`](#LaMEM.LaMEM_Model.rm_phase!-Tuple{Model}) — Method

    rm_phase!(model::Model)

This removes all existing phases from `model`

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L98-L101)

[`LaMEM.LaMEM_Model.rm_vbox!`](#LaMEM.LaMEM_Model.rm_vbox!-Tuple{Model}) — Method

    rm_vbox!(model::Model)

This removes all existing velocity boxes from `model`

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L45-L48)

[`LaMEM.LaMEM_Model.set_air`](#LaMEM.LaMEM_Model.set_air-Tuple{}) — Method

    set_air(; Name="air", ID=0, rho=1, alpha=nothing, eta=1e17, G=nothing, nu=nothing, fr=nothing, ch=nothing, k=30,Cp=1000)

Sets an air phase, with high conductivity

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L355-L358)

[`LaMEM.LaMEM_Model.set_geom!`](#LaMEM.LaMEM_Model.set_geom!-Tuple{Model,%20GeomSphere}) — Method

This sets the geometry

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L235-L238)

[`LaMEM.LaMEM_Model.stress_strainrate_0D`](#LaMEM.LaMEM_Model.stress_strainrate_0D-Tuple{Any,%20Vector}) — Method

    τ = stress_strainrate_0D(rheology, ε_vec::Vector; n=8, T=700, nstep_max=2, clean=true)

Computes the stress for a given strain rate and 0D rheology setup, for viscous creep rheologies. `n` is the resolution in `x,z`, `T` the temperature, `nstep_max` the number of time steps, `ε_vec` the strainrate vector (in 1/s).

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Utils.jl#L442-L449)

[`LaMEM.LaMEM_Model.within_bounds`](#LaMEM.LaMEM_Model.within_bounds-Tuple{Model,%20GeophysicalModelGenerator.CartData}) — Method

    within_bounds(model::Model, topography::CartData)

Verifies that the bounds of the topography grid are larger than that of the model

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ErrorChecking.jl#L45-L49)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile) — Function

    write_LaMEM_inputFile(d::Model,fname::String; dir=pwd())

Writes a LaMEM input file based on the data stored in Model

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Model.jl#L151-L155)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile-Tuple{Any,%20BoundaryConditions}) — Method

    write_LaMEM_inputFile(io, d::BoundaryConditions)

Writes the boundary conditions related parameters to file

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/BoundaryConditions.jl#L487-L490)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile-Tuple{Any,%20FreeSurface}) — Method

    write_LaMEM_inputFile(io, d::FreeSurface)

Writes the free surface related parameters to file

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/FreeSurface.jl#L116-L119)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile-Tuple{Any,%20GeomBox}) — Method

    write_LaMEM_inputFile(io, d::GeomBox)

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ModelSetup.jl#L240-L242)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile-Tuple{Any,%20GeomCylinder}) — Method

    write_LaMEM_inputFile(io, d::GeomCylinder)

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ModelSetup.jl#L459-L461)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile-Tuple{Any,%20GeomEllipsoid}) — Method

    write_LaMEM_inputFile(io, d::GeomEllipsoid)

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ModelSetup.jl#L187-L189)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile-Tuple{Any,%20GeomHex}) — Method

    write_LaMEM_inputFile(io, d::GeomHex)

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ModelSetup.jl#L343-L345)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile-Tuple{Any,%20GeomLayer}) — Method

    write_LaMEM_inputFile(io, d::GeomLayer)

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ModelSetup.jl#L408-L410)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile-Tuple{Any,%20GeomRidgeSeg}) — Method

    write_LaMEM_inputFile(io, d::GeomRidgeSeg)

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ModelSetup.jl#L302-L304)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile-Tuple{Any,%20GeomSphere}) — Method

    write_LaMEM_inputFile(io, d::GeomSphere)

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ModelSetup.jl#L138-L141)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile-Tuple{Any,%20Grid}) — Method

    write_LaMEM_inputFile(io, d::Grid)

This writes grid info to a LaMEM input file

**Example**

    julia> d=LaMEM.Grid(coord_x=[0.0, 0.7, 0.8, 1.0], bias_x=[0.3,1.0,3.0], nel_x=[10,4,2])
    julia> io = open("test.dat","w")
    julia> LaMEM.write_LaMEM_inputFile(io, d)
    julia> close(io)

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Grid.jl#L234-L248)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile-Tuple{Any,%20Materials}) — Method

    write_LaMEM_inputFile(io, d::Output)

Writes the free surface related parameters to file

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Materials.jl#L788-L791)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile-Tuple{Any,%20ModelSetup}) — Method

    write_LaMEM_inputFile(io, d::ModelSetup)

Writes options related to the Model Setup to disk

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/ModelSetup.jl#L482-L485)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile-Tuple{Any,%20Output}) — Method

    write_LaMEM_inputFile(io, d::Output)

Writes the free surface related parameters to file

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Output.jl#L203-L206)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile-Tuple{Any,%20PassiveTracers}) — Method

    write_LaMEM_inputFile(io, d::PassiveTracers)

Writes the boundary conditions related parameters to file

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/PassiveTracers.jl#L78-L81)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile-Tuple{Any,%20SolutionParams}) — Method

    write_LaMEM_inputFile(io, d::SolutionParams)

Writes the boundary conditions related parameters to file

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/SolutionParams.jl#L169-L172)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile-Tuple{Any,%20Solver}) — Method

    write_LaMEM_inputFile(io, d::Solver)

Writes the free surface related parameters to file

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Solver.jl#L73-L76)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile-Tuple{Any,%20Time}) — Method

Writes the Time related parameters to file

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Time.jl#L88-L90)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile-Tuple{Any,%20VelocityBox}) — Method

    write_LaMEM_inputFile(io, d::GeomSphere)

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/BoundaryConditions.jl#L77-L80)

[`LaMEM.LaMEM_Model.write_LaMEM_inputFile_PETSc`](#LaMEM.LaMEM_Model.write_LaMEM_inputFile_PETSc-Tuple{Any,%20Solver}) — Method

    write_LaMEM_inputFile_PETSc(io, d::Solver)

Writes the (optional) PETSc options to file

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Solver.jl#L109-L112)

[`LaMEM.Run.run_lamem`](#LaMEM.Run.run_lamem) — Function

    run_lamem(model::Model, cores::Int64=1, args:String=""; wait=true)

Performs a LaMEM run for the parameters that are specified in `model`

[source](https://github.com/JuliaGeodynamics/LaMEM.jl/blob/74478e1d388185f1d7175728c594033dbb087686/src/LaMEM_ModelGeneration/Model.jl#L190-L194)

[« Notebooks](../juliasetup_pluto/)[Run LaMEM »](../runlamem/)
