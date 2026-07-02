// Prose adapted closely from documentation/document.pdf. Keep this separate
// from generated parser metadata so the manual remains the editorial source.
const commandExplanations = {
  NODE: [
    "Nodes are the geometric reference points of a FEMaster model. Each node has an ID and three global coordinates. Elements, loads, supports, constraints, and features refer to nodes either directly or through node sets.",
    "The optional NSET keyword adds the nodes to a named set and defaults to NALL. A node may have up to six mechanical degrees of freedom: three translations and three rotations. Which of these are active depends on the connected elements and model objects. Solid elements usually need only translations, while shells, beams, point-mass springs, connectors, and constraints may require rotations.",
    "Node IDs must be non-negative integers; ID 0 is valid. Node coordinates are always global coordinates. Local coordinate systems do not change node positions: they only change how vector components are interpreted by commands such as *SUPPORT, *CLOAD, *DLOAD, *VLOAD, *SOLIDSECTION, *SHELLSECTION, and *CONNECTOR.",
    "The presence of a node alone does not imply six active unknowns. FEMaster builds the active degree-of-freedom set from the model objects that reference each node when the model is assembled. This keeps purely solid models smaller while allowing mixed solid-shell-beam models to use the same nodal component convention."
  ],
  ELEMENT: [
    "Elements connect nodes and define the interpolation used over a line, surface, or volume region. They are the fundamental discretization objects of the finite-element model: the selected element types determine how nodal degrees of freedom are converted into strains, curvatures, stresses, internal forces, mass contributions, geometric stiffness terms, and recoverable result fields.",
    "An element alone does not define a complete physical model. The element type defines the kinematic assumption and numerical integration domain; sections and materials define how that kinematic state becomes stiffness, mass, stress, or section force. Loads, supports, coordinate systems, and analysis controls complete the model.",
    "TYPE is required. ELSET is optional and defaults to EALL. Element IDs should be non-negative and unique because they are referenced by element sets, sections, loads, result fields, diagnostics, and postprocessing operations. Connectivity order matters: it defines local reference orientation, face numbering, integration-point layout, and, for structural elements, the local axes used for section behaviour.",
    "FEMaster provides solid, shell, beam, truss, and specialized panel element families. Solids represent three-dimensional continua, shells represent thin-walled structures, and line elements reduce slender members to section-based one-dimensional models. Element selection should follow the physical idealization; mesh quality, distortion, warping, and inappropriate element choices can reduce accuracy even when the solver converges."
  ],
  SURFACE: [
    "Surfaces describe loadable or connectable element sides. They are used for pressure loads, distributed surface loads, tie constraints, contact, and surface couplings. A surface is a selection of element sides; it is not a material or section property.",
    "A surface can be defined from one element side or generated for a selected side of every element in an element set. The side can be written as a number or as a token such as S1 or S2. For shells, SPOS and SNEG select the positive or negative side of the shell midsurface.",
    "Surface loads and interactions depend on the selected side and its normal direction. For shells, positive and negative sides share the same geometric midsurface but have opposite normals. Pressure sign, traction direction, contact gap, and interface behaviour therefore depend on element orientation and the side stored in the surface definition."
  ],
  NSET: [
    "A node set groups node IDs under a name. Sets are one of the main readability tools in a FEMaster deck because higher-level commands can target a modelling region instead of repeating raw mesh IDs. With GENERATE, the dataline is interpreted as a start ID, end ID, and increment; without GENERATE, the following datalines contain explicit node IDs.",
    "Node sets are used by supports, concentrated loads, couplings, connector definitions, point masses, diagnostic output, and result-extraction workflows. The set should describe the physical or modelling role of its nodes, while the active degrees of freedom still depend on the elements and objects connected to those nodes.",
    "A node may belong to several sets. Overlap is valid when the sets serve different purposes; it becomes relevant only when commands using those sets create contradictory constraints or duplicate loads. Because IDs can change after remeshing, set membership should be treated as mesh data and checked whenever the mesh changes."
  ],
  ELSET: [
    "An element set groups element IDs. With GENERATE, FEMaster creates a compact regular range; without it, datalines list explicit IDs. Element sets are central to the physical model because section definitions assign materials and geometric properties to named element regions.",
    "Solid, shell, beam, and truss sections, volume and inertial loads, topology fields, RBM constraints, and diagnostic workflows all use element sets. A set should normally represent a physical region or modelling role such as a solid core, shell skin, beam spar, topology domain, or gravity region.",
    "One element may belong to several sets when those sets are used for compatible purposes. The main risk is unintended multiple physical assignment, such as assigning incompatible sections to the same elements. Generated ranges are convenient for controlled meshes but can silently refer to the wrong region after remeshing, so production sets should preserve geometric intent rather than accidental numbering."
  ],
  SFSET: [
    "A surface set groups surface IDs and names a boundary or interface region above the raw element topology. With GENERATE, FEMaster creates a regular surface-ID range; otherwise the datalines list explicit surface IDs.",
    "Surface sets are used by pressure and distributed loads, ties, contact, and surface couplings. They allow these commands to target meaningful boundary names rather than individual faces. The orientation carried by the underlying surfaces remains part of the definition and controls pressure sign and interface direction.",
    "Surface sets can overlap, but operational overlap should be reviewed because the same region can otherwise receive duplicate loads or appear unintentionally on both sides of an interaction. Like other ID-based sets, surface sets must be regenerated or checked after the mesh changes."
  ],
  MATERIAL: [
    "Material definitions provide density, thermal expansion, and elastic constitutive behaviour. Sections reference materials by name; elements receive stiffness and mass only after a compatible section has assigned that material to an element set. A material can be defined once and reused by several shell, solid, beam, or truss sections.",
    "*MATERIAL activates or creates a named material entry. The name is supplied through MATERIAL or its NAME alias. The following *ELASTIC, *DENSITY, and *THERMALEXPANSION commands apply to the currently active material, so their position inside the material context determines which material receives the data."
  ],
  ELASTIC: [
    "*ELASTIC assigns the constitutive stiffness of the active material. The TYPE keyword selects the expected material model and therefore the shape of the datalines. FEMaster supports isotropic, generalised isotropic, orthotropic, and directly supplied ABD stiffness data.",
    "ISOTROPIC uses Young's modulus E and Poisson's ratio nu; the shear modulus follows from the isotropic relation. GENISO keeps the normal isotropic terms from E and nu but accepts an independent shear modulus G. This is useful when shear stiffness should be tuned independently, but it is only physically isotropic when G satisfies the isotropic relation.",
    "ORTHOTROPIC defines three material axes with separate Young's moduli, shear moduli, and Poisson ratios. It is intended for direction-dependent materials such as laminates, wood-like materials, and unidirectional composites. The material orientation supplied by the section is therefore part of the physical definition.",
    "TYPE=ABD reads a homogenized shell or laminate stiffness representation. The first 36 values form a 6 by 6 ABD matrix in row-major order and the last four values form a 2 by 2 transverse-shear matrix. This form is appropriate when laminate stiffness has already been homogenized outside FEMaster."
  ],
  DENSITY: [
    "*DENSITY sets the mass density of the active material. Density is used when FEMaster builds consistent or lumped mass contributions, inertia loads, inertia-relief mass properties, eigenfrequency systems, and transient equations.",
    "For continuum elements, density is integrated over the element volume. For shells and beams, section thickness or section area converts density into mass per area or mass per length. Density must use the same length, force, mass, and time convention as the rest of the deck."
  ],
  THERMALEXPANSION: [
    "*THERMALEXPANSION sets the scalar thermal expansion coefficient of the active material. Thermal loads combine this coefficient with the difference between the nodal temperature field and the reference temperature to form thermal strain.",
    "The equivalent mechanical load is assembled through the material constitutive law. A temperature field can exist without thermal expansion data, but it does not produce a meaningful elastic thermal-strain contribution for a material whose expansion coefficient is not defined."
  ],
  SOLIDSECTION: [
    "Sections connect element sets to physical data that the topology alone does not contain. A solid section assigns material behaviour to a set of three-dimensional continuum elements. It does not define thickness or area because these quantities are already represented by the volume mesh.",
    "MATERIAL supplies elastic stiffness, density, thermal expansion, and other constitutive data used by the selected analysis. The assignment affects stiffness and mass assembly, stress recovery, and every loadcase that depends on the material response.",
    "ORIENTATION is optional. It usually has no practical effect for an isotropic material, but it places the material axes for orthotropic or otherwise direction-dependent behaviour. The referenced element set should contain compatible solid elements and cover the intended physical region."
  ],
  SHELLSECTION: [
    "A shell section assigns material and thickness to a set of shell elements. Shells replace a thin three-dimensional body with a two-dimensional midsurface, so the section carries physical information that would otherwise be represented by the volume mesh.",
    "Thickness affects membrane stiffness, bending stiffness, transverse-shear behaviour, mass, stress recovery, and shell resultants. MATERIAL provides the constitutive response, while the optional ORIENTATION defines local material or resultant axes for direction-dependent shell behaviour.",
    "The referenced element set should contain shell elements with consistent normals and intended physical thickness. Section data and shell orientation together determine how local membrane, bending, shear, and resultant quantities are interpreted."
  ],
  PROFILE: [
    "A beam profile stores reusable cross-section constants for beam elements. It describes the geometry that a line element does not contain, including area, bending inertias, torsional constant, and any additional constants expected by the beam formulation.",
    "*BEAMSECTION references the profile by name and combines it with a material and local section direction. The resulting section defines axial, bending, and torsional stiffness as well as section-force recovery for the beam elements.",
    "Profile constants use the model's unit system. Area has units of length squared, bending inertias and torsional constants use higher powers of length, so an inconsistent unit conversion can change beam stiffness by orders of magnitude."
  ],
  BEAMSECTION: [
    "A beam section assigns a material and profile to a set of beam elements. Element connectivity defines the beam line; the section defines how that line behaves as a physical member. Together, material and profile determine axial, bending, and torsional stiffness, mass contribution where density is available, and section-force recovery.",
    "The optional dataline defines a local section direction. It is used to construct the beam coordinate system and determines which profile inertia acts about which physical direction. The direction must not be parallel to the beam axis, because a nearly parallel definition makes local-axis construction ill-conditioned.",
    "Beam orientation is part of the model rather than a display setting. Non-axisymmetric profiles have different stiffness about their local axes. If one element set contains members with incompatible intended orientations, the members should be separated or supplied with section directions that preserve the intended engineering axes."
  ],
  TRUSSSECTION: [
    "A truss section assigns material and cross-sectional area to truss elements. Truss elements carry axial force only, so the material provides axial constitutive response and AREA scales axial stiffness, mass contribution where density exists, and stress recovery.",
    "The section does not define bending inertia, torsional stiffness, offsets, or rotational behaviour. If the physical member transfers moment, resists bending, or has important eccentricity, a beam, shell, or solid idealization is required instead.",
    "AREA must be consistent with the model units and the intended idealization. It may be the actual area of a rod or an effective value chosen to reproduce axial stiffness, but the interpretation directly affects stress because truss stress is based on axial force divided by area."
  ],
  FIELD: [
    "*FIELD defines generic tabular data at nodes, elements, or element integration points. Fields can be defined at root level or inside a loadcase block depending on their use, and other commands reference them by name.",
    "TYPE selects the domain. NODE fields provide one row per node ID and are used for temperatures, initial velocities, or other nodal data. ELEMENT fields provide one row per element ID and are used for topology density and orientation. ELEMENT_IP stores data associated with element integration points.",
    "COLS specifies the component count and must be greater than zero; FEMaster limits fields to 64 components. FILL controls initialization of unspecified values and defaults to ZERO, while NAN initializes them with NaN. Datalines may omit individual values without overwriting existing entries."
  ],
  ORIENTATION: [
    "*ORIENTATION defines a named local coordinate system. Sections, loads, supports, and connectors reference these names when values are easier or physically more meaningful in a local basis. FEMaster transforms those components into the global basis used by the assembled equations.",
    "RECTANGULAR defines a Cartesian local basis from one, two, or three direction vectors. FEMaster normalizes and completes the basis. Two or three non-collinear directions provide the clearest frame; ambiguous or nearly parallel vectors do not define stable local axes.",
    "CYLINDRICAL defines a system from an origin, axis direction, and reference direction. The local radial and tangential directions depend on the evaluation point. Near the cylindrical axis the radial direction is not uniquely defined, so the reference direction establishes the angular frame.",
    "The same orientation can be used to interpret support components, nodal forces and moments, surface tractions, body forces, solid and shell material axes, and connector directions. In every case the local directions are transformed into global equation components before assembly."
  ],
  POINTMASS: [
    "Features add concentrated properties or auxiliary behaviour without defining element connectivity. *POINTMASS assigns concentrated mass, rotational inertia, and optional spring stiffness to the nodes of a node set. It is intended for equipment masses, simplified components, sensors, local inertias, or concentrated replacement springs.",
    "The dataline contains scalar mass, three rotational inertias, three translational spring stiffnesses, and three rotational spring stiffnesses. Missing values are interpreted as zero. The feature can therefore represent pure lumped mass, mass with rotary inertia, springs, or a combination of these properties.",
    "Point masses contribute to mass and inertia in eigenfrequency, transient, inertial-load, and inertia-relief workflows. In a static loadcase, mass alone does not create weight; gravity or acceleration must be introduced through a load definition."
  ],
  SUPPORT: [
    "Loads and supports are first assigned to named collectors. A support collector is a reusable list of support entries, and a loadcase later selects the collectors that participate in that analysis with *SUPPORTS. This allows the same model to combine fixture, symmetry, or temporary stabilization conditions differently in separate loadcases.",
    "*SUPPORT adds an entry to the collector named by SUPPORT_COLLECTOR. The collector is created when its name is first used. Each dataline targets a node set or a single node ID and prescribes values for the six nodal components: three translations followed by three rotations. Values are usually zero, but nonzero enforced displacement or rotation is also supported.",
    "Without ORIENTATION, the six values are interpreted in global coordinates. With an orientation, translational and rotational components are interpreted in that local system and transformed to global components before the constraint equations are built. Multiple support entries in the active collectors are assembled into the loadcase's global constraint set."
  ],
  CLOAD: [
    "*CLOAD adds concentrated nodal forces and moments to the collector named by LOAD_COLLECTOR. The target can be a node set or a single node ID. The six numeric components are generalized forces corresponding to the six nodal components: the first three are forces and the final three are moments.",
    "When ORIENTATION is present, force and moment components are entered in that local system and transformed to global components during assembly. AMPLITUDE references a scalar time function that multiplies the load entry in transient analysis. The load contributes only when its collector is activated by *LOADS in a loadcase."
  ],
  DLOAD: [
    "*DLOAD defines a vector traction on surfaces and stores it in a named load collector. The dataline targets a surface ID or surface set and supplies three traction components. During assembly, FEMaster integrates the traction with the surface shape functions to obtain equivalent nodal forces.",
    "Use DLOAD when the traction direction is known as a vector. With ORIENTATION, its components are interpreted in a local basis and transformed into global coordinates. Without an orientation they are global components. AMPLITUDE can scale the spatial traction pattern during transient analysis.",
    "The surface definition determines where the load is integrated, while the vector determines its direction. Unlike pressure, DLOAD does not derive its direction from the surface normal. It is assembled only when the selected LOAD_COLLECTOR is listed by *LOADS in the active loadcase."
  ],
  PLOAD: [
    "*PLOAD defines scalar pressure on a surface ID or surface set and stores the entry in a named load collector. FEMaster converts the scalar pressure into a traction using the normal of the selected surface, then integrates that traction into equivalent nodal forces.",
    "The selected element side and its normal direction determine the pressure sign convention. Opposite shell sides share the same midsurface but have opposite normals. PLOAD is therefore appropriate when the load should follow the surface normal; use *DLOAD when the load direction should be supplied as an explicit vector.",
    "AMPLITUDE can scale the pressure magnitude during transient analysis. As with every collector-based load, the pressure has no effect until its LOAD_COLLECTOR is activated by *LOADS inside a loadcase."
  ],
  VLOAD: [
    "*VLOAD defines a body-force density over an element region. The target is an element ID or element set, and the three dataline components form the volume-load vector that FEMaster integrates over the element domain to obtain equivalent nodal forces.",
    "The supplied values are force density per unit volume in the model's unit convention. ORIENTATION allows the vector to be entered in a local system, and AMPLITUDE allows time-dependent scaling. If gravitational acceleration is intended, density and unit conventions must be reflected consistently in the supplied load formulation."
  ],
  TLOAD: [
    "*TLOAD adds thermal strain loading to a load collector. It references a nodal temperature field and a reference temperature. The temperature difference is combined with the thermal expansion coefficient of each material to form thermal strain.",
    "FEMaster assembles the corresponding equivalent mechanical load through the element constitutive law. A meaningful thermal load therefore requires a compatible temperature field, materials with *THERMALEXPANSION data, and activation of the selected load collector through *LOADS."
  ],
  INERTIALOAD: [
    "*INERTIALOAD creates equivalent nodal loads from a prescribed rigid-body acceleration field over an element set. Its dataline supplies a reference point, translational acceleration at that point, angular velocity, and angular acceleration. These quantities define the acceleration at every material point in the selected region.",
    "FEMaster combines that acceleration field with the model mass and integrates the equivalent inertial load over the selected elements. CONSIDER_POINT_MASSES controls whether *POINTMASS features add concentrated inertia contributions. The resulting entry is stored in a load collector and must be activated by *LOADS."
  ],
  AMPLITUDE: [
    "An amplitude is a reusable scalar time function. Loads that reference an amplitude retain their spatial distribution while their magnitude is multiplied by the amplitude value during transient analysis.",
    "Each dataline supplies a time-value pair. LINEAR interpolates between adjacent points, STEP keeps the previous value until the next point, and NEAREST uses the value of the closest tabulated time. This separates the temporal history from the load collector that defines where and in which direction the load acts."
  ],
  RBM: [
    "Constraints define kinematic equations between degrees of freedom. *RBM generates up to six homogeneous equations for the nodes belonging to a selected element set. Its purpose is to suppress free rigid-body modes without fixing an arbitrary node, which is useful for stabilization, diagnostic models, and inertia-relief workflows.",
    "FEMaster collects the nodes used by the selected structural elements and determines the center of the region. Three equations suppress average translation, while three further equations suppress average rigid rotation from the first moment of the displacement field about that center.",
    "Only active translational degrees of freedom participate. ELSET also accepts SET as an alias and defaults to EALL. The generated equations enter the same global constraint set as supports, couplings, ties, and connectors."
  ],
  COUPLING: [
    "*COUPLING connects one master node set to either a slave node set or a slave surface set. Exactly one of SLAVE and SFSET is supplied. The six dataline flags select which translational and rotational components participate, and values greater than zero activate the corresponding component.",
    "KINEMATIC coupling generates compatibility equations. Slave translation follows master translation plus the small-rotation contribution produced by the offset between master and slave, while enabled rotations are tied directly. Equations are generated only for components whose degrees of freedom exist.",
    "STRUCTURAL coupling does not impose kinematic compatibility. It distributes a generalized master force and moment to slave translational forces in a least-norm, work-equivalent form so that the slave forces reproduce the master resultant. The master set normally represents one reference node because larger master sets can change the intended kinematic interpretation."
  ],
  TIE: [
    "*TIE binds slave nodes to a master surface set or master line set. The slave region can be a node set or a surface set. For each slave node, FEMaster searches for the closest master geometry within DISTANCE and projects the slave position onto that geometry.",
    "At the projected point, the master shape functions interpolate master motion. For each compatible active degree of freedom, FEMaster generates an equation that ties the slave value to this interpolated master value. Surface and line masters use their respective interpolation functions.",
    "If no master candidate is found within the search distance, no equation is generated for that slave node. ADJUST=YES moves the slave coordinate to the projected master position before equations are assembled; ADJUST=NO leaves the original geometry unchanged."
  ],
  CONTACT: [
    "*CONTACT defines frictionless node-to-surface penalty contact. The master side must be a surface set. The slave side can be a node set or surface set; for a slave surface set, FEMaster collects the unique nodes belonging to those surfaces.",
    "For each slave node, FEMaster searches nearby master surfaces, selects the closest projected point, and evaluates the normal gap using the master normal and CLEARANCE. FLIP reverses the master normal. Contact becomes active when the gap indicates penetration, and PENALTY scales the normal restoring force.",
    "The opposite contact force is distributed to master nodes with the surface shape functions. Only translational degrees of freedom are used. CONTACT is not a linear constraint equation: it contributes nonlinear internal force and tangent stiffness only in NONLINEARSTATIC, and linear stiffness assembly rejects models containing contact definitions."
  ],
  CONNECTOR: [
    "*CONNECTOR creates idealized relative-motion constraints between two node sets in a named coordinate system. Supported connector types map to masks that determine which of the six local relative translations and rotations are constrained.",
    "For a constrained translation, FEMaster transforms the selected local axis into global coordinates and constrains the relative displacement of the two node sets along that direction. Rotational components use the same directional transformation for relative rotations.",
    "The connector therefore acts along local axes rather than simply tying global components. COORDINATESYSTEM is part of the physical definition and controls which directions remain free or constrained for BEAM, HINGE, CYLINDRICAL, TRANSLATOR, JOIN, and JOINRX connector types."
  ],
  LOADCASE: [
    "A loadcase defines which analysis is executed on the previously defined model. The model describes what exists: nodes, elements, surfaces, sets, materials, sections, features, loads, supports, constraints, and fields. The loadcase describes which equation is assembled, which collectors are active, how constraints are handled, which solver is used, and which analysis-specific parameters are applied.",
    "TYPE is required and selects LINEARSTATIC, LINEARSTATICTOPO, NONLINEARSTATIC, EIGENFREQ, LINEARBUCKLING, or LINEARTRANSIENT. NAME is optional and identifies the loadcase in the deck and result output. Loadcases execute in input order and result blocks are written in the same order.",
    "Linear static analysis solves small-displacement equilibrium. Topology static analysis modifies element stiffness with density and orientation fields. Eigenfrequency analysis computes natural frequencies and mode shapes. Linear buckling performs a preload solve and then a geometric-stiffness eigenproblem. Nonlinear static analysis applies load incrementally with tangent iterations and is the only loadcase that assembles contact. Linear transient analysis solves linear dynamic response over time with implicit Newmark integration.",
    "Commands inside the loadcase activate collectors and configure the selected equation. The common operators are stiffness, mass, damping, and geometric stiffness. Supports and constraints form an equation set that is either reduced with a Nullspace basis or, where supported, introduced with Lagrange multipliers."
  ],
  SUPPORTS: [
    "*SUPPORTS activates one or more named support collectors that were populated by *SUPPORT. Entries in those collectors contribute prescribed-motion equations to the active loadcase's global constraint set.",
    "Support collectors are normally selected explicitly. If no collector is listed and the special all-support collector exists, FEMaster may use that collector as the active support set. Different loadcases can select different combinations without redefining the underlying support entries."
  ],
  LOADS: [
    "*LOADS activates named load collectors populated by *CLOAD, *DLOAD, *PLOAD, *VLOAD, *TLOAD, or *INERTIALOAD. Only collectors listed here contribute to the loadcase right-hand side; defining a load collector at root level does not assemble it automatically.",
    "Several collectors can be combined in one loadcase, while another loadcase can reuse a different combination. In transient analysis, amplitudes attached to entries in the active collectors are evaluated during time integration so the assembled load vector varies with time."
  ],
  SOLVER: [
    "*SOLVER selects the execution device and the strategy used for linear systems in the active loadcase. DEVICE=CPU is the standard path; DEVICE=GPU uses supported GPU solver routines and data structures when available.",
    "METHOD=DIRECT solves by matrix factorization followed by triangular solves or an equivalent direct procedure. Direct solvers are robust and predictable for small and medium models, with memory growth as their main limitation for large three-dimensional systems.",
    "METHOD=INDIRECT uses an iterative path that improves an approximate solution until its residual is sufficiently small. It can reduce factorization memory but is more sensitive to conditioning, scaling, matrix definiteness, and constraint handling. In FEMaster it is intended for NULLSPACE-reduced systems, not the saddle-point system created by Lagrange multipliers."
  ],
  CONSTRAINTMETHOD: [
    "Active supports and constraints define an algebraic equation set. *CONSTRAINTMETHOD selects how those equations are introduced into LINEARSTATIC and LINEARSTATICTOPO analyses. NONLINEARSTATIC uses the Nullspace path, while EIGENFREQ, LINEARBUCKLING, and LINEARTRANSIENT apply Nullspace constraints internally and do not accept this command.",
    "NULLSPACE parameterizes all admissible displacements as a particular constrained displacement plus a basis of independent coordinates. The stiffness, mass, and damping operators are projected into that reduced space, constrained degrees of freedom are removed from the solve, and the full displacement is recovered afterward.",
    "LAGRANGE retains the original displacement vector and adds one multiplier for every constraint equation. This preserves the full displacement variables and provides generalized constraint reactions, but creates an indefinite saddle-point system.",
    "In FEMaster, LAGRANGE requires a supported direct solver combination for linear static analysis. NULLSPACE is compatible with direct and iterative solution paths and is also the representation used by eigenvalue, transient, and current nonlinear workflows."
  ],
  NONLINEAR: [
    "*NONLINEAR configures the incremental Newton solve used by NONLINEARSTATIC and is valid only inside that loadcase type. It controls the initial increment count, maximum iterations, convergence tolerance, optional regularization of zero tangent rows, and the regularization scale.",
    "INCREMENTS sets the initial load increment as the reciprocal of the supplied count. FEMaster advances a scalar load factor from zero to one. In every increment it assembles current internal force and tangent stiffness, projects the residual into Nullspace coordinates, solves the reduced Newton equation, and applies the correction.",
    "An increment converges when the normalized reduced residual is below TOL. If convergence is not reached within MAXITER, FEMaster restores the last accepted state and cuts the load increment back. Fast convergence permits later increments to grow.",
    "REGULARIZE_ZERO_ROWS can add small diagonal stiffness to tangent rows whose norm is effectively zero, scaled by REGULARIZATION_ALPHA. The current nonlinear path requires a DIRECT solver and NULLSPACE constraints; it has no line search, arc-length control, contact damping, or automatic penalty scaling."
  ],
  INERTIARELIEF: [
    "Inertia relief is used for static analysis of free or nearly free structures. A free body under an unbalanced external load cannot have static equilibrium in an inertial frame because it accelerates as a rigid body. FEMaster computes an equivalent rigid-body acceleration and adds inertia loads that balance the external resultant force and moment.",
    "The remaining static solution describes elastic deformation relative to that accelerating rigid-body motion. This is appropriate when the body is intentionally free, such as a component under manoeuvre load, rather than when physical bearings, fixtures, bolts, or symmetry supports are merely missing from the model.",
    "CONSIDER_POINT_MASSES controls whether *POINTMASS features participate in the mass properties and inertial balance. It should include lumped masses that represent physical equipment and exclude features that are not intended to affect rigid-body inertia."
  ],
  REBALANCELOADS: [
    "*REBALANCELOADS modifies the active external loads so that their resultant force and moment vanish. It is a numerical balancing operation intended for load patterns that should already be self-equilibrated but retain small residuals from discretization, surface selection, rounding, or simplified load introduction.",
    "This differs from inertia relief. Inertia relief represents physical rigid-body acceleration of a free structure through inertia loads. Rebalancing assumes that the applied load itself should have zero resultant and corrects its residual imbalance before the static solve."
  ],
  REQUESTSTIFFNESS: [
    "*REQUESTSTIFFNESS requests stiffness-related matrix output for diagnostics or external comparison; it does not change the physical model. In LINEARSTATIC and LINEARSTATICTOPO it writes active and reduced stiffness data, and the linear static path also writes the reduced right-hand side.",
    "In LINEARBUCKLING it writes preload and reduced stiffness matrices. In NONLINEARSTATIC it writes the final tangent stiffness and final reduced tangent. The output is intended for assembly validation, reduced-matrix comparison, and investigation of singular or unexpectedly constrained models."
  ],
  REQUESTSTGEOM: [
    "*REQUESTSTGEOM requests geometric stiffness matrix output in LINEARBUCKLING. Geometric stiffness is assembled from the preload stress state and participates in the generalized eigenproblem that produces buckling factors and modes.",
    "The command is a diagnostic output control and does not change the model or the buckling equation. It is useful when validating preload stress effects, comparing buckling operators externally, or investigating unexpected buckling factors."
  ],
  CONSTRAINTSUMMARY: [
    "*CONSTRAINTSUMMARY requests a diagnostic listing of the internally generated constraint equations for the active loadcase. Output is grouped by origin, including supports, RBM equations, couplings, ties, and connectors.",
    "It reports the equations that FEMaster generated without changing them. The listing can be used to inspect which supports and interactions entered a loadcase, how many tie projections were found, and how the assembled model acquired or lost independent degrees of freedom."
  ],
  NUMEIGENVALUES: [
    "*NUMEIGENVALUES sets the number of eigenpairs requested by EIGENFREQ or LINEARBUCKLING. The count must be positive and is limited by the size of the reduced system available to the eigensolver.",
    "In eigenfrequency analysis, each pair contains a natural frequency and mode shape. In linear buckling, each pair contains an ideal buckling factor and corresponding mode around the preloaded reference state."
  ],
  SIGMA: [
    "*SIGMA sets the eigenvalue shift used by LINEARBUCKLING. After the static preload, FEMaster reduces elastic and geometric stiffness and solves a generalized eigenproblem for buckling factors. The shift guides the eigensolver toward the part of that spectrum that should be returned.",
    "A poorly selected shift can direct the solver to an unintended spectral region or miss factors of engineering interest. Reported buckling factors are sorted after the eigensolve.",
    "When SIGMA is zero, FEMaster estimates a shift from the pre-buckling displacement with a Rayleigh quotient. If that estimate is not finite and positive, it uses a diagonal-ratio estimate; if no valid estimate exists, the fallback shift is 1e-12."
  ],
  TOPODENSITY: [
    "*TOPODENSITY selects the element field that supplies density-like values for LINEARSTATICTOPO. The referenced field must use the ELEMENT domain with one component.",
    "During stiffness assembly, the element density is combined with the penalization exponent selected by *TOPOEXPONENT. Low density reduces an element's effective stiffness contribution, while high density approaches the full material and section stiffness."
  ],
  TOPOORIENT: [
    "*TOPOORIENT selects the element field that supplies orientation-like values for LINEARSTATICTOPO. The referenced field must use the ELEMENT domain with three components.",
    "The field is available to element and material formulations that evaluate local orientation data. It allows a topology analysis to vary directional stiffness information over the element domain independently of the fixed mesh connectivity."
  ],
  TOPOEXPONENT: [
    "*TOPOEXPONENT sets the density penalization exponent used by LINEARSTATICTOPO. FEMaster scales an element's full stiffness contribution with its density-like field value raised to this exponent.",
    "An exponent of one gives linear density scaling. Larger values penalize intermediate densities more strongly, so low-density elements lose stiffness faster while density values near one remain closer to the full element stiffness."
  ],
  TIME: [
    "*TIME defines the time interval and fixed time step for LINEARTRANSIENT. It accepts either start time, end time, and step size, or only end time and step size; in the two-value form the start time is zero.",
    "The interval length and time step determine the number of integration steps. A smaller step resolves the temporal response more finely but increases the number of effective-system solves and result states. The step should be chosen together with the relevant frequencies and load amplitudes."
  ],
  NEWMARK: [
    "*NEWMARK sets the beta and gamma parameters of the implicit Newmark integration used by LINEARTRANSIENT. FEMaster applies the method to the reduced mass, damping, and stiffness system after Nullspace constraint projection.",
    "The common values beta=0.25 and gamma=0.5 give the average-acceleration method for linear systems. Other values alter stability, numerical damping, and integration accuracy and therefore change the computed time response."
  ],
  DAMPING: [
    "*DAMPING defines damping for LINEARTRANSIENT. The current implementation supports Rayleigh damping, where the damping matrix is a linear combination of the mass and stiffness matrices.",
    "The alpha coefficient multiplies mass and has stronger influence at low frequencies. The beta coefficient multiplies stiffness and has stronger influence at high frequencies. FEMaster forms the same combination after reducing mass and stiffness into the constrained coordinate space."
  ],
  WRITEEVERY: [
    "*WRITEEVERY controls how often LINEARTRANSIENT writes result states. With TYPE=STEPS, the value is rounded to an integer and FEMaster writes every N integration steps.",
    "With TYPE=TIME, the value is a physical output interval. FEMaster converts it to a step stride using the active analysis time step. The final time state is written even when it does not lie exactly on the selected stride."
  ],
  INITIALVELOCITY: [
    "*INITIALVELOCITY selects a nodal field containing six velocity components for LINEARTRANSIENT. FEMaster maps the field into the active degree-of-freedom vector and then projects it into the constrained coordinate space as the initial reduced velocity.",
    "Without this command, initial reduced velocity is zero. Initial acceleration is not read from a field; it is computed internally from the initial displacement and velocity, the active reduced matrices, and the load at the start time."
  ],
  OVERVIEW: [
    "*OVERVIEW prints a summary of the model entities registered at the point where the command is read. It is a diagnostic command for checking the current deck structure and does not create, modify, or activate analysis data.",
    "The overview reports counts and registered objects that help confirm whether geometry, sets, properties, loads, supports, constraints, and loadcases have been read as intended. It can be placed at root level as an inspection checkpoint before solving the model."
  ]
};
