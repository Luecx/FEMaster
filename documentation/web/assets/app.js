const FIG = "../pages/figures/element_schematics/";

const navGroups = [
  {
    title: "Start",
    items: [
      ["overview", "Overview"],
      ["quick-start", "Quick start"],
      ["dsl-basics", "DSL basics"]
    ]
  },
  {
    title: "Model",
    items: [
      ["nodes", "Nodes"],
      ["elements", "Elements"],
      ["surfaces", "Surfaces"],
      ["sets", "Sets"],
      ["materials", "Materials"],
      ["sections", "Sections"],
      ["fields", "Fields"]
    ]
  },
  {
    title: "Commands",
    items: [
      ["commands", "All commands"],
      ["loads-supports", "Loads and supports"],
      ["constraints", "Constraints"],
      ["analysis-controls", "Analysis controls"]
    ]
  },
  {
    title: "Analysis",
    items: [
      ["loadcases", "Loadcases"],
      ["linear-static", "Linear static"],
      ["nonlinear-static", "Nonlinear static"],
      ["modal-buckling", "Modal and buckling"],
      ["transient", "Transient"]
    ]
  },
  {
    title: "Output",
    items: [
      ["results", "Result format"],
      ["examples", "Examples"],
      ["implementation-notes", "Implementation notes"]
    ]
  }
];

const commandCategories = [
  "Geometry",
  "Properties",
  "Loads",
  "Constraints",
  "Loadcase",
  "Transient",
  "Topology",
  "Diagnostics"
];

const commands = [
  {
    id: "NODE",
    name: "NODE",
    category: "Geometry",
    purpose: "Define nodal coordinates and optionally add nodes to a node set.",
    validIn: "root",
    syntax: "*NODE, NSET=<optional-set>\n<id>, <x>, <y>, <z>",
    details: [
      "Node IDs are non-negative row IDs used directly by fields, connectivity, loads, supports and result output.",
      "If NSET is omitted, nodes are added to NALL. Missing coordinates are interpreted as zero, but production decks should write all three coordinates explicitly.",
      "A node does not automatically have six active unknowns. Active DOFs are created by connected elements, constraints, connectors, supports and features."
    ],
    links: ["nodes", "sets"]
  },
  {
    id: "ELEMENT",
    name: "ELEMENT",
    category: "Geometry",
    purpose: "Define element topology, type and optional element-set membership.",
    validIn: "root",
    syntax: "*ELEMENT, TYPE=<type>, ELSET=<optional-set>\n<element-id>, <node-1>, <node-2>, ...",
    details: [
      "TYPE is required. ELSET defaults to EALL.",
      "Connectivity order defines local topology, face numbering, interpolation orientation and, for structural elements, local axes.",
      "Registered types include C3D4, C3D5, C3D6, C3D8, C3D10, C3D13, C3D15, C3D20, C3D20R, S3, S4, MITC4, S6, S8, QSPT, B33 and T33."
    ],
    links: ["elements", "sections"]
  },
  {
    id: "SURFACE",
    name: "SURFACE",
    category: "Geometry",
    purpose: "Define element sides as named surfaces for loads, ties, contact and couplings.",
    validIn: "root",
    syntax: "*SURFACE, NAME=<surface-set>\n<surface-id>, <element-id>, <side>\n\n*SURFACE, NAME=<surface-set>, TYPE=ELEMENT\n<element-set-or-id>, <side>",
    details: [
      "A surface is a geometric side selection, not a material property.",
      "Shell sides can use SPOS and SNEG. Surface normal direction matters for pressure loads and contact gap sign.",
      "TYPE=ELEMENT creates surfaces for the selected side of every element in a set or for one element ID."
    ],
    links: ["surfaces", "constraints"]
  },
  {
    id: "NSET",
    name: "NSET",
    category: "Geometry",
    purpose: "Define a named node set.",
    validIn: "root",
    syntax: "*NSET, NAME=<name>\n<node-id-1>, <node-id-2>, ...",
    details: [
      "Node sets keep decks readable and are accepted by supports, concentrated loads, couplings, contact slaves and point masses.",
      "The same node may belong to several sets."
    ],
    links: ["sets", "nodes"]
  },
  {
    id: "ELSET",
    name: "ELSET",
    category: "Geometry",
    purpose: "Define a named element set.",
    validIn: "root",
    syntax: "*ELSET, NAME=<name>\n<element-id-1>, <element-id-2>, ...",
    details: [
      "Element sets are used by sections, volume loads, inertia relief, RBM constraints, topology fields and surface generation.",
      "Names should reflect modelling intent, for example SOLID_CORE, SPAR_BEAMS or TOPO_DOMAIN."
    ],
    links: ["sets", "elements"]
  },
  {
    id: "SFSET",
    name: "SFSET",
    category: "Geometry",
    purpose: "Define a named surface set.",
    validIn: "root",
    syntax: "*SFSET, NAME=<name>\n<surface-id-1>, <surface-id-2>, ...",
    details: [
      "Surface sets are used by pressure loads, distributed loads, tie constraints, contact and surface couplings.",
      "Keep surface orientation in mind. A correctly named set can still have the opposite normal direction."
    ],
    links: ["sets", "surfaces"]
  },
  {
    id: "MATERIAL",
    name: "MATERIAL",
    category: "Properties",
    purpose: "Create or activate a named material context.",
    validIn: "root",
    syntax: "*MATERIAL, NAME=<name>",
    details: [
      "Material subcommands such as ELASTIC, DENSITY and THERMALEXPANSION apply to the active material.",
      "A material only affects elements after it is referenced by a section command."
    ],
    links: ["materials", "sections"]
  },
  {
    id: "ELASTIC",
    name: "ELASTIC",
    category: "Properties",
    purpose: "Assign elastic material behaviour to the active material.",
    validIn: "material context",
    syntax: "*ELASTIC, TYPE=<ISOTROPIC|ORTHOTROPIC|...>\n<elastic constants>",
    details: [
      "The accepted data line depends on the selected material model.",
      "For common isotropic elasticity, supply Young's modulus and Poisson ratio in a consistent unit system."
    ],
    links: ["materials"]
  },
  {
    id: "DENSITY",
    name: "DENSITY",
    category: "Properties",
    purpose: "Assign mass density to the active material.",
    validIn: "material context",
    syntax: "*DENSITY\n<rho>",
    details: [
      "Density contributes to mass matrices, point-mass combined inertia workflows and transient/modal analyses.",
      "Use units consistent with length, force and time in the deck."
    ],
    links: ["materials", "transient"]
  },
  {
    id: "THERMALEXPANSION",
    name: "THERMALEXPANSION",
    category: "Properties",
    purpose: "Assign thermal expansion data to the active material.",
    validIn: "material context",
    syntax: "*THERMALEXPANSION\n<alpha>",
    details: [
      "Thermal expansion is used by thermal-load workflows that reference temperature fields.",
      "The command requires an active material context."
    ],
    links: ["materials", "fields"]
  },
  {
    id: "SOLIDSECTION",
    name: "SOLIDSECTION",
    category: "Properties",
    purpose: "Assign material and optional orientation to solid elements.",
    validIn: "root",
    syntax: "*SOLIDSECTION, ELSET=<set>, MATERIAL=<material>, ORIENTATION=<optional>",
    details: [
      "Solid elements require a solid section before analysis.",
      "The optional orientation is used by material models that need local axes."
    ],
    links: ["sections", "elements"]
  },
  {
    id: "SHELLSECTION",
    name: "SHELLSECTION",
    category: "Properties",
    purpose: "Assign material, thickness and optional orientation to shell elements.",
    validIn: "root",
    syntax: "*SHELLSECTION, ELSET=<set>, MATERIAL=<material>, THICKNESS=<t>, ORIENTATION=<optional>",
    details: [
      "Shell thickness controls membrane, bending and transverse shear stiffness.",
      "Use shell sides and normals consistently for pressure loads and contact-like workflows."
    ],
    links: ["sections", "elements"]
  },
  {
    id: "BEAMSECTION",
    name: "BEAMSECTION",
    category: "Properties",
    purpose: "Assign material, profile and local-axis data to beam elements.",
    validIn: "root",
    syntax: "*BEAMSECTION, ELSET=<set>, MATERIAL=<material>, PROFILE=<profile>\n<n1x>, <n1y>, <n1z>",
    details: [
      "Beam sections need profile constants and a material.",
      "The local axis vector controls section orientation and therefore bending axes."
    ],
    links: ["sections", "elements"]
  },
  {
    id: "TRUSSSECTION",
    name: "TRUSSSECTION",
    category: "Properties",
    purpose: "Assign material and cross-sectional area to truss elements.",
    validIn: "root",
    syntax: "*TRUSSSECTION, ELSET=<set>, MATERIAL=<material>, AREA=<area>",
    details: [
      "Trusses carry axial force only.",
      "Use beam elements when bending, torsion or rotations are relevant."
    ],
    links: ["sections", "elements"]
  },
  {
    id: "PROFILE",
    name: "PROFILE",
    category: "Properties",
    purpose: "Define beam profile constants.",
    validIn: "root",
    syntax: "*PROFILE, NAME=<name>\n<area>, <Iy>, <Iz>, <J>, ...",
    details: [
      "Beam sections reference profiles by name.",
      "The exact constants must match the beam formulation and unit system."
    ],
    links: ["sections"]
  },
  {
    id: "FIELD",
    name: "FIELD",
    category: "Properties",
    purpose: "Define node, element or integration-point data fields.",
    validIn: "root",
    syntax: "*FIELD, NAME=<name>, DOMAIN=NODE|ELEMENT|ELEMENT_IP, COMPONENTS=<n>\n<row-id>, <value-1>, ...",
    details: [
      "Fields store structured model data such as topology density, orientation, initial velocity or temperature.",
      "Commands that consume fields validate domain and component count."
    ],
    links: ["fields", "analysis-controls"]
  },
  {
    id: "ORIENTATION",
    name: "ORIENTATION",
    category: "Properties",
    purpose: "Define rectangular or cylindrical local coordinate systems.",
    validIn: "root",
    syntax: "*ORIENTATION, NAME=<name>, TYPE=RECTANGULAR|CYLINDRICAL\n<data...>",
    details: [
      "Orientations transform local load, support, material and connector components into global coordinates.",
      "Node coordinates themselves always remain global."
    ],
    links: ["materials", "loads-supports"]
  },
  {
    id: "POINTMASS",
    name: "POINTMASS",
    category: "Properties",
    purpose: "Add concentrated mass, rotary inertia and optional spring stiffness to node sets.",
    validIn: "root",
    syntax: "*POINTMASS, NSET=<set>, MASS=<m>, ...",
    details: [
      "Point masses can contribute translational mass, rotational inertia and optional linear springs.",
      "Inertia relief can optionally include point-mass features."
    ],
    links: ["fields", "analysis-controls"]
  },
  {
    id: "SUPPORT",
    name: "SUPPORT",
    category: "Loads",
    purpose: "Add support entries to a support collector.",
    validIn: "root",
    syntax: "*SUPPORT, SUPPORT_COLLECTOR=<name>, ORIENTATION=<optional>\n<node-set-or-id>, <ux>, <uy>, <uz>, <rx>, <ry>, <rz>",
    details: [
      "Values prescribe kinematic displacement or rotation components. Most fixtures use zero values.",
      "If ORIENTATION is supplied, the six values are interpreted in local axes and transformed to global constraint equations."
    ],
    links: ["loads-supports", "analysis-controls"]
  },
  {
    id: "CLOAD",
    name: "CLOAD",
    category: "Loads",
    purpose: "Add concentrated nodal forces and moments to a load collector.",
    validIn: "root",
    syntax: "*CLOAD, LOAD_COLLECTOR=<name>, ORIENTATION=<optional>, AMPLITUDE=<optional>\n<node-set-or-id>, <Fx>, <Fy>, <Fz>, <Mx>, <My>, <Mz>",
    details: [
      "Load components are conjugate to Ux, Uy, Uz, Rx, Ry, Rz.",
      "Amplitudes are evaluated by transient analyses; static analyses use the defined value."
    ],
    links: ["loads-supports", "transient"]
  },
  {
    id: "DLOAD",
    name: "DLOAD",
    category: "Loads",
    purpose: "Add vector surface tractions to a load collector.",
    validIn: "root",
    syntax: "*DLOAD, LOAD_COLLECTOR=<name>, ORIENTATION=<optional>, AMPLITUDE=<optional>\n<surface-set-or-id>, <Fx>, <Fy>, <Fz>",
    details: [
      "Use DLOAD when the traction direction is known as a vector.",
      "If ORIENTATION is supplied, the vector is interpreted locally and transformed to global coordinates."
    ],
    links: ["loads-supports", "surfaces"]
  },
  {
    id: "PLOAD",
    name: "PLOAD",
    category: "Loads",
    purpose: "Add pressure loads to a load collector.",
    validIn: "root",
    syntax: "*PLOAD, LOAD_COLLECTOR=<name>, AMPLITUDE=<optional>\n<surface-set-or-id>, <pressure>",
    details: [
      "Pressure is converted to traction with the selected surface normal.",
      "If the sign is unexpected, inspect surface side, element orientation and shell SPOS/SNEG selection."
    ],
    links: ["loads-supports", "surfaces"]
  },
  {
    id: "VLOAD",
    name: "VLOAD",
    category: "Loads",
    purpose: "Add volume loads over element regions.",
    validIn: "root",
    syntax: "*VLOAD, LOAD_COLLECTOR=<name>, ORIENTATION=<optional>, AMPLITUDE=<optional>\n<element-set-or-id>, <Fx>, <Fy>, <Fz>",
    details: [
      "VLOAD represents a body-force density over structural elements.",
      "Use consistent units with the model geometry and material density."
    ],
    links: ["loads-supports", "elements"]
  },
  {
    id: "TLOAD",
    name: "TLOAD",
    category: "Loads",
    purpose: "Add thermal loads from a temperature field.",
    validIn: "root",
    syntax: "*TLOAD, LOAD_COLLECTOR=<name>, FIELD=<temperature-field>, REFERENCE=<Tref>",
    details: [
      "Thermal loads require material thermal expansion data and a compatible temperature field.",
      "The reference temperature defines the zero-strain thermal state."
    ],
    links: ["loads-supports", "fields"]
  },
  {
    id: "INERTIALOAD",
    name: "INERTIALOAD",
    category: "Loads",
    purpose: "Add rigid-body inertial loads to a load collector.",
    validIn: "root",
    syntax: "*INERTIALOAD, LOAD_COLLECTOR=<name>, ELSET=<set>, CONSIDER_POINT_MASSES=0|1\n<center>, <acceleration>, <angular terms>",
    details: [
      "The load describes acceleration and angular velocity/acceleration effects over an element set.",
      "It is a load definition, different from the INERTIARELIEF loadcase control."
    ],
    links: ["loads-supports", "analysis-controls"]
  },
  {
    id: "AMPLITUDE",
    name: "AMPLITUDE",
    category: "Loads",
    purpose: "Define a time function for transient load scaling.",
    validIn: "root",
    syntax: "*AMPLITUDE, NAME=<name>, INTERPOLATION=LINEAR|STEP\n<time>, <value>",
    details: [
      "Loads can reference amplitudes. Linear transient analysis evaluates them during time integration.",
      "Without an amplitude, a load uses a constant scale of one."
    ],
    links: ["loads-supports", "transient"]
  },
  {
    id: "COUPLING",
    name: "COUPLING",
    category: "Constraints",
    purpose: "Define kinematic or structural couplings between a master node set and slave nodes or surfaces.",
    validIn: "root",
    syntax: "*COUPLING, MASTER=<master-nset>, TYPE=KINEMATIC|STRUCTURAL, SLAVE=<slave-nset>\n<Ux>, <Uy>, <Uz>, <Rx>, <Ry>, <Rz>\n\n*COUPLING, MASTER=<master-nset>, TYPE=KINEMATIC|STRUCTURAL, SFSET=<surface-set>\n<Ux>, <Uy>, <Uz>, <Rx>, <Ry>, <Rz>",
    details: [
      "KINEMATIC creates constraint equations based on small-rotation rigid-body motion.",
      "STRUCTURAL distributes generalized master loads to slave translational forces in a work-equivalent way.",
      "The master set should normally contain one reference node."
    ],
    links: ["constraints"]
  },
  {
    id: "TIE",
    name: "TIE",
    category: "Constraints",
    purpose: "Bind slave nodes to master surfaces or lines with interpolation equations.",
    validIn: "root",
    syntax: "*TIE, MASTER=<master>, SLAVE=<slave>, ADJUST=NO|YES, DISTANCE=<search-distance>",
    details: [
      "The slave region may be a node set or a surface set.",
      "For each slave node, FEMaster searches the closest master geometry within DISTANCE. If no candidate is found, no equation is generated for that slave node.",
      "ADJUST=YES moves slave coordinates to the projected master position before equations are assembled."
    ],
    links: ["constraints", "surfaces"]
  },
  {
    id: "CONTACT",
    name: "CONTACT",
    category: "Constraints",
    purpose: "Define frictionless node-to-surface penalty contact for nonlinear static analysis.",
    validIn: "root, assembled only by NONLINEARSTATIC",
    syntax: "*CONTACT, MASTER=<master-surface-set>, SLAVE=<slave-node-or-surface-set>, DISTANCE=<search-distance>, PENALTY=<normal-stiffness>, CLEARANCE=<clearance>, FLIP=NO|YES",
    details: [
      "MASTER must be a surface set. SLAVE may be a node set or a surface set; surface slaves are expanded to unique slave nodes.",
      "The normal gap is g = (x_s - x_m) dot n - clearance. Contact activates only for g < 0.",
      "The current implementation is frictionless, translational-only, penalty-based and has no contact pressure result output. The tangent keeps closest-point coordinates and normal fixed."
    ],
    links: ["constraints", "nonlinear-static"]
  },
  {
    id: "CONNECTOR",
    name: "CONNECTOR",
    category: "Constraints",
    purpose: "Define idealized relative-motion constraints between two node sets.",
    validIn: "root",
    syntax: "*CONNECTOR, TYPE=BEAM|HINGE|CYLINDRICAL|TRANSLATOR|JOIN|JOINRX, NSET1=<set>, NSET2=<set>, COORDINATESYSTEM=<orientation>",
    details: [
      "Connector type chooses which local translations and rotations are constrained.",
      "The coordinate system controls the local axes used by the connector equations."
    ],
    links: ["constraints"]
  },
  {
    id: "RBM",
    name: "RBM",
    category: "Constraints",
    purpose: "Generate rigid-body-motion suppression equations over an element set.",
    validIn: "root",
    syntax: "*RBM, ELSET=<elset>",
    details: [
      "RBM creates up to six homogeneous equations for participating nodes.",
      "It is useful for free-body stabilization, diagnostics and inertia relief workflows."
    ],
    links: ["constraints", "analysis-controls"]
  },
  {
    id: "LOADCASE",
    name: "LOADCASE",
    category: "Loadcase",
    purpose: "Start an analysis block and execute it when the block closes.",
    validIn: "root",
    syntax: "*LOADCASE, TYPE=LINEARSTATIC|LINEARSTATICTOPO|NONLINEARSTATIC|EIGENFREQ|LINEARBUCKLING|LINEARTRANSIENT, NAME=<optional-name>\n  <loadcase commands>",
    details: [
      "Loadcases run in input order. Result blocks are written in the same order.",
      "Nested loadcase blocks are not supported."
    ],
    links: ["loadcases"]
  },
  {
    id: "SUPPORTS",
    name: "SUPPORTS",
    category: "Loadcase",
    purpose: "Activate support collectors for the active loadcase.",
    validIn: "LOADCASE",
    syntax: "*SUPPORTS\n<support-collector-1>, <support-collector-2>, ...",
    details: [
      "If no support collector is listed, the solver may use the special all-support collector if it exists.",
      "Active support equations are assembled together with other linear constraint equations."
    ],
    links: ["analysis-controls", "loads-supports"]
  },
  {
    id: "LOADS",
    name: "LOADS",
    category: "Loadcase",
    purpose: "Activate load collectors for the active loadcase.",
    validIn: "LOADCASE",
    syntax: "*LOADS\n<load-collector-1>, <load-collector-2>, ...",
    details: [
      "Load collectors are assembled only when listed.",
      "In transient analysis, amplitudes make the active load vector time-dependent."
    ],
    links: ["analysis-controls", "loads-supports"]
  },
  {
    id: "SOLVER",
    name: "SOLVER",
    category: "Loadcase",
    purpose: "Select solver device and method.",
    validIn: "LOADCASE",
    syntax: "*SOLVER, DEVICE=CPU|GPU, METHOD=DIRECT|INDIRECT",
    details: [
      "DIRECT factorizes the linearized system. INDIRECT uses an iterative path and is most appropriate with NULLSPACE constraints.",
      "NONLINEARSTATIC currently requires METHOD=DIRECT."
    ],
    links: ["analysis-controls"]
  },
  {
    id: "CONSTRAINTMETHOD",
    name: "CONSTRAINTMETHOD",
    category: "Loadcase",
    purpose: "Select NULLSPACE or LAGRANGE constraint handling for supported static loadcases.",
    validIn: "LOADCASE for LINEARSTATIC, LINEARSTATICTOPO and parser-level NONLINEARSTATIC",
    syntax: "*CONSTRAINTMETHOD, TYPE=NULLSPACE|LAGRANGE",
    details: [
      "LINEARSTATIC and LINEARSTATICTOPO can use NULLSPACE or LAGRANGE, with LAGRANGE requiring DIRECT solving.",
      "NONLINEARSTATIC currently runs only with NULLSPACE; selecting LAGRANGE is rejected when the solver starts.",
      "EIGENFREQ, LINEARBUCKLING and LINEARTRANSIENT use NULLSPACE internally and do not accept this command."
    ],
    links: ["analysis-controls", "constraints"]
  },
  {
    id: "NONLINEAR",
    name: "NONLINEAR",
    category: "Loadcase",
    purpose: "Configure NONLINEARSTATIC increment and iteration controls.",
    validIn: "LOADCASE, only NONLINEARSTATIC",
    syntax: "*NONLINEAR, INCREMENTS=<n>, MAXITER=<n>, TOL=<tol>, REGULARIZE_ZERO_ROWS=ON|OFF, REGULARIZATION_ALPHA=<alpha>",
    details: [
      "Defaults are INCREMENTS=10, MAXITER=20, TOL=1e-8, REGULARIZE_ZERO_ROWS=ON and REGULARIZATION_ALPHA=1e-4.",
      "The solver uses full Newton corrections, adaptive load increment growth/cutback and a reduced residual norm.",
      "There is currently no line search, arc-length control, contact damping or automatic penalty scaling."
    ],
    links: ["analysis-controls", "nonlinear-static"]
  },
  {
    id: "INERTIARELIEF",
    name: "INERTIARELIEF",
    category: "Loadcase",
    purpose: "Enable inertia relief for free or nearly free linear static models.",
    validIn: "LOADCASE, LINEARSTATIC and derived static topology path",
    syntax: "*INERTIARELIEF, CONSIDER_POINT_MASSES=0|1",
    details: [
      "Inertia relief cannot be used with active support collectors in the same loadcase.",
      "It adjusts external loads and adds temporary RBM suppression equations for the solve."
    ],
    links: ["analysis-controls"]
  },
  {
    id: "REBALANCELOADS",
    name: "REBALANCELOADS",
    category: "Loadcase",
    purpose: "Balance residual resultant force and moment in linear static loadcases.",
    validIn: "LOADCASE, LINEARSTATIC and derived static topology path",
    syntax: "*REBALANCELOADS",
    details: [
      "This is a numerical load-balancing tool, not a physical inertia model.",
      "Like inertia relief, it requires no active support collectors in the loadcase."
    ],
    links: ["analysis-controls"]
  },
  {
    id: "REQUESTSTIFFNESS",
    name: "REQUESTSTIFFNESS",
    category: "Diagnostics",
    purpose: "Request stiffness matrix output for supported loadcases.",
    validIn: "LOADCASE",
    syntax: "*REQUESTSTIFFNESS, FILE=<base-file>",
    details: [
      "Linear static writes K, A and in the standard static path b.",
      "Buckling writes preload K and A. Nonlinear static writes final Kt and final reduced tangent."
    ],
    links: ["analysis-controls", "results"]
  },
  {
    id: "REQUESTSTGEOM",
    name: "REQUESTSTGEOM",
    category: "Diagnostics",
    purpose: "Request geometric stiffness output for linear buckling.",
    validIn: "LOADCASE, only LINEARBUCKLING",
    syntax: "*REQUESTSTGEOM, FILE=<base-file>",
    details: [
      "Writes geometric stiffness Kg and reduced geometric operator B for buckling diagnostics.",
      "This command is rejected outside LINEARBUCKLING."
    ],
    links: ["analysis-controls", "modal-buckling"]
  },
  {
    id: "CONSTRAINTSUMMARY",
    name: "CONSTRAINTSUMMARY",
    category: "Diagnostics",
    purpose: "Request a generated-constraint diagnostic report.",
    validIn: "LOADCASE",
    syntax: "*CONSTRAINTSUMMARY",
    details: [
      "The report groups supports, RBM equations, couplings, ties, connectors and other equations.",
      "It is diagnostic output only and does not change the model."
    ],
    links: ["constraints", "analysis-controls"]
  },
  {
    id: "NUMEIGENVALUES",
    name: "NUMEIGENVALUES",
    category: "Loadcase",
    purpose: "Set the number of requested eigenpairs.",
    validIn: "LOADCASE, EIGENFREQ or LINEARBUCKLING",
    syntax: "*NUMEIGENVALUES\n<count>",
    details: [
      "The count must be positive.",
      "The solver clamps the requested count to the available reduced system size."
    ],
    links: ["modal-buckling"]
  },
  {
    id: "SIGMA",
    name: "SIGMA",
    category: "Loadcase",
    purpose: "Set the buckling eigenvalue shift.",
    validIn: "LOADCASE, only LINEARBUCKLING",
    syntax: "*SIGMA\n<shift>",
    details: [
      "If SIGMA is zero or omitted, FEMaster estimates a shift from the preload response or falls back to a conservative value.",
      "A poor shift can guide the eigensolver to an unintended part of the spectrum."
    ],
    links: ["modal-buckling", "analysis-controls"]
  },
  {
    id: "TOPODENSITY",
    name: "TOPODENSITY",
    category: "Topology",
    purpose: "Select the element density field for topology static analysis.",
    validIn: "LOADCASE, only LINEARSTATICTOPO",
    syntax: "*TOPODENSITY, FIELD=<element-field>",
    details: [
      "The field must use ELEMENT domain and one component.",
      "The solver applies rho^p stiffness scaling using TOPOEXPONENT."
    ],
    links: ["loadcases", "fields"]
  },
  {
    id: "TOPOORIENT",
    name: "TOPOORIENT",
    category: "Topology",
    purpose: "Select the element orientation field for topology static analysis.",
    validIn: "LOADCASE, only LINEARSTATICTOPO",
    syntax: "*TOPOORIENT, FIELD=<element-field>",
    details: [
      "The field must use ELEMENT domain and three components.",
      "It is passed to the model as material orientation data during stiffness assembly."
    ],
    links: ["loadcases", "fields"]
  },
  {
    id: "TOPOEXPONENT",
    name: "TOPOEXPONENT",
    category: "Topology",
    purpose: "Set SIMP-like density penalization exponent.",
    validIn: "LOADCASE, only LINEARSTATICTOPO",
    syntax: "*TOPOEXPONENT\n<p>",
    details: [
      "The active stiffness scale is rho^p.",
      "Values greater than one penalize intermediate densities."
    ],
    links: ["loadcases"]
  },
  {
    id: "TIME",
    name: "TIME",
    category: "Transient",
    purpose: "Set transient time interval and fixed step size.",
    validIn: "LOADCASE, only LINEARTRANSIENT",
    syntax: "*TIME\n<t_start>, <t_end>, <dt>\n\n*TIME\n<t_end>, <dt>",
    details: [
      "Two values use t_start=0.",
      "The solver writes the final state even if output cadence does not fall exactly on it."
    ],
    links: ["transient", "analysis-controls"]
  },
  {
    id: "NEWMARK",
    name: "NEWMARK",
    category: "Transient",
    purpose: "Set Newmark-beta integration parameters.",
    validIn: "LOADCASE, only LINEARTRANSIENT",
    syntax: "*NEWMARK\n<beta>, <gamma>",
    details: [
      "Defaults in the transient loadcase are beta=0.25 and gamma=0.5.",
      "Other values change numerical damping, accuracy and stability."
    ],
    links: ["transient", "analysis-controls"]
  },
  {
    id: "DAMPING",
    name: "DAMPING",
    category: "Transient",
    purpose: "Set Rayleigh damping for transient analysis.",
    validIn: "LOADCASE, only LINEARTRANSIENT",
    syntax: "*DAMPING, TYPE=RAYLEIGH\n<alpha>, <beta>",
    details: [
      "Rayleigh damping is assembled in reduced space as Cr = alpha Mr + beta A.",
      "alpha is mass-proportional; beta is stiffness-proportional."
    ],
    links: ["transient", "analysis-controls"]
  },
  {
    id: "WRITEEVERY",
    name: "WRITEEVERY",
    category: "Transient",
    purpose: "Control transient result cadence.",
    validIn: "LOADCASE, only LINEARTRANSIENT",
    syntax: "*WRITEEVERY, TYPE=STEPS|TIME\n<value>",
    details: [
      "TYPE=STEPS writes every rounded integer number of steps.",
      "TYPE=TIME converts a physical output interval to a step stride."
    ],
    links: ["transient", "results"]
  },
  {
    id: "INITIALVELOCITY",
    name: "INITIALVELOCITY",
    category: "Transient",
    purpose: "Set transient initial velocity from a node field.",
    validIn: "LOADCASE, only LINEARTRANSIENT",
    syntax: "*INITIALVELOCITY, FIELD=<node-field>",
    details: [
      "The referenced field must have NODE domain and exactly six components.",
      "The field is reduced into active coordinates and used as the initial velocity."
    ],
    links: ["transient", "fields"]
  },
  {
    id: "OVERVIEW",
    name: "OVERVIEW",
    category: "Diagnostics",
    purpose: "Print a model overview.",
    validIn: "root",
    syntax: "*OVERVIEW",
    details: [
      "Useful for checking model sizes, registered entities and high-level deck state before loadcases run.",
      "It is a diagnostic command and does not modify the model."
    ],
    links: ["overview"]
  }
];

const elements = [
  ["C3D4", "4-node tetrahedral solid", "Linear tetrahedra are robust for complex meshes but coarse for bending and stress gradients.", "c3d4_schematic.png"],
  ["C3D5", "5-node pyramid input", "Pyramid-style transition input, mapped internally to a degenerated hexahedral representation.", null],
  ["C3D6", "6-node wedge solid", "Useful for swept triangular meshes, layers and transition regions.", "c3d6_schematic.png"],
  ["C3D8", "8-node hexahedral solid", "Good default for structured solid meshes with regular block-like topology.", "c3d8_schematic.png"],
  ["C3D10", "10-node quadratic tetrahedral solid", "Preferred tetrahedral choice for stress work on complex geometry.", "c3d10_schematic.png"],
  ["C3D13", "13-node quadratic pyramid solid", "Higher-order transition element for mixed hex/tet meshes.", "c3d13_schematic.png"],
  ["C3D15", "15-node quadratic wedge solid", "Higher-order wedge for swept and layered solid regions.", "c3d15_schematic.png"],
  ["C3D20", "20-node quadratic hexahedral solid", "Higher-order hexahedron with full higher-order integration.", "c3d20_schematic.png"],
  ["C3D20R", "20-node reduced-integration hexahedral solid", "Quadratic hexahedron with reduced stiffness integration.", "c3d20r_schematic.png"],
  ["S3", "3-node triangular shell", "Low-order triangular shell for simple shell regions and transitions.", "s3_schematic.png"],
  ["S4", "4-node quadrilateral shell", "Standard low-order quadrilateral shell.", "s4_schematic.png"],
  ["MITC4", "4-node MITC shell", "Quadrilateral shell with MITC-style tied shear interpolation.", "mitc4_schematic.png"],
  ["S6", "6-node quadratic triangular shell", "Higher-order triangular shell.", "s6_schematic.png"],
  ["S8", "8-node quadratic quadrilateral shell", "Higher-order quadrilateral shell.", "s8_schematic.png"],
  ["QSPT", "4-node shear-panel type", "Specialized quadrilateral shear-panel formulation.", "qspt_schematic.png"],
  ["B33", "3D beam", "Beam element for axial, bending, shear and torsional response.", "b33_schematic.png"],
  ["T33", "3D truss", "Axial-only line element.", "t33_schematic.png"]
];

const pages = {
  "overview": {
    section: "Start",
    title: "FEMaster Documentation",
    subtitle: "A navigable web reference for the FEMaster input deck, model objects, solver controls and result fields.",
    render: () => `
      <div class="grid three">
        ${card("Command catalog", "Every registered command has syntax, scope, notes and related links.", "commands")}
        ${card("Loadcase solver map", "Static, nonlinear, modal, buckling and transient paths explained against the implementation.", "loadcases")}
        ${card("Element gallery", "Solid, shell, beam and truss families with the existing schematic images.", "elements")}
      </div>
      <div class="callout">
        This web documentation is intentionally static. It can be opened directly from <code>documentation/web/index.html</code> and does not require Node, a generator or a server.
      </div>
      <h2>What the web docs add</h2>
      <div class="grid two">
        ${plainCard("Fast navigation", "The left tree mirrors the manual structure but keeps command detail pages one click away. Hash routes make links stable.")}
        ${plainCard("Command-first workflow", "The catalog is generated from command metadata in this file, so search, cards and detail pages stay consistent.")}
        ${plainCard("Solver caveats close to controls", "Important implementation limits like nonlinear DIRECT/NULLSPACE requirements are visible next to the relevant commands.")}
        ${plainCard("Image reuse", "The element pages reference the existing PNG schematics used by the LaTeX manual.")}
      </div>
    `
  },
  "quick-start": {
    section: "Start",
    title: "Quick start",
    subtitle: "A minimal path from input deck concepts to a solvable static model.",
    render: () => `
      <h2>Small static deck skeleton</h2>
      <pre><code>*NODE, NSET=NALL
0, 0, 0, 0
1, 1, 0, 0

*ELEMENT, TYPE=T33, ELSET=BAR
0, 0, 1

*MATERIAL, NAME=STEEL
  *ELASTIC, TYPE=ISOTROPIC
  210000, 0.3

*TRUSSSECTION, ELSET=BAR, MATERIAL=STEEL, AREA=1.0

*SUPPORT, SUPPORT_COLLECTOR=FIXED
0, 0, 0, 0, , ,

*CLOAD, LOAD_COLLECTOR=FORCE
1, 1000, 0, 0, 0, 0, 0

*LOADCASE, TYPE=LINEARSTATIC, NAME=bar-static
  *SUPPORTS
  FIXED
  *LOADS
  FORCE
  *SOLVER, DEVICE=CPU, METHOD=DIRECT
  *CONSTRAINTMETHOD, TYPE=NULLSPACE</code></pre>
      <h2>Model order</h2>
      <ol>
        <li>Define geometry: nodes, elements, surfaces and sets.</li>
        <li>Define properties: materials, sections, profiles, orientations and fields.</li>
        <li>Define loads, supports and constraints into collectors.</li>
        <li>Run one or more loadcases that activate collectors and solver controls.</li>
      </ol>
      <div class="link-list">
        ${pageLink("nodes", "Nodes")}
        ${pageLink("commands", "All commands")}
        ${pageLink("linear-static", "Linear static analysis")}
      </div>
    `
  },
  "dsl-basics": {
    section: "Start",
    title: "DSL basics",
    subtitle: "Command scopes, keywords, data lines and naming conventions.",
    render: () => `
      <div class="grid two">
        ${plainCard("Commands", "Commands start with an asterisk. Keywords are written on the command line, data follows on one or more lines. Command and keyword names are normalized by the parser.")}
        ${plainCard("Scopes", "MATERIAL opens a material context. LOADCASE opens an analysis context. Commands such as SUPPORTS, LOADS, SOLVER and supported CONSTRAINTMETHOD variants configure the active loadcase.")}
        ${plainCard("Targets", "Many commands accept either a set name or a numeric ID, for example CLOAD node targets, DLOAD surface targets and VLOAD element targets.")}
        ${plainCard("Component order", "Mechanical nodal vectors use Ux, Uy, Uz, Rx, Ry, Rz. Not every element contributes all six components.")}
      </div>
    `
  },
  "nodes": topicPage("Model", "Nodes", "Nodes are geometric reference points with global coordinates and active DOFs derived from connected model objects.", `
    <pre><code>*NODE, NSET=&lt;optional-set&gt;
&lt;id&gt;, &lt;x&gt;, &lt;y&gt;, &lt;z&gt;</code></pre>
    <p>Node IDs are non-negative integer row identifiers. They are used directly by connectivity, fields, supports, loads and result output.</p>
    <p>A node may have up to six mechanical components, <code>Ux, Uy, Uz, Rx, Ry, Rz</code>. Solid-only nodes usually activate only translations. Shells, beams, connectors, couplings, point-mass springs and supports can make rotations relevant.</p>
  `),
  "surfaces": topicPage("Model", "Surfaces", "Surfaces are element-side selections used by loads, ties, contact and surface couplings.", `
    <pre><code>*SURFACE, NAME=&lt;surface-set&gt;
&lt;surface-id&gt;, &lt;element-id&gt;, &lt;side&gt;

*SURFACE, NAME=&lt;surface-set&gt;, TYPE=ELEMENT
&lt;element-set-or-id&gt;, &lt;side&gt;</code></pre>
    <p>Surface orientation matters. Pressure direction, contact gap sign and shell positive/negative sides all depend on the selected element side and normal.</p>
  `),
  "sets": topicPage("Model", "Sets", "Sets group nodes, elements and surfaces so decks remain readable and loadcases can activate meaningful regions.", `
    <div class="table-wrap"><table><thead><tr><th>Set</th><th>Command</th><th>Typical use</th></tr></thead><tbody>
      <tr><td>Node set</td><td><code>*NSET</code></td><td>Supports, concentrated loads, couplings, slave nodes, point masses.</td></tr>
      <tr><td>Element set</td><td><code>*ELSET</code></td><td>Sections, volume loads, RBM, topology fields, generated surfaces.</td></tr>
      <tr><td>Surface set</td><td><code>*SFSET</code></td><td>Pressure, traction, ties, contact, surface coupling.</td></tr>
    </tbody></table></div>
    <p>One entity may belong to several sets. That is normal. Problems arise only when later commands impose contradictory supports, duplicate loads or unintended overlapping constraints.</p>
  `),
  "materials": topicPage("Model", "Materials", "Materials collect elastic, density and thermal data before sections assign them to elements.", `
    <pre><code>*MATERIAL, NAME=&lt;name&gt;
  *ELASTIC, TYPE=ISOTROPIC
  &lt;E&gt;, &lt;nu&gt;
  *DENSITY
  &lt;rho&gt;</code></pre>
    <p>A material alone does not affect the model. It must be referenced by a section command such as <code>*SOLIDSECTION</code>, <code>*SHELLSECTION</code>, <code>*BEAMSECTION</code> or <code>*TRUSSSECTION</code>.</p>
  `),
  "sections": topicPage("Model", "Sections", "Sections assign material and geometric properties to element sets.", `
    <div class="grid two">
      ${plainCard("SOLIDSECTION", "Assigns material and optional orientation to solid elements.")}
      ${plainCard("SHELLSECTION", "Assigns material, thickness and optional orientation to shell elements.")}
      ${plainCard("BEAMSECTION", "Assigns material, profile and local axis data to beams.")}
      ${plainCard("TRUSSSECTION", "Assigns material and cross-sectional area to trusses.")}
    </div>
    <p>Elements without a compatible section do not have enough physical data for structural assembly.</p>
  `),
  "fields": topicPage("Model", "Fields", "Fields store node, element or integration-point data consumed by loadcase controls and model features.", `
    <div class="table-wrap"><table><thead><tr><th>Use</th><th>Domain</th><th>Components</th></tr></thead><tbody>
      <tr><td>Topology density</td><td>ELEMENT</td><td>1</td></tr>
      <tr><td>Topology orientation</td><td>ELEMENT</td><td>3</td></tr>
      <tr><td>Initial velocity</td><td>NODE</td><td>6</td></tr>
      <tr><td>Temperature</td><td>NODE or element-related</td><td>model dependent</td></tr>
    </tbody></table></div>
  `),
  "commands": {
    section: "Commands",
    title: "All commands",
    subtitle: "Search, filter and open detailed command pages.",
    render: renderCommandIndex
  },
  "loads-supports": topicPage("Commands", "Loads and supports", "Loads and supports are collected first and activated later by loadcases.", `
    <p>Load collectors are never assembled unless listed in <code>*LOADS</code>. Support collectors are normally selected with <code>*SUPPORTS</code>; if no support collector is listed, the solver may use the special all-support collector if present.</p>
    <div class="grid two">
      ${commandSummary("SUPPORT")}
      ${commandSummary("CLOAD")}
      ${commandSummary("DLOAD")}
      ${commandSummary("PLOAD")}
      ${commandSummary("VLOAD")}
      ${commandSummary("TLOAD")}
    </div>
  `),
  "constraints": topicPage("Commands", "Constraints and interactions", "Linear constraints generate equations. Contact contributes nonlinear force and tangent stiffness.", `
    <div class="grid two">
      ${commandSummary("RBM")}
      ${commandSummary("COUPLING")}
      ${commandSummary("TIE")}
      ${commandSummary("CONTACT")}
      ${commandSummary("CONNECTOR")}
    </div>
    <div class="callout warning">Contact is not a linear constraint equation. It is assembled only by <code>NONLINEARSTATIC</code> and is sensitive to normal orientation, penalty stiffness and load increment size.</div>
  `),
  "analysis-controls": topicPage("Commands", "Analysis controls", "Loadcase commands select collectors, solver method, constraint method and analysis-specific options.", `
    <div class="grid two">
      ${commandSummary("SUPPORTS")}
      ${commandSummary("LOADS")}
      ${commandSummary("SOLVER")}
      ${commandSummary("CONSTRAINTMETHOD")}
      ${commandSummary("NONLINEAR")}
      ${commandSummary("REQUESTSTIFFNESS")}
    </div>
  `),
  "elements": {
    section: "Model",
    title: "Elements",
    subtitle: "Registered element types with schematic images reused from the LaTeX manual.",
    render: renderElements
  },
  "loadcases": topicPage("Analysis", "Loadcases", "A loadcase decides which equation is assembled, which collectors are active and how constraints are handled.", `
    <div class="grid two">
      ${loadcaseCard("LINEARSTATIC", "Small-displacement static equilibrium with optional NULLSPACE or LAGRANGE constraints.", "linear-static")}
      ${loadcaseCard("LINEARSTATICTOPO", "Linear static solve with density and orientation fields scaling element stiffness.", "loadcases")}
      ${loadcaseCard("NONLINEARSTATIC", "Incremental Newton solve with updated positions and penalty contact support.", "nonlinear-static")}
      ${loadcaseCard("EIGENFREQ", "Natural frequencies and mode shapes using reduced stiffness and mass operators.", "modal-buckling")}
      ${loadcaseCard("LINEARBUCKLING", "Static preload, geometric stiffness and generalized buckling eigenproblem.", "modal-buckling")}
      ${loadcaseCard("LINEARTRANSIENT", "Implicit Newmark time integration with optional Rayleigh damping.", "transient")}
    </div>
  `),
  "linear-static": topicPage("Analysis", "Linear static", "The standard path for displacement, reactions, stress and static load paths.", `
    <p>Physical equation: <code>K u = f</code>. With Nullspace constraints FEMaster solves <code>T^T K T q = T^T(f - K u_p)</code> and recovers <code>u = u_p + T q</code>.</p>
    <pre><code>*LOADCASE, TYPE=LINEARSTATIC, NAME=static-load
  *SUPPORTS
  FIXED
  *LOADS
  FORCE
  *SOLVER, DEVICE=CPU, METHOD=DIRECT
  *CONSTRAINTMETHOD, TYPE=NULLSPACE</code></pre>
    <p>LAGRANGE is available for linear static direct solves. INDIRECT should be used with NULLSPACE.</p>
  `),
  "nonlinear-static": topicPage("Analysis", "Nonlinear static", "Incremental residual iterations with tangent stiffness, updated coordinates and contact assembly.", `
    <p>At load factor <code>lambda</code>, the residual is <code>lambda f_total - f_int(u)</code>. The reduced Newton equation is <code>T^T K_t T dq = T^T r</code>.</p>
    <pre><code>*LOADCASE, TYPE=NONLINEARSTATIC, NAME=contact-step
  *SUPPORTS
  FIXED
  *LOADS
  PUSH
  *SOLVER, DEVICE=CPU, METHOD=DIRECT
  *CONSTRAINTMETHOD, TYPE=NULLSPACE
  *NONLINEAR, INCREMENTS=50, MAXITER=30, TOL=1e-8</code></pre>
    <div class="callout warning">Current limits: DIRECT solver only, NULLSPACE constraints only, full Newton corrections, no line search or arc-length control.</div>
    ${commandSummary("CONTACT")}
  `),
  "modal-buckling": topicPage("Analysis", "Modal and buckling", "Eigenfrequency and linear buckling both use Nullspace constraints internally.", `
    <h2>Eigenfrequency</h2>
    <p>Solves <code>T^T K T psi = omega^2 T^T M T psi</code>. Loads are not used.</p>
    <h2>Linear buckling</h2>
    <p>Solves a static preload, builds geometric stiffness <code>K_G</code>, then solves <code>A phi = lambda (-B) phi</code> with <code>A=T^TKT</code> and <code>B=T^TK_G T</code>.</p>
    <div class="grid two">
      ${commandSummary("NUMEIGENVALUES")}
      ${commandSummary("SIGMA")}
      ${commandSummary("REQUESTSTGEOM")}
    </div>
  `),
  "transient": topicPage("Analysis", "Linear transient", "Implicit Newmark integration over reduced stiffness, mass and damping operators.", `
    <p>Reduced equation: <code>M_r qddot + C_r qdot + K_r q = r(t)</code>. Loads with amplitudes are evaluated over time.</p>
    <pre><code>*LOADCASE, TYPE=LINEARTRANSIENT, NAME=transient
  *SUPPORTS
  FIXED
  *LOADS
  PULSE
  *TIME
  0, 1.0, 0.001
  *NEWMARK
  0.25, 0.5
  *DAMPING, TYPE=RAYLEIGH
  0.0, 1e-4
  *WRITEEVERY, TYPE=STEPS
  10</code></pre>
    <div class="grid two">${commandSummary("TIME")}${commandSummary("NEWMARK")}${commandSummary("DAMPING")}${commandSummary("INITIALVELOCITY")}</div>
  `),
  "results": topicPage("Output", "Result format", "FEMaster writes text-based .res files organized into loadcases and fields.", `
    <pre><code>LOADCASE &lt;id&gt; &lt;name-or-type&gt;
FIELD &lt;name&gt; TYPE=&lt;type&gt; ROWS=&lt;rows&gt; COLS=&lt;cols&gt;
&lt;values...&gt;</code></pre>
    <div class="table-wrap"><table><thead><tr><th>Loadcase</th><th>Typical fields</th></tr></thead><tbody>
      <tr><td>Linear Static</td><td>DISPLACEMENT, STRAIN, STRESS, STRESS_TOP, STRESS_BOT, SHELL_RESULTANTS, EXTERNAL_FORCES, REACTION_FORCES, LOCAL_SECTION_FORCES, optional SHEAR_FLOW.</td></tr>
      <tr><td>Nonlinear Static</td><td>DISPLACEMENT, STRAIN, STRESS, STRESS_TOP, STRESS_BOT, EXTERNAL_FORCES, INTERNAL_FORCES, REACTION_FORCES.</td></tr>
      <tr><td>Eigenfrequency</td><td>MODE_SHAPE_i, PARTICIPATION_i, EIGENVALUES, EIGENFREQUENCIES, FREQUENCIES.</td></tr>
      <tr><td>Buckling</td><td>BUCKLING_MODE_i, BUCKLING_FACTORS.</td></tr>
      <tr><td>Transient</td><td>DISPLACEMENT_k, VELOCITY_k, ACCELERATION_k.</td></tr>
    </tbody></table></div>
  `),
  "examples": topicPage("Output", "Examples", "Starter patterns for common analyses.", `
    <div class="grid two">
      ${plainCard("Static sanity check", "Start with LINEARSTATIC, DIRECT, NULLSPACE and request stiffness if the model is singular or over-constrained.")}
      ${plainCard("Contact debug step", "Use NONLINEARSTATIC with small INCREMENTS, moderate PENALTY and explicit CONTACT normal verification.")}
      ${plainCard("Buckling preload", "Validate the static preload first, then run LINEARBUCKLING with NUMEIGENVALUES and optional SIGMA.")}
      ${plainCard("Transient response", "Define AMPLITUDE curves, TIME, NEWMARK and WRITEEVERY before increasing model size.")}
    </div>
  `),
  "implementation-notes": topicPage("Output", "Implementation notes", "Important solver behaviours that are easy to miss from syntax alone.", `
    <div class="grid two">
      ${plainCard("Active DOFs", "The active index matrix is built from model objects. Nodes alone do not create all six DOFs.")}
      ${plainCard("Constraint reports", "CONSTRAINTSUMMARY only prints generated equations. It does not change the equations or solve.")}
      ${plainCard("Nonlinear positions", "NONLINEARSTATIC restores model positions after the run. Current positions are updated during assembly.")}
      ${plainCard("Contact tangent", "Current HEAD contact tangent is a symmetric penalty outer product with fixed closest-point coordinates and normal.")}
    </div>
  `)
};

function topicPage(section, title, subtitle, body) {
  return { section, title, subtitle, render: () => body };
}

function card(title, text, pageId) {
  return `<a class="card" href="#/${pageId}"><h3>${title}</h3><p class="muted">${text}</p></a>`;
}

function plainCard(title, text) {
  return `<div class="card"><h3>${title}</h3><p class="muted">${text}</p></div>`;
}

function loadcaseCard(name, text, pageId) {
  return `<a class="card" href="#/${pageId}"><span class="pill">${name}</span><h3>${name}</h3><p class="muted">${text}</p></a>`;
}

function pageLink(pageId, text) {
  return `<a href="#/${pageId}">${text}</a>`;
}

function commandSummary(id) {
  const cmd = commands.find((item) => item.id === id);
  if (!cmd) return "";
  return `<a class="card" href="#/command/${cmd.id}"><span class="command-name">*${cmd.name}</span><p>${cmd.purpose}</p><p class="muted">${cmd.validIn}</p></a>`;
}

function renderCommandIndex() {
  const activeCategory = state.commandCategory;
  const categories = ["All", ...commandCategories];
  const filtered = commands.filter((cmd) => activeCategory === "All" || cmd.category === activeCategory);
  return `
    <div class="command-toolbar">
      ${categories.map((cat) => `<button class="filter-button ${cat === activeCategory ? "active" : ""}" type="button" data-category="${cat}">${cat}</button>`).join("")}
    </div>
    <div class="command-grid">
      ${filtered.map(renderCommandCard).join("")}
    </div>
  `;
}

function renderCommandCard(cmd) {
  return `
    <a class="command-card" href="#/command/${cmd.id}">
      <div class="pill-row"><span class="pill">${cmd.category}</span><span class="pill warn">${cmd.validIn}</span></div>
      <h3 class="command-name">*${cmd.name}</h3>
      <p>${cmd.purpose}</p>
    </a>
  `;
}

function renderCommandDetail(cmd) {
  return `
    <div class="hero">
      <div class="eyebrow">Command / ${cmd.category}</div>
      <h1>*${cmd.name}</h1>
      <p class="lead">${cmd.purpose}</p>
      <div class="pill-row"><span class="pill">${cmd.category}</span><span class="pill warn">${cmd.validIn}</span></div>
    </div>
    <h2>Syntax</h2>
    <pre><code>${escapeHtml(cmd.syntax)}</code></pre>
    <h2>Details</h2>
    <ul>${cmd.details.map((item) => `<li>${item}</li>`).join("")}</ul>
    <h2>Related topics</h2>
    <div class="link-list">${cmd.links.map((id) => pageLink(id, pages[id]?.title || id)).join("")}</div>
  `;
}

function renderElements() {
  return `
    <div class="element-grid">
      ${elements.map(([name, title, text, image]) => `
        <article class="element-card">
          ${image ? `<img src="${FIG}${image}" alt="${name} schematic">` : `<div class="empty">No schematic yet</div>`}
          <div class="element-card-body">
            <span class="pill">${name}</span>
            <h3>${title}</h3>
            <p class="muted">${text}</p>
            <pre><code>*ELEMENT, TYPE=${name}, ELSET=&lt;set&gt;
&lt;id&gt;, &lt;connectivity...&gt;</code></pre>
          </div>
        </article>
      `).join("")}
    </div>
  `;
}

function renderSearch(query) {
  const q = query.trim().toLowerCase();
  const pageHits = Object.entries(pages)
    .filter(([id, page]) => `${id} ${page.title} ${page.subtitle || ""}`.toLowerCase().includes(q));
  const commandHits = commands
    .filter((cmd) => `${cmd.name} ${cmd.category} ${cmd.purpose} ${cmd.details.join(" ")}`.toLowerCase().includes(q));
  const elementHits = elements
    .filter((element) => element.join(" ").toLowerCase().includes(q));

  return `
    <div class="hero">
      <div class="eyebrow">Search</div>
      <h1>Search results</h1>
      <p class="lead">Matches for <code>${escapeHtml(query)}</code>.</p>
    </div>
    <h2>Commands</h2>
    <div class="command-grid">${commandHits.length ? commandHits.map(renderCommandCard).join("") : `<div class="empty">No command matches.</div>`}</div>
    <h2>Topics</h2>
    <div class="grid two">${pageHits.length ? pageHits.map(([id, page]) => card(page.title, page.subtitle || "", id)).join("") : `<div class="empty">No topic matches.</div>`}</div>
    <h2>Elements</h2>
    <div class="grid two">${elementHits.length ? elementHits.map(([name, title, text]) => plainCard(name, `${title}. ${text}`)).join("") : `<div class="empty">No element matches.</div>`}</div>
  `;
}

function escapeHtml(value) {
  return value.replace(/[&<>"']/g, (char) => ({
    "&": "&amp;",
    "<": "&lt;",
    ">": "&gt;",
    '"': "&quot;",
    "'": "&#039;"
  }[char]));
}

const state = {
  commandCategory: "All",
  search: ""
};

function buildSidebar() {
  const nav = document.getElementById("tree-nav");
  nav.innerHTML = navGroups.map((group, index) => `
    <section class="tree-section" data-open="${index < 3 ? "true" : "false"}">
      <button class="tree-heading" type="button">${group.title}</button>
      <div class="tree-items">
        ${group.items.map(([id, label]) => `<a class="tree-link" href="#/${id}" data-page-id="${id}"><span>${label}</span><span class="tree-kicker">${pages[id]?.section || ""}</span></a>`).join("")}
      </div>
    </section>
  `).join("");

  nav.querySelectorAll(".tree-heading").forEach((button) => {
    button.addEventListener("click", () => {
      const section = button.closest(".tree-section");
      section.dataset.open = section.dataset.open === "true" ? "false" : "true";
    });
  });
}

function parseRoute() {
  const raw = location.hash.replace(/^#\/?/, "");
  return raw || "overview";
}

function setActiveNav(route) {
  document.querySelectorAll(".tree-link").forEach((link) => {
    link.classList.toggle("active", link.dataset.pageId === route);
  });
}

function render() {
  const pageEl = document.getElementById("page");
  const breadcrumbs = document.getElementById("breadcrumbs");

  if (state.search.trim()) {
    pageEl.innerHTML = renderSearch(state.search);
    breadcrumbs.textContent = "Search";
    setActiveNav("");
    pageEl.focus();
    return;
  }

  const route = parseRoute();
  if (route.startsWith("command/")) {
    const id = route.split("/")[1]?.toUpperCase();
    const cmd = commands.find((item) => item.id === id);
    if (!cmd) {
      pageEl.innerHTML = `<div class="empty">Unknown command.</div>`;
      breadcrumbs.textContent = "Command";
      return;
    }
    pageEl.innerHTML = renderCommandDetail(cmd);
    breadcrumbs.textContent = `Commands / *${cmd.name}`;
    setActiveNav("commands");
    pageEl.focus();
    return;
  }

  const page = pages[route] || pages.overview;
  pageEl.innerHTML = `
    <div class="hero">
      <div class="eyebrow">${page.section}</div>
      <h1>${page.title}</h1>
      ${page.subtitle ? `<p class="lead">${page.subtitle}</p>` : ""}
    </div>
    ${page.render()}
  `;
  breadcrumbs.textContent = `${page.section} / ${page.title}`;
  setActiveNav(route);
  pageEl.focus();
}

function initEvents() {
  document.getElementById("global-search").addEventListener("input", (event) => {
    state.search = event.target.value;
    render();
  });

  document.getElementById("nav-toggle").addEventListener("click", () => {
    document.body.classList.toggle("nav-open");
  });

  document.addEventListener("click", (event) => {
    if (event.target.closest("a[href^='#/']")) {
      document.body.classList.remove("nav-open");
    }
    const filter = event.target.closest("[data-category]");
    if (filter) {
      state.commandCategory = filter.dataset.category;
      render();
    }
  });

  window.addEventListener("hashchange", () => {
    state.search = "";
    document.getElementById("global-search").value = "";
    render();
  });
}

buildSidebar();
initEvents();
render();
