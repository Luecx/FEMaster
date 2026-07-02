// ─────────────────────────────────────────────────────────────────────────────
// FEMaster Documentation — app.js
// Combines the original data + fully rewritten rendering/nav/theme logic
// ─────────────────────────────────────────────────────────────────────────────
const sitePages = {
  "overview": { title: "FEMaster Documentation", section: "Start", file: "index.html" },
  "installation": { title: "Installation guide", section: "Start", file: "pages/installation.html" },
  "quick-start": { title: "Quick start", section: "Start", file: "pages/quick-start.html" },
  "dsl": { title: "DSL", section: "DSL", file: "pages/dsl/index.html" },
  "dsl-model": { title: "Model", section: "DSL", file: "pages/dsl/model/index.html" },
  "dsl-properties": { title: "Properties", section: "DSL", file: "pages/dsl/properties/index.html" },
  "dsl-loads": { title: "Loads and supports", section: "DSL", file: "pages/dsl/loads-and-supports/index.html" },
  "dsl-constraints": { title: "Constraints and interactions", section: "DSL", file: "pages/dsl/constraints/index.html" },
  "dsl-analysis": { title: "Analysis", section: "DSL", file: "pages/dsl/analysis/index.html" },
  "results": { title: "Result format", section: "Reference", file: "pages/results.html" },
  "examples": { title: "Examples", section: "Reference", file: "pages/examples.html" },
  "implementation-notes": { title: "Implementation notes", section: "Reference", file: "pages/implementation-notes.html" }
};

const commandGroups = [
  {
    id: "dsl-model",
    title: "Model",
    folder: "model",
    commandIds: ["NODE", "ELEMENT", "SURFACE", "NSET", "ELSET", "SFSET"]
  },
  {
    id: "dsl-properties",
    title: "Properties",
    folder: "properties",
    commandIds: ["MATERIAL", "ELASTIC", "DENSITY", "THERMALEXPANSION", "SOLIDSECTION", "SHELLSECTION", "BEAMSECTION", "TRUSSSECTION", "PROFILE", "FIELD", "ORIENTATION", "POINTMASS"]
  },
  {
    id: "dsl-loads",
    title: "Loads and supports",
    folder: "loads-and-supports",
    commandIds: ["SUPPORT", "CLOAD", "DLOAD", "PLOAD", "VLOAD", "TLOAD", "INERTIALOAD", "AMPLITUDE"]
  },
  {
    id: "dsl-constraints",
    title: "Constraints",
    folder: "constraints",
    commandIds: ["COUPLING", "TIE", "CONTACT", "CONNECTOR", "RBM"]
  },
  {
    id: "dsl-analysis",
    title: "Analysis",
    folder: "analysis",
    commandIds: ["LOADCASE", "SUPPORTS", "LOADS", "SOLVER", "CONSTRAINTMETHOD", "NONLINEAR", "INERTIARELIEF", "REBALANCELOADS", "REQUESTSTIFFNESS", "REQUESTSTGEOM", "CONSTRAINTSUMMARY", "NUMEIGENVALUES", "SIGMA", "TOPODENSITY", "TOPOORIENT", "TOPOEXPONENT", "TIME", "NEWMARK", "DAMPING", "WRITEEVERY", "INITIALVELOCITY", "OVERVIEW"]
  }
];

const navGroups = [
  {
    title: "Start",
    children: [
      { id: "overview", label: "Overview" },
      { id: "installation", label: "Installation guide" },
      { id: "quick-start", label: "Quick start" }
    ]
  },
  {
    title: "DSL",
    children: [
      {
        id: "dsl",
        label: "DSL",
        children: commandGroups.map((group) => ({
          id: group.id,
          label: group.title,
          children: group.commandIds.map((commandId) => ({
            commandId,
            label: `*${commandId}`
          }))
        }))
      }
    ]
  },
  {
    title: "Reference",
    children: [
      { id: "results", label: "Result format" },
      { id: "examples", label: "Examples" },
      { id: "implementation-notes", label: "Implementation notes" }
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

const commandCopy = {
  "NODE": "Use *NODE to create the coordinate rows that every later model definition refers to. A node is only a geometric point until elements, loads, supports, constraints, or features make degrees of freedom active on it.",
  "ELEMENT": "Use *ELEMENT to turn node connectivity into finite elements. The TYPE keyword chooses the topology and therefore controls how many node ids each dataline must provide.",
  "SURFACE": "Use *SURFACE to name element faces or line sides so later commands can apply pressure, traction, ties, contact, or couplings to a geometric boundary instead of to raw element ids.",
  "NSET": "Use *NSET to group nodes under a readable name. This keeps load, support, coupling, and feature definitions independent from long lists of node ids.",
  "ELSET": "Use *ELSET to group elements into modelling regions. Sections, volume loads, topology controls, and generated surfaces usually target element sets instead of individual element ids.",
  "SFSET": "Use *SFSET to group surfaces into reusable boundary regions. Pressure loads, distributed loads, ties, contact, and surface couplings can then refer to the surface set name.",
  "MATERIAL": "Use *MATERIAL to create or select a material context. The following material subcommands add elastic, density, and thermal data to that active material.",
  "ELASTIC": "Use *ELASTIC inside a material to define its stiffness law. The TYPE keyword selects the expected dataline shape, from simple isotropic constants to orthotropic or ABD data.",
  "DENSITY": "Use *DENSITY inside a material when mass matters. Modal, transient, inertia, and point-mass workflows depend on consistent density units.",
  "THERMALEXPANSION": "Use *THERMALEXPANSION inside a material when temperature changes should generate strain. Thermal loads only become meaningful when the referenced material has expansion data.",
  "SOLIDSECTION": "Use *SOLIDSECTION to assign a material to solid elements. Without a section, solid connectivity has geometry but no structural stiffness data.",
  "SHELLSECTION": "Use *SHELLSECTION to assign material and thickness to shell elements. Thickness controls membrane, bending, and shear stiffness for the selected shell set.",
  "BEAMSECTION": "Use *BEAMSECTION to assign material, profile, and local section orientation to beam elements. The optional dataline defines the local section direction used for bending axes.",
  "TRUSSSECTION": "Use *TRUSSSECTION to assign material and area to truss elements. Trusses carry axial response, so use beam elements when bending or torsion are required.",
  "PROFILE": "Use *PROFILE to define reusable beam cross-section constants. Beam sections reference this profile by name instead of repeating area and inertia values.",
  "FIELD": "Use *FIELD to attach structured numeric data to nodes, elements, or integration points. Fields drive topology density, orientation, initial velocity, temperature-like workflows, and other model data.",
  "ORIENTATION": "Use *ORIENTATION to define local coordinate systems. Materials, loads, supports, connectors, and section axes can then interpret components in local rather than global directions.",
  "POINTMASS": "Use *POINTMASS to add concentrated mass, rotary inertia, and optional spring stiffness to nodes. It is useful for lumped hardware, fixtures, payloads, and simplified dynamic models.",
  "SUPPORT": "Use *SUPPORT to collect prescribed nodal components into a named support collector. The collector is only active in an analysis when a loadcase references it with *SUPPORTS.",
  "CLOAD": "Use *CLOAD for concentrated nodal forces and moments. Each dataline targets a node or node set and contributes six mechanical components to a named load collector.",
  "DLOAD": "Use *DLOAD for vector traction on surfaces. The vector can be global or interpreted through an orientation, and all records accumulate in the selected load collector.",
  "PLOAD": "Use *PLOAD for scalar pressure on surfaces. The actual force direction follows the selected surface normal, so element side and shell side choices matter.",
  "VLOAD": "Use *VLOAD for body-force style loading over elements or element sets. Typical uses include gravity-like acceleration fields or distributed volumetric forcing.",
  "TLOAD": "Use *TLOAD to create thermal loading from a temperature field and a reference temperature. It needs compatible material thermal expansion data to produce strain loads.",
  "INERTIALOAD": "Use *INERTIALOAD to add acceleration and angular inertial effects as loads. This is a load definition, separate from the loadcase-level inertia relief option.",
  "AMPLITUDE": "Use *AMPLITUDE to define a reusable time scale curve. Transient loads can reference it so the same load collector changes magnitude over time.",
  "COUPLING": "Use *COUPLING to relate a master node set to slave nodes or a slave surface. Kinematic coupling enforces motion compatibility; structural coupling distributes generalized work.",
  "TIE": "Use *TIE to bind slave nodes or surfaces to master surfaces or lines through interpolation equations. Search distance, adjustment, and surface normals determine which slave points are tied.",
  "CONTACT": "Use *CONTACT for frictionless penalty contact in nonlinear static analysis. It contributes nonlinear force and tangent terms rather than ordinary linear constraint equations.",
  "CONNECTOR": "Use *CONNECTOR to impose idealized relative-motion constraints between two node sets. The connector type and coordinate system define which local movements are locked.",
  "RBM": "Use *RBM to add rigid-body suppression equations for an element set. It is mostly a stabilization and diagnostic tool for free or nearly free models.",
  "LOADCASE": "Use *LOADCASE to start an analysis block. Commands inside the block choose collectors, solver settings, constraint handling, and analysis-specific controls.",
  "SUPPORTS": "Use *SUPPORTS inside a loadcase to activate one or more support collectors. Only listed collectors participate in that analysis step.",
  "LOADS": "Use *LOADS inside a loadcase to activate one or more load collectors. Unlisted load collectors remain defined in the model but are not assembled for that step.",
  "SOLVER": "Use *SOLVER inside a loadcase to choose device and linear solution strategy. The available combinations depend on constraint handling and analysis type.",
  "CONSTRAINTMETHOD": "Use *CONSTRAINTMETHOD to choose how supported static loadcases enforce constraints. NULLSPACE reduces the system; LAGRANGE augments it with multipliers and has stricter solver limits.",
  "NONLINEAR": "Use *NONLINEAR to configure increment and Newton iteration controls for nonlinear static analysis. It controls convergence tolerance, iteration count, and regularization behaviour.",
  "INERTIARELIEF": "Use *INERTIARELIEF in static free-body workflows where supports are intentionally absent. It balances external loads by adding inertial effects and temporary stabilization.",
  "REBALANCELOADS": "Use *REBALANCELOADS to remove residual resultant force and moment from static loads. It is a numerical balancing tool, not a replacement for physical inertia modelling.",
  "REQUESTSTIFFNESS": "Use *REQUESTSTIFFNESS to write stiffness-related matrices for inspection or downstream tooling. It is diagnostic output and does not alter the solve itself.",
  "REQUESTSTGEOM": "Use *REQUESTSTGEOM in buckling analysis to write geometric stiffness data. This is mainly useful when validating preload and buckling operators.",
  "CONSTRAINTSUMMARY": "Use *CONSTRAINTSUMMARY to print the generated constraint equations for a loadcase. It helps diagnose supports, couplings, ties, connectors, and rigid-body equations.",
  "NUMEIGENVALUES": "Use *NUMEIGENVALUES to request how many eigenpairs a modal or buckling loadcase should compute. The solver may still limit the count to the available reduced system size.",
  "SIGMA": "Use *SIGMA to set the spectral shift for linear buckling. A useful shift steers the eigensolver toward the expected buckling factor range.",
  "TOPODENSITY": "Use *TOPODENSITY to select the element density field for topology static analysis. The field values scale element stiffness through the configured penalization exponent.",
  "TOPOORIENT": "Use *TOPOORIENT to select an element orientation field for topology static analysis. The field supplies orientation data used during stiffness assembly.",
  "TOPOEXPONENT": "Use *TOPOEXPONENT to set the density penalization exponent for topology static analysis. Larger values penalize intermediate densities more strongly.",
  "TIME": "Use *TIME to define the transient analysis interval and fixed time increment. The accepted dataline can include an explicit start time or just end time and step size.",
  "NEWMARK": "Use *NEWMARK to set the beta and gamma parameters for implicit transient integration. These parameters control stability, numerical damping, and integration accuracy.",
  "DAMPING": "Use *DAMPING to add proportional damping to transient analysis. Currently this is Rayleigh damping with mass- and stiffness-proportional coefficients.",
  "WRITEEVERY": "Use *WRITEEVERY to control transient result cadence. It can write every N steps or at an approximate physical time interval.",
  "INITIALVELOCITY": "Use *INITIALVELOCITY to initialize transient velocities from a node field. The field must match the active nodal component convention expected by the solver.",
  "OVERVIEW": "Use *OVERVIEW as an in-deck diagnostic checkpoint. It prints model counts and registered entities without changing the model state."
};

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


const relatedTopicMap = {
  nodes: "dsl-model",
  surfaces: "dsl-model",
  sets: "dsl-model",
  materials: "dsl-properties",
  sections: "dsl-properties",
  fields: "dsl-properties",
  "loads-supports": "dsl-loads",
  constraints: "dsl-constraints",
  "analysis-controls": "dsl-analysis",
  loadcases: "dsl-analysis",
  "linear-static": "dsl-analysis",
  "nonlinear-static": "dsl-analysis",
  "modal-buckling": "dsl-analysis",
  transient: "dsl-analysis",
  results: "results",
  examples: "examples",
  overview: "overview"
};

const commandGroupById = Object.fromEntries(
  commandGroups.flatMap((group) => group.commandIds.map((commandId) => [commandId, group]))
);

// ─────────────────────────────────────────────────────────────────────────────
// THEME TOGGLE
// ─────────────────────────────────────────────────────────────────────────────
(function initTheme() {
  const root = document.documentElement;
  const storageKey = 'femaster-docs-theme';
  const system = matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
  const initialTheme = root.getAttribute('data-theme');
  let current = initialTheme === 'light' || initialTheme === 'dark' ? initialTheme : system;

  try {
    const savedTheme = localStorage.getItem(storageKey);
    if (savedTheme === 'light' || savedTheme === 'dark') current = savedTheme;
  } catch (_) {
    // localStorage may be unavailable when the docs are opened via file://.
  }

  root.setAttribute('data-theme', current);

  function syncIcon(theme) {
    const icon = theme === 'dark'
      ? '<circle cx="12" cy="12" r="5"/><path d="M12 1v2M12 21v2M4.22 4.22l1.42 1.42M18.36 18.36l1.42 1.42M1 12h2M21 12h2M4.22 19.78l1.42-1.42M18.36 5.64l1.42-1.42"/>'
      : '<path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/>';
    const label = theme === 'dark' ? 'Switch to light mode' : 'Switch to dark mode';

    document.querySelectorAll('[data-theme-toggle]').forEach((button) => {
      button.setAttribute('aria-label', label);
      button.setAttribute('title', label);
      button.innerHTML = `<svg data-theme-icon width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" aria-hidden="true">${icon}</svg>`;
    });
  }
  syncIcon(current);

  document.querySelectorAll('[data-theme-toggle]').forEach(btn => {
    btn.addEventListener('click', () => {
      current = current === 'dark' ? 'light' : 'dark';
      root.setAttribute('data-theme', current);
      try { localStorage.setItem(storageKey, current); } catch (_) {}
      syncIcon(current);
    });
  });
})();

// ─────────────────────────────────────────────────────────────────────────────
// HELPERS
// ─────────────────────────────────────────────────────────────────────────────
function currentPageId() {
  return document.body.dataset.pageId || 'overview';
}

function currentRootPrefix() {
  const marker = '/pages/';
  const index = location.pathname.indexOf(marker);
  if (index === -1) return '';
  const tail = location.pathname.slice(index + marker.length);
  const depth = Math.max(1, tail.split('/').filter(Boolean).length);
  return '../'.repeat(depth);
}

function toPageHref(pageId) {
  const p = sitePages[pageId] || sitePages[relatedTopicMap[pageId]] || sitePages.overview;
  return `${currentRootPrefix()}${p.file}`;
}

function commandFileName(id) { return `${id.toLowerCase()}.html`; }

function commandHref(id) {
  const group = commandGroupById[id];
  if (!group) return toPageHref('dsl');
  return `${currentRootPrefix()}pages/dsl/${group.folder}/${commandFileName(id)}`;
}

function escapeHtml(v) {
  return String(v).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#039;'}[c]));
}

function inlineCommandLinks(value) {
  return escapeHtml(String(value)).replace(/\*([A-Z0-9_]+)/g, (match, id) => {
    const target = commands.find(cmd => cmd.id === id);
    return target ? `<a href="${commandHref(id)}"><code>${match}</code></a>` : `<code>${match}</code>`;
  });
}

// ─────────────────────────────────────────────────────────────────────────────
// SIDEBAR NAV
// ─────────────────────────────────────────────────────────────────────────────
const SVG_CHEVRON = `<svg class="tree-toggle-icon" viewBox="0 0 10 10" fill="none" stroke="currentColor" stroke-width="1.5" aria-hidden="true"><path d="M3 2l4 3-4 3"/></svg>`;
const SVG_CHEVRON_HEADING = `<svg class="tree-heading-chevron" viewBox="0 0 12 12" fill="none" stroke="currentColor" stroke-width="1.5" aria-hidden="true"><path d="M4 2l4 4-4 4"/></svg>`;

function buildSidebar() {
  const nav = document.getElementById('tree-nav');
  if (!nav) return;
  nav.innerHTML = navGroups.map((group, idx) => `
    <section class="tree-section" data-open="${idx < 2 ? 'true' : 'false'}">
      <button class="tree-heading" type="button">${escapeHtml(group.title)}${SVG_CHEVRON_HEADING}</button>
      <div class="tree-items">${renderTreeNodes(group.children || [], 0)}</div>
    </section>
  `).join('');

  nav.querySelectorAll('.tree-heading').forEach(btn => {
    btn.addEventListener('click', () => {
      const section = btn.closest('.tree-section');
      section.dataset.open = section.dataset.open === 'true' ? 'false' : 'true';
    });
  });
  nav.querySelectorAll('.tree-toggle').forEach(btn => {
    btn.addEventListener('click', () => {
      const node = btn.closest('.tree-node');
      node.dataset.open = node.dataset.open === 'true' ? 'false' : 'true';
    });
  });
}

function renderTreeNodes(items, level) {
  return items.map(item => {
    const hasChildren = Boolean(item.children?.length);
    const commandId = item.commandId || '';
    const pageId = item.id || '';
    const href = commandId ? commandHref(commandId) : toPageHref(pageId);
    const activeAttr = commandId ? `data-command-id="${commandId}"` : `data-page-id="${pageId}"`;
    const isActive = commandId
      ? commandId === (document.body.dataset.commandId || '')
      : pageId === currentPageId();
    const isOpen = hasChildren && (level === 0 || treeItemContainsActive(item));

    return `
      <div class="tree-node" data-open="${isOpen ? 'true' : 'false'}">
        <div class="tree-row" style="--level:${level}">
          ${hasChildren
            ? `<button class="tree-toggle" type="button" aria-label="Expand ${item.label}">${SVG_CHEVRON}</button>`
            : `<span class="tree-spacer"></span>`}
          <a class="tree-link${isActive ? ' active' : ''}" href="${href}" ${activeAttr}>${escapeHtml(item.label)}</a>
        </div>
        ${hasChildren ? `<div class="tree-children">${renderTreeNodes(item.children, level + 1)}</div>` : ''}
      </div>`;
  }).join('');
}

function treeItemContainsActive(item) {
  if (item.id && item.id === currentPageId()) return true;
  const bodyCmd = document.body.dataset.commandId || '';
  if (item.commandId && item.commandId === bodyCmd) return true;
  return Boolean(item.children?.some(treeItemContainsActive));
}

function setActiveNav() {
  const bodyCmd = document.body.dataset.commandId || '';
  document.querySelectorAll('.tree-link').forEach(link => {
    const active = (link.dataset.pageId && link.dataset.pageId === currentPageId() && !bodyCmd)
      || (link.dataset.commandId && link.dataset.commandId === bodyCmd);
    link.classList.toggle('active', active);
  });
}

// ─────────────────────────────────────────────────────────────────────────────
// BREADCRUMBS
// ─────────────────────────────────────────────────────────────────────────────
function updateBreadcrumbs(label) {
  const el = document.getElementById('breadcrumbs');
  if (!el) return;
  const sep = `<span class="breadcrumb-sep">/</span>`;
  if (label) { el.innerHTML = escapeHtml(label); return; }
  const bodyCmd = document.body.dataset.commandId || '';
  if (bodyCmd) {
    const cmd = commands.find(c => c.id === bodyCmd);
    const group = commandGroupById[bodyCmd];
    if (cmd && group) {
      el.innerHTML = `<a href="${toPageHref('dsl')}">DSL</a>${sep}<a href="${toPageHref(group.id)}">${escapeHtml(group.title)}</a>${sep}<span>*${escapeHtml(cmd.name)}</span>`;
      return;
    }
  }
  const p = sitePages[currentPageId()] || sitePages.overview;
  if (currentPageId().startsWith('dsl-') && currentPageId() !== 'dsl') {
    el.innerHTML = `<a href="${toPageHref('dsl')}">DSL</a>${sep}<span>${escapeHtml(p.title)}</span>`;
    return;
  }
  el.innerHTML = p.section !== p.title
    ? `<span>${escapeHtml(p.section)}</span>${sep}<span>${escapeHtml(p.title)}</span>`
    : `<span>${escapeHtml(p.title)}</span>`;
}

// ─────────────────────────────────────────────────────────────────────────────
// COMMAND CARD RENDERING
// ─────────────────────────────────────────────────────────────────────────────
function getCommandDoc(cmd) {
  const generated = (typeof commandDocs !== 'undefined' && commandDocs[cmd.id]) || {};
  const manual = (typeof commandManualDocs !== 'undefined' && commandManualDocs[cmd.id]) || {};
  const group = commandGroupById[cmd.id];
  const summary = manual.summary || commandCopy[cmd.id] || generated.description || cmd.purpose;
  return {
    id: cmd.id, name: cmd.name,
    summary, category: manual.category || group?.title || cmd.category,
    group, cmd, generated,
    admittedIn: admittedUnder(cmd),
    explanation: (typeof commandExplanations !== 'undefined' && commandExplanations[cmd.id]) || [summary],
    syntax: cmd.id === 'ELEMENT' ? [elementSyntax()] : manual.syntax || deriveSyntax(cmd, generated),
    examples: cmd.id === 'ELEMENT' ? [] : manual.examples || deriveExamples(cmd, generated),
    keywordsManual: manual.keywordsManual || null,
    keywords: generated.keywords || [],
    datalinesManual: manual.datalinesManual || null,
    variants: generated.variants || [],
    validation: manual.validation || deriveValidation(cmd, generated),
    related: manual.related || deriveRelated(cmd)
  };
}

function elementSyntax() {
  return {
    title: 'Basic form',
    code: '*ELEMENT, TYPE=<type>, ELSET=<optional-set>\n<element-id>, <node-1>, ..., <node-n>',
    description: ''
  };
}

function renderCommandCard(cmd) {
  const doc = getCommandDoc(cmd);
  return `
    <a class="command-card" href="${commandHref(cmd.id)}">
      <span class="command-name">*${cmd.name}</span>
      <div class="command-card-desc">${escapeHtml(doc.summary)}</div>
      <div class="command-card-context"><span>Valid in</span> ${doc.admittedIn.map(value => `<code>${escapeHtml(value)}</code>`).join(' ')}</div>
    </a>`;
}

function renderCommandIndex() {
  const activeCategory = state.commandCategory;
  const categories = ['All', ...commandCategories];
  const filtered = commands.filter(cmd => activeCategory === 'All' || cmd.category === activeCategory);
  return `
    <div class="command-toolbar">
      ${categories.map(cat => `<button class="filter-button${cat === activeCategory ? ' active' : ''}" type="button" data-category="${cat}">${cat}</button>`).join('')}
    </div>
    <div class="command-grid">${filtered.map(renderCommandCard).join('')}</div>`;
}

function renderCommandGroupGrid(groupId) {
  const group = commandGroups.find(g => g.id === groupId);
  if (!group) return '';
  const cmds = group.commandIds.map(id => commands.find(c => c.id === id)).filter(Boolean);
  return `<div class="command-grid">${cmds.map(renderCommandCard).join('')}</div>`;
}

// ─────────────────────────────────────────────────────────────────────────────
// COMMAND DETAIL PAGE
// ─────────────────────────────────────────────────────────────────────────────
function renderCommandDetail(cmd) {
  const doc = getCommandDoc(cmd);
  return `
    <article class="command-page">
      ${renderCommandHero(doc)}
      ${renderCommandExplanation(doc)}
      ${renderCommandOverview(doc)}
      ${renderExamples(doc)}
      ${renderKeywordReference(doc)}
      ${renderDatalineReference(doc)}
      ${renderListBlock('Validation rules', doc.validation, 'validation-list')}
      ${cmd.id === 'ELEMENT' ? renderElementLibraryBlock() : ''}
      ${renderRelatedCommands(doc)}
    </article>`;
}

function renderCommandExplanation(doc) {
  return `<section class="command-explanation">
    <h2>Overview</h2>
    ${doc.explanation.map(paragraph => `<p>${inlineCommandLinks(paragraph)}</p>`).join('')}
  </section>`;
}

function renderCommandHero(doc) {
  return `
    <header class="command-header">
      <nav class="page-breadcrumbs" aria-label="Breadcrumb">
        <a href="${toPageHref('dsl')}">DSL</a>
        <span class="meta-separator">/</span>
        ${doc.group ? `<a href="${toPageHref(doc.group.id)}">${escapeHtml(doc.group.title)}</a><span class="meta-separator">/</span>` : ''}
        <span>*${escapeHtml(doc.name)}</span>
      </nav>
      <div class="command-titlebar">
        <div>
          <div class="eyebrow">DSL command</div>
          <h1 class="command-title"><span class="command-prefix">*</span>${escapeHtml(doc.name)}</h1>
          <p class="command-desc">${escapeHtml(doc.summary)}</p>
        </div>
        <div class="command-context"><span>Valid in</span>${doc.admittedIn.map(value => `<code>${escapeHtml(value)}</code>`).join('')}</div>
      </div>
    </header>`;
}

function renderCommandOverview(doc) {
  const syntaxHtml = doc.syntax.map(item => `
    <div class="syntax-variant">
      ${item.title ? `<h3>${escapeHtml(item.title)}</h3>` : ''}
      <pre><code>${escapeHtml(item.code || '')}</code></pre>
      ${item.description ? `<p>${escapeHtml(item.description)}</p>` : ''}
    </div>`).join('');

  return `<section class="syntax-card">
    <div class="syntax-card-header">Syntax</div>
    <div class="syntax-card-body">${syntaxHtml}</div>
  </section>`;
}

function renderExamples(doc) {
  if (!doc.examples.length) return '';
  return `
    <section class="doc-block">
      <div class="doc-block-header"><h2>Minimal example</h2></div>
      ${doc.examples.map(ex => `
        <div class="example-card">
          ${ex.title ? `<h3>${escapeHtml(ex.title)}</h3>` : ''}
          <pre><code>${escapeHtml(ex.code)}</code></pre>
          ${ex.description ? `<p>${escapeHtml(ex.description)}</p>` : ''}
        </div>`).join('')}
    </section>`;
}

function renderKeywordReference(doc) {
  const manual = doc.keywordsManual;
  if (manual) {
    return `<section class="doc-block keyword-reference">
      <div class="doc-block-header"><h2>Keyword reference</h2></div>
      ${renderKeywordTable('Required keywords', manual.required || [])}
      ${renderKeywordTable('Optional keywords', manual.optional || [])}
      ${(manual.groups || []).map(renderKeywordGroup).join('')}
    </section>`;
  }
  const required = doc.keywords.filter(kw => kw.status === 'required');
  const optional = doc.keywords.filter(kw => kw.status !== 'required');
  if (!required.length && !optional.length) return renderEmptyBlock('Keyword reference', 'This command has no command-line keywords.');
  return `<section class="doc-block keyword-reference">
    <div class="doc-block-header"><h2>Keyword reference</h2></div>
    ${renderKeywordTable('Required keywords', required.map(keywordToRow))}
    ${renderKeywordTable('Optional keywords', optional.map(keywordToRow))}
  </section>`;
}

function renderKeywordTable(title, rows) {
  if (!rows.length) return '';
  return `<section class="reference-subsection"><h3>${escapeHtml(title)}</h3>
    <div class="table-wrap"><table class="reference-table">
      <thead><tr><th>Keyword</th><th>Allowed values</th><th>Default</th><th>Description</th></tr></thead>
      <tbody>${rows.map(row => `<tr>
        <td><code>${escapeHtml(row.name)}</code></td>
        <td>${formatCellCode(row.values || row.allowed || '—')}</td>
        <td>${formatCellCode(row.default || '—')}</td>
        <td>${escapeHtml(row.description || '')}</td>
      </tr>`).join('')}</tbody>
    </table></div>
  </section>`;
}

function renderKeywordGroup(group) {
  return `<section class="reference-subsection"><h3>${escapeHtml(group.title)}</h3>
    <div class="table-wrap"><table class="reference-table">
      <thead><tr><th>Keyword</th><th>Relation</th><th>Description</th></tr></thead>
      <tbody>${(group.rows || []).map(row => `<tr>
        <td><code>${escapeHtml(row.name)}</code></td>
        <td>${escapeHtml(row.relation)}</td>
        <td>${escapeHtml(row.description)}</td>
      </tr>`).join('')}</tbody>
    </table></div>
  </section>`;
}

function renderDatalineReference(doc) {
  const datalines = doc.id === 'ELEMENT'
    ? [elementConnectivityDataline()]
    : doc.datalinesManual || generatedDatalines(doc.variants);
  if (!datalines.length) return renderEmptyBlock('Dataline reference', 'This command does not require datalines.');
  return `<section class="doc-block">
    <div class="doc-block-header"><h2>Dataline reference</h2></div>
    <div class="dataline-list">${datalines.map(renderDatalineCard).join('')}</div>
  </section>`;
}

function elementConnectivityDataline() {
  return {
    title: 'Element connectivity',
    signature: '<element-id>, <node-1>, ..., <node-n>',
    rows: 'one or more',
    description: 'The number and order of node ids are defined by TYPE.',
    fields: [
      { name: 'element-id', meaning: 'Unique element id', type: 'int' },
      { name: 'node-1 ... node-n', meaning: 'Ordered element connectivity', type: 'int' }
    ]
  };
}

function renderDatalineCard(line) {
  return `<article class="dataline-card">
    <div class="dataline-card-header">
      <h3>${escapeHtml(line.title || 'Dataline')}</h3>
      ${line.rows ? `<span>${escapeHtml(line.rows)}</span>` : ''}
    </div>
    ${line.signature ? `<pre><code>${escapeHtml(line.signature)}</code></pre>` : ''}
    ${line.description ? `<p>${escapeHtml(line.description)}</p>` : ''}
    ${(line.fields || []).length ? `
      <div class="table-wrap"><table class="reference-table">
        <thead><tr><th>Field</th><th>Meaning</th><th>Type</th></tr></thead>
        <tbody>${line.fields.map(f => `<tr>
          <td><code>${escapeHtml(f.name)}</code></td>
          <td>${escapeHtml(f.meaning || '')}</td>
          <td>${escapeHtml(f.type || '')}</td>
        </tr>`).join('')}</tbody>
      </table></div>` : ''}
  </article>`;
}

function renderListBlock(title, items, className) {
  if (!items?.length) return '';
  return `<section class="doc-block">
    <div class="doc-block-header"><h2>${escapeHtml(title)}</h2></div>
    <ul class="${className}">${items.map(i => `<li>${escapeHtml(i)}</li>`).join('')}</ul>
  </section>`;
}

function renderRelatedCommands(doc) {
  if (!doc.related.length) return '';
  return `<section class="doc-block">
    <div class="doc-block-header"><h2>Related commands</h2></div>
    <div class="related-command-grid">${doc.related.map(item => {
      const href = commands.find(c => c.id === item.command) ? commandHref(item.command) : '#';
      return `<a class="related-command-card" href="${href}"><code>*${escapeHtml(item.command)}</code><span class="muted">${escapeHtml(item.description || '')}</span></a>`;
    }).join('')}</div>
  </section>`;
}

function renderElementLibraryBlock() {
  return `<section class="doc-block"><div class="doc-block-header"><h2>Element library</h2></div>${renderElements()}</section>`;
}

function renderEmptyBlock(title, text) {
  return `<section class="doc-block">
    <div class="doc-block-header"><h2>${escapeHtml(title)}</h2></div>
    <div class="empty-state">${escapeHtml(text)}</div>
  </section>`;
}

function renderElements() {
  const figureBase = `${currentRootPrefix()}../pages/figures/element_schematics/`;
  return `<div class="element-grid">${elements.map(([name, title, text, image]) => `
    <article class="element-card">
      ${image ? `<img src="${figureBase}${image}" alt="${name} element schematic" loading="lazy">`
               : `<div class="empty-state">No schematic</div>`}
      <div class="element-card-body">
        <code class="element-type">${name}</code>
        <h3>${escapeHtml(title)}</h3>
        <p class="muted">${escapeHtml(text)}</p>
      </div>
    </article>`).join('')}</div>`;
}

// ─────────────────────────────────────────────────────────────────────────────
// SEARCH
// ─────────────────────────────────────────────────────────────────────────────
function renderSearch(query) {
  const q = query.trim().toLowerCase();
  const commandHits = commands.filter(cmd =>
    `${cmd.name} ${cmd.category} ${cmd.purpose} ${cmd.details.join(' ')}`.toLowerCase().includes(q));
  const pageHits = Object.entries(sitePages).filter(([id, p]) =>
    `${id} ${p.title} ${p.section}`.toLowerCase().includes(q));

  return `
    <div class="page-hero">
      <div class="eyebrow">Search</div>
      <h1>Results for <code>${escapeHtml(query)}</code></h1>
    </div>
    <h2>Commands</h2>
    <div class="command-grid">${commandHits.length ? commandHits.map(renderCommandCard).join('') : `<div class="empty-state">No command matches.</div>`}</div>
    <h2>Topics</h2>
    <div class="grid two">${pageHits.length ? pageHits.map(([id, p]) => `<a class="card" href="${toPageHref(id)}"><h3>${escapeHtml(p.title)}</h3><p>${escapeHtml(p.section)}</p></a>`).join('') : `<div class="empty-state">No topic matches.</div>`}</div>`;
}

// ─────────────────────────────────────────────────────────────────────────────
// DERIVE HELPERS (from generated data)
// ─────────────────────────────────────────────────────────────────────────────
function deriveSyntax(cmd, generated) {
  if (cmd.id === 'ELEMENT') {
    return [{ title: 'Basic form', code: cmd.syntax, description: '' }];
  }
  const req = (generated.keywords || []).filter(kw => kw.status === 'required').map(kw => `${kw.name}=<${kw.name.toLowerCase()}>`);
  const opt = (generated.keywords || []).filter(kw => kw.status !== 'required').map(kw => `${kw.name}=<optional>`);
  const first = generatedDatalines(generated.variants || [])[0];
  return [{ title: 'Basic form',
    code: `*${cmd.name}${req.length ? `, ${req.join(', ')}` : ''}${opt.length ? `, ${opt.join(', ')}` : ''}${first?.signature ? `\n${first.signature}` : ''}`,
    description: 'Generated from the command keywords and dataline layout.' }];
}

function deriveExamples(cmd, generated) {
  const code = deriveSyntax(cmd, generated)[0]?.code;
  return code ? [{ title: 'Minimal form', code, description: 'Replace placeholder values with model-specific names and numbers.' }] : [];
}

function deriveValidation(cmd, generated) {
  return (generated.keywords || []).filter(kw => kw.status === 'required').map(kw => `${kw.name} must be provided.`);
}

function deriveRelated(cmd) {
  return (cmd.links || []).map(id => {
    const pageId = relatedTopicMap[id] || id;
    const target = commands.find(c => c.id === pageId.toUpperCase());
    return target ? { command: target.id, description: 'Related DSL command.' } : null;
  }).filter(Boolean);
}

function keywordToRow(keyword) {
  return {
    name: keyword.name,
    values: keyword.allowed?.length ? keyword.allowed.join(', ') : 'valid value',
    default: keyword.default || '—',
    description: keyword.description || (keyword.status === 'required' ? 'Required.' : 'Optional.')
  };
}

function generatedDatalines(variants) {
  const out = [];
  for (const variant of variants || []) {
    for (const layout of variant.layouts || []) {
      out.push({
        title: variantTitle(variant),
        signature: layout.fields,
        rows: formatLineRange(layout.range),
        description: variantDescription(variant),
        fields: (layout.entries || []).map(e => ({
          name: e.name,
          meaning: stripType(e.description),
          type: extractType(e.description)
        }))
      });
    }
  }
  return out;
}

function stripType(v) { return (v || '').replace(/\s*\[[^\]]+\]\s*$/, ''); }
function extractType(v) { const m = (v || '').match(/\[([^\]]+)\]\s*$/); return m ? m[1] : ''; }

function formatCellCode(value) {
  if (!value || value === '—') return '—';
  return String(value).split(/,\s*/).map(p => p.trim()).filter(Boolean).map(p => `<code>${escapeHtml(p)}</code>`).join(' ');
}

function variantTitle(variant) {
  const c = variantCondition(variant.when);
  if (c === 'always') return 'Default dataline';
  if (c.includes('TYPE=')) return c.replace('TYPE=', 'TYPE ');
  return `Dataline form ${variant.number}`;
}

function variantCondition(v) {
  return formatWhen(v).replace(/([A-Z0-9_]+) in \{ ([^}]*) \}/g, '$1=$2').replaceAll(', ', ' | ');
}

function variantDescription(variant) {
  const c = variantCondition(variant.when);
  return c === 'always' ? 'Default dataline format.' : `Use when ${c}.`;
}

function formatLineRange(r) {
  return String(r || '').replace('0..∞','zero or more').replace('1..∞','one or more').replace('1..1','exactly one').replace('0..1','optional');
}
function formatLayoutKind(k) { return k === 'single-line' ? 'repeatable row' : 'multi-row block'; }
function formatWhen(v) {
  if (!v || v === '(none)') return 'always';
  return v.replaceAll('self.keys["','').replaceAll('"]','');
}

function admittedUnder(cmd) {
  const exact = {
    ELASTIC:['*MATERIAL'],DENSITY:['*MATERIAL'],THERMALEXPANSION:['*MATERIAL'],
    SUPPORTS:['*LOADCASE'],LOADS:['*LOADCASE'],SOLVER:['*LOADCASE'],
    CONSTRAINTMETHOD:['*LOADCASE, TYPE=LINEARSTATIC','*LOADCASE, TYPE=NONLINEARSTATIC'],
    NONLINEAR:['*LOADCASE, TYPE=NONLINEARSTATIC'],
    INERTIARELIEF:['*LOADCASE, TYPE=LINEARSTATIC'],
    REBALANCELOADS:['*LOADCASE, TYPE=LINEARSTATIC'],
    REQUESTSTIFFNESS:['*LOADCASE'],REQUESTSTGEOM:['*LOADCASE, TYPE=LINEARBUCKLING'],
    CONSTRAINTSUMMARY:['*LOADCASE'],
    NUMEIGENVALUES:['*LOADCASE, TYPE=EIGENFREQ','*LOADCASE, TYPE=LINEARBUCKLING'],
    SIGMA:['*LOADCASE, TYPE=LINEARBUCKLING'],
    TOPODENSITY:['*LOADCASE, TYPE=LINEARSTATICTOPO'],
    TOPOORIENT:['*LOADCASE, TYPE=LINEARSTATICTOPO'],
    TOPOEXPONENT:['*LOADCASE, TYPE=LINEARSTATICTOPO'],
    TIME:['*LOADCASE, TYPE=LINEARTRANSIENT'],NEWMARK:['*LOADCASE, TYPE=LINEARTRANSIENT'],
    DAMPING:['*LOADCASE, TYPE=LINEARTRANSIENT'],WRITEEVERY:['*LOADCASE, TYPE=LINEARTRANSIENT'],
    INITIALVELOCITY:['*LOADCASE, TYPE=LINEARTRANSIENT']
  };
  return exact[cmd.id] || ['Root level'];
}

// ─────────────────────────────────────────────────────────────────────────────
// STATE & RENDER ORCHESTRATION
// ─────────────────────────────────────────────────────────────────────────────
const state = { commandCategory: 'All', search: '' };
let originalPageHtml = null;

function hydrateCommandDetail() {
  const node = document.querySelector('[data-command-detail]');
  if (!node) return false;
  const id = node.dataset.commandDetail.toUpperCase();
  const cmd = commands.find(c => c.id === id);
  const pageEl = node.closest('#page');
  if (pageEl) pageEl.innerHTML = cmd ? renderCommandDetail(cmd) : `<div class="empty-state">Unknown command: ${escapeHtml(id)}</div>`;
  return true;
}

function hydrateDynamicBlocks() {
  if (hydrateCommandDetail()) return;
  document.querySelectorAll('[data-command-index]').forEach(node => { node.innerHTML = renderCommandIndex(); });
  document.querySelectorAll('[data-command-group]').forEach(node => { node.innerHTML = renderCommandGroupGrid(node.dataset.commandGroup); });
  document.querySelectorAll('[data-elements]').forEach(node => { node.innerHTML = renderElements(); });
  document.querySelectorAll('[data-command-list]').forEach(node => {
    const ids = (node.dataset.commandList || '').split(',').map(id => id.trim().toUpperCase());
    node.innerHTML = `<div class="command-grid">${ids.map(id => { const cmd = commands.find(c => c.id === id); return cmd ? renderCommandCard(cmd) : ''; }).join('')}</div>`;
  });
}

function render() {
  const pageEl = document.getElementById('page');
  if (!pageEl) return;
  if (originalPageHtml === null) originalPageHtml = pageEl.innerHTML;

  if (state.search.trim()) {
    pageEl.innerHTML = renderSearch(state.search);
    updateBreadcrumbs('Search');
    setActiveNav();
    return;
  }

  if (pageEl.innerHTML !== originalPageHtml) pageEl.innerHTML = originalPageHtml;
  hydrateDynamicBlocks();
  updateBreadcrumbs();
  setActiveNav();
}

function initEvents() {
  const search = document.getElementById('global-search');
  if (search) {
    search.addEventListener('input', e => { state.search = e.target.value; render(); });
  }

  const navToggle = document.getElementById('nav-toggle');
  if (navToggle) {
    navToggle.addEventListener('click', () => document.body.classList.toggle('nav-open'));
  }

  document.addEventListener('click', e => {
    if (e.target.closest('.tree-link')) document.body.classList.remove('nav-open');
    if (document.body.classList.contains('nav-open') && !e.target.closest('.sidebar') && !e.target.closest('#nav-toggle')) {
      document.body.classList.remove('nav-open');
    }
    const filter = e.target.closest('[data-category]');
    if (filter) { state.commandCategory = filter.dataset.category; hydrateDynamicBlocks(); }
  });
}

// ─────────────────────────────────────────────────────────────────────────────
// INIT
// ─────────────────────────────────────────────────────────────────────────────
buildSidebar();
initEvents();
render();
