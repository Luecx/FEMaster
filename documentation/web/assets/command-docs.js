const commandDocs = {"AMPLITUDE":{"description":"Define a reusable scalar time history. Each data line specifies a time-value pair. The TYPE keyword controls interpolation (STEP, NEAREST, LINEAR). Loads referencing the amplitude will scale their components by the interpolated value at the current analysis time.","keywords":[{"name":"NAME","status":"required","default":"","allowed":[],"aliases":[],"description":"Amplitude identifier"},{"name":"TYPE","status":"optional","default":"LINEAR","allowed":["STEP","NEAREST","LINEAR"],"aliases":[],"description":"Interpolation scheme: STEP, NEAREST, or LINEAR (default)"}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..∞","kind":"single-line","fields":"TIME, VALUE","entries":[{"name":"TIME","description":"Time coordinate [double]"},{"name":"VALUE","description":"Amplitude value at TIME [double]"}]}]}]},"BEAMSECTION":{"description":"Assign a beam section with a local orientation. The target element set must contain only beam elements. The provided direction n1 is normalized and orthogonalized internally; if n1 is (near-)collinear with the element axis, a stable fallback is applied.","keywords":[{"name":"ELSET","status":"required","default":"","allowed":[],"aliases":[],"description":"Target element set; must contain only beam elements."},{"name":"MATERIAL","status":"required","default":"","allowed":[],"aliases":["MAT"],"description":"Material name (must exist)."},{"name":"PROFILE","status":"required","default":"","allowed":[],"aliases":[],"description":"Section/profile identifier, e.g., IPE80, RECT_20x10, CHS60x3."}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"Optional data line with the section direction n1 = (N1_x, N1_y, N1_z).","layouts":[{"range":"0..1","kind":"single-line","fields":"N11, N12, N13","entries":[{"name":"N11 - N13","description":"Section orientation vector components: N1_x, N1_y, N1_z. [double]"}]}]}]},"CLOAD":{"description":"Apply concentrated nodal forces and moments. Each data line supplies six components (Fx,Fy,Fz,Mx,My,Mz) for a target node or node set. Values accumulate when multiple records address the same entity. Optional ORIENTATION lets you specify the components in a local coordinate system, while AMPLITUDE references a time history used to scale the load during transient analyses.","keywords":[{"name":"AMPLITUDE","status":"optional","default":"","allowed":[],"aliases":[],"description":"Optional time amplitude that scales the load components."},{"name":"LOAD_COLLECTOR","status":"required","default":"","allowed":[],"aliases":[],"description":"Name of the active load collector that receives the contributions."},{"name":"ORIENTATION","status":"optional","default":"","allowed":[],"aliases":[],"description":"Optional coordinate system used to interpret the six load components."}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"Each data line defines a target node set/id followed by Fx, Fy, Fz, Mx, My, Mz (local or global depending on ORIENTATION).","layouts":[{"range":"1..∞","kind":"single-line","fields":"TARGET, LOAD1, LOAD2, LOAD3, LOAD4, LOAD5, LOAD6","entries":[{"name":"TARGET","description":"Node set name or single node id (integer). [string]"},{"name":"LOAD1 - LOAD6","description":"Fx, Fy, Fz, Mx, My, Mz. Uses the ORIENTATION basis when provided, otherwise global axes. [double]"}]}]}]},"CONNECTOR":{"description":"Define a connector between two node sets in a given coordinate system. Both node sets must exist and be non-empty. Their composition should fit the connector’s semantics. The coordinate system defines the local axes.","keywords":[{"name":"COORDINATESYSTEM","status":"required","default":"","allowed":[],"aliases":[],"description":"Name of a coordinate system."},{"name":"NSET1","status":"required","default":"","allowed":[],"aliases":[],"description":"First node set."},{"name":"NSET2","status":"required","default":"","allowed":[],"aliases":[],"description":"Second node set."},{"name":"TYPE","status":"required","default":"","allowed":[],"aliases":[],"description":"Connector type identifier."}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"CONSTRAINTMETHOD":{"description":"Select constraint backend for LINEARSTATIC/LINEARSTATICTOPO: NULLSPACE or LAGRANGE. Constraint | Backend | DIRECT | INDIRECT NULLSPACE | CPU MKL | Yes | Yes NULLSPACE | CPU Eigen | Yes | Yes NULLSPACE | GPU | Yes | Yes NULLSPACE | GPU cuDSS | Yes | Yes LAGRANGE | CPU MKL | Yes | No LAGRANGE | CPU Eigen | Limited | No LAGRANGE | GPU | No | No LAGRANGE | GPU cuDSS | Yes | No","keywords":[{"name":"TYPE","status":"required","default":"","allowed":["NULLSPACE","LAGRANGE"],"aliases":[],"description":""}],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"CONSTRAINTSUMMARY":{"description":"Enable constraint summary output for the active loadcase.","keywords":[],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"CONTACT":{"description":"Define frictionless node-to-surface penalty contact. MASTER must be a surface set; SLAVE may be a node set or surface set. Contact contributes only in NONLINEARSTATIC.","keywords":[{"name":"CLEARANCE","status":"optional","default":"0","allowed":[],"aliases":[],"description":"Allowed normal clearance before contact activates"},{"name":"DISTANCE","status":"required","default":"","allowed":[],"aliases":[],"description":"Search distance for candidate master surfaces"},{"name":"FLIP","status":"optional","default":"NO","allowed":["NO","YES"],"aliases":[],"description":"Flip master surface normals"},{"name":"MASTER","status":"required","default":"","allowed":[],"aliases":[],"description":"Master surface set"},{"name":"PENALTY","status":"required","default":"","allowed":[],"aliases":[],"description":"Normal penalty stiffness"},{"name":"SLAVE","status":"required","default":"","allowed":[],"aliases":[],"description":"Slave node set or slave surface set"}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"COUPLING":{"description":"Create kinematic or structural couplings between a MASTER and a SLAVE (node set) or SFSET (surface). Provide exactly one of SLAVE or SFSET. DOF entries > 0 enable coupling for [Ux, Uy, Uz, Rx, Ry, Rz]. KINEMATIC makes slaves follow the master; STRUCTURAL distributes/collects loads in a work-equivalent way. Avoid near-collinear axis/geometry definitions as they can cause unrealistic results.","keywords":[{"name":"MASTER","status":"required","default":"","allowed":[],"aliases":[],"description":"Master node set (commonly exactly one node)."},{"name":"SFSET","status":"optional","default":"","allowed":[],"aliases":[],"description":"Slave surface set."},{"name":"SLAVE","status":"optional","default":"","allowed":[],"aliases":[],"description":"Slave node set."},{"name":"TYPE","status":"required","default":"","allowed":["KINEMATIC","STRUCTURAL"],"aliases":[],"description":""}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"One data line: six DOF flags (1/0) for [Ux, Uy, Uz, Rx, Ry, Rz]. Values > 0 enable coupling.","layouts":[{"range":"1..1","kind":"single-line","fields":"DOF1, DOF2, DOF3, DOF4, DOF5, DOF6","entries":[{"name":"DOF1 - DOF6","description":"Coupled degrees of freedom (1=on, 0=off). [double]"}]}]}]},"DAMPING":{"description":"Assign damping for the active loadcase. Currently supported: Rayleigh proportional damping (C = α M + β K).","keywords":[{"name":"TYPE","status":"required","default":"","allowed":["RAYLEIGH"],"aliases":[],"description":"Damping model type. Only RAYLEIGH is supported."}],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"One data line: α, β coefficients for Rayleigh damping.","layouts":[{"range":"1..1","kind":"single-line","fields":"RAYLEIGH1, RAYLEIGH2","entries":[{"name":"RAYLEIGH1 - RAYLEIGH2","description":"α, β (mass- and stiffness-proportional coefficients). [double]"}]}]}]},"DENSITY":{"description":"Assign density to the active material.","keywords":[],"admittedRaw":["parent.command in { MATERIAL }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..1","kind":"single-line","fields":"RHO","entries":[{"name":"RHO","description":"Density value [double]"}]}]}]},"DLOAD":{"description":"Apply distributed surface tractions to surfaces or surface sets. Each line provides a target followed by three components (Fx,Fy,Fz). When ORIENTATION is specified the components are measured in that local frame; otherwise they are interpreted in the global Cartesian basis. AMPLITUDE references a time history that scales the traction. Contributions from multiple records add linearly in the active load collector.","keywords":[{"name":"AMPLITUDE","status":"optional","default":"","allowed":[],"aliases":[],"description":"Optional time amplitude that scales the traction vector"},{"name":"LOAD_COLLECTOR","status":"required","default":"","allowed":[],"aliases":[],"description":"Name of the load collector that aggregates the distributed loads"},{"name":"ORIENTATION","status":"optional","default":"","allowed":[],"aliases":[],"description":"Optional coordinate system describing the local traction directions"}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..∞","kind":"single-line","fields":"TARGET, LOAD1, LOAD2, LOAD3","entries":[{"name":"TARGET","description":"Surface set or surface id [string]"},{"name":"LOAD1 - LOAD3","description":"Fx,Fy,Fz components, local if ORIENTATION is set [double]"}]}]}]},"ELASTIC":{"description":"Assign elastic properties to the active material.","keywords":[{"name":"TYPE","status":"required","default":"","allowed":[],"aliases":[],"description":"GENISO, ORTHO, ORTHOTROPIC, ABD} Elasticity formulation (ISO/GENERALISED ISOTROPIC/ORTHO/ABD)"}],"admittedRaw":["parent.command in { MATERIAL }"],"variants":[{"number":1,"when":"self.keys[\"TYPE\"] in { ISO, ISOTROPIC }","description":"","layouts":[{"range":"1..1","kind":"single-line","fields":"E, NU","entries":[{"name":"E","description":"Young's modulus [double]"},{"name":"NU","description":"Poisson ratio [double]"}]}]},{"number":2,"when":"self.keys[\"TYPE\"] in { GENERALISEDISOTROPIC, GENERALISED_ISOTROPIC,","description":"","layouts":[{"range":"1..1","kind":"single-line","fields":"E, NU, G","entries":[{"name":"E","description":"Young's modulus [double]"},{"name":"NU","description":"Poisson ratio [double]"},{"name":"G","description":"Independent shear modulus [double]"}]}]},{"number":3,"when":"self.keys[\"TYPE\"] in { ORTHO, ORTHOTROPIC }","description":"","layouts":[{"range":"1..1","kind":"single-line","fields":"E1, E2, E3, G23, G13, G12, NU23, NU13, NU12","entries":[{"name":"E1","description":"Young's modulus in 1-direction [double]"},{"name":"E2","description":"Young's modulus in 2-direction [double]"},{"name":"E3","description":"Young's modulus in 3-direction [double]"},{"name":"G23","description":"Shear modulus G23 [double]"},{"name":"G13","description":"Shear modulus G13 [double]"},{"name":"G12","description":"Shear modulus G12 [double]"},{"name":"NU23","description":"Poisson ratio nu23 [double]"},{"name":"NU13","description":"Poisson ratio nu13 [double]"},{"name":"NU12","description":"Poisson ratio nu12 [double]"}]}]},{"number":4,"when":"self.keys[\"TYPE\"] in { ABD }","description":"","layouts":[{"range":"1..8","kind":"multiline","fields":"DATA1, DATA2, DATA3, DATA4, DATA5, DATA6, DATA7, DATA8, DATA9, DATA10, DATA11, DATA12, DATA13, DATA14, DATA15, DATA16, DATA17, DATA18, DATA19, DATA20, DATA21, DATA22, DATA23, DATA24, DATA25, DATA26, DATA27, DATA28, DATA29, DATA30, DATA31, DATA32, DATA33, DATA34, DATA35, DATA36, DATA37, DATA38, DATA39, DATA40","entries":[{"name":"DATA1 - DATA40","description":"ABD row-major (36 values), then shear row-major (4 values) [double]"}]}]}]},"ELEMENT":{"description":"Define finite elements by id and connectivity.","keywords":[{"name":"ELSET","status":"optional","default":"EALL","allowed":[],"aliases":[],"description":"Element set to activate while reading"},{"name":"TYPE","status":"required","default":"","allowed":[],"aliases":[],"description":"S3, S4, MITC4, MITC4FRT, S6, S8, QSPT} Element topology/type"}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"self.keys[\"TYPE\"] in { C3D4 }","description":"","layouts":[{"range":"1..∞","kind":"multiline","fields":"ID, N1, N2, N3, N4","entries":[{"name":"ID","description":"Element id [int]"},{"name":"N1 - N4","description":"C3D4 connectivity (4 nodes) [int]"}]}]},{"number":2,"when":"self.keys[\"TYPE\"] in { C3D6 }","description":"","layouts":[{"range":"1..∞","kind":"multiline","fields":"ID, N1, N2, N3, N4, N5, N6","entries":[{"name":"ID","description":"Element id [int]"},{"name":"N1 - N6","description":"C3D6 connectivity (6 nodes) [int]"}]}]},{"number":3,"when":"self.keys[\"TYPE\"] in { C3D8 }","description":"","layouts":[{"range":"1..∞","kind":"multiline","fields":"ID, N1, N2, N3, N4, N5, N6, N7, N8","entries":[{"name":"ID","description":"Element id [int]"},{"name":"N1 - N8","description":"C3D8 connectivity (8 nodes) [int]"}]}]},{"number":4,"when":"self.keys[\"TYPE\"] in { C3D10 }","description":"","layouts":[{"range":"1..∞","kind":"multiline","fields":"ID, N1, N2, N3, N4, N5, N6, N7, N8, N9, N10","entries":[{"name":"ID","description":"Element id [int]"},{"name":"N1 - N10","description":"C3D10 connectivity (10 nodes) [int]"}]}]},{"number":5,"when":"self.keys[\"TYPE\"] in { C3D15 }","description":"","layouts":[{"range":"1..∞","kind":"multiline","fields":"ID, N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12, N13, N14, N15","entries":[{"name":"ID","description":"Element id [int]"},{"name":"N1 - N15","description":"C3D15 connectivity (15 nodes) [int]"}]}]},{"number":6,"when":"self.keys[\"TYPE\"] in { C3D20 }","description":"","layouts":[{"range":"1..∞","kind":"multiline","fields":"ID, N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12, N13, N14, N15, N16, N17, N18, N19, N20","entries":[{"name":"ID","description":"Element id [int]"},{"name":"N1 - N20","description":"C3D20 connectivity (20 nodes) [int]"}]}]},{"number":7,"when":"self.keys[\"TYPE\"] in { C3D20R }","description":"","layouts":[{"range":"1..∞","kind":"multiline","fields":"ID, N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12, N13, N14, N15, N16, N17, N18, N19, N20","entries":[{"name":"ID","description":"Element id [int]"},{"name":"N1 - N20","description":"C3D20R connectivity (20 nodes) [int]"}]}]},{"number":8,"when":"self.keys[\"TYPE\"] in { T3 }","description":"","layouts":[{"range":"1..∞","kind":"multiline","fields":"ID, N1, N2","entries":[{"name":"ID","description":"Element id [int]"},{"name":"N1 - N2","description":"T3 connectivity (2 nodes) [int]"}]}]},{"number":9,"when":"self.keys[\"TYPE\"] in { S3 }","description":"","layouts":[{"range":"1..∞","kind":"multiline","fields":"ID, N1, N2, N3","entries":[{"name":"ID","description":"Element id [int]"},{"name":"N1 - N3","description":"S3 connectivity (3 nodes) [int]"}]}]},{"number":10,"when":"self.keys[\"TYPE\"] in { S4 }","description":"","layouts":[{"range":"1..∞","kind":"multiline","fields":"ID, N1, N2, N3, N4","entries":[{"name":"ID","description":"Element id [int]"},{"name":"N1 - N4","description":"S4 connectivity (4 nodes) [int]"}]}]},{"number":11,"when":"self.keys[\"TYPE\"] in { QSPT }","description":"","layouts":[{"range":"1..∞","kind":"multiline","fields":"ID, N1, N2, N3, N4","entries":[{"name":"ID","description":"Element id [int]"},{"name":"N1 - N4","description":"QSPT connectivity (4 nodes) [int]"}]}]},{"number":12,"when":"self.keys[\"TYPE\"] in { MITC4 }","description":"","layouts":[{"range":"1..∞","kind":"multiline","fields":"ID, N1, N2, N3, N4","entries":[{"name":"ID","description":"Element id [int]"},{"name":"N1 - N4","description":"MITC4 connectivity (4 nodes) [int]"}]}]},{"number":13,"when":"self.keys[\"TYPE\"] in { MITC4FRT }","description":"","layouts":[{"range":"1..∞","kind":"multiline","fields":"ID, N1, N2, N3, N4","entries":[{"name":"ID","description":"Element id [int]"},{"name":"N1 - N4","description":"MITC4FRT connectivity (4 nodes) [int]"}]}]},{"number":14,"when":"self.keys[\"TYPE\"] in { S6 }","description":"","layouts":[{"range":"1..∞","kind":"multiline","fields":"ID, N1, N2, N3, N4, N5, N6","entries":[{"name":"ID","description":"Element id [int]"},{"name":"N1 - N6","description":"S6 connectivity (6 nodes) [int]"}]}]},{"number":15,"when":"self.keys[\"TYPE\"] in { S8 }","description":"","layouts":[{"range":"1..∞","kind":"multiline","fields":"ID, N1, N2, N3, N4, N5, N6, N7, N8","entries":[{"name":"ID","description":"Element id [int]"},{"name":"N1 - N8","description":"S8 connectivity (8 nodes) [int]"}]}]},{"number":16,"when":"self.keys[\"TYPE\"] in { B33 }","description":"","layouts":[{"range":"1..∞","kind":"single-line","fields":"ID, N1, N2, N3","entries":[{"name":"ID","description":"Element id [int]"},{"name":"N1 - N3","description":"B33 connectivity (2 nodes plus optional orientation node) [int]"}]}]},{"number":17,"when":"self.keys[\"TYPE\"] in { C3D5 }","description":"","layouts":[{"range":"1..∞","kind":"multiline","fields":"ID, N1, N2, N3, N4, N5","entries":[{"name":"ID","description":"Element id [int]"},{"name":"N1 - N5","description":"C3D5 connectivity (5 nodes) [int]"}]}]}]},"ELSET":{"description":"Define named element sets.","keywords":[{"name":"ELSET","status":"required","default":"","allowed":[],"aliases":["NAME"],"description":"Name of the element set to activate or create"},{"name":"GENERATE","status":"flag","default":"","allowed":[],"aliases":[],"description":"Interpret rows as start,end[,increment] ranges"}],"admittedRaw":["(none)"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..∞","kind":"single-line","fields":"START, END, INC","entries":[{"name":"START","description":"First element id [int]"},{"name":"END","description":"Last element id [int]"},{"name":"INC","description":"Increment [int]"}]}]},{"number":2,"when":"NOT(self.has_key(\"GENERATE\")) OR","description":"","layouts":[{"range":"0..∞","kind":"single-line","fields":"ID1, ID2, ID3, ID4, ID5, ID6, ID7, ID8, ID9, ID10, ID11, ID12, ID13, ID14, ID15, ID16, ID17, ID18, ID19, ID20, ID21, ID22, ID23, ID24, ID25, ID26, ID27, ID28, ID29, ID30, ID31, ID32","entries":[{"name":"ID1 - ID32","description":"Element ids (up to 32 per line) [int]"}]}]}]},"FIELD":{"description":"Create or populate a generic field (NODE/ELEMENT/ELEMENT_NODAL/ELEMENT_IP).","keywords":[{"name":"COLS","status":"required","default":"","allowed":[],"aliases":[],"description":"Number of components"},{"name":"FILL","status":"optional","default":"ZERO","allowed":["ZERO","NAN"],"aliases":[],"description":"Initialization (ZERO default)"},{"name":"NAME","status":"required","default":"","allowed":[],"aliases":[],"description":"Field name"},{"name":"TYPE","status":"required","default":"","allowed":["NODE","ELEMENT","ELEMENT_NODAL","ELEMENT_IP","IP"],"aliases":[],"description":"Field domain"}],"admittedRaw":["parent.command in { ROOT } OR parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"0..∞","kind":"single-line","fields":"ID, V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, V19, V20, V21, V22, V23, V24, V25, V26, V27, V28, V29, V30, V31, V32, V33, V34, V35, V36, V37, V38, V39, V40, V41, V42, V43, V44, V45, V46, V47, V48, V49, V50, V51, V52, V53, V54, V55, V56, V57, V58, V59, V60, V61, V62, V63, V64","entries":[{"name":"ID","description":"Row id [int]"},{"name":"V1 - V64","description":"Values [string]"}]}]}]},"INERTIALOAD":{"description":"Rigid-body inertial load on element sets. Provide CENTER, CENTER_ACC, OMEGA, ALPHA. TARGET must be an element set. Coriolis term is not included. Optional keyword CONSIDER_POINT_MASSES includes all point-mass features.","keywords":[{"name":"CONSIDER_POINT_MASSES","status":"optional","default":"0","allowed":[],"aliases":[],"description":"If true, include all POINTMASS features in this inertial load"},{"name":"LOAD_COLLECTOR","status":"required","default":"","allowed":[],"aliases":[],"description":"Load collector that stores the inertial loads"}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..∞","kind":"single-line","fields":"TARGET, CENTER1, CENTER2, CENTER3, CENTER_ACC1, CENTER_ACC2, CENTER_ACC3, OMEGA1, OMEGA2, OMEGA3, ALPHA1, ALPHA2, ALPHA3","entries":[{"name":"TARGET","description":"Element set name [string]"},{"name":"CENTER1 - CENTER3","description":"Center point x,y,z [double]"},{"name":"CENTER_ACC1 - CENTER_ACC3","description":"Center linear acceleration ax,ay,az [double]"},{"name":"OMEGA1 - OMEGA3","description":"Angular velocity wx,wy,wz [double]"},{"name":"ALPHA1 - ALPHA3","description":"Angular acceleration ax,ay,az [double]"}]}]}]},"INERTIARELIEF":{"description":"Enable inertia relief for the active linear static load case. CONSIDER_POINT_MASSES controls whether POINTMASS features are included.","keywords":[{"name":"CONSIDER_POINT_MASSES","status":"optional","default":"1","allowed":[],"aliases":[],"description":"If true, inertia relief includes all POINTMASS features"}],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"INITIALVELOCITY":{"description":"Set initial velocity for transient analysis from a node FIELD. Usage: *INITIALVELOCITY, FIELD=NAME","keywords":[{"name":"FIELD","status":"required","default":"","allowed":[],"aliases":[],"description":"Name of a node field with 6 components (ux,uy,uz,rx,ry,rz)"}],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"LOADCASE":{"description":"Begin a load case definition block.","keywords":[{"name":"NAME","status":"optional","default":"","allowed":[],"aliases":[],"description":""},{"name":"TYPE","status":"required","default":"","allowed":[],"aliases":[],"description":"LINEARTRANSIENT, NONLINEARSTATIC}"}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"LOADS":{"description":"Assign load collectors to the active loadcase.","keywords":[],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..∞","kind":"single-line","fields":"LOAD1, LOAD2, LOAD3, LOAD4, LOAD5, LOAD6, LOAD7, LOAD8, LOAD9, LOAD10, LOAD11, LOAD12, LOAD13, LOAD14, LOAD15, LOAD16","entries":[{"name":"LOAD1 - LOAD16","description":"Load collector names [string]"}]}]}]},"MATERIAL":{"description":"Activate or create a material definition.","keywords":[{"name":"MATERIAL","status":"required","default":"","allowed":[],"aliases":["NAME"],"description":"Identifier of the material entry"}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"NEWMARK":{"description":"Set Newmark-β integration parameters (β, γ). Defaults are 0.25, 0.5.","keywords":[],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"One data line: β, γ","layouts":[{"range":"1..1","kind":"single-line","fields":"NEWMARK1, NEWMARK2","entries":[{"name":"NEWMARK1 - NEWMARK2","description":"β, γ parameters for Newmark-β. [double]"}]}]}]},"NODE":{"description":"Define nodes with optional coordinates and assign them to a node set.","keywords":[{"name":"NSET","status":"optional","default":"NALL","allowed":[],"aliases":[],"description":"Target node set to populate (defaults to NALL)"}],"admittedRaw":["(none)"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"0..∞","kind":"single-line","fields":"ID, X, Y, Z","entries":[{"name":"ID","description":"Node identifier [int]"},{"name":"X","description":"X-coordinate [double]"},{"name":"Y","description":"Y-coordinate [double]"},{"name":"Z","description":"Z-coordinate [double]"}]}]}]},"NONLINEAR":{"description":"Configure NONLINEARSTATIC increment and iteration controls.","keywords":[{"name":"INCREMENTS","status":"optional","default":"10","allowed":[],"aliases":[],"description":""},{"name":"MAXITER","status":"optional","default":"20","allowed":[],"aliases":[],"description":""},{"name":"REGULARIZATION_ALPHA","status":"optional","default":"1e-4","allowed":[],"aliases":[],"description":""},{"name":"REGULARIZE_ZERO_ROWS","status":"optional","default":"ON","allowed":[],"aliases":[],"description":""},{"name":"TOL","status":"optional","default":"1e-8","allowed":[],"aliases":[],"description":""}],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"NSET":{"description":"Define named node sets.","keywords":[{"name":"GENERATE","status":"flag","default":"","allowed":[],"aliases":[],"description":"Interpret rows as start,end[,increment] ranges"},{"name":"NSET","status":"required","default":"","allowed":[],"aliases":["NAME"],"description":"Name of the node set to activate or create"}],"admittedRaw":["(none)"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..∞","kind":"single-line","fields":"START, END, INC","entries":[{"name":"START","description":"First node id [int]"},{"name":"END","description":"Last node id [int]"},{"name":"INC","description":"Increment [int]"}]}]},{"number":2,"when":"NOT(self.has_key(\"GENERATE\")) OR","description":"","layouts":[{"range":"0..∞","kind":"single-line","fields":"ID1, ID2, ID3, ID4, ID5, ID6, ID7, ID8, ID9, ID10, ID11, ID12, ID13, ID14, ID15, ID16, ID17, ID18, ID19, ID20, ID21, ID22, ID23, ID24, ID25, ID26, ID27, ID28, ID29, ID30, ID31, ID32","entries":[{"name":"ID1 - ID32","description":"Node ids (up to 32 per line) [int]"}]}]}]},"NUMEIGENVALUES":{"description":"Set number of eigenvalues for buckling/eigenfrequency loadcases.","keywords":[],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..1","kind":"single-line","fields":"COUNT","entries":[{"name":"COUNT","description":"Number of eigenvalues [int]"}]}]}]},"ORIENTATION":{"description":"Define coordinate systems (rectangular or cylindrical).","keywords":[{"name":"DEFINITION","status":"optional","default":"","allowed":[],"aliases":[],"description":"Unused legacy parameter"},{"name":"NAME","status":"required","default":"","allowed":[],"aliases":[],"description":"Coordinate system identifier"},{"name":"TYPE","status":"required","default":"","allowed":["RECTANGULAR","CYLINDRICAL"],"aliases":[],"description":""}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"self.keys[\"TYPE\"] in { RECTANGULAR }","description":"","layouts":[{"range":"1..3","kind":"multiline","fields":"DATA1, DATA2, DATA3, DATA4, DATA5, DATA6, DATA7, DATA8, DATA9","entries":[{"name":"DATA1 - DATA9","description":"Rectangular system vectors [double]"}]}]},{"number":2,"when":"self.keys[\"TYPE\"] in { RECTANGULAR }","description":"","layouts":[{"range":"1..2","kind":"multiline","fields":"DATA1, DATA2, DATA3, DATA4, DATA5, DATA6","entries":[{"name":"DATA1 - DATA6","description":"Rectangular system vectors (two vectors) [double]"}]}]},{"number":3,"when":"self.keys[\"TYPE\"] in { RECTANGULAR }","description":"","layouts":[{"range":"1..1","kind":"single-line","fields":"DATA1, DATA2, DATA3","entries":[{"name":"DATA1 - DATA3","description":"Rectangular system vector [double]"}]}]},{"number":4,"when":"self.keys[\"TYPE\"] in { CYLINDRICAL }","description":"","layouts":[{"range":"1..3","kind":"multiline","fields":"DATA1, DATA2, DATA3, DATA4, DATA5, DATA6, DATA7, DATA8, DATA9","entries":[{"name":"DATA1 - DATA9","description":"Cylindrical system vectors [double]"}]}]}]},"OVERVIEW":{"description":"Print a summary/overview of the current model (nodes, elements, sets, materials, sections, couplings).","keywords":[],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"PLOAD":{"description":"Pressure loads on surfaces. AMPLITUDE may reference a time history that scales the pressure magnitude.","keywords":[{"name":"AMPLITUDE","status":"optional","default":"","allowed":[],"aliases":[],"description":"Optional time amplitude used to scale the pressure value"},{"name":"LOAD_COLLECTOR","status":"required","default":"","allowed":[],"aliases":[],"description":"Target load collector that groups the loads"}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..∞","kind":"single-line","fields":"TARGET, P","entries":[{"name":"TARGET","description":"Surface set or surface id [string]"},{"name":"P","description":"Pressure value [double]"}]}]}]},"POINTMASS":{"description":"Assign point-mass feature to a node set.","keywords":[{"name":"NSET","status":"required","default":"","allowed":[],"aliases":[],"description":"Target node set"}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..∞","kind":"multiline","fields":"MASS, INERTIA1, INERTIA2, INERTIA3, SPRING1, SPRING2, SPRING3, ROTSPRING1, ROTSPRING2, ROTSPRING3","entries":[{"name":"MASS","description":"Point mass [double]"},{"name":"INERTIA1 - INERTIA3","description":"Rotary inertia Ix,Iy,Iz [double]"},{"name":"SPRING1 - SPRING3","description":"Translational spring constants [double]"},{"name":"ROTSPRING1 - ROTSPRING3","description":"Rotational spring constants [double]"}]}]}]},"PROFILE":{"description":"Define beam profile properties in this order: A, Iy, Iz, Jt, Iyz, ey, ez, refy, refz. Only the first 4 are required. Convention: Iyz = integral_A(y*z*dA), i.e. without a leading minus sign.","keywords":[{"name":"PROFILE","status":"required","default":"","allowed":[],"aliases":["NAME"],"description":"Identifier of the profile"}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..1","kind":"single-line","fields":"A, IY, IZ, JT, IYZ, EY, EZ, REFY, REFZ","entries":[{"name":"A","description":"Cross-section area A [double]"},{"name":"IY","description":"Second moment of area about local y-axis (Iy) [double]"},{"name":"IZ","description":"Second moment of area about local z-axis (Iz) [double]"},{"name":"JT","description":"Torsional constant (Jt) [double]"},{"name":"IYZ","description":"Product of inertia: Iyz = integral_A(y*z*dA), no minus sign [double]"},{"name":"EY","description":"Offset in local y: ey = y(SP) - y(SMP) [double]"},{"name":"EZ","description":"Offset in local z: ez = z(SP) - z(SMP) [double]"},{"name":"REFY","description":"Reference-line offset in local y: refy = y(REF) - y(SMP) [double]"},{"name":"REFZ","description":"Reference-line offset in local z: refz = z(REF) - z(SMP) [double]"}]}]}]},"RBM":{"description":"Add rigid-body-motion suppression equations for an element set. ELSET defaults to EALL.","keywords":[{"name":"ELSET","status":"optional","default":"EALL","allowed":[],"aliases":["SET"],"description":"Target element set used to build RBM equations"}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"REBALANCELOADS":{"description":"Enable rigid-body rebalancing of external loads (sum F=M=0) for the active linear static load case.","keywords":[],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"REQUESTSTGEOM":{"description":"Request geometric stiffness output for LINEARBUCKLING loadcases.","keywords":[{"name":"FILE","status":"optional","default":"","allowed":[],"aliases":[],"description":""}],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"REQUESTSTIFFNESS":{"description":"Request stiffness matrix output for supported loadcases.","keywords":[{"name":"FILE","status":"optional","default":"","allowed":[],"aliases":[],"description":""}],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"SFSET":{"description":"Create or activate a surface set and add surface IDs (lists or ranges).","keywords":[{"name":"GENERATE","status":"flag","default":"","allowed":[],"aliases":[],"description":"Interpret data lines as ranges start,end[,inc]."},{"name":"SFSET","status":"optional","default":"SFALL","allowed":[],"aliases":["NAME"],"description":"Surface set name (default: SFALL)."}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"Range rows with explicit increment: start, end, inc.","layouts":[{"range":"0..∞","kind":"single-line","fields":"RANGE31, RANGE32, RANGE33","entries":[{"name":"RANGE31 - RANGE33","description":"start, end, inc [int]"}]}]},{"number":2,"when":"(none)","description":"Range rows with implicit increment: start, end (inc=1).","layouts":[{"range":"0..∞","kind":"single-line","fields":"RANGE21, RANGE22","entries":[{"name":"RANGE21 - RANGE22","description":"start, end (inc=1) [int]"}]}]},{"number":3,"when":"NOT(self.has_key(\"GENERATE\")) OR","description":"Explicit surface IDs (one or more per line).","layouts":[{"range":"0..∞","kind":"single-line","fields":"IDS1, IDS2, IDS3, IDS4, IDS5, IDS6, IDS7, IDS8, IDS9, IDS10, IDS11, IDS12, IDS13, IDS14, IDS15, IDS16, IDS17, IDS18, IDS19, IDS20, IDS21, IDS22, IDS23, IDS24, IDS25, IDS26, IDS27, IDS28, IDS29, IDS30, IDS31, IDS32","entries":[{"name":"IDS1 - IDS32","description":"Surface IDs (up to 32 per line). [int]"}]}]}]},"SHELLSECTION":{"description":"Assign a shell section thickness to an element set.","keywords":[{"name":"ELSET","status":"required","default":"","allowed":[],"aliases":[],"description":"Target element set"},{"name":"MATERIAL","status":"required","default":"","allowed":[],"aliases":["MAT"],"description":"Material name"},{"name":"ORIENTATION","status":"optional","default":"","allowed":[],"aliases":[],"description":"Optional coordinate system for shell n1/n2 material/resultant axes"}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..1","kind":"single-line","fields":"THICKNESS","entries":[{"name":"THICKNESS","description":"Shell thickness [double]"}]}]}]},"SIGMA":{"description":"Set eigen-shift parameter for LINEARBUCKLING loadcases.","keywords":[],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..1","kind":"single-line","fields":"SIGMA","entries":[{"name":"SIGMA","description":"Shift parameter [double]"}]}]}]},"SOLIDSECTION":{"description":"Assign a material to a solid element set.","keywords":[{"name":"ELSET","status":"required","default":"","allowed":[],"aliases":[],"description":"Target element set"},{"name":"MATERIAL","status":"required","default":"","allowed":[],"aliases":["MAT"],"description":"Material name"},{"name":"ORIENTATION","status":"optional","default":"","allowed":[],"aliases":[],"description":"Optional coordinate system for solid material axes"}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"SOLVER":{"description":"Configure solver options for the active loadcase. Constraint | Backend | DIRECT | INDIRECT NULLSPACE | CPU MKL | Yes | Yes NULLSPACE | CPU Eigen | Yes | Yes NULLSPACE | GPU | Yes | Yes NULLSPACE | GPU cuDSS | Yes | Yes LAGRANGE | CPU MKL | Yes | No LAGRANGE | CPU Eigen | Limited | No LAGRANGE | GPU | No | No LAGRANGE | GPU cuDSS | Yes | No","keywords":[{"name":"DEVICE","status":"optional","default":"CPU","allowed":["CPU","GPU"],"aliases":[],"description":""},{"name":"METHOD","status":"optional","default":"DIRECT","allowed":["DIRECT","INDIRECT"],"aliases":[],"description":""}],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"SUPPORT":{"description":"Define nodal supports via support collectors.","keywords":[{"name":"ORIENTATION","status":"optional","default":"","allowed":[],"aliases":[],"description":"Optional orientation coordinate system"},{"name":"SUPPORT_COLLECTOR","status":"required","default":"","allowed":[],"aliases":[],"description":"Support collector name"}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..∞","kind":"single-line","fields":"TARGET, DOF1, DOF2, DOF3, DOF4, DOF5, DOF6","entries":[{"name":"TARGET","description":"Node set or id [string]"},{"name":"DOF1 - DOF6","description":"Support values for ux,uy,uz,rx,ry,rz [double]"}]}]}]},"SUPPORTS":{"description":"Assign support collectors to the active loadcase.","keywords":[],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..∞","kind":"single-line","fields":"SUPP1, SUPP2, SUPP3, SUPP4, SUPP5, SUPP6, SUPP7, SUPP8, SUPP9, SUPP10, SUPP11, SUPP12, SUPP13, SUPP14, SUPP15, SUPP16","entries":[{"name":"SUPP1 - SUPP16","description":"Support collector names [string]"}]}]}]},"SURFACE":{"description":"Define surfaces from element faces (by element id) or entire element sets. TYPE must be ELEMENT. Use SFSET/NAME to select the active surface set.","keywords":[{"name":"SFSET","status":"optional","default":"SFALL","allowed":[],"aliases":["NAME"],"description":"Set name to activate/create (default: SFALL). Elements with 2D faces contribute surfaces; beam elements contribute 1D lines under the same set name."}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"Three fields per line: ID, ELEM_ID, SIDE (SIDE accepts 'S#' or integer).","layouts":[{"range":"1..∞","kind":"single-line","fields":"ID, ELEM_ID, SIDE","entries":[{"name":"ID","description":"Surface id [int]"},{"name":"ELEM_ID","description":"Element id [int]"},{"name":"SIDE","description":"Face id, e.g. S1 or 1 [string]"}]}]},{"number":2,"when":"(none)","description":"Two fields per line: TARGET, SIDE (TARGET = ELSET name or single element id).","layouts":[{"range":"1..∞","kind":"single-line","fields":"TARGET, SIDE","entries":[{"name":"TARGET","description":"ELSET name or element id (int) [string]"},{"name":"SIDE","description":"Face id, e.g. S1 or 1 [string]"}]}]}]},"THERMALEXPANSION":{"description":"Assign thermal expansion coefficient to the active material.","keywords":[],"admittedRaw":["parent.command in { MATERIAL }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..1","kind":"single-line","fields":"ALPHA","entries":[{"name":"ALPHA","description":"Thermal expansion coefficient [double]"}]}]}]},"TIE":{"description":"Bind a slave set (node or surface set) to a master surface/line set (tie constraint).","keywords":[{"name":"ADJUST","status":"optional","default":"NO","allowed":["NO","YES"],"aliases":[],"description":"Determines projection of slave to master"},{"name":"DISTANCE","status":"required","default":"","allowed":[],"aliases":[],"description":"Search distance"},{"name":"MASTER","status":"required","default":"","allowed":[],"aliases":[],"description":"Master surface set or line set"},{"name":"SLAVE","status":"required","default":"","allowed":[],"aliases":[],"description":"Slave node set or slave surface set"}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"TIME":{"description":"Set time window and step for transient analysis. Accepts either: (t_start, t_end, dt) or (t_end, dt).","keywords":[],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"One data line: t_start, t_end, dt","layouts":[{"range":"1..1","kind":"single-line","fields":"TIME1, TIME2, TIME3","entries":[{"name":"TIME1 - TIME3","description":"t_start, t_end, dt [double]"}]}]},{"number":2,"when":"(none)","description":"One data line: t_end, dt","layouts":[{"range":"1..1","kind":"single-line","fields":"TIME1, TIME2","entries":[{"name":"TIME1 - TIME2","description":"t_end, dt [double]"}]}]}]},"TLOAD":{"description":"Thermal load referencing a temperature field.","keywords":[{"name":"LOAD_COLLECTOR","status":"required","default":"","allowed":[],"aliases":[],"description":"Target load collector that groups the loads"},{"name":"REFERENCETEMPERATURE","status":"required","default":"","allowed":[],"aliases":[],"description":"Reference temperature value"},{"name":"TEMPERATUREFIELD","status":"required","default":"","allowed":[],"aliases":[],"description":"Name of the temperature field to apply"}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"TOPODENSITY":{"description":"Select the element density field for LINEARSTATICTOPO loadcases.","keywords":[{"name":"FIELD","status":"required","default":"","allowed":[],"aliases":["NAME"],"description":"Density field name (ELEMENT, 1 component)"}],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"TOPOEXPONENT":{"description":"Set penalization exponent for LINEARSTATICTOPO loadcases.","keywords":[],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..1","kind":"single-line","fields":"EXPONENT","entries":[{"name":"EXPONENT","description":"Penalization exponent [double]"}]}]}]},"TOPOORIENT":{"description":"Select the orientation field for LINEARSTATICTOPO loadcases.","keywords":[{"name":"FIELD","status":"required","default":"","allowed":[],"aliases":["NAME"],"description":"Orientation field name (ELEMENT, 3 components)"}],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[]}]},"TRUSSSECTION":{"description":"Assign a truss section area to an element set.","keywords":[{"name":"ELSET","status":"required","default":"","allowed":[],"aliases":[],"description":"Target truss element set"},{"name":"MATERIAL","status":"required","default":"","allowed":[],"aliases":["MAT"],"description":"Material name"}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..1","kind":"single-line","fields":"AREA","entries":[{"name":"AREA","description":"Truss cross-sectional area [double]"}]}]}]},"VLOAD":{"description":"Apply body-force vectors to element sets or individual elements. Each record supplies Fx, Fy, Fz components that are distributed according to the element formulation. Specify ORIENTATION to express the components in a local coordinate system (e.g., gravity in cylindrical axes); otherwise the global basis is used. AMPLITUDE allows scaling the body force with a time history.","keywords":[{"name":"AMPLITUDE","status":"optional","default":"","allowed":[],"aliases":[],"description":"Optional time amplitude used to scale the body-force components"},{"name":"LOAD_COLLECTOR","status":"required","default":"","allowed":[],"aliases":[],"description":"Load collector that stores the volumetric loads"},{"name":"ORIENTATION","status":"optional","default":"","allowed":[],"aliases":[],"description":"Optional coordinate system for interpreting the body-force vector"}],"admittedRaw":["parent.command in { ROOT }"],"variants":[{"number":1,"when":"(none)","description":"","layouts":[{"range":"1..∞","kind":"single-line","fields":"TARGET, LOAD1, LOAD2, LOAD3","entries":[{"name":"TARGET","description":"Element set or element id [string]"},{"name":"LOAD1 - LOAD3","description":"Fx,Fy,Fz components, local if ORIENTATION is set [double]"}]}]}]},"WRITEEVERY":{"description":"Control output frequency during transient analysis. Use TYPE=STEPS with an integer N (write every N steps), or TYPE=TIME with a positive Δt_write in seconds (overrides steps).","keywords":[{"name":"TYPE","status":"optional","default":"STEPS","allowed":["STEPS","TIME"],"aliases":[],"description":""}],"admittedRaw":["parent.command in { LOADCASE }"],"variants":[{"number":1,"when":"(none)","description":"One data line: either integer steps (TYPE=STEPS) or Δt_write seconds (TYPE=TIME).","layouts":[{"range":"1..1","kind":"single-line","fields":"VALUE","entries":[{"name":"VALUE","description":"N (steps) or Δt_write (seconds), depending on TYPE. [double]"}]}]}]}};


const commandManualDocs = {
  "DLOAD": {
    "title": "Distributed surface load",
    "summary": "Apply distributed vector tractions to surfaces.",
    "category": "Loads and supports",
    "badges": [
      "Loads",
      "Root command",
      "Collector-based",
      "Surface target",
      "Supports orientation",
      "Supports amplitude"
    ],
    "quickFacts": [
      [
        "Valid in",
        "Root level"
      ],
      [
        "Creates",
        "Distributed load entries in a load collector"
      ],
      [
        "Activated by",
        "*LOADS inside *LOADCASE"
      ],
      [
        "Targets",
        "Surface id or surface set"
      ],
      [
        "Components",
        "load_1, load_2, load_3"
      ],
      [
        "Internal effect",
        "Integrated into the global external force vector during assembly"
      ]
    ],
    "syntax": [
      {
        "title": "Global surface traction",
        "code": "*DLOAD, LOAD_COLLECTOR=<collector>\n<surface>, <Fx>, <Fy>, <Fz>",
        "description": "Applies a distributed vector traction in global coordinates."
      },
      {
        "title": "Local surface traction",
        "code": "*DLOAD, LOAD_COLLECTOR=<collector>, ORIENTATION=<orientation>\n<surface>, <F1>, <F2>, <F3>",
        "description": "Interprets the components in a local coordinate system before transforming them to global coordinates."
      },
      {
        "title": "Transient-scaled traction",
        "code": "*DLOAD, LOAD_COLLECTOR=<collector>, AMPLITUDE=<amplitude>\n<surface>, <Fx>, <Fy>, <Fz>",
        "description": "References an amplitude curve so the load can be scaled during transient analysis."
      }
    ],
    "examples": [
      {
        "title": "Minimal distributed surface traction",
        "code": "*DLOAD, LOAD_COLLECTOR=LOADS\nSURF_TOP, 0.0, 0.0, -1.0\n\n*LOADCASE, TYPE=LINEARSTATIC\n  *LOADS\n  LOADS\n  *SUPPORTS\n  BCS",
        "description": "This stores a distributed vector traction on SURF_TOP in the LOADS collector. It is assembled only when that collector is activated in the loadcase."
      }
    ],
    "keywordsManual": {
      "required": [
        {
          "name": "LOAD_COLLECTOR",
          "values": "collector name",
          "default": "\u2014",
          "description": "Name of the load collector that receives the distributed load entries."
        }
      ],
      "optional": [
        {
          "name": "ORIENTATION",
          "values": "orientation name",
          "default": "global basis",
          "description": "Interprets load components in a local coordinate system before assembly."
        },
        {
          "name": "AMPLITUDE",
          "values": "amplitude name",
          "default": "constant factor 1",
          "description": "Scales the load during transient or incremental analyses."
        }
      ],
      "groups": []
    },
    "datalinesManual": [
      {
        "title": "Surface traction entry",
        "signature": "<target>, <load_1>, <load_2>, <load_3>",
        "rows": "one or more",
        "fields": [
          {
            "name": "target",
            "meaning": "Surface id or surface set",
            "type": "string/int",
            "notes": "Region where the traction is applied."
          },
          {
            "name": "load_1",
            "meaning": "First traction component",
            "type": "float",
            "notes": "Global Fx or local F1 when ORIENTATION is set."
          },
          {
            "name": "load_2",
            "meaning": "Second traction component",
            "type": "float",
            "notes": "Global Fy or local F2 when ORIENTATION is set."
          },
          {
            "name": "load_3",
            "meaning": "Third traction component",
            "type": "float",
            "notes": "Global Fz or local F3 when ORIENTATION is set."
          }
        ]
      }
    ],
    "semantics": [
      "FEMaster stores each dataline as a distributed load entry inside the selected load collector.",
      "The command does not immediately assemble a global force vector. During a loadcase, activated load collectors are evaluated and integrated over the referenced surfaces.",
      "If ORIENTATION is provided, components are interpreted in that local coordinate system and transformed to global coordinates before assembly."
    ],
    "validation": [
      "LOAD_COLLECTOR must be provided.",
      "The target surface id or surface set must exist.",
      "ORIENTATION must reference an existing orientation if provided.",
      "AMPLITUDE must reference an existing amplitude if provided.",
      "Each dataline must contain one target and three numeric components."
    ],
    "mistakes": [
      "Defining the load collector but never activating it with *LOADS.",
      "Using *DLOAD for normal pressure. Use *PLOAD when the load should follow the surface normal.",
      "Forgetting that ORIENTATION changes the meaning of load_1, load_2 and load_3.",
      "Applying the load to a surface with an unexpected side or normal definition."
    ],
    "related": [
      {
        "command": "SURFACE",
        "description": "Defines element faces that can receive surface loads."
      },
      {
        "command": "SFSET",
        "description": "Groups surfaces for reuse as load targets."
      },
      {
        "command": "LOADS",
        "description": "Activates load collectors inside a loadcase."
      },
      {
        "command": "LOADCASE",
        "description": "Runs the analysis where the load is assembled."
      },
      {
        "command": "ORIENTATION",
        "description": "Defines optional local load directions."
      },
      {
        "command": "PLOAD",
        "description": "Use for scalar pressure normal to a surface."
      }
    ]
  },
  "PLOAD": {
    "title": "Pressure load",
    "summary": "Apply scalar pressure to surfaces using the selected surface normal.",
    "category": "Loads and supports",
    "badges": [
      "Loads",
      "Root command",
      "Collector-based",
      "Surface target",
      "Normal pressure"
    ],
    "quickFacts": [
      [
        "Valid in",
        "Root level"
      ],
      [
        "Creates",
        "Pressure load entries in a load collector"
      ],
      [
        "Activated by",
        "*LOADS inside *LOADCASE"
      ],
      [
        "Targets",
        "Surface id or surface set"
      ],
      [
        "Components",
        "pressure magnitude"
      ],
      [
        "Internal effect",
        "Converted to equivalent nodal forces during surface integration"
      ]
    ],
    "syntax": [
      {
        "title": "Surface pressure",
        "code": "*PLOAD, LOAD_COLLECTOR=<collector>\n<surface>, <pressure>",
        "description": "Applies pressure using the stored surface normal."
      },
      {
        "title": "Amplitude-scaled pressure",
        "code": "*PLOAD, LOAD_COLLECTOR=<collector>, AMPLITUDE=<amplitude>\n<surface>, <pressure>",
        "description": "Scales the pressure magnitude by an amplitude curve during transient analysis."
      }
    ],
    "examples": [
      {
        "title": "Minimal pressure load",
        "code": "*PLOAD, LOAD_COLLECTOR=PRESSURE\nSKIN_OUTER, 2.5\n\n*LOADCASE, TYPE=LINEARSTATIC\n  *LOADS\n  PRESSURE",
        "description": "This defines a pressure on SKIN_OUTER and activates the pressure collector in a static loadcase."
      }
    ],
    "keywordsManual": {
      "required": [
        {
          "name": "LOAD_COLLECTOR",
          "values": "collector name",
          "default": "\u2014",
          "description": "Load collector that receives the pressure entries."
        }
      ],
      "optional": [
        {
          "name": "AMPLITUDE",
          "values": "amplitude name",
          "default": "constant factor 1",
          "description": "Scales the pressure in transient or incremental analyses."
        }
      ],
      "groups": []
    },
    "datalinesManual": [
      {
        "title": "Pressure entry",
        "signature": "<target>, <pressure>",
        "rows": "one or more",
        "fields": [
          {
            "name": "target",
            "meaning": "Surface id or surface set",
            "type": "string/int",
            "notes": "Region where the pressure is applied."
          },
          {
            "name": "pressure",
            "meaning": "Pressure magnitude",
            "type": "float",
            "notes": "Sign follows FEMaster's surface normal convention."
          }
        ]
      }
    ],
    "semantics": [
      "FEMaster stores pressure entries in a load collector and converts them into equivalent nodal forces when the collector is assembled.",
      "The pressure direction is determined from the target surface normal, not from a user-supplied vector."
    ],
    "validation": [
      "LOAD_COLLECTOR must be provided.",
      "The target surface id or surface set must exist.",
      "AMPLITUDE must reference an existing amplitude if provided.",
      "Each dataline must contain one target and one numeric pressure value."
    ],
    "mistakes": [
      "Using *PLOAD when a fixed global vector is intended. Use *DLOAD for vector tractions.",
      "Unexpected pressure sign caused by element side or shell side orientation.",
      "Defining a pressure collector without activating it in a loadcase."
    ],
    "related": [
      {
        "command": "SURFACE",
        "description": "Defines target faces and their side labels."
      },
      {
        "command": "SFSET",
        "description": "Groups surfaces for pressure application."
      },
      {
        "command": "DLOAD",
        "description": "Use for vector surface tractions."
      },
      {
        "command": "LOADS",
        "description": "Activates the load collector."
      }
    ]
  },
  "CLOAD": {
    "title": "Concentrated nodal load",
    "summary": "Apply nodal forces and moments to nodes or node sets.",
    "category": "Loads and supports",
    "badges": [
      "Loads",
      "Root command",
      "Collector-based",
      "Node target",
      "Supports orientation",
      "Supports amplitude"
    ],
    "quickFacts": [
      [
        "Valid in",
        "Root level"
      ],
      [
        "Creates",
        "Concentrated nodal load entries"
      ],
      [
        "Activated by",
        "*LOADS inside *LOADCASE"
      ],
      [
        "Targets",
        "Node id or node set"
      ],
      [
        "Components",
        "Fx, Fy, Fz, Mx, My, Mz"
      ],
      [
        "Internal effect",
        "Adds directly to the external force vector for active DOFs"
      ]
    ],
    "syntax": [
      {
        "title": "Global nodal force/moment",
        "code": "*CLOAD, LOAD_COLLECTOR=<collector>\n<node-or-nset>, <Fx>, <Fy>, <Fz>, <Mx>, <My>, <Mz>",
        "description": "Applies six global force and moment components to a node or node set."
      },
      {
        "title": "Local nodal load",
        "code": "*CLOAD, LOAD_COLLECTOR=<collector>, ORIENTATION=<orientation>\n<node-or-nset>, <F1>, <F2>, <F3>, <M1>, <M2>, <M3>",
        "description": "Interprets components in the given local coordinate system."
      }
    ],
    "examples": [
      {
        "title": "Pull a node set in x direction",
        "code": "*CLOAD, LOAD_COLLECTOR=FORCE\nTIP_NODES, 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0\n\n*LOADCASE, TYPE=LINEARSTATIC\n  *LOADS\n  FORCE",
        "description": "The load is stored in FORCE and applied only by the loadcase that activates FORCE."
      }
    ],
    "keywordsManual": {
      "required": [
        {
          "name": "LOAD_COLLECTOR",
          "values": "collector name",
          "default": "\u2014",
          "description": "Collector receiving the nodal load entries."
        }
      ],
      "optional": [
        {
          "name": "ORIENTATION",
          "values": "orientation name",
          "default": "global basis",
          "description": "Interprets force and moment components locally."
        },
        {
          "name": "AMPLITUDE",
          "values": "amplitude name",
          "default": "constant factor 1",
          "description": "Scales the load in transient analyses."
        }
      ],
      "groups": []
    },
    "datalinesManual": [
      {
        "title": "Nodal load entry",
        "signature": "<target>, <Fx>, <Fy>, <Fz>, <Mx>, <My>, <Mz>",
        "rows": "one or more",
        "fields": [
          {
            "name": "target",
            "meaning": "Node id or node set",
            "type": "string/int",
            "notes": "All addressed nodes receive the listed components."
          },
          {
            "name": "Fx, Fy, Fz",
            "meaning": "Force components",
            "type": "float",
            "notes": "Global or local translations depending on ORIENTATION."
          },
          {
            "name": "Mx, My, Mz",
            "meaning": "Moment components",
            "type": "float",
            "notes": "Only active rotational DOFs can receive moments."
          }
        ]
      }
    ],
    "semantics": [
      "FEMaster accumulates concentrated load records into a load collector.",
      "During assembly, active collectors contribute to the global external force vector. Repeated records on the same node or set add together."
    ],
    "validation": [
      "LOAD_COLLECTOR must be provided.",
      "The target node or node set must exist.",
      "Each dataline must contain one target and six numeric components.",
      "ORIENTATION and AMPLITUDE must reference existing objects if provided."
    ],
    "mistakes": [
      "Applying moments to solid-only nodes without active rotational DOFs.",
      "Forgetting to activate the load collector in *LOADS.",
      "Using local components without defining ORIENTATION."
    ],
    "related": [
      {
        "command": "NSET",
        "description": "Defines reusable node targets."
      },
      {
        "command": "LOADS",
        "description": "Activates load collectors."
      },
      {
        "command": "ORIENTATION",
        "description": "Defines local load directions."
      },
      {
        "command": "AMPLITUDE",
        "description": "Defines time-dependent scale factors."
      }
    ]
  },
  "VLOAD": {
    "title": "Volume load",
    "summary": "Apply body-force vectors to element regions.",
    "category": "Loads and supports",
    "badges": [
      "Loads",
      "Root command",
      "Collector-based",
      "Element target",
      "Supports orientation",
      "Supports amplitude"
    ],
    "quickFacts": [
      [
        "Valid in",
        "Root level"
      ],
      [
        "Creates",
        "Volume load entries"
      ],
      [
        "Activated by",
        "*LOADS inside *LOADCASE"
      ],
      [
        "Targets",
        "Element id or element set"
      ],
      [
        "Components",
        "Fx, Fy, Fz"
      ],
      [
        "Internal effect",
        "Integrated over element volume or formulation domain"
      ]
    ],
    "syntax": [
      {
        "title": "Global body force",
        "code": "*VLOAD, LOAD_COLLECTOR=<collector>\n<element-or-elset>, <Fx>, <Fy>, <Fz>",
        "description": "Applies a body-force vector in global coordinates."
      },
      {
        "title": "Local body force",
        "code": "*VLOAD, LOAD_COLLECTOR=<collector>, ORIENTATION=<orientation>\n<element-or-elset>, <F1>, <F2>, <F3>",
        "description": "Interprets the vector in a local coordinate system."
      }
    ],
    "examples": [
      {
        "title": "Gravity-like volume load",
        "code": "*VLOAD, LOAD_COLLECTOR=GRAVITY\nALL_SOLIDS, 0.0, 0.0, -9.81\n\n*LOADCASE, TYPE=LINEARSTATIC\n  *LOADS\n  GRAVITY",
        "description": "The body force is stored in GRAVITY and assembled when GRAVITY is activated."
      }
    ],
    "keywordsManual": {
      "required": [
        {
          "name": "LOAD_COLLECTOR",
          "values": "collector name",
          "default": "\u2014",
          "description": "Collector receiving volume load entries."
        }
      ],
      "optional": [
        {
          "name": "ORIENTATION",
          "values": "orientation name",
          "default": "global basis",
          "description": "Interprets body-force components locally."
        },
        {
          "name": "AMPLITUDE",
          "values": "amplitude name",
          "default": "constant factor 1",
          "description": "Scales the body force in transient analyses."
        }
      ],
      "groups": []
    },
    "datalinesManual": [
      {
        "title": "Volume load entry",
        "signature": "<target>, <load_1>, <load_2>, <load_3>",
        "rows": "one or more",
        "fields": [
          {
            "name": "target",
            "meaning": "Element id or element set",
            "type": "string/int",
            "notes": "Region receiving the body force."
          },
          {
            "name": "load_1..load_3",
            "meaning": "Body-force vector components",
            "type": "float",
            "notes": "Global Fx/Fy/Fz or local F1/F2/F3."
          }
        ]
      }
    ],
    "semantics": [
      "FEMaster stores body-force records in a load collector and evaluates them over the referenced element region during load assembly."
    ],
    "validation": [
      "LOAD_COLLECTOR must be provided.",
      "The target element or element set must exist.",
      "Each dataline must contain one target and three numeric components."
    ],
    "mistakes": [
      "Confusing acceleration units with force-density units in a chosen unit system.",
      "Targeting shell or beam regions where the intended volume interpretation is not appropriate.",
      "Forgetting to activate the collector with *LOADS."
    ],
    "related": [
      {
        "command": "ELSET",
        "description": "Defines element regions."
      },
      {
        "command": "LOADS",
        "description": "Activates volume load collectors."
      },
      {
        "command": "INERTIALOAD",
        "description": "Use for rigid-body acceleration effects."
      }
    ]
  },
  "SUPPORT": {
    "title": "Support definition",
    "summary": "Store prescribed nodal components in a support collector.",
    "category": "Loads and supports",
    "badges": [
      "Boundary conditions",
      "Root command",
      "Collector-based",
      "Node target",
      "Supports orientation"
    ],
    "quickFacts": [
      [
        "Valid in",
        "Root level"
      ],
      [
        "Creates",
        "Support entries in a support collector"
      ],
      [
        "Activated by",
        "*SUPPORTS inside *LOADCASE"
      ],
      [
        "Targets",
        "Node id or node set"
      ],
      [
        "Components",
        "Ux, Uy, Uz, Rx, Ry, Rz"
      ],
      [
        "Internal effect",
        "Creates prescribed DOF constraints"
      ]
    ],
    "syntax": [
      {
        "title": "Global support values",
        "code": "*SUPPORT, SUPPORT_COLLECTOR=<collector>\n<node-or-nset>, <Ux>, <Uy>, <Uz>, <Rx>, <Ry>, <Rz>",
        "description": "Prescribes components in the global coordinate system."
      },
      {
        "title": "Local support values",
        "code": "*SUPPORT, SUPPORT_COLLECTOR=<collector>, ORIENTATION=<orientation>\n<node-or-nset>, <U1>, <U2>, <U3>, <R1>, <R2>, <R3>",
        "description": "Interprets prescribed components in a local coordinate system."
      }
    ],
    "examples": [
      {
        "title": "Fix a node set",
        "code": "*SUPPORT, SUPPORT_COLLECTOR=BCS\nFIXED_END, 0, 0, 0, 0, 0, 0\n\n*LOADCASE, TYPE=LINEARSTATIC\n  *SUPPORTS\n  BCS",
        "description": "The support is stored in BCS and active only in loadcases that list BCS under *SUPPORTS."
      }
    ],
    "keywordsManual": {
      "required": [
        {
          "name": "SUPPORT_COLLECTOR",
          "values": "collector name",
          "default": "\u2014",
          "description": "Support collector receiving the boundary condition entries."
        }
      ],
      "optional": [
        {
          "name": "ORIENTATION",
          "values": "orientation name",
          "default": "global basis",
          "description": "Interprets support components in local axes."
        }
      ],
      "groups": []
    },
    "datalinesManual": [
      {
        "title": "Support entry",
        "signature": "<target>, <Ux>, <Uy>, <Uz>, <Rx>, <Ry>, <Rz>",
        "rows": "one or more",
        "fields": [
          {
            "name": "target",
            "meaning": "Node id or node set",
            "type": "string/int",
            "notes": "Nodes whose components are prescribed."
          },
          {
            "name": "Ux, Uy, Uz",
            "meaning": "Prescribed translations",
            "type": "float/blank",
            "notes": "Use zero for fixed translations."
          },
          {
            "name": "Rx, Ry, Rz",
            "meaning": "Prescribed rotations",
            "type": "float/blank",
            "notes": "Relevant only where rotational DOFs exist."
          }
        ]
      }
    ],
    "semantics": [
      "FEMaster stores supports as collector entries. When a loadcase activates the collector, the entries become prescribed DOF constraints in the solver system."
    ],
    "validation": [
      "SUPPORT_COLLECTOR must be provided.",
      "The target node or node set must exist.",
      "Each dataline must contain one target and six component fields.",
      "ORIENTATION must exist if provided."
    ],
    "mistakes": [
      "Defining supports but not activating the support collector in *SUPPORTS.",
      "Constraining rotations on nodes that never receive rotational DOFs.",
      "Leaving a model under-constrained while expecting a static solve to be regular."
    ],
    "related": [
      {
        "command": "NSET",
        "description": "Defines support targets."
      },
      {
        "command": "SUPPORTS",
        "description": "Activates support collectors in a loadcase."
      },
      {
        "command": "CONSTRAINTMETHOD",
        "description": "Controls how constraints are handled by supported solvers."
      }
    ]
  },
  "LOADCASE": {
    "title": "Loadcase block",
    "summary": "Start an analysis block and execute it when the block closes.",
    "category": "Analysis",
    "badges": [
      "Analysis",
      "Context command",
      "Executes solve"
    ],
    "quickFacts": [
      [
        "Valid in",
        "Root level"
      ],
      [
        "Creates",
        "An analysis/loadcase context"
      ],
      [
        "Contains",
        "Loads, supports, solver, and analysis controls"
      ],
      [
        "Depends on",
        "Previously defined model data"
      ],
      [
        "Internal effect",
        "Triggers assembly and solve for the selected type"
      ]
    ],
    "syntax": [
      {
        "title": "Static loadcase",
        "code": "*LOADCASE, TYPE=LINEARSTATIC, NAME=<name>\n  *SUPPORTS\n  <support-collector>\n  *LOADS\n  <load-collector>\n  *SOLVER, DEVICE=CPU, METHOD=DIRECT",
        "description": "Runs a linear static analysis using selected collectors."
      }
    ],
    "examples": [
      {
        "title": "Minimal linear static loadcase",
        "code": "*LOADCASE, TYPE=LINEARSTATIC, NAME=static\n  *SUPPORTS\n  BCS\n  *LOADS\n  LOADS\n  *SOLVER, DEVICE=CPU, METHOD=DIRECT\n  *CONSTRAINTMETHOD, TYPE=NULLSPACE",
        "description": "The loadcase activates support and load collectors, configures the solver, and then executes the analysis."
      }
    ],
    "keywordsManual": {
      "required": [
        {
          "name": "TYPE",
          "values": "analysis type",
          "default": "\u2014",
          "description": "Selects the analysis implementation."
        }
      ],
      "optional": [
        {
          "name": "NAME",
          "values": "loadcase name",
          "default": "analysis type",
          "description": "Human-readable result/loadcase name."
        }
      ],
      "groups": []
    },
    "datalinesManual": [],
    "semantics": [
      "A loadcase opens an analysis context. Commands inside it select collectors and configure the solve.",
      "When the block is complete, FEMaster assembles active model data for the selected analysis type and writes result fields."
    ],
    "validation": [
      "TYPE must name a supported loadcase type.",
      "Referenced collectors and fields must exist.",
      "Solver and analysis-specific commands must be compatible with the selected TYPE."
    ],
    "mistakes": [
      "Defining load collectors but not listing them under *LOADS.",
      "Using transient controls in a static loadcase.",
      "Forgetting supports in a constrained static model."
    ],
    "related": [
      {
        "command": "LOADS",
        "description": "Activates load collectors."
      },
      {
        "command": "SUPPORTS",
        "description": "Activates support collectors."
      },
      {
        "command": "SOLVER",
        "description": "Selects solution strategy."
      },
      {
        "command": "CONSTRAINTMETHOD",
        "description": "Selects constraint handling for supported static loadcases."
      }
    ]
  },
  "LOADS": {
    "title": "Active load collectors",
    "summary": "Activate load collectors for the current loadcase.",
    "category": "Analysis",
    "badges": [
      "Loadcase command",
      "Collector activation"
    ],
    "quickFacts": [
      [
        "Valid in",
        "*LOADCASE"
      ],
      [
        "Creates",
        "No new loads"
      ],
      [
        "Modifies",
        "Active load collector list"
      ],
      [
        "Targets",
        "Load collector names"
      ],
      [
        "Internal effect",
        "Controls which load entries are assembled"
      ]
    ],
    "syntax": [
      {
        "title": "Activate collectors",
        "code": "*LOADS\n<collector-1>, <collector-2>",
        "description": "Lists one or more load collectors to assemble in the current loadcase."
      }
    ],
    "examples": [
      {
        "title": "Activate one load collector",
        "code": "*LOADCASE, TYPE=LINEARSTATIC\n  *LOADS\n  FORCE",
        "description": "Only loads stored in FORCE participate in this loadcase."
      }
    ],
    "keywordsManual": {
      "required": [],
      "optional": [],
      "groups": []
    },
    "datalinesManual": [
      {
        "title": "Collector list",
        "signature": "<load_collector_1>, <load_collector_2>, ...",
        "rows": "one or more",
        "fields": [
          {
            "name": "load_collector_i",
            "meaning": "Load collector name",
            "type": "string",
            "notes": "Collector must have been created by load commands."
          }
        ]
      }
    ],
    "semantics": [
      "FEMaster reads the listed collector names into the active loadcase. During assembly, only these collectors contribute external loads."
    ],
    "validation": [
      "Each listed load collector must exist.",
      "At least one collector should be listed when the analysis needs external loading."
    ],
    "mistakes": [
      "Assuming all defined loads are active automatically.",
      "Misspelling a collector name."
    ],
    "related": [
      {
        "command": "CLOAD",
        "description": "Creates nodal load collectors."
      },
      {
        "command": "DLOAD",
        "description": "Creates distributed load entries."
      },
      {
        "command": "LOADCASE",
        "description": "Provides the active analysis context."
      }
    ]
  },
  "SUPPORTS": {
    "title": "Active support collectors",
    "summary": "Activate support collectors for the current loadcase.",
    "category": "Analysis",
    "badges": [
      "Loadcase command",
      "Boundary activation"
    ],
    "quickFacts": [
      [
        "Valid in",
        "*LOADCASE"
      ],
      [
        "Creates",
        "No new supports"
      ],
      [
        "Modifies",
        "Active support collector list"
      ],
      [
        "Targets",
        "Support collector names"
      ],
      [
        "Internal effect",
        "Controls which prescribed DOFs are constrained"
      ]
    ],
    "syntax": [
      {
        "title": "Activate support collectors",
        "code": "*SUPPORTS\n<support-collector-1>, <support-collector-2>",
        "description": "Lists support collectors that constrain the current loadcase."
      }
    ],
    "examples": [
      {
        "title": "Activate boundary conditions",
        "code": "*LOADCASE, TYPE=LINEARSTATIC\n  *SUPPORTS\n  BCS",
        "description": "Only supports stored in BCS participate in this loadcase."
      }
    ],
    "keywordsManual": {
      "required": [],
      "optional": [],
      "groups": []
    },
    "datalinesManual": [
      {
        "title": "Collector list",
        "signature": "<support_collector_1>, <support_collector_2>, ...",
        "rows": "one or more",
        "fields": [
          {
            "name": "support_collector_i",
            "meaning": "Support collector name",
            "type": "string",
            "notes": "Collector must have been created by *SUPPORT."
          }
        ]
      }
    ],
    "semantics": [
      "The command selects which support collectors become active constraints in the loadcase."
    ],
    "validation": [
      "Each listed support collector must exist.",
      "Static models generally need enough active supports or other constraints to remove rigid body motion."
    ],
    "mistakes": [
      "Defining *SUPPORT entries but never activating them.",
      "Activating a wrong collector name and leaving the model unconstrained."
    ],
    "related": [
      {
        "command": "SUPPORT",
        "description": "Creates support collector entries."
      },
      {
        "command": "LOADCASE",
        "description": "Provides the active analysis context."
      },
      {
        "command": "CONSTRAINTMETHOD",
        "description": "Controls constraint handling."
      }
    ]
  },
  "COUPLING": {
    "title": "Coupling constraint",
    "summary": "Create coupling equations between a master reference node and a slave region.",
    "category": "Constraints",
    "badges": [
      "Constraint",
      "Root command",
      "Creates equations",
      "Node or surface slave"
    ],
    "quickFacts": [
      [
        "Valid in",
        "Root level"
      ],
      [
        "Creates",
        "Constraint equations"
      ],
      [
        "Master",
        "Reference node set"
      ],
      [
        "Slave",
        "Node set or surface set"
      ],
      [
        "Types",
        "KINEMATIC, STRUCTURAL"
      ],
      [
        "DOFs",
        "Ux, Uy, Uz, Rx, Ry, Rz"
      ]
    ],
    "syntax": [
      {
        "title": "Node-set slave",
        "code": "*COUPLING, TYPE=KINEMATIC, MASTER=<master-nset>, SLAVE=<slave-nset>\n<Ux>, <Uy>, <Uz>, <Rx>, <Ry>, <Rz>",
        "description": "Couples a slave node set to a master reference node set."
      },
      {
        "title": "Surface slave",
        "code": "*COUPLING, TYPE=STRUCTURAL, MASTER=<master-nset>, SFSET=<slave-surface-set>\n<Ux>, <Uy>, <Uz>, <Rx>, <Ry>, <Rz>",
        "description": "Uses a surface set as the slave region."
      }
    ],
    "examples": [
      {
        "title": "Kinematic coupling to reference node",
        "code": "*NSET, NAME=REF_NODE\n1000\n\n*NSET, NAME=SLAVE_NODES\n1, 2, 3, 4\n\n*COUPLING, TYPE=KINEMATIC, MASTER=REF_NODE, SLAVE=SLAVE_NODES\n1, 1, 1, 1, 1, 1",
        "description": "All translational and rotational coupling flags are enabled between the reference node and slave nodes."
      }
    ],
    "keywordsManual": {
      "required": [
        {
          "name": "TYPE",
          "values": "KINEMATIC, STRUCTURAL",
          "default": "\u2014",
          "description": "Selects the coupling formulation."
        },
        {
          "name": "MASTER",
          "values": "node set name",
          "default": "\u2014",
          "description": "Master/reference node set, typically containing one node."
        }
      ],
      "optional": [],
      "groups": [
        {
          "title": "Slave selection",
          "rows": [
            {
              "name": "SLAVE",
              "relation": "exactly one of SLAVE or SFSET",
              "description": "Slave node set."
            },
            {
              "name": "SFSET",
              "relation": "exactly one of SLAVE or SFSET",
              "description": "Slave surface set."
            }
          ]
        }
      ]
    },
    "datalinesManual": [
      {
        "title": "DOF flags",
        "signature": "<Ux>, <Uy>, <Uz>, <Rx>, <Ry>, <Rz>",
        "rows": "exactly one",
        "fields": [
          {
            "name": "Ux",
            "meaning": "Translation x",
            "type": "0/1",
            "notes": "Couple x-translation if active."
          },
          {
            "name": "Uy",
            "meaning": "Translation y",
            "type": "0/1",
            "notes": "Couple y-translation if active."
          },
          {
            "name": "Uz",
            "meaning": "Translation z",
            "type": "0/1",
            "notes": "Couple z-translation if active."
          },
          {
            "name": "Rx",
            "meaning": "Rotation x",
            "type": "0/1",
            "notes": "Couple rotation around x if meaningful."
          },
          {
            "name": "Ry",
            "meaning": "Rotation y",
            "type": "0/1",
            "notes": "Couple rotation around y if meaningful."
          },
          {
            "name": "Rz",
            "meaning": "Rotation z",
            "type": "0/1",
            "notes": "Couple rotation around z if meaningful."
          }
        ]
      }
    ],
    "semantics": [
      "FEMaster generates constraint equations between the master node and slave region.",
      "KINEMATIC couplings make slave motion follow rigid-body motion of the master. STRUCTURAL couplings distribute generalized master forces and moments to the slave region in a work-equivalent manner."
    ],
    "validation": [
      "MASTER must reference an existing node set.",
      "MASTER should normally contain exactly one reference node.",
      "Exactly one of SLAVE or SFSET must be provided.",
      "TYPE must be KINEMATIC or STRUCTURAL.",
      "The dataline must contain six DOF flags and at least one should be active."
    ],
    "mistakes": [
      "MASTER contains multiple nodes although one reference node was intended.",
      "Both SLAVE and SFSET are provided.",
      "Rotational DOFs are enabled where they are not meaningful.",
      "The slave region is geometrically degenerate."
    ],
    "related": [
      {
        "command": "NSET",
        "description": "Defines master and slave node sets."
      },
      {
        "command": "SFSET",
        "description": "Defines slave surface regions."
      },
      {
        "command": "TIE",
        "description": "Alternative for interpolation-based binding."
      },
      {
        "command": "CONTACT",
        "description": "Alternative for unilateral interaction."
      },
      {
        "command": "CONSTRAINTMETHOD",
        "description": "Controls how constraints are handled by the solver."
      }
    ]
  },
  "TIE": {
    "title": "Tie constraint",
    "summary": "Bind slave nodes or surfaces to a master surface or line.",
    "category": "Constraints",
    "badges": [
      "Constraint",
      "Root command",
      "Creates equations",
      "Surface search"
    ],
    "quickFacts": [
      [
        "Valid in",
        "Root level"
      ],
      [
        "Creates",
        "Interpolation constraint equations"
      ],
      [
        "Targets",
        "Master surface/line and slave region"
      ],
      [
        "Internal effect",
        "Adds tie equations to the constraint system"
      ]
    ],
    "syntax": [
      {
        "title": "Tie slave nodes to a master surface",
        "code": "*TIE, MASTER=MASTER_SURF, SLAVE=SLAVE_NODES, DISTANCE=0.05\n",
        "description": "Searches master geometry for each slave node within DISTANCE."
      }
    ],
    "examples": [
      {
        "title": "Minimal tie",
        "code": "*TIE, MASTER=MASTER_SURF, SLAVE=SLAVE_NODES, DISTANCE=0.05",
        "description": "The slave nodes are tied to the closest matching master geometry."
      }
    ],
    "semantics": [
      "FEMaster searches master geometry for slave nodes and creates interpolation equations for found candidates.",
      "ADJUST can move slave coordinates to the projected master location before equations are assembled."
    ],
    "validation": [
      "MASTER and SLAVE must exist.",
      "DISTANCE must be provided.",
      "ADJUST must be YES or NO if present."
    ],
    "mistakes": [
      "Using too small a DISTANCE and silently missing slave nodes.",
      "Tying surfaces with unexpected normals or side selections."
    ],
    "related": [
      {
        "command": "SURFACE",
        "description": "Defines master geometry."
      },
      {
        "command": "NSET",
        "description": "Defines slave node sets."
      },
      {
        "command": "COUPLING",
        "description": "Alternative reference-node constraint."
      },
      {
        "command": "CONTACT",
        "description": "Alternative unilateral interaction."
      }
    ]
  },
  "CONTACT": {
    "title": "Penalty contact",
    "summary": "Define frictionless node-to-surface contact for nonlinear static analysis.",
    "category": "Constraints",
    "badges": [
      "Constraint",
      "Root command",
      "Nonlinear only",
      "Penalty method"
    ],
    "quickFacts": [
      [
        "Valid in",
        "Root level"
      ],
      [
        "Creates",
        "Contact interaction definition"
      ],
      [
        "Assembled by",
        "NONLINEARSTATIC"
      ],
      [
        "Targets",
        "Master surface set and slave node/surface set"
      ],
      [
        "Internal effect",
        "Adds nonlinear contact force and tangent terms"
      ]
    ],
    "syntax": [
      {
        "title": "Frictionless contact",
        "code": "*CONTACT, MASTER=MASTER_SURF, SLAVE=SLAVE_NODES, DISTANCE=0.1, PENALTY=1e6\n",
        "description": "Defines contact candidates and normal penalty stiffness."
      }
    ],
    "examples": [
      {
        "title": "Minimal contact",
        "code": "*CONTACT, MASTER=TOOL_SURF, SLAVE=PART_NODES, DISTANCE=0.1, PENALTY=1e6",
        "description": "Contact is evaluated only by a nonlinear static loadcase."
      }
    ],
    "semantics": [
      "Contact is not a linear constraint equation. FEMaster evaluates gap and penalty force during nonlinear static iterations.",
      "MASTER normals control the sign of the normal gap; FLIP can reverse them."
    ],
    "validation": [
      "MASTER must be a surface set.",
      "SLAVE must reference a node set or surface set.",
      "DISTANCE and PENALTY must be provided.",
      "Contact is only assembled by NONLINEARSTATIC."
    ],
    "mistakes": [
      "Trying to use contact in linear static analysis.",
      "Penalty stiffness too small or too large.",
      "Unexpected normals causing contact to open instead of close."
    ],
    "related": [
      {
        "command": "SURFACE",
        "description": "Defines master surfaces."
      },
      {
        "command": "NONLINEAR",
        "description": "Controls nonlinear iterations."
      },
      {
        "command": "LOADCASE",
        "description": "Selects NONLINEARSTATIC."
      },
      {
        "command": "TIE",
        "description": "Use for permanent binding instead of contact."
      }
    ]
  },
  "MATERIAL": {
    "title": "Material context",
    "summary": "Create or select a material and attach material subcommands to it.",
    "category": "Properties",
    "badges": [
      "Material",
      "Root command",
      "Context command"
    ],
    "quickFacts": [
      [
        "Valid in",
        "Root level"
      ],
      [
        "Creates",
        "Named material context"
      ],
      [
        "Modified by",
        "*ELASTIC, *DENSITY, *THERMALEXPANSION"
      ],
      [
        "Used by",
        "Section commands"
      ]
    ],
    "syntax": [
      {
        "title": "Material with isotropic elasticity",
        "code": "*MATERIAL, NAME=STEEL\n  *ELASTIC, TYPE=ISOTROPIC\n  210000, 0.3\n  *DENSITY\n  7.85e-9",
        "description": "Material subcommands apply to the active material."
      }
    ],
    "examples": [
      {
        "title": "Minimal material",
        "code": "*MATERIAL, NAME=STEEL\n  *ELASTIC, TYPE=ISOTROPIC\n  210000, 0.3",
        "description": "Defines a material that can be referenced by a section."
      }
    ],
    "semantics": [
      "FEMaster stores material properties under the active material name. Sections later bind those properties to element sets."
    ],
    "validation": [
      "NAME must be provided.",
      "Material names referenced by sections must exist.",
      "Required material subproperties depend on the analysis and element type."
    ],
    "mistakes": [
      "Defining a material but never assigning it through a section.",
      "Mixing inconsistent units for stiffness and density."
    ],
    "related": [
      {
        "command": "ELASTIC",
        "description": "Defines stiffness."
      },
      {
        "command": "DENSITY",
        "description": "Defines mass density."
      },
      {
        "command": "SOLIDSECTION",
        "description": "Assigns material to solid elements."
      },
      {
        "command": "SHELLSECTION",
        "description": "Assigns material to shell elements."
      }
    ]
  },
  "ELASTIC": {
    "title": "Elastic material law",
    "summary": "Assign elastic stiffness data to the active material.",
    "category": "Properties",
    "badges": [
      "Material",
      "Material subcommand",
      "Variant by TYPE"
    ],
    "quickFacts": [
      [
        "Valid in",
        "*MATERIAL"
      ],
      [
        "Creates",
        "Elastic law data"
      ],
      [
        "Depends on",
        "Active material context"
      ],
      [
        "Variants",
        "Isotropic, generalized isotropic, orthotropic, ABD"
      ]
    ],
    "syntax": [
      {
        "title": "Isotropic elasticity",
        "code": "*ELASTIC, TYPE=ISOTROPIC\n<E>, <nu>",
        "description": "Defines Young's modulus and Poisson ratio for the active material."
      }
    ],
    "examples": [
      {
        "title": "Steel-like isotropic material",
        "code": "*MATERIAL, NAME=STEEL\n  *ELASTIC, TYPE=ISOTROPIC\n  210000, 0.3",
        "description": "The elastic data belongs to STEEL."
      }
    ],
    "semantics": [
      "FEMaster attaches the elastic model to the currently active material. The TYPE keyword selects the expected dataline and material formulation."
    ],
    "validation": [
      "Must appear inside a MATERIAL context.",
      "TYPE must be supported.",
      "Dataline values must match the selected TYPE."
    ],
    "mistakes": [
      "Placing *ELASTIC before *MATERIAL.",
      "Supplying isotropic two-value data while TYPE requests orthotropic data."
    ],
    "related": [
      {
        "command": "MATERIAL",
        "description": "Provides the active material context."
      },
      {
        "command": "DENSITY",
        "description": "Adds mass properties."
      },
      {
        "command": "SHELLSECTION",
        "description": "Uses material stiffness in shells."
      }
    ]
  },
  "DENSITY": {
    "title": "Material density",
    "summary": "Assign mass density to the active material.",
    "category": "Properties",
    "badges": [
      "Material",
      "Material subcommand",
      "Mass property"
    ],
    "quickFacts": [
      [
        "Valid in",
        "*MATERIAL"
      ],
      [
        "Creates",
        "Density value on active material"
      ],
      [
        "Used by",
        "Mass matrix, modal, transient, inertial workflows"
      ]
    ],
    "syntax": [
      {
        "title": "Density row",
        "code": "*DENSITY\n<rho>",
        "description": "Assigns density to the active material."
      }
    ],
    "examples": [
      {
        "title": "Material with density",
        "code": "*MATERIAL, NAME=STEEL\n  *DENSITY\n  7.85e-9",
        "description": "Density units must match the model unit system."
      }
    ],
    "semantics": [
      "FEMaster stores density on the active material and uses it when assembling mass-related operators."
    ],
    "validation": [
      "Must appear inside a MATERIAL context.",
      "Dataline must contain exactly one numeric density."
    ],
    "mistakes": [
      "Using density units inconsistent with force and length units.",
      "Expecting density to matter in a stiffness-only static check."
    ],
    "related": [
      {
        "command": "MATERIAL",
        "description": "Provides the active material context."
      },
      {
        "command": "DAMPING",
        "description": "Uses mass-proportional damping in transient analysis."
      },
      {
        "command": "INERTIARELIEF",
        "description": "May depend on mass properties."
      }
    ]
  },
  "NSET": {
    "title": "Node set",
    "summary": "Create a named reusable set of nodes.",
    "category": "Model",
    "badges": [
      "Model",
      "Root command",
      "Set definition"
    ],
    "quickFacts": [
      [
        "Valid in",
        "Root level"
      ],
      [
        "Creates",
        "Node set"
      ],
      [
        "Targets",
        "Node ids"
      ],
      [
        "Used by",
        "Supports, nodal loads, couplings, point masses"
      ]
    ],
    "syntax": [
      {
        "title": "Explicit node set",
        "code": "*NSET, NAME=FIXED\n1, 2, 3, 4",
        "description": "Adds listed nodes to FIXED."
      }
    ],
    "examples": [
      {
        "title": "Boundary node set",
        "code": "*NSET, NAME=FIXED\n1, 2, 3, 4",
        "description": "FIXED can now be referenced by supports or loads."
      }
    ],
    "semantics": [
      "FEMaster stores node set membership by name. Later commands resolve the name to the contained node ids."
    ],
    "validation": [
      "The set name must be provided.",
      "Referenced node ids should exist before use.",
      "GENERATE rows must contain valid integer ranges."
    ],
    "mistakes": [
      "Creating a set with the wrong NAME keyword spelling.",
      "Using a set before the corresponding nodes are defined or imported."
    ],
    "related": [
      {
        "command": "NODE",
        "description": "Defines node ids."
      },
      {
        "command": "SUPPORT",
        "description": "Uses node sets as targets."
      },
      {
        "command": "CLOAD",
        "description": "Uses node sets as targets."
      },
      {
        "command": "COUPLING",
        "description": "Uses node sets for master/slave regions."
      }
    ]
  },
  "ELSET": {
    "title": "Element set",
    "summary": "Create a named reusable set of elements.",
    "category": "Model",
    "badges": [
      "Model",
      "Root command",
      "Set definition"
    ],
    "quickFacts": [
      [
        "Valid in",
        "Root level"
      ],
      [
        "Creates",
        "Element set"
      ],
      [
        "Targets",
        "Element ids"
      ],
      [
        "Used by",
        "Sections, volume loads, topology controls"
      ]
    ],
    "syntax": [
      {
        "title": "Explicit element set",
        "code": "*ELSET, NAME=SOLIDS\n1, 2, 3, 4",
        "description": "Adds listed elements to SOLIDS."
      }
    ],
    "examples": [
      {
        "title": "Solid region",
        "code": "*ELSET, NAME=SOLIDS\n1, 2, 3, 4\n\n*SOLIDSECTION, ELSET=SOLIDS, MATERIAL=STEEL",
        "description": "The element set is later used by a section."
      }
    ],
    "semantics": [
      "FEMaster stores element membership by name. Sections and element-based loads resolve the name during model setup or assembly."
    ],
    "validation": [
      "The set name must be provided.",
      "Referenced element ids should exist before use.",
      "GENERATE ranges must be valid."
    ],
    "mistakes": [
      "Assigning a shell section to an element set that contains solids.",
      "Forgetting that sections target element sets, not individual nodes."
    ],
    "related": [
      {
        "command": "ELEMENT",
        "description": "Defines element ids."
      },
      {
        "command": "SOLIDSECTION",
        "description": "Assigns material to solid sets."
      },
      {
        "command": "VLOAD",
        "description": "Targets element sets."
      },
      {
        "command": "SURFACE",
        "description": "Can generate surfaces from element sets."
      }
    ]
  },
  "SURFACE": {
    "title": "Surface definition",
    "summary": "Create named element-side surface entries for loads and interactions.",
    "category": "Model",
    "badges": [
      "Model",
      "Root command",
      "Surface target"
    ],
    "quickFacts": [
      [
        "Valid in",
        "Root level"
      ],
      [
        "Creates",
        "Surface entries in a surface set"
      ],
      [
        "Targets",
        "Element faces or element sets"
      ],
      [
        "Used by",
        "DLOAD, PLOAD, TIE, CONTACT, COUPLING"
      ]
    ],
    "syntax": [
      {
        "title": "Element face surface",
        "code": "*SURFACE, NAME=SURF_TOP\n1, 10, S3",
        "description": "Creates surface id 1 from side S3 of element 10."
      }
    ],
    "examples": [
      {
        "title": "Surface for pressure",
        "code": "*SURFACE, NAME=SURF_TOP, TYPE=ELEMENT\nTOP_ELEMENTS, S3",
        "description": "Generates one surface entry for side S3 of each element in TOP_ELEMENTS."
      }
    ],
    "semantics": [
      "FEMaster stores surface entries as references to element sides. Later commands use those entries as geometric integration or search targets."
    ],
    "validation": [
      "The referenced element or element set must exist.",
      "SIDE must be valid for the element topology.",
      "Surface orientation must be considered for pressure and contact."
    ],
    "mistakes": [
      "Using the wrong side id and applying loads on the opposite face.",
      "Assuming pressure sign is independent of surface normal."
    ],
    "related": [
      {
        "command": "ELEMENT",
        "description": "Defines element topology and side numbering."
      },
      {
        "command": "SFSET",
        "description": "Groups surface ids."
      },
      {
        "command": "DLOAD",
        "description": "Applies vector traction."
      },
      {
        "command": "PLOAD",
        "description": "Applies pressure."
      },
      {
        "command": "CONTACT",
        "description": "Uses master surfaces."
      }
    ]
  },
  "SFSET": {
    "title": "Surface set",
    "summary": "Create a named reusable set of surface ids.",
    "category": "Model",
    "badges": [
      "Model",
      "Root command",
      "Set definition"
    ],
    "quickFacts": [
      [
        "Valid in",
        "Root level"
      ],
      [
        "Creates",
        "Surface set"
      ],
      [
        "Targets",
        "Surface ids"
      ],
      [
        "Used by",
        "Surface loads, tie, contact, surface coupling"
      ]
    ],
    "syntax": [
      {
        "title": "Explicit surface set",
        "code": "*SFSET, NAME=LOAD_SURF\n1, 2, 3, 4",
        "description": "Groups surface ids for reuse."
      }
    ],
    "examples": [
      {
        "title": "Reusable load surface",
        "code": "*SFSET, NAME=LOAD_SURF\n1, 2, 3, 4\n\n*PLOAD, LOAD_COLLECTOR=P\nLOAD_SURF, 1.0",
        "description": "The pressure can target the surface set name."
      }
    ],
    "semantics": [
      "FEMaster stores surface id membership by name. Loads and interactions resolve the set to the contained surface entries."
    ],
    "validation": [
      "The surface set name should be provided.",
      "Referenced surface ids should exist.",
      "GENERATE rows must contain valid ranges."
    ],
    "mistakes": [
      "Confusing *SURFACE, which creates surface entries, with *SFSET, which groups existing surface ids."
    ],
    "related": [
      {
        "command": "SURFACE",
        "description": "Creates surface ids."
      },
      {
        "command": "DLOAD",
        "description": "Targets surface sets."
      },
      {
        "command": "PLOAD",
        "description": "Targets surface sets."
      },
      {
        "command": "TIE",
        "description": "Uses surface sets."
      },
      {
        "command": "CONTACT",
        "description": "Uses surface sets."
      }
    ]
  },
  "ELEMENT": {
    "title": "Finite element definition",
    "summary": "Create elements by assigning topology and node connectivity.",
    "category": "Model",
    "badges": [
      "Model",
      "Root command",
      "Topology",
      "Has element library"
    ],
    "quickFacts": [
      [
        "Valid in",
        "Root level"
      ],
      [
        "Creates",
        "Finite elements"
      ],
      [
        "Depends on",
        "Node ids"
      ],
      [
        "Targets",
        "Optional ELSET"
      ],
      [
        "Datalines",
        "One element id plus topology-specific connectivity"
      ]
    ],
    "syntax": [
      {
        "title": "Hex element",
        "code": "*ELEMENT, TYPE=C3D8, ELSET=SOLIDS\n1, 1, 2, 3, 4, 5, 6, 7, 8",
        "description": "Creates one C3D8 element and adds it to SOLIDS."
      }
    ],
    "examples": [
      {
        "title": "Minimal element row",
        "code": "*ELEMENT, TYPE=T33, ELSET=BARS\n1, 1, 2",
        "description": "The TYPE keyword determines how many connectivity entries are expected."
      }
    ],
    "semantics": [
      "FEMaster creates element objects from the selected TYPE and connectivity row. The element type controls interpolation, side numbering, active DOFs, integration, and compatible sections."
    ],
    "validation": [
      "TYPE must be provided.",
      "Node ids in connectivity must exist.",
      "The dataline length must match the selected TYPE.",
      "ELSET defaults to EALL if omitted."
    ],
    "mistakes": [
      "Using the wrong node ordering and getting flipped normals or bad element orientation.",
      "Assigning a section incompatible with the element family.",
      "Forgetting to define the referenced nodes first."
    ],
    "related": [
      {
        "command": "NODE",
        "description": "Defines node coordinates."
      },
      {
        "command": "ELSET",
        "description": "Groups elements."
      },
      {
        "command": "SURFACE",
        "description": "References element sides."
      },
      {
        "command": "SOLIDSECTION",
        "description": "Assigns material to solid elements."
      },
      {
        "command": "SHELLSECTION",
        "description": "Assigns material and thickness to shells."
      }
    ]
  },
  "NODE": {
    "title": "Node definition",
    "summary": "Create nodal coordinates and optional node-set membership.",
    "category": "Model",
    "badges": [
      "Model",
      "Root command",
      "Geometry"
    ],
    "quickFacts": [
      [
        "Valid in",
        "Root level"
      ],
      [
        "Creates",
        "Nodes"
      ],
      [
        "Targets",
        "Optional NSET"
      ],
      [
        "Datalines",
        "Node id plus x, y, z coordinates"
      ],
      [
        "Internal effect",
        "Coordinates are referenced by elements and geometry operations"
      ]
    ],
    "syntax": [
      {
        "title": "Node coordinates",
        "code": "*NODE, NSET=NALL\n1, 0.0, 0.0, 0.0\n2, 1.0, 0.0, 0.0",
        "description": "Creates two nodes and adds them to NALL."
      }
    ],
    "examples": [
      {
        "title": "Two nodes",
        "code": "*NODE, NSET=NALL\n1, 0.0, 0.0, 0.0\n2, 1.0, 0.0, 0.0",
        "description": "Nodes become useful when elements, loads, supports, or constraints reference them."
      }
    ],
    "semantics": [
      "FEMaster stores nodes as coordinate rows. Active DOFs are not created by nodes alone; they arise from connected model objects and constraints."
    ],
    "validation": [
      "Node ids must be integers.",
      "Each dataline should provide id, x, y, z.",
      "Duplicate node ids should be avoided."
    ],
    "mistakes": [
      "Expecting every node to have all six DOFs automatically.",
      "Using inconsistent coordinate units."
    ],
    "related": [
      {
        "command": "ELEMENT",
        "description": "Connects nodes into elements."
      },
      {
        "command": "NSET",
        "description": "Groups nodes."
      },
      {
        "command": "SUPPORT",
        "description": "Constrains node components."
      },
      {
        "command": "CLOAD",
        "description": "Applies nodal loads."
      }
    ]
  },
  "SOLVER": {
    "title": "Solver selection",
    "summary": "Choose solver device and method for the active loadcase.",
    "category": "Analysis",
    "badges": [
      "Loadcase command",
      "Solver control"
    ],
    "quickFacts": [
      [
        "Valid in",
        "*LOADCASE"
      ],
      [
        "Creates",
        "Solver configuration"
      ],
      [
        "Depends on",
        "Analysis type and constraint method"
      ],
      [
        "Options",
        "DEVICE, METHOD"
      ]
    ],
    "syntax": [
      {
        "title": "CPU direct solver",
        "code": "*SOLVER, DEVICE=CPU, METHOD=DIRECT",
        "description": "Selects direct CPU solve for the active loadcase."
      }
    ],
    "examples": [
      {
        "title": "Static solve setup",
        "code": "*LOADCASE, TYPE=LINEARSTATIC\n  *SOLVER, DEVICE=CPU, METHOD=DIRECT\n  *CONSTRAINTMETHOD, TYPE=NULLSPACE",
        "description": "Solver settings apply only inside the current loadcase."
      }
    ],
    "semantics": [
      "FEMaster stores solver options on the active loadcase. Compatibility is checked when the analysis starts."
    ],
    "validation": [
      "SOLVER must appear inside a loadcase.",
      "DEVICE and METHOD values must be supported.",
      "Some constraint methods require direct solving."
    ],
    "mistakes": [
      "Selecting LAGRANGE constraints with an incompatible indirect solver.",
      "Assuming GPU is available in a CPU-only build."
    ],
    "related": [
      {
        "command": "LOADCASE",
        "description": "Provides the analysis context."
      },
      {
        "command": "CONSTRAINTMETHOD",
        "description": "Controls constraint handling."
      },
      {
        "command": "NONLINEAR",
        "description": "Controls nonlinear static iterations."
      }
    ]
  },
  "ORIENTATION": {
    "title": "Local coordinate system",
    "summary": "Define a named rectangular or cylindrical coordinate system.",
    "category": "Properties",
    "badges": [
      "Model data",
      "Root command",
      "Local axes"
    ],
    "quickFacts": [
      [
        "Valid in",
        "Root level"
      ],
      [
        "Creates",
        "Named coordinate system"
      ],
      [
        "Used by",
        "Loads, supports, sections, connectors"
      ],
      [
        "Types",
        "RECTANGULAR, CYLINDRICAL"
      ]
    ],
    "syntax": [
      {
        "title": "Rectangular orientation",
        "code": "*ORIENTATION, NAME=LOCAL_XY, TYPE=RECTANGULAR\n1, 0, 0\n0, 1, 0",
        "description": "Defines local axes from direction vectors."
      }
    ],
    "examples": [
      {
        "title": "Local basis for loads",
        "code": "*ORIENTATION, NAME=LOCAL_LOAD, TYPE=RECTANGULAR\n1, 0, 0\n0, 1, 0",
        "description": "Other commands can reference LOCAL_LOAD."
      }
    ],
    "semantics": [
      "FEMaster stores the coordinate system and uses it to transform local components into global coordinates when referenced."
    ],
    "validation": [
      "NAME and TYPE must be provided.",
      "The vector data must define a usable basis.",
      "Referenced orientations must exist before use."
    ],
    "mistakes": [
      "Using nearly collinear direction vectors.",
      "Forgetting that local load/support components are transformed to global coordinates."
    ],
    "related": [
      {
        "command": "DLOAD",
        "description": "Can use local traction components."
      },
      {
        "command": "CLOAD",
        "description": "Can use local force components."
      },
      {
        "command": "SUPPORT",
        "description": "Can use local prescribed components."
      },
      {
        "command": "BEAMSECTION",
        "description": "Uses local section orientation."
      }
    ]
  },
  "AMPLITUDE": {
    "title": "Amplitude curve",
    "summary": "Define a reusable scalar time-history for load scaling.",
    "category": "Loads and supports",
    "badges": [
      "Root command",
      "Time function",
      "Load scaling"
    ],
    "quickFacts": [
      [
        "Valid in",
        "Root level"
      ],
      [
        "Creates",
        "Named amplitude curve"
      ],
      [
        "Used by",
        "CLOAD, DLOAD, PLOAD, VLOAD"
      ],
      [
        "Datalines",
        "time, value pairs"
      ]
    ],
    "syntax": [
      {
        "title": "Linear amplitude",
        "code": "*AMPLITUDE, NAME=PULSE, TYPE=LINEAR\n0.0, 0.0\n1.0, 1.0",
        "description": "Defines a piecewise-linear scale factor."
      }
    ],
    "examples": [
      {
        "title": "Load scaling curve",
        "code": "*AMPLITUDE, NAME=RAMP\n0.0, 0.0\n1.0, 1.0",
        "description": "Loads can reference RAMP by name."
      }
    ],
    "semantics": [
      "FEMaster stores ordered time-value pairs. During transient analysis, referenced loads are scaled by the interpolated amplitude value."
    ],
    "validation": [
      "NAME must be provided.",
      "TYPE must be STEP, NEAREST, or LINEAR if present.",
      "Each dataline must contain time and value."
    ],
    "mistakes": [
      "Referencing an amplitude name that was never defined.",
      "Using unsorted or duplicate times in a way that makes interpolation ambiguous."
    ],
    "related": [
      {
        "command": "CLOAD",
        "description": "Can reference amplitudes."
      },
      {
        "command": "DLOAD",
        "description": "Can reference amplitudes."
      },
      {
        "command": "TIME",
        "description": "Defines transient time range."
      },
      {
        "command": "LOADS",
        "description": "Activates scaled loads."
      }
    ]
  }
};
