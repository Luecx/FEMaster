# 103 mortar parallel-block benchmark

This benchmark isolates surface-to-surface mortar contact between two matching
C3D8 faces. The lower block is fixed. The upper block has a 0.5 mm initial gap,
is prevented from moving laterally and receives a prescribed top displacement
of -1.0 mm at lambda = 1.

## Geometry and material

- contact area: A = 1 mm^2
- upper-block height: H = 100 mm
- initial gap: g0 = 0.5 mm
- prescribed top displacement: u_top(lambda) = -lambda mm
- Young's modulus: E = 1000 N/mm^2
- Poisson ratio: nu = 0
- axial stiffness of upper block: k = E A / H = 10 N/mm
- mortar penalty: epsilon = 1e5 N/mm^3

Surface-to-surface contact has no search-distance parameter. Every master
surface is tested directly during the current common-plane segmentation.

## Exact pre-contact response

Before contact the upper block can translate rigidly in z and therefore remains
stress free. Its lower face follows the prescribed top displacement exactly:

    u_bottom = -lambda mm
    g(lambda) = g0 - lambda mm = 0.5 - lambda mm

Therefore contact must remain completely inactive for lambda < 0.5.
At lambda = 0.5 the faces touch with zero pressure.

| lambda | free gap [mm] | expected contact |
|------:|--------------:|------------------|
| 0.1 | +0.4 | open |
| 0.2 | +0.3 | open |
| 0.3 | +0.2 | open |
| 0.4 | +0.1 | open |
| 0.5 |  0.0 | just touching, zero pressure |

For every one of these increments the contact force must remain zero. Projected
segments and mortar quadrature points may already exist because projected overlap
is independent of normal closure.

## Post-contact response

For lambda > 0.5 the upper block compresses axially against the fixed lower
block. With a pure penalty law and small-strain axial stiffness k, let delta_p
be the positive penetration. Equilibrium gives

    epsilon A delta_p = k (lambda - 0.5 - delta_p)

and therefore

    delta_p = k (lambda - 0.5) / (epsilon A + k)
    F       = epsilon A delta_p

With k = 10 N/mm and epsilon A = 100000 N/mm:

| lambda | imposed overtravel [mm] | penalty penetration [mm] | reaction magnitude [N] |
|------:|--------------------------:|-------------------------:|-----------------------:|
| 0.6 | 0.1 | 9.9990001e-6 | 0.9999000 |
| 0.7 | 0.2 | 1.9998000e-5 | 1.9998000 |
| 0.8 | 0.3 | 2.9997000e-5 | 2.9997000 |
| 0.9 | 0.4 | 3.9996000e-5 | 3.9996000 |
| 1.0 | 0.5 | 4.9995000e-5 | 4.9995000 |

The nonlinear solid formulation introduces only a very small correction because
the maximum axial strain is about 0.5%. The contact onset, gap sign, overlap
area and force direction are exact geometric checks independent of that small
constitutive/geometric correction.

Augmented-Lagrange updates may reduce the penetration below the pure-penalty
values. They must not change the onset of contact or create force for lambda <=
0.5.

## Mortar invariants to inspect

For this matching one-face/one-face case:

- master normal on C3D8 S2: +z
- slave normal on C3D8 S1: -z
- normal dot product: -1
- projected overlap area: 1 mm^2
- exactly one visible master layer everywhere
- no self-contact rejection
- no hidden-layer ambiguity
- all four slave mortar nodes are geometrically equivalent
- total master and slave contact forces must be equal and opposite

This case should be the first regression benchmark for every change to the
surface-to-surface mortar geometry or active-set formulation.
