# Rotations.jl

*3D rotations made easy in Julia*

[![Build Status](https://travis-ci.org/FugroRoames/Rotations.jl.svg?branch=static_arrays)](https://travis-ci.org/FugroRoames/Rotations.jl)

This package implements various 3D rotation parameterizations and defines
conversions between them. At their heart, each rotation parameterization *is*
a 3×3 unitary (orthogonal) matrix (based on the [StaticArrays.jl package](https://github,com/andyferris/StaticArrays.jl)),
and acts to rotate a 3-vector about the origin through matrix-vector multiplication.

While the `RotMatrix` type is a dense representation of a `3×3` matrix, we also
have sparse (or computed, rather) representations such as quaternions,
angle-axis parameterizations, and Euler angles.

For composing rotations about the origin with other transformations, see also
the [CoordinateTransformations.jl](https://github.com/FugroRoames/CoordinateTransformations.jl)
or [AffineTransforms.jl](https://github.com/timholy/AffineTransforms.jl) packages.

### Example Usage

```julia
using Rotations, StaticArrays

# create the null rotation (identity matrix)
id = eye(RotMatrix{Float64})

# create a random rotation matrix (uniformly distributed over all 3D rotations)
r = rand(RotMatrix) # uses Float64 by default

# create a point
p = SVector(1.0, 2.0, 3.0) # from StaticArrays.jl, but could use any AbstractVector...

# convert to a quaternion (Quat) and rotate the point
q = Quat(r)
p_rotated = q * p

# Compose rotations
q2 = rand(Quat)
q3 = q * q2

# Take the inverse (equivalent to transpose)
q_inv = q'
q_inv == inv(q)
p ≈ q_inv * (q * p)

# convert to a Stereographic quaternion projection (recommended for applications with differentiation)
spq = SPQuat(r)

# convert to the Angle-axis parameterization, or related Rodrigues vector
aa = AngleAxis(r)
rv = RodriguesVec(r)

# convert to Euler angles, composed of X/Y/Z axis rotations (Z applied first)
# (all combinations of "RotABC" are defined)
r_xyz = RotXYZ(r)

# Rotation about the X axis by 0.1 radians
r_x = RotX(0.1)

# Composing axis rotations together automatically results in Euler parameterization
RotX(0.1) * RotY(0.2) * RotZ(0.3) === RotXYZ(0.1, 0.2, 0.3)
```

### Rotation Parameterizations

1. **Rotation Matrix** `RotMatrix{T}`

    A 3 x 3 rotation matrix storing the rotation.  This is a simple wrapper for
    a [StaticArrays](https://github.com/andyferris/StaticArrays.jl) `SMatrix{3,3,T}`.
    A rotation matrix `R` should have the property `I = R * R'`, but this isn't
    enforced by the constructor.


2. **Arbitrary Axis Rotation** `AngleAxis{T}`

    A 3D rotation with fields `theta`, `axis_x`, `axis_y`, and
    `axis_z` to store the rotation angle and axis of the rotation.
    Like all other types in this package, once it is constructed it acts and
    behaves as a 3×3 `AbstractMatrix`. The axis will be automatically
    renormalized by the constructor to be a unit vector, so that `theta` always
    represents the rotation angle in radians.

3. **Quaternions** `Quat{T}`

    A 3D rotation parameterized by a unit quaternion. Note that the constructor
    will renormalize the quaternion to be a unit quaternion, and that although
    they follow the same multiplicative *algebra* as quaternions, it is better
    to think of `Quat` as a 3×3 matrix rather than as a quaternion *number*.

4. **Rodrigues Vector** `RodriguesVec{T}`

    A 3D rotation encoded by an angle-axis representation as `angle * axis`.  
    This type is used in packages such as [OpenCV](http://docs.opencv.org/2.4/modules/calib3d/doc/camera_calibration_and_3d_reconstruction.html#void%20Rodrigues%28InputArray%20src,%20OutputArray%20dst,%20OutputArray%20jacobian%29).

    Note: If you're differentiating a Rodrigues Vector check the result is what
    you expect at theta = 0.  The first derivative of the rotation *should*
    behave, but higher-order derivatives of it (as well as parameterization
    conversions) should be tested.  The Stereographic Quaternion Projection is
    the recommended three parameter format for differentiation.

5. **Stereographic Projection of a unit Quaternion** `SPQuat{T}`

    A 3D rotation encoded by the stereographic projection of a unit quaternion.  This projection can be visualized as a pin hole camera, with the pin hole matching the quaternion `[-1,0,0,0]` and the image plane containing the origin and having normal direction `[1,0,0,0]`.  The "null rotation" `Quaternion(1.0,0,0,0)` then maps to the `SPQuat(0,0,0)`

    These are similar to the Rodrigues vector in that the axis direction is stored in an unnormalized form, and the rotation angle is encoded in the length of the axis.  This type has the nice property that the derivatives of the rotation matrix w.r.t. the `SPQuat` parameters are rational functions, making the `SPQuat` type a good choice for differentiation / optimization.

6. **Cardinal axis rotations** `RotX{T}`, `RotY{T}`, `RotZ{T}`

    Sparse representations of 3D rotations about the X, Y, or Z axis, respectively.

7. **Two-axis rotations** `RotXY{T}`, etc

    Conceptually, these are compositions of two of the cardinal axis rotations above,
    so that `RotXY(x, y) == RotX(x) * RotY(y)` (note that the order of application to
    a vector is right-to-left, as-in matrix-matrix-vector multiplication, `Mx * My * v`).

8. **Euler Angles - Three-axis rotations** `RotXYZ{T}`, `RotXYX{T}`, etc

    A composition of 3 cardinal axis rotations is typically known as a Euler
    angle parameterization of a 3D rotation. The rotations with 3 unique axes,
    such as `RotXYZ`, are said to follow the [**Tait Byran**](https://en.wikipedia.org/wiki/Euler_angles#Tait.E2.80.93Bryan_angles) angle ordering,
    while those which repeat (e.g. `EulerXYX`) are said to use [**Proper Euler**](https://en.wikipedia.org/wiki/Euler_angles#Conventions) angle ordering.

    Like the two-angle versions, read the application of the rotations to a
    vector from right-to-left, so that `RotXYZ(x, y, z) * v = RotX(x) * RotY(x) * RotZ(z) * v`.

### Import / Export

All parameterizations can be converted to and from (mutable or immutable)
3×3 matrices, e.g.

```julia
using StaticArrays, Rotations

# export
q = Quat(1.0,0,0,0)
matrix_mutable = Array(q)
matrix_immutable = SMatrix{3,3}(q)

# import
q2 = Quaternion(matrix_mutable)
q3 = Quaternion(matrix_immutable)
```

### Notes

This package assumes [active (right handed) rotations](https://en.wikipedia.org/wiki/Active_and_passive_transformation) where applicable.


### Why use immutables / StaticArrays?

They're faster (Julia's `Array` and BLAS aren't great for 3x3 matrices) and
don't need preallocating or garbage collection. For example, see this (now
outdated) benchmark case:

```julia
cd(Pkg.dir("Rotations") * "/test")
include("benchmark.jl")
BenchMarkRotations.benchmark_mutable()
# Rotating using mutables
#   0.077123 seconds (3 allocations: 208 bytes)
#   Rotating using immutables
#   0.003569 seconds (4 allocations: 160 bytes)
```
