# Rotations.jl

This package implements various rotation parameterizations and defines conversions between them.  All rotation variables are stored as as immutable types.

This package assumes [active (right handed) rotations](https://en.wikipedia.org/wiki/Active_and_passive_transformation) wheere applicable.

### Rotation Parameterizations

1. **Rotation Matrix** `RotMatrix{T <: AbstractFloat}`

    A 3 x 3 rotation matrix storing the rotation.  This is a typealias of a [FixedSizeArrays](https://github.com/SimonDanisch/FixedSizeArrays.jl) `Mat{3,3,T}`.  A rotation matrix `R` should have the property `I = R * R<sup>T</sup>`.



2. **Arbitrary Axis Rotation** `AngleAxis{T <: AbstractFloat}`

    A 4 element immutable array with fields `theta`, `x`, `y`, and `z` to store the rotation angle and axis of the rotation.  This is a typealias of a [FixedSizeArrays](https://github.com/SimonDanisch/FixedSizeArrays.jl) `FixedVectorNoTuple{4, T}`.



3. **EulerAngles** `EulerAngles{Order <: TaitByranOrder, T <: AbstractFloat}`

    A 3 element immutable array which stores the Euler angles with [**Tait Byran**](https://en.wikipedia.org/wiki/Euler_angles#Tait.E2.80.93Bryan_angles) angle ordering (e.g. `EulerXYZ`) encoded by the template parameter `Order`.   This is a typealias of a [FixedSizeArrays](https://github.com/SimonDanisch/FixedSizeArrays.jl) `FixedVectorNoTuple{3, T}`.



4. **ProperEulerAngles** `ProperEulerAngles{Order <: ProperEulerOrder, T <: AbstractFloat}`

    A 3 element immutable array which stores the Euler angles with [**Proper Euler**](https://en.wikipedia.org/wiki/Euler_angles#Conventions) angle ordering (e.g. `EulerXYX`) encoded by the template paramter `Order`.   This is a typealias of a [FixedSizeArrays](https://github.com/SimonDanisch/FixedSizeArrays.jl) `FixedVectorNoTuple{3, T}`.



5. **Quaternions** `Quaternion{T <: Real}`

    A 4 element immutable array containing the Quaternion representation of the rotation.  This uses the [Quaternions](https://github.com/JuliaGeometry/Quaternions.jl) package.



6. **Stereographic Quaternion Projection** `SpQuat{T <: AbstractFloat}`

    A 3 element immutable array containing the stereographic projection of a unit Quaternion.  This gives a compact three element representation of the rotation, but unlike EulerAngles the derivitives of the rotation matrix w.r.t. the `SpQuat` parameters are rational functions.  This makes the `SpQuat` type a good choice for use in optimization problems.

    This projection can be visualized as a pin hole camera, with the pin hole matching the Quaternion `[-1,0,0,0]` and the image plane containing the origin and having normal direction `[1,0,0,0]`.  The "no rotation" `Quaternion(1.0,0,0,0)` then maps to the `SpQuat(0,0,0)`




### Rotation Conversions

The `convert_rotations` function and copy constructor are used to transform between different rotation parameterizations, e.g.

```julia

    # create a matrix
    R = RotMatrix(eye(3))

    # convert to Euler Angles (using the default ordering scheme)
    ea = convert_rotation(EulerAngles, R)

    # convert to Euler Angles specifying the order
    ea = convert_rotation(EulerAngles{Rotations.EulerXYZ}, R)

    # a Quaternion
    q = convert_rotation(Quaternion, R)

    # or equivalently
    q = Quaternion(R)

    # change the element type while we're at it
    q = convert_rotation(Quaternion{Float32}, R)


```

The `convert` function is used to convert element types for the same parameterization, e.g.

```julia

    # create a matrix
    Rf64 = RotMatrix(eye(3))

    # convert to a different element type
    Rf32 = convert(RotMatrix{Float32}, Rf64)

```



### Import / Export

All parameterizations can be converted to and from mutable / immutable vectors, e.g.

```julia

    using FixedSizeArrays

    # export
    q = Quaternion(1.0,0,0,0)
    v_mutable = Vector(q)
    v_immutable = Vec(q)
    
    # import
    q2 = Quaternion(v_mutable)
    q2 = Quaternion(v_immutable)

```

