# Rotations

This package implements various rotation parameterisations and defined convewrsaions between them.  The design is to store rotation variables as immutable types where allowed by dependency.

This package assumes [active (right handed) rotations](https://en.wikipedia.org/wiki/Active_and_passive_transformation) wheere applicable.

### Rotation parameterisations

1. **Rotation Matrix** `RotMatrix{T <: AbstractFloat}`

    A 3 x 3 rotation matrix stroing the rotation.  This is a typealias of a [FixedSizeArrays](https://github.com/SimonDanisch/FixedSizeArrays.jl) `Mat{3,3,T}`.  A rotation matrix `R` should have the property `I = R * R<sup>T</sup>`.


2. **Arbitrary Axis Rotation** `AngleAxis{T <: AbstractFloat}`

    A 4 element immutable array with fields `theta`, `x`, `y`, and `z` to store the rotation angle and axis of the rotation.  This is a typealias of a [FixedSizeArrays](https://github.com/SimonDanisch/FixedSizeArrays.jl) `FixedVectorNoTuple{4, T}`


3. **EulerAngles** `EulerAngles{Order, T <: AbstractFloat}`

    A 3 element immutable array which stores the Euler angles with [**Tait Byran** angle ordering](https://en.wikipedia.org/wiki/Euler_angles#Tait.E2.80.93Bryan_angles) (e.g. `EulerXYZ`) encoded by the template paramter `Order`.   This is a typealias of a [FixedSizeArrays](https://github.com/SimonDanisch/FixedSizeArrays.jl) `FixedVectorNoTuple{3, T}`


3. **EulerAngles** `EulerAngles{Order <: EulerOrder, T <: AbstractFloat}`

    A 3 element immutable array which stores the Euler angles with [**Tait Byran** angle ordering](https://en.wikipedia.org/wiki/Euler_angles#Tait.E2.80.93Bryan_angles) (e.g. `EulerXYZ`) encoded by the template paramter `Order`.   This is a typealias of a [FixedSizeArrays](https://github.com/SimonDanisch/FixedSizeArrays.jl) `FixedVectorNoTuple{3, T}`



### Rotation conversions

The `convert_rotations` function and copy constructor methodology is used to transform between different rotation parameterizations.

The `convert` function is used to convert 
