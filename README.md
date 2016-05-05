# Rotations.jl

### Example Usage

```julia

    using FixedSizeArrays
    using Rotations

    # create a rotation matrix
    R = eye(RotMatrix{Float64})

    # create a point
    X = Vec(1.0, 2.0, 3.0)

    # convert to Euler Angles (using the default ordering scheme)
    ea = EulerAngles(R)
    Xo = rotate(ea, X)

    # convert to Euler Angles specifying the angle order
    ea = EulerAngles{Rotations.EulerXYZ}(R)
    Xo = rotate(ea, X)

    # convert to proper Euler Angles (using the default ordering scheme)
    ea = ProperEulerAngles(R)
    Xo = rotate(ea, X)

    # convert to proper Euler Angles specifying the angle order
    ea = ProperEulerAngles{Rotations.EulerXYX}(R)
    Xo = rotate(ea, X)

    # a Quaternion
    q = Quaternion(R)
    Xo = rotate(ea, X)

    # Stereo graphic projection of a quaternion (recommended for optimization problems)
    spq = SpQuat(R)
    Xo = rotate(spq, X)

```

This package implements various 3D rotation parameterizations and defines conversions between them.  All rotation variables are stored as as immutable types.

This package assumes [active (right handed) rotations](https://en.wikipedia.org/wiki/Active_and_passive_transformation) where applicable.

### Rotation Parameterizations

1. **Rotation Matrix** `RotMatrix{T <: AbstractFloat}`

    A 3 x 3 rotation matrix storing the rotation.  This is a typealias of a [FixedSizeArrays](https://github.com/SimonDanisch/FixedSizeArrays.jl) `Mat{3,3,T}`.  A rotation matrix `R` should have the property `I = R * R'`.



2. **Arbitrary Axis Rotation** `AngleAxis{T <: AbstractFloat}`

    A 4 element immutable array with fields `theta`, `x`, `y`, and `z` to store the rotation angle and axis of the rotation.



3. **EulerAngles** `EulerAngles{Order <: TaitByranOrder, T <: AbstractFloat}`

    A 3 element immutable array which stores the Euler angles with [**Tait Byran**](https://en.wikipedia.org/wiki/Euler_angles#Tait.E2.80.93Bryan_angles) angle ordering (e.g. `EulerXYZ`) encoded by the template parameter `Order`.



4. **ProperEulerAngles** `ProperEulerAngles{Order <: ProperEulerOrder, T <: AbstractFloat}`

    A 3 element immutable array which stores the Euler angles with [**Proper Euler**](https://en.wikipedia.org/wiki/Euler_angles#Conventions) angle ordering (e.g. `EulerXYX`) encoded by the template paramter `Order`.



5. **Quaternions** `Quaternion{T <: Real}`

    A 4 element immutable array containing the Quaternion representation of the rotation.  This uses the [Quaternions](https://github.com/JuliaGeometry/Quaternions.jl) package.



6. **Stereographic Quaternion Projection** `SpQuat{T <: AbstractFloat}`

    A 3 element immutable array containing the stereographic projection of a unit Quaternion.  This gives a compact three element representation of the rotation, but unlike EulerAngles the derivitives of the rotation matrix w.r.t. the `SpQuat` parameters are rational functions.  This makes the `SpQuat` type a good choice for use in optimization problems.

    This projection can be visualized as a pin hole camera, with the pin hole matching the Quaternion `[-1,0,0,0]` and the image plane containing the origin and having normal direction `[1,0,0,0]`.  The "no rotation" `Quaternion(1.0,0,0,0)` then maps to the `SpQuat(0,0,0)`






The `convert` function is used to convert element types for the same parameterization, e.g.

```julia

    spq_f64 = SpQuat(0.0,0.0,0.0)
    spq_32 = convert(SpQuat{Float32}, spq_f64)

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


### Why use immutables / FixedSizeArrays?

They're faster (BLAS isn't great for 3x3 matrices) and don't need preallocating:

```julia

    function benchmark(n::Int=1_000_000)

        function rotate_mutable(R, X, n)
            Xb, Xo = zeros(3), zeros(3)
            for i = 1:n
                A_mul_B!(Xb, R, X)
                # @inbounds Xo[1] += Xb[1]; @inbounds Xo[2] += Xb[2]; @inbounds Xo[3] += Xb[3];
            end
            return Xo
        end

        function rotate_immutable(R, X, n)
            Xo = Vec(0.0,0,0)
            for i = 1:n
                Xb = R * X
                # Xo += Xb
            end
            return Xo
        end

        # Initialise
        R_mute, R_immute = eye(3), eye(RotMatrix)
        X_mute, X_immute = zeros(3), Vec(0.0, 0.0, 0.0)

        # and test
        rotate_mutable(R_mute, X_mute, 1)
        println("Rotating using mutables")
        @time Xo = rotate_mutable(R_mute, X_mute, n)

        rotate_immutable(R_immute, X_immute, 1)
        println("Rotating using immutables")
        @time Xo = rotate_immutable(R_immute, X_immute, n)

    end

    benchmark()
    # Rotating using mutables
    #   0.077123 seconds (3 allocations: 208 bytes)
    #   Rotating using immutables
    #   0.003569 seconds (4 allocations: 160 bytes)

```

