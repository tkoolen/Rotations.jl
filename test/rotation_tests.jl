# function to perform tests of the rotation functions in the Rotations module
using FixedSizeArrays
using Quaternions
using Rotations


# reset this for testing
srand(0)

numel(X) = (isa(X, Vector) || isa(X, Vec) || isa(X, Matrix) || isa(X, Mat)) ? length(X) : Rotations.nvars(X)

# a macro to test if the contents of two containers are approximatley equal
macro contents_approx_eq(a, b)
    quote
        n = numel($(esc(a)))
        @test n == numel($(esc(b)))
        for i = 1:n
            ai, bi = getindex($(esc(a)), i), getindex($(esc(b)), i)
            @test typeof(ai) == typeof(bi)
            @test_approx_eq ai bi
        end
    end
end

# a macro to test if tow types are approximatey equal
macro types_approx_eq(a, b)
    quote
        @test typeof($(esc(a))) == typeof($(esc(b)))
        @contents_approx_eq($(esc(a)), $(esc(b)))
    end
end

# a macro to test if the contents of two containers are approximatley equal, without examining the element types
macro contents_approx_eq_notype(a, b)
    quote
        n = numel($(esc(a)))
        @test n == numel($(esc(b)))
        for i = 1:n
            ai, bi = getindex($(esc(a)), i), getindex($(esc(b)), i)
            @test_approx_eq ai bi
        end
    end
end


# build a full list of rotation types including the different ordering schemas
rot_types = Vector{Any}(0)
for rt in Rotations.RotTypeList
    if (Rotations.n_params(rt) == 2)

        # get the super type for the order parameter
        order_type = super(Rotations.default_params(rt)[1])
        for order in subtypes(order_type)
            push!(rot_types, rt{order})
        end
    end

    # ordered ones should have defaults so leave them in
    push!(rot_types, rt)
end

# define null rotations for conveniences
null_rotation(::Type{RotMatrix}) = RotMatrix(eye(3))
null_rotation(::Type{Quaternion}) = Quaternion(1.0,0,0,0)
null_rotation(::Type{AngleAxis}) = AngleAxis(0.0, 1.0, 0, 0)
null_rotation(::Type{SpQuat}) = SpQuat(0.0,0.0,0.0)
null_rotation(::Type{EulerAngles}) = EulerAngles(0.0,0.0,0.0)
null_rotation(::Type{ProperEulerAngles}) = ProperEulerAngles(0.0,0.0,0.0)
null_rotation{ORD}(::Type{EulerAngles{ORD}}) = EulerAngles{ORD, Float64}(0.0,0.0,0.0)
null_rotation{ORD}(::Type{ProperEulerAngles{ORD}}) = ProperEulerAngles{ORD, Float64}(0.0,0.0,0.0)


#########################################################################
# Check fixed relationships
#########################################################################

# Do no rotation
R = RotMatrix(eye(3))
# println("********************************\nIndentity checks\n********************************\n")
for rt in rot_types
    # println(rt)
    rot_var = rt(R)
    null_var = null_rotation(rt)
    @types_approx_eq(rot_var, null_var)
end

#########################################################################
# Check conversion to and from mutable and immutable vectors,
# as well as element conversions
#########################################################################


# println("\n\n\n********************************\nVector conversion checks\n********************************\n")
R = RotMatrix(eye(3))
eltypes = subtypes(AbstractFloat)  # only abstarct floats are supported by all
for rt in rot_types

    # println("$(rt)")
    rot_var = rt(R)
    
    # export to immutable
    ivu = Vec(rot_var)
    @contents_approx_eq(rot_var, ivu)

    # export to mutable
    mvu = Vector(rot_var)
    @contents_approx_eq(rot_var, mvu)

    # import from immutable
    rot_ivu = rt(ivu)
    @types_approx_eq(rot_var, rot_ivu)

    # import from mutable
    rot_ivu = rt(mvu)
    @types_approx_eq(rot_var, rot_ivu)

    # test typed stuff
    for eT in eltypes
        
        # export to immutable
        ivt = Vec{Rotations.nvars(rt), eT}(rot_var)
        ivc = convert(Vec{Rotations.nvars(rt), eT}, ivu)
        @types_approx_eq(ivt, ivc)
    
        # export to mutable
        mvt = Vector{eT}(rot_var)
        mvc = convert(Vector{eT}, mvu)
        @types_approx_eq(mvt, mvc)

        # import from immutable
        rot_ivt = rt(ivt)
        @contents_approx_eq_notype(rot_var, rot_ivt)

        # import from mutable
        rot_mvt = rt(mvt)
        @contents_approx_eq_notype(rot_var, rot_mvt)

        # and test the element conversion on the rotation parameterization directly
        rot_c = convert(Rotations.add_params(rt, eT), rot_var)
        @contents_approx_eq(rot_c, rot_mvt)

    end
end


#########################################################################
# Test conversions between rotation types
#########################################################################

# println("\n\n\n********************************\nRotation conversion checks\n********************************\n")

# test random round trip conversions
repeats = 1000
thresh = 1e-6
eye3 = @fsa([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
for rt_in in rot_types

    for rt_out in rot_types

        # println("$(rt_in) - > $(rt_out)")

        # and each test
        #fcount = 0
        for i = 1:repeats
            
            # start with a random quaternion
            q = nquatrand()
            X = convert_rotation(rt_in, q)

            # round trip conversion
            Xd = convert_rotation(rt_in, convert_rotation(rt_out, X))

            # compare rotations before and after the round trip
            Rout = RotMatrix(X) * RotMatrix(Xd)'  # should be the ident
            rd = vecnorm(eye3 - Rout)
            @test rd <= thresh
            #fcount += (rd > thresh)
        end
        #if (fcount > 0)
        #    warn("Failed $(fcount) / $(repeats): $(rt_in) - > $(rt_out)")
        #end
    end
end



###
# For bench marking the set_subnormal method of dealing with subnormal numbers
# the result of this is that playing with the subnormal mode is slow
function testsubnormal(n::Int=10000000)

    function rottoeuler2{T}(::Type{ProperEulerAngles{Rotations.EulerXYX}}, R::RotMatrix{T})

        # temporarily dis allow sub normal numbers
        sn = get_zero_subnormals();
        set_zero_subnormals(true);
	    t1 = atan2(R[2, 1], -R[3, 1])  
        set_zero_subnormals(sn);
	    
        ct1, st1 = cos(t1), sin(t1)

	    ProperEulerAngles{Rotations.EulerXYX,T}(
		    T(t1),
		    T(atan2((R[1, 2] * R[1, 2] + R[1, 3] * R[1, 3])^(1/2), R[1, 1])),
		    T(atan2(- R[2, 3]*ct1 - R[3, 3]*st1, R[2, 2]*ct1 + R[3, 2]*st1))
		    )
    end

    function calc1(R, n)
        for i = 1:n
            ea = Rotations.rottoeuler(ProperEulerAngles{Rotations.EulerXYX}, R)
        end
    end

    function calc2(R, n)
        for i = 1:n
            ea = rottoeuler2(ProperEulerAngles{Rotations.EulerXYX}, R)
        end
    end


    # a random rotation for testing
    q = nquatrand()
    R = RotMatrix(q)


    # and test
    println("DONT set subnormals to zero")
    @time calc1(R, n)

    println("DO set subnormals to zero")
    @time calc2(R, n)

end













