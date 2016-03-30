"""

Function to convert a covariance matrix containing one rotation parameterization to another:

Inputs:

    tfunc     - transformation function, should take a Vector with length size(cX, 1) as an input and return the transformed state vector

    cX        - covariance matrix

    xbar      - mean X in the input space (set = nothing to enforce zero mean in input and output)

    oT_cov    - output type for the covariance matrix (default Matrix)

    oT_mean   - output type for the covariance matrix (default Vector)


Experimental - This is not a fast implementation

"""
function convert_covariance{T, U, V}(tfunc, cX, xbar::V=zeros(size(cX,1)), ::Type{T}=Any, ::Type{U}=Any)

    # generate the sigma points
    (sigmas, weights) = sigma_points(cX, xbar)

    # get the output size
    siz = ((T == Any) ? length(tfunc(sigmas[:,1])) : size(T,1))::Int

    # and convert
    sigma_trans = zeros(siz, size(sigmas, 2))
    for i = 1:size(sigmas,2)
        sigma_trans[:,i] = tfunc(sigmas[:,i])     # type it as the input type before converting to the output type
    end

    # and build the output mean
    xbar_out = (V == Void) ? nothing : vec(mean(sigma_trans,2))

    # convert the sigma points back into a covariance
    cX = sigma_point_cov(sigma_trans, weights, xbar_out)
    if (T != Any)
        cX = T(cX)::T
    end
    if (V != Void) && (U != Any)
        xbar_out = U(xbar_out)::U
    end

    return cX, xbar_out  # return the mean as well
end



# generate sigma points for a specifified covaraince
function sigma_points{V}(cX, xbar::V=nothing)

    if !isa(cX, AbstractMatrix)
        cX = Matrix(cX)  # chol works with Matrix not other types
    end

    # save some time and use and existing ukf approach instead of actual maths
    n = size(cX,1)

    W0 = -0.2
    R = full(chol((n/(1-W0)) * cX, Val{:L}))::Matrix{Float64}

    # weights
    weights = zeros(2*n+1)
    weights[1] = W0
    weights[2:end] = (1.0-W0)/(2.0*n)

    # sigma points
    if (V != Void)
        xbar = Vector(xbar)
        sigma_points = [xbar  (xbar .+  R)   xbar .- R]
    else
        sigma_points = [zeros(size(cX,1))  R   -R]
    end

    #TODO: Need to put in special catch if angle sigma points get too big so they don't wrap around
    return sigma_points, weights

end


# generate covariance from sigma points
function sigma_point_cov(sigma_points, weights, xbar=nothing)

    # recover the covariance
    #cX_out = (sigmaPointWeightsObs'.* sigmaPoints) * sigmaPoints'
    if xbar != nothing
        sigma_points = sigma_points .- xbar
    end
    cov = sigma_points * (sigma_points' .* weights)  # maybe better than the above (less copy but mult across rows)

    #This is poop!! Force covariance matrix symmetric.
    cov = 0.5*(cov + cov')
    return cov

end



#=
# investigate the distribution of SpQuat vars relative to Gaussian Euler Angles
# TLDR it appears pretty Gaussian itself
using PyPlot
function test_dist(std=10.0, n::Int=100_000)

    # generate the Euler angles
    cX_euler = diagm(fill(std,3)*pi/180) .^ 2  # keep it small angle
    sX_euler = chol(cX_euler, Val{:L})
    ea_vec = [EulerAngles(sX_euler * randn(size(sX_euler, 1))) for i in 1:n]

    # convert them
    x, y, z = zeros(n), zeros(n), zeros(n)
    for (i, ea) in enumerate(ea_vec)
        spq = SpQuat(ea)
        x[i], y[i], z[i] = (spq...)
    end

    # show the results
    figure();
    plt[:hist](x, div(n, 50)); title("x (std=$(std))")
    figure();
    plt[:hist](y, div(n, 50)); title("y (std=$(std))")
    figure();
    plt[:hist](z, div(n, 50)); title("z (std=$(std))")
    return (x,y,z)

end
=#

