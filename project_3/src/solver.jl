
#=
As opposed to Python, methods in Julia work as function(object, args...)
instead of object.function(args...).
As a result, all the functions relevant for this solver,
including all the difference schemes, time step finder, and hydro-solver (step function)
are simply functions specialized to get a solver of type Star as input.
=#

### this is like a class in Python, except it just holds the variables, and not functions
struct Star
    g::Float64
    dx::Float64
    dy::Float64
    M0::Float64
    R0::Float64
    T0::Float64
    P0::Float64
    ρ::Array{Float64}
    u::Array{Float64}
    w::Array{Float64}
    e::Array{Float64}
    T::Array{Float64}
    P::Array{Float64}
end


function perturbation!(T::Array{Float64}; σ=15, A=10000)
    """
    Function for adding gaussing noise to an array T,
    with standard deviation σ and amplitude A.
    This is an in-place function (hence the !).
    """

    ny, nx = size(T)

    ### position of maximum amplitude ###
    y0 = Int(ny/2 - 0.25*ny) #
    x0 = Int(nx/2)

    ### mask ###
    x = reshape(collect(LinRange(0, 1, nx)), 1, nx)
    y = reshape(collect(LinRange(0, 1, ny)), ny, 1)

    μy = y[y0]
    μx = x[x0]

    T .+= A*exp.(-0.5*σ^2*((x .- μx).^2 .+ (y .- μy).^2))
end


function initialize(M0::Float64=1.989e30, R0::Float64=6.96e8,
                    T0::Float64=5778.0,   P0::Float64=1.8e8;
                    ∇::Float64=0.402,
                    nx::Int64=300, ny::Int64=100,
                    disturbe::Bool=false)

    global mᵤ     = 1.660539066e-27 # atomic mass [kg]
    global k_B    = 1.3806e-23      # Boltzmann constant [J/K]
    global μ      = 0.61            # mean molecular weight
    global γ      = 5/3             # ratio between specific heat capacities (ideal gas)
    global G      = 6.674e-11       # gravitational constant [m³·kg⁻¹·s⁻²]
    global height = 4.0e6           # height of the computational box (y-range) [m]
    global width  = 1.2e7           # width of the computational box (x-range) [m]

    dx = width/nx   # step size in x-direction
    dy = height/ny  # step in in y-direction
    g = -G*M0/R0^2  # gravitational acceleration

    T = zeros(Float64, ny, nx) # temperature array
    P = zeros(Float64, ny, nx) # pressure array

    y_range = LinRange((ny - 1)*dy, 0, ny) # discrete y-values in the box (like np.linspace)

    ### initialize T and P ###
    A = ∇*g*μ*mᵤ/k_B

    @inbounds for (j, y) in enumerate(y_range)
        T[j,:] .= T0 - A*y
    end

    P[:,:] = (T/T0).^(1/∇)*P0

    ### add gaussian perturbation if requested ###
    disturbe ? perturbation!(T) : nothing

    ### initialize energy and density ###
    e = P/(γ - 1)
    ρ = (μ*mᵤ/k_B)*P ./ T

    ### initialize horizontal and vertical velocity ###
    u = zeros(Float64, ny, nx)
    w = zeros(Float64, ny, nx)

    ### initialise and return star object ###
    Star(g, dx, dy, M0, R0, T0, P0, ρ, u, w, e, T, P)
end

function equilibrium(model::Star)
    """ Check whether star is in hydrostatic equilibrium """

    println("\nChecking hydrostatic equilibrium criterium\n")

    dPdy = central_y(model.P, model.dx)
    gρ = model.g*model.ρ

    diff = abs.((gρ - dPdy))[2:end-1,:]
    meandiff = sum(diff)/length(diff)
    maxdiff = maximum(diff)

    println("Mean absolute difference: $(round(meandiff, digits=2))")
    println("Maximum absolute difference: $(round(maxdiff, digits=2))\n")

    return diff
end

function central_x(func::Array{Float64}, dx::Float64)
    """
    Central difference scheme in x-direction
    with periodic boundary conditions
    """

    nx = size(func)[2]
    deriv = zeros(Float64, size(func))

    ### inner points ###
    @inbounds for i = 2:nx-1
        deriv[:,i] .= func[:,i+1] - func[:,i-1]
    end

    ### boundaries ###
    deriv[:,1] = func[:,2] - func[:,end]
    deriv[:,end] = func[:,1] - func[:,end-1]

    return deriv*0.5/dx

end

function upwind_x(func::Array{Float64}, vel::Array{Float64}, dx::Float64)

    """
    upwind difference scheme in x-direction
    """

    ny, nx = size(func)
    deriv = zeros(size(func))

    ### inner points ###
    @inbounds for j = 1:ny
        @inbounds for i = 2:nx-1
            if vel[j, i] >= 0
                deriv[j, i] = func[j,i] - func[j,i-1]
            else
                deriv[j, i] = func[j,i+1] - func[j,i]
            end
        end
    end

    ### boundaries ###
    @inbounds for j = 1:ny
        if vel[j, 1] >= 0
            deriv[j, 1] = func[j, 1] - func[j, end]
        else
            deriv[j, 1] = func[j,2] - func[j,1]
        end

        if vel[j, end] >= 0
            deriv[j, end] = func[j, end] - func[j, end-1]
        else
            deriv[j, end] = func[j,1] - func[j,end]
        end

    end

    return deriv/dx
end

function central_y(func, dy)
    """
    Central difference scheme in y-direction (inner points)
    """

    ny = size(func)[1]
    deriv = zeros(Float64, size(func))

     @inbounds for j = 2:ny-1
        deriv[j,:] = func[j+1,:] - func[j-1,:]
    end

    return deriv*0.5/dy
end

function upwind_y(func::Array{Float64}, vel::Array{Float64}, dy::Float64)

    """
    upwind difference scheme in x-direction
    """

    ny, nx = size(func)
    deriv = zeros(size(func))

    @inbounds for j = 2:ny-1
        @inbounds for i = 1:nx
            if vel[j, i] >= 0
                deriv[j, i] = func[j,i] - func[j-1,i]
            else
                deriv[j, i] = func[j+1,i] - func[j,i]
            end
        end
    end

    return deriv/dy

end

function boundary_conditions(model::Star)
    """ Sets the vertical boundary conditions for the model """

    ### vertical velocity ###
    model.w[1,:] .= 0
    model.w[end,:] .= 0

    ### horizontal velocity ###
    model.u[1,:] = (4*model.u[2,:] .- model.u[3,:])/3.0
    model.u[end,:] = (4*model.u[end-1,:] .- model.u[end-2,:])/3.0

    ### energy ###
    C1 = -model.g*μ*mᵤ/k_B
    model.e[1,:] = (4*model.e[2,:] .- model.e[3,:]) ./ (3 .- 2*model.dy*C1 ./ model.T[1,:])
    model.e[end,:] = (4*model.e[end-1,:] .- model.e[end-2,:]) ./ (3 .+ 2*model.dy*C1 ./ model.T[end,:])

    ### density ###
    C2 = (γ - 1)*μ*mᵤ/k_B
    model.ρ[1,:] = C2*model.e[1,:] ./ model.T[1,:]
    model.ρ[end,:] = C2*model.e[end,:] ./ model.T[end,:]

end

function continuity_equation(model::Star)
    du_dx = central_x(model.u, model.dx)
    dw_dy = central_y(model.w, model.dy)

    dρ_dx = upwind_x(model.ρ, model.u, model.dx)
    dρ_dy = upwind_y(model.ρ, model.w, model.dy)

    dρ_dt =  - model.ρ .* (du_dx .+ dw_dy) .-
               model.u .* dρ_dx .-
               model.w .* dρ_dy

    return dρ_dt

end

function momentum_equation_horizontal(model::Star)
    du_dx = upwind_x(model.u, model.u, model.dx)
    dw_dy = upwind_y(model.u, model.u, model.dy)

    ρu = model.ρ .* model.u

    dρu_dx = upwind_x(ρu, model.u, model.dx)
    dρu_dy = upwind_y(ρu, model.w, model.dy)

    dP_dx = central_x(model.P, model.dx)

    dρu_dt = - ρu .* (du_dx .+ dw_dy) .-
               model.u .* dρu_dx .-
               model.w .* dρu_dy .-
               dP_dx

    return dρu_dt, ρu

end

function momentum_equation_vertical(model::Star)
    du_dx = upwind_x(model.u, model.w, model.dx)
    dw_dy = upwind_y(model.u, model.w, model.dy)

    ρw = model.ρ .* model.w

    dρw_dx = upwind_x(ρw, model.u, model.dx)
    dρw_dy = upwind_y(ρw, model.w, model.dy)

    dP_dy = central_y(model.P, model.dy)

    dρw_dt =  - ρw .* (du_dx .+ dw_dy) .-
                model.u .* dρw_dx .-
                model.w .* dρw_dy .-
                dP_dy .+ model.g*model.ρ

    return dρw_dt, ρw

end

function energy_equation(model::Star)
    du_dx = central_x(model.u, model.dx)
    dw_dy = central_y(model.w, model.dy)

    de_dx = upwind_x(model.e, model.u, model.dx)
    de_dy = upwind_y(model.e, model.w, model.dy)

    de_dt = - model.e .* (du_dx .+ dw_dy) .-
              model.u .* de_dx .-
              model.w .* de_dy .-
              model.P .* (du_dx .+ dw_dy)

    return de_dt

end

function timestep(model::Star, dρ_dt, dρu_dt, dρw_dt, de_dt, ρw, ρu)
    """ Find variable time step length """

    ### velocities close to zero ###
    u_zero = abs.(model.u) .< 1e-6
    w_zero = abs.(model.w) .< 1e-6

    ### if all velocities are zero ###
    if all(u_zero) && all(w_zero)
        return 1e-2

    ### if just all u-velocities are zero ###
    elseif all(u_zero) && !all(w_zero)
        w_ignore = w_zero .== false

        max_ρ = maximum(abs.(dρ_dt ./ model.ρ))
        max_w = maximum(abs.(dρw_dt[w_ignore] ./ ρw[w_ignore]))
        max_e = maximum(abs.(de_dt ./ model.e))
        max_x = maximum(abs.(model.w[w_ignore] ./ model.dy))

        δ = max(max_ρ, max_w, max_e, max_x)

    ### if just all w-velocities are zero ###
    elseif all(w_zero) && !all(u_zero)
        u_ignore = u_zero .== false

        max_ρ = maximum(abs.(dρ_dt ./ model.ρ))
        max_u = maximum(abs.(dρu_dt[u_ignore] ./ ρu[u_ignore]))
        max_e = maximum(abs.(de_dt ./ model.e))
        max_y = maximum(abs.(model.u[u_ignore] ./ model.dx))

        δ = max(max_ρ, max_u, max_e, max_x)

    ### if both velocities have non-zero values ###
    else
        w_ignore = w_zero .== false
        u_ignore = u_zero .== false

        max_ρ = maximum(abs.(dρ_dt ./ model.ρ))
        max_w = maximum(abs.(dρw_dt[w_ignore] ./ ρw[w_ignore]))
        max_u = maximum(abs.(dρu_dt[u_ignore] ./ ρu[u_ignore]))
        max_e = maximum(abs.(de_dt ./ model.e))
        max_x = maximum(abs.(model.w[w_ignore] ./ model.dy))
        max_y = maximum(abs.(model.u[u_ignore] ./ model.dx))

        δ = max(max_ρ, max_w, max_u, max_e, max_x, max_y)

    end

    Δt = 0.1/δ

    # return 1e-2 if time step is smaller than 1e-2, to speed up
    Δt < 1e-2 ? Δt = 1e-2 : nothing

    return Δt

end

function hydro_solver(model::Star)

    ### continuity equation ###
    dρ_dt = continuity_equation(model)

    ### momentum equation ###
    dρu_dt, ρu = momentum_equation_horizontal(model)
    dρw_dt, ρw = momentum_equation_vertical(model)

    ### energy equation ###
    de_dt = energy_equation(model)

    ### compute time step ###
    Δt = timestep(model, dρ_dt, dρu_dt, dρw_dt, de_dt, ρw, ρu)

    ### step inner points ###
    model.ρ[2:end-1,:] = (model.ρ .+ dρ_dt*Δt)[2:end-1,:]
    model.u[2:end-1,:] = ((ρu .+ dρu_dt*Δt) ./ model.ρ)[2:end-1,:]
    model.w[2:end-1,:] = ((ρw .+ dρw_dt*Δt) ./ model.ρ)[2:end-1,:]
    model.e[2:end-1,:] = (model.e .+ de_dt*Δt)[2:end-1,:]

    ### set boundary conditions ###
    boundary_conditions(model)

    ### update P and T ###
    model.P[:,:] = (γ - 1)*model.e
    model.T[:,:] = model.e ./ model.ρ * (γ - 1)*(μ*mᵤ)/k_B

    return Δt

end
