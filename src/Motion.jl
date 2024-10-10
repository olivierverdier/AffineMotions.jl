@doc raw"""
General interface for maps of the form ``φ \colon M \to \mathfrak{g}``,
where ``M`` is a manifold, and ``\mathfrak{g}`` is the Lie algebra
of a Lie group ``G`` action on the manifold ``M``.

```julia
G = SpecialOrthogonal(3)
ξ = rand(TangentSpace(G, identity_element(G)))
M = Sphere(2)
A = RotationAction(M, G)
# use a rigid motion for the sake of the example
φ = RigidMotion(A, ξ)
x = rand(M)
φ(x) # ξ
φ'(x) # linear operator Alg(G) → Alg(G)
integrate(φ, x) # solution of x'=φ(x)⋅x
0.1*φ # rescaled motion
```
"""

abstract type AbstractMotion{TA<:AbstractGroupAction{LeftAction}} end


"""
    get_action(m::AbstractMotion) :: AbstractGroupAction

The action for the given motion.
"""
get_action(m::AbstractMotion) = m.A

#--------------------------------
# AbstractMotion Interface
#--------------------------------
"""
    get_dynamics(::AbstractMotion, x) :: Alg(G)

The actual motion x -> φ(x) in Alg(G).

**Shortcut**: `m(x)` ≡ `get_dynamics(m, x)`
"""
function get_dynamics end

@doc raw"""
    get_lin_at(φ::AbstractMotion, ::Any) :: Function

The linear operator ``ξ ↦ ⟨dφ,ξ⟩_x+[φ(x),ξ]`` at the point ``x`` for motion ``φ``.
"""
function get_lin_at end
#--------------------------------


(m::AbstractMotion)(x) = get_dynamics(m, x)

Base.adjoint(m::AbstractMotion) = x -> get_lin_at(m, x)

abstract type AbstractAffineMotion{TA} <: AbstractMotion{TA} end
abstract type CompositeAffineMotion{TA} <: AbstractAffineMotion{TA} end
abstract type SimpleAffineMotion{TA} <: AbstractAffineMotion{TA} end


@doc raw"""
    rescale_motion(s::Number, φ::Motion) :: Typeof(φ)

Rescale the motion by ``s``.
This corresponds to the new motion ``(sφ)(x) := s(φ(x))``.
The motion `φ` is conveniently rescaled with the syntax `s*φ`.
"""
function rescale_motion end

get_lin_at(m::AbstractAffineMotion, ::Any) = get_lin(m)

Base.:*(s::Number, m::AbstractAffineMotion) = rescale_motion(s, m)

"""
Compute the tangent vector corresponding to the motion
at the point x on the manifold.
"""
Base.:*(m::AbstractMotion, x) = apply_diff_group(get_action(m), Identity(base_group(get_action(m))), m(x), x)

# include("Motion/AffineMotion.jl")


Base.:-(m::AbstractAffineMotion) = -1*m



# Integration

function _prepare_problem(
    G, # Group
    dyn,
    method,
    )
    u0 = identity_element(G)
    diffop = ManifoldDiffEq.LieManifoldDiffEqOperator{Float64}(dyn)
    prob = ManifoldDiffEq.ManifoldODEProblem(diffop, u0, (0, 1.0), G)
    # new:
    A = GroupOperationAction(G)
    algo = method(G, Manifolds.GroupExponentialRetraction(), A)
    return prob, algo
end

# function integrate_problem(prob, algo, dt, method)
#     # alg = ManifoldDiffEq.ManifoldLieEuler(M, ExponentialRetraction(), action)
#     sol = OrdinaryDiffEq.solve(prob, algo, dt=dt)
#     return sol
# end


"""
    integrate_lift(m::AbstractMotion, x0) :: AbstractODESolution

Integrate the lifted motion in the acting group.
"""
function integrate_lift(
    m::AbstractMotion,
    x0;
    dt=0.1,
    method=ManifoldDiffEq.RKMK4
    )
    action = get_action(m)
    dyn(χ,p,t) = m(apply(action, χ, x0))
    # prob, algo = _prepare_problem(action, dyn, method)
    prob, algo = _prepare_problem(base_group(action), dyn, method)
    # return integrate_problem(prob, algo, dt, method)
    return OrdinaryDiffEq.solve(prob, algo, dt=dt)
end

"""
    integrate(motion, x0::TM) :: TM

Integrate the motion `motion` from the point `x0` on the manifold.
"""
function integrate(
    motion::AbstractMotion,
    x0;
    dt=0.1
)
    action = get_action(motion)
    # manifold = group_manifold(action)
    # @assert is_point(manifold, x0)

    sol = integrate_lift(motion, x0; dt)
    χ = last(sol.u)
    x = apply(action, χ, x0)
    return x
end


"""
 Compute matrix for Dφ in the given basis
"""
function get_lin_mat(m::AbstractMotion, x, B::AbstractBasis)
    G = base_group(get_action(m))
    op = get_lin_at(m, x)
    return GU.matrix_from_lin_endomorphism(G, op, B)
end




@doc raw"""
     morphism(φ::Motion, x, B::AbstractBasis)

Return the Lie algebra morphism ``\exp(Dφ)``,
where ``Dφ`` is the linear part of the motion ``φ``.
"""
morphism(motion, x, B) = begin
    mm = get_lin_mat(motion, x, B)
    return exp(mm)
end

include("Motion/Simple/Rigid.jl")
include("Motion/Simple/Translation.jl")
include("Motion/Simple/FlatAffine.jl")
include("Motion/Sum.jl")
include("Motion/AffineMotion.jl")
include("Motion/GroupMotion.jl")




