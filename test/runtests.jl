include("header.jl")

using AffineMotions

using ManifoldGroupUtils

using Manifolds

import MultiAffine: MultiDisplacementGroup, MultiAffineAction
import LinearAlgebra

using Random
rng = Random.default_rng()

"""
Compose the adjoint action with the matrix M, i.e.,
given M ∈ Alg(G) ⊗ Alg(G)*,
the new operator is defined by
ξ ↦ χ (Mξ) χ⁻¹
and the function returns the corresponding matrix.
"""
function compose_adjoint(
    G,
    elt,
    morph_mat,
    B::AbstractBasis
    ) 
    op(ξ) = adjoint_action(G, elt, ξ)
    return ManifoldGroupUtils.compose_lie_matrix_op(G, op, morph_mat, B)
end

@doc raw"""
     compute_morphism(φ::Motion, x, B::AbstractBasis)

Integrate the lift of the motion ``φ`` at the point ``x``,
this gives a group element ``χ``.
This allows to compute the associate morphism, i.e., the operator
```math
ξ ↦ χ^{-1} (\exp(Dφ)ξ) χ
```
"""
function compute_morphism(motion, x, B; dt=0.1)
    action = AffineMotions.get_action(motion)

    sol = AffineMotions.integrate_lift(x, motion; dt)
    χ = last(sol.u)


    G = base_group(action)
    mm = AffineMotions.get_lin_mat(motion, x, B)
    morph = compose_adjoint(G, inv(G, χ), exp(mm), B)

    return χ, morph
end


# Base.length(x::ProductRepr) = sum(map(Base.length, submanifold_components(x)))
# Base.convert(::Type{ArrayPartition}, x::ProductRepr) = ArrayPartition(submanifold_components(x)...)
# Base.convert(::Type{ProductRepr}, x::ArrayPartition) = ProductRepr(submanifold_components(x)...)


check_zero_motion(rm::RigidMotion) = begin
    A = AffineMotions.get_action(rm)
    rm0 = rm + ZeroMotion(A)
    return (rm0 isa RigidMotion) && (rm0 ≈ rm)
end

@testset "Zero Motion" begin
    G = SpecialOrthogonal(3)
    M = Sphere(2)
    A = RotationAction(M, G)
    m = ZeroMotion(A)
    x = [1.0, 0, 0]
    @test m(x) ≈ zero_vector(G, identity_element(G))
    ξ = rand_lie(rng, G)
    @test check_zero_motion(RigidMotion(A, ξ))
    @test isapprox(G, m'(x)(ξ), zero_vector(G, identity_element(G)))
    @test RigidMotion(A) isa ZeroMotion
end

@testset "RigidMotion Exceptions" begin
    G = SpecialOrthogonal(3)
    M = Sphere(2)
    A = RotationAction(M, G)
    Test.@test_throws ErrorException RigidMotion(A, 0)
    Test.@test_throws TypeError RigidMotion(GroupOperationAction(G,(RightAction(), RightSide())), rand_lie(rng, G))
end

@testset "Motion Composition" begin
    G = MultiDisplacementGroup(4,2)
    M = randn(rng, 2, 2)
    action = GroupOperationAction(G)
    vel = rand_lie(rng, G)
    motions = [RigidMotion(action, vel),
               TranslationMotion(G,vel,RightSide()),
               AdjointLinearMotion(G, rand(rng, 2,2), LeftSide()),
               FlatAffineMotion(M, zeros(2)),
               ZeroMotion(action),
               RigidMotion(action, vel) + TranslationMotion(G, vel, LeftSide()),
               ]
    Test.@test_throws MethodError RigidMotion(action, vel) + TranslationMotion(G, vel, RightSide())
    @testset "Sum/Rescale $m" for m in motions
        @test m ≈ m
        @test 0.5 * m isa typeof(m)
        @test 2 * (0.5 * m) ≈ m
        Test.@test m + (-m) ≈ 0*m broken = isa(m, AffineMotions.AffineMotionSum)
        Test.@test 2 * m ≈ m + m broken = isa(m, AffineMotions.AffineMotionSum)
        # if m is AffineMotionSum{TA, TV}, the sum is AffineMotionSum{TA, TV'} with another TV, hence the following two cases:
        if m isa AffineMotions.AffineMotionSum
            @test m + m isa AffineMotions.AffineMotionSum
        else
            @test m+m isa typeof(m)
        end
    end
end

@testset "Motion Sum" begin
    G = MultiDisplacementGroup(4,2)
    m1 = AdjointLinearMotion(G, ones(2,2), LeftSide())
    ξ = rand_lie(rng, G)

    motions = (
        rm=RigidMotion(GroupOperationAction(G), ξ),
        tm=TranslationMotion(G, ξ, LeftSide()),
        lm=AdjointLinearMotion(G, randn(rng, 2, 2), LeftSide()),
    )

    @testset "Sum type" for m in motions
        m + m isa typeof(m)
        m + m + m isa typeof(m)
    end

    rm, tm, lm = motions

    @test (rm + tm) + lm ≈ rm + (tm + lm)

    @test (rm + tm) + (rm + tm) isa AffineMotions.AffineMotionSum

    S = rm + tm + lm
    χ = rand(rng, G)
    @test isapprox(algebra(G), S(χ), rm(χ) + tm(χ) + lm(χ))
    @test isapprox(algebra(G), S'(χ)(ξ), rm'(χ)(ξ) + tm'(χ)(ξ) + lm'(χ)(ξ))
    actions = [AffineMotions.get_action(m) for m in motions]
    Sa = AffineMotions.get_action(S)
    @test all(actions) do a
        return Sa == a
    end
end

@testset "Rigid Motion Morphism" begin
    G = SpecialOrthogonal(3)
    M = Sphere(2)
    x = [1,0,0]
    A = RotationAction(M, G)
    ξ = rand_lie(rng, G)
    rm = RigidMotion(A, ξ)
    M = last(compute_morphism(rm, x, DefaultOrthogonalBasis()))
    @test isapprox(M, LinearAlgebra.I)
end

test_files = [
    "test_integ.jl",
    "test_swap.jl",
    "test_flat.jl",
    "test_utils.jl",
]

include_tests(test_files)

