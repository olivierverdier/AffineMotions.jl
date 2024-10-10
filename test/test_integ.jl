using Manifolds

import LinearAlgebra
import Random

rng = Random.default_rng()

@doc raw"""
If ``F`` denotes the flow of the motion ``φ``,
and if ``E = \exp(Dφ)``
then
```math
F(\exp(ξ) ⋅ x_0) = \exp(Eξ) ⋅ F(x_0)
```
"""
check_geodesic_preservation(m::AbstractAffineMotion, x0, ξ, B=DefaultOrthogonalBasis()) = begin
    A = AffineMotions.get_action(m)
    G = base_group(A)

    # compute F(exp(ξ) ⋅ x0)
    χ1 = exp_lie(base_group(A), ξ)
    x0_ = apply(A, χ1, x0)
    left = integrate(m, x0_)

    # compute E = exp(Dφ) (a matrix)
    E = AffineMotions.morphism(m, x0, B)

    # compute exp(Eξ) ⋅ F(x0)
    X = get_coordinates_lie(G, ξ, B)
    Eξ = get_vector_lie(G, E * X, B)
    χ = exp_lie(G, Eξ)
    x1 = integrate(m, x0)
    right = apply(A, χ, x1)

    return isapprox(G, left, right)
end

@testset "geodesic" for G in [
    MultiDisplacementGroup(3, 2),
]
    ξ = rand_lie(rng, G)
    ξ0 = rand_lie(rng, G)

    A = GroupOperationAction(G, (LeftAction(), LeftSide()))

    @testset "m" for m in [
        RigidMotion(A, ξ0),
        TranslationMotion(G, ξ0, LeftSide()),
        TranslationMotion(G, ξ0, RightSide()),
    ]
        A = AffineMotions.get_action(m)
        x0 = rand(rng, group_manifold(A))
        @test check_geodesic_preservation(m, x0, ξ)
    end
end

@testset "∫ $name" for (name, G) in [
    "SO4" => SpecialOrthogonal(4),
    "MD4" => MultiDisplacementGroup(4),
    "MD43" => MultiDisplacementGroup(4, 2),
    "SE4" => SpecialEuclidean(4),
]
    action = GroupOperationAction(G)
    vel = rand_lie(rng, G)
    # submanifold_component(vel, 1) .= 0
    # rm = make_rigid_motion(action, vel)
    rm = RigidMotion(action, vel)


    sol = AffineMotions.integrate_lift(1.0*rm, identity_element(G); dt=.01)
    # test that sol(1) ≡ exp(ξ), for a rigid motion ξ
    expected = exp_lie(G, vel)

    @testset "Integration" begin
        computed = last(sol.u)
        # @error G vel computed expected
        @test isapprox(G, computed, expected)

        rm_ = RigidMotion(action, -vel)
        sol_ = AffineMotions.integrate_lift(rm+rm_, identity_element(G);dt=.01)
        id = last(sol_.u)
        @test isapprox(G, id, identity_element(G))

        tm = TranslationMotion(G, vel, LeftSide())
        sol = AffineMotions.integrate_lift(tm, identity_element(G); dt=.01)
        @test isapprox(G, last(sol.u), exp_lie(G, -vel))
    end

    # vel_ = rand_lie(rng, G)
    # v3 = lie_bracket(G, vel, vel_)
    # @show vel_
    # @show adjoint_action(G, inv(G, χ), v3)

    @testset "compose adjoint" begin
        B = DefaultOrthogonalBasis()
        # mm = get_adjoint_matrix(G, vel, B)
        mm = AffineMotions.get_lin_mat(rm, identity_element(G), B)
        res_mat = compose_adjoint(G, inv(G, expected), exp(mm), B)
        # display(res_mat)
        @test isapprox(res_mat, LinearAlgebra.I)
    end



end
