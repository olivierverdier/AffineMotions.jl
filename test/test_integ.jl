using Manifolds

import LinearAlgebra
import Random

rng = Random.default_rng()


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


    sol = AffineMotions.integrate_lift(1.0*rm, identity_element(G), .01)
    # test that sol(1) ≡ exp(ξ), for a rigid motion ξ
    expected = exp_lie(G, vel)

    @testset "Integration" begin
        computed = last(sol.u)
        # @error G vel computed expected
        @test isapprox(G, computed, expected)

        rm_ = RigidMotion(action, -vel)
        sol_ = AffineMotions.integrate_lift(rm+rm_, identity_element(G), .01)
        id = last(sol_.u)
        @test isapprox(G, id, identity_element(G))

        tm = TranslationMotion(G, vel, LeftSide())
        sol = AffineMotions.integrate_lift(tm, identity_element(G), .01)
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
        res_mat = AffineMotions.compose_adjoint(G, inv(G, expected), exp(mm), B)
        # display(res_mat)
        @test isapprox(res_mat, LinearAlgebra.I)
    end



end
