rng = Random.default_rng()

@testset for G in [
    SpecialOrthogonal(3),
    ], side in [
        LeftSide(),
        RightSide(),
    ]
    ξ = rand_lie(rng, G)
    χ = rand(rng, G)
    χ_ = apply(GroupOperationAction(G, (LeftAction(), side)), exp_lie(G, ξ), χ)
    v = last(first(AffineMotions.pose_velocities(G, [χ, χ_])))
    m = AffineMotions.rigid_motion(G, χ, v, side)
    computed = integrate(χ, m)
    @test isapprox(G, computed, χ_)
end
