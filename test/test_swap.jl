using Manifolds

function test_constant_dict(d::Dict, comp)
    if length(d) == 0
        return
    end
    v_ = first(d).second
    @testset "Testing $k" for (k, v) in d
        @test comp(v_, v)
    end
end

@testset "Swap Group Motion $G" for G in [
    SpecialOrthogonal(3),
    # SpecialEuclidean(3),
    MultiDisplacementGroup(3, 1),
    MultiDisplacementGroup(3, 2),
    ]
    ξ = rand_lie(rng, G)

    # x0 = identity_element(G)
    x0 = rand(rng, G)

    @testset "Rigid/Translation" for side in
        [
            LeftSide(),
            RightSide(),
         ]

        R1 = RigidMotion(GroupOperationAction(G, (LeftAction(), side)), ξ)
        T1 = TranslationMotion(G, ξ, switch_side(side))
        M1 = Dict(
            :R1 => R1,
            :R1_ => swap_group_motion(R1),
            :R1__ => AffineMotions._swap_group_motion(R1),
            :T1 => T1,
            :T1_ => swap_group_motion(T1),
            :T1__ => AffineMotions._swap_group_motion(T1),
        )

        rM1 = Dict([(s => integrate(x0, v)) for (s,v) in M1]...)
        # for v in values(rM1)
        #     @test isapprox(G, rM1[:R1], v)
        # end
        test_constant_dict(rM1, (a,b)->isapprox(G, a, b))

    end

end


@testset "Swap AdjointLinearMotion" begin
    G = MultiDisplacementGroup(3,2)
    x0 = rand(rng, G)
    m1 = AdjointLinearMotion(G, [1.0 0;0 0], LeftSide())
    m2 = swap_group_motion(m1)
    MAM = Dict(
        :m1 => integrate(x0, m1),
        :m2 => integrate(x0, m2),
    )
    test_constant_dict(MAM, (a,b) -> isapprox(G, a, b))
end

@testset "Swap AdjointLinear sum" begin
    G = MultiDisplacementGroup(3,2)
    ξ = rand_lie(rng, G)
    m1 = AdjointLinearMotion(G, [1.0 0;0 0], LeftSide())
    rm = RigidMotion(GroupOperationAction(G), ξ)
    @test swap_group_motion(m1+rm) isa AffineMotions.AffineMotionSum
end
