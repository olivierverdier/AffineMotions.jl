using Test
using AffineMotions

import ManifoldGroupUtils: rand_lie

using Manifolds

import MultiAffine: MultiDisplacement, MultiAffineAction
import LinearAlgebra

using Random
rng = Random.default_rng()


# Base.length(x::ProductRepr) = sum(map(Base.length, submanifold_components(x)))
# Base.convert(::Type{ArrayPartition}, x::ProductRepr) = ArrayPartition(submanifold_components(x)...)
# Base.convert(::Type{ProductRepr}, x::ArrayPartition) = ProductRepr(submanifold_components(x)...)


@testset "Zero Motion" begin
    G = SpecialOrthogonal(3)
    M = Sphere(2)
    A = RotationAction(M, G)
    m = ZeroMotion(A)
    x = [1., 0, 0]
    @test m(x) ≈ zero_vector(G, identity_element(G))
    ξ = rand_lie(rng, G)
    @test isapprox(G, m'(x)(ξ), zero_vector(G, identity_element(G)))
    @test RigidMotion(A) isa ZeroMotion
end

@testset "Motion Composition" begin
    G = MultiDisplacement(4,2)
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
    @test_throws MethodError RigidMotion(action, vel) + TranslationMotion(G, vel, RightSide())
    @testset "Sum/Rescale $m" for m in motions
        @test m ≈ m
        @test .5*m isa typeof(m)
        @test 2*(.5*m) ≈ m
        @test 2*m ≈ m+m broken=isa(m, AffineMotions.AffineMotionSum)
        # if m is AffineMotionSum{TA, TV}, the sum is AffineMotionSum{TA, TV'} with another TV, hence the following two cases:
        if m isa AffineMotions.AffineMotionSum
            @test m+m isa AffineMotions.AffineMotionSum
        else
            @test m+m isa typeof(m)
        end
    end
end

@testset "Motion Sum" begin
    G = MultiDisplacement(4,2)
    m1 = AdjointLinearMotion(G, ones(2,2), LeftSide())
    ξ = rand_lie(rng, (G))

    motions = (
    rm = RigidMotion(GroupOperationAction(G), ξ),
    tm = TranslationMotion(G, ξ, LeftSide()),
    lm = AdjointLinearMotion(G, randn(rng, 2,2), LeftSide()),
    )

    @testset "Sum type" for m in motions
        m+m isa typeof(m)
        m+m+m isa typeof(m)
    end

    rm,tm,lm = motions

    @test (rm + tm) + lm ≈ rm + (tm + lm)

    @test (rm + tm) + (rm + tm) isa AffineMotions.AffineMotionSum
end



include("test_integ.jl")

include("test_swap.jl")

include("test_flat.jl")