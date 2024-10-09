module AffineMotions

using Manifolds

import ManifoldDiffEq
import OrdinaryDiffEq

using ManifoldGroupUtils
import ManifoldGroupUtils as GU

export AbstractMotion, AbstractAffineMotion,
    RigidMotion, TranslationMotion,
    FlatAffineMotion,
    ZeroMotion,
    get_flat_action,
    integrate,
    # compose_adjoint,
    swap_group_motion

include("Motion.jl")
include("Utils.jl")

include("../ext/AdjointLinearExt/AdjointLinearExt.jl")

end
