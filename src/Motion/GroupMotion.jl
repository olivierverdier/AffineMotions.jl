
"""
    swap_group_motion(::Motion)

Swap from standard to dual group action, while preserving the same dynamics.
"""
swap_group_motion(m::AbstractRigidMotion{<:GroupOperationAction{LeftAction, LeftSide}}) = TranslationMotion(base_group(get_action(m)), get_velocity(m), RightSide())
swap_group_motion(m::TranslationMotion{RightSide}) = RigidMotion(GroupOperationAction(m.G), m.vel)
swap_group_motion(m::AbstractRigidMotion{<:GroupOperationAction{LeftAction, RightSide}}) = TranslationMotion(base_group(get_action(m)), get_velocity(m), LeftSide())
swap_group_motion(m::TranslationMotion{LeftSide}) = RigidMotion(GroupOperationAction(m.G, (LeftAction(), RightSide())), m.vel)

# _swap_inv(::GroupOperationAction, G, χ) = inv(G,χ)
# _swap_inv(::GroupOperationAction{LeftAction, RightSide}, G, χ) = χ
_swap_adjoint_action(::GroupOperationAction{LeftAction, LeftSide}, G, χ, ξ) = adjoint_action(G, χ, ξ, RightAction())
_swap_adjoint_action(::GroupOperationAction{LeftAction, RightSide}, G, χ, ξ) = adjoint_action(G, χ, ξ)
_swap_group_action(A::GroupOperationAction{LeftAction, LeftSide}) = GroupOperationAction(base_group(A), (LeftAction(), RightSide()))
_swap_group_action(A::GroupOperationAction{LeftAction, RightSide}) = GroupOperationAction(base_group(A))

