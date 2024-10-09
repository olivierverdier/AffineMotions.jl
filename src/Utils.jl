
@doc raw"""
    pose_velocities(G, poses)

Take a list ``χ_i`` and returns a list ``(χ_i, v_i)``
 such that ``χ_{i+1} = \exp_{χ_{i}}(v_i)``
"""
pose_velocities(G, poses) = [(poses[i], log(G, poses[i], poses[i+1])) for i in 1:length(poses)-1]

_switch_sign(::LeftSide, ξ) = ξ
_switch_sign(::RightSide, ξ) = -ξ

"""
    rigid_motion(G, χ, v, side)

The rigid motion starting at ``χ`` in the direction
``v`` (tangent to ``χ``), using either the standard
left action or the dual one.
"""
rigid_motion(G, χ, v, side) = begin
    ξ = translate_to_id(G, χ, v, switch_side(side))
    return RigidMotion(GroupOperationAction(G, (LeftAction(), side)), _switch_sign(side, ξ))
end

rigid_motions(G, vels, side) = map(vels) do (χ, v)
    return rigid_motion(G, χ, v, side)
end
