struct AgentBasedModel{A, S, F, P}
    agents::A
    space::S
    scheduler::F
    properties::P
end

function AgentBasedModel(a::Vector{A},s::S) where {A<:AbstractAgentM, S<:AbstractSpacesTuple}
    AgentBasedModel(a,s,nothing, nothing)
end

function AgentBasedModel(N <:Int, Max <: Int ,s::S) where {A<:AbstractAgentM, S<:AbstractSpacesTuple}
    # To be implemented
end
