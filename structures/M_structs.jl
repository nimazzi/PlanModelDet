module M_structs

using JuMP

export SPs,MLo,MUo,RMPs,LSs

mutable struct SPs{AF1<:Array{Float64,1}}
    m::JuMP.Model
    x::AF1
    c::AF1
    θ::Float64
    λ::AF1
    ϕ::AF1
end

mutable struct MLo
    m::JuMP.Model
    c::Array{JuMP.ConstraintRef,2}
end

mutable struct MUo
    m::JuMP.Model
    v::Array{JuMP.VariableRef,2}
end

mutable struct LSs{AF1<:Array{Float64,1},AF2<:Array{Float64,2}}
    G::AF2
   Bi::AF2
   Bo::AF2
    R::AF2
    E::AF1
end

mutable struct RMPs{I<:Int64,F<:Float64,AF1<:Array{Float64,1},AF2<:Array{Float64,2}}
   mL::MLo
   mU::MUo
   mB::MUo
   sp::SPs
   xm::AF1
   cm::AF1
   nx::I
   nc::I
    x::AF2
    h::AF2
    c::AF2
   Io::UnitRange{Int64}
   βL::AF1
   βU::AF1
   Lb::F
   Ub::F
    δ::F
   LB::AF1
   UB::AF1
    Δ::AF1
   it::I
    ϵ::F
    J::I
   ls::LSs
end


end
