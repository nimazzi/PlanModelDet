module O_structs

export ms_struct,mp_struct,pp_struct,ps_struct,hc_struct,sd_struct,al_params

mutable struct ms_struct{UI<:UnitRange{Int64}}
      P::UI
     Ii::UI
     Io::UI
end

mutable struct mp_struct{AF1<:Array{Float64,1},AF2<:Array{Float64,2}}
      xh::AF2
      xm::AF1
      ci::AF2
      cf::AF1
     map::Dict{Tuple{Int64,Int64},Array{Int64,1}}
end


mutable struct pp_struct{AF1<:Array{Float64,1},AF2<:Array{Float64,2},AF3<:Array{Float64,3},F<:Float64}
   cvg::AF1
   cfg::AF1
    eg::AF1
    ηg::AF1
    rg::AF1
    ηb::AF1
    pb::AF1
    pr::AF3
    pd::AF2
    cs::F
end

mutable struct ps_struct{UI<:UnitRange{Int64},I<:Int64}
    S::UI
    H::UI
    P::UI
    G::UI
    B::UI
    R::UI
   b0::I
   r0::I
end

mutable struct hc_struct{I<:Int64,AF2<:Array{Float64,2}}
    h::AF2
    c::AF2
   nx::I
   nc::I
  nx0::I
   ny::I
end

mutable struct al_params
   ϵ::Float64
   J::Int64
end

mutable struct sd_struct{AI1<:Array{Int64,1}}
   yrs::AI1
   yri::AI1
    bi::AI1
 techs::Array{String,1}
end

end
