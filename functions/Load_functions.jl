function circ(h::Int64,H::UnitRange{Int64})::Int64
    (h==minimum(H)) ? idx = maximum(H) : idx = h-1
    return idx
end

function RMP0!(m::JuMP.Model,
              ms::ms_struct,
              mp::mp_struct,
              hc::hc_struct)::JuMP.Model

    # */ ------------ variables ------------ /* #
    @variable(m, f )
    @variable(m, xi[ms.P,ms.Ii] >= .0 )
    @variable(m, xo[1:hc.nx,ms.Io])
    # */ ----------------------------------- /* #

    # */ ----------- obj function ---------- /* #
    @objective(m, Min, f )
    # */ ----------------------------------- /* #

    # */ ----------- constraints ----------- /* #
    @constraint(m, c1, f >= sum(mp.ci[p,ii]*xi[p,ii] for p in ms.P for ii in ms.Ii) + sum(mp.cf[p]*xo[p,io] for p in ms.P for io in ms.Io))
    @constraint(m, c2[p=ms.P,io=ms.Io], xo[p,io] == mp.xh[io,p] + sum(xi[p,ii] for ii in mp.map[p,io]) )
    @constraint(m, c3[p=ms.P,io=ms.Io], xo[p,io] >= .0 )
    @constraint(m, c4[p=ms.P,io=ms.Io], xo[p,io] <= mp.xm[p] )
    @constraint(m, c6[       io=ms.Io], xo[hc.nx0+1,io] == hc.h[io,1] )
    # */ ----------------------------------- /* #

    return m
end

function LM!(m::JuMP.Model,
            ms::ms_struct,
            mp::mp_struct,
            hc::hc_struct)::JuMP.Model

    RMP0!(m,ms,mp,hc)
    # */ ------------ variables ------------ /* #
    @variable(m, β[ms.Io] >= 0., container=Array)
    @variable(m, ν[ms.Io] >= 0., container=Array)
    @variable(m, ϕ[ms.Io] >= 0., container=Array)
    # */ ----------- obj function ---------- /* #
    @objective(m, Min, objective_function(m) + sum(β[io] for io in ms.Io))
    # */ --------- constraints ---------- /* #
    @constraint(m, cl01[io=ms.Io], ν[io] + ϕ[io]*hc.c[io,1] <= β[io] )

    return m
end

function UM!(m::JuMP.Model,
            ms::ms_struct,
            mp::mp_struct,
            hc::hc_struct)::JuMP.Model

    RMP0!(m,ms,mp,hc)
    # */ ------------ variables ------------ /* #
    @variable(m, β[ms.Io] >= 0. )
    # */ ----------- obj function ---------- /* #
    @objective(m, Min, objective_function(m) + sum(β[io] for io in ms.Io))
    # */ --------- constraints ---------- /* #
    @constraint(m, cu01[io=ms.Io], 0. <= β[io] )
    @constraint(m, cu02[io=ms.Io,ix=1:hc.nx], 0. <= m[:xo][ix,io] )
    @constraint(m, cu03[io=ms.Io], 0. == 1. )

    return m
end

# */ ---------------------------------------------------------------------- /* #

function SP!(m::JuMP.Model,
            ps::ps_struct,
            pp::pp_struct,
            hc::hc_struct)::JuMP.Model

    # */ ------------ variables ------------ /* #
    @variable(m, ϕ[1:hc.nc] )
    @variable(m, c0 )
    @variable(m, yG[ps.G,ps.S,ps.H] >= 0  )
    @variable(m, yI[ps.B,ps.S,ps.H] >= 0  )
    @variable(m, yO[ps.B,ps.S,ps.H] >= 0  )
    @variable(m, yL[ps.B,ps.S,ps.H] >= 0  )
    @variable(m, yR[ps.R,ps.S,ps.H] >= 0  )
    @variable(m, yS[     ps.S,ps.H] >= 0  )
    @variable(m, x[1:hc.nx] )
    # */ ----------- obj function ---------- /* #
    @objective(m, Min, c0 )
    # */ --------- constraints ---------- /* #
    @constraint(m, c01, c0 == sum(sum(sum((pp.cvg[g]+pp.cfg[g]/pp.ηg[g])*yG[g,s,h] for g in ps.G) + pp.cs*yS[s,h] for h in ps.H) for s in ps.S))
    @constraint(m, c02, ϕ[1] == sum(sum(sum(pp.eg[g]/pp.ηg[g]*yG[g,s,h] for g in ps.G) for h in ps.H) for s in ps.S))
    @constraint(m, c03[b=ps.B,s=ps.S,h=ps.H], yL[b,s,h] - yL[b,s,circ(h,ps.H)] == pp.ηb[b]*yI[b,s,h] - yO[b,s,h])
    @constraint(m, c04[g=ps.G,s=ps.S,h=ps.H], yG[g,s,h] - yG[g,s,circ(h,ps.H)] <= pp.rg[g]*x[g])
    @constraint(m, c05[g=ps.G,s=ps.S,h=ps.H], yG[g,s,circ(h,ps.H)] - yG[g,s,h] <= pp.rg[g]*x[g])
    @constraint(m, c06[g=ps.G,s=ps.S,h=ps.H], yG[g,s,h] <= x[g])
    @constraint(m, c07[b=ps.B,s=ps.S,h=ps.H], yI[b,s,h] <= pp.pb[b]*x[ps.b0+b])
    @constraint(m, c08[b=ps.B,s=ps.S,h=ps.H], yO[b,s,h] <= pp.pb[b]*x[ps.b0+b])
    @constraint(m, c09[b=ps.B,s=ps.S,h=ps.H], yL[b,s,h] <= x[ps.b0+b])
    @constraint(m, c10[r=ps.R,s=ps.S,h=ps.H], yR[r,s,h] <= pp.pr[r,s,h]*x[ps.r0+r])
    @constraint(m, c11[s=ps.S,h=ps.H], sum(yG[g,s,h] for g in ps.G) + sum(yO[b,s,h]-yI[b,s,h] for b in ps.B) + sum(yR[r,s,h] for r in ps.R) + yS[s,h] == -x[hc.nx0+1]*pp.pd[s,h] )

    return m
end

function create_SPo(ps::ps_struct,
                    pp::pp_struct,
                    hc::hc_struct)::SPs

    gurobi_env = @suppress Gurobi.Env()
    m = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,gurobi_env,OutputFlag=0))
    SP!(m,ps,pp,hc)
    x,λ = zeros(hc.nx),zeros(hc.nx)
    c,ϕ = zeros(hc.nc),zeros(hc.nc)
    θ   = 0.
    return SPs(m,x,c,θ,λ,ϕ)

end

function update_and_solve_SPo!(spo::SPs,
                                 x::Array{Float64,1},
                                 c::Array{Float64,1})::SPs

    spo.x .= x
    spo.c .= c
    fix.(spo.m[:x],spo.x;force=true)
    @objective(spo.m, Min, spo.m[:c0] + sum(spo.c[ic]*spo.m[:ϕ][ic] for ic in 1:length(spo.c)))
    optimize!(spo.m)
    spo.θ = objective_value(spo.m)
    spo.λ .= dual.(FixRef.(spo.m[:x]))
    spo.ϕ .= value.(spo.m[:ϕ])

    return spo

end

# */ ---------------------------------------------------------------------- /* #

function create_RMPo(ms::ms_struct,
                     mp::mp_struct,
                     ps::ps_struct,
                     pp::pp_struct,
                     hc::hc_struct,
                     al::al_params)::RMPs

    gurobi_env = @suppress Gurobi.Env()
    rmL = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,gurobi_env,OutputFlag=0))
    rmU = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,gurobi_env,OutputFlag=0))
    rmB = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer,gurobi_env,OutputFlag=0))
    LM!(rmL,ms,mp,hc)
    UM!(rmU,ms,mp,hc)
    UM!(rmB,ms,mp,hc)
    rmLo = MLo(rmL,Array{ConstraintRef,2}(undef,hc.ny,al.J+1))
    rmUo = MUo(rmU,Array{VariableRef  ,2}(undef,hc.ny,al.J+1))
    rmBo = MUo(rmB,Array{VariableRef  ,2}(undef,hc.ny,al.J+1))
    sp = create_SPo(ps,pp,hc)
    xm,cm = vcat(minimum(mp.xh,dims=1)[:],[minimum(hc.h)]),[minimum(hc.c)]
    ls = LSs(zeros(ms.Io[end],ps.G[end]),zeros(ms.Io[end],ps.B[end]),zeros(ms.Io[end],ps.B[end]),zeros(ms.Io[end],ps.R[end]),zeros(ms.Io[end]))
    return RMPs(rmLo,rmUo,rmBo,sp,xm,cm,hc.nx,hc.nc,zeros(hc.nx,hc.ny),hc.h,hc.c,ms.Io,zeros(hc.ny),zeros(hc.ny),-Inf,Inf,Inf,zeros(0),zeros(0),zeros(0),0,al.ϵ,al.J,ls)

end

function build_cuts_RMP!(rmp::RMPs)::RMPs

    j = rmp.it+1
    for i in rmp.Io
        rmp.mL.c[i,j] = @constraint(rmp.mL.m, rmp.mL.m[:ν][i] + sum(rmp.mL.m[:ϕ][i,ic]*rmp.sp.c[ic] for ic in 1:rmp.nc) >= rmp.sp.θ + sum(rmp.sp.λ[ix]*(rmp.mL.m[:xo][ix,i]-rmp.sp.x[ix]) for ix in 1:rmp.nx) )
        rmp.mU.v[i,j] = @variable(rmp.mU.m,base_name="μ[$(i),$(j)]",lower_bound=0.)
        set_normalized_coefficient(rmp.mU.m[:cu01][i],rmp.mU.v[i,j],rmp.sp.θ + sum(rmp.sp.ϕ[ic]*(rmp.c[i,ic]-rmp.sp.c[ic]) for ic in 1:rmp.nc))
        for ix in 1:rmp.nx set_normalized_coefficient(rmp.mU.m[:cu02][i,ix],rmp.mU.v[i,j],rmp.sp.x[ix]) end
        set_normalized_coefficient(rmp.mU.m[:cu03][i],rmp.mU.v[i,j],1.)
        rmp.mB.v[i,j] = @variable(rmp.mB.m,base_name="μ[$(i),$(j)]",lower_bound=0.)
        set_normalized_coefficient(rmp.mB.m[:cu01][i],rmp.mB.v[i,j],rmp.sp.θ + sum(rmp.sp.ϕ[ic]*(rmp.c[i,ic]-rmp.sp.c[ic]) for ic in 1:rmp.nc))
        for ix in 1:rmp.nx set_normalized_coefficient(rmp.mB.m[:cu02][i,ix],rmp.mB.v[i,j],rmp.sp.x[ix]) end
        set_normalized_coefficient(rmp.mB.m[:cu03][i],rmp.mB.v[i,j],1.)
    end
    return rmp

end

function do_step0_RMP!(rmp::RMPs)::RMPs

    println("")
    println(" PlanModel algorithm started \r")
    println("")
    update_and_solve_SPo!(rmp.sp,rmp.xm,rmp.cm)
    build_cuts_RMP!(rmp)
    return rmp

end

function do_step_RMP!(rmp::RMPs)::RMPs

    rmp.it += 1
    optimize!(rmp.mL.m)
    rmp.x .= value.(rmp.mL.m[:xo])
    rmp.Lb = objective_value(rmp.mL.m)
    for ix in 1:rmp.nx for iy in rmp.Io fix(rmp.mB.m[:xo][ix,iy],rmp.x[ix,iy];force=true) end end
    optimize!(rmp.mB.m)
    rmp.βL .= value.(rmp.mL.m[:β])
    rmp.βU .= value.(rmp.mB.m[:β])
    i0 = findmax(rmp.βU-rmp.βL)[2]
    update_and_solve_SPo!(rmp.sp,rmp.x[:,i0],rmp.c[i0,:])
    build_cuts_RMP!(rmp)
    optimize!(rmp.mU.m)
    rmp.Ub = objective_value(rmp.mU.m)
    rmp.δ = (rmp.Ub-rmp.Lb)/rmp.Ub*100
    println(" j = $(rmp.it), Δ = $(@sprintf("%.2e",rmp.δ))% \r")
    push!(rmp.LB,rmp.Lb)
    push!(rmp.UB,rmp.Ub)
    push!(rmp.Δ,rmp.δ)
    return rmp

end

function solve_RMP!(rmp::RMPs)::RMPs

    do_step0_RMP!(rmp)
    for j in 1:rmp.J
        do_step_RMP!(rmp)
        (rmp.δ <= rmp.ϵ) ? break : nothing
    end
    return rmp

end

function last_step_RMP!(rmp::RMPs,
                         hc::hc_struct,
                         ps::ps_struct)::RMPs

    xo,co =zeros(hc.nx),zeros(hc.nc)
    for io in RMP.Io
        xo .= value.(RMP.mU.m[:xo][:,io])
        co .= RMP.c[io]
        update_and_solve_SPo!(RMP.sp,xo,co)
        for g in ps.G
            RMP.ls.G[io,g] = sum(value.(RMP.sp.m[:yG][g,:,:]))/1.0e6
        end
        for b in ps.B
            RMP.ls.Bi[io,b] = sum(value.(RMP.sp.m[:yI][b,:,:]))/1.0e6
            RMP.ls.Bo[io,b] = sum(value.(RMP.sp.m[:yO][b,:,:]))/1.0e6
        end
        for r in ps.R
            RMP.ls.R[io,r] = sum(value.(RMP.sp.m[:yR][r,:,:]))/1.0e6
        end
        RMP.ls.E[io] = value(RMP.sp.m[:ϕ][1])/1.0e6
    end
    return rmp

end

# */ ---------------------------------------------------------------------- /* #

function write_output(rmp::RMPs,
                       ms::ms_struct,
                       ps::ps_struct,
                       sd::sd_struct)

    XLSX.openxlsx("$(pwd())/output/results.xlsx", mode="w") do xf

        sheet = xf[1]
        XLSX.rename!(sheet, "Investment")
        sheet["A1"] = "Investment in new capacity"
        sheet["A3"] = "Years"
        sheet["A5", dim=1] = sd.yrs
        sheet["B3", dim=2] = sd.techs
        sheet["B4", dim=2] = repeat(["(MW)"],length(sd.techs))
        xi = zeros(ms.P[end],ms.Ii[end])
        xi .= value.(rmp.mU.m[:xi])
        iO = (sd.bi.*(1:length(sd.bi)))[Bool.(sd.bi)]
        Xi = zeros(ms.Io[end],ms.P[end])
        for ii in ms.Ii Xi[iO[ii],:] .= xi[:,ii] end
        sheet["B5"] = Xi
        sheet[3,ms.P[end]+3] = "Total_cost"
        sheet[4,ms.P[end]+3] = "(10^9 £)"
        sheet[5,ms.P[end]+3] = objective_value(rmp.mU.m)*1.0e-9
        sheet[3,ms.P[end]+4] = "Investment_cost"
        sheet[4,ms.P[end]+4] = "(10^9 £)"
        sheet[5,ms.P[end]+4] = sum(mp.ci[p,ii]*xi[p,ii] for p in ms.P for ii in ms.Ii)*1.0e-9

        XLSX.addsheet!(xf)
        sheet = xf[2]
        XLSX.rename!(sheet, "Operation")
        sheet["A1"] = "Total yearly available capacity"
        sheet["A3"] = "Years"
        sheet["A5", dim=1] = sd.yrs
        sheet["B3", dim=2] = sd.techs
        sheet["B4", dim=2] = repeat(["(MW)"],length(sd.techs))
        Xo = zeros(ms.Io[end],ms.P[end])
        for p in ms.P for io in ms.Io Xo[io,p] = value(rmp.mU.m[:xo][p,io])  end end
        sheet["B5"] = Xo
        sheet[3,ms.P[end]+3] = "FixO&M"
        sheet[4,ms.P[end]+3] = "(10^9 £)"
        sheet[5,ms.P[end]+3] = sum(mp.cf[p]*Xo[io,p] for p in ms.P for io in ms.Io)*1.0e-9
        sheet[3,ms.P[end]+4] = "Operation cost"
        sheet[4,ms.P[end]+4] = "(10^9 £)"
        sheet[5,ms.P[end]+4] = sum(value.(rmp.mU.m[:β]))*1.0e-9

        XLSX.addsheet!(xf)
        sheet = xf[3]
        XLSX.rename!(sheet, "Usage_Conventional")
        sheet["A1"] = "Total yearly usage of conventional technologies"
        sheet["A3"] = "Years"
        sheet["A5", dim=1] = sd.yrs
        sheet["B3", dim=2] = sd.techs[ps.G]
        sheet["B4", dim=2] = repeat(["(TWh)"],ps.G[end])
        sheet["B5"] = RMP.ls.G

        XLSX.addsheet!(xf)
        sheet = xf[4]
        XLSX.rename!(sheet, "Usage_Storage")
        sheet["A1"] = "Total yearly usage of storage technologies"
        sheet["A3"] = "Years"
        sheet["A5", dim=1] = sd.yrs
        sheet["B3", dim=2] = repeat(sd.techs[ps.b0.+ps.B],2) .* vcat(repeat([" (charge)"],ps.B[end]),repeat([" (discharge)"],ps.B[end]))
        sheet["B4", dim=2] = repeat(["(TWh)"],2*ps.B[end])
        sheet[XLSX.encode_column_number(2).*"5"]           = RMP.ls.Bi
        sheet[XLSX.encode_column_number(2+ps.B[end]).*"5"] = RMP.ls.Bo

        XLSX.addsheet!(xf)
        sheet = xf[5]
        XLSX.rename!(sheet, "Usage_Renewable")
        sheet["A1"] = "Total yearly usage of renewable technologies"
        sheet["A3"] = "Years"
        sheet["A5", dim=1] = sd.yrs
        sheet["B3", dim=2] = sd.techs[ps.r0.+ps.R]
        sheet["B4", dim=2] = repeat(["(TWh)"],ps.R[end])
        sheet["B5"] = RMP.ls.R

        XLSX.addsheet!(xf)
        sheet = xf[6]
        XLSX.rename!(sheet, "Emissions")
        sheet["A1"] = "Total yearly CO2 emissions"
        sheet["A3"] = "Years"
        sheet["A5", dim=1] = sd.yrs
        sheet["B3"] = "CO2 emissions"
        sheet["B4"] = "(MtCO2)"
        sheet["B5", dim=1] = RMP.ls.E

    end

end

# */ ---------------------------------------------------------------------- /* #
