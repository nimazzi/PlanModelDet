function load_data()::Tuple{ms_struct,mp_struct,ps_struct,pp_struct,hc_struct,sd_struct,al_params}
    XF = XLSX.readxlsx("$(pwd())/input/data.xlsx")

    # */ --- Conventional --------------------------------------------------*/ #
    GDF = convert(Array{Float64,2},XF["Conventional!B3:F8"])
    cvg = GDF[:,1]
    cfg = GDF[:,2]
    eg  = GDF[:,3]
    ηg  = GDF[:,4]
    rg  = GDF[:,5]
    # */ -------------------------------------------------------------------*/ #

    # */ --- Storage -------------------------------------------------------*/ #
    BDF = convert(Array{Float64,2},XF["Storage!B3:C5"])
    ηb = BDF[:,1]
    pb = BDF[:,2]
    # */ -------------------------------------------------------------------*/ #

    # */ --- Renewables ----------------------------------------------------*/ #
    RDF = convert(Array{Float64,2},XF["Renewables!B3:M2192"])
    pr = zeros(3,4,2190)
    for r in 1:3 for s in 1:4 pr[r,s,:] = RDF[:,(r-1)*4+s] end end
    # */ -------------------------------------------------------------------*/ #

    # */ --- Renewables ----------------------------------------------------*/ #
    DDF = convert(Array{Float64,2},XF["Demand!B2:E2191"])
    pd = zeros(4,2190)
    for s in 1:4 pd[s,:] = DDF[:,s] end
    # */ -------------------------------------------------------------------*/ #

    # */ --- Others --------------------------------------------------------*/ #
    ODF = convert(Array{Float64,2},XF["Others!A2:A2"])
    cs = ODF[1,1]
    # */ -------------------------------------------------------------------*/ #

    # */ --- Investment ----------------------------------------------------*/ #
    IDF = convert(Array{Float64,2},XF["Investment!B3:F14"])
    ci = IDF[:,1]
    cf = IDF[:,2]
    lf = convert(Array{Int64,1},IDF[:,3])
    tc = convert(Array{Int64,1},IDF[:,4])
    xm = IDF[:,5]
    # */ -------------------------------------------------------------------*/ #

    # */ --- TimeLine ------------------------------------------------------*/ #
    ny = length(XF["TimelineData!A:A"])-2
    yrs = convert(Array{Int64,2},XF["TimelineData!A3:A$(ny+2)"])[:]
    TDF1 = convert(Array{Float64,2},XF["TimelineData!B3:D$(ny+2)"])
    νD  = TDF1[:,1]
    co  = TDF1[:,2]
    ib  = convert(Array{Int64,1},TDF1[:,3])
    yri = yrs[Bool.(ib)]
    TDF2 = convert(Array{Float64,2},XF["TimelineData!F3:Q$(ny+2)"])
    xh = TDF2[:,:]
    # */ -------------------------------------------------------------------*/ #

    # */ --- ms_struct -----------------------------------------------------*/ #
    P  = 1:length(ci)
    Ii = 1:sum(ib)
    Io = 1:length(ib)

    Ms = ms_struct(P,Ii,Io)
    # */ -------------------------------------------------------------------*/ #

    # */ --- mp_struct -----------------------------------------------------*/ #
    cinv = zeros(P[end],Ii[end])
    for p in P for i in Ii
        mlt = yrs[end]-(yri[i]+tc[p])+1
        (mlt >= lf[p]) ? cinv[p,i] = ci[p] : cinv[p,i] = ci[p]*(mlt/lf[p])
    end end
    map = Dict{Tuple{Int64,Int64},Array{Int64,1}}()
    for p in P for io in Io
        ayoi = zeros(Int64,0)
        for ii in Ii
            (tc[p] <= yrs[io]-yri[ii] <= lf[p]) ? push!(ayoi,ii) : nothing
        end
        map[p,io] = ayoi
    end end

    Mp = mp_struct(xh,xm,cinv,cf,map)
    # */ -------------------------------------------------------------------*/ #

    # */ --- ps_struct -----------------------------------------------------*/ #
    S = 1:4
    H = 1:2190
    G = 1:6
    B = 1:3
    R = 1:3
    b0 = 6
    r0 = 9

    Ps = ps_struct(S,H,P,G,B,R,b0,r0)
    # */ -------------------------------------------------------------------*/ #

    # */ --- pp_struct -----------------------------------------------------*/ #

    Pp = pp_struct(cvg,cfg,eg,ηg,rg,ηb,pb,pr,pd,cs)
    # */ -------------------------------------------------------------------*/ #

    # */ --- hc_struct -----------------------------------------------------*/ #
    h,c = zeros((ny,1)),zeros((ny,1))
    h[:,1] = -νD
    c[:,1] = co
    nx0 = P[end]

    HC = hc_struct(h,c,nx0+1,1,nx0,ny)
    # */ -------------------------------------------------------------------*/ #

    # */ --- sd_struct -----------------------------------------------------*/ #
    ny = length(XF["TimelineData!A:A"])-2
    yrs = convert(Array{Int64,2},XF["TimelineData!A3:A$(ny+2)"])[:]
    TDF1 = convert(Array{Float64,2},XF["TimelineData!B3:D$(ny+2)"])
    ib  = convert(Array{Int64,1},TDF1[:,3])
    yri = yrs[Bool.(ib)]
    techs = convert(Array{String,1},XF["TimelineData!F1:Q1"][:])
    SD = sd_struct(yrs,yri,ib,techs)
    # */ -------------------------------------------------------------------*/ #

    # */ --- al_params -----------------------------------------------------*/ #
    prms = XF["Solution_Params!C1:C2"]
    J = Int64(prms[1])
    ϵ = Float64(prms[2])
    AL = al_params(ϵ,J)
    # */ -------------------------------------------------------------------*/ #

    return Ms,Mp,Ps,Pp,HC,SD,AL

end

ms,mp,ps,pp,hc,sd,al = load_data()
