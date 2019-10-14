cd(dirname(@__FILE__))
include("./functions/Load_all_stuff.jl")

RMP = create_RMPo(ms,mp,ps,pp,hc,al)
solve_RMP!(RMP)
write_output(RMP,ms)
println("")
println(" finished - optimal solution found")
# */------------------------------------------------------------------------/* #
