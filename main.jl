cd(dirname(@__FILE__))
include("./functions/Load_all_stuff.jl")

# */------------------------------------------------------------------------/* #
RMP = create_RMPo(ms,mp,ps,pp,hc,al)
solve_RMP!(RMP)
last_step_RMP!(RMP,hc,ps)
write_output(RMP,ms,ps,sd)
println("")
println(" finished - optimal solution found")
# */------------------------------------------------------------------------/* #
