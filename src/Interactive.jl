module Interactive

using Pkg
Pkg.add("DifferentialEquations")
Pkg.add("Interact")
Pkg.add("Blink")
Pkg.add("Plots")

using DifferentialEquations
using Interact
using Blink
using Plots


function seirhcd_ode(du,u,p,t)
    S,E,I,H,C,R,D = u
    α,β,γ,δ,η,λ,p_a,p_c,p_f = p
    du[1] = -β*S*I+δ*R
    du[2] = β*S*I-α*E
    du[3] = α*E-γ*I
    du[4] = γ*(1-p_a)*I+λ*(1-p_f)*C-η*H
    du[5] = η*p_c*H-λ*C
    du[6] = γ*p_a*I+η*(1-p_c)*H-δ*R
    du[7] = λ*p_f*C
end
parms = [0.3,0.05,0.01,0.15,0.5,0.5,0.79,0.12,0.33]
init = [0.9,0.1,0,0.0,0,0,0]
tspan = (0.0,1000.0)
seirhcd_prob = ODEProblem(seirhcd_ode,init,tspan,parms)
seirhcd_sol = solve(seirhcd_prob,saveat = 0.1)
p = plot(seirhcd_sol[1,:], title = " Simulation of SEIHCRD-Model",label = ["S" "E" "I" "H" "C" "R" "D"],linewidth = 2,legend = :outertopright,dpi = 200)


function solve_seirhcd(α,β,γ,δ,η,λ,p_a,p_c,p_f,b)
    parms = [α,β,γ,δ,η,λ, p_a,p_c,p_f]
    init = [0.9,0.1,0.0,0.0,0,0,0]
    tspan = (0.0,1000.0)
    seirhcd_prob = ODEProblem(seirhcd_ode,init,tspan,parms)
    seirhcd_sol = solve(seirhcd_prob,saveat = 0.1)
    f = findall(x -> x == 1, b)
    return transpose(seirhcd_sol[f,:])
end

function interactive_simulation()
    update = button("Update")
    # reset = button("Reset")
    s_alpha = slider(0:0.01:1, label = "α", value = 0.3)
    s_beta = slider(0:0.01:1, label = "β", value = 0.05)
    s_gamma= slider(0:0.01:1, label = "γ", value = 0.01)
    s_delta = slider(0:0.01:1, label = "δ", value = 0.0)
    s_eta = slider(0:0.01:1, label = "η", value = 0.5)
    s_lambda = slider(0:0.01:1, label = "λ", value = 0.5)
    s_pa = slider(0:0.01:1, label = "p_a", value = 0.79)
    s_pc = slider(0:0.01:1, label = "p_c", value = 0.12)
    s_pf = slider(0:0.01:1, label = "p_f", value = 0.33)
    S = toggle(label = "S")
    S[] = true
    E = toggle(label = "E")
    E[] = true
    I = toggle(label = "I")
    I[] = true
    H = toggle(label = "H")
    H[] = true
    C = toggle(label = "C")
    C[] = true
    R = toggle(label = "R")
    R[] = true
    D = toggle(label = "D")
    D[] = true
    l = ["S" "E" "I" "H" "C" "R" "D"]
    output = Interact.@map (&update; solve_seirhcd(s_alpha[],s_beta[],s_gamma[],s_delta[],s_eta[],s_lambda[],s_pa[],s_pc[],s_pf[],[S[],E[],I[],H[],C[],R[],D[]]))
    plt = Interact.@map plot(&output, size=(750,750),title = " Simulation of SEIHCRD-Model",label = permutedims(l[findall(x -> x == 1, [S[],E[],I[],H[],C[],R[],D[]])]),linewidth = 2,legend = :outertopright,dpi = 200)
    wdg = Widget(["update" => update, "s_alpha" => s_alpha, "s_beta" => s_beta, "s_gamma" => s_gamma, "s_delta" => s_delta,"s_eta" => s_eta, "s_lambda" => s_lambda,"s_pa" => s_pa,"s_pc" => s_pc, "s_pf" => s_pf, "S" => S, "E" => E, "I" => I, "H" => H, "C" => C, "R" => R, "D" => D ],output = output)
    @layout! wdg hbox(plt, vbox(:s_alpha, :s_beta, :s_gamma, :s_delta,:s_eta,:s_lambda,:s_pa,:s_pc,:s_pf, :update), vbox(:S,:E,:I,:H,:C,:R,:D))
    w = Window()
    body!(w, wdg)
end

interactive_simulation()

end
