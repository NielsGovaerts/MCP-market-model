using JuMP
using ExcelReaders
using DataFrames
using Iterators
using Complementarity
using Gadfly


function to_dict(data, sets)
    d = Dict()
    for (n, idx) in enumerate(product(sets...))
        d[idx...] = D_data[n]
    end
    return d
end

function next_idx(index, idx)
    pos = findin(index,[idx])
    return index[pos+1][1]

end

function prev_idx(index, idx)
    pos = findin(index,[idx])
    return index[pos-1][1]

end




### Data
################################################################################

tar = 0.15


xls_file = joinpath(@__DIR__, "Data_1.xlsx")

data_1 = readxlsheet(DataFrame, xls_file, "Demand")
data_2 = readxlsheet(DataFrame, xls_file, "Tech")
data_3 = readxlsheet(DataFrame, xls_file, "Users")
data_4 = readxlsheet(DataFrame, xls_file, "Plants")
data_5 = readxlsheet(DataFrame, xls_file, "Res")

D_data = [data_1[:Demand_1] data_1[:Demand_2]]
Di_data = data_1[:Demand_i]
lambda_data = data_1[:lambda]

dr = data_2[:Discount_rate][1]
pvi = data_2[:PV_installation][1]
pvl = data_2[:PV_lifetime][1]
si = data_2[:S_installation][1]
sl = data_2[:S_lifetime][1]
K = data_2[:K][1]
eff_CH = data_2[:eff][1]
eff_DC = data_2[:eff][1]
ICW_year = data_2[:ICW_year][1]


share_users_data = data_3[:Share_users]
numb_users_data = data_3[:Numb_users]
PVi_max_data = data_3[:PV_I_max]
Si_max_data = data_3[:S_I_max]

P_max_data = data_4[:Pmax]
AVC_data = data_4[:AVC]
ICP_year_data = data_4[:ICP_year]


LFS_data = data_5[:LFS]
LFW_data = data_5[:LFW]

numb_timesteps = length(Di_data)
numb_user_types = length(share_users_data)
numb_plants = length(AVC_data)

C_PV_annual = pvi/((1-(1/(1+dr))^pvl)/dr)
C_PV = C_PV_annual/365
C_S_annual = si/((1-(1/(1+dr))^sl)/dr)
C_S = C_S_annual/365

ICW = ICW_year/365
ICP_data = ICP_year_data/365


timesteps = [@sprintf "t%03d" t for t in 1:numb_timesteps]
user_types = [@sprintf "i%02d" i for i in 1:numb_user_types]
plants = [@sprintf "j%02d" j for j in 1:numb_plants]

D = Dict()
Di = Dict()
lambda_e = Dict()
share_users = Dict()
numb_users = Dict()
PVi_max = Dict()
Si_max = Dict()
P_max = Dict()
AVC = Dict()
ICP = Dict()
LFS = Dict()
LFW = Dict()
for t in 1:length(timesteps), i in 1:length(user_types), j in 1:length(plants)
    D[timesteps[t], user_types[i]] = D_data[t,i]
    Di[timesteps[t]] = Di_data[t]
    lambda_e[timesteps[t]] = lambda_data[t]
    share_users[user_types[i]] = share_users_data[i]
    numb_users[user_types[i]] = numb_users_data[i]
    PVi_max[user_types[i]] = PVi_max_data[i]
    Si_max[user_types[i]] = Si_max_data[i]
    P_max[plants[j]] = P_max_data[j]
    AVC[plants[j]] = AVC_data[j]
    ICP[plants[j]] = ICP_data[j]
    LFS[timesteps[t]] = LFS_data[t]
    LFW[timesteps[t]] = LFW_data[t]
end

# D = to_dict(D_data, [timesteps, user_types])












### MCP model
################################################################################

m = MCPModel()

# Primale variabelen
@variable(m, W[t in timesteps, i in user_types])
@variable(m, E[t in timesteps, i in user_types] >= 0)
@variable(m, CH[t in timesteps, i in user_types] >= 0)
@variable(m, DC[t in timesteps, i in user_types] >= 0)
@variable(m, PVi[i in user_types] >= 0)
@variable(m, Si[i in user_types] >= 0)
@variable(m, Wtp[i in user_types] >= 0)

@variable(m, P[t in timesteps, j in plants] >= 0)
@variable(m, PW[t in timesteps] >= 0)
@variable(m, CP[j in plants] >= 0)
@variable(m, CW >= 0)

# Duale variabelen
@variable(m, alpha1[t in timesteps, i in user_types] >= 0)
@variable(m, alpha2[t in timesteps, i in user_types] >= 0)
@variable(m, beta1[t in timesteps, i in user_types] >= 0)
@variable(m, beta2[t in timesteps, i in user_types] >= 0)
@variable(m, beta_2_1[i in user_types] >= 0)
@variable(m, beta_2_2[i in user_types] >= 0)
@variable(m, gamma[t in timesteps, i in user_types] >= 0)
@variable(m, delta[t in timesteps, i in user_types] >= 0)
@variable(m, epsilon[t in timesteps, i in user_types] >= 0)
@variable(m, zeta[i in user_types] >= 0)
@variable(m, eta[i in user_types] >= 0)
@variable(m, theta[i in user_types] >= 0)


@variable(m, mu[t in timesteps, j in plants] >= 0)
@variable(m, nu[t in timesteps] >= 0)

@variable(m, lambda1[t in timesteps] >= 0)
@variable(m, lambda2[t in timesteps] >= 0)


@NLexpression(m, alpha[t in timesteps, i in user_types], alpha1[t,i] - alpha2[t,i])
@NLexpression(m, beta[t in timesteps, i in user_types], beta1[t,i] - beta2[t,i])
@NLexpression(m, beta_2[i in user_types], beta_2_1[i] - beta_2_2[i])
@NLexpression(m, lambda[t in timesteps], lambda1[t] - lambda2[t])


### Consumers
################################################################################

# Gelijkheidsbeperkingen
@mapping(m, pb1[t in timesteps, i in user_types],
                W[t,i] + DC[t,i] + LFS[t]*PVi[i] - D[t,i] - CH[t,i])

@mapping(m, pb2[t in timesteps, i in user_types],
                -W[t,i] - DC[t,i] - LFS[t]*PVi[i] + D[t,i] + CH[t,i])

@mapping(m, soc1[t in timesteps[2:end],i in user_types],
                E[t,i] - E[prev_idx(timesteps,t),i]
                - CH[t,i]*eff_CH + DC[t,i]/eff_DC)

@mapping(m, soc2[t in timesteps[2:end],i in user_types],
                -E[t,i] + E[prev_idx(timesteps,t),i]
                + CH[t,i]*eff_CH - DC[t,i]/eff_DC)

@mapping(m, soc01[i in user_types],
                E[timesteps[1],i] - Si[i]/2 - CH[timesteps[1],i]*eff_CH
                + DC[timesteps[1],i]/eff_DC)

@mapping(m, soc02[i in user_types],
                -E[timesteps[1],i] + Si[i]/2 + CH[timesteps[1],i]*eff_CH
                - DC[timesteps[1],i]/eff_DC)

@mapping(m, socn1[i in user_types], E[timesteps[end],i] - Si[i]/2)

@mapping(m, socn2[i in user_types], -E[timesteps[end],i] + Si[i]/2)

# Ongelijkheidsbeperkingen
@mapping(m, e_max[t in timesteps,i in user_types], Si[i] - E[t,i])

@mapping(m, ch_max[t in timesteps,i in user_types], K*Si[i] - CH[t,i])

@mapping(m, dc_max[t in timesteps,i in user_types], K*Si[i] - DC[t,i])

@mapping(m, si_max[i in user_types], Si_max[i] - Si[i])

@mapping(m, pvi_max[i in user_types],
    PVi_max[i] - PVi[i])

@mapping(m, wtp[i in user_types], Wtp[i] - sum(W[t,i] for t in timesteps))

# Stationariteitsvoorwaarden
@mapping(m, l_w[t in timesteps,i in user_types],
                lambda[t] - alpha[t,i] + theta[i])

@mapping(m, l_e[t in timesteps[1:end-1],i in user_types],
                -beta[t,i] + beta[next_idx(timesteps,t),i] + gamma[t,i])

@mapping(m, l_en[i in user_types],
                -beta[timesteps[end],i] - beta_2[i] + gamma[timesteps[end],i])

@mapping(m, l_ch[t in timesteps,i in user_types],
                alpha[t,i] + beta[t,i]*eff_CH + delta[t,i])

@mapping(m, l_dc[t in timesteps,i in user_types],
                -alpha[t,i] - beta[t,i]/eff_DC + epsilon[t,i])

@mapping(m, l_pvi[i in user_types],
                C_PV - sum(alpha[t,i]*LFS[t] for t in timesteps) + zeta[i])

@mapping(m, l_si[i in user_types],
                C_S - sum(gamma[t,i] + delta[t,i]*K
                + epsilon[t,i]*K for t in timesteps)
                + eta[i] + beta[timesteps[1],i]/2 + beta_2[i]/2)

@mapping(m, l_wtp[i in user_types], tar - theta[i])




### Producers
################################################################################

# Ongelijkheidsbeperkingen
@mapping(m, p_max[t in timesteps,j in plants], CP[j] - P[t,j])

@mapping(m, pw_max[t in timesteps], CW*LFW[t] - PW[t])

# Stationariteitsvoorwaarden
@mapping(m, l_p[t in timesteps, j in plants], AVC[j] - lambda[t] + mu[t,j])

@mapping(m, l_pw[t in timesteps], -lambda[t] + nu[t])

@mapping(m, l_cp[j in plants], ICP[j] - sum(mu[t,j] for t in timesteps))

@mapping(m, l_cw, ICW - sum(nu[t]*LFW[t] for t in timesteps))





### Market clearing
################################################################################

@mapping(m, mc1[t in timesteps],
            sum(P[t,j] for j in plants) + PW[t]
            - Di[t] - sum(numb_users[i]*W[t,i] for i in user_types))

@mapping(m, mc2[t in timesteps],
            -sum(P[t,j] for j in plants) - PW[t]
            + Di[t] + sum(numb_users[i]*W[t,i] for i in user_types))



### Complementarity
################################################################################

# Consumers
@complementarity(m,pb1,alpha1)
@complementarity(m,pb2,alpha2)
for i in user_types
    @complementarity(m,soc01[i],beta1[timesteps[1],i])
    @complementarity(m,soc02[i],beta2[timesteps[1],i])
end
for i in user_types, t in timesteps[2:end]
    @complementarity(m,soc1[t,i],beta1[t,i])
    @complementarity(m,soc2[t,i],beta2[t,i])
end
@complementarity(m,socn1,beta_2_1)
@complementarity(m,socn2,beta_2_2)
@complementarity(m,e_max,gamma)
@complementarity(m,ch_max,delta)
@complementarity(m,dc_max,epsilon)
@complementarity(m,pvi_max,zeta)
@complementarity(m,si_max,eta)
@complementarity(m,wtp,theta)
@complementarity(m,l_w,W)

for i in user_types, t in timesteps[1:end-1]
    @complementarity(m,l_e[t,i],E[t,i])
end

for i in user_types
    @complementarity(m,l_en[i],E[timesteps[end],i])
end

@complementarity(m,l_ch,CH)
@complementarity(m,l_dc,DC)
@complementarity(m,l_si,Si)
@complementarity(m,l_pvi,PVi)
@complementarity(m,l_wtp,Wtp)



# Producers
@complementarity(m,p_max,mu)
@complementarity(m,pw_max,nu)
@complementarity(m,l_p,P)
@complementarity(m,l_pw,PW)
@complementarity(m,l_cp,CP)
@complementarity(m,l_cw,CW)



# Market clearing
@complementarity(m,mc1,lambda1)
@complementarity(m,mc2,lambda2)














































### Solve
################################################################################

# PATHSolver.options(
#         convergence_tolerance=1e-8, output=:yes, time_limit=3600,
#         cumulative_iteration_limit = 50000,
#         lemke_rank_deficiency_iterations = 1000,
#         major_iteration_limit = 1000, minor_iteration_limit = 1000,
#         proximal_perturbation = 0, nms_initial_reference_factor = 10,
#         nms_memory_size = 2, nms_mstep_frequency = 1, lemke_start_type = "slack",
#         crash_perturb = "no", crash_method = "none")

PATHSolver.options(
        convergence_tolerance=1e-4, output=:yes, time_limit=3600,
        cumulative_iteration_limit = 500000,
        lemke_rank_deficiency_iterations = 100,
        major_iteration_limit = 10000, minor_iteration_limit = 10000)



status = solveMCP(m)

@show status



### Results
################################################################################

@show getvalue(Si)
@show getvalue(PVi)


Dist_cost_PV = getvalue(tar*Wtp[user_types[1]])*365
Dist_cost_active = getvalue(tar*Wtp[user_types[2]])*365
Costs_rec_day = getvalue(tar*sum(share_users[i]*Wtp[i] for i in user_types))
Costs_rec_year = 365*Costs_rec_day

lambda_r = zeros(numb_timesteps)
W_r = zeros(numb_timesteps,numb_user_types)
CH_r = zeros(numb_timesteps, numb_user_types)
DC_r = zeros(numb_timesteps, numb_user_types)
E_r = zeros(numb_timesteps, numb_user_types)
P_r = zeros(numb_timesteps, numb_plants)
Wtp_r = zeros(numb_user_types)
PVi_r = zeros(numb_user_types)
Si_r = zeros(numb_user_types)
CP_r = zeros(numb_plants)
for t in 1:numb_timesteps, i in 1:numb_user_types, j in 1:numb_plants
    lambda_r[t] = getvalue(lambda[timesteps[t]])
    W_r[t,i] = getvalue(W[timesteps[t], user_types[i]])
    CH_r[t,i] = getvalue(CH[timesteps[t], user_types[i]])
    DC_r[t,i] = getvalue(DC[timesteps[t], user_types[i]])
    E_r[t,i] = getvalue(E[timesteps[t], user_types[i]])
    P_r[t,j] = getvalue(P[timesteps[t], plants[j]])
    Wtp_r[i] = getvalue(Wtp[user_types[i]])
    PVi_r[i] = getvalue(PVi[user_types[i]])
    Si_r[i] = getvalue(Si[user_types[i]])
    CP_r[j] = getvalue(CP[plants[j]])
end


d_cons_tot = zeros(numb_timesteps)
for t = 1:numb_timesteps
    d_cons_tot[t] = sum(numb_users[i]*D[timesteps[t],i] for i in user_types)
end

d_tot = zeros(numb_timesteps)
for t = 1:numb_timesteps
    d_tot[t] = d_cons_tot[t] + Di[timesteps[t]]
end

p_plants = zeros(numb_timesteps)
for t=1:numb_timesteps
    p_plants[t] = sum(P_r[t,j] for j=1:numb_plants)
end

w_cons_tot = zeros(numb_timesteps)
for t = 1:numb_timesteps
    w_cons_tot[t] = getvalue(sum(numb_users[i]*(
    D[timesteps[t],i] - PVi[i]*LFS[timesteps[t]] + CH[timesteps[t],i]
    - DC[timesteps[t],i]) for i in user_types))
end

d_pros = zeros(numb_timesteps)
for t = 1:numb_timesteps
    d_pros[t] = numb_users[user_types[1]]*D[timesteps[t],user_types[1]]
end

w_pros = zeros(numb_timesteps)
for t = 1:numb_timesteps
    w_pros[t] = getvalue(numb_users[user_types[1]]*(
    D[timesteps[t],user_types[1]] - PVi[user_types[1]]*LFS[timesteps[t]]))
end

production_cons = zeros(numb_timesteps)
for t = 1:numb_timesteps
    production_cons[t] =
    getvalue(sum(numb_users[i]*PVi[i]*LFS[timesteps[t]] for i in user_types))
end

system_costs = getvalue(sum(sum(AVC[j]*P[t,j] for t in timesteps)
            + ICP[j]*CP[j] for j in plants) + ICW*CW)

producer_surplus = getvalue(sum(
        sum((lambda_r[t] - AVC[j])*P[timesteps[t],j] for t=1:numb_timesteps)
        + ICP[j]*CP[j] for j in plants))

consumer_surplus = sum(numb_users_data[i]*sum(                                  # VOLL is 1000 €/MWh
            1*D_data[t,i] - lambda_r[t]*W_r[t,i] for t=1:numb_timesteps)
            - tar*Wtp_r[i] - C_PV*PVi_r[i] - C_S*Si_r[i] for i=1:numb_user_types)

active_surplus = numb_users_data[1]*sum(
            1*D_data[t,1] - lambda_r[t]*W_r[t,1] for t=1:numb_timesteps)
            - tar*Wtp_r[1] - C_PV*PVi_r[1] - C_S*Si_r[1]

passive_surplus = numb_users_data[2]*sum(
                        1*D_data[t,2] - lambda_r[t]*W_r[t,2] for t=1:numb_timesteps)
                        - tar*Wtp_r[2] - C_PV*PVi_r[2] - C_S*Si_r[2]

@show system_costs
@show producer_surplus
@show active_surplus
@show passive_surplus
@show consumer_surplus
@show Dist_cost_active
@show Dist_cost_passive
@show Costs_rec_year

### Plots
################################################################################

df1 = DataFrame(t=1:numb_timesteps)
df1[:lambda] = lambda_r
df1[:CH] = CH_r[:,1]
df1[:DC] = DC_r[:,1]
df1[:E] = E_r[:,1]
df1[:W] = W_r[:,1]
df1[:wcons] = w_cons_tot
df1[:dcons] = d_cons_tot
df1[:wpros] = w_pros
df1[:dpros] = d_pros
push!(df1, @data([df1[end,:t]+1, df1[end,:lambda], df1[end,:CH], df1[end,:DC],
                df1[end,:E], df1[end,:W], df1[end,:wcons], df1[end,:dcons],
                df1[end,:wpros], df1[end,:dpros]]))
df1[:t] -= 1


# df2 = DataFrame(t=1:numb_timesteps)
# df2[:lambda] = lambda_r
# df2[:CH] = CH_r[:,2]
# df2[:DC] = DC_r[:,2]
# df2[:E] = E_r[:,2]
# df2[:W] = W_r[:,2]


df3 = DataFrame(t=1:numb_timesteps)
# df3[:P3] = P_r[:,3]
df3[:P2] = P_r[:,2]
df3[:P1] = P_r[:,1]
df3[:t] -= 0.5


dfs = stack(df3)





plot_lambda = plot(df1, x=:t, y=:lambda, Geom.step,
                Guide.xlabel("Time [h]"),Guide.ylabel("Price [€/kWh]"))

plot_ch = plot(df1, x=:t, y=:CH, Geom.step,
                Guide.xlabel("Time [h]"),Guide.ylabel("Charge [kW]"))

plot_dc = plot(df1, x=:t, y=:DC, Geom.step,
                Guide.xlabel("Time [h]"),Guide.ylabel("Discharge [kW]"))

filename = joinpath(@__DIR__, "Storage_test.pdf")
draw(PDF(filename, 27cm, 21cm), vstack(plot_lambda, plot_ch, plot_dc))

plot_cons = plot(
        layer(df1, x=:t, y=:dcons, Geom.step, Theme(default_color=colorant"red")),
        layer(df1, x=:t, y=:wcons, Geom.step, Theme(default_color=colorant"white")),
                Guide.xlabel("Time [h]"),Guide.ylabel("Offtake [kWh]"))

plot_pros = plot(
        layer(df1, x=:t, y=:dpros, Geom.step, Theme(default_color=colorant"red")),
        layer(df1, x=:t, y=:wpros, Geom.step, Theme(default_color=colorant"white")),
                Guide.xlabel("Time [h]"),Guide.ylabel("Offtake [kWh]"))

# plot_p = plot(dfs, x=:t, y=:value, color=:variable, Geom.bar(position=:stack))
