using DifferentialEquations; using CSV; using DataFrames;using StatsPlots; using DiffEqFlux;
using Statistics; using Distributions;
using Turing; using JLD2

healthy_p1= load("h_p1.jld2")["data"];
healthy_p2= load("h_p2.jld2")["data"];
healthy_p3= load("h_p3.jld2")["data"];
healthy_p4= load("h_p4.jld2")["data"];

ifg_p1= load("ifg_ml_p1.jld2")["data"];
ifg_p2= load("ifg_ml_p2.jld2")["data"];
ifg_p3= load("ifg_ml_p3.jld2")["data"];
ifg_p4= load("ifg_ml_p4.jld2")["data"];

igt_p1=load("igt_ml_p1.jld2")["data"];
igt_p2=load("igt_ml_p2.jld2")["data"];
igt_p3=load("igt_ml_p3.jld2")["data"];
igt_p4=load("igt_ml_p4.jld2")["data"];

ifg_igt_p1= load("ifg_igt_ml_p1.jld2")["data"];
ifg_igt_p2= load("ifg_igt_ml_p2.jld2")["data"];
ifg_igt_p3= load("ifg_igt_ml_p3.jld2")["data"];
ifg_igt_p4= load("ifg_igt_ml_p4.jld2")["data"];

t2dm_p1= load("t2dm_ml_p1.jld2")["data"];
t2dm_p2= load("t2dm_ml_p2.jld2")["data"];
t2dm_p3= load("t2dm_ml_p3.jld2")["data"];
t2dm_p4= load("t2dm_ml_p4.jld2")["data"];
#= healthy_pvals=  load("healthy_pvals.jld")["data"]
ifg_pvals=      load("ifg_pvals.jld")["data"]
ifg_igt_pvals=  load("ifg_igt_pvals.jld")["data"]
igt_pvals=      load("igt_pvals.jld")["data"]
t2dm_pvals=     load("t2dm_pvals.jld")["data"] =#

#Alpha and Beta values for Healthy individuals
healthy_alpha_p1=   mean((healthy_p1))^2/var((healthy_p1))
healthy_beta_p1=    mean((healthy_p1))/var((healthy_p1))

healthy_alpha_p2=   mean((healthy_p2))^2/var((healthy_p2))
healthy_beta_p2=    mean((healthy_p2))/var((healthy_p2))

healthy_alpha_p3=   mean((healthy_p3))^2/var((healthy_p3))
healthy_beta_p3=    mean((healthy_p3))/var((healthy_p3))

healthy_alpha_p4=   mean((healthy_p4))^2/var((healthy_p4))
healthy_beta_p4=    mean((healthy_p4))/var((healthy_p4))


#Alpha and Beta values for IFG individuals
ifg_alpha_p1=       mean((ifg_p1))^2/var((ifg_p1))
ifg_beta_p1=        mean((ifg_p1))/var((ifg_p1))

ifg_alpha_p2=       mean((ifg_p2))^2/var((ifg_p2))
ifg_beta_p2=        mean((ifg_p2))/var((ifg_p2))

ifg_alpha_p3=       mean((ifg_p3))^2/var((ifg_p3))
ifg_beta_p3=        mean((ifg_p3))/var((ifg_p3))

ifg_alpha_p4=       mean((ifg_p4))^2/var((ifg_p4))
ifg_beta_p4=        mean((ifg_p4))/var((ifg_p4))


#Alpha and Beta values for IGT individuals
mean(igt_p1)^2/var(igt_p1)
igt_alpha_p1=       mean((igt_p1))^2/var((igt_p1))
igt_beta_p1=        mean(igt_p1)/var(igt_p1)

igt_alpha_p2=       mean(igt_p2)^2/var(igt_p2)
igt_beta_p2=        mean(igt_p2)/var(igt_p2)

igt_alpha_p3=       mean(igt_p3)^2/var(igt_p3)
igt_beta_p3=        mean(igt_p3)/var(igt_p3)

igt_alpha_p4=       mean(igt_p4)^2/var(igt_p4)
igt_beta_p4=        mean(igt_p4)/var(igt_p4)

#Alpha and Beta values for IFG_IGT individuals
ifg_igt_alpha_p1=   mean(ifg_igt_p1)^2/var(ifg_igt_p1)
ifg_igt_beta_p1=    mean(ifg_igt_p1)/var(ifg_igt_p1)

ifg_igt_alpha_p2=   mean(ifg_igt_p2)^2/var(ifg_igt_p2)
ifg_igt_beta_p2=    mean(ifg_igt_p2)/var(ifg_igt_p2)

ifg_igt_alpha_p3=   mean(ifg_igt_p3)^2/var(ifg_igt_p3)
ifg_igt_beta_p3=    mean(ifg_igt_p3)/var(ifg_igt_p3)

ifg_igt_alpha_p4=   mean(ifg_igt_p4)^2/var(ifg_igt_p4)
ifg_igt_beta_p4=    mean(ifg_igt_p4)/var(ifg_igt_p4)

#Alpha and Beta values for T2DM individuals
t2dm_alpha_p1=   mean((t2dm_p1))^2/var((t2dm_p1))
t2dm_beta_p1=    mean((t2dm_p1))/var((t2dm_p1))

t2dm_alpha_p2=   mean((t2dm_p2))^2/var((t2dm_p2))
t2dm_beta_p2=    mean((t2dm_p2))/var((t2dm_p2))
mean((t2dm_p2))
t2dm_alpha_p3=   mean((t2dm_p3))^2/var((t2dm_p3))
t2dm_beta_p3=    mean((t2dm_p3))/var((t2dm_p3))

t2dm_alpha_p4=   mean((t2dm_p4))^2/var((t2dm_p4))
t2dm_beta_p4=    mean((t2dm_p4))/var((t2dm_p4))


 plot(Gamma(healthy_alpha_p1,1/healthy_beta_p1),xlabel="Sample",ylabel="Density",title="PDF of Healthy k1")
plot(Gamma(healthy_alpha_p2,1/healthy_beta_p2),xlabel="Sample",ylabel="Density",title="PDF of Healthy k2")
plot(Gamma(healthy_alpha_p3,1/healthy_beta_p3),xlabel="Sample",ylabel="Density",title="PDF of Healthy k3")
plot(Gamma(healthy_alpha_p4,1/healthy_beta_p4),xlabel="Sample",ylabel="Density",title="PDF of Healthy k4")

plot(Gamma(ifg_alpha_p1,1/ifg_beta_p1),title="PDF of IFG k1",xlabel="Sample",ylabel="Density")
plot(Gamma(ifg_alpha_p2,1/ifg_beta_p2),title="PDF of IFG k2",xlabel="Sample",ylabel="Density")
plot(Gamma(ifg_alpha_p3,1/ifg_beta_p3),title="PDF of IFG k3",xlabel="Sample",ylabel="Density")
plot(Gamma(ifg_alpha_p4,1/ifg_beta_p4),title="PDF of IFG k4",xlabel="Sample",ylabel="Density")

plot(Gamma(igt_alpha_p1,1/igt_beta_p1),title="PDF of IGT k1",xlabel="Sample",ylabel="Density")
plot(Gamma(igt_alpha_p2,1/igt_beta_p2),title="PDF of IGT k2",xlabel="Sample",ylabel="Density")
plot(Gamma(igt_alpha_p3,1/igt_beta_p3),title="PDF of IGT k3",xlabel="Sample",ylabel="Density")
plot(Gamma(igt_alpha_p4,1/igt_beta_p4),title="PDF of IGT k4",xlabel="Sample",ylabel="Density")

plot(Gamma(ifg_igt_alpha_p1,1/ifg_igt_beta_p1),title="PDF of IFG_IGT k1",xlabel="Sample",ylabel="Density")
plot(Gamma(ifg_igt_alpha_p2,1/ifg_igt_beta_p2),title="PDF of IFG_IGT k2",xlabel="Sample",ylabel="Density")
plot(Gamma(ifg_igt_alpha_p3,1/ifg_igt_beta_p3),title="PDF of IFG_IGT k3",xlabel="Sample",ylabel="Density")
plot(Gamma(ifg_igt_alpha_p4,1/ifg_igt_beta_p4),title="PDF of IFG_IGT k4",xlabel="Sample",ylabel="Density")

plot(Gamma(t2dm_alpha_p1,1/t2dm_beta_p1),title="PDF of T2DM k1",xlabel="Sample",ylabel="Density")
plot(Gamma(t2dm_alpha_p2,1/t2dm_beta_p2),title="PDF of T2DM k2",xlabel="Sample",ylabel="Density")
plot(Gamma(t2dm_alpha_p3,1/t2dm_beta_p3),title="PDF of T2DM k3",xlabel="Sample",ylabel="Density")
plot(Gamma(t2dm_alpha_p4,1/t2dm_beta_p4),title="PDF of T2DM k4",xlabel="Sample",ylabel="Density")
 =#

function diffeq_with_cb(du,u,p,t)
    c = [
        0.00551, #f index 1
        17/70, # vg index 2
        0.043, # gbliv index 3
        1, # beta index 4
        31, # taui index 5
        3, # taud index 6
        13/70, # vi index 7
        9, # Gthpl index 8
        30, # t_integralwindow index 9
        0.1, # c1 index 10
        0.043*(p[11]+p[12])/p[12]-p[5]*1*p[13], # c2 index 11
        p[7]*p[12]/(1*31*p[13])*30 # c3 index 12
    ] 
    Mg,Gpl,Ipl,Int=u

    if t > input[3]
        mgmeal = p[10].*p[1].^p[10].*(t-input[3]).^(p[10]-1) .* exp(-(p[1].*(t-input[3])).^p[10]) .* input[1]
    else
        mgmeal=0;
    end
    mgpl = p[2] .* Mg
    du[1] = mgmeal - mgpl
    
    # glucose in plasma
    gliv    = c[3] - p[3] .* (Gpl-p[12])-p[4] .* c[4] .* (Ipl-p[13])

    ggut    = p[2] .* (c[1]/(c[2]*input[2])) .* Mg;
    gnonit  = c[11] .* (Gpl./(p[11]+Gpl));
    git     = p[5] .* c[4] .* Ipl .* (Gpl./(p[11]+Gpl))
    if Gpl > c[8]
        gren  = c[10] ./ (c[2]*input[2]) .* (Gpl - c[8])
    else
        gren = 0
    end
    du[2] = gliv + ggut - gnonit - git - gren

    
    #Insulin in plasma
    #u[3] represents Ipl
    #u[4] represents Gint
   
    Gpl_lowerbound = Gpl_saved[1]
    du[4] = (Gpl-p[12]) - (Gpl_lowerbound-p[12]); #dGint/dt
    ipnc = (c[4].^-1) .* (p[6].*(Gpl-p[12]) +#=(p[7]/c[5]).*Int+=# (p[7]/c[5]) .*p[12]* c[9] + (p[8].*c[6]) .* du[2])

    iliv    = c[12].*Ipl;
    iif     = p[9].*(Ipl-p[13])
    
    du[3] = ipnc - iliv - iif
    
end
@model function fitdiffeq(data,prob1)

    phi ~ InverseGamma(2,3)
    #alpha ~ truncated(Normal((0.013371373352062), (0.01)),(0.0001),(0.1))
    #alpha ~ Beta(0.03,1)
    alpha ~ Gamma(ifg_alpha_p1,1/ifg_beta_p1) 
    #beta ~ Beta(0.02,1)
    beta ~ Gamma(ifg_alpha_p2,1/ifg_beta_p2)
    #beta ~ truncated(Normal((.0025428), (0.0003)),(0.0001), (0.05))
    #gamma ~ truncated(Normal((0.49758), (0.5)),(0.1),(10))
    gamma ~ Gamma(ifg_alpha_p3, 1/ifg_beta_p3)
    #delta ~ truncated(Normal((2.48325), (0.1)),(0.583886959936111*0.001), (10))
    delta ~ Gamma(ifg_alpha_p4,1/ifg_beta_p4)
    p = [
	alpha  #k1
    6.33e-1  #k2
    5.00e-5  #k3
    1.00e-3  #k4
    beta  #k5
    gamma  #k6
    2.20e-2  #k7
    delta     #k8
    1.08e-2  #k9
	#2.6     REMOVED
    1.35     #sigma /index 10
    0.63     #Km /index 11
    true_data[z,1]     #Gb /index 12
    true_data[z,6]     #Ib /index 13
    ]    
    
    predicted= solve(prob1, alg_hints=[:stiff];
                    p=p,   
                    saveat=[0,30,60,90,120],maxiters=10000,save_idxs=[2,3])
    predicted[2,:]=predicted[2,:]/10
    Threads.@threads for i in 1:length(predicted)
        data[:,i] ~ MvNormal(predicted[i],phi)
    end

    return predicted
end
#odedata= Array(Matlab_Data)

ifg_test= CSV.read("IFG_df_test.csv",DataFrame)
data= CSV.read("example_data.csv",DataFrame)
true_data= Array(ifg_test);
z=1;
true_glu= true_data[z,3:7]
true_ins= true_data[z,9:13];
actual_data=hcat(true_glu,true_ins)';

#sampled_sol= Array( true_sol)+ 0.8 * randn(size(Array(true_sol)))
#sample_solution=sampled_sol[2:3,:]
p = [
	1.35e-2  #k1
    6.33e-1  #k2
    5.00e-5  #k3
    1.00e-3  #k4
    3.80e-3  #k5
    5.82e-1  #k6
    2.20e-2  #k7
    4.71     #k8
    1.08e-2  #k9
	#2.6     REMOVED
    1.35     #sigma /index 10
    0.63     #Km /index 11
    true_data[z,3]      #Gb /index 12
    true_data[z,9]     #Ib /index 13
]    

u0=[0.0,true_glu[1], true_ins[1],0.0]
Gpl_saved=u0[2];
tspan=(0.0,120.0);
prob= ODEProblem(diffeq_with_cb,u0,tspan,p);

actual_data[2,:]=actual_data[2,:]/10
model= fitdiffeq(actual_data, prob);

    # INPUTS
input = [
        75000, # D
        70, # Mb
        0 # t_meal_start
    ]
chain = sample(model, NUTS(0.65), MCMCThreads() ,1000,1 ; progress=true)
p[1]=mean(chain["alpha"])
p[5]= mean(chain["beta"])
p[6]= mean(chain["gamma"])
p[8]= mean(chain["delta"])

prob= ODEProblem(diffeq_with_cb,u0,tspan,p)
solution= solve(prob, alg= Rosenbrock23(), saveat=1)

sim_glu= plot(solution, vars=2, xlabel="time",ylabel="glucose, mmol/L", labels="Combination sampling");
sim_glu= scatter!([0,30,60,90,120],true_data[z,3:7], labels="Measured data", legend= :outertopright);
#sim_glu=plot!(Vector(glucose_data[z,:]),labels="Matlab");
sim_ins= plot(solution, vars=3, xlabel="time",ylabel="insulin, mmol/L", labels="Combination sampling");
sim_ins= scatter!([0,30,60,90,120],true_data[z,9:13], labels="Measured data", legend= :outertopright);
#sim_ins= plot!(Vector(insulin_data[z,:]),labels="Matlab");
x=plot(plot(sim_glu,sim_ins),layout=(1,1),legendfontsize=8,legend=:outertopleft,size=(1600,900));
display(x) 