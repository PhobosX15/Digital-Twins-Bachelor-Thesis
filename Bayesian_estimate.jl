#Load packages
#using Pkg  #(uncomment this if you need to add packages for the first time)
using DifferentialEquations; using CSV; using DataFrames; using Statistics; 
using DiffEqFlux; using Optim;
using Distributions; using Turing;
using StatsPlots, Measures;

# Load sample scenarios
example_data= Array(CSV.read("example_data.csv", DataFrame));

#  parameter values
global p=[
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
    missing     #Gb /index 12 (missing parametes depend on scenarios, will be updated below)
    missing    #Ib /index 13
];

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
    missing, # c2 index 11
    missing # c3 index 12
    ]  ;

# Input Meal
input = [
        75000, # D
        70, # Mb
        0 # t_meal_start
];

# Pre-existing data to estimate priors
data_p=CSV.read("ps.csv",DataFrame); # Contains all classes of individuals

general_alpha_p1= mean((data_p.p1))^2/var((data_p.p1))
general_beta_p1=  mean((data_p.p1))/var((data_p.p1))

general_alpha_p2= mean((data_p.p2))^2/var((data_p.p2))
general_beta_p2=  mean((data_p.p2))/var((data_p.p2))
julia_p3_values=[0.666,0.93315,0.3771,0.40025,0.2935,0.3155]
p3_var= var(julia_p3_values)
mean(julia_p3_values)
julia_p4_values= [4.723,2.4060,3.3602,2.0495,2.2768,0.0840]
p4_var= var(julia_p4_values)
general_alpha_p3= mean(julia_p3_values)^2/p3_var
general_beta_p3=  mean((julia_p3_values))/p3_var

general_alpha_p4= mean((julia_p4_values))^2/p4_var
general_beta_p4=  mean((julia_p4_values))/p4_var

# Adapted E-DES model
# input 1) Initial Parameter values, 2) Initial conditions and 3) time-span
function diffeq_bayes(du,u,p,t)
    
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
    ipnc = (c[4].^-1) .* (p[6].*(Gpl-p[12]) + (p[7]/c[5]) .*p[12]* c[9] + (p[8].*c[6]) .* du[2])

    iliv    = c[12].*Ipl;
    iif     = p[9].*(Ipl-p[13])
    
    du[3] = ipnc - iliv - iif
end


#set tspan
tspan=(0,120)
# Estimates value of all scenarios
for a in 1:7
    
    true_glu= example_data[a,1:5]
    true_ins= example_data[a,6:10]
    global p[12]= true_glu[1]
    global p[13]= true_ins[1]
    
    c[11]=0.043*(p[11]+p[12])/p[12]-p[5]*1*p[13]
    c[12]=p[7]*p[12]/(1*31*p[13])*30
    global Gpl_saved= p[12]
    u0=[0.0,true_glu[1], true_ins[1],0.0]
    prob=ODEProblem(diffeq_bayes,u0,tspan, p)
    @model function fitdiffeq(data,prob1)

        phi ~ InverseGamma(2,3) # Random noise
       
        alpha ~ Gamma(general_alpha_p1,1/general_beta_p1) 
        
        beta ~ Gamma(general_alpha_p2,1/general_beta_p2)
    
        gamma ~ Gamma(general_alpha_p3, 1/general_beta_p3)
       
        delta ~ Gamma(general_alpha_p4,1/general_beta_p4)
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
        true_glu[1]     #Gb /index 12
        true_ins[1]    #Ib /index 13
        ]   
        
        predicted= solve(prob1, alg=QNDF();
                        p=p,   
                        saveat=[0,30,60,90,120],maxiters=10000,save_idxs=[2,3])
        predicted[2,:]= predicted[2,:]/10
        for i in 1:length(predicted)
            data[:,i] ~ MvNormal(predicted[i],phi)
        end
    
        return predicted
    end
    
    
    actual_data=hcat(true_glu,true_ins)' # combining glucose and insulin as true_solution
    actual_data[2,:] = actual_data[2,:]/10 # normalizing insulin by a factor of 0.1
    model= fitdiffeq(actual_data, prob);

    chain = sample(model, NUTS(0.65),MCMCSerial() ,1000,1 ; progress=true) # Sampling
    
    # updating p values to Mean of chain
   
    
    global p=[
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
        true_glu[1]     #Gb /index 12 (missing parametes depend on scenarios, will be updated below)
        true_ins[1]    #Ib /index 13
    ];
    global p[1]=mean(chain["alpha"]) 
    global p[5]=mean(chain["beta"])
    global p[6]=mean(chain["gamma"])
    global p[8]=mean(chain["delta"])
    # remake problem
    prob1=ODEProblem(diffeq_bayes,u0,tspan, p)
    sol1= solve(prob1,alg=QNDF(),saveat =1,save_idxs=[2,3])
    
    #Plotting Glucose
    glu=  plot(sol1,vars=1,xlabel="time",ylabel="glucose, mmol/L",linecolor="green", lw=3,label= "Julia Bayesian Glucose")
    glu= scatter!([0,30,60,90,120],example_data[a,1:5], markersize=5 ,label="True Val");
    
    
    #Plotting Insulin
    ins= plot(sol1,vars=2,xlabel="time",ylabel="glucose, mmol/L",linecolor="blue", lw=3,label="Julia Bayesian Insulin")
    ins= scatter!([0,30,60,90,120],xlabel="time",ylabel="insulin, mU/L",example_data[a,6:10], markersize=5,label="True Val");
    x=plot(plot(glu,ins),layout=(1,2),legend=:outertop, title="Scenario "*string(a))
    display(x)

    # Resetting p values to avoid bad initial guess
    global p=[
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
    missing     #Gb /index 12 (missing parametes depend on scenarios, will be updated below)
    missing    #Ib /index 13
];
end