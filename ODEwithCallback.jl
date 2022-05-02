using Pkg;using CSV; using DataFrames; using Dierckx;using DifferentialEquations;using  LSODA;using IfElse;
using DiffEqDevTools
using FileIO;using JLD2;
using Plots
data= CSV.read("example_data.csv",DataFrame)
est_p_values= CSV.read("p_values.csv", DataFrame)
glucose_data=CSV.read("glucose_data.csv", DataFrame)
insulin_data=CSV.read("insulin_data.csv", DataFrame)

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
    5.0      #Gb /index 12
    8.58     #Ib /index 13
]    

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
    
    # INPUTS
input = [
        75000, # D
        70, # Mb
        0 # t_meal_start
    ]

function diffeq_with_cb(du,u,p,t)
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
    t_lowerbound= t- c[9];

    if (t > c[9]) && (length(t_saved)>1) && (length(t_saved) == length(Gpl_saved))
        spl= Spline1D(t_saved,Gpl_saved)
        Gpl_lowerbound= spl(t_lowerbound)
    else
        Gpl_lowerbound = Gpl_saved[1]  #is called when t < t_integralwindow, or if there is no saved step yet (steps are only saved at pre-defined time points)
    end
    #Gpl_lowerbound = Gpl_saved[1]
    du[4] = (Gpl-p[12]) - (Gpl_lowerbound-p[12]); #dGint/dt
    ipnc = (c[4].^-1) .* (p[6].*(Gpl-p[12]) +#=(p[7]/c[5]).*Int+=# (p[7]/c[5]) .*p[12]* c[9] + (p[8].*c[6]) .* du[2])

    iliv    = c[12].*Ipl;
    iif     = p[9].*(Ipl-p[13])
    
    du[3] = ipnc - iliv - iif
    
end

function differentialequations_1(du,u,p,t)
    Mg,Gpl,Ipl,Int=u
    
   
    # INPUTS
    

   
    mgmeal = IfElse.ifelse(t> input[3], p[10].*p[1].^p[10].*(t-input[3]).^(p[10]-1) .* exp(-(p[1].*(t-input[3])).^p[10]) .* input[1],0)
    
    mgpl = p[2] .* Mg
    du[1] = mgmeal - mgpl
    
    # glucose in plasma
    gliv    = c[3] - p[3] .* (Gpl-p[12])-p[4] .* c[4] .* (Ipl-p[13])

    ggut    = p[2] .* (c[1]/(c[2]*input[2])) .* Mg;
    gnonit  = c[11] .* (Gpl./(p[11]+Gpl));
    git     = p[5] .* c[4] .* Ipl .* (Gpl./(p[11]+Gpl))
    gren= IfElse.ifelse(Gpl > c[8], c[10] ./ (c[2]*input[2]) .* (Gpl - c[8]), 0 )
    du[2] = gliv + ggut - gnonit - git - gren

    
    #Insulin in plasma
    #u[3] represents Ipl
    #u[4] represents Gint
    t_lowerbound= t- c[9];

    #=if (t > c[9]) && (length(t_saved)>1) && (length(t_saved) == length(Gpl_saved))
        spl=Spline1D(t_saved,Gpl_saved)
        Gpl_lowerbound= spl(t_lowerbound);
    else
        Gpl_lowerbound = Gpl_saved[1];  #is called when t < t_integralwindow, or if there is no saved step yet (steps are only saved at pre-defined time points)
    end =#
    Gpl_lowerbound = Gpl_saved[1]
    du[4] = (Gpl-p[12]) - (Gpl_lowerbound-p[12]); #dGint/dt
    ipnc = (c[4].^-1) .* (p[6].*(Gpl-p[12]) #=+(p[7]/c[5]).*Int=#+ (p[7]/c[5]) .*p[12]* c[9] + (p[8].*c[6]) .* du[2])

    iliv    = c[12].*Ipl;
    iif     = p[9].*(Ipl-p[13])
    
    du[3] = ipnc - iliv - iif

end
#= setups = [Dict(:alg=>Rodas4())
          Dict(:alg=>Rosenbrock23())
          Dict(:alg=>Rodas5())
          Dict(:alg=>QNDF())
          ]
labels = [
    "Julia: Rodas4"
    "Julia: Rosenbrock23"
    "Julia: Rodas5"
    "Julia: QNDF"
]
abstols = 1.0 ./ 10.0 .^ (6:13)
reltols = 1.0 ./ 10.0 .^ (3:10)
wp = WorkPrecisionSet(prob,abstols,reltols,setups;
                      names = labels,print_names = true,
                      appxsol=sol,dense=false,
                      save_everystep=false,numruns=100,maxiters=10000000,
                      timeseries_errors=false,verbose=false)
        
plot(wp,title="Stiff ODE Problem",legend=:outertopleft,
                      color=permutedims([repeat([:LightGreen],1)...,repeat([:DarkGreen],1)...,
                      :Red,repeat([:Orange],1)...,repeat([:Yellow],3)...,
                      repeat([:Blue],2)...,:Purple]),size = (800,350),
                      xticks = 10.0 .^ (-12:1:5),
                      yticks = 10.0 .^ (-6:0.5:5),
                      bottom_margin=5Plots.mm)
                    
length(sol)
aproxsol= zeros()
 Vector(glucose_data[3,:])
ndims(sol) =#

#= p[1]=0.009332432;
p[5]= 0.001325156;
p[6]= 0.360093938;
p[8]=3.36292971;

preset_times=[0:1:120;]
tspan=(0,120)
u0=[0.0,5.0,11.5,0.0]
Gpl_saved=[]
Ipl_saved=[]
t_saved=[]
affect!(integrator)= push!(Gpl_saved,integrator.u[2]) , push!(t_saved,integrator.t), push!(Ipl_saved, integrator.u[3])
cb= PresetTimeCallback(preset_times,affect!)
prob= ODEProblem(differentialequations_1,u0,tspan,p)
sol= solve(prob,callback=cb,saveat=preset_times)
x=2 =#

 #=
simglu=plot(sol, vars=2, xlabel="time",ylabel="glucose, mmol/L", labels="QNDF,cb")
simins=plot(sol, vars=3, xlabel="time",ylabel="glucose, mmol/L", labels="QNDF,cb")
x=plot(plot(simglu,simins),layout=(1,1),legendfontsize=8,legend=:outertopleft,size=(1400,1000));
display(x)

setups = [Dict(:alg=>Rodas4())
          Dict(:alg=>Rodas4(),:abstol=>1e-8,:reltol=>1e-4)
          Dict(:alg=>Rodas4(),:abstol=>1e-4,:reltol=>1e-4)
          Dict(:alg=>Rodas4(),:abstol=>1e-14,:reltol=>1e-14)
          Dict(:alg=>Rosenbrock23())
          Dict(:alg=>Rosenbrock23(),:abstol=>1e-8,:reltol=>1e-4)
          Dict(:alg=>Rosenbrock23(),:abstol=>1e-4,:reltol=>1e-4)
          Dict(:alg=>Rosenbrock23(),:abstol=>1e-14,:reltol=>1e-14)
          Dict(:alg=>Rodas5())
          Dict(:alg=>Rodas5(),:abstol=>1e-8,:reltol=>1e-4)
          Dict(:alg=>Rodas5(),:abstol=>1e-4,:reltol=>1e-4)
          Dict(:alg=>Rodas5(),:abstol=>1e-14,:reltol=>1e-14)
          Dict(:alg=>QNDF())
          Dict(:alg=>QNDF(),:abstol=>1e-8,:reltol=>1e-4)
          Dict(:alg=>QNDF(),:abstol=>1e-4,:reltol=>1e-4)
          Dict(:alg=>QNDF(),:abstol=>1e-14,:reltol=>1e-14)
          ]

labels = [
            "Rodas4"
            "Rodas4 a 1e-8,r 1e-4"
            "Rodas4 a 1e-4,r 1e-4"
            "Rodas4	a 1e-14,r 1e-14"
            "Rosenbrock23"
            "Rosenbrock23 a 1e-8,r 1e-4"
            "Rosenbrock23 a 1e-4,r 1e-4"
            "Rosenbrock23 a 1e-14,r 1e-14"
            "Rodas5"
            "Rodas5 a 1e-8,r 1e-4"
            "Rodas5 a 1e-4,r 1e-4"
            "Rodas5	a 1e-14,r 1e-14"
            "QNDF"
            "QNDF a 1e-8,r 1e-4"
            "QNDF a 1e-4,r 1e-4"
            "QNDF a 1e-14,r 1e-14"
            
        ]


labels=["Rodas","Rosenbrock23","Rodas5","QNDF"] =#
Gpl_saved=[]
Ipl_saved=[]
t_saved=[]
Errors_mse= Dict()
Errors_mae= Dict()
for i in 1:length(data[!,1])
    filename= "QNDF_nocb_scen"*string(i)*".png"
    preset_times=[0:1:120;]
    tspan=(0,120)
    u0=[0.0,data[i,1],data[i,6],0.0]
    p[1]= est_p_values[i,1];
    p[5]= est_p_values[i,2];
    p[6]= est_p_values[i,3];
    p[8]= est_p_values[i,4];
    Gpl_saved=[]
    Ipl_saved=[]
    t_saved=[]
    affect!(integrator)= push!(Gpl_saved,integrator.u[2]) , push!(t_saved,integrator.t), push!(Ipl_saved, integrator.u[3])
    cb= PresetTimeCallback(preset_times,affect!)
    prob= ODEProblem(diffeq_with_cb,u0,tspan,p)
    sol1= solve(prob,QNDF(),callback=cb,saveat=1)
    Gpl_saved=[]
    Ipl_saved=[]
    t_saved=[]
    sol2= solve(prob,QNDF(), abstol=1e-8, reltol=1e-4,callback=cb,saveat=1)
    Gpl_saved=[]
    Ipl_saved=[]
    t_saved=[]
    sol3= solve(prob,QNDF(),abstol=1e-4, reltol=1e-4,callback=cb,saveat=1)
    Gpl_saved=[]
    Ipl_saved=[]
    t_saved=[]
    sol4= solve(prob,QNDF(), abstol=1e-14, reltol=1e-14,callback=cb,saveat=1)

    msescore1=0
    msescore2=0
    msescore3=0
    msescore4=0

    msescore1=(((sol1(30)[2]- data[i,2]))^2+((sol1(60)[2]- data[i,3]))^2+((sol1(90)[2]- data[i,4]))^2+((sol1(120)[2]- data[i,5]))^2)/4
    Errors_mse[filename*"glu_mse_autotol"]= msescore1
    msescore2=(((sol2(30)[2]- (data[i,2]))^2)+((sol2(60)[2]- (data[i,3]))^2)+((sol2(90)[2]- (data[i,4]))^2)+((sol2(120)[2]- (data[i,5]))^2))/4
    Errors_mse[filename*"glu_mse_r-8,a-4"]= msescore2
    msescore3=(((sol3(30)[2]- (data[i,2]))^2)+((sol3(60)[2]- (data[i,3]))^2)+((sol3(90)[2]- (data[i,4]))^2)+((sol3(120)[2]- (data[i,5]))^2))/4
    Errors_mse[filename*"glu_mse_r-4,a-4"]= msescore3
    msescore4=(((sol4(30)[2]- (data[i,2]))^2)+((sol4(60)[2]- (data[i,3]))^2)+((sol4(90)[2]- (data[i,4]))^2)+((sol4(120)[2]- (data[i,5]))^2))/4
    Errors_mse[filename*"glu_mse_r-14,a-14"]=msescore4

    maescore1=0
    maescore2=0
    maescore3=0
    maescore4=0

    maescore1= (abs((sol1(30)[2]- (data[i,2])))+ abs((sol1(60)[2]- (data[i,3])))+abs((sol1(90)[2]- (data[i,4])))+abs((sol1(120)[2]- (data[i,5]))))/4
    Errors_mae[filename*"glu_mae_autotol"]= maescore1
    maescore2= (abs((sol2(30)[2]- (data[i,2])))+ abs((sol2(60)[2]- (data[i,3])))+abs((sol2(90)[2]- (data[i,4])))+abs((sol2(120)[2]- (data[i,5]))))/4
    Errors_mae[filename*"glu_mae_r-8,a-4"]= maescore2
    maescore3= (abs((sol3(30)[2]- (data[i,2])))+ abs((sol3(60)[2]- (data[i,3])))+abs((sol3(90)[2]- (data[i,4])))+abs((sol3(120)[2]- (data[i,5]))))/4
    Errors_mae[filename*"glu_mae_r-4,a-4"]= maescore3
    maescore4= (abs((sol4(30)[2]- (data[i,2])))+ abs((sol4(60)[2]- (data[i,3])))+abs((sol4(90)[2]- (data[i,4])))+abs((sol4(120)[2]- (data[i,5]))))/4
    Errors_mae[filename*"glu_mae_r-14,a-14"]= maescore4

    sim_glu= plot(sol1, vars=2, xlabel="time",ylabel="glucose, mmol/L", labels="QNDF,nocb");
    sim_glu= plot!(sol2, vars=2, labels="QNDF, a-8,r-4,_nocb");
    sim_glu= plot!(sol3, vars=2, labels="QNDF, a-4,r-4,_nocb");
    sim_glu= plot!(sol4, vars=2, labels="QNDF, a-14,r-14,_nocb");
    sim_glu= plot!(Vector(glucose_data[i,:]), labels="MATLAB");
    sim_glu= scatter!([0,30,60,90,120],Vector(data[i,1:5]), labels="Measured data", legend= :outertopright);

    sim_ins= plot(sol1, vars=3, xlabel="time",ylabel="insulin, mmol/L", labels="QNDF,nocb");
    sim_ins= plot!(sol2, vars=3, labels="QNDF, a-8,r-4,_nocb");
    sim_ins= plot!(sol3, vars=3, labels="QNDF, a-4,r-4,_nocb");
    sim_ins= plot!(sol4, vars=3, labels="QNDF, a-14,r-14,_nocb");
    sim_ins= plot!(Vector(insulin_data[i,:]), labels="MATLAB");
    sim_ins= scatter!([0,30,60,90,120],Vector(data[i,6:10]), labels="Measured data", legend= :outertopright);

    x=plot(plot(sim_glu,sim_ins),layout=(1,1),legendfontsize=8,legend=:outertopleft,size=(1600,800));
    display(x)
    png(filename)
end

sort(Errors_mae; byvalue=true)
sort(Errors_mse; byvalue=true)
save("errors_mse_QNDF_nocb.jld2",Errors_mse)
save("errors_mae_QNDF_nocb.jld2", Errors_mae)

QNDF_nocb_mse= load("errors_mse_QNDF_nocb.jld2")
QNDF_nocb_mae=load("errors_mae_QNDF_nocb.jld2")
QNDF_cb_mse=load("errors_mse_QNDF_cb.jld2")
QNDF_cb_mae=load("errors_mae_QNDF_cb.jld2")
Rodas4_nocb_mse= load("errors_mse_Rodas4_nocb.jld2")
Rodas4_nocb_mae= load("errors_mae_Rodas4_nocb.jld2")
Rodas5_nocb_mse= load("errors_mse_Rodas5_nocb.jld2")
Rodas5_nocb_mae= load("errors_mae_Rodas5_nocb.jld2")

merged_mse_errors= merge(QNDF_nocb_mse,Rodas4_nocb_mse,Rodas5_nocb_mse,QNDF_cb_mse)
merged_mae_errors= merge(QNDF_nocb_mae,QNDF_cb_mae,Rodas4_nocb_mae,Rodas5_nocb_mae)

merged_mse_errors=sort(merged_mse_errors, byvalue= true)
merged_mae_errors=sort(merged_mae_errors, byvalue= true)

CSV.write("merged_mae_errors.xlsx",merged_mae_errors)
CSV.write("merged_mse_errors.csv",merged_mse_errors)
merged_mse_errors.keys
length(glucose_data[!,1])
glucose_data[1,121]
matlab_dict=Dict()
for i in 1:length(glucose_data[!,1])
    filename= "Scen"* string(i)
    msescore=0
    msescore=(((glucose_data[i,30]- data[i,2]))^2+((glucose_data[i,60]- data[i,3]))^2+((glucose_data[i,90]- data[i,4]))^2+((glucose_data[i,120]- data[i,5]))^2)/4
    maescore=0
    maescore= (abs((glucose_data[i,30]- (data[i,2])))+ abs((glucose_data[i,60]- (data[i,3])))+abs((glucose_data[i,90]- (data[i,4])))+abs((glucose_data[i,120]- (data[i,5]))))/4
    mse_ins=0
    mse_ins=(((insulin_data[i,30]- data[i,7]))^2+((insulin_data[i,60]- data[i,8]))^2+((insulin_data[i,90]- data[i,9]))^2+((insulin_data[i,120]- data[i,10]))^2)/4
    mae_ins=0
    mae_ins= (abs((insulin_data[i,30]- (data[i,7])))+ abs((insulin_data[i,60]- (data[i,8])))+abs((insulin_data[i,90]- (data[i,9])))+abs((insulin_data[i,120]- (data[i,10]))))/4
    matlab_dict[filename*"mse_glu"]= msescore
    matlab_dict[filename*"mae_glu"]= maescore
    matlab_dict[filename*"mse_ins"]= mse_ins
    matlab_dict[filename*"mae_ins"]= mae_ins

end
data[5,2]
matlab_dict
CSV.write("Matlab_errors.csv", matlab_dict)

for i in 1:length(data[!,1])
    filename= "Best_models "*string(i)*".png"
    preset_times=[0:1:120;]
    tspan=(0,120)
    u0=[0.0,data[i,1],data[i,6],0.0]
    p[1]= est_p_values[i,1];
    p[5]= est_p_values[i,2];
    p[6]= est_p_values[i,3];
    p[8]= est_p_values[i,4];
    Gpl_saved=[]
    Ipl_saved=[]
    t_saved=[]
    affect!(integrator)= push!(Gpl_saved,integrator.u[2]) , push!(t_saved,integrator.t), push!(Ipl_saved, integrator.u[3])
    cb= PresetTimeCallback(preset_times,affect!)
    prob= ODEProblem(diffeq_with_cb,u0,tspan,p)
    sol1= solve(prob,QNDF(),abstol=1e-8, reltol=1e-4,callback=cb,saveat=1)
    Gpl_saved=[]
    Ipl_saved=[]
    t_saved=[]
    prob2= ODEProblem(differentialequations_1,u0,tspan,p)
    sol2= solve(prob2,QNDF(),callback=cb,saveat=1)
    Gpl_saved=[]
    Ipl_saved=[]
    t_saved=[]
    sol3= solve(prob2,Rodas5(),abstol=1e-4, reltol=1e-8,callback=cb,saveat=1)
    
    sim_glu= plot(sol1, vars=2, xlabel="time",ylabel="glucose, mmol/L", labels="QNDF,cb,r-14,a-14");
    sim_glu= plot!(sol2, vars=2, labels="QNDF, autotol,_nocb");
    sim_glu= plot!(sol3, vars=2, labels="Rodas5, a-4,r-8,_nocb");
    sim_glu= scatter!([0,30,60,90,120],Vector(data[i,1:5]), labels="Measured data", legend= :outertopright);

    sim_ins= plot(sol1, vars=3, xlabel="time",ylabel="insulin, mmol/L", labels="QNDF,cb,r-14,a-14");
    sim_ins= plot!(sol2, vars=3, labels="QNDF, autotol,_nocb");
    sim_ins= plot!(sol3, vars=3, labels="Rodas5, a-4,r-4,_nocb");
    sim_ins= scatter!([0,30,60,90,120],Vector(data[i,6:10]), labels="Measured data", legend= :outertopright);
    x=plot(plot(sim_glu,sim_ins),layout=(1,1),legendfontsize=8,legend=:outertopleft,size=(1600,900));
    display(x)
    png(filename)
end