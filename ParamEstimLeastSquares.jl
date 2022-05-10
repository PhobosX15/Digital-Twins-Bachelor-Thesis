using DifferentialEquations
using Pkg

using  CSV
using DataFrames
using Plots

data= CSV.read("example_data.csv",DataFrame)
p = [
	1.35e-2  #k1
    
    3.80e-3  #k5
    5.82e-1  #k6
   
    4.71     #k8
]    


p_fixed = [
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
    data[1,1]     #Gb /index 12
    data[1,6]    #Ib /index 13
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
        0.043*(p_fixed[11]+p_fixed[12])/p_fixed[12]-p[2]*1*p_fixed[13], # c2 index 11
        p_fixed[7]*p_fixed[12]/(1*31*p_fixed[13])*30 # c3 index 12
    ] 
    
    # INPUTS
input = [
        75000, # D
        70, # Mb
        0 # t_meal_start
    ]

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
        0.043*(p_fixed[11]+p_fixed[12])/p_fixed[12]-p[2]#=p[5]=#*1*p_fixed[13], # c2 index 11
        p_fixed[7]*p_fixed[12]/(1*31*p_fixed[13])*30 # c3 index 12
    ] 
    Mg,Gpl,Ipl,Int=u

    if t > input[3]
        mgmeal = p_fixed[10].*p[1].^p_fixed[10].*(t-input[3]).^(p_fixed[10]-1) .* exp(-(p[1].*(t-input[3])).^p_fixed[10]) .* input[1]
    else
        mgmeal=0;
    end
    mgpl = p_fixed[2] .* Mg
    du[1] = mgmeal - mgpl
    
    # glucose in plasma
    gliv    = c[3] - p_fixed[3] .* (Gpl-p_fixed[12])-p_fixed[4] .* c[4] .* (Ipl-p_fixed[13])

    ggut    = p_fixed[2] .* (c[1]/(c[2]*input[2])) .* Mg;
    gnonit  = c[11] .* (Gpl./(p_fixed[11]+Gpl));
    git     = p[2]#=p5=# .* c[4] .* Ipl .* (Gpl./(p_fixed[11]+Gpl))
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

    #=if (t > c[9]) && (length(t_saved)>1) && (length(t_saved) == length(Gpl_saved))
        spl= Spline1D(t_saved,Gpl_saved)
        Gpl_lowerbound= spl(t_lowerbound)
    else
        Gpl_lowerbound = Gpl_saved[1]  #is called when t < t_integralwindow, or if there is no saved step yet (steps are only saved at pre-defined time points)
    end=#
    Gpl_lowerbound = Gpl_saved[1]
    du[4] = (Gpl-p_fixed[12]) - (Gpl_lowerbound-p_fixed[12]); #dGint/dt
    ipnc = (c[4].^-1) .* (p[3]#=p6=#.*(Gpl-p_fixed[12]) +#=(p[7]/c[5]).*Int+=# (p_fixed[7]/c[5]) .*p_fixed[12]* c[9] + (p[4]#=p8=#.*c[6]) .* du[2])

    iliv    = c[12].*Ipl;
    iif     = p_fixed[9].*(Ipl-p_fixed[13])
    
    du[3] = ipnc - iliv - iif
    
end
u0=[0.0,data[1,1],data[1,6],0.0]
tspan=(0.0,120.0)
Gpl_saved=[5.0]
Ipl_saved=[11.5]

t=[0,30,60,90,120]
prob= ODEProblem(diffeq_with_cb,u0,tspan,p)
#= sol = solve(prob,QNDF(),p=p,saveat =[0.0,30.0,60.0,90.0,120.0])
sol
    solution= sol[2:3,:]
    solution=solution[2,:]./10
    solution =#
#Using Flux
true_data= Array(data)
true_glu= true_data[1,1:5]
true_ins= true_data[1,6:10]
actual_data=hcat(true_glu,true_ins)'
normalized_data= actual_data
normalized_data[2,:]= normalized_data[2,:]/10
#= sampled_sol= Array( sol)+ 0.8 * randn(size(Array(sol)))
sampled_sol= sampled_sol[2:3,:]
sampled_sol[2,:]=sampled_sol[2,:]./10
sampled_sol

actual_data =#
using DiffEqFlux
abs2
function loss(p)
    sol = solve(prob,QNDF(),p=p,saveat =[0.0,30.0,60.0,90.0,120.0],save_idxs=[2,3])
    sol[2,:]= sol[2,:]/10
    
    loss= sum(abs2, sol.-normalized_data)
    return loss,sol
    
end

callback = function (p,l,pred)
    display(l)
    plt= plot(pred)
    display(plt)
    return false
end

result_ode= DiffEqFlux.sciml_train(loss, p,cb=callback,maxiters=100; 
                                    lower_bounds=[(0.013511373352062*0.001),(0.003806228136097*0.01),(0.583886959936111*0.001),(0.583886959936111*0.001)], 
                                    upper_bounds=[(0.013511373352062*10),(0.003806228136097*10),(0.583886959936111*10),(0.583886959936111*10)]
                                    )
p=result_ode.u

sampled_prob= ODEProblem(diffeq_with_cb,u0,tspan,p)
sampled_sol= solve(sampled_prob,QNDF(),saveat=1,save_idxs=[2,3])

sim_glu= plot(sampled_sol, vars=1, xlabel="time",ylabel="glucose, mmol/L", labels="Simulation with uniform sampled params");
sim_glu= scatter!([0,30,60,90,120],actual_data[1,:], labels="Measured data", legend= :outertopright);

sim_ins= plot(sampled_sol, vars=2, xlabel="time",ylabel="insulin, mmol/L", labels="Simulation with uniform sampled params");
sim_ins= scatter!([0,30,60,90,120],actual_data[2,:], labels="Measured data", legend= :outertopright);

x=plot(plot(sim_glu,sim_ins),layout=(1,1),legendfontsize=8,legend=:outertopleft,size=(1600,900));
display(x)
