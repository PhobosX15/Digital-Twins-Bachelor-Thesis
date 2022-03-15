using Pkg
using PyCall
pygui(:tk)
using CSV; using DataFrames; using Dierckx;using DifferentialEquations;using PyPlot;

data = CSV.read("example_data.csv", DataFrame)
est_p_values= CSV.read("p_values.csv", DataFrame)
sim_MG=0
sim_Gint=0;
data_meal= 75e3; #in milligrams
meal_start_time=0;
bodymass= 70; # in kg
timespan = range(0,120,121);
time_intervals = [0,30,60,90,120]
p_values= zeros(length(data[!,1]),4)
residuals= zeros(length(data[!,1]),10)
sim_glu= zeros(length(data[!,1]),121);
sim_ins= zeros(length(data[!,1]),121);






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

function differentialequations_1(du,u,p,t)
    Mg,Gpl,Ipl,Int=u

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
        spl=Spline1D(t_saved,Gpl_saved)
        Gpl_lowerbound= spl(t_lowerbound);
    else
        Gpl_lowerbound = Gpl_saved[1];  #is called when t < t_integralwindow, or if there is no saved step yet (steps are only saved at pre-defined time points)
    end
    
    du[4] = (Gpl-p[12]) - (Gpl_lowerbound-p[12]); #dGint/dt
    ipnc = (c[4].^-1) .* (p[6].*(Gpl-p[12]) +(p[7]/c[5]).*Int+ (p[7]/c[5]) .*p[12]* c[9] + (p[8].*c[6]) .* du[2])

    iliv    = c[12].*Ipl;
    iif     = p[9].*(Ipl-p[13])
    
    du[3] = ipnc - iliv - iif

end
        
function compute(u0, Gpl_saved,Insulin_saved)
    
    t=0
    global t_saved=[0]
    

    while   t < 120
        tspan = (t, t+1)
        
        prob = ODEProblem(differentialequations_1,  u0,tspan,p)
        sol=solve(prob, DP5())
         t=  t+1

        push!(t_saved, t)
        push!(Gpl_saved, last(sol[2,:]))
        push!(Insulin_saved,last(sol[3,:]))
        u0=last(sol)
        
        
    end
    
    return Gpl_saved, Insulin_saved, u0
end 
    
#Gpl_saved, Insulin_saved, u=compute()

using Plots
pyplot()
Plots.PyPlotBackend()

function plotting(Gpl_saved, Insulin_saved)
    pygui(true)
    for i in 1: length(Gpl_saved)
        
        g =plot(Gpl_saved[i], vars=2, xlabel="time", ylabel="glucose, mmol/L", xlims=(0, 120),labels="Estimated")
        g=scatter!([0,30,60,90,120],Vector(data[i,1:5]),labels="Actual")
	    h = plot(Insulin_saved[i], vars=3, xlabel="time", ylabel="insulin, mU/L", xlims=(0, 120), labels="Estimated");
        h= scatter!([0,30,60,90,120], Vector(data[i,6:10]),labels="Actual");
        p=plot(plot(g,h), layout=(1,1),legend= :outertopleft )
        display(p)
        
    end

end



function scenarios()
    Gpl_all = []
    Insulin_all=[]
    for i in 1:length(data[!,1])
        u0=[0.0,data[i,1],data[i,6],0.0]
        global Gpl_saved= [data[i,1]]
        global Insulin_saved=[data[i,6]]
        p[1]= est_p_values[i,1];
        p[5]= est_p_values[i,2];
        p[6]= est_p_values[i,3];
        p[8]= est_p_values[i,4];
        Gpl_saved, Insulin_saved, u= compute(u0, Gpl_saved,Insulin_saved)
        push!(Gpl_all, Gpl_saved)
        push!(Insulin_all,Insulin_saved)
        
    end

    return Gpl_all, Insulin_all
end




Gpl_all, Insulin_all=scenarios()
length(data[!,1])
plotting(Gpl_all,Insulin_all)


