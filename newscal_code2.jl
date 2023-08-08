#import necessary libraries
using DifferentialEquations
using PyPlot
using DelimitedFiles
using Random

#define function which returns solution to set of differential equations
function photnumbevo(T_end, K, NCl, ν, σv, i, γ, g)

#define fixed parameters
κ=1
Δ=0
   
function G(x)
    λ=1
    k=(2*pi)/λ
    return g*cos(k*x)
end
    
#define atom trajectories
x0=[j/NCl for j=1:NCl]
Random.seed!(i)
v0=[randn()*σv for j=1:NCl]

function cavityall(du, u, p, t)  
          
xt=[x0[i]+t*v0[i] for i=1:NCl]
#abbreviations
aσparr=[G(x0[i]+t*v0[i])*u[i+3] for i=1:NCl]
aσp=sum(aσparr)
        
aσmarr=[G(x0[i]+t*v0[i])*u[NCl+i+3] for i=1:NCl]
aσm=sum(aσmarr)
        
σpσmarr=[G(x0[i]+t*v0[i])*u[2*NCl+i+3] for i=1:NCl]
σpσm=sum(σpσmarr)
        
σmarr=[G(x0[i]+t*v0[i])*u[3*NCl+i+3] for i=1:NCl]
σm=sum(σmarr)
pσm=[σm-σmarr[i] for i=1:NCl]
  
                
# write down manually derived equations in mixed order
# for mean field we need completely different equations      
du[1]=-(κ/2-im*Δ)*u[1]  - im*K*σm #for <a>
du[2]= -(κ-2*im*Δ)*u[2] + im*K*(aσp - conj(aσp)) #for <a^†a>
du[3]=-(κ-2*im*Δ)*u[3] - 2*im*K*aσm #for <aa>
        
#use loop for each cluster 
for i=4:(3+NCl) #eq for <aσ+>
du[i]=  -(κ/2+γ/2+ν/2-im*Δ)*u[i] + im*G(xt[i-3])*u[2] -im*G(xt[i-3])*u[i+2*NCl]+
        -2*im*G(xt[i-3])*(u[2]*u[i+2*NCl]+u[1]*conj(u[i+4*NCl])+
        +u[i+4*NCl]*conj(u[1])-2*u[1]*conj(u[1])*u[i+2*NCl]) +
        -im*K*conj(u[i+3*NCl])*pσm[i-3] +
        -im*G(xt[i-3])*(K-1)*conj(u[i+3*NCl])*u[i+3*NCl]-im*G(xt[i-3])*u[i+2*NCl]
end   

for i=(4+NCl):(3+2*NCl) #eq for <aσ->
du[i]=  -(κ/2+γ/2+ν/2-im*Δ)*u[i]- im*G(xt[i-3-NCl])*u[3]+
        +2*im*G(xt[i-3-NCl])*(u[3]*u[i+NCl]+2*u[1]*u[i+3*NCl]-2*u[1]*u[1]*u[i+NCl]) +
        - im*K*u[i+2*NCl]*pσm[i-3-NCl]-im*G(xt[i-3-NCl])*(K-1)*u[i+2*NCl]*u[i+2*NCl]
            
end
            
for i=(4+2*NCl):(3+3*NCl) #eq for <σee>
du[i]=  - (γ+ν)*u[i] + ν+im*G(xt[i-3-2*NCl])*(conj(u[i-2*NCl]) - u[i-2*NCl]) 
end

for i=(4+3*NCl):(3+4*NCl) #eq for <σ->
du[i]=  -(γ+ν)/2*u[i]-im*G(xt[i-3-3*NCl])*u[1]+ 2*im*G(xt[i-3-3*NCl])*u[i+NCl]
end    

for i=(4+4*NCl):(3+5*NCl)  #eq for <aσee>        
du[i]=  -(γ+ν+κ/2-im*Δ)*u[i] + ν*u[1]- im*G(xt[i-3-4*NCl])*(u[3]*conj(u[i-NCl])+
        +2*u[1]*u[i-4*NCl]-2*u[1]*u[1]*conj(u[i-NCl])) +
        +im*G(xt[i-3-4*NCl])*(u[2]*u[i-NCl]+conj(u[i-4*NCl])*u[1]+
        +u[i-3*NCl]*conj(u[1])-2*conj(u[1])*u[1]*u[i-NCl]) +
        -im*K*u[i-2*NCl]*pσm[i-3-4*NCl]-im*G(xt[i-3-4*NCl])*(K-1)*u[i-2*NCl]*u[i-NCl]
end                 
end        
    
# initial conditions        
u0=[0.01+0.0*im for i=1:(3+5*NCl)]
    
tspan=(0.0, T_end)
prob=ODEProblem(cavityall, u0, tspan)
sol=solve(prob, saveat=0.1);
    
#pick solutions for interesting field operators
return [sol[1, :], sol[2, :], sol[3, :], sol[4, :], sol[2*NCl+2, :], sol[3*NCl+2, :]]  
end

#define function calculating for different seeds and save the data in an array
function calcdata(T_end, K, NCl, ν, σv, g, samplesize)

    #save data in txt file with parameters in name
    open("inhom/qmixed_Tend$T_end.K$K.NCl$NCl.nu$ν.sv$σv.g$g.txt", "w") do io
        ini=1
        for j=(ini):(ini+samplesize-1)   
        writedlm(io, [abs.(photnumbevo(T_end, K, NCl, ν, σv, j, 10^(-8), g)[2])])
        end
    end
end

tableK=[2500, 2000, 1500, 1000, 500]
scanv1=collect(0.0:0.01:0.1)
scanv2=collect(0.12:0.02:0.7)
scanv=0.5*vcat(scanv1, scanv2);

using ProgressMeter
prog=Progress(length(tableK)*length(scanv))

#calculate for all parameters using a loop
samplesize=50;
for i in 1:length(tableK)
Threads.@threads for j in 1:length(scanv)
            calcdata(400, tableK[i], 400, 0.5, scanv[j],0.00136, samplesize)
            next!(prog)
        end
end

#reads out data from txt file
function readout_te(T_end, K, NCl, ν, σv, g)
    yprob=readdlm("qmixed_Tend$T_end.K$K.NCl$NCl.nu$ν.sv$σv.g$g.txt");
    te=sum(yprob, dims=1) 
    return transpose(te)/size(yprob, 1) 
    #averaging over all trajectories with different random distributions
end

function readout_av(T_end, K, NCl, ν, σv, g)
    yprob=readout_te(T_end, K, NCl, ν, σv, g)
    #average over each trajectory
    sum=0;
    start=2000; #at this time the quasi-steady-state is achieved
    size=length(yprob)-start+1;
    for i=start:length(yprob)
        sum+=yprob[i]
    end
    return sum/size
end 

#reads out average quasi-steady state values and plots them
function readout_graph(T_end, K, NCl, ν, g)
    graph=zeros(length(scanv))
    graph=readout_av.(T_end, K, NCl, ν, scanv, g)
 plot(scanv, graph) 
return graph 
end
