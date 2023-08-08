using DifferentialEquations
using DelimitedFiles
#using Statistics
#using PyPlot
using Random
cd(@__DIR__)

#Cavity with moving atoms

function photnumbevo(N, K, σv, i)

κ=1
γ=0
Δ=0
g=2*1.36*0.001
ν=0.5;  
    
function G(x)
    λ=1
    k=(2*pi)/λ
    return g*cos(k*x)
end
    
   
x0=[j/N for j=0:(N-1)]
Random.seed!(i)
v0=[randn()*σv for j=1:N]

#the DiffEq slightly differs from the paper
#as we use different time scales
#and the equation for <σ^z> instead of <σ^ee>
function cavityall(du, u, p, t)  
    GU = 0
    for i=1:N
        GU+=G(x0[i]+t*v0[i])*u[i+1]
    end
    #<a>
    du[1]=-κ*u[1]+  K * GU
    #<σ->
    for i=2:(N+1)
        du[i]=(im*Δ-γ-ν)*u[i]+G(x0[i-1]+t*v0[i-1])*u[N+i]*u[1]
    end
    #<σ^z>
    for i=(N+2):(2*N+1)
        du[i]=-2*(γ+ν)*u[i]-2*G(x0[i-1-N]+t*v0[i-1-N])*(conj(u[1])*u[i-N]+conj(u[i-N])*u[1])+2*(ν-γ)
    end
end

s0=[0.0+0.0*im for i=1:N]
z0=[-1.0+0.0*im for i=1:N];
u0=[0.01+0.0*im; s0; z0];

tspan=(0.0, 400.0)
prob=ODEProblem(cavityall, u0, tspan)
sol=solve(prob, reltol=1e-10, abstol=1e-10, saveat=0.2);

L=length(sol.t)
expn=abs2.(sol[1, :])

return expn
end


function calcdata(N, K, σv, samplesize)
    open("data/data_N$N.sc$K.sv$σv.txt", "w") do io
        ini=0
        for j=(ini+1):(ini+samplesize)   
            writedlm(io, [photnumbevo.(N, K, σv, j)])
        end
    end
end

tablesc=[500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500]

samplesize=50;

scanv1=collect(0:0.01:0.1)
scanv2=collect(0.12:0.02:0.7)

scanv=vcat(scanv1, scanv2);

#using ProgressMeter
#prog=Progress(length(tablesc)*length(scanv))

Threads.@threads for j in 1:length(scanv)
	for i in 1:length(tablesc)
            calcdata(400, tablesc[i], scanv[j], samplesize)
        #next!(prog)
    end
end

open("done.txt", "w") do io
	writedlm(io, 0)
end

