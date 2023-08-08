using DifferentialEquations
using DelimitedFiles
using Random
cd(@__DIR__)


#Cavity with moving atoms

function photnumbevo(T_end, scaling, N, ν, σv, i, γ, g, veldistcut)
    
    
κ=1
γ=0
Δ=0
#g=sqrt(0.0015/200)    #6.2*0.000001
#ν=0.5;  
  
#N=400    
   
function G(x)
    λ=1
    k=(2*pi)/λ
    return g*cos(k*x)
end
    
    
function myfilter(b1, cut)
 return filter(x -> abs(x)<cut, b1)
end
 
x0=[j/N for j=0:(N-1)]
Random.seed!(i)
v0_ini=[randn() for j=1:N]
v0=myfilter(v0_ini, veldistcut)*σv  
N=length(v0)    

function cavityall(du, u, p, t)  
    GU = 0
    for i=1:N
        GU+=G(x0[i]+t*v0[i])*u[i+1]
    end
                                 
  du[1]=-κ*u[1]+  scaling * GU
  for i=2:(N+1)
  	du[i]=(im*Δ-γ-ν)*u[i]+G(x0[i-1]+t*v0[i-1])*u[N+i]*u[1]
  end
            
  for i=(2+N):(2*N+1)
  	du[i]=-2*(γ+ν)*u[i]-2*G(x0[i-1-N]+t*v0[i-1-N])*(conj(u[1])*u[i-N]+conj(u[i-N])*u[1])+2*(ν-γ)
  end
end

s0=[0.0+0.0*im for i=1:N]
z0=[-1.0+0.0*im for i=1:N];
u0=[0.01+0.0*im; s0; z0];

tspan=(0.0, T_end)
prob=ODEProblem(cavityall, u0, tspan)
sol=solve(prob, abstol=1e-10, reltol=1e-10, saveat=0.1);

L=length(sol.t)
expn=abs2.(sol[1, :])

return expn
end

function calcdata(T_end, scalesize, N, ν, σv, g, samplesize, veldistcut)

    open("data/veldistcut_mean_Tend$T_end.scsi$scalesize.N$N.nu$ν.sv$σv.g$g.veldistcut$veldistcut.txt", "w") do io
        ini=1
        #ini=get_samplesize(N, σv)
        for j=(ini):(ini+samplesize-1)   
        writedlm(io, [photnumbevo.(T_end, scalesize, N, ν, σv, j, 0, g, veldistcut)])
        end
    end
end

tablescsi=[2500];
scanv1=collect(0.01:0.01:0.1)
scanv2=collect(0.12:0.02:0.7)
scanv=vcat(scanv1, scanv2);
samplesize=50;
tablecut=[1, 0.5];


#using ProgressMeter
#prog=Progress(length(tablescsi)*length(scanv)*length(tablecut))
Threads.@threads for j in 1:length(scanv)
	for i in 1:length(tablescsi)
		for k in 1:length(tablecut)
			calcdata(200, tablescsi[i], 400, 0.5, scanv[j], 0.00272, samplesize, tablecut[k])
			#next!(prog)
		end
	end
end


open("data/done.txt", "w") do io
        writedlm(io, 0)
end
