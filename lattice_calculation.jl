#import necessary libraries
using DifferentialEquations
using DelimitedFiles
using Random
cd(@__DIR__)

function integratefrom(x, y, tstart)
    int=0
    istart=0
    for i=1:length(x)
        if x[i]<=tstart<x[i+1]
            istart=i
        end
    end
    ibx=[tstart; [x[i] for i=istart:length(x)]]
    iby=[y[istart]+(y[istart+1]-y[istart])*(tstart-x[istart]); [y[i] for i=istart:length(x)]]
    for i=1:(length(ibx)-1)
    int+=(iby[i]+iby[i+1])*(ibx[i+1]-ibx[i])*1/2
    #println("sum: $int \n")
    end
    return int
end

#define function which returns solution to set of differential equations
function photnumbevo(T_end, K, NCl, ν, σv, i, γ, g, ΔE, λlatt)

    κ=1
    Δ=0

    function G(x)
        λ=1
        k=(2*pi)/λ
        return g*cos(k*x)
    end
    #define atom trajectories
    #x0=[j/NCl for j=1:NCl]
    Random.seed!(i)
    x0=[j*λlatt for j=1:NCl]
    varray=randn(NCl)*σv
    function cavityall(du, u, p, t)
        GU = 0
        for i=1:NCl
            GU+=G(u[i+1+2*NCl])*u[i+1]
        end

        du[1]=-κ*u[1]+  K* GU

        for i=2:(1+NCl)
            du[i]=(im*Δ-γ-ν)*u[i]+G(u[i+2*NCl])*u[NCl+i]*u[1]
        end

        for i=(2+NCl):(1+2*NCl)
                du[i]=-2*(γ+ν)*u[i]-2*G(u[i+NCl])*(conj(u[1])*u[i-NCl]+conj(u[i-NCl])*u[1])+2*(ν-γ)
        end
        #now the atom movement
        for i=(2+2*NCl):(1+3*NCl)
            du[i]=u[i+NCl]
        end
        #atom velocities
        for i=(2+3*NCl):(1+4*NCl)
            du[i]=-ΔE*2*π*sin(2*π/λlatt*u[i-NCl])
        end
    end

    s00=[0.0+0.0*im for i=1:NCl];
    u0=[0.1+0.0*im; s00; s00; x0; varray];

    tspan=(0.0, T_end)
    prob=ODEProblem(cavityall, u0, tspan)
    sol=solve(prob, abstol=1e-12, reltol=1e-12);

    #pick solutions for interesting field operators
    return sol
    #return ibres
end


scanv1=collect(0.0:0.01:0.1)
scanv2=collect(0.12:0.02:0.7)
scanv3=collect(0.75:0.05:1.5)
scanv=2*vcat(scanv1, scanv2, scanv3);

samplesize=10

tableΔE=[0.0, 0.00001, 0.0001, 0.001, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3]
tableλlatt=[0.5, 0.99]

function findssval(T_end, K, NCl, ν, σv, γ, g, ΔE, λlatt, samplesize)
    finres=0
    for i=1:samplesize
        res=photnumbevo(T_end, K, NCl, ν, σv, i, γ, g, ΔE, λlatt)
        ibres=integratefrom(res.t, abs2.(res[1,:]), 200)/200
        finres+=ibres
    end
    return finres/samplesize
end

curveres=zeros(length(tableΔE), length(scanv))

#using ProgressMeter
#prog=Progress(length(tableΔE)*length(scanv)*length(tableλlatt))
for k=1:length(tableλlatt)
    for i=1:length(tableΔE)
        curveres=zeros(length(scanv))
        Threads.@threads for j=1:length(scanv)
            curveres[j]=findssval(400, 10000, 100, 0.5, scanv[j], 0, 2*0.00136, tableΔE[i], tableλlatt[k], samplesize)
            #next!(prog)
        end
        open("data/nmcurve.deltaE$(tableΔE[i]).lambdalatt$(tableλlatt[k]).txt", "w") do io
            writedlm(io, curveres)
        end
    end
end


open("data/done.txt", "w") do io
    writedlm(io, 0)
end