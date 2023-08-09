using DelimitedFiles
cd(@__DIR__)
using PyPlot
pygui(true)

scanv1=collect(0:0.01:0.1)
scanv2=collect(0.12:0.02:0.7)
tarray=collect(0:0.1:200);

scanv=vcat(scanv1, scanv2);


#reading data from file and processing

function readout_te(scaling, σv)
    yprob=readdlm("data/data_N400.sc$scaling.sv$σv.txt");
    te=sum(yprob, dims=1) 
    return transpose(te)/size(yprob, 1)
end

function readout_av(scaling, σv)
    yprob=readout_te(scaling, σv)
    sum=0;
    start=1000;
    size=length(yprob)-start+1;
    for i=start:length(yprob)
        sum+=yprob[i]
    end
    
    return sum/size
end 

function readout_graph(scaling) #returns the average photon number at N over σv
    graph=zeros(length(scanv))
    graph=readout_av.(scaling, scanv)
return graph 
end 

function findroot(x, f, c)
    L=length(x)
    i=1
    r=0
    for i=1:L
        if f[i]<c
        r=(f[i-1]-c)*(x[i]-x[i-1])/(f[i-1]-f[i])+x[i-1];break
        end
    end
    return r
end

function readout_root(scaling)
    graph=zeros(length(scanv))
    graph=readout_av.(scaling, scanv)
    
    r=findroot(scanv, graph/maximum(graph), 0.5) #missing: graph/Natoms, limit 0.001
    return r
end

##############################################################################################################
# (5) plots
tablesc=[500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500]

yval=[0.5 for i=2:10]

graph1000=readout_graph(2500)
graph900=readout_graph(2250)
graph800=readout_graph(2000)
graph700=readout_graph(1750)
graph600=readout_graph(1500)
graph500=readout_graph(1250)
graph400=readout_graph(1000)
graph300=readout_graph(750)
graph200=readout_graph(500);

rc("text",  usetex="True")
rc("font", family="serif") #sans-serif
rc("font", size=18)

#paper plot
close("fixg"); figure("fixg")#, figsize=(12,5))
subplot(131)
scale=1000;
plot(0.5*scanv, graph1000/scale, label=L"N=1 \cdot 10^6")
plot(0.5*scanv, graph900/scale, label=L"N=9 \cdot 10^5")
plot(0.5*scanv, graph800/scale, label=L"N=8 \cdot 10^5")
plot(0.5*scanv, graph700/scale, label=L"N=7 \cdot 10^5")
plot(0.5*scanv, graph600/scale, label=L"N=6 \cdot 10^5")
plot(0.5*scanv, graph500/scale, label=L"N=5 \cdot 10^5")
plot(0.5*scanv, graph400/scale, label=L"N=4 \cdot 10^5")
plot(0.5*scanv, graph300/scale, label=L"N=3 \cdot 10^5")
plot(0.5*scanv, graph200/scale, label=L"N=2 \cdot 10^5")


#title(L"$γ=10^{-8}, Δ=0, g=1.36 \cdot 10^{-3}, ν=0.5, κ=1$, 400 clusters")
annotate("(a)", xy=(0.01, 0.06), xycoords="axes fraction")
xlabel(L"\sigma_\mathrm{v}/(\lambda \kappa)")
ylabel(L"\mathrm{photon\; number}/10^3")
legend()
grid()

subplot(132)
plot(0.5*scanv, graph1000/graph1000[1], label=L"N=1 \cdot 10^6")
plot(0.5*scanv, graph900/graph900[1], label=L"N=9 \cdot 10^5")
plot(0.5*scanv, graph800/graph800[1], label=L"N=8 \cdot 10^5")
plot(0.5*scanv, graph700/graph700[1], label=L"N=7 \cdot 10^5")
plot(0.5*scanv, graph600/graph600[1], label=L"N=6 \cdot 10^5")
plot(0.5*scanv, graph500/graph500[1], label=L"N=5 \cdot 10^5")
plot(0.5*scanv, graph400/graph400[1], label=L"N=4 \cdot 10^5")
plot(0.5*scanv, graph300/graph300[1], label=L"N=3 \cdot 10^5")
plot(0.5*scanv, graph200/graph200[1], label=L"N=2 \cdot 10^5")

scatter(0.5*readout_root.(tablesc), yval, color="r")
#title(L"\gamma=10^{-8}, \Delta=0, g=1.36 \cdot 10^{-3}, \; \nu=0.5, \; \kappa=1 ,\; \mathrm{400 \; clusters}")
annotate("(b)", xy=(0.01, 0.06), xycoords="axes fraction")
xlabel(L"\sigma_\mathrm{v}/(\lambda \kappa)")
ylabel("photon number (normalized)")
#legend()
grid()


#root
subplot(133)
pygui(true)
scatter(tablesc*400/1000000, (0.5*readout_root.(tablesc)).^2, color="r", label="average")
annotate("(c)", xy=(0.01, 0.06), xycoords="axes fraction")
ylabel(L"\sigma_\mathrm{v}^2/(\lambda \kappa)^2 \sim T \mathrm{-half-maximum}")
xlabel(L"N/10^6")
grid()