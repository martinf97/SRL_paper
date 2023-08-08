using PyPlot
pygui(true)
using DelimitedFiles
cd(@__DIR__)


scanv1=collect(0.01:0.01:0.1)
scanv2=collect(0.12:0.02:0.7)
scanv=vcat(scanv1, scanv2);


#reading data from file and processing


function mean_te(T_end, scalesize, N, ν, σv, g)
 yprob=readdlm("data/mean_Tend$T_end.scsi$scalesize.N$N.nu$ν.sv$σv.g$g.txt");
 
 te=sum(yprob, dims=1) 
return transpose(te)/size(yprob, 1)
end



function mean_av(T_end, scalesize, N, ν, σv, g)
    yprob=mean_te(T_end, scalesize, N, ν, σv, g)
    sum=0;
    start=1000;
    size=length(yprob)-start+1;
    for i=start:length(yprob)
        sum+=yprob[i]
    end
    
    return sum/size
end 


function mean_graph(T_end, scalesize, N, ν, g, norm=1) #plots the average photon number at N over σv
    graph=zeros(length(scanv))
    graph=mean_av.(T_end, scalesize, N, ν, scanv, g)
    
    res=0
    if norm==0
        res=graph
    end
    if norm==1
        res=graph/maximum(graph)
        end
#plot(scanv, graph/graph[1], label="$scalesize clustersize") #missing: graph/Natoms /graph[1]
return res 
end 

function veldist_mean_te(T_end, scalesize, N, ν, σv, g, veldistcut)
 yprob=readdlm("data/veldistcut_mean_Tend$T_end.scsi$scalesize.N$N.nu$ν.sv$σv.g$g.veldistcut$veldistcut.txt");
 
 te=sum(yprob, dims=1) 
return transpose(te)/size(yprob, 1)
end



function veldist_mean_av(T_end, scalesize, N, ν, σv, g, veldistcut)
    yprob=veldist_mean_te(T_end, scalesize, N, ν, σv, g, veldistcut)
    sum=0;
    start=1000;
    size=length(yprob)-start+1;
    for i=start:length(yprob)
        sum+=yprob[i]
    end
    
    return sum/size
end 


function veldist_mean_graph(T_end, scalesize, N, ν, g, veldistcut, norm=1) #plots the average photon number at N over σv
    graph=zeros(length(scanv))
    graph=veldist_mean_av.(T_end, scalesize, N, ν, scanv, g, veldistcut)
    
    res=0
    if norm==0
        res=graph
    end
    if norm==1
        res=graph/maximum(graph)
        end
#plot(scanv, graph/graph[1], label="$scalesize clustersize") #missing: graph/Natoms /graph[1]
return res 
end 


rc("text",  usetex="True")
rc("font", family="serif")
rc("font", size=18)

#################################################################################################################
# (10) plots

#paper plot
close("veldistcut"); figure("veldistcut", figsize=(12, 5))
subplot(121)
scale=1000;
#plot(0.5*scanv, mean_graph(200, 2500, 400, 0.5, 0.0027, 0)/scale, label=L"\mathrm{no \;cut}")
plot(0.5*scanv, veldist_mean_graph(200, 2500, 400, 0.5, 0.00272, 1.0, 0)/scale, label=L"\mathrm{cut\; at\;} \sigma_\mathrm{v}")
plot(0.5*scanv, veldist_mean_graph(200, 2500, 400, 0.5, 0.00272, 0.5, 0)/scale, label=L"\mathrm{cut\; at\; } \sigma_\mathrm{v}/2")
#plot(0.5*scanv, veldist_mixed_graph(200, 2500, 400, 0.5, 0.0027, 0.5, 0)/scale, linestyle="dotted")
#plot(0.5*scanv, veldist_mixed_graph(200, 2500, 400, 0.5, 0.0027, 1.0, 0)/scale, linestyle="dotted")

grid()
legend()
annotate("(a)", xy=(0.05, 0.08), xycoords="axes fraction")
xlabel(L"\sigma_\mathrm{v}/(\lambda \kappa)")
ylabel(L"\mathrm{photon\; number}/10^3")
subplot(122)
#plot(0.5*scanv, mean_graph(200, 2500, 400, 0.5, 0.00272), label="no cut")
plot(0.5*scanv, veldist_mean_graph(200, 2500, 400, 0.5, 0.00272, 1.0), label=L"\mathrm{cut\; at\; 1.0\;} \sigma_\mathrm{v}/(\lambda \kappa)")
plot(0.5*scanv, veldist_mean_graph(200, 2500, 400, 0.5, 0.00272, 0.5), label=L"\mathrm{cut\; at\; 0.5\;} \sigma_\mathrm{v}/(\lambda \kappa)")
grid()

#red_patch = Patch(color="red", label="The red data")
#legend()
annotate("(b)", xy=(0.05, 0.08), xycoords="axes fraction")

#title(L"$γ=10^{-8}, Δ=0, g=1.35 \cdot 10^{-3}, ν=0.5, κ=1$, 400 clusters")
xlabel(L"\sigma_\mathrm{v}/(\lambda \kappa)")
ylabel("photon number (normalized)")
#save veldistcut_norm_up.pdf
#tight_layout()

pygui(true)
scale=1000;
plot(0.5*scanv, mean_graph(200, 2500, 400, 0.5, 0.0027, 0)/scale, label=L"\mathrm{no \;cut}")
plot(0.5*scanv, veldist_mean_graph(200, 2500, 400, 0.5, 0.0027, 1.0, 0)/scale, label=L"\mathrm{cut\; at\; 1.0\;} \sigma_\mathrm{v}/(\lambda \kappa)")
plot(0.5*scanv, veldist_mean_graph(200, 2500, 400, 0.5, 0.0027, 0.5, 0)/scale, label=L"\mathrm{cut\; at\; 0.5\;} \sigma_\mathrm{v}/(\lambda \kappa)")
plot(0.5*scanv, veldist_mixed_graph(200, 2500, 400, 0.5, 0.0027, 0.5, 0)/scale, linestyle="dotted")
plot(0.5*scanv, veldist_mixed_graph(200, 2500, 400, 0.5, 0.0027, 1.0, 0)/scale, linestyle="dotted")

grid()
legend()
annotate("(a)", xy=(0.05, 0.08), xycoords="axes fraction")
xlabel(L"\sigma_\mathrm{v}/(\lambda \kappa)")
ylabel(L"\mathrm{photon\; number}/10^3")

#title(L"$γ=10^{-8}, Δ=0, g=1.35 \cdot 10^{-3}, ν=0.5, κ=1$, 400 clusters")

#save aveldistcut.pdf

pygui(true)
plot(0.5*scanv, mean_graph(200, 2500, 400, 0.5, 0.0027), label="no cut")
plot(0.5*scanv, veldist_mean_graph(200, 2500, 400, 0.5, 0.0027, 1.0), label=L"\mathrm{cut\; at\; 1.0\;} \sigma_\mathrm{v}/(\lambda \kappa)")
plot(0.5*scanv, veldist_mean_graph(200, 2500, 400, 0.5, 0.0027, 0.5), label=L"\mathrm{cut\; at\; 0.5\;} \sigma_\mathrm{v}/(\lambda \kappa)")
grid()

#red_patch = Patch(color="red", label="The red data")
legend()
annotate("(b)", xy=(0.05, 0.08), xycoords="axes fraction")

#title(L"$γ=10^{-8}, Δ=0, g=1.35 \cdot 10^{-3}, ν=0.5, κ=1$, 400 clusters")
xlabel(L"\sigma_\mathrm{v}/(\lambda \kappa)")
ylabel("photon number (normalized)")
#save veldistcut_norm_up.pdf

################################################################################################################

pygui(false)
plot(scanv, mean_graph(200, 2500, 400, 0.5, 0.0027, 0)/10^6, label="no cut")
plot(scanv, veldist_mean_graph(200, 2500, 400, 0.5, 0.0027, 1.0, 0)/norm10, label="cut at 1.0 σv")
plot(scanv, veldist_mean_graph(200, 2500, 400, 0.5, 0.0027, 0.5, 0)/norm05, label="cut at 0.5 σv")
grid()
legend()
xlabel("σv")
ylabel("photons per atom")

pygui(true)
plot(scanv, mean_graph(200, 2500, 400, 0.5, 0.0027, 0)/10^12, label="no cut")
plot(scanv, veldist_mean_graph(200, 2500, 400, 0.5, 0.0027, 1.0, 0)/norm10^2, label="cut at 1.0 σv")
plot(scanv, veldist_mean_graph(200, 2500, 400, 0.5, 0.0027, 0.5, 0)/norm05^2, label="cut at 0.5 σv")
grid()
legend()
xlabel("σv")
ylabel("photons per atomnumber^2")

function getnofatoms(i, veldistcut)
Random.seed!(i)
v0_ini=[randn() for j=1:400]
v0=myfilter(v0_ini, veldistcut)  
return length(v0)
end

function myfilter(b1, cut)
 return filter(x -> abs(x)<cut, b1)
end

arr10=[getnofatoms(i, 1.0) for i=1:50]
arr05=[getnofatoms(i, 0.5) for i=1:50];

norm10=sum(arr10)/50*2500

norm05=sum(arr05)/50*2500

b1=randn(10)

myfilter(b1, 1)

xarr=[1 2 3; 2 3 4; 3 4 5]

contourplot()

xarr=[5+i for i=1:10]
yarr=[20+3*i for i=1:15];

z=Array[]

zarr=zeros(length(yarr), length(xarr));

for i=1:15
    for j=1:10
        zarr[i, j]= xarr[j]*yarr[i]
    end
end

contourf(yarr, xarr, transpose(zarr))
legend()

zarr

transpose(zarr)

plot_surface(xarr, yarr, trazarr)


