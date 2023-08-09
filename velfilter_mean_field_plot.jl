using PyPlot
using DelimitedFiles
cd(@__DIR__)

scanv1=collect(0.01:0.01:0.1)
scanv2=collect(0.12:0.02:0.7)
scanv=vcat(scanv1, scanv2);

function mean_te(T_end, scalesize, N, ν, σv, g)
 yprob=readdlm("veldistcut/mean_Tend$T_end.scsi$scalesize.N$N.nu$ν.sv$σv.g$g.txt");
 
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
 yprob=readdlm("veldistcut/veldistcut_mean_Tend$T_end.scsi$scalesize.N$N.nu$ν.sv$σv.g$g.veldistcut$veldistcut.txt");
 
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
pygui(true);

#################################################################################################################
# (10) plots

#paper plot
close("veldistcut"); figure("veldistcut", figsize=(12, 5))
subplot(121)
scale=1000;
plot(0.5*scanv, mean_graph(200, 2500, 400, 0.5, 0.0027, 0)/scale, label=L"\mathrm{no \;cut}")
plot(0.5*scanv, veldist_mean_graph(200, 2500, 400, 0.5, 0.0027, 1.0, 0)/scale, label=L"\mathrm{cut\; at\;} \sigma_\mathrm{v}")
plot(0.5*scanv, veldist_mean_graph(200, 2500, 400, 0.5, 0.0027, 0.5, 0)/scale, label=L"\mathrm{cut\; at\; } \sigma_\mathrm{v}/2")
#plot(0.5*scanv, veldist_mixed_graph(200, 2500, 400, 0.5, 0.0027, 0.5, 0)/scale, linestyle="dotted")
#plot(0.5*scanv, veldist_mixed_graph(200, 2500, 400, 0.5, 0.0027, 1.0, 0)/scale, linestyle="dotted")

grid()
legend()
annotate("(a)", xy=(0.05, 0.08), xycoords="axes fraction")
xlabel(L"\sigma_\mathrm{v}/(\lambda \kappa)")
ylabel(L"\mathrm{photon\; number}/10^3")
subplot(122)
plot(0.5*scanv, mean_graph(200, 2500, 400, 0.5, 0.0027), label="no cut")
plot(0.5*scanv, veldist_mean_graph(200, 2500, 400, 0.5, 0.0027, 1.0), label=L"\mathrm{cut\; at\; 1.0\;} \sigma_\mathrm{v}/(\lambda \kappa)")
plot(0.5*scanv, veldist_mean_graph(200, 2500, 400, 0.5, 0.0027, 0.5), label=L"\mathrm{cut\; at\; 0.5\;} \sigma_\mathrm{v}/(\lambda \kappa)")
grid()
annotate("(b)", xy=(0.05, 0.08), xycoords="axes fraction")
xlabel(L"\sigma_\mathrm{v}/(\lambda \kappa)")
ylabel("photon number (normalized)")

