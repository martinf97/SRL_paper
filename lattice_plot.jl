using DelimitedFiles
using PyPlot
pygui(true)
cd(@__DIR__)

scanv1=collect(0.0:0.01:0.1)
scanv2=collect(0.12:0.02:0.7)
scanv3=collect(0.75:0.05:1.5)
scanv=2*vcat(scanv1, scanv2, scanv3);



include("base_functions.jl")
function percentage_of_trapped_atoms(σv, ΔE)
    if σv==0
        if ΔE>0
            return 1
        elseif ΔE==0
            return 0
        end
    end
    xarr=linplace(0, 80*σv, 50000)
    res=gauss.(xarr, 0, σv)
    return res=1-2*integratefrom(xarr, res, sqrt(2*ΔE))
end

function find_fractions(xarr, yarr, frac, ΔE)
    trapfracarr=percentage_of_trapped_atoms.(xarr, ΔE)
    index, xval=findyintersect(xarr, trapfracarr, frac)
    yval=yarr[index]+(yarr[index+1]-yarr[index])/(xarr[index+1]-xarr[index])*(xval-xarr[index])
    return xval, yval
end


##########################################################################
##########################################################################
#paper plot
hbar=1.05457182*1e-34
uconv=1.66053906660*1e-27
#calcium
λconv=657*1e-9
κconv=7.3*1e5
mconv=40.078*uconv
#conversion factors
vconv=λconv^2*κconv*mconv/(hbar * 2* pi)
Econv=vconv^2


rc("text",  usetex="True")
rc("font", family="serif")
rc("font", size=30)
figure("lattice", figsize=(18, 8))
subplot(121)
    whichlambda=0.5
    #curve0000001=readdlm("data/nmcurve.deltaE1.0e-6.lambdalatt$whichlambda.txt")
    curve000001=readdlm("data/nmcurve.deltaE1.0e-5.lambdalatt$whichlambda.txt")
    curve00001=readdlm("data/nmcurve.deltaE0.0001.lambdalatt$whichlambda.txt")
    curve0001=readdlm("data/nmcurve.deltaE0.001.lambdalatt$whichlambda.txt")
    curve001=readdlm("data/nmcurve.deltaE0.01.lambdalatt$whichlambda.txt")
    curve002=readdlm("data/nmcurve.deltaE0.02.lambdalatt$whichlambda.txt")
    curve005=readdlm("data/nmcurve.deltaE0.05.lambdalatt$whichlambda.txt")
    curve01=readdlm("data/nmcurve.deltaE0.1.lambdalatt$whichlambda.txt")
    curve02=readdlm("data/nmcurve.deltaE0.2.lambdalatt$whichlambda.txt")
    curve03=readdlm("data/nmcurve.deltaE0.3.lambdalatt$whichlambda.txt")
    curve00=readdlm("data/nmcurve.deltaE0.0.lambdalatt$whichlambda.txt")
    title(L"$\lambda_\mathrm{latt}=%$whichlambda$") #works
    color="grey"
    lw=1
    alpha=0.5
    #ax1.title(L"$\lambda_\mathrm{latt}=%$whichlambda$") #works
    plot(0.5*scanv*vconv, curve00/1e3, color=color, linewidth=lw, alpha=alpha)
    #plot(0.5*scanv*vconv, curve00001/1e3, color=color, linewidth=lw, alpha=alpha)
    plot(0.5*scanv*vconv, curve0001/1e3, color=color, linewidth=lw, alpha=alpha)
    plot(0.5*scanv*vconv, curve001/1e3, color=color, linewidth=lw, alpha=alpha)
    plot(0.5*scanv*vconv, curve002/1e3, color=color, linewidth=lw, alpha=alpha)
    plot(0.5*scanv*vconv, curve005/1e3, color=color, linewidth=lw, alpha=alpha)
    plot(0.5*scanv*vconv, curve01/1e3, color=color, linewidth=lw, alpha=alpha)
    plot(0.5*scanv*vconv, curve02/1e3, color=color, linewidth=lw, alpha=alpha)
    plot(0.5*scanv*vconv, curve03/1e3, color=color, linewidth=lw, alpha=alpha)
    cmap_str="viridis_r"
    size_p=150
    scatter(0.5*scanv*vconv, curve03/1e3, marker="|", label=L"\Delta E/E_\mathrm{rec}=3 \cdot 10^2" , s=size_p, c=percentage_of_trapped_atoms.(0.5*scanv*vconv, vconv^2*0.3), vmin=0, vmax=1, cmap=cmap_str)
    scatter(0.5*scanv*vconv, curve02/1e3, marker=7, label=L"\Delta E/E_\mathrm{rec}=2 \cdot 10^2" , s=size_p, c=percentage_of_trapped_atoms.(0.5*scanv*vconv, vconv^2*0.2), vmin=0, vmax=1, cmap=cmap_str)
    scatter(0.5*scanv*vconv, curve01/1e3, marker="x", label=L"\Delta E/E_\mathrm{rec}=1 \cdot 10^2" , s=size_p, c=percentage_of_trapped_atoms.(0.5*scanv*vconv, vconv^2*0.1), vmin=0, vmax=1, cmap=cmap_str)
    scatter(0.5*scanv*vconv, curve005/1e3, marker="v", label=L"\Delta E/E_\mathrm{rec}=5 \cdot 10^1", s=size_p, c=percentage_of_trapped_atoms.(0.5*scanv*vconv, vconv^2*0.05), vmin=0, vmax=1, cmap=cmap_str)
    scatter(0.5*scanv*vconv, curve002/1e3, marker=".", label=L"\Delta E/E_\mathrm{rec}=2 \cdot 10^1", s=size_p, c=percentage_of_trapped_atoms.(0.5*scanv*vconv, vconv^2*0.02), vmin=0, vmax=1, cmap=cmap_str)
    scatter(0.5*scanv*vconv, curve001/1e3, marker="s", label=L"\Delta E/E_\mathrm{rec}=1 \cdot 10^1", s=size_p, c=percentage_of_trapped_atoms.(0.5*scanv*vconv, vconv^2*0.01), vmin=0, vmax=1, cmap=cmap_str)
    scatter(0.5*scanv*vconv, curve0001/1e3, marker="*", label=L"\Delta E/E_\mathrm{rec}=1", s=size_p, c=percentage_of_trapped_atoms.(0.5*scanv*vconv, vconv^2*0.001), vmin=0, vmax=1, cmap=cmap_str)
    #scatter(0.5*scanv*vconv, curve00001/1e3, marker="+", label=L"\Delta E/E_\mathrm{rec}=1", s=size_p, c=percentage_of_trapped_atoms.(0.5*scanv*vconv, vconv^2*0.0001), vmin=0, vmax=1, cmap=cmap_str)
    scatter(0.5*scanv*vconv, curve00/1e3, marker="o", label=L"\Delta E/E_\mathrm{rec}=0", s=size_p, c=zeros(length(scanv*vconv)), vmin=0, vmax=1, cmap=cmap_str)
    ylim(-10,250)
    xlabel(L"\sigma_\mathrm{v}/v_\mathrm{rec}")
    ylabel(L"\mathrm{photon\;number}/10^{3}")
    grid()
    annotate("(a)", xy=(0.01, 0.06), xycoords="axes fraction")
subplot(122)
    whichlambda=0.99
    #curve0000001=readdlm("data/nmcurve.deltaE1.0e-6.lambdalatt$whichlambda.txt")
    curve000001=readdlm("data/nmcurve.deltaE1.0e-5.lambdalatt$whichlambda.txt")
    curve00001=readdlm("data/nmcurve.deltaE0.0001.lambdalatt$whichlambda.txt")
    curve0001=readdlm("data/nmcurve.deltaE0.001.lambdalatt$whichlambda.txt")
    curve001=readdlm("data/nmcurve.deltaE0.01.lambdalatt$whichlambda.txt")
    curve002=readdlm("data/nmcurve.deltaE0.02.lambdalatt$whichlambda.txt")
    curve005=readdlm("data/nmcurve.deltaE0.05.lambdalatt$whichlambda.txt")
    curve01=readdlm("data/nmcurve.deltaE0.1.lambdalatt$whichlambda.txt")
    curve02=readdlm("data/nmcurve.deltaE0.2.lambdalatt$whichlambda.txt")
    curve03=readdlm("data/nmcurve.deltaE0.3.lambdalatt$whichlambda.txt")
    curve00=readdlm("data/nmcurve.deltaE0.0.lambdalatt$whichlambda.txt")
    title(L"$\lambda_\mathrm{latt}=%$whichlambda$") #works
    plot(0.5*scanv*vconv, curve00/1e3, color=color, linewidth=lw, alpha=alpha)
    #plot(0.5*scanv*vconv, curve00001/1e3, color=color, linewidth=lw, alpha=alpha)
    plot(0.5*scanv*vconv, curve0001/1e3, color=color, linewidth=lw, alpha=alpha)
    plot(0.5*scanv*vconv, curve001/1e3, color=color, linewidth=lw, alpha=alpha)
    plot(0.5*scanv*vconv, curve002/1e3, color=color, linewidth=lw, alpha=alpha)
    plot(0.5*scanv*vconv, curve005/1e3, color=color, linewidth=lw, alpha=alpha)
    plot(0.5*scanv*vconv, curve01/1e3, color=color, linewidth=lw, alpha=alpha)
    plot(0.5*scanv*vconv, curve02/1e3, color=color, linewidth=lw, alpha=alpha)
    plot(0.5*scanv*vconv, curve03/1e3, color=color, linewidth=lw, alpha=alpha)
    #cmap_str="viridis_r"
    scatter(0.5*scanv*vconv, curve03/1e3, marker="|", label=L"\Delta E/E_\mathrm{rec}=3 \cdot 10^2" , s=size_p, c=percentage_of_trapped_atoms.(0.5*scanv*vconv, vconv^2*0.3), vmin=0, vmax=1, cmap=cmap_str)
    scatter(0.5*scanv*vconv, curve02/1e3, marker=7, label=L"\Delta E/E_\mathrm{rec}=2 \cdot 10^2" , s=size_p, c=percentage_of_trapped_atoms.(0.5*scanv*vconv, vconv^2*0.2), vmin=0, vmax=1, cmap=cmap_str)
    scatter(0.5*scanv*vconv, curve01/1e3, marker="x", label=L"\Delta E/E_\mathrm{rec}=1 \cdot 10^2" , s=size_p, c=percentage_of_trapped_atoms.(0.5*scanv*vconv, vconv^2*0.1), vmin=0, vmax=1, cmap=cmap_str)
    scatter(0.5*scanv*vconv, curve005/1e3, marker="v", label=L"\Delta E/E_\mathrm{rec}=5 \cdot 10^1", s=size_p, c=percentage_of_trapped_atoms.(0.5*scanv*vconv, vconv^2*0.05), vmin=0, vmax=1, cmap=cmap_str)
    scatter(0.5*scanv*vconv, curve002/1e3, marker=".", label=L"\Delta E/E_\mathrm{rec}=2 \cdot 10^1", s=size_p, c=percentage_of_trapped_atoms.(0.5*scanv*vconv, vconv^2*0.02), vmin=0, vmax=1, cmap=cmap_str)
    scatter(0.5*scanv*vconv, curve001/1e3, marker="s", label=L"\Delta E/E_\mathrm{rec}=1 \cdot 10^1", s=size_p, c=percentage_of_trapped_atoms.(0.5*scanv*vconv, vconv^2*0.01), vmin=0, vmax=1, cmap=cmap_str)
    scatter(0.5*scanv*vconv, curve0001/1e3, marker="*", label=L"\Delta E/E_\mathrm{rec}=1", s=size_p, c=percentage_of_trapped_atoms.(0.5*scanv*vconv, vconv^2*0.001), vmin=0, vmax=1, cmap=cmap_str)
    #scatter(0.5*scanv*vconv, curve00001/1e3, marker="+", label=L"\Delta E/E_\mathrm{rec}=1", s=size_p, c=percentage_of_trapped_atoms.(0.5*scanv*vconv, vconv^2*0.0001), vmin=0, vmax=1, cmap=cmap_str)
    scatter(0.5*scanv*vconv, curve00/1e3, marker="o", label=L"\Delta E/E_\mathrm{rec}=0", s=size_p, c=zeros(length(scanv*vconv)), vmin=0, vmax=1, cmap=cmap_str)
    legend()
    ylim(-10,250)
    xlabel(L"\sigma_\mathrm{v}/v_\mathrm{rec}")
    grid()
    annotate("(b)", xy=(0.01, 0.06), xycoords="axes fraction")
    colorbar(label="fraction of trapped atoms")
    #colorbar()
    tight_layout()
