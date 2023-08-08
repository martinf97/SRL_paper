
function integrate(x, y)
    int=0
    for i=1:(length(x)-1)
    int+=(y[i]+y[i+1])*(x[i+1]-x[i])*1/2
    #println("sum: $int \n")
    end
    return int
end



function integratetill(x, y, tend)
    if last(x)==tend
        return integrate(x, y)
    end
    int=0
    iend=0
    for i=1:length(x)
        if x[i]<=tend<x[i+1]
            iend=i
        end
    end
    ibx=[[x[i] for i=1:iend]; tend]
    iby=[[y[i] for i=1:iend]; y[iend]+(y[iend+1]-y[iend])*(tend-x[iend])]
    for i=1:(length(ibx)-1)
    int+=(iby[i]+iby[i+1])*(ibx[i+1]-ibx[i])*1/2
    #println("sum: $int \n")
    end
    return int
end

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
    end
    return int
end

function findintersect(xarr, yarr, xval)
    ib=0
    res=0
    for i=1:length(xarr)
        if xarr[i]==xval
            return yarr[i]
        elseif xarr[i]<xval && xarr[i+1]>xval
            ib=i
        end
    end
    res=yarr[ib]+(yarr[ib+1]-yarr[ib])/(xarr[ib+1]-xarr[ib])*(xval-xarr[ib])
    return res
end

function findyintersect(xarr, yarr, yval)
    ib=0
    res=0
    for i=1:length(yarr)
        if yarr[i]==yval
            return xarr[i], yarr[i]
        elseif yarr[i]>yval && yarr[i+1]<yval
            ib=i
        end
    end
    
    res=xarr[ib]+(xarr[ib+1]-xarr[ib])/(yarr[ib+1]-yarr[ib])*(yval-yarr[ib])
    
    return ib, res
end



#######################################################################
function dimred(scanfracp, scanpl, res)
    xarr=[0.0]
    yarr=[0.0]
    for i=1:length(scanfracp)
        for j=1:length(scanpl)
            push!(xarr, scanfracp[i]*scanpl[j])
            push!(yarr, res[j, i])
        end
    end
    return xarr, yarr
end

function zoom(res, istart, jstart, length)
    ibres=zeros(length, length)
    for i=1:length
        for j=1:length
            ibres[i, j]=res[istart+i-1, jstart+j-1]
        end
    end
    return ibres
end

function findzeros(res)
    int=0
    for i=1:size(res, 1)
        for j=1:size(res, 2)
            if res[i, j]!=0
                int+=1
            end
        end
    end
    return int/length(res), int/length(res)+(size(res, 1)+size(res, 2)-1)/length(res)
end

function findzerosoph(res, tol=0)
    for i=1:size(res, 1)
        for j=1:size(res, 2)
            if res[i, j]<=tol
                println("0 at pos $i , $j")
            end
        end
    end
    return
end

function cutarr(x, a=1, b=length(x))
    res=[x[i] for i=a:b]
    return res
end

function cap(arr, a)
    arr<a ? ret=arr : ret=a
    return ret
end

function gauss(x, x0, σx)
    return 1/sqrt(2*pi*σx^2)*exp(-(x-x0)^2/(2*σx^2))
end

function find_cross(x1, x2, y1, y2, c)
    x=x1+(x2-x1)/(y2-y1)*(c-y1)
    return x
end

function find_lw(ω, spec)
    spec ./=maximum(spec)
    len=length(ω)
    x1=0
    x2=0
    for i1=1:len
        if spec[i1]>0.5
            x1=i1
            break
        end
    end
    i2_s=argmax(spec)
    for i2=i2_s:len
        if spec[i2]<0.5
            x2=i2
            break
        end
    end
    xl=find_cross(ω[x1-1],ω[x1],spec[x1-1],spec[x1],0.5)
    xr=find_cross(ω[x2-1],ω[x2],spec[x2-1],spec[x2],0.5)
    return xr-xl
end

path=pwd()
function save(xarr, name)
    open("$path/$name.txt", "w") do io
            writedlm(io, xarr)
    end
    return 0
end
function load(name)
    yprob=readdlm("$path/data/$name.txt")
    return yprob
end
function logplace(a, b, n)
    α=(n-1)/log(b/a)
    res=a*exp.(([i for i=1:n].-1)/α)
    res[n]=b
    return res
end
function linplace(a, b, c)
    k=(b-a)/(c-1)
    return collect(a:k:b)
end


string="done"
