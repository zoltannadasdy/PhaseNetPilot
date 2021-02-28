#cd("C:\\Users\\znada\\projects\\PhaseNetJulia")
cd("C:\\Users\\User\\home\\projects\\PhaseLearn\\PhaseNetJulia")
## parameter setting
n=64; m=64; # defione the dimensions of n x m cellular matrix
SurrDim=3; # define the size of the surrounding space
SurRad=Int(floor(SurrDim/2)) # half SurrDim as radius of surround
RunCycles=10; # define the number of time steps to simulate
M=rand(n,m)*2*pi.-pi; # initialize cellular matrix M with radom phase
# With the iPhasGrad you can tweak how long you want the network remember the input
iPhaseGrad=10 # inverse phase gradient must be a positive integer [1 ...n]. The larger the iPhasGrad the smaller the phase influence is
#M=round.(rand(n,m))*pi.-pi/2; # initialize M with radom phase
StimDim=10;
fontsize=12 # to write the time counter on the plots
#StimStart=[80 160]
#StimStp=[120 200]
StimStart=[20 120]
StimStop=[40 140]
l=StimStart[1];
#StimStart=[40]
#StimStop=[60]
# allocation of arrays and matrices
#STIM=round.(rand(StimDim,StimDim))*pi.-(pi/2).+0.2; # define a random pattern that serves as a stimulus
STIM=round.(rand(StimDim,StimDim))*2pi.-(pi/2); # define a random pattern that serves as a stimulus
TM=zeros(n,m); # a temporary n x m container
TimeCourseM=zeros(RunCycles,n,m); # the history of states of the cellular matrix
TimeCourseXC=zeros(RunCycles,n,m); # the history of states of the cellular matrix
TimeCourseXXC=zeros(RunCycles,n,m); # the history of states of the cellular matrix
N=zeros(SurrDim,SurrDim);
NM=zeros(SurrDim,SurrDim);
PhaseD=zeros(SurrDim,SurrDim);
MM=zeros(SurrDim,SurrDim);

# defining the layout
using Plots
gr()
using PerceptualColourMaps

# Functions used
# compute the angulare fifference between two phases
function angdiff(a, b)
    r = (b - a) % pi
    if r ≥ pi
        r -= 2*pi
    end
    return r
end

function wrapToPi(a)
r=a;
    if a > pi
        r = -pi+(a % pi)
    elseif a < (-pi)
        r = pi-(-a % pi)
    end
    return r
end

function learning_rule(PhaseD)
#   dPhase=(PhaseD.*log10.(PhaseD.^2).*exp(-1*PhaseD.^2))/16
   dPhase=(PhaseD.*log10.(PhaseD.^2*0.1).*exp(-0.2*PhaseD.^2))/16
   return dPhase
end

# This function implements the plasticity (learning rule)
function plasticity(N,PhaseD,CntrPhase)
    UpdateMat=zeros(SurrDim,SurrDim)
    #####
# #    UpdateMat=wrapToPi.(N.-(PhaseD.*log10.(PhaseD.^2).*exp(-1*PhaseD.^2))/8)
    UpdateMat=wrapToPi.(N.-learning_rule(PhaseD))

#    UpdateMat=wrapToPi.(N.-((PhaseD.+1).*log10.((PhaseD.+1).^2*0.1).*exp(-0.2*PhaseD.^2))/32)
#    UpdateMat=wrapToPi.(N.-((PhaseD.+1).*log10.((PhaseD.+1).^2*0.1).*exp(-0.2*PhaseD.^2))/32) # STDP function positive first
#    UpdateMat=wrapToPi.(((1/(.8*sqrt(2pi))*exp(-(1/2).*((PhaseD.+1)./.8)^2))./32).-0.004) # symmetrical inverted U
    #####
    UpdateMat[SurRad+1,SurRad+1]=CntrPhase   # reload the original center value which became NaN after applying the STDP function
    replace!(UpdateMat, NaN=>0.0)
    replace!(UpdateMat, Inf=>0.0)
    return UpdateMat
end

## Display learning rule function
function disp_funct()
pi_frac=2*pi/1000
x=zeros(1,1000)
y=zeros(1,1000)

    for i in 1:1000
        PhaseD=angdiff.(0,i*pi_frac)
        x[1,i]=i*pi_frac+(-1*pi)
        X=x[1,i]
        y[1,i]=learning_rule(X)  # asymmetric
    end
    replace!(y, NaN=>0.0)
    replace!(y, Inf=>0.0)
    plot(x[1,:],y[1,:])
    current()
    # Saving the SVG for futire editing
    savefig("STDP-1-fn.svg")
end


## Main loop
# construct the cell grid and run the state machine
anim = @animate for k in 0:RunCycles
if k==0
    annotate!(50,  50, text(string(k), :left, fontsize, color="white"))
    heatmap!(M,show = true, c = cmap("C8"), colorbar = true, clims=(-pi, pi))
    annotate!(50,  50, text(string(k), :left, fontsize, color="black"))  # put the counter on
else
    for i in SurRad+1:n-SurRad
        for j in SurRad+1:m-SurRad
            N=M[i-SurRad:i+SurRad,j-SurRad:j+SurRad]       # N is a copy of a sliding subsmatrix of M
            MM=zeros(SurrDim,SurrDim)  # MM will represent the effect of the central cell on the surround
            MM=MM.+M[i,j]               # load the center to the surround matrix
#            PhaseD=angdiff.(N[:,:],MM[:,:])       # we compute the phase difference
            PhaseD=angdiff.(N[:,:],MM[:,:])       # we compute the phase difference
#            NM=wrapToPi.(N.-(PhaseD./(iPhaseGrad*pi))) # NM is original content of the submatrix updated by the phase difference effect


#            NM=wrapToPi.(N.-(PhaseD.*log10.(PhaseD.^2).*exp(-1*PhaseD.^2))/8) # STDP function stabilize and freeze (hasa long lasting memory)

#            NM=wrapToPi.(N.-(PhaseD.*log10.(PhaseD.^2*0.1).*exp(-0.2*PhaseD.^2))/32) # STDP function positive first
#             NM=wrapToPi.(N.-((PhaseD.+1).*log10.((PhaseD.+1).^2*0.1).*exp(-0.2*PhaseD.^2))/32) # STDP function positive first

#             NM=wrapToPi.(1/(2*sqrt(2pi))*exp(-(1/2)).*((PhaseD.-0)./2)^2) # symmetrical U
#             NM=wrapToPi.(((1/(.8*sqrt(2pi))*exp(-(1/2).*((PhaseD.+1)./.8)^2))./32).-0.004) # symmetrical inverted U

#            NM=wrapToPi.(N.-(-PhaseD.*log10.(PhaseD.^2*0.1).*exp(-0.2*PhaseD.^2))/32) # STDP function negative first

#            NM=wrapToPi.(N.-(PhaseD.*log10.((PhaseD.^2).*0.00001).*exp(-1*PhaseD.^2))/80) # STDP function
#            NM[SurRad+1,SurRad+1]=M[i,j]   # reload the original center value which became NaN after applying the STDP function
#            replace!(NM, NaN=>0.0)
#            replace!(NM, Inf=>0.0)
            NM=plasticity(N,PhaseD,M[i,j])
            TM[i-SurRad:i+SurRad,j-SurRad:j+SurRad]=NM
            M[i-SurRad:i+SurRad,j-SurRad:j+SurRad]=NM    # temporarily storing the new values of the matrix
        end
    end
# injecting a zero phase block when 80<k<100 and make it drifting across the cell grid.
for ii in 1:size(StimStart,2)
    if k>StimStart[ii] && k<StimStop[ii]
#        M[2+k-80:2+StimDim-1+k-80,2+k-80:2+StimDim-1+k-80]=STIM; # introducing the stimulus
    if ii < 2
        M[2+l:2+l+StimDim-1,(2+l):(2+StimDim-1+l)]=STIM[1:StimDim,1:StimDim];
#        M[2+l:2+l+StimDim-1,m-(2+StimDim-1+l):n-(2+l)]=STIM[l:l+StimDim-1,l:l+StimDim-1]; # drifting from lower left corner to the center
#    else
#        M[2+l:2+l+StimDim-1,(2+l):(2+StimDim-1+l)]=STIM[1:StimDim,1:StimDim]; # drifting from lower left corner to the center
    end
#        M[60-(2+StimDim-1+k-80):60-(2+k-80),2+k-80:2+StimDim-1+k-80]=STIM; # drifting from top left to lowe right
    end
end

    #counter (white to erase the trace of previous number)
    annotate!(50,  50, text(string(k-1), :left, fontsize, color="white"))
    # --------------------
    # Store the state in a 3-dim matrix, time (1-dim) and cells (2-dim)
    TimeCourseM[k,:,:]=M
    # --------------------

    #    heatmap!(p,M,show = true, c = :hsv, colorbar = true, clims=(-pi, pi), subplot = 1)
    heatmap!(M,show = true, c = cmap("C8"), colorbar = true, clims=(-pi, pi))
    annotate!(50,  50, text(string(k), :left, fontsize, color="black"))  # put the counter on
    # cmap needs to be precompiled
end
    current()
    println("t = $k")
end
# save animation to a GIF
gif(anim, "out4slow10x-neighbor3_slow8x_3stim-STDP6-invU_static_stim.gif", fps=15)

current()
