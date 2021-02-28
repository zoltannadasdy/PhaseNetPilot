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

pi_frac=2*pi/1000
x=zeros(1,1000)
y=zeros(1,1000)


for i in 1:1000
#PhaseD=angdiff.(0,i*pi_frac+(-1*pi))
#PhaseD=angdiff.(0,i*pi_frac)-pi
#PhaseD=angdiff.(0,i*pi_frac)
#PhaseD=1*exp(angdiff.(0,i*pi_frac))/pi/2
PhaseD=angdiff.(0,i*pi_frac)
x[1,i]=i*pi_frac+(-1*pi)
X=x[1,i]

#y[1,i]=(X-1)*log10((X-1)^2*.1).*exp(-.2*X^2)/32  # asymmetric
#y[1,i]=-X*log10((X-1)^2*.1).*exp(-.2*(X-1)^2)/32
#y[1,i]=-X*log10(abs(X)^2*.1)
y[1,i]=((1/(2*sqrt(2pi))*exp(-(1/2)*((X+1)/2)^2))/32)-0.004  # symmetrical U



#y[1,i]=-X*log10(X^2*.1)*exp(-.2*X^2) / 32
#y[1,i]=X*log10(X^2) / 2
end
replace!(y, NaN=>0.0)
replace!(y, Inf=>0.0)

#for i in 501:1000
##PhaseD=angdiff.(0,i*pi_frac+(-1*pi))
##PhaseD=angdiff.(0,i*pi_frac)-pi
#PhaseD=angdiff.(0,i*pi_frac)
#PhaseD=-1*exp(pi-angdiff.(0,i*pi_frac))/pi/2
#x[1,i]=i*pi_frac+(-1*pi)
#y[1,i]=PhaseD
#end


using Plots

plot(x[1,:],y[1,:])

current()
# Saving the SVG for futire editing
savefig("STDP-1-fn.svg")
