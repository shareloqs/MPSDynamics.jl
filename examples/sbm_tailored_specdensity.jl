#=
    Example of a zero-temperature Spin-Boson Model with a tailored spectral density J(ω). This spectral density can be either defined with discrete
    modes or continuous formula.

    The dynamics is simulated using the T-TEDOPA method that maps the normal modes environment into a non-uniform tight-binding chain.

    H = \\frac{ω_0}{2} σ_z + Δ σ_x + c_0 σ_x(b_0^\\dagger + b_0) + \\sum_{i=0}^{N-1} t_i (b_{i+1}^\\dagger b_i +h.c.) + \\sum_{i=0}^{N-1} ϵ_i b_i^\\dagger b_i 

    Two variants of the one-site Time Dependent Variational Principal (TDVP) are presented for the time evolution of the quantum state.
=#

using MPSDynamics, Plots, LaTeXStrings

const ∞  = Inf
#----------------------------
# Physical parameters
#----------------------------

ω0 = 0.2 # TLS gap

Δ = 0.0 # tunneling 

#-----------------------
# Options and parameters for the Spectral Density
# Definition of the frequency ranges and J(ω) formula for the continuous case
#-----------------------

#Jω_type=:Discrete

Jω_type=:Continuous

if Jω_type==:Discrete
    
    freqs = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50] # Frequency value
    
    Jω = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10] # Value of the spectral density at the respective frequency
    
    N = length(freqs) # The number of mode in the chain is set with the number of discrete frequencies
    
    d = 6 # number of Fock states of the chain modes

elseif Jω_type==:Continuous 

    N = 30 # length of the chain

    d = 6 # number of Fock states of the chain modes

    α = 0.1 # coupling strength
    
    s = 1 # ohmicity
    
    ωc = 1.0 # Cut-off of the spectral density J(ω) 
    
    #β = 100 # Thermalized environment
    β = ∞ # Case zero temperature T=0, β → ∞

    # AB segments correspond to the different frequency ranges where the spectral density is defined.
    # It corresponds to i=1 ; i=2 ; i=3 and i=4 of the Jω_fct(ω,i) function.

    # The y formula and the frequency ranges can be changed while keeping the extra terms for thermalized environments.
    # This example reproduces an Ohmic type spectral density. 

    AB = [[-Inf -ωc];[-ωc 0];[0 ωc];[ωc Inf]]

    function Jω_fct(ω,i)
        if i==1 # First segment of AB (here [-Inf -ωc])
            y = 0
        elseif i==2 # Second segment of AB (here [-ωc 0])
            if β == ∞
                y = 0 # zero for negative frequencies when β == Inf   
            else
                y = -2*α*abs.(ω).^s / ωc^(s-1) .* (coth.((β/2).*ω) .+ 1) ./2 # .* (coth.((β/2).*ω) .+ 1) ./2 and the abs.(ω) have to be added when β != Inf (T-TEDOPA formulation) 
            end
        elseif i==3 # Third segment of AB (here [0 ωc])
            if β == ∞
                y = 2*α*ω.^s / ωc^(s-1) 
            else
                y = 2*α*ω.^s / ωc^(s-1) .* (coth.((β/2).*ω) .+ 1) ./2  # .* (coth.((β/2).*ω) .+ 1) ./2 has to be added when β != Inf (T-TEDOPA formulation) 
            end       
        elseif i==4 # Fourth segment of AB (here [ωc Inf])
            y = 0
        end
        return y
    end
    
else

    throw(ErrorException("The spectral density has either to be Continuous or Discrete"))

end

#-----------------------
# Calculation of chain parameters
#-----------------------

if Jω_type==:Discrete

    cpars = chaincoeffs_finiteT_discrete(β, freqs, Jω; procedure=:Lanczos, Mmax=5000, save=false)  # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0

#= If cpars is stored in "../ChainOhmT/discreteT" 
    curdir = @__DIR__
    dir_chaincoeff = abspath(joinpath(curdir, "../ChainOhmT/discreteT"))
    cpars  = readchaincoeffs("$dir_chaincoeff/chaincoeffs.h5", N, β) # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0
=#

elseif Jω_type==:Continuous

    cpars = chaincoeffs_finiteT(N, β, false; α=α, s=s, J=Jω_fct, ωc=ωc, mc=size(AB)[1], mp=0, AB=AB, iq=1, idelta=2, procedure=:Lanczos, Mmax=5000, save=false)  # chain parameters, i.e. on-site energies ϵ_i, hopping energies t_i, and system-chain coupling c_0

end


#-----------------------
# Simulation parameters
#-----------------------

dt = 0.5 # time step

tfinal = 30.0 # simulation time

method = :TDVP1 # Regular one-site TDVP (fixed bond dimension)

# method = :DTDVP # Adaptive one-site TDVP (dynamically updating bond dimension)

convparams = [2,4,6] # MPS bond dimension (1TDVP)

# convparams = [1e-2, 1e-3, 1e-4] # threshold value of the projection error (DTDVP)

#---------------------------
# MPO and initial state MPS
#---------------------------

H = spinbosonmpo(ω0, Δ, d, N, cpars) # MPO representation of the Hamiltonian

ψ = unitcol(1,2) # Initial up-z system state 

A = productstatemps(physdims(H), state=[ψ, fill(unitcol(1,d), N)...]) # MPS representation of |ψ>|Vacuum>

#---------------------------
# Definition of observables
#---------------------------

ob1 = OneSiteObservable("sz", sz, 1)

ob2 = OneSiteObservable("chain mode occupation", numb(d), (2,N+1))

ob3 = TwoSiteObservable("SXdisp", sx, MPSDynamics.disp(d), [1], collect(2:N+1))

#-------------
# Simulation
#------------

A, dat = runsim(dt, tfinal, A, H;
                name = "ohmic spin boson model",
                method = method,
                obs = [ob2,ob3],
                convobs = [ob1],
                params = @LogParams(Δ, ω0, N, Jω_type),
                convparams = convparams,
                verbose = false,
                savebonddims = true, # this keyword argument enables the bond dimension at each time step to be saved when using DTDVP
                save = true,
                plot = true,
                );

#----------
# Plots
#----------

display(method == :TDVP1 && plot(dat["data/times"], dat["convdata/sz"], label=["Dmax = 2" "Dmax = 4" "Dmax = 6"], xlabel=L"t",ylabel=L"\sigma_z"))

method == :DTDVP && plot(dat["data/times"], dat["convdata/sz"], label=["p = 1e-2" "p = 1e-3" "p = 1e-4"], xlabel=L"t",ylabel=L"\sigma_z") 

method == :DTDVP && heatmap(dat["data/times"], collect(0:N+1), dat["data/bonddims"], xlabel=L"t",ylabel="bond index")

heatmap(dat["data/times"], collect(1:N), abs.(dat["data/SXdisp"][1,:,:]), xlabel=L"t",ylabel="chain mode")
