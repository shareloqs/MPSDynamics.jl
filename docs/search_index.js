var documenterSearchIndex = {"docs":
[{"location":"#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"The MPSDynamics.jl package provides an easy to use interface for performing tensor network simulations on matrix product states (MPS) and tree tensor network (TTN) states. Written in the Julia programming language, MPSDynamics.jl is a versatile open-source package providing a choice of several variants of the Time-Dependent Variational Principle (TDVP) method for time evolution.  The package also provides strong support for the measurement of observables, as well as the storing and logging of data, which makes it a useful tool for the study of many-body physics.  The package has been developed with the aim of studying non-Markovian open system dynamics at finite temperature using the state-of-the-art numerically exact Thermalized-Time Evolving Density operator with Orthonormal Polynomials Algorithm (T-TEDOPA) based on environment chain mapping. However the methods implemented can equally be applied to other areas of physics.","category":"page"},{"location":"#Table-of-Contents","page":"Introduction","title":"Table of Contents","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"","category":"page"},{"location":"#Installation","page":"Introduction","title":"Installation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"The package may be installed by typing the following into a Julia REPL","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"    ] add https://github.com/shareloqs/MPSDynamics.git","category":"page"},{"location":"#Functions","page":"Introduction","title":"Functions","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Modules = [MPSDynamics]","category":"page"},{"location":"#MPSDynamics.OneSiteObservable-Tuple{Any, Any, Any}","page":"Introduction","title":"MPSDynamics.OneSiteObservable","text":"OneSiteObservable(name,op,sites)\n\nComputes the local expectation value of the one-site operator op on the specified sites. Used to define one-site observables that are obs and convobs parameters for the runsim function.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.MPOtoVector-Tuple{ITensors.MPO}","page":"Introduction","title":"MPSDynamics.MPOtoVector","text":"MPOtoVector(mpo::MPO)\n\nConvert an ITensors chain MPO into a form compatible with MPSDynamics\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.addchild!-Tuple{MPSDynamics.Tree, Int64}","page":"Introduction","title":"MPSDynamics.addchild!","text":"addchild!(tree::Tree, id::Int)\n\nAdd child to node id of tree.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.addchildren!-Tuple{MPSDynamics.Tree, Int64, Int64}","page":"Introduction","title":"MPSDynamics.addchildren!","text":"addchildren!(tree::Tree, id::Int, n::Int)\n\nAdd n children to node id of tree.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.chaincoeffs_ohmic-Tuple{Any, Any, Any}","page":"Introduction","title":"MPSDynamics.chaincoeffs_ohmic","text":"chaincoeffs_ohmic(N, α, s; ωc=1, soft=false)\n\nGenerate chain coefficients ϵ_0ϵ_1t_0t_1c_0 for an Harmonic bath at zero temperature with a power law spectral density given by:\n\nsoft cutoff: J(ω) = 2παω_c (fracωω_c)^s exp(-ωω_c) \n\nhard cutoff: J(ω) = 2παω_c (fracωω_c)^s θ(ω-ω_c)\n\nThe coefficients parameterise the chain Hamiltonian\n\nH = H_S + c_0 A_SB_0+sum_i=0^N-1t_i (b_i+1^dagger b_i +hc) + sum_i=0^N ϵ_ib_i^dagger b_i\n\nwhich is unitarily equivalent (before the truncation to N sites) to\n\nH = H_S + A_Sint_0^dωsqrtfracJ(ω)πB_ω + int_0^dωωb_ω^dagger b_ω\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.chainmps-Tuple{Int64, Int64, Int64}","page":"Introduction","title":"MPSDynamics.chainmps","text":"chainmps(N::Int, site::Int, numex::Int)\n\nGenerate an MPS with numex excitations on site\n\nThe returned MPS will have bond-dimensions and physical dimensions numex+1\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.chainmps-Tuple{Int64, Vector{Int64}, Int64}","page":"Introduction","title":"MPSDynamics.chainmps","text":"chainmps(N::Int, sites::Vector{Int}, numex::Int)\n\nGenerate an MPS with numex excitations of an equal super-position over sites\n\nThe returned MPS will have bond-dimensions and physical dimensions numex+1\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.chainprop-Tuple{Any, Any}","page":"Introduction","title":"MPSDynamics.chainprop","text":"chainprop(t, cparams...)\n\nPropagate an excitation placed initially on the first site of a tight-binding chain with parameters given by cparams for a time t and return occupation expectation for each site.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.dynamap-NTuple{4, Any}","page":"Introduction","title":"MPSDynamics.dynamap","text":"dynamap(ps1,ps2,ps3,ps4)\n\nCalulate complete dynamical map to time step at which ps1, ps2, ps3 and ps4 are specified.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.electron2kmps","page":"Introduction","title":"MPSDynamics.electron2kmps","text":"electronkmps(N::Int, k::Vector{Int}, spin=:Up, chainparams=[fill(1.0,N), fill(1.0,N-1)])\n\nGenerate an MPS with 2 electrons in k-states k1 and k2.\n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.electronkmps","page":"Introduction","title":"MPSDynamics.electronkmps","text":"electronkmps(N::Int, k::Int, spin=:Up, chainparams=[fill(1.0,N), fill(1.0,N-1)])\n\nGenerate an MPS for an electron with momentum k.\n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.elementmpo-Tuple{Any, Vararg{Any}}","page":"Introduction","title":"MPSDynamics.elementmpo","text":"elementmpo(M, el...)\n\nReturn the element of the MPO M for the set of physical states el...\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.elementmps-Tuple{Any, Vararg{Any}}","page":"Introduction","title":"MPSDynamics.elementmps","text":"elementmps(A, el...)\n\nReturn the element of the MPS A for the set of physical states el...\n\nExamples\n\njulia> A = chainmps(6, [2,4], 1);\n\njulia> elementmps(A, 1, 2, 1, 1, 1, 1)\n0.7071067811865475\n\njulia> elementmps(A, 1, 1, 1, 2, 1, 1)\n0.7071067811865475\n\njulia> elementmps(A, 1, 2, 1, 2, 1, 1)\n0.0\n\njulia> elementmps(A, 1, 1, 1, 1, 1, 1)\n0.0\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.entanglemententropy-Tuple{Any}","page":"Introduction","title":"MPSDynamics.entanglemententropy","text":"entanglemententropy(A)\n\nFor a list of tensors A representing a right orthonormalized MPS, compute the entanglement entropy for a bipartite cut for every bond.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.findchainlength-Tuple{Any, Any}","page":"Introduction","title":"MPSDynamics.findchainlength","text":"findchainlength(T, cparams...; eps=10^-6)\n\nEstimate length of chain required for a particular set of chain parameters by calulating how long an excitation on the first site takes to reach the end. The chain length is given as the length required for the excitation to have just reached the last site after time T.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.findchild-Tuple{MPSDynamics.TreeNode, Int64}","page":"Introduction","title":"MPSDynamics.findchild","text":"findchild(node::TreeNode, id::Int)\n\nReturn integer corresponding to the which number child site id is of node.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.hbathchain-Tuple{Int64, Int64, Any, Vararg{Any}}","page":"Introduction","title":"MPSDynamics.hbathchain","text":"hbathchain(N::Int, d::Int, chainparams, longrangecc...; tree=false, reverse=false, coupletox=false)\n\nCreate an MPO representing a tight-binding chain of N oscillators with d Fock states each. Chain parameters are supplied in the standard form: chainparams =ϵ_0ϵ_1t_0t_1c_0. The output does not itself represent a complete MPO but will possess an end which is open and should be attached to another tensor site, usually representing the system.\n\nArguments\n\nreverse: If reverse=true create a chain were the last (i.e. Nth) site is the site which couples to the system\ncoupletox: Used to choose the form of the system coupling. coupletox=true gives a non-number conserving coupling of the form H_textI= A_textS(b_0^dagger + b_0) where A_textS is a system operator, while coupletox=false gives the number-converving coupling H_textI=(A_textS b_0^dagger + A_textS^dagger b_0)\ntree: If true the resulting chain will be of type TreeNetwork; useful for construcing tree-MPOs \n\nExample\n\nOne can constuct a system site tensor to couple to a chain by using the function up to populate the tensor. For example, to construct a system site with Hamiltonian Hs and coupling operator As, the system tensor M is constructed as follows for a non-number conserving interaction:\n\nu = one(Hs) # system identity\nM = zeros(1,3,2,2)\nM[1, :, :, :] = up(Hs, As, u)\n\nThe full MPO can then be constructed with:\n\nHmpo = [M, hbathchain(N, d, chainparams, coupletox=true)...]\n\nSimilarly for a number conserving interaction the site tensor would look like:\n\nu = one(Hs) # system identity\nM = zeros(1,4,2,2)\nM[1, :, :, :] = up(Hs, As, As', u)\n\nAnd the full MPO would be\n\nHmpo = [M, hbathchain(N, d, chainparams; coupletox=false)...]\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.isingmpo-Tuple{Int64}","page":"Introduction","title":"MPSDynamics.isingmpo","text":"isingmpo(N; J=1.0, h=1.0)\n\nReturn the MPO representation of a N-spin 1D Ising model with external field vech = (00h).\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.measure-Tuple{Any, OneSiteObservable}","page":"Introduction","title":"MPSDynamics.measure","text":"measure(A, O; kwargs...)\n\nMeasure observable O on mps state A\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.measure1siteoperator-Tuple{Vector, Any, Vector{Int64}}","page":"Introduction","title":"MPSDynamics.measure1siteoperator","text":"measure1siteoperator(A::Vector, O, sites::Vector{Int})\n\nFor a list of tensors A representing a right orthonormalized MPS, compute the local expectation value of a one-site operator O for every site or just one if it is specified.\n\nFor calculating operators on single sites this will be more efficient if the site is on the left of the mps.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.measure1siteoperator-Tuple{Vector, Any}","page":"Introduction","title":"MPSDynamics.measure1siteoperator","text":"measure1siteoperator(A::Vector, O)\n\nFor a list of tensors A representing a right orthonormalized MPS, compute the local expectation value of a one-site operator O for every site.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.measure2siteoperator-Tuple{Vector, Any, Any, Int64, Int64}","page":"Introduction","title":"MPSDynamics.measure2siteoperator","text":" measure2siteoperator(A::Vector, M1, M2, j1, j2)\n\nCaculate expectation of M1*M2 where M1 acts on site j1 and M2 acts on site j2, assumes A is right normalised.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.modemps","page":"Introduction","title":"MPSDynamics.modemps","text":"modemps(N::Int, k::Vector{Int}, numex::Int, chainparams=[fill(1.0,N), fill(1.0,N-1)])\n\nGenerate an MPS with numex excitations of an equal superposition of modes k of a bosonic tight-binding chain.\n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.modemps-2","page":"Introduction","title":"MPSDynamics.modemps","text":"modemps(N::Int, k::Int, numex::Int, chainparams=[fill(1.0,N), fill(1.0,N-1)])\n\nGenerate an MPS with numex excitations of mode k of a bosonic tight-binding chain. \n\nchainparams takes the form [e::Vector, t::Vector] where e are the on-site energies and t are the hoppping parameters.\n\nThe returned MPS will have bond-dimensions and physical dimensions numex+1\n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.mpsembed!-Tuple{Vector, Int64}","page":"Introduction","title":"MPSDynamics.mpsembed!","text":"mpsembed(A::Vector, Dmax::Int)\n\nEmbed MPS A in manifold of max bond-dimension Dmax\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.mpsleftnorm!","page":"Introduction","title":"MPSDynamics.mpsleftnorm!","text":"mpsleftnorm!(A::Vector, jq::Int=length(A))\n\nLeft orthoganalise MPS A up to site jq.\n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.mpsmixednorm!-Tuple{MPSDynamics.TreeNetwork, Int64}","page":"Introduction","title":"MPSDynamics.mpsmixednorm!","text":"mpsmixednorm!(A::TreeNetwork, id::Int)\n\nNormalise tree-MPS A such that orthogonality centre is on site id.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.mpsmixednorm!-Tuple{Vector, Int64}","page":"Introduction","title":"MPSDynamics.mpsmixednorm!","text":"mpsmixednorm!(A::Vector, OC::Int)\n\nPut MPS A into mixed canonical form with orthogonality centre on site OC.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.mpsmoveoc!-Tuple{MPSDynamics.TreeNetwork, Int64}","page":"Introduction","title":"MPSDynamics.mpsmoveoc!","text":"mpsmoveoc!(A::TreeNetwork, id::Int)\n\nMove the orthogonality centre of right normalised tree-MPS A to site id.\n\nThis function will be more efficient than using mpsmixednorm! if the tree-MPS is already right-normalised.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.mpsrightnorm!","page":"Introduction","title":"MPSDynamics.mpsrightnorm!","text":"mpsrightnorm!(A::Vector, jq::Int=1)\n\nRight orthoganalise MPS A up to site jq.\n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.mpsrightnorm!-Tuple{MPSDynamics.TreeNetwork}","page":"Introduction","title":"MPSDynamics.mpsrightnorm!","text":"mpsrightnorm!(A::TreeNetwork)\n\nWhen applied to a tree-MPS, right normalise towards head-node.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.mpsshiftoc!-Tuple{MPSDynamics.TreeNetwork, Int64}","page":"Introduction","title":"MPSDynamics.mpsshiftoc!","text":"mpsshiftoc!(A::TreeNetwork, newhd::Int)\n\nShift the orthogonality centre by one site, setting new head-node newhd.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.normmps-Tuple{MPSDynamics.TreeNetwork}","page":"Introduction","title":"MPSDynamics.normmps","text":"normmps(net::TreeNetwork; mpsorthog=:None)\n\nWhen applied to a tree-MPS mpsorthog=:Left is not defined.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.normmps-Tuple{Vector}","page":"Introduction","title":"MPSDynamics.normmps","text":"normmps(A::Vector; mpsorthog=:None)\n\nCalculate norm of MPS A.\n\nSetting mpsorthog=:Right/:Left will calculate the norm assuming right/left canonical form. Setting mpsorthog=OC::Int will cause the norm to be calculated assuming the orthoganility center is on site OC. If mpsorthog is :None the norm will be calculated as an MPS-MPS product.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.orthcentersmps-Tuple{Vector}","page":"Introduction","title":"MPSDynamics.orthcentersmps","text":"orthcentersmps(A)\n\nCompute the orthoganality centres of MPS A.\n\nReturn value is a list in which each element is the corresponding site tensor of A with the orthoganility centre on that site. Assumes A is right normalised.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.physdims-Tuple{Vector}","page":"Introduction","title":"MPSDynamics.physdims","text":"physdims(M)\n\nReturn the physical dimensions of an MPS or MPO M.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.productstatemps","page":"Introduction","title":"MPSDynamics.productstatemps","text":"productstatemps(physdims::Dims, Dmax=1; state=:Vacuum, mpsorthog=:Right)\n\nReturn an MPS representing a product state with local Hilbert space dimensions given by physdims.\n\nBy default all bond-dimensions will be 1 since the state is a product state. However, to embed the product state in a manifold of greater bond-dimension, Dmax can be set accordingly.\n\nThe indvidual states of the MPS sites can be provided by setting state to a list of column vectors. Setting state=:Vacuum will produce an MPS in the vacuum state (where the state of each site is represented by a column vector with a 1 in the first row and zeros elsewhere). Setting state=:FullOccupy will produce an MPS in which each site is fully occupied (ie. a column vector with a 1 in the last row and zeros elsewhere).\n\nThe argument mpsorthog can be used to set the gauge of the resulting MPS.\n\nExample\n\njulia> ψ = unitcol(1,2); d = 6; N = 30; α = 0.1; Δ = 0.0; ω0 = 0.2; s = 1\n\njulia> cpars = chaincoeffs_ohmic(N, α, s)\n\njulia> H = spinbosonmpo(ω0, Δ, d, N, cpars)\n\njulia> A = productstatemps(physdims(H), state=[ψ, fill(unitcol(1,d), N)...]) # MPS representation of |ψ>|Vacuum>\n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.productstatemps-2","page":"Introduction","title":"MPSDynamics.productstatemps","text":"productstatemps(N::Int, d::Int, Dmax=1; state=:Vacuum, mpsorthog=:Right)\n\nReturn an N-site MPS with all local Hilbert space dimensions given by d. \n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.randisometry-Tuple{Type, Int64, Int64}","page":"Introduction","title":"MPSDynamics.randisometry","text":"randisometry([T=Float64], dims...)\n\nConstruct a random isometry\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.randmps","page":"Introduction","title":"MPSDynamics.randmps","text":"randmps(N::Int, d::Int, Dmax::Int, T=Float64)\n\nConstruct a random, N-site, right-normalised MPS with all local Hilbert space dimensions given by d.\n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.randmps-2","page":"Introduction","title":"MPSDynamics.randmps","text":"randmps(tree::Tree, physdims, Dmax::Int, T::Type{<:Number} = Float64)\n\nConstruct a random, right-normalised, tree-MPS, with structure given by tree and max bond-dimension given by Dmax.\n\nThe local Hilbert space dimensions are specified by physdims which can either be of type Dims{length(tree)}, specifying the dimension of each site, or of type Int, in which case the same local dimension is used for every site.\n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.randmps-Union{Tuple{N}, Tuple{Tuple{Vararg{Int64, N}}, Int64}, Tuple{Tuple{Vararg{Int64, N}}, Int64, Type{<:Number}}} where N","page":"Introduction","title":"MPSDynamics.randmps","text":"randmps(physdims::Dims{N}, Dmax::Int, T::Type{<:Number} = Float64) where {N}\n\nConstruct a random, right-normalised MPS with local Hilbert space dimensions given by physdims and max bond-dimension given by Dmax. \n\nT specifies the element type, eg. use T=ComplexF64 for a complex valued MPS.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.randtree-Tuple{Int64, Int64}","page":"Introduction","title":"MPSDynamics.randtree","text":"randtree(numnodes::Int, maxdegree::Int)\n\nConstruct a random tree with nummodes modes and max degree maxdegree.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.rmsd-Tuple{Any, Any}","page":"Introduction","title":"MPSDynamics.rmsd","text":"rmsd(dat1::Vector{Float64}, dat2::Vector{Float64})\n\nCalculate the root mean squared difference between two measurements of an observable over the same time period.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.spinbosonmpo-NTuple{5, Any}","page":"Introduction","title":"MPSDynamics.spinbosonmpo","text":"spinbosonmpo(ω0, Δ, d, N, chainparams; rwa=false, tree=false)\n\nGenerate MPO for a spin-1/2 coupled to a chain of harmonic oscillators, defined by the Hamiltonian\n\nH = fracω_02σ_z + Δσ_x + c_0σ_x(b_0^dagger+b_0) + sum_i=0^N-1 t_i (b_i+1^dagger b_i +hc) + sum_i=0^N ϵ_ib_i^dagger b_i.\n\nThe spin is on site 1 of the MPS and the bath modes are to the right.\n\nThis Hamiltonain is unitarily equivalent (before the truncation to N sites) to the spin-boson Hamiltonian defined by\n\nH =  fracω_02σ_z + Δσ_x + σ_xint_0^ dωsqrtfracJ(ω)π(b_ω^dagger+b_ω) + int_0^ dω ωb_ω^dagger b_ω.\n\nThe chain parameters, supplied by chainparams=ϵ_0ϵ_1t_0t_1c_0, can be chosen to represent any arbitrary spectral density J(ω) at any temperature.\n\nThe rotating wave approximation can be made by setting rwa=true.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.svdmps-Tuple{Any}","page":"Introduction","title":"MPSDynamics.svdmps","text":"svdmps(A)\n\nFor a right normalised mps A compute the full svd spectrum for a bipartition at every bond.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.svdtrunc-Tuple{Any}","page":"Introduction","title":"MPSDynamics.svdtrunc","text":"U, S, Vd = svdtrunc(A; truncdim = max(size(A)...), truncerr = 0.)\n\nPerform a truncated SVD, with maximum number of singular values to keep equal to truncdim or truncating any singular values smaller than truncerr. If both options are provided, the smallest number of singular values will be kept. Unlike the SVD in Julia, this returns matrix U, a diagonal matrix (not a vector) S, and Vt such that A ≈ U * S * Vt\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.xyzmpo-Tuple{Int64}","page":"Introduction","title":"MPSDynamics.xyzmpo","text":"xyzmpo(N::Int; Jx=1.0, Jy=Jx, Jz=Jx, hx=0., hz=0.)\n\nReturn the MPO representation of the N-spin XYZ Hamiltonian with external field vech=(h_x 0 h_z).\n\nH = sum_n=1^N-1 -J_x σ_x^n σ_x^n+1 - J_y σ_y^n σ_y^n+1 - J_z σ_z^n σ_z^n+1 + sum_n=1^N(- h_x σ_x^n - h_z σ_z^n).\n\n\n\n\n\n","category":"method"}]
}
