var documenterSearchIndex = {"docs":
[{"location":"#MPSDynamics.jl-Documentation","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.jl Documentation","text":"","category":"section"},{"location":"","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.jl Documentation","text":"Modules = [MPSDynamics]","category":"page"},{"location":"#MPSDynamics.MPOtoVector-Tuple{ITensors.MPO}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.MPOtoVector","text":"MPOtoVector(mpo::MPO)\n\nConvert an ITensors chain MPO into a form compatible with MPSDynamics\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.addchild!-Tuple{MPSDynamics.Tree,Int64}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.addchild!","text":"addchild!(tree::Tree, id::Int)\n\nAdd child to node id of tree.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.addchildren!-Tuple{MPSDynamics.Tree,Int64,Int64}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.addchildren!","text":"addchildren!(tree::Tree, id::Int, n::Int)\n\nAdd n children to node id of tree.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.chaincoeffs_ohmic-Tuple{Any,Any,Any}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.chaincoeffs_ohmic","text":"chaincoeffs_ohmic(nummodes, α, s, beta=\"inf\"; wc=1, soft=false)\n\nGenerate chain coefficients for an Harmonic bath at zero temperature with a power law spectral density given by: \n\nsoft cutoff: J(ω) = 2παω_c (fracωω_c)^s exp(-ωω_c) \n\nhard cutoff: J(ω) = 2παω_c (fracωω_c)^s θ(ω-ω_c)\n\nThe Hamiltonian is given by:\n\nH = fracω_02σ_z + Δσ_x + σ_xsum_kg_k(b_k^dagger+b_k) + sum_kω_kb_k^dagger b_k\n\nAnd the spectral density is defined by:\n\nJ(ω)  πsum_kg_k^2δ(ω-ω_k)\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.chainmps-Tuple{Int64,Array{Int64,1},Int64}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.chainmps","text":"chainmps(N::Int, sites::Vector{Int}, numex::Int)\n\nGenerate an MPS with numex excitations of an equal super-position over sites\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.chainmps-Tuple{Int64,Int64,Int64}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.chainmps","text":"chainmps(N::Int, site::Int, numex::Int)\n\nGenerate an MPS with numex excitations on site\n\nThe returned MPS will have bond-dimensions and physical dimensions numex+1\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.chainprop-Tuple{Any,Any}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.chainprop","text":"chainprop(t, cparams...)\n\nPropagate an excitation placed initially on the first site of a tight-binding chain with parameters given by cparams for a time t and return occupation expectation for each site.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.dynamap-NTuple{4,Any}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.dynamap","text":"dynamap(ps1,ps2,ps3,ps4)\n\nCalulate complete dynamical map to time step at which ps1, ps2, ps3 and ps4 are specified.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.electron2kmps","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.electron2kmps","text":"electronkmps(N::Int, k::Vector{Int}, spin=:Up, chainparams=[fill(1.0,N), fill(1.0,N-1)])\n\nGenerate an MPS with 2 electrons in k-states k1 and k2.\n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.electronkmps","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.electronkmps","text":"electronkmps(N::Int, k::Int, spin=:Up, chainparams=[fill(1.0,N), fill(1.0,N-1)])\n\nGenerate an MPS for an electron with momentum k.\n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.elementmpo-Tuple{Any,Vararg{Any,N} where N}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.elementmpo","text":"elementmpo(M, el...)\n\nReturn the element of the MPO M for the set of physical states el...\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.elementmps-Tuple{Any,Vararg{Any,N} where N}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.elementmps","text":"elementmps(A, el...)\n\nReturn the element of the MPS A for the set of physical states el...\n\nExamples\n\njulia> A = chainmps(6, [2,4], 1);\n\njulia> elementmps(A, 1, 2, 1, 1, 1, 1)\n0.7071067811865475\n\njulia> elementmps(A, 1, 1, 1, 2, 1, 1)\n0.7071067811865475\n\njulia> elementmps(A, 1, 2, 1, 2, 1, 1)\n0.0\n\njulia> elementmps(A, 1, 1, 1, 1, 1, 1)\n0.0\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.entanglemententropy-Tuple{Any}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.entanglemententropy","text":"entanglemententropy(A)\n\nFor a list of tensors A representing a right orthonormalized MPS, compute the entanglement entropy for a bipartite cut for every bond.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.findchainlength-Tuple{Any,Any}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.findchainlength","text":"findchainlength(T, cparams...; eps=10^-6)\n\nEstimate length of chain required for a particular set of chain parameters by calulating how long an excitation on the first site takes to reach the end. The chain length is given as the length required for the excitation to have just reached the last site after time T.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.findchild-Tuple{MPSDynamics.TreeNode,Int64}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.findchild","text":"findchild(node::TreeNode, id::Int)\n\nReturn integer corresponding to the which number child site id is of node.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.measure-Tuple{Any,OneSiteObservable}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.measure","text":"measure(A, O; kwargs...)\n\nmeasure observable O on mps state A\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.measure1siteoperator-Tuple{Array{T,1} where T,Any,Array{Int64,1}}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.measure1siteoperator","text":"measure1siteoperator(A, O)\n\nFor a list of tensors A representing a right orthonormalized MPS, compute the local expectation value of a one-site operator O for every site or just one if i is specified.\n\nFor calculating operators on single sites this will be more efficient if the site is on the left of the mps.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.measure2siteoperator-Tuple{Array{T,1} where T,Any,Any,Int64,Int64}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.measure2siteoperator","text":" measure2siteoperator(A::Vector, M1, M2, j1, j2)\n\nCaculate expectation of M1*M2 where M1 acts on site j1 and M2 acts on site j2, assumes A is right normalised.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.modemps","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.modemps","text":"modemps(N::Int, k::Vector{Int}, numex::Int, chainparams=[fill(1.0,N), fill(1.0,N-1)])\n\nGenerate an MPS with numex excitations of an equal superposition of modes k of a bosonic tight-binding chain.\n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.modemps-2","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.modemps","text":"modemps(N::Int, k::Int, numex::Int, chainparams=[fill(1.0,N), fill(1.0,N-1)])\n\nGenerate an MPS with numex excitations of mode k of a bosonic tight-binding chain. \n\nchainparams takes the form [e::Vector, t::Vector] where e are the on-site energies and t are the hoppping parameters.\n\nThe returned MPS will have bond-dimensions and physical dimensions numex+1\n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.mpsembed!-Tuple{Array{T,1} where T,Int64}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.mpsembed!","text":"mpsembed(A::Vector, Dmax::Int)\n\nEmbed MPS A in manifold of max bond-dimension Dmax\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.mpsleftnorm!","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.mpsleftnorm!","text":"mpsleftnorm!(A::Vector, jq::Int=length(A))\n\nLeft orthoganalise MPS A up to site jq.\n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.mpsmixednorm!-Tuple{Array{T,1} where T,Int64}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.mpsmixednorm!","text":"mpsmixednorm!(A::Vector, OC::Int)\n\nPut MPS A into mixed canonical form with orthogonality centre on site OC.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.mpsmixednorm!-Tuple{MPSDynamics.TreeNetwork,Int64}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.mpsmixednorm!","text":"mpsmixednorm!(A::TreeNetwork, id::Int)\n\nNormalise tree-MPS A such that orthogonality centre is on site id.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.mpsmoveoc!-Tuple{MPSDynamics.TreeNetwork,Int64}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.mpsmoveoc!","text":"mpsmoveoc!(A::TreeNetwork, id::Int)\n\nMove the orthogonality centre of right normalised tree-MPS A to site id.\n\nThis function will be more efficient than using mpsmixednorm! if the tree-MPS is already right-normalised.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.mpsrightnorm!","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.mpsrightnorm!","text":"mpsrightnorm!(A::Vector, jq::Int=1)\n\nRight orthoganalise MPS A up to site jq.\n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.mpsrightnorm!-Tuple{MPSDynamics.TreeNetwork}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.mpsrightnorm!","text":"mpsrightnorm!(A::TreeNetwork)\n\nWhen applied to a tree-MPS, right normalise towards head-node.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.mpsshiftoc!-Tuple{MPSDynamics.TreeNetwork,Int64}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.mpsshiftoc!","text":"mpsshiftoc!(A::TreeNetwork, newhd::Int)\n\nShift the orthogonality centre by one site, setting new head-node newhd.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.normmps-Tuple{Array{T,1} where T}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.normmps","text":"normmps(A::Vector; mpsorthog=:None)\n\nCalculate norm of MPS A.\n\nSetting mpsorthog=:Right/:Left will calculate the norm assuming right/left canonical form. Setting mpsorthog=OC::Int will cause the norm to be calculated assuming the orthoganility center is on site OC. If mpsorthog is :None the norm will be calculated as an MPS-MPS product.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.normmps-Tuple{MPSDynamics.TreeNetwork}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.normmps","text":"normmps(net::TreeNetwork; mpsorthog=:None)\n\nWhen applied to a tree-MPS mpsorthog=:Left is not defined.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.orthcentersmps-Tuple{Array{T,1} where T}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.orthcentersmps","text":"orthcentersmps(A)\n\nCompute the orthoganality centres of MPS A.\n\nReturn value is a list in which each element is the corresponding site tensor of A with the orthoganility centre on that site. Assumes A is right normalised.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.physdims-Tuple{Array{T,1} where T}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.physdims","text":"physdims(M)\n\nReturn the physical dimensions of an MPS or MPO M.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.productstatemps","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.productstatemps","text":"productstatemps(N::Int, d::Int, Dmax=1; state=:Vacuum, mpsorthog=:Right)\n\nReturn an N-site MPS with all local Hilbert space dimensions given by d. \n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.productstatemps-2","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.productstatemps","text":"productstatemps(physdims::Dims, Dmax=1; state=:Vacuum, mpsorthog=:Right)\n\nReturn an MPS representing a product state with local Hilbert space dimensions given by physdims.\n\nBy default all bond-dimensions will be 1 since the state is a product state. However, to embed the product state in a manifold of greater bond-dimension, Dmax can be set accordingly.\n\nThe indvidual states of the MPS sites can be provdided by setting state to a list of column vectors. Setting state=:Vacuum will produce an MPS in the vacuum state (where the state of each site is represented by a column vector with a 1 in the first row and zeros elsewhere). Setting state=:FullOccupy will produce an MPS in which each site is fully occupied (ie. a column vector with a 1 in the last row and zeros elsewhere).\n\nThe argument mpsorthog can be used to set the gauge of the resulting MPS.\n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.randisometry-Tuple{Type,Int64,Int64}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.randisometry","text":"randisometry([T=Float64], dims...)\n\nConstruct a random isometry\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.randmps","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.randmps","text":"randmps(N::Int, d::Int, Dmax::Int, T=Float64)\n\nConstruct a random, N-site, right-normalised MPS with all local Hilbert space dimensions given by d.\n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.randmps-2","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.randmps","text":"randmps(tree::Tree, physdims, Dmax::Int, T::Type{<:Number} = Float64)\n\nConstruct a random, right-normalised, tree-MPS, with structure given by tree and max bond-dimension given by Dmax.\n\nThe local Hilbert space dimensions are specified by physdims which can either be of type Dims{length(tree)}, specifying the dimension of each site, or of type Int, in which case the same local dimension is used for every site.\n\n\n\n\n\n","category":"function"},{"location":"#MPSDynamics.randmps-Union{Tuple{N}, Tuple{Tuple{Vararg{Int64,N}},Int64}, Tuple{Tuple{Vararg{Int64,N}},Int64,Type{#s131} where #s131<:Number}} where N","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.randmps","text":"randmps(physdims::Dims{N}, Dmax::Int, T::Type{<:Number} = Float64) where {N}\n\nConstruct a random, right-normalised MPS with local Hilbert space dimensions given by physdims and max bond-dimension given by Dmax. \n\nT specifies the element type, eg. use T=ComplexF64 for a complex valued MPS.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.randtree-Tuple{Int64,Int64}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.randtree","text":"randtree(numnodes::Int, maxdegree::Int)\n\nConstruct a random tree with nummodes modes and max degree maxdegree.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.rmsd-Tuple{Any,Any}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.rmsd","text":"rmsd(dat1::Vector{Float64}, dat2::Vector{Float64})\n\nCalculate the root mean squared difference between two measurements of an observable over the same time period.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.spinbosonmpo-NTuple{5,Any}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.spinbosonmpo","text":"spinbosonmpo(ω0, Δ, d, N, chainparams; rwa=false, tree=false)\n\nGenerate MPO for a spin-1/2 coupled to a chain of harmonic oscillators, defined by the Hamiltonian\n\nH = fracω_02σ_z + Δσ_x + c_0σ_x(b_k^dagger+b_k) + sum_i=0^N-1 t_i (b_i+1^dagger b_i +hc) + sum_i=1^N-1 ϵ_ib_i^dagger b_i.\n\nThis Hamiltonain is unitarily equivalent to the spin-boson Hamiltonian defined by\n\nH =  fracω_02σ_z + Δσ_x + σ_xint_0^ dωsqrtJ(ω)(b_ω^dagger+b_ω) + int_0^ωb_ω^dagger b_ω.\n\nThe chain parameters, supplied by chainparams=[[ϵ_0,ϵ_1,...],[t_0,t_1,...],c_0], can be chosen to represent any arbitrary spectral density J(ω) at any temperature.\n\nThe rotating wave approximation can be made by setting rwa=true.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.svdmps-Tuple{Any}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.svdmps","text":"svdmps(A)\n\nFor a right normalised mps A compute the full svd spectrum for a bipartition at every bond.\n\n\n\n\n\n","category":"method"},{"location":"#MPSDynamics.svdtrunc-Tuple{Any}","page":"MPSDynamics.jl Documentation","title":"MPSDynamics.svdtrunc","text":"U, S, Vd = svdtrunc(A; truncdim = max(size(A)...), truncerr = 0.)\n\nPerform a truncated SVD, with maximum number of singular values to keep equal to truncdim or truncating any singular values smaller than truncerr. If both options are provided, the smallest number of singular values will be kept. Unlike the SVD in Julia, this returns matrix U, a diagonal matrix (not a vector) S, and Vt such that A ≈ U * S * Vt\n\n\n\n\n\n","category":"method"}]
}
