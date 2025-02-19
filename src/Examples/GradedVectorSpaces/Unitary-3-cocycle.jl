#=------------------------------------------------
	Author: Liam Rogel
------------------------------------------------=#



function exponent_of_homology(G,k)
	torsion_coefficients=GAP.Globals.GroupHomology(G.X,k)
	torsion_coefficients=Vector{ZZRingElem}(GAP.gap_to_julia(torsion_coefficients))
	return length(torsion_coefficients) != 0 ? lcm(torsion_coefficients) : 1
end

function unitary_cocycle(G::Group, K::Field, k::Int, i::Int = 2) 
	GAP.Packages.install("HAP")
	GAP.Packages.load("HAP") 

	n = exponent_of_homology(G,k)
	if n == 1
		return trivial_cocycle(G,K,k)
	end

	ρ = root_of_unity(K,Int(n))
	
	x = gens(symmetric_group(Int(n)))[1]  #We create a cyclic group of order n
	B,g = sub(x)
	D = Dict{NTuple{k,elem_type(G)},elem_type(parent(ρ))}()

	A = GAP.Globals.TrivialGModuleAsGOuterGroup(G.X,B.X) #This is the cyclic group encoded as a trivial G-module;We need this weird G.X notation to get GAP elements; Need GAP.Globals because tehre is no wrapper for HAP
	R = GAP.Globals.ResolutionFiniteGroup(G.X,k+1)
	C = GAP.Globals.HomToGModule(R,A)
	CH = GAP.Globals.CohomologyModule(C,k)
	classes = GAP.Globals.Elements(GAP.Globals.ActedGroup(CH)) #This is the list of cohomoogy classes
	if length(classes)<i
		println("We have only $(length(classes)) many classes of cocycles")
		return D
	end
	get_representativeCocycle = GAP.evalstr("x -> x!.representativeCocycle")
    c = get_representativeCocycle(CH)(classes[i])
	f = GAP.Globals.Mapping(c)
	Elts = GAP.Globals.Elements(B.X)
	
	for g in Base.product([elements(G) for _ ∈ 1:k]...)
				#print(f(g.X,h.X,k.X))
				exponent = GAP.Globals.Position(Elts,f([h.X for h ∈ g]...))-1
				push!(D, reverse(g) => ρ^exponent)
				end
	return Cocycle(G,D)
end

""" 

	twisted_graded_vector_spaces(G::Group, i::Int)

Construct the category of twisted graded vectorspaces with the i-th 3-cocycle.
"""
function twisted_graded_vector_spaces(K::Field, G::Group, i::Int = 2, j::Int = 1) #Inputs finite group G and the number of Cocycle we want in the twisted_graded_vector_spaces
	GAP.Packages.install("HAP")
	GAP.Packages.load("HAP") 
    
    #K,ρ=cyclotomic_field(6*Int(n),"ρ") #for some reason cyclotomic_field wants Int64...
	ξ = unitary_cocycle(G,K,3,i) #One must know whetever we have i many classes
	#braid = unitary_cocycle(G,K,2,j)
    return graded_vector_spaces(K,G,ξ)
    #pentagon_axiom(VecGtw)
end


function twisted_graded_vector_spaces(G::Group, i::Int = 2, j::Int = 1)
	twisted_graded_vector_spaces(G, QQBarField(), i, j)
end