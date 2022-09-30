using Oscar
import GAP
GAP.Packages.install("HAP")
GAP.Packages.load("HAP") 

function exponent_of_third_homology(G)
	torsion_coefficients=GAP.Globals.GroupHomology(G.X,3)
	torsion_coefficients=Vector{fmpz}(GAP.gap_to_julia(torsion_coefficients))
	return lcm(torsion_coefficients)
end

function unitary_3_cocycle(G,n,ρ,i) 
	x=gens(symmetric_group(Int(n)))[1]  #We create a cyclic group of order n
	B,g=subgroup(x)
	D = Dict{NTuple{3,elem_type(G)},elem_type(parent(ρ))}()

	A=GAP.Globals.TrivialGModuleAsGOuterGroup(G.X,B.X) #This is the cyclic group encoded as a trivial G-module;We need this weird G.X notation to get GAP elements; Need GAP.Globals because tehre is no wrapper for HAP
	R=GAP.Globals.ResolutionFiniteGroup(G.X,4)
	C=GAP.Globals.HomToGModule(R,A)
	CH=GAP.Globals.CohomologyModule(C,3)
	classes=GAP.Globals.Elements(GAP.Globals.ActedGroup(CH)) #This is the list of cohomoogy classes
	if length(classes)<i
		println("We have only $(length(classes)) many classes of cocycles")
		return D
	end
	get_representativeCocycle = GAP.evalstr("x -> x!.representativeCocycle")
    c = get_representativeCocycle(CH)(classes[i])
	f=GAP.Globals.Mapping(c)
	Elts=GAP.Globals.Elements(B.X)
	
	for g in elements(G)
		for h in elements(G)
			for k in elements(G)
				#print(f(g.X,h.X,k.X))
				exponent=GAP.Globals.Position(Elts,f(g.X,h.X,k.X))-1
				push!(D,(g,h,k)=>ρ^exponent)
			end
		end
	end
	return Cocycle(G,D)
end

function twisted_graded_vector_spaces(G,i) #Inputs finite group G and the number of Cocycle we want in the twisted_graded_vector_spaces
    n=exponent_of_third_homology(G)
    K,ρ=CyclotomicField(Int(n),"ρ") #for some reason cyclotomic_field wants Int64...
    ξ=unitary_3_cocycle(G,n,ρ,i) #One must know whetever we have i many classes
    return GradedVectorSpaces(K,G,ξ)
    #pentagon_axiom(VecGtw)
end

