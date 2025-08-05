#=-------------------------------------------------
	Author: Liam Rogel
-------------------------------------------------=#

# #Calculates the first n quantum numbers where [2]=x
# function quantum(x,n)
	# 	eins=one(parent(x))
	# 	q=[eins,x]
	# 	for i in 2:n
		# 		push!(q,q[2]*q[i]-q[i-1])
	# 	end
	# 	return q
# end

# #Calculates Quantum factorial [n]!
# function fac(q,n)
       #        a=q[1]
       #        for i in 2:n 
       #        a*=q[i]
       #        end
       #        return a
# end

# #Calculates the Theta-net
# function Net(q,m,n,p)
       	#        	a=(-q[1])^(m+n+p)
	# 	a*=fac(q,m+n+p+1)*fac(q,n)*fac(q,m)*fac(q,p)//(fac(q,n+m)*fac(q,n+p)*fac(q,m+p))
       	#        	return a
# end

# #Calculates the Tetrahedron net following Kauffman-Lins
# function Tet(q,A,B,C,D,E,F)
	# 	a1=div(A+D+E,2)
	# 	a2=div(B+C+E,2)
	# 	a3=div(A+B+F,2)
	# 	a4=div(C+D+F,2)
	# 	b1=div(B+D+E+F,2)
	# 	b2=div(A+C+E+F,2)
	# 	b3=div(A+B+C+D,2)
	# 	m=max(a1,a2,a3,a4)
	# 	M=min(b1,b2,b3)
	# 	res=q[1]
	# 	res*=fac(q,b1-a1)*fac(q,b1-a2)*fac(q,b1-a3)*fac(q,b1-a4)*fac(q,b2-a1)*fac(q,b2-a2)*fac(q,b2-a3)*fac(q,b2-a4)*fac(q,b3-a1)*fac(q,b3-a2)*fac(q,b3-a3)*fac(q,b3-a4)//(fac(q,A)*fac(q,B)*fac(q,C)*fac(q,D)*fac(q,E)*fac(q,F))
	# 	sum=zero(parent(q[1]))
	# 	for i in m:M
	# 	sum+=(-1)^i*fac(q,i+1)//(fac(q,i-a1)*fac(q,i-a2)*fac(q,i-a3)*fac(q,i-a4)*fac(q,b1-i)*fac(q,b2-i)*fac(q,b3-i))
	# 	end
	# 	res*=sum
	# 	return res
# end

# #Gives the 6j-Symbols for TL-Algebra. Also called Recoupling in KL
# function tl_six_j_symbol(q,a,b,c,d,i,j)
	# 	res=one(parent(q[1]))
	# 	m1=div(a+d-i,2)
	# 	n1=div(a-d+i,2)
	# 	p1=div(d-a+i,2)
	# 	m2=div(b+c-i,2)
	# 	n2=div(c+i-b,2)
	# 	p2=div(i+b-c,2)
	# 	res*=Tet(q,a,b,c,d,i,j)*(-1)^i*q[i+1]//(Net(q,m1,n1,p1)*Net(q,m2,n2,p2))
	# 	return res
# end

#Verlinde admissible; returns i such that X_i occurs in X_a\otimes X_b in the verlinde category
#If a or b are greater or equal than m we should get an empty set; This functions does not do that
function ver_admissible(m,a,b)
	return [i for i in abs(a-b):2:min(a+b,2*m-a-b-2)]
end


#Verlinde fusion
function verlindefusionmultmatrix(m)
	M=zeros(Int,m,m,m)
	for i in 1:m, j in 1:m
		A=[c for c in abs(i-j)+1:2:min(2*m-i-j+2,i+j)]
		for k in A
			M[i,j,k]=1
		end
	end
	return M
end

function lambda(qq, a, b, c)
	return (-1)^(div(a+b-c, 2))*qq^(div(a*(a+2)+b*(b+2)-c*(c+2), 2))
end

#Creates the Wess-Zumino-Witten-Verlinde fusion Category
#It has fusion rules of the Dynking diagram of type A_n; 
#I.e. X_0 is the monoidal unit; X_1 generates the whole category; X_1\otimes X_i \simeq X_{i-1}\oplus X_{i+1}; exacept for X_0 and X_1\otimes X_{n-1}\simeq X_{n-2}
#construct the Verlinde type categorie with m many objects,  l-th associator and k-th braiding; Need (m+1, l) coprime
#furthermore for the braiding k we use even higher roots of unity,  but also need (m+1, k) coprime
#Formulas come from "Temperley-Lieb Recoupling theory"
@doc raw""" 

	verlinde_category(K::Ring, m::Int)

Return the Verlinde category specialized at a ``m``th rooth of unity.
"""
function verlinde_category(K::Ring, m::Int, l::Int = 1, t::Int = 1)
	#K=Oscar.QQBarField()
	z = root_of_unity(K, 2*m+2)
   	q = quantum(z^l+z^-l, 2*m+2)
	
   	#println(q)

   	M = verlindefusionmultmatrix(m)
   	C = six_j_category(K, ["X$i" for i in 0:m-1])
   	set_tensor_product!(C, M)
	C.ass = Array{MatElem}(undef,m,m,m,m)

	a = (i,j,k,l) -> verlinde_six_j_symbol(K,q,i,j,k,l,m)
	set_attribute!(C, :six_j_symbol, a)

   	# for lx in 0:m-1,  ly in 0:m-1,  lz in 0:m-1,  lw in 0:m-1
	#    #println("$lx,  $ly,  $lz,  $lw")
	#    li = intersect(ver_admissible(m, lx, ly), ver_admissible(m, lw, lz))
	#    lj = intersect(ver_admissible(m, lx, lw), ver_admissible(m, ly, lz))
	#    gr = length(li)
	#    C.ass[lx+1, ly+1, lz+1, lw+1] = matrix(K, gr, gr, [tl_six_j_symbol(q, ly, lx, lw, lz, j, i) for i in li,  j in lj])	
   	# end


	#zz = root_of_unity(K, 4*m+4)
	has_braiding = K == QQBarField() || typeof(K) == CalciumField || is_square(z)
	if has_braiding
		zz = sqrt(z)
		braid = Array{MatElem, 3}(undef,  m, m, m) #slight mistake,  need to check formulas again
		set_braiding!(C, braid)

		b = (i,j,k) -> verlinde_r_symbol(K,i,j,k,m,zz,t)
		set_attribute!(C, :r_symbol, b)
		# for a in 0:m-1
		# 	for b in 0:m-1
		# 		M=ver_admissible(m, a, b)
		# 		#println(a,b,M)
		# 		for i in 0:m-1
		# 			if i in M
		# 				#@show braid[a+1, b+1, i+1] = matrix(K,  1,  1,  [(-1)^(i)*q[i+1]*lambda(zz^(k+m+1), a, b, i)//Net(q, div(a+b-i, 2), div(a-b+i, 2), div(-a+b+i, 2))])
		# 				braid[a+1, b+1, i+1] = matrix(K,  1,  1,  [lambda(zz^(k),a,b,i)])
		# 				#braid[a+1, b+1, i+1] = matrix(K,  1,  1,  [(-1)^(div(a+b+i,2))*zz^(div(i*(i+2)-a*(a+2)-b*(b+2),2))])
		# 			else 
		# 				braid[a+1, b+1, i+1] = matrix(K,  0,  0,  [])
		# 			end
		# 		end
		# 	end
		# end

		# 
	end



   	set_one!(C,  [mod(i, m) == 1 ? 1 : 0 for i ∈ 1:m])
	if has_braiding 
	   	set_name!(C,  "Verlinde Category $m with associator $l and braiding $t")
	else
		set_name!(C,  "Verlinde Category $m with associator $l")
	end
   	return C
end

function verlinde_six_j_symbol(K::Ring, q, i::Int, j::Int, k::Int, l::Int, m::Int)
	lx,ly,lz,lw = i-1,j-1,k-1,l-1
	li = intersect(ver_admissible(m, lx, ly), ver_admissible(m, lw, lz))
	lj = intersect(ver_admissible(m, lx, lw), ver_admissible(m, ly, lz))
	gr = length(li)
	
	matrix(K, gr, gr, [tl_six_j_symbol(q, ly, lx, lw, lz, j, i) for i in li,  j in lj])	
end

function verlinde_r_symbol(K::Ring, i::Int, j::Int, k::Int, m::Int, zz::RingElem, t::Int = 1)
	a,b = i-1,j-1
	M = ver_admissible(m, a, b)
	#println(a,b,M)
	if k-1 in M
		#@show braid[a+1, b+1, i+1] = matrix(K,  1,  1,  [(-1)^(i)*q[i+1]*lambda(zz^(k+m+1), a, b, i)//Net(q, div(a+b-i, 2), div(a-b+i, 2), div(-a+b+i, 2))])
		return matrix(K,  1,  1,  [lambda(zz^(t),a,b,k-1)])
		#braid[a+1, b+1, i+1] = matrix(K,  1,  1,  [(-1)^(div(a+b+i,2))*zz^(div(i*(i+2)-a*(a+2)-b*(b+2),2))])
	else 
		return matrix(K,  0,  0,  [])
	end
end

@doc raw""" 

	verlinde_category(m::Int)

Return the Verlinde category specialized at a ``m``th rooth of unity.
"""
verlinde_category(m::Int, l::Int = 1, k::Int = 1) = verlinde_category(cyclotomic_field(4*m+4)[1],m, l, k)