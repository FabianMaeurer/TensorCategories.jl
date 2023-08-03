#=-------------------------------------------------
	Author: Liam Rogel
-------------------------------------------------=#

#Calculates the first n quantum numbers where [2]=x
function quantum(x,n)
	eins=one(parent(x))
	q=[eins,x]
	for i in 2:n
		push!(q,q[2]*q[i]-q[i-1])
	end
	return q
end

#Calculates Quantum factorial [n]!
function fac(q,n)
       a=q[1]
       for i in 2:n 
       a*=q[i]
       end
       return a
end

#Calculates the Theta-net
function Net(q,m,n,p)
       	a=(-q[1])^(m+n+p)
	a*=fac(q,m+n+p+1)*fac(q,n)*fac(q,m)*fac(q,p)//(fac(q,n+m)*fac(q,n+p)*fac(q,m+p))
       	return a
end

#Calculates the Tetrahedron net following Kauffman-Lins
function Tet(q,A,B,C,D,E,F)
       a1=div(A+D+E,2)
       a2=div(B+C+E,2)
       a3=div(A+B+F,2)
       a4=div(C+D+F,2)
       b1=div(B+D+E+F,2)
       b2=div(A+C+E+F,2)
       b3=div(A+B+C+D,2)
       m=max(a1,a2,a3,a4)
       M=min(b1,b2,b3)
       res=q[1]
	res*=fac(q,b1-a1)*fac(q,b1-a2)*fac(q,b1-a3)*fac(q,b1-a4)*fac(q,b2-a1)*fac(q,b2-a2)*fac(q,b2-a3)*fac(q,b2-a4)*fac(q,b3-a1)*fac(q,b3-a2)*fac(q,b3-a3)*fac(q,b3-a4)//(fac(q,A)*fac(q,B)*fac(q,C)*fac(q,D)*fac(q,E)*fac(q,F))
       	sum=zero(parent(q[1]))
       	for i in m:M
		sum+=(-1)^i*fac(q,i+1)//(fac(q,i-a1)*fac(q,i-a2)*fac(q,i-a3)*fac(q,i-a4)*fac(q,b1-i)*fac(q,b2-i)*fac(q,b3-i))
       	end
       	res*=sum
       	return res
end

#Gives the 6j-Symbols for TL-Algebra. Also called Recoupling in KL
function SixJCategory(q,a,b,c,d,i,j)
       res=one(parent(q[1]))
       m1=div(a+d-i,2)
       n1=div(a-d+i,2)
       p1=div(d-a+i,2)
       m2=div(b+c-i,2)
       n2=div(c+i-b,2)
       p2=div(i+b-c,2)
	res*=Tet(q,a,b,c,d,i,j)*(-1)^i*q[i+1]//(Net(q,m1,n1,p1)*Net(q,m2,n2,p2))
       return res
end

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

#Creates the Wess-Zumino-Witten-Verlinde fusion Category
#It has fusion rules of the Dynking diagram of type A_n; 
#I.e. X_0 is the monoidal unit; X_1 generates the whole category; X_1\otimes X_i \simeq X_{i-1}\oplus X_{i+1}; exacept for X_0 and X_1\otimes X_{n-1}\simeq X_{n-2}
function I2subcategory(m)
	n=div(m,2)
	F,z=cyclotomic_field(2*m)
	q=quantum(z+z^-1,2*m-2)

	M=fusionmultmatrix(m)
	C=SixJCategory(F,[join(["Bs",repeat("ts",i)]) for i in 0:n-1])
	set_tensor_product!(C,M)

	#Now set associators using SixJCategory symbols above
	#The associator of (XY)Z-W over J and I cooresponds to SixJCategory(Y,X,W,Z,J,I)

	for lx in 1:n-1, ly in 1:n-1, lz in 1:n-1, lw in 1:n-1
		li=intersect(m_admissible(m-1,2*lx,2*ly),m_admissible(m-1,2*lz,2*lw))
		lj=intersect(m_admissible(m-1,2*lx,2*lw),m_admissible(m-1,2*ly,2*lz))
		gr=length(li)
		C.ass[lx+1,ly+1,lz+1,lw+1]=matrix(F,gr,gr,[SixJCategory(q,2*ly,2*lx,2*lw,2*lz,j,i) for i in li, j in lj])	
	end
	
	A=[0 for s in simples(C)]::Vector{Int64}
	A[1]=1
	set_one!(C,A)
	set_spherical!(C, [F(1) for s ∈ simples(C)])
	set_name!(C, "Fusion subcategory of I₂($m)")
	return C
end

function Verlinde(m)
	 #construct a Verlinde type categorie
	 K=QQBar
	 z=root_of_unity(K,2*m+2)
	q=quantum(z+z^-1,2*m+2)
	#println(q)

	M=verlindefusionmultmatrix(m)
	C=SixJCategory(K,["X$i" for i in 0:m-1])
	set_tensor_product!(C,M)

	for lx in 0:m-1, ly in 0:m-1, lz in 0:m-1, lw in 0:m-1
		#println("$lx, $ly, $lz, $lw")
		li=intersect(ver_admissible(m,lx,ly),ver_admissible(m,lw,lz))
		lj=intersect(ver_admissible(m,lx,lw),ver_admissible(m,ly,lz))
		gr=length(li)
		C.ass[lx+1,ly+1,lz+1,lw+1]=matrix(K,gr,gr,[SixJCategory(q,ly,lx,lw,lz,j,i) for i in li, j in lj])	
	end
	set_one!(C, [mod(i,m) == 1 ? 1 : 0 for i ∈ 1:m])
	TensorCategories.set_name!(C, "Verlinde Category $m")
	return C
end