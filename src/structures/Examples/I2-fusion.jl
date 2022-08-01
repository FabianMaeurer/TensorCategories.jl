#=-------------------------------------------------
	Examples by Liam Rogel
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
function SixJ(q,a,b,c,d,i,j)
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

#Return the m-admissible pairings to the tuple (a,b) (See 7.1)
function m_admissible(m,a,b)
	mini=abs(a-b)
	maxi=min(2*m-2-a-b,a+b)
	A=[i for i in mini:2:maxi]
	return A
end

#Returns a 2(m-1)^3 Matrix describing the monoidal product in the multifusionversion of I2(m)
function multifusionmultmatrix(m)
	#N gives multtable for "half" the objects if for example Bs\otimes Bts would exist	
	n=m-1
	N=zeros(Int,n,n,n)
	for i in 1:n, j in 1:n
		A=m_admissible(n,i-1,j-1)
		for k in A
			N[i,j,k+1]=1
		end
	end

	#Now we modify to get the whole table
	M=zeros(Int,2*n,2*n,2*n)	
	for i in 1:n, j in 1:n, k in 1:n
		if rem(i,2)==1
			M[i,j,k]=N[i,j,k]
			M[i+n,j+n,k+n]=N[i,j,k]
		else
			M[i,j+n,k]=N[i,j,k]
			M[i+n,j,k+n]=N[i,j,k]
		end
	end
	return M
end

#Returns a div(m,2)^3 Matrix describing the monoidal product in the fusion subcat of I2(m
#More precise: It is a sixtyfourth of the above matrix, zoomed in onto the objects (Bs,Bsts,Bststs,...)
function fusionmultmatrix(m)
	n=div(m,2)
	M=zeros(Int,n,n,n)
	for i in 1:n, j in 1:n
		A=m_admissible(m-1,2*(i-1),2*(j-1))
		for k in A
			M[i,j,div(k,2)+1]=1
		end
	end
	return M
end

#Creates the categorification of biggest cell in I2(m). So far only the associators work. ev/coev needs more work
function I2(m)
	n=m-1
	F,z=cyclotomic_field(2*m)
	q=quantum(z+z^-1,2*n)
	#println(q)

	M=multifusionmultmatrix(m)
	C=RingCategory(F,vcat([if rem(i,2)==1 join(["Bs",repeat("ts",div(i,2))]) else join(["B",repeat("st",div(i,2))]) end for i in 1:n],[if rem(i,2)==1 join(["Bt",repeat("st",div(i,2))]) else join(["B",repeat("ts",div(i,2))]) end for i in 1:n]))
	set_tensor_product!(C,M)

	#Now set associators using SixJ symbols above
	#Note, that Bs, Bt always have trivial associators (prove this)
	#Further note, that the order of elements is different. 
	#The associator of (XY)Z-W over J and I cooresponds to SixJ(Y,X,W,Z,J,I)

	for lx in 1:n-1, ly in 1:n-1, lz in 1:n-1, lw in 1:n-1
		li=intersect(m_admissible(n,lx,ly),m_admissible(n,lz,lw))
		lj=intersect(m_admissible(n,lx,lw),m_admissible(n,ly,lz))
		gr=length(li)
		#println(lx,ly,lz,lw,li,lj)
		if rem(lx,2)==1 && rem(ly,2)==1
	 		C.ass[lx+1,ly+1+n,lz+1,lw+1]=matrix(F,gr,gr,[SixJ(q,ly,lx,lw,lz,j,i) for i in li, j in lj])	
			C.ass[lx+1+n,ly+1,lz+1+n,lw+1+n]=matrix(F,gr,gr,[SixJ(q,ly,lx,lw,lz,j,i) for i in li, j in lj])
		elseif rem(lx,2)==1 && rem(ly,2)==0
			C.ass[lx+1,ly+1+n,lz+1+n,lw+1]=matrix(F,gr,gr,[SixJ(q,ly,lx,lw,lz,j,i) for i in li, j in lj])	
			C.ass[lx+1+n,ly+1,lz+1,lw+1+n]=matrix(F,gr,gr,[SixJ(q,ly,lx,lw,lz,j,i) for i in li, j in lj])
		elseif rem(lx,2)==0 && rem(ly,2)==1
			C.ass[lx+1,ly+1,lz+1+n,lw+1]=matrix(F,gr,gr,[SixJ(q,ly,lx,lw,lz,j,i) for i in li, j in lj])	
			C.ass[lx+1+n,ly+1+n,lz+1,lw+1+n]=matrix(F,gr,gr,[SixJ(q,ly,lx,lw,lz,j,i) for i in li, j in lj])
		else
			C.ass[lx+1,ly+1,lz+1,lw+1]=matrix(F,gr,gr,[SixJ(q,ly,lx,lw,lz,j,i) for i in li, j in lj])	
			C.ass[lx+1+n,ly+1+n,lz+1+n,lw+1+n]=matrix(F,gr,gr,[SixJ(q,ly,lx,lw,lz,j,i) for i in li, j in lj])
		end
	end
	
	return C
end

#Creates the fusion subcategory version. Objects Bs, Bsts, Bststs,...
function I2subcategory(m)
	n=div(m,2)
	F,z=cyclotomic_field(2*m)
	q=quantum(z+z^-1,2*m-2)

	M=fusionmultmatrix(m)
	C=RingCategory(F,[join(["Bs",repeat("ts",i)]) for i in 0:n-1])
	set_tensor_product!(C,M)

	#Now set associators using SixJ symbols above
	#Note, that Bs, Bt always have trivial associators (prove this)
	#Further note, that the order of elements is different. 
	#The associator of (XY)Z-W over J and I cooresponds to SixJ(Y,X,W,Z,J,I)

	for lx in 1:n-1, ly in 1:n-1, lz in 1:n-1, lw in 1:n-1
		li=intersect(m_admissible(m-1,2*lx,2*ly),m_admissible(m-1,2*lz,2*lw))
		lj=intersect(m_admissible(m-1,2*lx,2*lw),m_admissible(m-1,2*ly,2*lz))
		gr=length(li)
		C.ass[lx+1,ly+1,lz+1,lw+1]=matrix(F,gr,gr,[SixJ(q,2*ly,2*lx,2*lw,2*lz,j,i) for i in li, j in lj])	
	end
	
	A=[0 for s in simples(C)]::Vector{Int64}
	A[1]=1
	set_one!(C,A)
	set_spherical!(C, [F(1) for s ∈ simples(C)])
	
	return C
end
