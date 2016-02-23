

# additional methods for fixed size arrays

# make an immutable vector one longer
macro fsa_extend_tuple(vec, N, X)
	cexpr = Expr(:tuple)
	for i = 1:N
		push!(cexpr.args, :($vec[$i]))
	end
	push!(cexpr.args, :($X))
	return esc(Expr(:call, Vec, cexpr))		

end

# macro to extend the size of an immutable matrix by 1
macro fsa_projective(Tmat, N, M, T)

	mat_expr = Expr(:tuple)

	# add a zero to each column
	for i = 1:M
		vec_expr = Expr(:tuple)
		for j = 1:N
			push!(vec_expr.args, :($Tmat[$j, $i]))
		end
		push!(vec_expr.args, :($T(0)))
	    push!(mat_expr.args, vec_expr)
	end

	# add an extra column
	vec_expr = Expr(:tuple)
	for j = 1:N
		push!(vec_expr.args, :($T(0)))
	end
	push!(vec_expr.args, :($T(1)))
	push!(mat_expr.args, vec_expr)


	# combine all the tuples
	return esc(Expr(:call, :Mat, mat_expr))
end

# macro to test whihch elements of a fixed array are null
macro fsa_isnan(Tmat, N, M)

	mat_expr = Expr(:tuple)

	# add a zero to each column
	for i = 1:M
		vec_expr = Expr(:tuple)
		for j = 1:N
			push!(vec_expr.args, :(isnan($Tmat[$j, $i])))
		end
	    push!(mat_expr.args, vec_expr)
	end

	# combine all the tuples
	return esc(Expr(:call, :Mat, mat_expr))
end

macro fsa_isnan_vec(V, N)
	vec_expr = Expr(:tuple)
	for i = 1:N
		push!(vec_expr.args, :(isnan($V[$i])))
	end
	return esc(Expr(:call, :Vec, vec_expr))

end
