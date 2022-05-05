struct Dual
    v
    δ
end

Base.:+(a::Dual, b::Dual) = Dual(a.v + b.v, a.δ + b.δ)
Base.:*(a::Dual, b::Dual) = Dual(a.v * b.v, a.v * b.δ + a.δ * b.v)
Base.log(a::Dual) = Dual(log(a.v), a.δ / a.v)
function Base.max(a::Dual, b::Dual)
    v = max(a.v, b.v)
    δ = a.v > b.v ? a.δ : a.v < b.v ? b.v : b.δ : NaN
    return Dual(v, δ)
end
function Base.max(a::Dual, b::Int)
    v = max(a.v, b)
    δ = a.v > b ? a.δ : a.v < b ? 0 : NaN
    return Dual(v, δ)
end

b = Dual(2, 0)
a = Dual(3, 1)

# ln(ab + max(a, 2))
c1 = b * a
c2 = max(a, 2)
c3 = c1 + c2
c4 = log(c3)

println(c4.v, " ", c4.δ)
