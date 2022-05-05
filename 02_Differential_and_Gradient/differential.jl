diff_forward(f, x, h=sqrt(eps(Float64))) = (f(x + h) - f(x)) / h
diff_backward(f, x, h=sqrt(eps(Float64))) = (f(x) - f(x - h)) / h
diff_central(f, x, h=cbrt(eps(Float64))) = (f(x + h) - f(x - h)) / (2 * h)
diff_complex(f, x; h = 1e-20) = imag(f(x + h * im)) / h

f(x) = 3x^2 + 2x
f_prime(x) = 6x + 2

println("f_prime(2) = ", f_prime(2))
println("diff_forward = ", diff_forward(f, 2))
println("diff_backward = ", diff_backward(f, 2))
println("diff_central = ", diff_central(f, 2))
println("diff_complex = ", diff_complex(f, 2))
