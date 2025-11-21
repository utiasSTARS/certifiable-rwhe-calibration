
# Short script demonstrating the necessity of the positive orthant constraint for permutation matrices
N = 1001
a = collect(range(-1.0, stop=1.0, length=N))
b = a
res = zeros(N, N)

function residual(a, b)
    a^2 + b^2 - a - b + a*b
end

for i in 1:lastindex(a) 
    for j in 1:lastindex(b) 
        res[i, j] = residual(a[i], b[j])
    end
end


findall(abs.(res) .< 0.0001)