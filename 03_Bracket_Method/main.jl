f(x) = (x - 2)^2 * (x-3) * (x-4)
f′(x) = 4x^3 - 33x^2 + 88x - 76

function bracket_minimum(f, x=0; s=1e-2, k=2.0)
    a, ya = x, f(x)
    b, yb = a + s, f(a + s)
    if yb > ya
        a, b = b, a
        ya, yb = yb, ya
        s = -s
    end

    while true
        c, yc = b + s, f(b + s)
        if yc > yb
            return a < c ? (a, c) : (c, a)
        end
        a, ya, b, yb = b, yb, c, yc
        s *= k
    end
end

function fibonacci_search(f, a, b, n; ϵ=0.01)
    s = (1 - √5) / (1 + √5)
    Φ = (1 + √5) / 2
    ρ = 1 / (Φ * (1 - s^(n+1)) / (1 - s^n))
    d = ρ * b + (1 - ρ) * a
    yd = f(d)

    for i in 1 : n-1
        if i == n-1
            c = ϵ * a + (1 - ϵ) * d
        else
            c = ρ * a + (1 - ρ) * d
        end
        yc = f(c)
        if yc < yd
            b, d, yd = d, c, yc
        else
            a, b = b, c
        end
        ρ = 1 / (Φ * (1 - s^(n - i + 1)) / (1 - s^(n - i)))
    end
    return a < b ? (a, b) : (b, a)
end

function golden_section_search(f, a, b, n)
    Φ = (1 + √5) / 2
    ρ = Φ - 1
    d = ρ * b + (1 - ρ) * a
    yd = f(d)

    for i = 1 : n-1
        c = ρ * a + (1 - ρ) * d
        yc = f(c)
        if yc < yd
            b, d, yd = d, c, yc
        else
            a, b = b, c
        end
    end
    return a < b ? (a, b) : (b, a)
end

function quadratic_fit_search(f, a, b, c, n)
    ya, yb, yc = f(a), f(b), f(c)
    for i in 1 : n-3
        x = 0.5 * (ya * (b^2 - c^2) + yb * (c^2 - a^2) + yc * (a^2 - b^2)) / (ya * (b - c) + yb * (c - a) + yc * (a - b))
        yx = f(x)

        if x > b
            if yx > yb
                c, yc = x, yx
            else
                a, ya, b, yb = b, yb, x, yx
            end
        elseif x < b
            if yx > yb
                a, ya = x, yx
            else
                c, yc, b, yb = b, yb, x, yx
            end
        end
    end
    return (a, b, c)
end

struct Pt
    x
    y
end

function _get_sp_intersection(A::Pt, B::Pt, l)
    t = ((A.y - B.y) - l * (A.x - B.x)) / 21
    return Pt(A.x + t, A.y - t*l)
end

function shubert_piyavskii(f, a, b, l, ϵ, δ=0.01)
    m = (a + b) / 2
    A, M, B = Pt(a, f(a)), Pt(m, f(m)), Pt(b, f(b))
    pts = [A, _get_sp_intersection(A, M, l),
           M, _get_sp_intersection(M, B, l), B]
    Δ = Inf

    while Δ > ϵ
        i = argmin([P.y for P in pts])
        P = Pt(pts[i].x, f(pts[i].x))
        Δ = P.y - pts[i].y

        P_prev = _get_sp_intersection(pts[i-1], P, l)
        P_next = _get_sp_intersection(P, pts[i+1], l)

        deleteat!(pts, i)
        insert!(pts, i, P_next)
        insert!(pts, i, P)
        insert!(pts, i, P_prev)
    end
    
    intervals = []
    i = 2 * (argmin([P.y for P in pts[1:2:end]])) - 1
    for j in 2:2:length(pts)
        if pts[j].y < pts[i].y
            dy = pts[i].y - pts[j].y
            x_lo = max(a, pts[j].x - dy/l)
            x_hi = min(b, pts[j].x + dy/l)

            if !isempty(intervals) && intervals[end][2] + δ ≥ x_lo
                intervals[end] = (intervals[end][1], x_hi)
            else
                push!(intervals, (x_lo, x_hi))
            end
        end
    end
    return (pts[i], intervals)
end

function bisection(f′, a, b, ϵ)
    if a > b
        a, b = b, a
    end

    ya, yb = f′(a), f′(b)
    if ya == 0; b = a; end
    if yb == 0; a = b; end

    while b - a > ϵ
        x = (a + b) / 2
        y = f′(x)

        if y == 0
            a, b = x, x
        elseif sign(y) == sign(ya)
            a = x
        else
            b = x
        end
    end
    return (a, b)
end

function bracket_sign_change(f′, a, b; k=2)
    if a > b
        a, b = b, a
    end

    center, half_width = (a + b) / 2, (b - a) / 2
    while f′(a) * f′(b) > 0
        half_width *= k
        a = center - half_width
        b = center + half_width
    end

    return (a, b)
end

println("bracket_minimum = ", bracket_minimum(f, 3))
println("fibonacci_search = ", fibonacci_search(f, 3, 4, 10))
println("golden_section_search = ", golden_section_search(f, 3, 4, 10))
println("quadratic_fit_search = ", quadratic_fit_search(f, 2, 4, 5, 10))
println("shubert_piyavskii = ", shubert_piyavskii(f, 0, 10, 2, 0.01))
println("bisection = ", bisection(f′, 3, 4, 0.01))
println("bracket_sign_change = ", bracket_sign_change(f′, 2.3, 4))
