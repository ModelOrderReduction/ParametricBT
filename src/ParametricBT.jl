module ParametricBT

using FiniteDifferences, LinearAlgebra, MatrixEquations

function compute_taylor_coeffs(A, x0, order)
    @assert length(x0) == 1 "Only implemented for scalar x0"
    @assert order >= 0 "Order must be non-negative"
    Acoeffs = Vector{typeof(A(x0))}(undef, 0)
    for i = 0:order
        Ai = zeros(size(A(x0))...)
        for j in eachindex(Ai)
            Ai[j] = central_fdm(5, i)(s -> A(s)[j], x0[1])
        end
        Ai ./= factorial(i)
        push!(Acoeffs, Ai)
    end
    return Acoeffs
end

function balancing(A, B, C, order, x0)
    Ac = compute_taylor_coeffs(A, x0, order)
    Bc = compute_taylor_coeffs(B, x0, order)
    Cc = compute_taylor_coeffs(C, x0, order)
    P0, R = symlyapc(Ac[1], Bc[1])
    Q0, L = symlyapc(Ac[1]', Cc[1]')
    U, Sig, V = svd(R'*L)
    Tl = L * V * diagm(0 => 1 ./ sqrt.(Sig))
    Tr = R * U * diagm(0 => 1 ./ sqrt.(Sig))
    Tr = Tr'
    for l in 1:(order+1)
    end
    return nothing
end

function symlyapc(A, B)
    P = plyapc(A, B)
    Pc = P*P'
    return Pc, P
end
function symlyapc(A::Adjoint, B::Adjoint)
    P = plyapc(A, B)
    Pc = P'*P
    return Pc, P'
end


end
