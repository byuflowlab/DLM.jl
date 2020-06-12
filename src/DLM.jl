module DLM

export influence_matrix, influence_matrix!
export pressure_coefficients
export AbstractKernel, ParabolicKernel, QuarticKernel

# coefficients for 12-term approximation of kernel integrals
const a = (0.000319759140, -0.000055461471, 0.002726074362, 0.005749551566,
    0.031455895072, 0.106031126212, 0.406838011567, 0.798112357155,
    -0.417749229098, 0.077480713894, -0.012677284771, 0.001787032960)
const b = 0.009054814793

# kernel approximations
"""
    AbstractKernel()

Abstract type for approximation of the planar and nonplanar influence terms in
the numerator of the doublet lattice method kernel
"""
abstract type AbstractKernel end

"""
    ParabolicKernel()

Parabolic approximation for the planar and nonplanar influence terms in the
numerator of the doublet lattice method kernel.
"""
struct ParabolicKernel <: AbstractKernel end

"""
    QuarticKernel()

Quartic approximation for the planar and nonplanar influence terms in the
numerator of the doublet lattice method kernel.
"""
struct QuarticKernel <: AbstractKernel end

"""
    approximate_kernel_integrals(u, k)

Solves for the integrals I0 and J0 using a 12-term exponential approximation
given by Desmarais in "AN ACCURATE METHOD FOR EVALUATING THE KERNEL OF THE
INTEGRAL EQUATION RELATING LIFT TO DOWNWASH IN UNSTEADY POTENTIAL FLOW".

# Arguments
- u = (M*R - x0)/(β^2*x0)
- k = ω*r/U
"""
function approximate_kernel_integrals(u, k)

    TF = promote_type(typeof(u), typeof(k))

    I0 = zero(complex(TF))
    J0 = zero(complex(TF))
    @inbounds for j = 1:12
        p = b*2^j
        pk = complex(p, k)
        tmp = a[j]*exp(-pk*u)/pk
        I0 += tmp
        J0 += tmp*(u+1/pk)
    end

    return I0, J0
end

"""
    kernel_coefficients(x0, y0, z0, ω, U, M)

Solves for the kernel coefficients K1 and K2 and their steady state counterparts
K10 and K20.

See "A Doublet-Lattice Method for Calculating Lift Distributions
on Oscillating Surfaces in Subsonic Flows" by Edward Albano and William P Rodden.

# Arguments
- x0 = x-ξ
- y0 = y-η
- z0 = z-ζ
- ω = frequency
- U = Velocity
- M = Mach Number
"""
function kernel_coefficients(x0, y0, z0, ω, U, M)

    r = sqrt(y0^2 + z0^2)
    β = sqrt(1 - M^2)
    R = sqrt(x0^2+β^2*r^2)

    u = (M*R-x0)/(β^2*max(r, eps(typeof(r))))
    k = ω*r/U

    # use the symmetry properties of the integrand if u < 0
    if u < 0
        I0, J0 = approximate_kernel_integrals(0, k)
        I10 = 1 - 1im*k*I0
        I20 = 2 - 1im*k*I0 + k^2*J0
        u1pos = -u
        I0, J0 = approximate_kernel_integrals(u1pos, k)
        invR = 1/R
        invsqrt = 1/sqrt(1+u1pos^2)
        expk = exp(-1im*k*u1pos)
        I11 = (1-u1pos*invsqrt)*expk - 1im*k*I0
        I21 = ((2 + 1im*k*u1pos)*(1 - u1pos*invsqrt) -
            u1pos*invsqrt^3)*expk - 1im*k*I0 + J0*k^2
        I1 = complex(2*real(I10) - real(I11), imag(I11))
        I2 = complex(2*real(I20) - real(I21), imag(I21))
        expk = exp(-1im*k*u)
    else
        I0, J0 = approximate_kernel_integrals(u, k)
        invR = 1/R
        invsqrt = 1/sqrt(1+u^2)
        expk = exp(-1im*k*u)
        I1 = (1-u*invsqrt)*expk - 1im*k*I0
        I2 = ((2 + 1im*k*u)*(1 - u*invsqrt) -
            u*invsqrt^3)*expk - 1im*k*I0 + J0*k^2
    end

    K1 = I1 + M*r*invR*invsqrt*expk
    K2 = -I2 - 1im*k*invsqrt*expk*(M*r*invR)^2 -
       (M*r*invR)*((1 + u^2)*(β*r*invR)^2 +
       2 + M*r*u*invR)*expk*invsqrt^3

    K10 =  1 + x0/R
    K20 = -2 - x0/R*(2+β^2*r^2/R^2)

    return K1, K2, K10, K20
end

"""
    kernel_numerator(ω, U, M, x0, y0, z0, γr, γs)

Solves for the numerator terms of the kernel and subtracts out the steady-state
components of these terms.

See "Further Refinement of the Subsonic Doublet-Lattice Method" by William P. Rodden,
Paul F. Taylor and Samuel C. McIntosh Jr.
"""
function kernel_numerator(ω, U, M, x0, y0, z0, cγr, sγr, cγs, sγs)
    T1, T2 = cγr*cγs + sγr*sγs, (z0*cγr-y0*sγr)*(z0*cγs-y0*sγs)
    K1, K2, K10, K20 = kernel_coefficients(x0, y0, z0, ω, U, M)
    expn = exp(complex(0, -ω*x0/U))
    return (K1*expn-K10)*T1, (K2*expn-K20)*T2
end

"""
influence_coefficient(ω, U, M, (xr, yr, zr), (xi, yi, zi), (xo, yo, zo),
    chord, λ, γr, γs, kernel)

Solves for the influence coefficient of the sending panel on the receiving panel.

The planar and nonplanar influence terms in the numerator are approximated using
the approximation specified by `kernel`.

See "Further Refinement of the Subsonic Doublet-Lattice Method" by William P. Rodden,
Paul F. Taylor and Samuel C. McIntosh Jr.
"""
influence_coefficient

function influence_coefficient(ω, U, M, (xr, yr, zr),
    (xi, yi, zi), (xo, yo, zo), chord, cλ, cγr, sγr, cγs, sγs, ::ParabolicKernel)

    dx = xo-xi
    dy = yo-yi
    dz = zo-zi
    l = sqrt(dx^2 + dy^2 + dz^2)

    xm = (xo+xi)/2
    ym = (yo+yi)/2
    zm = (zo+zi)/2

    # Approximate numerator coefficients as parabolas
    x0i = xr - xi
    y0i = yr - yi
    z0i = zr - zi
    Ki1, Ki2 = kernel_numerator(ω, U, M, x0i, y0i, z0i, cγr, sγr, cγs, sγs)

    x0m = xr - xm
    y0m = yr - ym
    z0m = zr - zm

    Km1, Km2 = kernel_numerator(ω, U, M, x0m, y0m, z0m, cγr, sγr, cγs, sγs)

    x0o = xr - xo
    y0o = yr - yo
    z0o = zr - zo
    Ko1, Ko2 = kernel_numerator(ω, U, M, x0o, y0o, z0o, cγr, sγr, cγs, sγs)

    e = 1/2*l*cλ
    η = y0m*cγs + z0m*sγs
    ζ = -y0m*sγs + z0m*cγs

    if abs(ζ) == zero(ζ)
        F = 2*e/(η^2-e^2)
    else
        F =atan(2*e*abs(ζ),(η^2+ζ^2-e^2))/abs(ζ)
    end

    # steady term: from VLM
    d0 = steady_influence_coefficient(M, x0i, y0i, z0i, x0o, y0o, z0o, sγr, cγr)

    # planar term: ∫(A1*y^2 + B1*y + C1)/((η-y)^2+ζ^2)
    A1 = (Ki1 - 2*Km1 + Ko1)/(2*e^2)
    B1 = (Ko1 - Ki1)/(2*e)
    C1 = Km1

    d1 = ((η^2-ζ^2)*A1 + η*B1 + C1) * F +
        (1/2*B1 + η*A1) * log(((η-e)^2+ζ^2)/((η+e)^2 + ζ^2)) + 2*e*A1

    # nonplanar term: ∫(A1*y^2 + B1*y + C1)/((η-y)^2+ζ^2)^2
    if abs(ζ) == zero(ζ)
        D2 = zero(complex(typeof(d1)))
    else
        A2 = (Ki2 - 2*Km2 + Ko2)/(2*e^2)
        B2 = (Ko2 - Ki2)/(2*e)
        C2 = Km2

        α = (e/ζ)^2 * (1 - (η^2 + ζ^2 - e^2)/(2*e)*F)

        D2 = e/(η^2 + ζ^2 - e^2)*((
            (2*(η^2 + ζ^2 + e^2)*(e^2*A2 + C2) + 4*η*e^2*B2))/
            (((η + e)^2 + ζ^2)*((η - e)^2 + ζ^2)) -
             (α/e^2)*((η^2 + ζ^2)*A2 + η*B2 + C2))
    end

    return chord/(8*pi)*(d0 + d1 + D2)
end

function influence_coefficient(ω, U, M, (xr, yr, zr),
    (xi, yi, zi), (xo, yo, zo), chord, cλ, cγr, sγr, cγs, sγs, ::QuarticKernel)

    dx = xo-xi
    dy = yo-yi
    dz = zo-zi
    l = sqrt(dx^2 + dy^2 + dz^2)

    xm = (xi+xo)/2
    ym = (yi+yo)/2
    zm = (zi+zo)/2

    xmi = (xi+xm)/2
    ymi = (yi+ym)/2
    zmi = (zi+zm)/2

    xmo = (xm+xo)/2
    ymo = (ym+yo)/2
    zmo = (zm+zo)/2

    # Approximate numerator coefficients as quartic
    x0i = xr - xi
    y0i = yr - yi
    z0i = zr - zi
    Ki1, Ki2 = kernel_numerator(ω, U, M, x0i, y0i, z0i, cγr, sγr, cγs, sγs)

    x0mi = xr - xmi
    y0mi = yr - ymi
    z0mi = zr - zmi

    Kmi1, Kmi2 = kernel_numerator(ω, U, M, x0mi, y0mi, z0mi, cγr, sγr, cγs, sγs)

    x0m = xr - xm
    y0m = yr - ym
    z0m = zr - zm

    Km1, Km2 = kernel_numerator(ω, U, M, x0m, y0m, z0m, cγr, sγr, cγs, sγs)

    x0mo = xr - xmo
    y0mo = yr - ymo
    z0mo = zr - zmo

    Kmo1, Kmo2 = kernel_numerator(ω, U, M, x0mo, y0mo, z0mo, cγr, sγr, cγs, sγs)

    x0o = xr - xo
    y0o = yr - yo
    z0o = zr - zo
    Ko1, Ko2 = kernel_numerator(ω, U, M, x0o, y0o, z0o, cγr, sγr, cγs, sγs)

    e = 1/2*l*cλ
    η = y0m*cγs + z0m*sγs
    ζ = -y0m*sγs + z0m*cγs

    if abs(ζ) == zero(ζ)
        F = 2*e/(η^2-e^2)
    else
        F =atan(2*e*abs(ζ), (η^2+ζ^2-e^2))/abs(ζ)
    end

    # steady term: from VLM
    d0 = steady_influence_coefficient(M, x0i, y0i, z0i, x0o, y0o, z0o, sγr, cγr)

    # planar term: ∫(A1*y^4 + B1*y^3 + C1*y^2 + D1*y + E1)/((η-y)^2+ζ^2)
    A1 = -(1/(6*e^2))*(Ki1 - 16*Kmi1 + 30*Km1 - 16*Kmo1 + Ko1)
    B1 =  (1/(6*e  ))*(Ki1 -  8*Kmi1          +  8*Kmo1 - Ko1)
    C1 = Km1
    D1 = -(2/(3*e^3))*(Ki1 -  2*Kmi1          +  2*Kmo1 - Ko1)
    E1 =  (2/(3*e^4))*(Ki1 -  4*Kmi1 +  6*Km1 -  4*Kmo1 + Ko1)

    d1 = ((η^2-ζ^2)*A1 + η*B1 + C1 + η*(η^2 - 3*ζ^2)*D1 +
        (η^4 - 6*η^2*ζ^2 + ζ^4)*E1)*F + (η*A1 + 1/2*B1 +
        1/2*(3*η^2 - ζ^2)*D1 + 2*η*(η^2-ζ^2)*E1)*log(((η - e)^2 + ζ^2)/((η + e)^2 + ζ^2)) +
        2*e*(A1 + 2*η*D1 + (3*η^2 - ζ^2 + 1/3*e^2)*E1)

    # nonplanar term: ∫(A2*y^4 + B2*y^3 + C2*y^2 + D2*y + E2)/((η-y)^2+ζ^2)^2
    if abs(ζ) == zero(ζ)
        d2 = zero(complex(typeof(d1)))
    else
        A2 = -(1/(6*e^2))*(Ki2 - 16*Kmi2 + 30*Km2 - 16*Kmo2 + Ko2)
        B2 =  (1/(6*e  ))*(Ki2 -  8*Kmi2          +  8*Kmo2 - Ko2)
        C2 = Km2
        D2 = -(2/(3*e^3))*(Ki2 -  2*Kmi2          +  2*Kmo2 - Ko2)
        E2 =  (2/(3*e^4))*(Ki2 -  4*Kmi2 +  6*Km2 -  4*Kmo2 + Ko2)

        d2 = 1/(2*ζ^2)*(((η^2 + ζ^2)*A2 + η*B2 + C2 + η*(η^2 + 3*ζ^2)*D2 +
            (η^4 + 6*η^2*ζ^2 - 3*ζ^4)*E2)*F + 1/((η + e)^2 + ζ^2)*(((η^2 + ζ^2)*η +
            (η^2 - ζ^2)*e)*A2 + (η^2 + ζ^2 + η*e)*B2 + (η + e)*C2 +
            (η^4 - ζ^4 + (η^2 - 3*ζ^2)*η*e)*D2 + ((η^4 - 2*η^2*ζ^2 - 3*ζ^4)*η +
            (η^4 - 6*η^2*ζ^2 + ζ^4)*e)*E2) - 1/((η - e)^2 + ζ^2)*(((η^2 + ζ^2)*η -
            (η^2 - ζ^2)*e)*A2 + (η^2 + ζ^2 - η*e)*B2 + (η - e)*C2 +
            (η^4 - ζ^4 - (η^2 - 3*ζ^2)*η*e)*D2 + ((η^4 - 2*η^2*ζ^2 - 3*ζ^4)*η -
            (η^4 - 6*η^2*ζ^2 + ζ^4)*e)*E2) + (ζ^2*log(((η - e)^2 + ζ^2)/((η+e)^2 + ζ^2)))*D2 +
            4*ζ^2*(e + η*log(((η - e)^2 + ζ^2)/((η + e)^2 + ζ^2)))*E2)
    end

    return chord/(8*pi)*(d0 + d1 + d2)
end

"""
    steady_influence_coefficient(M, x0i, y0i, z0i, x0o, y0o, z0o, sγr, cγr)

Solves for the steady component of the influence of the sending panel on the
receiving panel using the vortex lattice method.
"""
function steady_influence_coefficient(M, x0i, y0i, z0i, x0o, y0o, z0o, sγr, cγr)

    # steady term (from VLM)
    ax = x0i/sqrt(1-M^2)
    ay = y0i
    az = z0i
    anrm = sqrt(ax^2 + ay^2 + az^2)
    ainv = 1/(anrm*(anrm - ax))

    bx = x0o/sqrt(1-M^2)
    by = y0o
    bz = z0o
    bnrm = sqrt(bx^2 + by^2 + bz^2)
    binv = 1/(bnrm*(bnrm - bx))

    vy =  az*ainv - bz*binv
    vz = -ay*ainv + by*binv

    ainv = 1/(anrm*bnrm*(anrm*bnrm + ax*bx + ay*by + az*bz))
    vy = vy + (az*bx - ax*bz)*(anrm + bnrm)*ainv
    vz = vz + (ax*by - ay*bx)*(anrm + bnrm)*ainv

    return -(sγr*vy - cγr*vz)
end

"""
    influence_matrix(ω, U, M, xyz::AbstractVector{<:AbstractArray{<:Number, 3}}, symmetric::AbstractVector; kernel=QuarticKernel())
    influence_matrix(ω, U, M, xyz::AbstractArray{<:Number, 3}, symmetric; kernel=QuarticKernel())
    influence_matrix(ω, U, M, xyzr::AbstractArray{<:Number, 3}, xyzs::AbstractArray{<:Number, 3}, symmetric; kernel=QuarticKernel())

Fills in aerodynamic influence coefficient matrix. Returns aerodynamic influence
coefficient matrix AIC with shape (nir, njr, nis, njs). For a non-allocating version,
see `influence_matrix!`.

Receiving and sending panel points must be aligned with the x-direction. The
freestream velocity must also be aligned in the x-direction.

# Arguments
- `ω`: oscillation frequency
- `U`: freestream velocity
- `M`: Mach number (for Prandtl-Glauert compressibility correction)
- `xyzr`: definition of receiving panel points with shape: (3, nir, njr)
- `xyzs`: definition of sending panel points with shape: (3, nis, njs)
- `symmetric`: flag indicating symmetry of each sending panel
- `kernel`: indicates which approximation of the kernel should be used
"""
influence_matrix

function influence_matrix(ω, U, M, xyz::AbstractVector{<:AbstractArray{<:Number, 3}},
        symmetric::AbstractVector; kernel=QuarticKernel())

    nx = sum(size.(xyz, 2))
    ny = sum(size.(xyz, 3))

    TF = promote_type(typeof(ω), typeof(U), typeof(M), eltype(eltype(xyz)))

    AIC = zeros(complex(TF), nx, ny, nx, ny)

    return influence_matrix!(AIC, ω, U, M, xyz, symmetric, kernel=kernel)
end

function influence_matrix(ω, U, M, xyzr::AbstractArray{<:Number, 3},
    xyzs::AbstractArray{<:Number, 3}, symmetric; kernel=QuarticKernel())

    nir = size(xyzr, 2)-1
    njr = size(xyzr, 3)-1
    nis = size(xyzs, 2)-1
    njs = size(xyzs, 3)-1

    TF = promote_type(typeof(ω), typeof(U), typeof(M), eltype(xyzr), eltype(xyzs))

    AIC = zeros(complex(TF), nir, njr, nis, njs)

    return influence_matrix!(AIC, ω, U, M, xyzr, xyzs, symmetric, kernel=kernel)
end

influence_matrix(ω, U, M, xyz::AbstractArray{<:Number, 3}, symmetric; kernel=QuarticKernel()) = influence_matrix(ω, U, M, xyz, xyz, symmetric, kernel=kernel)


"""
    influence_matrix!(AIC, ω, U, M, xyz::AbstractVector{<:AbstractArray{<:Number, 3}}, symmetric::AbstractVector; kernel=QuarticKernel())
    influence_matrix!(AIC, ω, U, M, xyz::AbstractArray{<:Number, 3}, symmetric; kernel=QuarticKernel())
    influence_matrix!(AIC, ω, U, M, xyzr::AbstractArray{<:Number, 3}, xyzs::AbstractArray{<:Number, 3}, symmetric; kernel=QuarticKernel())

Non-allocating version of [`influence_matrix`](@ref)
"""
influence_matrix!

function influence_matrix!(AIC, ω, U, M, xyz::AbstractVector{<:AbstractArray{<:Number, 3}}, symmetric::AbstractVector; kernel=QuarticKernel())

    nx = size(AIC, 1)
    ny = size(AIC, 2)
    nsurf = length(xyz)

    ir = 0
    jr = 0
    is = 0
    js = 0
    @inbounds for kr = 1:nsurf, ks = 1:nsurf
        nir = size(xyz[ir], 2)
        njr = size(xyz[ir], 3)
        nis = size(xyz[is], 2)
        njs = size(xyz[is], 3)

        ixr = ir+1:ir+nir
        iyr = jr+1:jr+njr
        ixs = is+1:is+nis
        iys = js+1:js+njs

        influence_matrix!(view(AIC, ixr, iyr, ixs, iys), ω, U, M, xyz[kr], xyz[ks], symmetric[ks], kernel=kernel)
    end

    return AIC
end

function influence_matrix!(AIC, ω, U, M, xyzr::AbstractArray{<:Number, 3}, xyzs::AbstractArray{<:Number, 3}, symmetric; kernel=QuarticKernel())
    nir = size(xyzr, 2)-1
    njr = size(xyzr, 3)-1
    nis = size(xyzs, 2)-1
    njs = size(xyzs, 3)-1

    @assert size(AIC) == (nir, njr, nis, njs)

    @inbounds for jr=1:njr, ir = 1:nir, js=1:njs, is=1:nis
        # receiving point at 3/4 chord
        xri = (1/4)*xyzr[1, ir, jr] + (3/4)*xyzr[1, ir+1, jr]
        yri = (1/4)*xyzr[2, ir, jr] + (3/4)*xyzr[2, ir+1, jr]
        zri = (1/4)*xyzr[3, ir, jr] + (3/4)*xyzr[3, ir+1, jr]
        xro = (1/4)*xyzr[1, ir, jr+1] + (3/4)*xyzr[1, ir+1, jr+1]
        yro = (1/4)*xyzr[2, ir, jr+1] + (3/4)*xyzr[2, ir+1, jr+1]
        zro = (1/4)*xyzr[3, ir, jr+1] + (3/4)*xyzr[3, ir+1, jr+1]
        xr = (xri+xro)/2
        yr = (yri+yro)/2
        zr = (zri+zro)/2
        # sending line at 1/4 chord
        xsi = (3/4)*xyzs[1, is, js] + (1/4)*xyzs[1, is+1, js]
        ysi = (3/4)*xyzs[2, is, js] + (1/4)*xyzs[2, is+1, js]
        zsi = (3/4)*xyzs[3, is, js] + (1/4)*xyzs[3, is+1, js]
        xso = (3/4)*xyzs[1, is, js+1] + (1/4)*xyzs[1, is+1, js+1]
        yso = (3/4)*xyzs[2, is, js+1] + (1/4)*xyzs[2, is+1, js+1]
        zso = (3/4)*xyzs[3, is, js+1] + (1/4)*xyzs[3, is+1, js+1]
        # average chord length of sending panel
        chord = (xyzs[1, is+1, js] - xyzs[1, is, js])/2 + (xyzs[1, is+1, js+1] - xyzs[1, is, js+1])/2
        # cos(sweep) of sending doublet line
        cλs = (yso-ysi)/(sqrt((xso-xsi)^2 + (yso-ysi)^2))
        # cos(dihedral) and sin(dihedral) of receiving point
        tmp = sqrt((yro-yri)^2+(zro-zri)^2)
        cγr = (yro-yri)/tmp
        sγr = (zro-zri)/tmp
        # cos(dihedral) and sin(dihedral) of sending line
        tmp = sqrt((yso-ysi)^2+(zso-zsi)^2)
        cγs = (yso-ysi)/tmp
        sγs = (zso-zsi)/tmp
        # compute influence of the sending panel on the receiving panel
        AIC[ir, jr, is, js] = influence_coefficient(ω, U, M, (xr, yr, zr),
            (xsi, ysi, zsi), (xso, yso, zso), chord, cλs, cγr, sγr, cγs, sγs, kernel)
        # add influence of reflected sending panel, if applicable
        if symmetric
            AIC[ir, jr, is, js] += influence_coefficient(ω, U, M, (xr, yr, zr),
                (xso, -yso, zso), (xsi, -ysi, zsi), chord, cλs, cγr, sγr, cγs, -sγs, kernel)
        end
    end

    return AIC
end

influence_matrix!(AIC, ω, U, M, xyz::AbstractArray{<:Number, 3}, symmetric; kernel=QuarticKernel()) = influence_matrix!(AIC, ω, U, M, xyz, xyz, symmetric, kernel=kernel)

"""
    pressure_coefficients(AIC, w)

Compute pressure coefficients `Cp` on each panel, given the downwash w.

# Arguments:
- `AIC`: aerodynamic influence coefficient matrix with shape (nir, njr, nis, njs)
- `w`: downwash on each panel, with shape (nir, njr)
"""
function pressure_coefficients(AIC, w)

    nir, njr, nis, njs = size(AIC)

    AIC_r = reshape(AIC, nir*njr, nis*njs)
    w_r = reshape(w, nir*njr)

    Cp = AIC_r\w_r

    return reshape(Cp, nir, njr)
end

# """
#     downwash!(w, xyz, Vinf, α, β)
#
# Compute the downwash on each of the panels in `xyz` given the freestream
# velocity `Vinf`, angle of attack `α`, sideslip angle `β`, and panel twist `θ`.
# Panel twist has the shape (ni, nj). Non-allocating.
# """
# function downwash!(w, xyz, Vinf, α, β, θ)
#
#     ca, sa = cos(α), sin(α)
#     cb, sb = cos(β), sin(β)
#
#     V = (Vinf*ca*cb, -Vinf*sb, Vinf*sa*cb)
#
#     ni = size(xyz, 2)-1
#     nj = size(xyz, 3)-1
#
#     @assert size(w) == (ni, nj)
#
#     @inbounds for j = 1:nj, i = 1:ni
#         # compute cos(dihedral) and sin(dihedral)
#         yri = (1/4)*xyz[2, i, j] + (3/4)*xyz[2, i+1, j]
#         zri = (1/4)*xyz[3, i, j] + (3/4)*xyz[3, i+1, j]
#         yro = (1/4)*xyz[2, i, j+1] + (3/4)*xyz[2, i+1, j+1]
#         zro = (1/4)*xyz[3, i, j+1] + (3/4)*xyz[3, i+1, j+1]
#         tmp = sqrt((yro-yri)^2+(zro-zri)^2)
#         cd = (yro-yri)/tmp
#         sd = (zro-zri)/tmp
#         # compute cos(theta) and sin(theta)
#         ct = cos(θ[i,j])
#         st = sin(θ[i,j])
#         # compute downwash on each panel
#         w[i, j] = V[1]*st + V[2]*(-ct*sd) + V[3]*(ct*cd)
#     end
#     return w
# end

# """
#     downwash(xyz, Vinf, α, β)
#
# Compute the downwash on each of the panels in `xyz` given the freestream
# velocity `Vinf`, angle of attack `α`, sideslip angle `β`, and panel twist `θ`.
# Panel twist has the shape (ni, nj). For a non-allocating version see `downwash!`.
# """
# function downwash(xyz, Vinf, α, β, θ)
#
#     TF = promote_type(eltype(xyz), typeof(Vinf), typeof(α), typeof(β), eltype(θ))
#
#     w = zeros(TF, size(xyz,2)-1, size(xyz,3)-1)
#
#     return downwash!(w, xyz, Vinf, α, β, θ)
# end

end # module
