using Test, DLM

@testset "Blair" begin

    nchord = 3
    nspan = 3

    # create geometry
    semispan = 12
    chord = 12
    symmetric = true
    xyz = zeros(3, nchord+1, nspan+1)
    for i = 1:nchord+1
        for j = 1:nspan+1
            xyz[1, i, j] = chord*(i-1)/nchord
            xyz[2, i, j] = semispan*(j-1)/nspan
        end
    end
    θ = zeros(nchord, nspan)

    # set flow conditions
    Vinf = 1
    α = 0
    β = 0
    b = 6
    kr = 1
    M = 0.5
    ω = Vinf*kr/b

    # get influence matrix
    AIC = influence_matrix(ω, Vinf, M, xyz, true)

    # set downwash
    w = downwash(xyz, Vinf, α, β, θ) .- 1im

    # solve
    Cp = pressure_coefficients(AIC, w)

    Cp_blair = [-5.4900e-01 + 6.2682e+00im -5.9144e-01 + 5.8092e+00im -5.8286e-01 + 4.5474e+00im;
                -3.8862e+00 + 2.4495e+00im -3.6405e+00 + 2.1530e+00im -2.8983e+00 + 1.4663e+00im;
                -3.8736e+00 + 1.1745e+00im -3.6234e+00 + 1.0281e+00im -2.8893e+00 + 7.1186e-01im]

    # test each panel for agreement
    for i = 1:length(Cp)
        @test isapprox(Cp[i], Cp_blair[i], atol=0.3)
    end
end

@testset "Rodden" begin

    nchord = 10
    nspan = 40

    AR = 20

    # create geometry
    b = 10
    semispan = b/2
    chord = b/AR
    symmetric = true

    area = semispan*chord

    xyz = zeros(3, nchord+1, nspan+1)
    for i = 1:nchord+1
        for j = 1:nspan+1
            xyz[1, i, j] = chord*(i-1)/nchord
            xyz[2, i, j] = semispan*(j-1)/nspan
        end
    end

    # set flow conditions
    Vinf = 1
    α = 0.5*pi/180
    β = 0
    M = 0.0

    # set frequency range
    nfreq = 25
    kr = range(0.0, 2.0, length=nfreq)

    # initialize matrices
    AIC = zeros(Complex{Float64}, nchord, nspan, nchord, nspan)
    w = zeros(Complex{Float64}, nchord, nspan)
    Cl = zeros(Complex{Float64}, nfreq)

    for k = 1:nfreq
        ω = Vinf*kr[k]/(0.5*chord)

        # get AIC
        @time influence_matrix!(AIC, ω, Vinf, M, xyz, true)

        # get downwash
        for i = 1:nchord
            for j = 1:nspan
                # get control point x-location
                x1 = (1/2)*xyz[1,i,j] + (1/2)*xyz[1,i,j+1]
                x2 = (1/2)*xyz[1,i+1,j] + (1/2)*xyz[1,i+1,j+1]
                x = (1/4)*x1 + (3/4)*x2
                # compute downwash
                # -1/U*(dh/dt + U*dh/dx), where h = (x - 0.25c)*e^{j*omega*t}
                w[i,j] = -1.0 - 1im*(ω/Vinf)*(-x-0.25*chord)
            end
        end

        # get pressure coefficients
        Cp = pressure_coefficients(AIC, w)

        # coe[fficient of pressure
        Cl[k] = sum(Cp)/(nchord*nspan)

        println("$k, $(kr[k]), $(real(Cl[k])), $(imag(Cl[k]))")
    end
end
