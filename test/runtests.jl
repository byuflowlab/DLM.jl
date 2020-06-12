using Test, DLM

# Source: "Robustness of the subsonic doublet lattice method" by van Zyl
@testset "AR 20 - 10 X 40 - Van Zyl" begin

    AR = 20

    nchord = 10
    nspan = 40

    # define geometry
    b = 10
    semispan = b/2
    chord = b/AR
    symmetric = true

    # generate geometry
    xyz = zeros(3, nchord+1, nspan+1)
    for i = 1:nchord+1
        for j = 1:nspan+1
            xyz[1, i, j] = chord*(i-1)/nchord
            xyz[2, i, j] = semispan*(j-1)/nspan
        end
    end

    # set flow conditions
    Vinf = 10
    M = 0.0

    # set frequency range
    nfreq = 12
    kr = [0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    # initialize matrices
    AIC = zeros(Complex{Float64}, nchord, nspan, nchord, nspan)
    w = zeros(Complex{Float64}, nchord, nspan)
    CL = zeros(Complex{Float64}, nfreq)

    for kernel in [ParabolicKernel(), QuarticKernel()]

        # loop through frequencies
        for k = 1:nfreq

            # reduced freqency is based on entire chord
            ω = Vinf*kr[k]/(chord)

            # get AIC
            influence_matrix!(AIC, ω, Vinf, M, xyz, true, kernel=kernel)

            # get downwash
            for i = 1:nchord
                for j = 1:nspan
                    # get control point x-location
                    x1 = (1/2)*xyz[1,i,j] + (1/2)*xyz[1,i,j+1]
                    x2 = (1/2)*xyz[1,i+1,j] + (1/2)*xyz[1,i+1,j+1]
                    x = (1/4)*x1 + (3/4)*x2
                    # compute downwash
                    # -1/U*(dh/dt + U*dh/dx), where h = (x - 0.5c)*e^{j*omega*t}
                    w[i,j] = -1.0 - 1im*(ω/Vinf)*(x-0.5*chord)
                end
            end

            # get pressure coefficients
            Cp = pressure_coefficients(AIC, w)

            # solve for lift coefficient
            CL[k] = sum(Cp)/(nspan*nchord)
        end

        # # 10 x 10 - extracted from plot
        # CL_van_Zyl = [5.5+0.0im, 4.4+0.2im, 3.9+1.6im, 3.5+4.1im, 3.2+6.5im, 3.0+8.5im,
        #     2.6+10.3im, 2.2+11.9im, 1.7+13.0im, 1.2+13.7im, 0.6+13.9im, -0.1+13.3im]

        # # 10 x 40 - extracted from plot
        # CL_van_Zyl = [5.5+0.0im, 4.3+0.2im, 3.8+1.5im, 3.4+4.0im, 3.2+6.2im, 2.9+8.2im,
        #     2.6+9.9im, 2.3+11.3im, 2.0+12.4im, 1.7+13.2im, 1.4+13.8im, 1.1+14.0im]

        # # extrapolated solution (grid independent) - extracted from plot
        # CL_van_Zyl = [5.45+0.0im, 4.27+0.2im, 3.80+1.6im, 3.54+4.1im, 3.46+6.4im,
        #     3.42+8.8im, 3.41+1.1im, 3.39+13.4im, 3.38+15.6im, 3.38+17.9im, 3.37+20.2im, 3.34+22.4im]

        if kernel === ParabolicKernel()
            # current outputs (validate new values before changing these)
            CL_test = [5.5+0.0im, 4.3+0.3im, 3.8+1.6im, 3.5+4.1im, 3.2+6.4im, 3.0+8.5im,
                2.7+10.3im, 2.4+11.8im, 2.0+13.1im, 1.6+14.0im, 1.2+14.6im, 0.8+14.7im]
        elseif kernel === QuarticKernel()
            # current outputs (identical to those produced by van Zyl)
            CL_test= [5.5+0.0im, 4.3+0.3im, 3.8+1.6im, 3.4+4.1im, 3.2+6.3im, 2.9+8.3im,
                2.6+10.0im, 2.3+11.4im, 2.0+12.5im, 1.7+13.3im, 1.4+13.8im, 1.1+14.1im]
        end

        # test that validated results have not changed
        for i = 1:nfreq
            @test isapprox(real(CL[i]), real(CL_test[i]), atol=0.1)
            @test isapprox(imag(CL[i]), imag(CL_test[i]), atol=0.1)
        end
    end
end
