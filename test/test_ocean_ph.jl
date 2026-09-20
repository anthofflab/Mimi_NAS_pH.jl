@testitem "Ocean pH" begin

    using Mimi

    # The five SSP CO₂ concentration scenarios the data file provides, in order
    # of increasing forcing.
    scenarios = ["SSP1-19", "SSP1-26", "SSP2-45", "SSP3-70", "SSP5-85"]

    pH_0 = 8.123 # the initial condition get_model sets

    final_pH = Dict{String, Float64}()

    for scenario in scenarios
        m = Mimi_NAS_pH.get_model(ssp_emissions_scenario = scenario)
        run(m)

        @test Mimi.time_labels(m) == collect(1950:2100)

        pH = m[:ocean_pH, :pH]

        @test length(pH) == 151
        @test all(isfinite, pH)

        # Timestep one is the initial condition, everything after it is computed.
        @test pH[1] == pH_0

        # Globally averaged ocean pH stays in a plausible range.
        @test all(p -> 7.5 < p < 8.2, pH[2:end])

        final_pH[scenario] = pH[end]
    end

    # More CO₂ means a more acidic ocean, so pH in 2100 should fall as the
    # scenarios get more forcing.
    @test issorted([final_pH[s] for s in scenarios], rev = true)
end

@testitem "pH approximation" begin

    # Equation 7 in Appendix F of the report, using the coefficients get_model
    # sets. The component deliberately lags the atmospheric CO₂ concentration by
    # one year, so pH in period t uses concentrations from period t-1.
    β1 = -0.3671
    β2 = 10.2328

    scenario = "SSP2-45"

    m = Mimi_NAS_pH.get_model(ssp_emissions_scenario = scenario)
    run(m)

    pH = m[:ocean_pH, :pH]
    atm_co2_conc = Mimi_NAS_pH.co2_data[:, scenario]

    @test pH[2:end] ≈ β1 .* log.(atm_co2_conc[1:end-1]) .+ β2
end
