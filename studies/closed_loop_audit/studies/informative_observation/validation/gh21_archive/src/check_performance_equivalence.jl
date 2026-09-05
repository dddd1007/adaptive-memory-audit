module Before
include("../validation/InformativeObservation_before_type_fix.jl")
end
include("InformativeObservation.jl")
using .InformativeObservation, TOML, Test
c=config(); s=settings(); maximum_difference=0.0
@testset "Type-stability changes preserve likelihood and fixed-data fits" begin
    for lambda in (0.0,1.0,2.0)
        old=Before.InformativeObservation
        d=simulate(c,320270906,lambda,1.0)
        od=old.Dataset(d.z0,d.plan,d.action,d.observed,d.y,d.lambda,d.alpha)
        x=[.25,-.55,.35,1.1]
        f,g,H=objective(x,d,c,gh(21)); of,og,oH=old.objective(x,od,c,old.gh(21))
        @test f==of
        @test g==og
        @test H==oH
    end
end
