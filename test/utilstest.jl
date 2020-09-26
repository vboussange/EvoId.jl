using DataFrames
""" function test_interpol()
function testing `interpolate_df(df,xlab,ylab,zlab)`
"""
function test_interpol()
    df = DataFrame(zeros(100^2,3),[:x,:y,:z])
    x = cumsum(rand(100));y = cumsum(rand(100))
    for (i,t) in enumerate(Iterators.product(x,y))
        df[ i, [:x,:y]] .= t
        df[i,:z] = sin(sum(t))
    end
    df = df[shuffle(1:100^2),:]
    fun = interpolate_df(df,:x,:y,:z)
    idx = 9999
    return (fun.itp(df.x[idx],df.y[idx]) ≈ df.z[idx])
end
@testset "Interpolations" begin
    @test test_interpol()
end
