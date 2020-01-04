using ForwardDiff


@testset "Test accuracy of FS gradients." begin
    fd_res = ForwardDiff.jacobian(x -> runCkine([100.0], x), rxntfR)

    x, dx = runCkineFS(100.0, rxntfR)

    println(fd_res .- dx)
end


@testset "Profile forward sensitivities." begin
    println("runCkineFS")
    @time runCkineFS(100.0, rxntfR)

    @profile runCkineFS(100.0, rxntfR)
    Profile.print(noisefloor = 2.0)
end