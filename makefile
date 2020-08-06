all: figureJ1.svg figureJ2.svg figureJ3.svg figureJ4.svg figureJ5.svg

figureJ1.svg:
	julia -e 'using Pkg; Pkg.activate("."); include("figures/figureJ1.jl"); figureJ1()'
    
figureJ2.svg:
	julia -e 'using Pkg; Pkg.activate("."); include("figures/figureJ2.jl"); figureJ2()'

figureJ3.svg:
	julia -e 'using Pkg; Pkg.activate("."); include("figures/figureJ3.jl"); figureJ3()'

figureJ4.svg:
	julia -e 'using Pkg; Pkg.activate("."); include("figures/figureJ4.jl"); figureJ4()'

figureJ5.svg:
	julia -e 'using Pkg; Pkg.activate("."); include("figures/figureJ5.jl"); figureJ5()'

figure%.svg:
	julia -e 'using Pkg; Pkg.activate("."); using gcSolver; gcSolver.figure$*()'

clean:
	rm -rf *.svg
