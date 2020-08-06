all: figureJ1.svg figureJ2.svg figureJ3.svg figureJ4.svg figureJ5.svg

figureJ1.svg:
	julia -e 'include("figures/figure1.jl"); figureJ1.jl()'
    
figureJ2.svg:
	julia -e 'include("figures/figure2.jl"); figureJ2.jl()'

figureJ3.svg:
	julia -e 'include("figures/figure3.jl"); figureJ3.jl()'

figureJ4.svg:
	julia -e 'include("figures/figure4.jl"); figureJ4.jl()'

figureJ5.svg:
	julia -e 'include("figures/figure5.jl"); figureJ5.jl()'

figure%.svg:
	julia -e 'using Pkg; Pkg.activate("."); using gcSolver; gcSolver.figure$*()'

clean:
	rm -rf *.svg
