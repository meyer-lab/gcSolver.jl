all: figureJ1.svg

figureJ1.svg:
	julia -e 'using Pkg; Pkg.activate("."); using gcSolver; gcSolver.figureAll()'

figure%.svg:
	julia -e 'using Pkg; Pkg.activate("."); using gcSolver; gcSolver.figure$*()'

clean:
	rm -rf *.svg
