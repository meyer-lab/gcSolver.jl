all: figureJ1.svg figureJ2.svg

figureJ1.svg:
	julia -e 'using Pkg; Pkg.activate("."); using gcSolver; gcSolver.figureJ1()'
    
figureJ2.svg:
	julia -e 'using Pkg; Pkg.activate("."); using gcSolver; gcSolver.figureJ2()'

figure%.svg:
	julia -e 'using Pkg; Pkg.activate("."); using gcSolver; gcSolver.figure$*()'

clean:
	rm -rf *.svg
