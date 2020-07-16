all: figureJ1.svg figureJ2.svg figureJ3.svg figureJ5.svg

figureJ1.svg:
	julia -e 'using Pkg; Pkg.activate("."); using gcSolver; gcSolver.figureJ1()'
    
figureJ2.svg:
	julia -e 'using Pkg; Pkg.activate("."); using gcSolver; gcSolver.figureJ2()'

figureJ3.svg:
	julia -e 'using Pkg; Pkg.activate("."); using gcSolver; gcSolver.figureJ3()'

figureJ5.svg:
	julia -e 'using Pkg; Pkg.activate("."); using gcSolver; gcSolver.figureJ5()'

figure%.svg:
	julia -e 'using Pkg; Pkg.activate("."); using gcSolver; gcSolver.figure$*()'

clean:
	rm -rf *.svg
