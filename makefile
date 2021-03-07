all: figureJ1.svg figureJ2.svg

figureJ1.svg figureJ2.svg:
	julia -e 'include("figures.jl")'

clean:
	rm -rf *.svg
