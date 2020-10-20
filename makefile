all: figureJ1.svg figureJ2.svg figureJ3.svg figureJ4.svg figureJ5.svg figureJ6.svg figureJ7.svg

figureJ1.svg figureJ2.svg figureJ3.svg figureJ4.svg figureJ5.svg figureJ6.svg figureJ7.svg:
	julia -e 'include("figures.jl")'

clean:
	rm -rf *.svg
