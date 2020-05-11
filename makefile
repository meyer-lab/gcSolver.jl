
all: figureJ1.svg
venv: venv/bin/activate

venv/bin/activate: requirements.txt
	test -d venv || virtualenv venv
	. venv/bin/activate && pip install -Uqr requirements.txt
	touch venv/bin/activate

figureJ1.svg:
	julia -e 'using Pkg; Pkg.activate("."); using gcSolver; gcSolver.figureAll()'

figure%.svg:
	julia -e 'using Pkg; Pkg.activate("."); using gcSolver; gcSolver.figure$*()'

coverage.cob:
	julia -e 'using Pkg; Pkg.add("Coverage"); using Coverage; Pkg.activate("."); Pkg.test("FcgR"; coverage=true); coverage = process_folder(); LCOV.writefile("coverage-lcov.info", coverage)'
	pip3 install --user lcov_cobertura
	python3 ~/.local/lib/python3.7/site-packages/lcov_cobertura.py coverage-lcov.info -o coverage.cob

clean:
	rm -rf *.svg venv output
