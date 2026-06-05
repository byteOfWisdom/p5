all:
	pdflatex -interaction=nonstopmode main.tex | tail
	biber main > /dev/null
	pdflatex -interaction=nonstopmode main.tex | tail

clean:
	rm *.log
	rm *.aux
	rm *.out
