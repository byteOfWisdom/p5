all:
	pdflatex -interaction=nonstopmode main.tex | tail
	biber main > /dev/null
	pdflatex -interaction=nonstopmode main.tex > /dev/null

clean:
	rm *.log
	rm *.aux
	rm *.out
