all:
	pdflatex -interaction=nonstopmode main.tex > /dev/null
	biber main > /dev/null
	pdflatex -interaction=nonstopmode main.tex | tail

clean:
	rm *.log
	rm *.aux
	rm *.out
	
