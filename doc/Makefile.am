
PDFTEX = pdflatex -shell-escape

lvxman_PDF = lavax-manual.pdf

# DOC_TARGETS = $(lvxman_PDF)

pdf_DATA = $(lvxman_PDF)

$(lvxman_PDF) : lavax-manual.tex
	$(PDFTEX) lavax-manual.tex
