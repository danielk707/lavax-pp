SUBDIRS = src

if HAVE_PDFLATEX

  lavax-manual : doc/lavax-manual.tex
	cd doc && $(PDFLATEX) lavax-manual.tex && cd ..

endif
