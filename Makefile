#
# Makefile for LaTeX project IEI_zh
#
# Reference:
# 	1. http://tex.stackexchange.com/questions/40738/
# 	2. Manual of latexmk
#
DOC = IEI_zh

.PHONY: $(DOC).pdf all clean

all: $(DOC).pdf

$(DOC).pdf: $(DOC).tex
	latexmk -xelatex -shell-escape $^

clean:
	latexmk -c
	-rm -r _minted-$(DOC) *.minted* *.pyg
