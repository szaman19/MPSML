# You want latexmk to *always* run, because make does not have all the info.
# Also, include non-file targets in .PHONY so they are run regardless of any
# file of the given name existing.


# The first rule in a Makefile is the one executed by default ("make"). It
# should always be the "all" rule, so that "make" and "make all" are identical.
all: manuscript supplementary cover_letter

# CUSTOM BUILD RULES

# In case you didn't know, '$@' is a variable holding the name of the target,
# and '$<' is a variable holding the (first) dependency of a rule.
# "raw2tex" and "dat2tex" are just placeholders for whatever custom steps
# you might have.


# MAIN LATEXMK RULE

# -pdf tells latexmk to generate PDF directly (instead of DVI).
# -pdflatex="" tells latexmk to call a specific backend with specific options.
# -use-make tells latexmk to call make for generating missing files.

# -interaction=nonstopmode keeps the pdflatex backend from stopping at a
# missing file reference and interactively asking you for an alternative.
# -pdf -pdflatex="pdflatex -interaction=nonstopmode"

manuscript: src/manuscript.tex 
	@latexmk -pdf -outdir=build -pdflatex="pdflatex -interaction=nonstopmode" src/manuscript.tex
	@mv build/manuscript.pdf bin/manuscript.pdf

supplementary: src/supplementary.tex
	@latexmk -pdf -outdir=build -pdflatex="pdflatex -interaction=nonstopmode" src/supplementary.tex
	@mv build/supplementary.pdf bin/supplementary.pdf

cover_letter: src/cover_letter.tex
	@latexmk -pdf -outdir=build -pdflatex="pdflatex -interaction=nonstopmode" src/cover_letter.tex
	@mv build/cover_letter.pdf bin/cover_letter.pdf

preview: src/manuscript.tex 
	@latexmk -pdf -outdir=build -pdflatex="pdflatex -interaction=nonstopmode" -pvc src/manuscript.tex
	@cp build/manuscript.pdf bin/manuscript.pdf

clean:
	@echo "Cleaning..."
	rm -vf build/*
	rm -vf bin/*


.PHONY: manuscript all clean preview supplementary cover_letter
