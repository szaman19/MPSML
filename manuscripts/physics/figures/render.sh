#!/bin/bash
cat > $1_converted.tex << EOF
\documentclass{standalone}
\usepackage{graphicx}
\usepackage{color}
\begin{document}
\def\svgwidth{3.2in}
\input{changeme.pdf_tex}
\end{document}
EOF

sed -i -e "s/changeme/$1/g" $1_converted.tex
pdflatex $1_converted.tex
rm $1_converted.aux $1_converted.log $1_converted.tex $1_converted.tex-e $1_converted.tex
open $1_converted.pdf
