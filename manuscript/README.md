# Integrated morphometric, molecular, and clinical characterization of Parkinson's disease pathology

This repository contains the LaTeX files and figures for the recent preprint "[Integrated morphometric, molecular, and clinical characterization of Parkinson's disease pathology]()."
You can find the full repository in support of the manuscript with code, data, and more on GitHub: [https://github.com/netneurolab/markello_ppmisnf](https://github.com/netneurolab/markello_ppmisnf).

Full citation:
> Markello, R.D., Shafiei, G., Tremblay, C., Postuma, R.B., Dagher, A., and Misic, B. (2020). Integrated morphometric, molecular, and clinical characterization of Parkinson's disease pathology. bioRxiv.

---

If you would (for whatever reason) like to re-build the PDF locally you can do so with `pdflatex` and `bibtex` (YMMV for other compilers):

```bash
pdflatex --interaction=nonstopmode main.tex
bibtex main.tex
pdflatex --interaction=nonstopmode main.tex
pdflatex --interaction=nonstopmode main.tex
```

PDF files in the [`figures`](./figures) directory were created directly from SVGs with [ImageMagick](https://imagemagick.org/index.php) using the following command:

```bash
convert -density 300 figname.svg figname.pdf
```
