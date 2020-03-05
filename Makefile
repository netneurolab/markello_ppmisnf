.PHONY: all help preprocess analysis results manuscript

PYTHON ?= python
JUPYTER ?= jupyter

all: analysis results manuscript

help:
	@echo "Please use 'make <target>' where <target> is one of:"
	@echo "  preprocess    to preprocess and parcellate neuroimaging data for analysis"
	@echo "  analysis      to run the primary computational analyses"
	@echo "  results       to run all results-generating code and generate PNG/SVG figures"
	@echo "  manuscript    to compile a PDF from the manuscript TeX files"
	@echo "  all           to run 'analysis', 'results', and 'manuscript'"

preprocess:
	@echo "Running data preprocessing\n"
	$(PYTHON) scripts/01_preprocess/01_parcellate_antslct.py
	$(PYTHON) scripts/01_preprocess/02_get_pdica_zscore.py

analysis:
	@echo "Running data analyses\n"
	$(PYTHON) scripts/02_analysis/01_prepare_snf_data.py
	$(PYTHON) scripts/02_analysis/02_snf_gridsearch.py

results:
	@echo "Executing all results code to generate figures\n"
	$(PYTHON) scripts/03_results/01_concatenation_vs_snf.py
	$(PYTHON) scripts/03_results/02_pd_patient_biotypes.py
	$(PYTHON) scripts/03_results/03_mri_contributions.py
	$(PYTHON) scripts/03_results/04_diffusion_embedding.py
	$(PYTHON) scripts/03_results/05_supplementary_results.py

manuscript:
	@echo "Generating PDF with pdflatex + bibtex"
	@cd manuscript && \
	 rm -f main.pdf && \
	 pdflatex --interaction=nonstopmode main > /dev/null && \
	 bibtex main > /dev/null && \
	 pdflatex --interaction=nonstopmode main > /dev/null && \
	 pdflatex --interaction=nonstopmode main > /dev/null && \
	 rm -f main.aux main.bbl main.blg main.log main.out mainNotes.bib main.synctex.gz
	@echo "Check ./manuscript/main.pdf for generated file"
