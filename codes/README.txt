To start an anylsis you need first to open de Terminal and type:

    jupyter notebook

A new tab will open in your browser. Then select the method you are interest on running:
    Functional_enrichment.ipynb
    GSEA.ipynb 

Walktrough: Functional_enrichment.ipynb



Walkthrough: GSEA.ipynb 

Two types of inputs are allowed:
	Two column Excel file with the gene and its ranking value (e.g. FC)
	(see ra)

	Excel matrix where each column is a condition and genes are in the rows, cells are expression values. A GSEA will be perfomed for each condition compared against all the other conditions (columns) in the file. 