Open a Terminal cd/ into codes/ folder and then type: 
    jupyter notebook  

A new tab will open in your browser. Then select the method you are interest on running:
    Functional_enrichment.ipynb
    GSEA.ipynb 

    Once you are in a notebook (the file you clicked in is a notebook)
    You need to run each code cell by clicking at it (inside where the text is) and next clicking the triangle run button


Walktrough: Functional_enrichment.ipynb
    Run first cell

    Run second cell and:
    Select the file you have your gene sets of interest (see example_groups_file.xls)
    Select a file with a list of terms of interest you do not want to display (optional filter, likely applied after you ran one analysis)
    Select a minimum of genes common in your gene set and the pathway to perform analyses (default 0)
    Select corrected p-value threshold (by default < 0.05)

    Run third cell
    If you want to save the heatmap just right click it as you would do with any image and Save as

Walkthrough: GSEA.ipynb 

Two types of inputs are allowed:
	Two column Excel file with the gene and its ranking value (e.g. FC).

	Excel matrix where each column is a condition and genes are in the rows, cells are expression values. A GSEA will be perfomed for each condition compared against all the other conditions (columns) in the file. 