# STDecon_benchmark
We conducted a comprehensive benchmarking of 16 existing methods, evaluating them with respect to accuracy, robustness and usability. We compared these methods under varying conditions, including different numbers of genes, spots, and cell subtypes:  
(1) selecting different numbers of genes in the mouse brain cortex (seqFISH+) dataset;  
(2) using smaller bin sizes on the mouse brain medial preoptic area (MERFISH) dataset to simulate more spatial spots;  
(3) dividing each of 5 broad cell type labels into 2, 4, and 8 cell subtype labels on the DV simulated dataset. 
## Implementation
In the ***/methods***, we showed the implementation of all 16 methods on the mouse brain medial preoptic area (MERFISH) dataset as an example. Please make sure the installation of all methods following these versions: CARD (v1.1), Cell2location(v0.1.3), CellDART(v0.1.2), DestVI(v1.0.4), DSTG(v0.0.1), ovoSpaRc(v0.4.4), RCTD(v2.0.1), SD<sup>2</sup>, Seurat(v5.0.1), SpatialDecon(v1.12.3), SpatialDWLS(v1.1.2), SPOTlight (v0.1.0), spSeudoMap(v1.1.0), STdeconvolve(v1.6.0), stereoscope(v.03), Tangram(v1.0.4).  
  
In the ***/evaluation***, we showed the approach used to qualify the performance of each method. In the ***/visualization***, we  showed the visualization method used in our manuscript.
