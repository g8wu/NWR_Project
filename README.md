# ConGen
### Using DAVID for GO Analysis
Used [DAVID 6.7](https://david-d.ncifcrf.gov) 
<br />Input lists of 100 top upregulated and downregulated genes are all in the DAVID folder.
<br /> Resulting webpages are also saved into the folder from Firefox browser.

1. Go to 'Functional Analysis'
2. Under List tab upload '100upreg.txt'/'100doreg.txt', select 'Gene list' > Submit
3. If DAVID gene 3000 threshold was hit, select under List tab 'Homo sapien'
   If DAVID 3000 gene threshold was still hit, select under Backgrounds > Affymetrix 3' IVT Backgrounds > Human Genome U133A Array
4. Select 'Clear All', under Gene_Ontology check 'GOTERM_BP_DIRECT' > Functional Annotation Clustering
5. Webpage htmls saved to DAVID folder

### Pairwise Analysis
<p float="left">
  <img src='/pictures/scatter_fibro_ipsc.png' width='350' />
  <img src='/pictures/ebd3_ipsc/scatter_ipsc_ebd31p.png' width='350' />
</p>

### 3D PCAs
<p float="left">
    <img src='/pictures/3dPCA(1).png' width='280' />
    <img src='/pictures/3dPCA_fullebvsipsc(1).png' width='500'>
</p>

### PCA Loading Score Heatmaps
<p float='left'
     <img src='/pictures/loadingscatter_fibro_ipsc.png' width='300' />
     <img src='/pictures/loadingscatter_full_ipsc_eb.png' width='300'>
</p>

### Heatmaps
<p float="left">
   <img src='/pictures/heatPCAload_fibro_ipsc.png' width='250' />
   <img src='/pictures/heatPCAload_full_ipsc_eb.png' width='250'>
</p>

##### Using viridis color package
<p float='left'>
   <img src='/pictures/pcaheat_fibro_ipsc_viridis.png' width='250' />
   <img src='/pictures/pcaheat_fibro_ipsc_magma.png' width='250'>
</p>

### UHCs
<p float="left">
   <img src='/pictures/uhc_fibro_ipsc.png' width='350' />
   <img src='/pictures/uhc_full_ipsc_eb.png' width='350' />
</p>
