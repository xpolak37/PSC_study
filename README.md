# PSC study



<div style="font-size: larger;">
Lukas Bajer <sup>1,2,*</sup>, Petra Polakovicova <sup>3,4,*</sup>, Marie Heczkova <sup>3</sup>, Kristian Holm <sup>5,6,7</sup>, Mikal J Hole <sup>5,6,7</sup>, Mojmir Hlavaty <sup>1,8</sup>, Alena Bohdanecka <sup>3,8</sup>, Pavel Drastich <sup>1</sup>, Filip Tichanek <sup>9</sup>, Malin H Meyer-Myklestad <sup>10,11</sup>, Asle W Medhus <sup>6,12</sup>, Dag Henrik Reikvam <sup>6,10</sup>, Kristin K Jørgensen <sup>5,6,13</sup>, Jan Brezina <sup>1</sup>, Peter Macinga <sup>1</sup>, Pavel Wohl <sup>1</sup>, Ondrej Fabian <sup>14</sup>, Johannes R Hov <sup>5,6,7,15,#</sup>, Monika Cahova <sup>3,#</sup>
</div>

<br>

<sup>*</sup> These authors have contributed equally to this work and share first authorship   
<sup>#</sup> These authors have contributed equally to this work and share last authorship   

<sup>1</sup> Institute for Clinical and Experimental Medicine, Department of Hepatogastroenterology, Prague, CR, Czech Republic
<sup>2</sup> Department of Internal Medicine, 2nd Faculty of Medicine, Charles University, Prague, CR, Czech Republic
<sup>3</sup> Institute for Clinical and Experimental Medicine, Center for Experimental Medicine, Prague, CR, Czech Republic
<sup>4</sup> Faculty of Science, Charles University, Prague, CR, Czech Republic
<sup>5</sup> Norwegian PSC Research Center, Department of Transplantation Medicine, Oslo University Hospital, Oslo, Norway
<sup>6</sup> Institute of Clinical Medicine, University of Oslo, Oslo, Norway
<sup>7</sup> Research Institute of Internal Medicine, Oslo University Hospital, Oslo, Norway
<sup>8</sup> First Faculty of Medicine, Charles University, Prague, CR, Czech Republic
<sup>9</sup> Institute for Clinical and Experimental Medicine, Department of Data Science, Prague, CR, Czech Republic
<sup>10</sup> Department of Infectious Diseases, Division of Medicine, Oslo University Hospital, Oslo, Norway
<sup>11</sup> Department of Microbiology, Division of Laboratory Medicine, Oslo University Hospital, Oslo, Norway
<sup>12</sup> Department of Gastroenterology, Division of Medicine, Oslo University Hospital, Oslo, Norway
<sup>13</sup> Department of Gastroenterology, Akershus University Hospital, Lørenskog, Norway
<sup>14</sup> Institute for Clinical and Experimental Medicine, Department of Pathology, Prague, CR, Czech Republic
<sup>15</sup> Section of Gastroenterology, Department of Transplantation Medicine, Oslo University Hospital, Oslo, Norway

---------------------------------------------------------------------------------------------------

**General information**

This repository provides a comprehensive report of the study **Geography-independent mucosal microbiota alterations in primary sclerosing cholangitis persist after liver transplantation**

All reported results can be reproduced using the code in this repository. Feel free to contact Petra Polakovicova by [petra.polakovicova@ikem.cz](petra.polakovicova@ikem.cz) if you have any questions about the computational part of the study.

📚 **Citation** 

If you find this code and report helpful, cite the original publication:

Bajer, L., Polakovicova, P., Heczkova, M., Holm, K., Hole, M. J., Hlavaty, M., Bohdanecka, A., Drastich, P., Tichanek, F., Meyer-Myklestad, M. H., Medhus, A. W., Reikvam, D. H., Jørgensen, K. K., Brezina, J., Macinga, P., Wohl, P., Fabian, O., Hov, J. R., & Cahova, M. (2025). Geography-independent mucosal microbiota alterations in primary sclerosing cholangitis persist after liver transplantation. JHEP Reports, 8(4), 101716. [https://doi.org/10.1016/j.jhepr.2025.101716](https://doi.org/10.1016/j.jhepr.2025.101716)

💾 **Data Availability**
- Czech cohort: 
> SRA PRJNA1250244
- Norwegian cohort:
> 10.1002/hep.32773

----------------------------------------------------------------------------------------------------

## Report info

This project analyses biopsy samples from two cohorts (Czech and Norwegian), namely data from amplicon sequencing of the 16S rRNA gene (region V3-V4) and related clinical parameters relevant to PSC disease. 370 subjects are included, from whom a total of 1083 samples were collected.

📁 **Project Structure**

Below is an overview of the folder structure:

- **scripts/** - Source code for bioinformatics processing of raw sequencing data to ASV taxonomic tables (mainly bash scripts)

- **analysis/**  
  - `scripts/`
	- `merged_sites/` – scripts for analyzing the data of *terminal ileum* and *colon* sites
		- `main_analysis/` – main part of analysis, whose results are directly reported in the publication
		- `supplementary_analysis/` – additional part of the analysis, where different methods or metrics were used
	- `split_sites/`  – scripts for analyzing the data of *terminal ileum*, *left_colon* and *right_colon* sites

  - `results/` – results generated directly via provided scripts


## Methodology

For detailed methodology, see the original publication. 

🖥️ **Bioinformatics processing** 

The Illumina paired-end reads were quality-checked, and after preprocessing, amplicon sequence variants (ASVs) were generated with Deblur in QIIME2 (2024.2) after trimming reads to 400 bp. The amplicon-region-specific Naive Bayes classifier was trained based on the SILVA Ref NR 99 database v 138.1 [Quast] via RESCRIPt QIIME 2 plugin. 

📊 **Statistical analysis**

All statistical analyses were performed in R (v4.3.1) on merged Czech and Norwegian microbiome datasets at the ASV level, which were subsequently analyzed separately for terminal ileum and colon samples. The analysis proceeded in two stages: first comparing post-LTx, pre-LTx, and healthy controls (HC), and second subdividing post-LTx into rPSC and non-rPSC groups, alongside HC, with an additional IBD versus non-IBD comparison within PSC patients. Alpha diversity (ASV richness and Shannon index) was calculated on rarefied ASV-level data, whereas beta diversity, differential abundance testing, and classification analyses were conducted at the genus level following filtering of low-depth samples and low-prevalence taxa. Beta diversity was assessed using robust Aitchison distance with permutational ANOVA, while differential abundance was determined using both linDA and MaAsLin2, retaining only concordant significant taxa across cohorts. Predictive modeling employed elastic net alongside random forest, gradient boosting, and k-nearest neighbors, with hyperparameter tuning, bootstrapped validation, and optimism-corrected AUC evaluation to assess classification performance and control overfitting. A microbial dysbiosis index (MDI) was calculated as the log ratio (implemented as a difference in clr-transformed data) between taxa increased and decreased in PSC, and its association with clinical parameters was evaluated using iterative Spearman correlations, with significance defined by consistent results across repeated subsampling.


## Results 

The code with reported results can be found:
- [Q1_analysis](https://xpolak37.github.io/PSC_study/analysis/scripts/merged_sites/main_analysis/psc_study_q1.html): reports the comparison of pre-LTx, post-LTx, and HC
- [Q2_analysis](https://xpolak37.github.io/PSC_study/analysis/scripts/merged_sites/main_analysis/psc_study_q2.html): reports the comparison of rPSC vs. non-rPSC
- [Q3_analysis](https://xpolak37.github.io/PSC_study/analysis/scripts/merged_sites/main_analysis/psc_study_q3.html): reports the comparison of IBD vs. non-IBD
- [Q5_analysis](https://xpolak37.github.io/PSC_study/analysis/scripts/merged_sites/main_analysis/psc_study_q5.html): reports the comparison between patients with low and high fecal calprotectin values
- [MDI_analysis](https://xpolak37.github.io/PSC_study/analysis/scripts/merged_sites/main_analysis/psc_study_mdi.html): calculation and exploration of the microbiome dysbiosis index
- [ALD_analysis](https://xpolak37.github.io/PSC_study/analysis/scripts/merged_sites/main_analysis/psc_study_ald.html): reports the comparison between Czech ALD and PSC samples
- [clinical_analysis](https://xpolak37.github.io/PSC_study/analysis/scripts/merged_sites/main_analysis/psc_study_clinical.html): explores associations between the microbiome and clinical features
- [ML_overfitting_check](https://xpolak37.github.io/PSC_study/analysis/scripts/merged_sites/main_analysis/psc_study_ml_overfitting_check.html): reports all models used in this study and shows the results with reshuffled labels, confirming that the original datasets perform well

---------------------------------------------------------------------------------------------------

## Acknowledgment

This study was supported by grant MH CR no. NU21J-06-00027, by the project National Institute for Research of Metabolic and Cardiovascular Diseases (Programme EXCELES, Project No. LX22NPO5104) - Funded by the European Union - Next Generation EU and by MH CR –DRO (Institute for Clinical and Experimental medicine –IKEM, IN 00023001”). JRH was funded by the European Research Council (StopAutoimmunity, no. 802544).


