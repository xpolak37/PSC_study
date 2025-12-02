# PSC study

---------------------------------------------------------------------------------------------------
**Authors and affiliations**

<div style="font-size: larger;">
Lukas Bajer<sup>1,*</sup>, Petra Polakovičova<sup>2,11*</sup>, Marie Heczkova<sup>2</sup>, Kristian Holm<sup>4,5,6</sup> Mikal J. Hole<sup>4,5,6</sup>, Mojmir Hlavaty<sup>1</sup>, Alena Bohdanecka<sup>2</sup>, Pavel Drastich<sup>1</sup>, Filip Tichanek<sup>3</sup>, Malin H. Meyer-Myklestad<sup>7</sup>, Asle W. Medhus<sup>5,8</sup>, Dag Henrik Reikvam<sup>5,7</sup>,  Kristin K. Jørgensen<sup>4,5,9</sup>, Jan Brezina<sup>1</sup>, Peter Macinga<sup>1</sup>, Pavel Wohl<sup>1</sup>, Johannes R. Hov<sup>4,5,6,10,#</sup>, Monika Cahova<sup>2,#</sup>
</div>

<br>

<sup>*</sup> These authors have contributed equally to this work and share first authorship   
<sup>#</sup> These authors have contributed equally to this work and share last authorship   

<sup>1</sup> Institute for Clinical and Experimental Medicine, Department of Hepatogastroenterology, Prague	
<sup>2</sup> Institute for Clinical and Experimental Medicine, Center for Experimental Medicine, Prague	
<sup>3</sup> Institute for Clinical and Experimental Medicine, Department of Data Science, Prague	
<sup>4</sup> Norwegian PSC Research Center, Department of Transplantation Medicine, Oslo University Hospital, Oslo, Norway	
<sup>5</sup> Institute of Clinical Medicine, University of Oslo, Oslo, Norway	
<sup>6</sup> Research Institute of Internal Medicine, Oslo University Hospital, Oslo, Norway	 
<sup>7</sup> Department of Infectious Diseases, Division of Medicine, Oslo University Hospital, Oslo, Norway	
<sup>8</sup> Department of Gastroenterology, Division of Medicine, Oslo University Hospital, Oslo, Norway	
<sup>9</sup> Department of Gastroenterology, Akershus University Hospital, Lørenskog, Norway	
<sup>10</sup> Section of Gastroenterology, Depatrment of Transplantation Medicine, Oslo University Hospital, Oslo, Norway	
<sup>11</sup> Faculty of Science, Charles University, Prague, Czech Republic

---------------------------------------------------------------------------------------------------

**General information**

This repository provides a comprehensive report of the study **Geography-independent mucosal microbiota alterations in primary sclerosing cholangitis persist after liver transplantation**

All reported results can be reproduced using the code in this repository. Feel free to contact Petra Polakovicova by [petra.polakovicova@ikem.cz](petra.polakovicova@ikem.cz) if you have any questions about the computational part of the study.

📚 **Citation** 

If you find this code and report helpful, cite the original publication:

> TO BE ADDED

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

🧬🖥️ **Bioinformatics processing** 

> TO BE ADDED

📊📈 **Statistical analysis**

> TO BE ADDED

📑✔️ ## Results

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


