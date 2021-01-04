# taser_metabolomics
Investigating metabolomic changes in RA patients treated with MTX over 18 months to direct predictive biomarker discovery

The previous repository- taser_initial - contains the R scripts from the initial exploratory work looking into the RA taser dataset from JC. 
This provided some insight into the data but the method was altered and ultimately cleaned up to be more reproducible (hopefully) when additional datasets are received. 

Predictive biomarker discovery will be the ultimate goal of my project, with this early part providing an exploratory approach, where the metabolites themselves will only be putatively identified and need further assessment (fragmentation possibly) to confirm their identities due to the limitations of the LC-MS approach.
Patients samples include those from baseline (prior to treatment) and 18 months following treatment initiation. In the taser trial, previous work looked at the 3 month changes but this current work will focus on the changes and associations over the 18 month period. This is due to the inclusion of the erosion imaging data alongside the disease activity data, the former not having been taken at 3 months, so metabolomic associations with erosion changes could not be achieved for this. 

Goals for this section of the analysis include the assessment of:
  - metabolites whose levels at baseline which are associated with clinical outcomes based on DAS44
  - metabolites whose levels change over the 18 months of treatment which are associated with the clinical outcomes
  - metabolites whose levels at baseline which are associated with erosion imaging (clinical outcomes from MRI and clinical outcomes from Xray)
  - metabolites whose levels change over 18 months which are associated with the erosion imaging outcomes
  
Data analysis will be done using R with the identification of metabolites which contribute towards the clinical outcomes being carried out using a machine learning approach, followed by correlating metabolite levels with the clinical outcomes and a differential analysis. These will be done to assess which metabolites are important in contributing towards the DAS44 and erosion imaging clinical outcomes (involving feature selection using ML), assessing which metabolites are associated with changes in the disease activity through the correlations and finally looking to see which metabolites change in their levels across baseline for respective clinical outcome subgroups, and which metabolites change over the 18 month period. This latter point will help in the discovery of predictive biomarkers and mechanism of action biomarkers, respectively. 

Files included:
20201127_resp_ints.csv - rows are samples and columns 14:1472 are peaks (peak ID shown only). Columns 1:13 show patient metadata and clinical disease measurements. Samples include A (baseline) and F (18 months post-treatment initiation). 
20201105_resp_diff.csv- rows are samples and columns 13:1471 are peaks. Columns 1:12 show clinical outcomes and categorised clinical outcomes for DAS44 changes. Only the changes in the sample measurements across baseline and 18 months are shown in this dataframe.
20200117_Taser_PeakIDs.csv- shows peak IDs from LC-MS analysis and annotations indicating possible metabolite names based on m/z ratio. These give basis for the putative identification of the metabolites in later analyses.
20200427_Taser_PeakMetadata.csv- provides a greater degree of information on each of the peaks from the LC-MS analysis. This includes the putative identification of the metabolites with KEGG and HMDB references. The m/z ratio, RT and metabolite type are listed for each. 
20200318_Taser_SampleSheet.csv- samples are rows (A and F) with the clinical outcome measurements for the two points in the study
20190713_Taser_PatientMetadata.csv- samples are rows. Patient metadata included, including age, sex and smoking status. Some missing data was shown.
20200430_Taser_NEG_PeakIntensities.csv- the adjusted negative ion mode peak data for A and F samples. Data was pre-processed in 20200430_Taser_Parsing_NEG.R to remove missing data and normalise the data. It was later batch corrected to provide the  20201127_resp_ints.csv and 20201105_resp_diff.csv files

