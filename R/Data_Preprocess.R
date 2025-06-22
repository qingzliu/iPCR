library(readxl)
library(dplyr)
library(readr)
library(ggplot2)
data_dir = "/Users/qingzliu/Documents/Thesis1/Data"
cell_line_exp_file = "CellignerData/CCLE_mat22.csv"
cell_line_GDSCanno_file = "TRANSACTData/model_list_20191104.csv"
cell_line_CCLEanno_file = "CellignerData/sample_info.csv"
cell_line_response_file1 = "TRANSACTData/GDSC1_fitted_dose_response_25Feb20.xlsx"
cell_line_response_file2 = "TRANSACTData/GDSC2_fitted_dose_response_25Feb20.xlsx"
mini_cancer_genes_file = "TRANSACTData/mini_cancer_lookup_genes.csv"

#### CCLE ####

## CCLE Gene Expression 
CCLE_exp =  readr::read_csv(file.path(data_dir, cell_line_exp_file)) %>% 
  as.data.frame() %>%
  tibble::column_to_rownames('...1') %>%
  as.matrix()
colnames(CCLE_exp) = stringr::str_match(colnames(CCLE_exp), '\\((.+)\\)')[,2]

## CCLE and GDSC cell line annotation
cell_line_GDSCanno = readr::read_csv(file.path(data_dir, cell_line_GDSCanno_file))
cell_line_CCLEanno = readr::read_csv(file.path(data_dir, cell_line_CCLEanno_file))
cell_line_CCLEanno = cell_line_CCLEanno[match(rownames(CCLE_exp), cell_line_CCLEanno$DepMap_ID), ]
cell_line_CCLEanno = as.data.frame(cell_line_CCLEanno)
rownames(cell_line_CCLEanno) = cell_line_CCLEanno$DepMap_ID

SIDM_ACH_match = cell_line_GDSCanno[,c("model_id", "BROAD_ID")]
SIDM_ACH_match_complete = SIDM_ACH_match[complete.cases(SIDM_ACH_match), ]
SIDM_ACH_match_complete = as.data.frame(SIDM_ACH_match_complete)

CCLE_SIDM_ACH_match = SIDM_ACH_match_complete[match(rownames(CCLE_exp), as.vector(SIDM_ACH_match_complete[,2])), ]
CCLE_SIDM_ACH_match = CCLE_SIDM_ACH_match[complete.cases(CCLE_SIDM_ACH_match), ]

## mini_cancer_genes expression & annotation
mini_cancer_genes <-  readr::read_csv(file.path(data_dir, mini_cancer_genes_file))
mini_cancer_genes = as.data.frame(mini_cancer_genes)
CCLE_exp = CCLE_exp[ ,intersect(colnames(CCLE_exp), mini_cancer_genes$ENSEMBL)]

hgnc_file = "CellignerData/hgnc_complete_set_7.24.2018.txt"
hgnc.complete.set = data.table::fread(file.path(data_dir, hgnc_file)) %>% as.data.frame()

mini_cancer_genes_anno = hgnc.complete.set[match(colnames(CCLE_exp), hgnc.complete.set$ensembl_gene_id), ]
mini_cancer_genes_anno = mini_cancer_genes_anno[match(colnames(CCLE_exp), mini_cancer_genes_anno$ensembl_gene_id), ]
rownames(mini_cancer_genes_anno) = mini_cancer_genes_anno$ensembl_gene_id
colnames(CCLE_exp) = mini_cancer_genes_anno$symbol

#### Pan-Cancer TCGA ####

## Pan-Cancer TCGA annotation download
panTCGA_anno_file = "XenaData/Survival_SupplementalTable_S1_20171025_xena_sp.txt"
panTCGA_anno = data.table::fread(file.path(data_dir, panTCGA_anno_file)) %>% as.data.frame()
panTCGA_anno = panTCGA_anno[!duplicated(panTCGA_anno$'_PATIENT'), ]
rownames(panTCGA_anno) = panTCGA_anno$'_PATIENT'

## Clean panTCGA, Treehouse, and target annotation data
Tumor_exp_file = "CellignerData/TCGA_mat_v11.tsv"
Tumor_anno_file = "CellignerData/clinical_TumorCompendium_v11_PolyA_2020-04-09.tsv"
Tumor_response_file = "TRANSACTData/response.csv"
Tumor_exp <-  readr::read_tsv(file.path(data_dir, Tumor_exp_file)) %>% 
  as.data.frame() %>%
  tibble::column_to_rownames('Gene') %>%
  as.matrix() %>% 
  t()
Tumor_anno = readr::read_tsv(file.path(data_dir, Tumor_anno_file)) %>% as.data.frame()
Tumor_exp = Tumor_exp[match(Tumor_anno$th_sampleid, rownames(Tumor_exp)), match(colnames(CCLE_exp), colnames(Tumor_exp))]
Tumor_anno$Cancer_Type = NA
Tumor_anno = Tumor_anno %>% relocate(Cancer_Type, .after = th_sampleid)

## Clean panTCGA Annotation
Tumor_anno$Cancer_Type = panTCGA_anno$`cancer type abbreviation`[match(Tumor_anno$site_donor_id, panTCGA_anno$'_PATIENT')] # add abbreviated cancer types
Tumor_anno[Tumor_anno$th_sampleid %in% c("TCGA-28-2510-01", "TCGA-2G-AAKO-01", "TCGA-2G-AAKO-05" , "TCGA-2G-AALF-01", "TCGA-2G-AALG-01", "TCGA-2G-AALN-01", "TCGA-2G-AALO-01", "TCGA-2G-AALQ-01", "TCGA-2G-AALR-01", "TCGA-2G-AALS-01", "TCGA-2G-AALT-01", "TCGA-2G-AALW-01", "TCGA-2G-AALX-01", "TCGA-2G-AALY-01", "TCGA-2G-AALZ-01","TCGA-2G-AAM2-01", "TCGA-2G-AAM3-01", "TCGA-2G-AAM4-01", "TCGA-5M-AAT5-01", "TCGA-5M-AATA-01", "TCGA-BH-A0B2-01", "TCGA-F5-6810-01", "TCGA-R8-A6YH-01"), 2] = c("GBM", rep("TGCT",17), "COAD", "COAD", "BRCA", "READ", "LGG")
Tumor_anno_raw = Tumor_anno
panTCGA_sampleid = Tumor_anno$th_sampleid[2942:12747]
table(substr(panTCGA_sampleid, 14, 15)) # examine the frequency of samples types
Tumor_anno = Tumor_anno[!Tumor_anno$th_sampleid %in% panTCGA_sampleid[which(substr(panTCGA_sampleid, 14, 15) == "02")][substr(panTCGA_sampleid[which(substr(panTCGA_sampleid, 14, 15) == "02")-1], 1, 12)==substr(panTCGA_sampleid[which(substr(panTCGA_sampleid, 14, 15) == "02")], 1, 12)], ] # remove "02" (Recurrent Solid Tumor) sample when "01" (Primary Solid Tumor) exists for a certain patient
Tumor_anno = Tumor_anno[!Tumor_anno$th_sampleid %in% panTCGA_sampleid[which(substr(panTCGA_sampleid, 14, 15) == "05" | substr(panTCGA_sampleid, 14, 15) == "07")], ] # remove "05" (Additional - New Primary) and "07" (Additional Metastatic)
panTCGA_sampleid = Tumor_anno$th_sampleid[2942:12747]
table(substr(panTCGA_sampleid, 14, 15)) # re-examine the frequency of samples types after removal
panTCGA_anno = Tumor_anno[2942:12702, ]
panTCGA_duplicated_anno = panTCGA_anno[panTCGA_anno$site_donor_id %in% panTCGA_anno$site_donor_id[duplicated(panTCGA_anno$site_donor_id)],  ] # 28 samples that have both "01" (Primary Solid Tumor) and "06" (Metastatic)

## Clean Treehouse & Target

# Among Treehouse and Target samples, only select those disease types (listed below) that appears in CCLE:
# 2 acute leukemia; 4 acute lymphoblastic leukemia; 5 acute megakaryoblastic leukemia
# 6 acute myeloid leukemia; 11 alveolar rhabdomyosarcoma; 14 atypical teratoid/rhabdoid tumor
# 15 cholangiocarcinoma; 17 chronic myelogenous leukemia (S02), acute lymphoblastic leukemia (S01) 
# 19 colon adenocarcinoma; 23 embryonal rhabdomyosarcoma; 25 endometrial stromal sarcoma
# 28 epithelioid sarcoma; 29 Ewing sarcoma; 36 glioblastoma multiforme; 37 glioma; 39 hepatoblastoma
# 40 hepatocellular carcinoma; 45 kidney clear cell carcinoma; 46 leukemia; 48 lung adenocarcinoma
# 50 malignant peripheral nerve sheath tumor; 51 medulloblastoma; 52 melanoma; 54 meningioma
# 61 neuroblastoma; 66 osteosarcoma; 67 ovarian serous cystadenocarcinoma; 73 retinoblastoma
# 74 rhabdoid tumor; 75 rhabdomyosarcoma; 77 sarcoma; 82 synovial sarcoma ; 83 teratoma
# 86 thyroid carcinoma; 88 undifferentiated pleomorphic sarcoma
TH_Target_anno = Tumor_anno[1:2941, ]
TH_Target_anno = TH_Target_anno[TH_Target_anno$disease %in% sort(unique(Tumor_anno[1:2941, ]$disease))[c(2,4,5,6,11,14,15,17,19,23,25,28,29,36,37,39,40,45,46,48,50,51,52,54,61,66,67,73,74,75,77,82,83,86,88)], ] # Among Treehouse and Target samples, only select those disease types that appears in CCLE 
TH_Target_anno = rbind(TH_Target_anno[1:1751,][!duplicated(substr(TH_Target_anno$th_sampleid[1:1751], 1, nchar(TH_Target_anno$th_sampleid[1:1751]) - 4)), ], TH_Target_anno[1752:2467,]) # delete duplicated samples in Treehouse 
Tumor_anno = rbind(TH_Target_anno, panTCGA_anno)
Tumor_exp = Tumor_exp[Tumor_anno$th_sampleid, ]

#### Normalize Tumor and Cell line Expression Data

## Mean and Scale Correction
CCLE_exp_norm = apply(CCLE_exp, 2, function(x){(x-mean(x))/sd(x)})
Tumor_exp_norm = apply(Tumor_exp, 2, function(x){(x-mean(x))/sd(x)})

#### Drug Response ####

## Download CCLE and panTumor Drug Response 
GDSC1_response_25Feb20 = read_excel(file.path(data_dir, cell_line_response_file1))
GDSC2_response_25Feb20 = read_excel(file.path(data_dir, cell_line_response_file2))
GDSC_response = rbind(GDSC1_response_25Feb20, GDSC2_response_25Feb20)
GDSC_response = GDSC_response[GDSC_response$SANGER_MODEL_ID%in%CCLE_SIDM_ACH_match$model_id, ]
TCGA_response_file = "TRANSACTData/response.csv"
TCGA_response = read_csv(file = file.path(data_dir, TCGA_response_file),
                         col_types = cols(...1 = col_skip()))

## Functions for Reading CCLE and panTumor Drug Response Data
read_GDSC_drug_data = function(Drug_name) {
  GDSC_response_drug = GDSC_response%>%filter(DRUG_NAME==Drug_name)%>%group_by(SANGER_MODEL_ID)%>%summarize(mean_AUC = mean(AUC, na.rm = TRUE))
  GDSC_response_drug = as.data.frame(GDSC_response_drug)
  rownames(GDSC_response_drug) = CCLE_SIDM_ACH_match$BROAD_ID[match(GDSC_response_drug$SANGER_MODEL_ID, CCLE_SIDM_ACH_match$model_id)]
  GDSC_response_drug = GDSC_response_drug[order(rownames(GDSC_response_drug)), ]
  return(GDSC_response_drug)
}
read_TCGA_drug_data = function(Drug_name) {
  TCGA_response_drug = TCGA_response%>%filter(drug_name == Drug_name)
  TCGA_response_drug = as.data.frame(TCGA_response_drug)
  TCGA_response_drug = TCGA_response_drug[TCGA_response_drug$bcr_patient_barcode%in%Tumor_anno$site_donor_id, ]
  TCGA_response_drug = TCGA_response_drug[!duplicated(TCGA_response_drug$bcr_patient_barcode), ]
  rownames(TCGA_response_drug) = TCGA_response_drug$bcr_patient_barcode
  return(TCGA_response_drug)
}

