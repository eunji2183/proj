
#clinical data
CLC_query <- GDCquery(project = cancer_type, 
                      data.category = "Clinical", 
                      file.type = "xml")
GDCdownload(CLC_query)
clinical <- GDCprepare_clinic(CLC_query, clinical.info = "patient")
write.csv(clinical,file = paste(cancer_type,"_clinical",".csv"))
#survival data
clinical_trait <- clinical  %>%
  dplyr::select(bcr_patient_barcode,gender,vital_status,
                age_at_initial_pathologic_diagnosis,
                days_to_death,days_to_last_followup,race_list) %>%
  distinct( bcr_patient_barcode, .keep_all = TRUE)
#organize dead data
dead_patient <- clinical_trait  %>%
  dplyr::filter(vital_status == 'Dead') %>%
  dplyr::select(-days_to_last_followup) %>%
  reshape::rename(c(bcr_patient_barcode = 'Barcode',
                    gender = 'Gender',
                    vital_status = 'OS',
                    days_to_death='OS.Time',
                    race_list = 'Race',
                    age_at_initial_pathologic_diagnosis = 'Age',
                    neoplasm_histologic_grade = 'Grade' )) %>%
  mutate(OS=ifelse(OS=='Dead',1,0))%>%
  mutate(OS.Time=OS.Time/365)
#organize the alive data
alive_patient <- clinical_trait %>%
  dplyr::filter(vital_status == 'Alive') %>%
  dplyr::select(-days_to_death) %>%
  reshape::rename(c(bcr_patient_barcode = 'Barcode',
                    gender = 'Gender',
                    vital_status = 'OS',
                    days_to_last_followup='OS.Time',
                    race_list = 'Race',
                    age_at_initial_pathologic_diagnosis = 'Age',
                    neoplasm_histologic_grade = 'Grade')) %>%
  mutate(OS=ifelse(OS=='Dead',1,0))%>%
  mutate(OS.Time=OS.Time/365)
#combine clincial data
survival_data <- rbind(dead_patient,alive_patient)
save(survival_data,file="/HDD8T/eunji/proj/lg/brcasurv.RData")
