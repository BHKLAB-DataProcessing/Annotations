library(tidyr)
library(sjmisc)
library(readxl)
library(CoreGx)
library(data.table)
library(BiocParallel)
library(stringr)
library(stringi)

setwd('../rawdata')
all_data <- read_excel("journal.pone.0140310.s012.XLSX",sheet=2)

modify_names <- function(data_table){
    colnames <- strsplit(names(data_table),"[...]")
    colnames <- unlist(lapply(colnames, function(z){ z[!is.na(z) & z != "" & is.na(as.numeric(z))]}))
    names(data_table) <- colnames
    return(data_table)
}


check_drug_w_column <- function(drug,annotation_col){
    annotation_dt_column <- annotation_col
    bad_chars <- c("-","/","_"," ")
    new_name <- ""
    if(tolower(drug)%in% tolower(annotation_dt_column)){
        new_name <-  annotation_dt_column[(tolower(annotation_dt_column)==tolower(drug))%in% TRUE]
        return(c(tolower(annotation_dt_column)==tolower(drug),new_name))
        }
    for (i in seq(1,length(bad_chars))){
        char <- bad_chars[i]
        if(str_replace_all(tolower(drug),char,"") %in% str_replace_all(tolower(annotation_dt_column),char,"")){
            new_name <- annotation_dt_column[(str_replace_all(tolower(annotation_dt_column),char,"") ==str_replace_all(tolower(drug),char,""))%in% TRUE]
            return(c(str_replace_all(tolower(annotation_dt_column),char,"") ==str_replace_all(tolower(drug),char,""),new_name))
            }
        }

    return(c(tolower(annotation_dt_column)==tolower(d1),new_name))
    }

check_drug<- function(drug,dt,dt_names){
    
    for (name in dt_names){
        results <-  check_drug_w_column(drug,dt[[name]])
        new_name <- tail(results,n=1)
        print(c(name,new_name))
        bools <- results[1:length(results)-1] %in% TRUE
        if(new_name != "" && !is.na(new_name)){
            return(dt$unique.drugid[bools])
        }
    }
    return("")
}

check_drug_wrapper<- function(drug,dt,dt_names){
    if(grepl("/",drug)){
        d1 <- strsplit(drug,"/")[[1]][1]
        d2 <-  strsplit(drug,"/")[[1]][2]

        new_combined_name <- check_drug(drug,dt,dt_names)
        new_d1_name <-  check_drug(d1,dt,dt_names)
        new_d2_name <-  check_drug(d2,dt,dt_names)
        if(new_combined_name!=""){return(new_combined_name)}
        else if(new_d1_name != ""){return(new_d1_name)}
        else if(new_d2_name !=""){return(new_d2_name)}
        else{return("")}
    }
    return(check_drug(drug,dt,dt_names))
}


viability_experiments <- modify_names(all_data[2:40])
viability_bliss <- modify_names(all_data[43:81])
cparp_experiments <- modify_names(all_data[84:122])
cparp_bliss <- modify_names(all_data[125:163])
drug_data <- read_excel("journal.pone.0140310.s010.XLSX",sheet=1)

drugs <- drug_data$Drug
annotPath <- "../rawdata"
drugAnnotations_DT <- fread(file.path(annotPath, 'drugs_with_ids.csv'))
`%notin%` <- Negate(`%in%`)
id_columns <- names(drugAnnotations_DT)[1:21]
new_drug_names <- lapply(drugs,check_drug_wrapper,dt=drugAnnotations_DT,dt_names=id_columns)
bad_drugs <- drugs[new_drug_names==""]


d1s <- unique(viability_experiments['Drug 1'])
d2s <- unique(viability_experiments['Drug 2'])
dataset_drugs <- sort(as.vector(union(unlist(d1s),unlist(d2s))))
annotation_drugs <- sort(as.vector(drug_data['Drug'][[1]]))
#Only one drug is not matching, so I just change it manually here
drug_data[drug_data["Drug"]=="Vismodegib/GDC0449","Drug"] = "vismodegib/GDC0449"
#One Drug Combo is misnamed, fixed here
viability_experiments[11556,'Drug Combo'] = paste0(viability_experiments[11556,'Drug 1'],paste0(" + ",viability_experiments[11556,'Drug Combo']))

drug_combos <- lapply(viability_experiments[3], function (x) strsplit(as.character(x), " + ", fixed=TRUE))[[1]]

