args=(commandArgs(TRUE))


library(dplyr)
library(edgeR)
library("genefilter")
library(readr)
library(stringr)
library(DMwR)


#path<-args[1]



projects<-list.dirs("/home/jgonzalezgom/references/TCGA_SEL") # donde estan
#projects<-list.dirs(path)
projects<-projects[2:length(projects)]


for(i in c(9)){
    
    name<-basename(projects[i])
    
    cat("loading files ",name,"\n")
 
    mirna<-read_tsv(paste0(projects[i],"/",name,".mirna.tsv"))
    
    cat("preprocessing data","\n")
    mirna<-as.data.frame(mirna)
    rownames(mirna)<-mirna$miRNA_ID
    mirna$miRNA_ID<-NULL
    counts<-read_tsv(paste0(projects[i],"/",name,".htseq_counts.tsv"))
    counts<-as.data.frame(counts)
    rownames(counts)<-counts$Ensembl_ID
    counts$Ensembl_ID<-NULL

    cat("\n","selection of patients with survival, miRNA and counts only in tumors","\n")
    headers<-colnames(counts)


    headers_samples<-paste(lapply(strsplit(paste(headers),"-"),"[", 4))
    tumor_samples<-which(as.numeric(substr(headers_samples,1,2))<10)
    peri_samples<-which(as.numeric(substr(headers_samples,1,2))>=11)
    
    RNA_tumor_names<-headers[c(tumor_samples, peri_samples)]
    intersect_selection<-intersect(RNA_tumor_names,colnames(mirna))
    intersect_selection_tumor<-intersect(intersect_selection,headers[tumor_samples])
    intersect_selection_peri<-intersect(intersect_selection,headers[peri_samples])

    if(length(intersect_selection_peri) <= 30){

        cat("This database has a minimum of 30 patients of the minority class","\n")
        y<-paste0(name,  " has less than 30 patients in the minority class") 
        write(y,file="/home/jgonzalezgom/02_AOVIEDO/TCGA_CONTROL_GENES_PATIENTS.txt",append=TRUE)

        next
    }
    

    cat(length(intersect_selection)," samples with paired information selected", "\n")

    mirna<-mirna[,intersect_selection]
    counts<-counts[,intersect_selection]
    counts_o<-counts
 
    dir.create(paste0("/home/jgonzalezgom/02_AOVIEDO/",name))

    cat("filtering low counts in RNAseq", "\n")

    selection_TUMOR<-intersect_selection_tumor
    selection_PERI<-intersect_selection_peri


    log_intensity_threshold<-5

    numsamples<-round(length(selection_TUMOR)/2)
    f1 <- kOverA(numsamples, log_intensity_threshold)
    ffun1 <- filterfun(f1)
    whichFilter.TUMOR<-genefilter(counts_o[,which(colnames(counts_o) %in% selection_TUMOR)], ffun1)

    numsamples<-round(length(selection_PERI)/2)
    f1 <- kOverA(numsamples, log_intensity_threshold)
    ffun1 <- filterfun(f1)
    whichFilter.PERI<-genefilter(counts_o[,which(colnames(counts_o) %in% selection_PERI)], ffun1)

    whichFilter <- (whichFilter.TUMOR | whichFilter.PERI )   

    selection_data_filter<-counts_o[whichFilter,]

    cat("DEG RNAseq", "\n")

    group<-c(rep("TUMOR", length(intersect_selection_tumor)),rep("PERI",length(intersect_selection_peri)))
    d0<-DGEList(selection_data_filter, group=group)
    d0 <- calcNormFactors(d0, norm.method="TMM", test.method="edger")
  
    # estimate dispersions
    d0 <- estimateCommonDisp(y = d0, verbose=T)
    d0 <- estimateTagwiseDisp(d0)  

    snames<-group
    design<-model.matrix(~0+snames)

    d1<-estimateGLMCommonDisp(d0, design)
    d1<- estimateGLMTrendedDisp(d1, design, method="power")
    d1<-estimateGLMTagwiseDisp(d1, design)


    fit <- glmFit(d1, design)

    #make contrasts

    PERIvsTUMOR<-glmLRT(fit, coef=ncol(fit$design), contrast=c(1,-1))
    
    PERIvsTUMOR<-topTags(PERIvsTUMOR, n=nrow(PERIvsTUMOR), adjust.method = "BH", sort.by = "PValue")
    
    PERIvsTUMOR<-PERIvsTUMOR$table

    write.table(PERIvsTUMOR, file=paste0("/home/jgonzalezgom/02_AOVIEDO/",name,"/",name,"_RNASEQ_DEG.txt"))


    GENES_DEG<-rownames(PERIvsTUMOR[which(PERIvsTUMOR$FDR < 0.05),])

    normalized_counts<-cpm(d1 ,normalized.lib.sizes = TRUE, log=T)
    normalized_counts_final<-normalized_counts[which(rownames(normalized_counts) %in% GENES_DEG),]

    cat("SAMPLING EQUILIBRADO PARA EL MODELO", "\n")

    if((length(intersect_selection_tumor)> 3*(length(intersect_selection_peri)) )){
        
        if(length(intersect_selection_tumor) > 500){
 
           patients<-sample(intersect_selection_tumor, round(length(intersect_selection_tumor)/2))
           patients_sel<-c(patients, intersect_selection_peri)

           class<-c(rep("TUMOR", round(length(intersect_selection_tumor)/2)),rep("PERI",length(intersect_selection_peri)))
        
        }else{

            patients<-sample(intersect_selection_tumor, length(intersect_selection_tumor))
            patients_sel<-c(patients, intersect_selection_peri)
            class<-c(rep("TUMOR", length(intersect_selection_tumor)),rep("PERI",length(intersect_selection_peri)))

        }

        normalized_counts_final_sel<-normalized_counts_final[,patients_sel]
        mirna_sel<-mirna[,patients_sel]

        ALL<-rbind(normalized_counts_final_sel, mirna_sel)

        normalized_counts_final_sel_t<-as.data.frame(t(ALL))
        #mirna_sel_t<-as.data.frame(t(mirna_sel))

        normalized_counts_final_sel_t$class<-as.factor(paste(class))
        #mirna_sel_t$class<-as.factor(paste(class))

        if(length(intersect_selection_peri < 50)){

            RNAseq_smote <- SMOTE(class ~ ., data = normalized_counts_final_sel_t, perc.over = 300, k = 5)
           # mirna_smote <- SMOTE(class ~ ., data = mirna_sel_t, perc.over = 300, k = 5)


        }else{

            RNAseq_smote <- SMOTE(class ~ ., data = normalized_counts_final_sel_t, perc.over = 200, k = 5)
           # mirna_smote <- SMOTE(class ~ ., data = mirna_sel_t, perc.over = 200, k = 5)

        }
        
        
        normalized_counts_final_sel_s<-RNAseq_smote
        #normalized_counts_final_sel_s$class<-NULL
        normalized_counts_final_sel_s<-t(normalized_counts_final_sel_s)

        #mirna_sel_s<-mirna_smote
        #mirna_sel_s$class<-NULL
        #mirna_sel_s<-t(mirna_sel_s)

    }else if(length(intersect_selection_tumor)> 2*(length(intersect_selection_peri))& (length(intersect_selection_tumor)< 3*(length(intersect_selection_peri)))){
        
        if(length(intersect_selection_tumor) > 500){
          
          patients<-sample(intersect_selection_tumor, round(length(intersect_selection_tumor)/2))
          patients_sel<-c(patients, intersect_selection_peri)

        class<-c(rep("TUMOR", round(length(intersect_selection_tumor)/2)),rep("PERI",length(intersect_selection_peri)))
        
        }else{
         
            patients<-sample(intersect_selection_tumor, length(intersect_selection_tumor))
            patients_sel<-c(patients, intersect_selection_peri)
            class<-c(rep("TUMOR", length(intersect_selection_tumor)),rep("PERI",length(intersect_selection_peri)))

        }
        normalized_counts_final_sel<-normalized_counts_final[,patients_sel]
        mirna_sel<-mirna[,patients_sel]

        ALL<-rbind(normalized_counts_final_sel, mirna_sel)

        normalized_counts_final_sel_t<-as.data.frame(t(ALL))
        #mirna_sel_t<-as.data.frame(t(mirna_sel))

        normalized_counts_final_sel_t$class<-as.factor(paste(class))
        #mirna_sel_t$class<-as.factor(paste(class))

        if(length(intersect_selection_peri < 50)){

            RNAseq_smote <- SMOTE(class ~ ., data = normalized_counts_final_sel_t, perc.over = 300, k = 5)
           # mirna_smote <- SMOTE(class ~ ., data = mirna_sel_t, perc.over = 300, k = 5)


        }else{

            RNAseq_smote <- SMOTE(class ~ ., data = normalized_counts_final_sel_t, perc.over = 200, k = 5)
           # mirna_smote <- SMOTE(class ~ ., data = mirna_sel_t, perc.over = 200, k = 5)

        }
        
        
        normalized_counts_final_sel_s<-RNAseq_smote
        #normalized_counts_final_sel_s$class<-NULL
        normalized_counts_final_sel_s<-t(normalized_counts_final_sel_s)
        
    }else if((length(intersect_selection_tumor)> length(intersect_selection_peri)) & (length(intersect_selection_tumor)< 2*(length(intersect_selection_peri)) )){
        
        patients_sel<-c(intersect_selection_tumor, intersect_selection_peri)
        normalized_counts_final_sel<-normalized_counts_final[,patients_sel]
        mirna_sel<-mirna[,patients_sel]

        class<-c(rep("TUMOR", length(intersect_selection_tumor)),rep("PERI",length(intersect_selection_peri)))
        
        #normalized_counts_final_sel_t<-as.data.frame(t(normalized_counts_final_sel))
        #mirna_sel_t<-as.data.frame(t(mirna_sel))

        ALL<-rbind(normalized_counts_final_sel, mirna_sel)

        normalized_counts_final_sel_t<-as.data.frame(t(ALL))
        #mirna_sel_t<-as.data.frame(t(mirna_sel))

        normalized_counts_final_sel_t$class<-as.factor(paste(class))
        #mirna_sel_t$class<-as.factor(paste(class))

        if(length(intersect_selection_peri < 50)){

            RNAseq_smote <- SMOTE(class ~ ., data = normalized_counts_final_sel_t, perc.over = 100, k = 5)
           # mirna_smote <- SMOTE(class ~ ., data = mirna_sel_t, perc.over = 300, k = 5)


        }else{

            RNAseq_smote <- SMOTE(class ~ ., data = normalized_counts_final_sel_t, perc.over = 100, k = 5)
           # mirna_smote <- SMOTE(class ~ ., data = mirna_sel_t, perc.over = 200, k = 5)

        }
        
        
        normalized_counts_final_sel_s<-RNAseq_smote
        #normalized_counts_final_sel_s$class<-NULL
        normalized_counts_final_sel_s<-t(normalized_counts_final_sel_s)
    
    }

    tabla<-table(RNAseq_smote$class)

    y<-paste0(name,  " GENES DEG= ",length(GENES_DEG)," ;PREVIOUS PAT T= " ,length(intersect_selection_tumor)," ;PREVIOUS PAT P= ",length(intersect_selection_peri) ," ;AFTER PAT T= ",tabla["TUMOR"]," ;AFTER PAT P= ",tabla["PERI"], ";FINAL= ",dim(normalized_counts_final_sel_s)[2])
    write(y,file="/home/jgonzalezgom/02_AOVIEDO/TCGA_CONTROL_GENES_PATIENTS.txt",append=TRUE)

    normalized_counts_final_t<-t(normalized_counts_final_sel_s)
    colnames(normalized_counts_final_t)<-sapply(strsplit(colnames(normalized_counts_final_t),"\\."),"[",1)
    colnames(normalized_counts_final_t)<-tolower(colnames(normalized_counts_final_t))
    names(normalized_counts_final_t) <- sapply(names(normalized_counts_final_t), function(x) {
      if(grepl("^ensg", x)) {
        return(toupper(x))
     } else {
       return(x)
      }
    })


    write.table(cbind("patient"=rownames(normalized_counts_final_t),normalized_counts_final_t), file=paste0("/home/jgonzalezgom/02_AOVIEDO/",name,"/",name,"_RNA_AND_mirRNA_Counts.txt"),sep="\t")

    # Remover las variables
    variables_remover <- setdiff(ls(), c("i","projects"))
    rm(list = variables_remover)
    gc()

   # i=i+1

}



