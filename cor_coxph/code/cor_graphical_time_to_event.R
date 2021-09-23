#Sys.setenv(TRIAL = "moderna_mock") # moderna_mock  janssen_pooled_real  janssen_pooled_mock  janssen_na_mock
#Sys.setenv(VERBOSE = 1) 
renv::activate(project = here::here(".."))    
    # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
    if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))
source(here::here("..", "_common.R"))
#-----------------------------------------------


library(kyotil) # p.adj.perm, getFormattedSummary

myprint(study_name)
myprint(verbose)

# COR defines the analysis to be done, e.g. D29, D57, D29start1
Args <- commandArgs(trailingOnly=TRUE)
if (length(Args)==0) Args=c(COR="D29D57") 
COR=Args[1]; myprint(COR)

# COR has a set of analysis-specific parameters defined in the config file
config.cor <- config::get(config = COR)
tpeak=as.integer(paste0(config.cor$tpeak))
tpeaklag=as.integer(paste0(config.cor$tpeaklag))
tinterm=as.integer(paste0(config.cor$tinterm))
myprint(tpeak, tpeaklag, tinterm)
if (length(tpeak)==0 | length(tpeaklag)==0 | length(tinterm)==0) stop("config "%.%COR%.%" misses some fields")


# path for figures and tables etc
save.results.to = here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(here::here("output"), "/", attr(config,"config"));
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR,"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))
    


dat.mock <- read.csv(here::here("..", "data_clean", data_name))

print(paste0("reading data from ",data_name))



mypdf (mfrow=c(1,3), file=paste0(save.results.to, "barplot_mixed"))     
    tmp.1=table(subset(dat.mock, ph1.intercurrent.cases       & Trt==1 & Bserostatus==0, "EventTimePrimaryD"%.%tinterm, drop=T))
    tmp.1=table(subset(dat.mock, ph1.intercurrent.cases       & Trt==1 & Bserostatus==0, "EventTimePrimaryD"%.%tinterm, drop=T))
    tmp.2=table(subset(dat.mock, dat.mock[["ph1.D"%.%tpeak]] & dat.mock[["EventIndPrimaryD"%.%tpeak]] & Trt==1 & Bserostatus==0, "EventTimePrimaryD"%.%tinterm, drop=T))
    tmp.3=table(subset(dat.mock, dat.mock[["ph1.D"%.%tpeak]] & dat.mock[["EventIndPrimaryD"%.%tpeak]] & Trt==1 & Bserostatus==0, "EventTimePrimaryD"%.%tpeak  , drop=T))
    
    tmp=cbinduneven(list(tmp.1, tmp.2, tmp.3))
    tmp=tmp[order(as.numeric(rownames(tmp))),]
    tmp.1=tmp[,1]; names(tmp.1)=rownames(tmp)
    tmp.2=tmp[,2]; names(tmp.2)=rownames(tmp)
    tmp.3=tmp[,3]; names(tmp.3)=rownames(tmp)
    
    barplot(tmp.1, main="D29 to COVID", xlab="Days", yaxt="n", xaxt="n"); title(line=3, main="Intercurrent Cases"); axis(2, at=0:10); axis(1, at=seq(0,200,by=10)); 
    barplot(tmp.2, main="D29 to COVID", xlab="Days", yaxt="n", xaxt="n"); title(line=3, main="Post Day 57 Cases");  axis(2, at=0:10); axis(1, at=seq(0,200,by=10)); 
    barplot(tmp.3, main="D57 to COVID", xlab="Days", yaxt="n", xaxt="n"); title(line=3, main="Post Day 57 Cases");  axis(2, at=0:10); axis(1, at=seq(0,200,by=10)); 
dev.off()  
