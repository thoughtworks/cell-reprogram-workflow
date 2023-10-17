computeMostSynergisticCores <- function(data=NULL, classes=NULL, infile=NULL, donorC=NULL, outDir="."){
  if(is.null(data)){ stop("data undefined") }
  if(is.null(classes)){ stop("target class undefined") } 
  if(is.null(donorC)){ stop("starting cell type undefined") }
  if(! file.exists(donorC) & ! file.exists((paste0("startpop/",donorC))) ){ stop("expression of starting cell type does not exist") }
  FDplus=1

  ## HCL
  data2 = t(data)
  data2 = data2[,colSums(data2)!=0]
  data2 <- scale(data2)
  rownames(data2) = gsub("^(.+)\\.\\d+","\\1",rownames(data2))
  cla = sort(unique(rownames(data2)))
  data2 = t(sapply(cla,function(x){
    colMeans(data2[rownames(data2) %in% x,,drop=FALSE])
  }))
  # Pearson distance
  com2 = combn(1:nrow(data2),2)
  dis = apply(com2,2,function(i){
    1 - cor(data2[i[1],], data2[i[2],])
  })
  names(dis) = apply(com2,2,function(i){ paste0(rownames(data2)[i[1]],"___vs___",rownames(data2)[i[2]]) })
  dis = as.matrix(dis)
  mdNa = do.call('rbind',strsplit(rownames(dis),"___vs___"))
  cla = rownames(data2)
  m = matrix(NA,nrow=length(cla),ncol=length(cla))
  rownames(m) = colnames(m) = cla
  naI = cbind(match(mdNa[,1],rownames(m)), match(mdNa[,2],colnames(m)))
  m[rbind(naI, cbind(naI[,2],naI[,1]))] = c(dis[,1], dis[,1])
  m = as.dist(m)
  colorLeafs <- function(x) {
    if (is.leaf(x) && attr(x, "label") %in% target.classes) {
      attr(x, "nodePar") <- list(lab.col="red", pch=NA)
      }
    return(x)
  }
  pdf(paste0(outDir,'/hc_Pearson_',geneexp,'.pdf'))
  hc1 <- hclust(m, method="average" )
  dd <- dendrapply(as.dendrogram(hc1), colorLeafs)
  plot(as.dendrogram(dd,hang=-1), cex=0.8)
  graphics.off()

    # sizeNormalization
    data = sweep(data,2,colSums(data),`/`)
    data[is.na(data)] = 0

    col.numbers = which(colnames(data) %in% classes)
    if(length(col.numbers)==0){ stop("target class does not exist") }
    # Die if target subpopulation does not have enough cells
    if(length(col.numbers) >3 ){
        a2 = as.matrix(data[,col.numbers, drop=FALSE])
    }else{
        stop("the number target cells less than 3")
    }
    a4 = log10(as.matrix(data[,-col.numbers]) + 1)

    ## Expressed TFs 
    expressed = rowSums(a2)!=0
    GEN = sort(names(which(expressed)))
    if(length(GEN) > 300){
        # Option 1: CV
        m = apply(a2,1,mean); s = apply(a2,1,sd)
        cv = (s/m)*100
        e = cv[names(cv) %in% GEN]
        GEN = names(e[order(e)])[1:300]

        # lower than background CV
#        mall = apply(data[,-col.numbers],1,mean); sall = apply(data[,-col.numbers],1,sd)
#        cvall = (sall/mall)*100
#        eall = cvall[ match(names(e), names(cvall))]
#        GEN = names(e[(e - eall) < 0])

        # random sampling from background
#        n = 50
#        E = matrix(NA,nrow=length(e),ncol=n)
#        nums = nums[!nums %in% col.numbers]
#        for(kk in 1:n){
#            dat = data[,sample(nums, round(length(nums)*0.7),replace=FALSE)]
#            mall = apply(dat,1,mean); sall = apply(dat,1,sd)
#            cvall = (sall/mall)*100
#            eall = cvall[ match(names(e), names(cvall))]
#            eall[is.na(eall)] = max(eall,na.rm=TRUE)
#            E[,kk] = eall
#        }        
#        pval = rep(NA,length(e))
#        for(kk in 1:length(e)){
#            pval[kk] = wilcox.test(E[kk,], mu=e[kk], alternative="greater")$p.value
#        }
#        names(pval) = names(e)
    }
    print("length GEN:")
    print(length(GEN))

    ## JSD
    pfs = read.delim(paste0("src/tool/code/PF_",species,".txt"),header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
    pfs = toupper(pfs[,1])
    print("Running C++ JSD...")
    data = log2(data+1)
    JSD_value = scJSD(data, classes)
    sJSD_expr_score = scJSD_score(JSD_value)
    JSD_rank_threshold = 10
    temp = sJSD_expr_score[labels(sJSD_expr_score) %in% GEN]
    temp = temp[1:JSD_rank_threshold]
    # just use PFs (either all expressed, or cv <150 too)
    PF = pfs[toupper(pfs) %in% toupper(GEN)]
    print("Printing PF")
    print(PF)
    nonJSD_GEN = PF[!toupper(PF) %in% toupper(labels(temp))]
    print(temp)
    print(nonJSD_GEN)
    GEN = labels(temp)
    all_GEN = c(GEN, nonJSD_GEN)# all_GEN is JSD_GEN + nonJSD_GEN

    ## entropy
    count = 1
    y2 = apply(a4,1,function(x){
        x = as.numeric(x)
        x[x<1e-6] = 0
        e = as.numeric(data[count,])
        e[e<1e-6] = 0      
        bre = seq(from=0, to=(round(max(e)+.05, 1)+0.1), length.out=(nclass.FD(e)+FDplus))
        count <<- count+1
        cut(x, breaks=bre, include.lowest=TRUE)
    })
    H = apply(y2,2,function(x){
        y1d = table(x)
        freqs = y1d/sum(y1d)
	h = -sum(ifelse(freqs > 0, freqs * log2(freqs), 0))	
	c(h, log2(prod(dim(y1d))), h/log2(prod(dim(y1d))))
    })
    H = as.matrix(t(H))
    H[is.nan(H[,3]),3] = 0
    colnames(H) = c('H','H.theoreticalMax','H.maxNormalized')
    H_all_GEN = H[rownames(H) %in% all_GEN,]  # for nonspecific
    H_all_GEN = H_all_GEN[match(all_GEN, rownames(H_all_GEN)),] # for nonspecific
    H = H[rownames(H) %in% GEN,]
    H = H[match(GEN, rownames(H)),]
    
    # Assign integers to symbols
    uniq = na.exclude(unique(c(y2)))
    y3 = apply(y2,2,function(x){match(x,uniq)})
    y3_all_GEN = y3[,match(all_GEN, colnames(y3))] # for nonspecific
    y3 = y3[,match(GEN, colnames(y3))]
    
    ## Compute all 3-TF combinations
    nmax = 3
    nr = choose(length(GEN), 2) + choose(length(GEN), 3)
    new.COM = matrix(nrow=nr, ncol=nmax)
    # 2
    C = t(combn(1:length(GEN), 2))
    new.COM[1:nrow(C),1:2] = C
    nr = nrow(C)+1
    # 3
    C = t(combn(1:length(GEN), 3))
    new.COM[nr:(nr+nrow(C)-1),1:3] = C
    
    ## C++
    print("Running C++ computeMMI...")
    outCombis = compute_MMI_JSD(new.COM-1, y3, H, GEN)
    outCombis$classes = classes
    outCombis$JSD = GEN
    a=do.call('rbind',outCombis$TCs)

    ##  if nonJSD_GEN are all NA
    flag = 0
    print(all(is.na(nonJSD_GEN)))
    print(nonJSD_GEN)
    if(all(is.na(nonJSD_GEN))){
	outCombisF = outCombis
        outCombisF$bestCombis = lapply(outCombisF$bestCombis,function(x){c(x,rep(NA,3))})
	flag = 1
    }
    if(flag==0){
        GEN_index = match(GEN,all_GEN)
    	if(length(nonJSD_GEN) < 3){
    	    nonJSD_GEN = c(nonJSD_GEN, rep(NA, 3-length(nonJSD_GEN)))
        }
        nonJSD_GEN_index = match(nonJSD_GEN,all_GEN)
    
        # compute all combinations of nonJSD_GEN
        tops = outCombis$final_set
        nmax = 3
        nr_new.COM = choose(length(nonJSD_GEN_index), nmax)

        # Split input for saving memory
        topMMI = topTC = 5 # reporting number of each run
        eachRun = 5000
        if(nr_new.COM < eachRun){
            inputIndex = t(as.matrix(c(1,nr_new.COM)))
        } else{
            tmp = seq(eachRun,nr_new.COM,eachRun)
            if(nr_new.COM %% eachRun ==0){
                inputIndex = cbind(seq(1,nr_new.COM,eachRun), c(tmp))
            }else{
                inputIndex = cbind(seq(1,nr_new.COM,eachRun), c(tmp, tmp[length(tmp)] + nr_new.COM %% eachRun))
            }
        }
        nrowF = topMMI * nrow(tops) * nrow(inputIndex)
        combisF_MMI = matrix(NA,nrow=nrowF,ncol=nmax+2)
        combisF_TC = matrix(NA,nrow=nrowF,ncol=nmax+2)
        nrF_MMI = nrF_TC = 0
        for(i in 1:nrow(tops)){
            new.COM = t(combn(c(nonJSD_GEN_index-1), nmax))
            #combinedM = cbind(matrix(tops[i,],nrow=nrow(new.COM),ncol=length(tops[i,]),byrow=TRUE),new.COM)
            for(j in 1:nrow(inputIndex)){
                combinedM = cbind(matrix(tops[i,],nrow=nrow(new.COM[inputIndex[j,1]:inputIndex[j,2],,drop=FALSE]),ncol=length(tops[i,]),byrow=TRUE),new.COM[inputIndex[j,1]:inputIndex[j,2],,drop=FALSE])
                combinedM = t(apply(combinedM,1,sort))
                outCombisF = compute_MMI_finalCombis2(combinedM, y3_all_GEN, H_all_GEN, all_GEN, length(tops[i,]))
                if(nrow(outCombisF$final_set) > 0){
                    # 1. MMI
                    o = order(outCombisF$TCs[,2])
                    m = cbind(i, outCombisF$TCs[o[1:topMMI],2], outCombisF$final_set[o[1:topMMI],])
                    #combisF_MMI[(1:topMMI + nrF_MMI),] = m
                    combisF_MMI[(1:topMMI + nrF_MMI),1:ncol(m)] = m

                    # 2. TC
                    o = rev(order(outCombisF$TCs[,1]))
                    m = cbind(i, outCombisF$TCs[o[1:topTC],1], outCombisF$final_set[o[1:topTC],])
                    combisF_TC[(1:topTC + nrF_TC),1:ncol(m)] = m
                }
                nrF_MMI = nrF_MMI + topMMI
                nrF_TC = nrF_TC + topTC
            }
            rm(new.COM)
        }
        # if no PFs are synergistic, just add them to the specific cores
        if( all(is.na(combisF_MMI)) ){
            outCombisF = outCombis
            outCombisF$bestCombis = lapply(outCombisF$bestCombis, function(x){c(x,nonJSD_GEN[1:3])})
        } else {
            ### Here only report the best one
            # take both top 10 MMI adn TC combis 
            topX = 10
            tTC = combisF_TC[order(-combisF_TC[,2]),,drop=FALSE][1:topX,,drop=FALSE]
            tTC = tTC[!is.na(tTC[,2]),,drop=FALSE]
            tMMI = combisF_MMI[order(combisF_MMI[,2]),,drop=FALSE][1:topX,,drop=FALSE]
            tMMI = tMMI[!is.na(tMMI[,2]),,drop=FALSE]

            # TCs
            TCs = rbind(cbind(tTC[,2],100), cbind(-100, tMMI[,2]))
            TCs = lapply(1:nrow(TCs),function(i){TCs[as.numeric(i),]})
            # final_set
            final_set = rbind(cbind(tops[tTC[,1],,drop=FALSE],tTC[,3:ncol(tTC),drop=FALSE]), cbind(tops[tMMI[,1],,drop=FALSE],tMMI[,3:ncol(tMMI),drop=FALSE]))
            # bestCombis
            bestCombis = lapply(1:nrow(final_set),function(i){ all_GEN[final_set[as.numeric(i),]+1] })
            outCombisF = list() 
            outCombisF$newCOM = NULL
            outCombisF$final_set = final_set
            outCombisF$bestCombis = bestCombis
            outCombisF$TCs = TCs
            outCombisF$JSD = GEN
            outCombisF$nonJSD = nonJSD_GEN
        }
    }
    

    a=do.call('rbind',outCombisF$TCs)
    if(!is.null(a)){
	core = outCombisF$bestCombis[order(a[,2])][1][[1]]
	print(core)
    }else{
        a=do.call('rbind',outCombis$TCs)
	core = outCombis$bestCombis[rev(order(a[,2]))][1][[1]]
        print(core)
    }
    
    # compute fold change of coreTFs
    spec = t(core[1:(length(core)-3)])
    nonspec = t(core[(length(core)-2):length(core)])

    # mean difference
    a2 = a2*1000000 #CPM
    specfc = rowMeans(log2(a2[c(spec),]+1),na.rm=TRUE) - donorE[match(c(toupper(spec)),toupper(names(donorE)))]
    names(specfc) = spec
    nonspec = nonspec[!is.na(nonspec)]
    nonspecfc = rowMeans(log2(a2[c(nonspec),,drop=FALSE]+1),na.rm=TRUE) - donorE[match(c(toupper(nonspec)),toupper(names(donorE)))]
    names(nonspecfc) = nonspec
    names(specfc)[specfc < 0 | is.na(specfc)] = tolower(names(specfc)[specfc < 0 | is.na(specfc)] )
    names(specfc)[is.na(specfc)] = toupper(names(specfc)[is.na(specfc)] )
    names(nonspecfc)[nonspecfc < 0 | is.na(nonspecfc)] = tolower(names(nonspecfc)[nonspecfc < 0 | is.na(nonspecfc)] )
    names(nonspecfc)[is.na(nonspecfc)] = toupper(names(nonspecfc)[is.na(nonspecfc)] )
    names(specfc)[specfc > 0 & !is.na(specfc)] = toupper(names(specfc)[specfc > 0 & !is.na(specfc)] )
    names(nonspecfc)[nonspecfc > 0 & !is.na(nonspecfc)] = toupper(names(nonspecfc)[nonspecfc > 0 & !is.na(nonspecfc)] )
    pecfc = as.matrix(round(specfc[rev(order(specfc))],digits=3))
    nonpecfc = as.matrix(round(nonspecfc[rev(order(nonspecfc))],digits=3))
    
    # print("printing non specific")
    # print(nonpecfc)

    ##output
    colnames(pecfc)=c("FoldChange")
    colnames(nonpecfc)=c("FoldChange")
    
    df1 = pecfc %>% as_tibble(rownames = "Gene") %>% mutate(core = "specific")
    df2 = nonpecfc %>% as_tibble(rownames = "Gene") %>% mutate(core = "pioneer")
    df3 = rbind(df1, df2) %>% mutate(Subpopulation = classes)
    write.table(df3, paste0(outDir,"/core_",classes,".txt"), sep=",", row.names = FALSE, quote = FALSE, col.names=TRUE)
}
