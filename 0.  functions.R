# Function to pull basic info from rawdata eg number of samples/bins/patients
# Requires: chr | start | stop | sample1CN | sample2CN
# Requires sample name in colname, in structure: patientID.sampleNumber
PullDataInfo <- function(rawdata) {
  
  # Dataframe identifying start and stop codon, and chromosome for each bin
  start.stop <- rawdata[,c(1:3)]
  start.stop$bin <- 1:nrow(start.stop)
  start.stop <- start.stop[, c(4, 1, 2, 3)]
  
  # Vector holding sample ID
  sampleIDs <- colnames(rawdata)[-c(1:3)] 
  
  # Vector holding patient identifiers, ie the string before the '.' in sample IDs
  patientIDs <- unique(sub("\\..*", "", colnames(rawdata)))[-c(1:3)] 
  
  # Number of samples
  noSamples <- length(sampleIDs)
  
  # Number of patients
  noPatients <- length(patientIDs)
  
  # list of number of samples per patient
  sampPerPatient <- list()
  for ( i in 1:noPatients ) {
    wd <- rawdata[, which(sub("\\..*", "", colnames(rawdata))==patientIDs[i])]
    sampPerPatient[[i]] <- ncol(wd)
  }
  
  # Number of bins
  noBins <- length(start.stop$bin)
  
  # Visualisation may require a vector to identify of chromosome ends and chromosome midpoints 
  chr.ends <- cumsum(table(start.stop$chr))
  list <- list()
  l <- 1
  for ( i in 1:length(chr.ends) ) {  
    if ( i == 1 ) { 
      list[[l]] <- chr.ends[i]/2 
      l <- l+1
    }
    else { 
      list[[l]] <- chr.ends[i-1] + ((chr.ends[i]-chr.ends[i-1])/2)
      l <- l+1
    }
  }
  chr.mid <- unlist(list)
  chr.ends <- data.frame(start=c(0,chr.ends[-22]), end=(chr.ends), col=c('G','W'))
  
  #average number of bins per patient
  binsPerPatient <- list()
  for (i in 1:length(patientIDs)) {
    wd <- data.frame(rawdata[,which(sub("\\..*", "", colnames(rawdata))==patientIDs[i])])
    binsPerPatient[[i]] <- round(mean(colSums(!is.na(wd))))
  }
  
  
  # Return data
  newData <- list('start.stop'=start.stop, 'sampleIDs'=sampleIDs, 'patientIDs'=patientIDs, 'noSamples'=noSamples, 'noPatients'=noPatients,
                  'sampPerPatient'=sampPerPatient, 'noBins'=noBins,'chr.mid'=chr.mid,  'chr.end'=chr.ends, 'binsPerPatient'=binsPerPatient)
  return(newData)
}

# Function to create 5 new data structures relating to clonality
# Requires: chr | start | stop | sample1CN | sample2CN
# dataInfo must be previously generated from rawdata using PullDataInfo
PullDataClonality <- function(rawdata, dataInfo) {
  
  # A dataframe with a bin per row and a patient per column, with values indicating clonality. 
  # 0=notCNA, 1=subclonalCNA, 2=clonalCNA
  clonal.data <- rawdata[,-c(1:3)]
  l <- 1
  clonal <- list()
  for ( k in 1:nrow(clonal.data) ) {
    for ( i in 1:length(dataInfo$patientIDs) ) {
      wd <- clonal.data[k, which(sub("\\..*", "", colnames(clonal.data))==dataInfo$patientIDs[i])]
      wd <- wd[, !is.na(wd)]
      
      if ( 1 %in% wd | 3 %in% wd ) { #if one of the samples has a mutation, proceed
        if ( length(unique(t(wd)))==1 ) { #if all the same then clonal
          clonal[[l]] <- 2
          l <- l + 1
        }
        else { #different = subclonal
          clonal[[l]] <- 1
          l <- l + 1
        }
      }
      else if ( length(wd)==0 ) { #all MR for that bit are NA
        clonal[[l]] <- NA
        l <- l + 1
      }
      else { #neither sample has a mutation
        clonal[[l]] <- 0 
        l <- l + 1
      }
    }
  }
  
  clonal.data <- data.frame(t(matrix(unlist(clonal), ncol=dataInfo$noBins)))
  colnames(clonal.data) <- dataInfo$patientIDs
  clonal.data[] <- lapply(clonal.data, factor, levels=unique(unlist(clonal.data)))
  
  
  # Dataframe detailing the counts of gains/losses and whether they are subclonal or clonal
  CNA.clo.counts <- data.frame(bin = dataInfo$start.stop$bin, chr1 = NA, chr2 = NA,
                               clonal.aneu = NA, subclonal.aneu = NA, gain = NA, loss = NA,
                               clonal.gain = NA, clonal.loss = NA, clonal.noCNA = NA, 
                               subclonal.gain = NA, subclonal.loss = NA)
  
  CNA.clo.counts$chr1 <- as.numeric(dataInfo$start.stop[match(CNA.clo.counts$bin, dataInfo$start.stop$bin), 2]) 
  CNA.clo.counts$chr2 <- as.numeric(dataInfo$start.stop[match(CNA.clo.counts$bin, dataInfo$start.stop$bin), 2]) 
  CNA.clo.counts$chr2 <- factor(CNA.clo.counts$chr2,levels = rev(seq(1:22))) 
  
  data <- rawdata[,-c(1:3)]
  for ( k in 1:nrow(data) ) { # for a bin
    clonal.all <- clonal.gain <- clonal.loss <- clonal.noCNA <- subclonal.all <- subclonal.gain <- subclonal.loss <- subclonal.noCNA <- 0
    
    for ( i in 1:dataInfo$noPatients ) { # for a patient
      
      # If its NA
      if ( is.na(clonal.data[k,i]) == TRUE ) {
        next
      }
      
      # If its clonal
      if ( clonal.data[k,i] == 2 ) {
        wd <- data[k, which(sub("\\..*", "", colnames(data))==dataInfo$patientIDs[i])]
        wd <- wd[, is.na(wd)!=TRUE]
        
        if ( 3 %in% wd ) { # if both are gains
          clonal.gain <- clonal.gain + 1
        }
        else if ( 1 %in% wd ) { # if both are losses
          clonal.loss <- clonal.loss + 1
        }
      }
      
      # if its subclonal
      else if ( clonal.data[k,i] == 1 ) {
        wd <- data[k, which(sub("\\..*", "", colnames(data))==dataInfo$patientIDs[i])]
        wd <- wd[, is.na(wd)!=TRUE]
        
        if ( 3 %in% wd ) { # if one is a gain
          subclonal.gain <- subclonal.gain + 1
        }
        if ( 1 %in% wd ) { # if one is a loss
          subclonal.loss <- subclonal.loss + 1
        }
      }
      
      # if its no CNA
      else if ( clonal.data[k,i] == 0 ) {
        clonal.noCNA <- clonal.noCNA + 1
      }
    }
    
    CNA.clo.counts$clonal.gain[k] <- clonal.gain
    CNA.clo.counts$clonal.loss[k] <- clonal.loss
    CNA.clo.counts$clonal.noCNA[k] <- clonal.noCNA
    CNA.clo.counts$subclonal.gain[k] <- subclonal.gain
    CNA.clo.counts$subclonal.loss[k] <- subclonal.loss
  }
  
  CNA.clo.counts$clonal.aneu <- CNA.clo.counts$clonal.gain + CNA.clo.counts$clonal.loss
  CNA.clo.counts$subclonal.aneu <- CNA.clo.counts$subclonal.gain + CNA.clo.counts$subclonal.loss
  CNA.clo.counts$gain <- CNA.clo.counts$clonal.gain + CNA.clo.counts$subclonal.gain
  CNA.clo.counts$loss <- CNA.clo.counts$clonal.loss + CNA.clo.counts$subclonal.loss
  
  # dataframe showing: bin | countGain/NoPatient | countLoss/NoPatient | 
  CloFreq <- cbind(bin=CNA.clo.counts$bin, gain=CNA.clo.counts$gain/dataInfo$noPatients, loss=CNA.clo.counts$loss/dataInfo$noPatients)
  
  # dataframe showing what percent of gain and loss are subclonal. On a patient basis: bin | gain | loss
  pcSubclonal <- data.frame(bin=1:dataInfo$noBins, gain=CNA.clo.counts$subclonal.gain / CNA.clo.counts$gain, loss=CNA.clo.counts$subclonal.loss / CNA.clo.counts$loss)
  
  # Count of noCNA, subclonalCNA, and clonalCNA by patient
  patientClo <- as.data.frame(t(sapply(clonal.data, table))) 
  patientClo <- patientClo[, c('0','1','2')]
  colnames(patientClo) <- c('noCNA','subclonal','clonal')
  patientClo$CNA <- patientClo$subclonal + patientClo$clonal
  patientClo$patient <- rownames(patientClo)
  
  # Return data
  newData <- list('clonal.data'=clonal.data, 'CNA.clo.counts'=CNA.clo.counts, 
                  'CloFreq'=CloFreq, 'pcSubclonal'=pcSubclonal, 'patientClo'=patientClo)
  return(newData)
}

# Function to calc pic score per patient across multiple samples
PIC <- function(data, sample.no, index) { 
  # based on table of unique patient nos
  # returns PIC per row (bin), calc across cols as given by index
  # PIC formula: 1- ((CN1/n)^2 + (CN2/n)^2 + (CN3/n)^2), where CN1 is no. of counts of copy number1
  PIC <- 1 - ((rowSums(data[,index]==1, na.rm = TRUE)/sample.no)^2 + 
                (rowSums(data[,index]==2, na.rm = TRUE)/sample.no)^2 + 
                (rowSums(data[,index]==3, na.rm = TRUE)/sample.no)^2)
  return(PIC)
}

# Function to create 4 new data structures relating to diversity
# Requires: chr | start | stop | sample1CN | sample2CN
# There must only be 2 multiregion samples per patient
# dataInfo must be previously generated from rawdata using PullDataInfo
PullDataDiversity <- function(rawdata, dataInfo) {
  
  # Cannot use averaged raw CN across the MR samples in case of variable number of samples per patient.
  # Will therefore can use a PIC.frac for each patient. 
  # This is a continuous, not binary, measure of clonality
  # Define the maximum pic score depending on the number of samples (up to 13)
  max.pics <- list()
  for ( d in 1:50 ) {
    if ( d %% 3 == 0 ) {
      n1 <- n2 <- n3 <- d/3
    }
    else if ( d %% 3 == 1 ) {
      n1 <- ((d-1)/3) + 1
      n2 <- ((d-1)/3)
      n3 <- ((d-1)/3)
    }
    else if ( d %% 3 == 2 ) {
      n1 <- ((d-2)/3) + 1
      n2 <- ((d-2)/3) + 1
      n3 <- ((d-2)/3)
    }
    max.pics[[d]] <- pic.score <- 1 - ( (n1/d)^2 + (n2/d)^2 + (n3/d)^2 )
  }
  
  # A dataframe with a bin per row and a patient per column, with values indicating pic score
  # And a dataframe showing pic.frac per patient
  pic <- list()
  pic.frac <- list()
  for ( i in 1:length(dataInfo$patientIDs)) {
    # Set working data as the cols holding the samples
    wd <- rawdata[, which(sub("\\..*", "", colnames(rawdata))==dataInfo$patientIDs[i])]
    
    # Record the number of samples
    upto <- ncol(wd)
    
    # Use PIC function on wd
    pic[[i]] <- PIC(wd, upto, c(1:upto))
    pic[[i]] <- na.omit(pic[[i]])
    
    # Define the max possible diversity given the number of sample, as maxPIC*number of bins
    # max.ITH <- max.pics[[upto]] * length(pic[[i]])
    max.ITH <- max.pics[[upto]] * dataInfo$binsPerPatient[[i]]
    
    # Store as a dataframe
    pic.frac[[i]] <- data.frame(pic.frac = sum(pic[[i]], na.rm = TRUE)/max.ITH)
  }
  
  pic.data <- as.data.frame(do.call('cbind',pic))
  colnames(pic.data) <- dataInfo$patientIDs
  pic.frac <- do.call('rbind',pic.frac)
  
  # Calulate average PIC per bin as a measure of average bin hetero across patients 
  ave.pic <- dataInfo$start.stop[,c(1:2)]
  ave.pic$avePic <- rowSums(pic.data)/(ncol(pic.data))
  
  # A dataframe of the propotion of genome gained/lost per sample, alongside ith
  pga <- data.frame(t(rawdata[,-c(1:3)]), check.names = FALSE)
  pga <- data.frame(prop.gain=apply(pga,1,function(x) sum(x == 3, na.rm = TRUE)/ncol(pga)),
                    prop.loss=apply(pga,1,function(x) sum(x == 1, na.rm = TRUE)/ncol(pga)))
  pga$prop.aneu <- pga$prop.gain + pga$prop.loss
  pga <- cbind(as.data.frame(lapply(pic.frac, rep, dataInfo$sampPerPatient)),
               pga)
  
  # Return data
  newData <- list('pic.data'=pic.data, 'pic.frac'=pic.frac, 'ave.pic'=ave.pic, 'pga'=pga)
  return(newData)
}

# Function to generate matrices for dip vs aneu, dip vs gain, dip vs loss, gain vs loss, and dip vs gain vs loss
genMatrices <- function(rawdata) {
  
  # Create list to store correlation matrices
  matrices.list <- rep(list(data.frame((rawdata[,-c(1:3)]))),5)
  names(matrices.list) = c('diploid.aneu', 'diploid.gain', 'diploid.loss', 'loss.gain', 'loss.dip.gain')
  
  # Convert to numeric matrices
  matrices.list <- lapply(matrices.list, function(x) {
    y <- data.frame(apply(x, 2, as.numeric))
    rownames(y) <- rownames(x)
    y
  })
  
  # In the diploid/aueploid matrix make 0=diploid, 1=aneuploid.
  matrices.list$diploid.aneu[matrices.list$diploid.aneu == 1 | matrices.list$diploid.aneu == 3] <- 1
  matrices.list$diploid.aneu[matrices.list$diploid.aneu == 2] <- 0
  
  # In the diploid/gain matrix make 0=diploid, 1=gain.
  matrices.list$diploid.gain[matrices.list$diploid.gain == 1] <- NA
  matrices.list$diploid.gain[matrices.list$diploid.gain == 2] <- 0
  matrices.list$diploid.gain[matrices.list$diploid.gain == 3] <- 1
  
  # In the diploid/loss matrix make 0=diploid, 1=loss. 
  matrices.list$diploid.loss[matrices.list$diploid.loss == 1] <- 1
  matrices.list$diploid.loss[matrices.list$diploid.loss == 2] <- 0
  matrices.list$diploid.loss[matrices.list$diploid.loss == 3] <- NA
  
  # In the loss/gain matrix make 0=loss, 1=gain. 
  matrices.list$loss.gain[matrices.list$loss.gain == 1] <- 0
  matrices.list$loss.gain[matrices.list$loss.gain == 2] <- NA
  matrices.list$loss.gain[matrices.list$loss.gain == 3] <- 1
  
  # # In the loss/dip/gain matrix make 1=loss, 2=diploid, 3=gain. Requires no changes
  
  return(matrices.list)
}

# Function to align bin with new bins and use weight average copy number
newAlignBins <- function(bins, cn.list) {
  # bins needs to be a dataframe holding: bin | chr | start | stop, for the bins to align to
  # cn.list is the output from ploidyRecentre with skipcol=3
  
  # Create dataframe for each patient and put into list
  cnBinned.list <- rep(list(bins), length(cn.list))
  
  # Convert to numeric
  cnBinned.list <- lapply(cnBinned.list, function (x) {
    x[] <- apply(x,2,as.numeric)
    x
  })
  
  # Add empty columns to hold output
  for ( i in 1:length(cnBinned.list) ) {
    sampleNo <- ncol(cn.list[[i]])-3 #
    cnBinned.list[[i]] <- do.call("cbind", list(cnBinned.list[[i]], rep(list(NA), sampleNo)))
    colnames(cnBinned.list[[i]])[-c(1:4)] <- colnames(cn.list[[i]])[-c(1:3)]
  }
  
  # Bin incoming dataset (cn.list)
  for ( p in 1:length(cnBinned.list) ) {
    print(paste(p,"/",length(cnBinned.list)))
    # per bin
    for ( b in 1:nrow(cnBinned.list[[p]]) ) {
      chr <- cnBinned.list[[p]]$chr[b]
      start <- cnBinned.list[[p]]$start[b] 
      stop <- cnBinned.list[[p]]$stop[b]
      bin <- cnBinned.list[[p]]$bin[b]
      
      wd <- cn.list[[p]][cn.list[[p]]$chr == chr,]
      wd <- wd[order(wd$start),]
      
      for ( r in 1:nrow(wd) ) {
        
        ## IF START IS BEFORE ROW R 
        if ( start < wd$start[r] ) {
          
          #1.3              # stop is also before start of first row of wd
          if ( stop < wd$start[r] ) {
            cnBinned.list[[p]][which(cnBinned.list[[p]]$bin==bin),-c(1:4)] <- "too early"
            break
          }
          #1.1              # if stop is in row r of wd, or predictor doesnt reach next bin, or this is the last row of wd...
          else if ( dplyr::between(stop, wd$start[r], wd$stop[r]) | stop < wd$start[r+1] | r == nrow(wd) ) { 
            cnBinned.list[[p]][which(cnBinned.list[[p]]$bin==bin),-c(1:4)] <- wd[r,-c(1:3)]
            break
          }
          
          #1.2              # ...or else stop is beyond row r of wd
          else if ( stop > wd$stop[r] ) {
            
            fraction <- list()
            f <- 1 
            
            # find number of extra rows required
            for ( e in 0:(nrow(wd)-r) ) {
              
              # if stop is before row r+e of wd, or stop doesnt reach next bin, or this is the last row of wd...
              if ( stop <= wd$stop[r+e] | stop < wd$start[r+e+1] | r+e == nrow(wd) ) {
                fraction[[f]] <- ( min(wd$stop[r+e], stop)-wd$start[r+e])/(stop-start) # we know it always covers wd$start
                break}
              
              else if (stop > wd$stop[r+e]) {
                fraction[[f]] <- ( min(wd$stop[r+e], stop)-wd$start[r+e])/(stop-start) # take fraction and move onto next
                f <- f + 1
                next}
            }
            
            # take weighted mean
            cnBinned.list[[p]][which(cnBinned.list[[p]]$bin==bin),-c(1:4)] <- colSums( (data.frame(wd[c(r:(r+e)),-c(1:3)]) *unlist(fraction)), na.rm=T) / sum(unlist(fraction))
            break
            
          }
          
        }
        
        ## IF START IS IN ROW R
        else if ( dplyr::between(start, wd$start[r], wd$stop[r]) ) { 
          
          #2.1              # if stop is in row r of wd, or predictor doesnt reach next bin, or this is the last row of wd...
          if ( dplyr::between(stop, wd$start[r], wd$stop[r]) | stop < wd$start[r+1] | r == nrow(wd) ) { 
            cnBinned.list[[p]][which(cnBinned.list[[p]]$bin==bin),-c(1:4)] <- wd[r,-c(1:3)]
            break
          }
          
          #2.2              # ...or else stop is beyond row r of wd
          else if ( stop > wd$stop[r] ) {
            
            fraction <- list()
            f <- 1 
            
            # find number of extra rows required
            for ( e in 0:(nrow(wd)-r) ) {
              
              # if stop is before row r+e of wd, or stop doesnt reach next bin, or this is the last row of wd...
              if ( stop <= wd$stop[r+e] | stop < wd$start[r+e+1] | r+e == nrow(wd) ) {
                fraction[[f]] <- ( min(wd$stop[r+e], stop) - max(wd$start[r+e], start))/(stop-start) # we dont know it always covers wd$start
                break}
              
              else if (stop > wd$stop[r+e]) {
                fraction[[f]] <- ( min(wd$stop[r+e], stop) - max(wd$start[r+e], start))/(stop-start) # take fraction and move onto next
                f <- f + 1
                next}
            }
            
            # take weighted mean
            cnBinned.list[[p]][which(cnBinned.list[[p]]$bin==bin),-c(1:4)] <- colSums( (data.frame(wd[c(r:(r+e)),-c(1:3)]) *unlist(fraction)), na.rm=T) / sum(unlist(fraction))
            break
            
          }
        }
        
        ## IF START IS IN A LATER ROW
        else if ( start > wd$stop[r] ){
          next
        }
        
      }
    }
  }
  return(cnBinned.list)
}

# not in function
'%!in%' <- function(x,y)!('%in%'(x,y))

# annotate genes
gene_anno <- function (top.bins, results) {
  #top.bins data must have bin=col1, chr=col2, start=col3, stop=col4
  #results (gene annotation) must have hgnc_symbol=col1, entrez=col2, chr=col3, start=col4, stop=col5, geneID=col6
  k <- 1
  i <- 1
  l <- 1
  g <- c <- sta <- sto <- id <- b <- NULL
  g.l <- c.l <- sta.l <- sto.l <- id.l <- b.l <- list() #to match genes to bins
  
  for ( k in 1:nrow(top.bins) ) {
    bin <- top.bins[k,1]
    chr <- top.bins[k,2]
    start <- top.bins[k,3]
    stop <- top.bins[k,4]
    
    wd <- results[which(results[,3] == chr),]
    
    for ( i in 1:nrow(wd) ) {
      if ( (chr == wd[i,3]) && ( between(wd[i,4], start, stop) || between(wd[i,5], start, stop)) ) {
        g <- append(g,wd[i,2])
        c <- append(c,wd[i,3])
        sta <- append(sta,wd[i,4])
        sto <- append(sto,wd[i,5])
        id <- append(id,wd[i,6])
        b <- append(b,bin)
      }
    }
    g.l[[l]] <- g
    c.l[[l]] <- c
    sta.l[[l]] <- sta
    sto.l[[l]] <- sto
    id.l[[l]] <- id
    b.l[[l]] <- b
    l <- l+1
    g <- c <- sta <- sto <- id <- b <- NULL
  }
  top.genes <- data.frame(bin=unlist(b.l), entrezgene=unlist(g.l), chr=unlist(c.l), gene.start=unlist(sta.l), gene.stop=unlist(sto.l), geneID=unlist(id.l))
  #top.genes[top.genes==""] <- NA
  top.genes <- na.omit(top.genes) #remove rows with blanks
  #colnames(top.genes) <- c("bin",colnames(wd)[2])
  #top.genes <- merge(top.genes, results, by = colnames(wd)[2]) #add in data from gene annotation data
  #colnames(top.genes)[4] <- "chr"
  return(top.genes)
}

# Function to assess if scalar
is.scalar <- function(x) is.atomic(x) && length(x) == 1L

# Function from Dijk et. al
dijk <- function ( seg_val ,seg_len, ploidy=NULL, purity=NULL ) {
  # check if input seg_val and seg_len are vectors of equal size
  if ( is.vector(seg_val) && is.vector(seg_len) && length(seg_val) != length(seg_len) )
    stop('Segment values (1st input argument) and segment lengths (second input argument) appear not to be column vectors of equal length')
  
  # check if ploidy is scalar argument or empty
  if ( (is.scalar(ploidy) && ploidy > 0) || is.null(ploidy) )
    warning('Ploidy is not a positive scalar or empty');
  
  # check if purity is a scalar between 0 and 1 argument or empty
  if ( (is.scalar(purity) && purity > 0 && purity <=1 ) || is.null(purity) )
    warning('Purity is not a positive scalar or empty');
  
  # specify default range of ploidy for grid search, if input ploidy is empty
  if ( is.null(ploidy) ) {
    ploidy = seq(1.5, 5, 0.01) 
  } 
  
  # specify default range of purity purity for grid search, if input purity is empty
  if ( is.null(purity) ) {
    purity = seq(0.2, 1, 0.01)
  }
  
  # get number of ploidies and purities for grid search
  Nploidy = length(ploidy);
  Npurity = length(purity);
  
  # initialize vectors a1 and a2 from all combinations of ploidy and purity, for the transformation of 
  # measured relative copy number profile (seg_val) to absolute values (q) using   
  #  q = seg_val*a1+a2. 
  a1 = matrix(0, nrow = Nploidy*Npurity, ncol = 1);
  a2 = matrix(0, nrow = Nploidy*Npurity, ncol = 1);
  purity_all = matrix(0, nrow = Nploidy*Npurity, ncol = 1); # vector that contains all purities used in 2D grid search
  ploidy_all = matrix(0, nrow = Nploidy*Npurity, ncol = 1); # vector that contains all ploidies used in 2D grid search
  
  for ( i in 1:Nploidy ) {
    a1[((i-1)*Npurity+1):(i*Npurity)] <- (purity*ploidy[i]+2*(1-purity))/purity;
    a2[((i-1)*Npurity+1):(i*Npurity)] <- -2*(1-purity)/purity;
    
    purity_all[((i-1)*Npurity+1):(i*Npurity)] = purity;
    ploidy_all[((i-1)*Npurity+1):(i*Npurity)] = ploidy[i];
  }
  
  # iniatilize output: CNH_out, ploidy_out and purity_out
  CNH_out = 1;
  purity_out = 0;
  ploidy_out = 0;
  
  # grid search over all ploidies and purities to infer CNH
  for (i in 1:Nploidy*Npurity) {
    # transform relative copy numbers to absolute numbers
    q <- a1[i]*seg_val+a2[i];
    
    # measure distance to closest integer value of each segment
    q_dist_down <- q %% 1; 
    q_dist_up <- 1 - (q %% 1);
    q_dist_min <- pmin(q_dist_up, q_dist_down);
    
    # calculate the mean distance of segments to integer values,
    # weighted by segment length
    CNHnew <- sum(q_dist_min*seg_len) / (sum(seg_len));
    
    # if the new weighted distance to integers is smaller than any
    # previously calculated, replace CNH_out, ploidy_out and purity_out with the new values.
    if (CNHnew < CNH_out) {
      CNH_out <- CNHnew;
      purity_out <- purity_all[i];
      ploidy_out <- ploidy_all[i];
    }
  }
  
  output <- list('CNH_out'=CNH_out, 'purity_out'=purity_out, 'ploidy_out'=ploidy_out)
  return(output)
}

IPH <- function(dataMatrix, dataInfo, nboot) { # need matrix with CN per patient (col) per bin (row)
  max.pics <- list()
  for ( d in 1:100 ) {
    if ( d %% 3 == 0 ) {
      n1 <- n2 <- n3 <- d/3
    }
    else if ( d %% 3 == 1 ) {
      n1 <- ((d-1)/3) + 1
      n2 <- ((d-1)/3)
      n3 <- ((d-1)/3)
    }
    else if ( d %% 3 == 2 ) {
      n1 <- ((d-2)/3) + 1
      n2 <- ((d-2)/3) + 1
      n3 <- ((d-2)/3)
    }
    
    max.pics[[d]] <- pic.score <- 1 - ( (n1/d)^2 + (n2/d)^2 + (n3/d)^2 )
  }
  
  # bootstrap selection of one sample per patient
  # bootstrap
  iph.frac <- list()
  for ( i in 1:nboot ) { 
    
    # select random sample per patient
    sample.selection <- list()
    for ( j in 1:dataInfo$noPatients ) { 
      wd <- dataMatrix[ ,which(sub("\\..*","",colnames(dataMatrix)) == dataInfo$patientIDs[j])]
      sample.selection[[j]] <- sample(wd, 1)
    }
    sample.selection <- do.call("cbind", sample.selection)
    
    # calculate pic per bin and sum across bins
    iph <- 1 - ( (rowSums(sample.selection==1, na.rm = TRUE) / dataInfo$noPatients )^2 + 
                   (rowSums(sample.selection==2, na.rm = TRUE) / dataInfo$noPatients )^2 + 
                   (rowSums(sample.selection==3, na.rm = TRUE) / dataInfo$noPatients )^2 )
    
    iph <- sum(iph, na.rm = TRUE)
    
    # divide by max possible pic
    max.ITH <- max.pics[[dataInfo$noPatients]] * dataInfo$noBins
    iph.frac[[i]] <- iph/max.ITH
  }
  
  return <- mean(unlist(iph.frac))
  
  return(return)
}

# get chromosome lengths and centromere positions
getChrLength <- function(genome = "BSgenome.Hsapiens.UCSC.hg19"){
  g <- BSgenome::getBSgenome(genome, masked=FALSE)
  data.frame(chrom=1:24, length=seqlengths(g)[1:24])
}
.chrAsNum <- function(tbl){
  tbl$chrom <- gsub("chr", "", tbl$chrom)
  tbl$chrom[tbl$chrom=="X"] <- 23
  tbl$chrom[tbl$chrom=="Y"] <- 24
  tbl$chrom <- as.numeric(tbl$chrom)
  tbl[order(tbl$chrom),]
}
getCentromeres <- function( genome="hg19" ){
  mySession <- try(browserSession("UCSC"), silent=TRUE)
  # In case of failure, try another mirror
  if(inherits(mySession, "try-error"))
    mySession <- rtracklayer::browserSession("UCSC",
                                url="http://genome-euro.ucsc.edu/cgi-bin/")
  genome(mySession) <- genome
  obj <- rtracklayer::ucscTableQuery(mySession, table="gap")
  tbl <- getTable(obj)
  tbl <- tbl[tbl$type=="centromere", c("chrom", "chromStart", "chromEnd")]
  colnames(tbl)[2:3] <- c("centromerStart", "centromerEnd")
  .chrAsNum(tbl)
}
makeHg19 <- function(){
  tbl <- merge(getChrLength(), getCentromeres(), by="chrom")
  cumlen <- c(0, cumsum(as.numeric(tbl$length))[-nrow(tbl)])
  cbind.data.frame(tbl, cumlen=cumlen)
}

# custome themes
theme_custom <- function(){
  theme_classic() %+replace%  
    theme(
      plot.margin=margin(t=0.2,r=0.5,b=0.5,l=0.5,"cm"),
      panel.background = element_blank(),
      plot.title = element_text(size=28, colour='black',face='bold', hjust=0),
      axis.line = element_line(size = 0.5, colour = "black"),
      axis.title = element_text(size=28, colour='black'),
      axis.text = element_text(size=28, colour='black'),
      axis.ticks.length=unit(0.2, "cm")
    )
}
theme_legend <- function(){
  theme_classic() %+replace%  
    theme(
      plot.margin = unit(c(t=-100,r=0,b=-100,l=0), "cm"),
      legend.position = "top",legend.direction="horizontal",
      legend.title = element_text(size=28, colour='black',face='bold'),
      legend.margin = margin(grid::unit(c(t=-100,r=0,b=-100,l=0),"cm")),
      legend.text = element_text(size=28, colour='black'),
      legend.key.height = grid::unit(0.8,"cm"),
      legend.key.width = grid::unit(1.4,"cm")
    )
}

# Function to plot frequency and clonality of gains and losses
cloFreqPlot <- function(clonalityData, dataInfo, annotation, title=NULL, ylab="Fraction of patients with aberration", xlab="Chromosome", colourChoice=c("#1B7837","#5AAE61","#A6DBA0","#D9F0D3","#E7D4E8","#C2A5CF","#9970AB","#762A83")) {
  
  # pull CloFreq
  CloFreq <- data.frame(clonalityData$CloFreq)
  
  # Make losses negative to create a mirror plot
  CloFreq$loss <- CloFreq$loss*-1 
  
  # Gather into long form
  CloFreq <- gather(CloFreq,CNA,freq,-bin) 
  
  # Look up the percent of the gain/loss that is subclonal
  CloFreq$pcSubclonal <- NA
  
  for ( i in 1: nrow(CloFreq) ) {
    bin <- CloFreq$bin[i]
    CNA <- CloFreq$CNA[i]
    
    if ( CNA == 'gain' ) {
      CloFreq$pcSubclonal[i] <- clonalityData$pcSubclonal$gain[bin] *100
    }
    if ( CNA == 'loss' ) {
      CloFreq$pcSubclonal[i] <- clonalityData$pcSubclonal$loss[bin] *100
    }
  }
  
  # Position label along y axis
  for ( i in 1:nrow(annotation) ) {
    bin <- annotation$x[i]
    data <- CloFreq[which(CloFreq$bin==bin),]
    annotation$y[i] <- data$freq[which.max(abs(data$freq))]
  }
  
  # Create plots
  fig <- ggplot() + 
    geom_rect(data = dataInfo$chr.end[which(dataInfo$chr.end$col=='W'),], 
              aes(NULL,NULL,xmin=start, xmax=end),
              fill = alpha("#CCCCCC", 0.1),
              ymin = -1,
              ymax = 1) +
    
    geom_bar(data=CloFreq, aes(fill=pcSubclonal, y=freq, x=bin),
             position="stack", stat="identity", width = 1) +
    scale_fill_gradientn(colours = colourChoice) +
    
    ggtitle(title) +
    scale_x_continuous(expand = c(0,0), name=xlab, breaks=dataInfo$chr.mid, labels = c(1:18,'\n19','20','\n21','22')) +
    scale_y_continuous(expand = c(0,0), name=ylab, limits = c(-1,1), breaks = c(seq(-1,1,0.2)), 
                       labels = c(1,0.8,0.6,0.4,0.2,0,0.2,0.4,0.6,0.8,1)) +
    
    geom_hline(yintercept = 0, size=0.3, color="black") +
    
    geom_point(data = annotation, aes(x = x, y = y),
               shape = 18, size = 0) +
    geom_text_repel(data = subset(annotation, y >= 0),
                    aes(x = x, y = y, label = label),
                    size = 7,
                    nudge_y = 0.2,direction = 'x',
                    arrow = arrow(length = unit(0.015, "npc"))) +
    geom_text_repel(data = subset(annotation, y < 0),
                    aes(x = x, y = y, label = label),
                    size = 7,
                    nudge_y = -0.2,direction = 'x',
                    arrow = arrow(length = unit(0.015, "npc"))) +
    annotate('text', y=c(0.9,-0.9), x=c(50,50), 
             label=c('italic(Gains)','italic(Losses)'), size=10, parse=TRUE, hjust = 0) +
    theme_custom() +
    theme(legend.position = "none", plot.title = element_text(hjust=0))
  
  # Create legend
  galo.legend <- as_ggplot(cowplot::get_legend(fig + 
                                                 guides(fill = guide_colorbar(title="% subclonal", label.position = "bottom",
                                                                              title.position = "left", title.vjust = 0.9)) +
                                                 theme_legend() ))
  
  # remove negative for loss to allow for calculation of correlation
  CloFreq$freq[which(CloFreq$freq<0)] <- CloFreq$freq[which(CloFreq$freq<0)]*-1
  return(list(CloFreqDF = CloFreq, plot = fig, legend=galo.legend))
}

# Function to plot AVERAGE PGA against actual CNA diversity
predictITHfromPGA <- function(diversityData, title=NULL, colourChoice="#273046", toplim, ylab="CNA diversity", xlab="PGA") {
  pga <- diversityData
  pga$sample <- rownames(pga)
  pga$patient <- sub('\\..*', '',rownames(pga))
  pga <- data.frame(cbind(aggregate(pga$prop.aneu, list(pga$patient), mean),
        aggregate(pga$pic.frac, list(pga$patient), mean)[2]))
  colnames(pga) <- c("patient","pga","pic.frac")
  #pga <- gather(pga, CNA, proportion, -sample, -patient, -pic.frac)
  
  fig <- 
    #ggplot(data=pga[which(pga$CNA=="prop.aneu"),], aes(y=pic.frac, x=proportion)) + 
    ggplot(data=pga, aes(y=pic.frac, x=pga)) +
    geom_smooth(method = "lm", formula = y ~ x, color="black", fill=colourChoice) +
    #geom_point(size=3, shape=21, color="black", fill=colourChoice) +
    geom_point(size=5, shape=21, color="black", fill=alpha(colourChoice, 0.9)) +
    ggpmisc::stat_poly_eq(formula = y ~ x, 
                          aes(label = paste(after_stat(adj.rr.label), "*\", \"*", after_stat(p.value.label), "*\"\"", sep = "")),
                          parse = TRUE, label.x = 0, label.y = 'top', size=10, ) +
    ggtitle(title) +
    coord_fixed(ratio=1) +
    scale_y_continuous(ylab, expand = c(0.01,0), limits = c(0,toplim)) +
    scale_x_continuous(xlab, expand = c(0.01,0), limits = c(0,toplim), breaks = c(seq(0,toplim,0.2)), labels = c(seq(0,toplim,0.2))) +
    theme_custom()
  
  #r2 <- lm(pic.frac ~ proportion, data=pga[which(pga$CNA=="prop.aneu"),])
  r2 <- lm(pic.frac ~ pga, data=pga)
  
  return(list(plot=fig, r2=r2))
}

# Function to plot heatmat using geom_tile using pariwise correlation matrix
cnaCoocHM <- function(corrMatrix, dataInfo) {
  # Melt into long form for geom_tile
  melted <- melt(corrMatrix, na.rm = FALSE)
  colnames(melted) <- c('first_bin','second_bin','correlation')
  
  # add 2x chr data for faceting
  melted <- cbind(melted, chr1 = as.numeric(dataInfo$start.stop[match(melted$first_bin, dataInfo$start.stop$bin), 2])) 
  melted <- cbind(melted, chr2 = as.factor(dataInfo$start.stop[match(melted$second_bin, dataInfo$start.stop$bin), 2]))
  melted$chr2 <- factor(melted$chr2, levels = rev(c(seq(1,22,1))) ) 

  
  # Prepare heatmaps 
  hm <- ggplot(melted, aes(first_bin, second_bin, fill = correlation)) +
    geom_tile() +
    scale_fill_distiller(palette ="BrBG", na.value="grey85", limits = c(-1,1), breaks=c(-1,0,1), guide=guide_colorbar()) +
    facet_grid(chr2 ~ chr1, scales='free', space = 'free') +
    
    scale_x_continuous(name=NULL, breaks = c(round(dataInfo$chr.mid, digits = 0)),  
                       labels = c(1:18,'\n19','20','\n21','22')) + 
    scale_y_continuous(name=NULL, breaks = c(round(dataInfo$chr.mid, digits = 0)),
                       labels = c(1:18,'19    ','20','21    ','22')) +
    
    theme_custom() +
    theme(
      legend.position = "none",
      panel.spacing = unit(0.2, "lines"),
      panel.background = element_blank(),
      panel.border=element_blank(),
      strip.text = element_blank()) 
  
  # create legend
  hm.legend <- as_ggplot(cowplot::get_legend(hm + 
                                               guides(fill = guide_colorbar(title="Correlation", label.position = "bottom",
                                                                            title.position = "left", title.vjust = 0.9)) +
                                               theme_legend()))   
  
  return(list(plot=hm, legend=hm.legend))
  
}
  
# Function to plot a split violin plot
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1, "group"]
  newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
      1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

# Function to predict using input beta model, and plot against actual
predictITHfromBeta <- function(rawData, dataInfo, predictors, actualITH = NULL, beta, colourChoice = "#273046", title=NULL, ylab = NULL, xlab = NULL) {
  multiReg.in <- t(rawData[,-c(1:3)]) 
  multiReg.in <- data.frame(apply(multiReg.in, 2, as.character), check.names = FALSE)
  multiReg.in <- data.frame(lapply(multiReg.in, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)
  rownames(multiReg.in) <- dataInfo$sampleIDs
  multiReg.in <- multiReg.in[,c(predictors$bin)]
  names(multiReg.in) = paste("bin_", names(multiReg.in), sep="")
  multiReg.in <- dummy_cols(multiReg.in, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
  predictedITH <- data.frame(predicted = predict(beta, multiReg.in))
  
  predictedITH <- data.frame(cbind(actualITH, predictedITH))
  
  if (!is.null(actualITH)) {
    plot <- ggplot(predictedITH, aes(x = actual, y = predicted)) + 
      geom_smooth(method = "lm", formula = y ~ x, color="black", fill=colourChoice) +
      geom_line(aes(group = patient), size=0.2, colour = colourChoice) +
      geom_point(shape=21, fill = alpha(colourChoice,0.9), size = 5) +
      geom_line(aes(x = actual, y = actual), linetype = "dashed") +
      ggtitle(title) +
      ylab(ylab) +
      xlab(xlab) +
      ggpmisc::stat_poly_eq(formula = y ~ x, 
                            aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                            parse = TRUE, label.x = 'left', label.y = 'top', size=10) +
      theme_custom() +
      theme(legend.position = "none", plot.title = element_text(hjust=0))
  } else {
    plot <- NULL
  }
  
  return(list(predictedITH=predictedITH, plot=plot))
}
  
  
