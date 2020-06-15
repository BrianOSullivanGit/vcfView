# vcfView
# Graphical tool to view, plot, process and re-filter a VCF.


# Make sure the required libraries are installed.

if (!requireNamespace("shiny", quietly = TRUE))
  BiocManager::install("shiny")
if (!requireNamespace("ggplot2", quietly = TRUE))
  BiocManager::install("ggplot2")
if (!requireNamespace("VariantAnnotation", quietly = TRUE))
  BiocManager::install("VariantAnnotation")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("annotate", quietly = TRUE))
  BiocManager::install("annotate")
if (!requireNamespace("reshape", quietly = TRUE))
  BiocManager::install("reshape")
if (!requireNamespace("pracma", quietly = TRUE))
  BiocManager::install("pracma")
if (!requireNamespace("neutralitytestr", quietly = TRUE))
  BiocManager::install("neutralitytestr")


# Load these libraries.
library(shiny)
library("ggplot2")
library(VariantAnnotation)
library('org.Hs.eg.db')
library(annotate)
library(reshape)  # For melt() function.
library("pracma") # For lsqnonneg()

library("ggrepel")
library(neutralitytestr)


# Track change in Max/Min AF freq range.
# We do this to avoid making the AF density plot reactive to these values.
previousMinF <<- 0
previousMaxF <<- 0

# Gene annotation dataframe
geneAnnotateDf <<- NULL
uniqueGeneList <<- NULL

# Data relating to current plot.
plot_data <- list(pathname="", selectedSample="", snvIdx = NULL, vcfFile = NULL, gg_b_object=NULL, result=NULL, dp_breaks=NULL)

# Current ggplot objects so we can ggsave them sgv file if required.
# TODO!!!! consolidate some of these globals into a common data structure &
# remove globals which have become redundant as a result to code development.
targetVcfGgplot <<- NULL
TopOneGgplot <<- NULL

# Cache the max depth / break of any density distribution we have plotted.
maxDepthInBreaks <<- 0



# Mutational Order
# The order in which they will be displayed in a plot.
SIGNATURE_ORDER = c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", 
                    "C[C>A]T", "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", 
                    "T[C>A]G", "T[C>A]T", "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "C[C>G]A", 
                    "C[C>G]C", "C[C>G]G", "C[C>G]T", "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", 
                    "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", 
                    "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "G[C>T]A", "G[C>T]C", 
                    "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T", "A[T>A]A", 
                    "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", 
                    "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", 
                    "T[T>A]T", "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", 
                    "C[T>C]G", "C[T>C]T", "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", 
                    "T[T>C]C", "T[T>C]G", "T[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", 
                    "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", 
                    "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")



# Function sets up filter option checkbox element and help text as tooltip.
filterCheckboxHtml <- function(inputId, value, helpText) {
  paste0("<br><div class='vcfViewTooltip'><label class='checkbox'><input type='checkbox' style='margin-left:-14px; margin-top:-1px;' name='", inputId, "' value='", value,"'>", value,"<span class='vcfViewTooltiptext vcfViewTooltip-left'>",helpText,"</span></label></div>")
}

# Make sure the choices of gene in the selectiveInput only contains the genes within the
# new subset the user has created (gChoices passed from the caller).

updateGeneSelectionList <- function(session, gChoices, currentGeneSelected)
{   
  
  if(any(is.na(gChoices)))
  {
    updateSelectizeInput(session,"geneList", choices = c("None Found"))
    return()
  }
  
  # Make sure the input is a gene referenced in our VCF before we go any further with
  # updating the selectizeInput list (to avoid triggering an update every time a key is pressed).
  if((any(currentGeneSelected == uniqueGeneList)  || currentGeneSelected == "None Found"))
  {
    
    # Make sure there is something to update it with.
    if(length(gChoices) > 0)
    {    
      # Make sure the current choice remains the current choice if still possible.
      # If not, let it take the default from the new list.
      if(any(gChoices == currentGeneSelected))
      {
        updateSelectizeInput(session,"geneList", choices = gChoices, selected = currentGeneSelected)
      }
      else
      {
        updateSelectizeInput(session,"geneList", choices = gChoices, selected = gChoices[1])
      }
      
    }
    
  }
}

# Function returns the Max bin depth that has been plotted on a density plot.
maxMedianBinDepth <- function(nBins, result)
{
  
  g=ggplot(result, aes(x = AF, y = ..density..)) + geom_histogram(bins=nBins)
  
  # Create gg build object to examine the details of this plot (bins/breaks etc.).
  gg_bd=ggplot_build(g)$data[[1]]
  
  dp_breaks = c()
  
  # Interpolate from allele frequency to depth..
  # Now create the median depth vector with entries to match each bin.
  # Make sure there are more than 2 bins!
  # (just in case we were passed something unreasonable).
  if (nBins >= 2){
    i = 1
    # Run through the breaks finding out what the median depth is within each one.
    while (i < length(gg_bd$xmin)) {
      dp_breaks[i] = median(result$DP[result$AF > gg_bd$xmin[i] & result$AF <= gg_bd$xmin[i+1]])
      # Check to see if this bin was empty.
      # Just set mean to zero in that case.
      if(is.na(dp_breaks[i])) {
        dp_breaks[i] = 0;
      }
      i = i+1
    }
    dp_breaks[i] = median(result$DP[result$AF > gg_bd$xmin[i]])
    # Check to see if this bin was empty.
    # Just set mean to zero in that case.
    # NA's will mess things up down the line.
    if(is.na(dp_breaks[i])) {
      dp_breaks[i] = 0;
    }
  }
  
  # Now create a colours vector to reflect the average depth at each bin.
  # Also we can't leave zeros in the vector passed to colorRampPalette() or it will skip them.
  # Just offset by the min value.
  # Make a second copy of the dp_breaks vector so we can display the actual depth on mouse over etc. (TODO!!!!)
  
  # There are many ways to setup the colour vector to get the visual resolution you may require.
  # Generally there will be a couple of bins sampled to very high depth but the average will be much lower than these.
  # Smooth it out by taking the log to highlight the detail at the lower depth bins.
  dp_colour_breaks = log(dp_breaks)
  dp_colour_breaks = (dp_colour_breaks * 100000) + 1
  
  maxDp = max(dp_colour_breaks)
  
  if(maxDp > maxDepthInBreaks) {
    maxDepthInBreaks <<- maxDp
  }
  else {
    maxDp = maxDepthInBreaks
  }
  
  #####################################################################    
  # Record where R puts the breaks.
  # Cache to ggplot build object global so it can be accessed on mouse click etc.
  # Also cache the result so it can be used by the "TopOne" inset plot.
  plot_data$gg_b_object <<- gg_bd
  plot_data$result <<- result
  plot_data$dp_breaks <<- dp_breaks
  
  
  
  # Create an interpolation pallet to reflect this.
  colours <- colorRampPalette(c("red", "green"))(maxDp)[dp_colour_breaks]
  
  return(colours)
  
}

# processVcf()
# Load the VCF.
# Convert to the required datastructure.
# Load annotation.

processVcf <- function(csvPath, session) {
  
  if(is.null(csvPath))
  {
    return(NULL)
  }
  
  # Parse the filter names and descriptions from the VCF metadata.
  
  filetype = summary(file(csvPath))$class
  
  con = NULL
  
  if(filetype == "gzfile") {
    con = gzfile(csvPath,"r")
  }
  else {
    con = file(csvPath, "r")
  }
  
  ID = c()
  Description = c()
  genome_str = c()
  tumourSlotNum = 0
  normalSlotNum = 0
  normalSampleName = ""
  tumourSampleName = ""
  
  
  while ( TRUE ) {
    line = readLines(con, n = 1)
    # Keep reading, as long as we are reading metadata lines.
    if ( length(line) == 0 || !grepl("^#",line)) {
      break
    }
    
    # Parce out the filter descriptions for the tooltip.
    if(grepl("^##FILTER=<ID=",line)) {
      line = gsub("^##FILTER=<ID=", "", line)
      line = gsub("Description=\"", "", line)
      line = gsub("\".*", "", line)
      fentry = unlist(strsplit(line,','))
      ID = c(ID, fentry[1])
      Description = c(Description, fentry[2])
    }
    
    # Try and deduce which version of the genome we are working with
    # from the GATK commandline so we can load the correct annotation etc.
    if(grepl("^##GATKCommandLine.*ID=[mM]u[tT]ect2",line)) {
      line = gsub("^##GATKCommandLine.*--reference[ ][ ]*.*GRCh", "", line)
      
      # Alternative reference sequence may be identified as (for MuTect2 GATK3),
      line = gsub("^##GATKCommandLine.*reference_sequence=.*GRCh", "", line)
      
      genome_str =  gsub("[^0-9].*$", "", line)
      
      # Sub 37 with 19 to conform to UCSC naming convention in their
      # annotation libraries which we load later on.
      genome_str =  gsub("^37$", "19", genome_str)
      
    }
    
    # Find out which slot is allocated to tumour/normal.
    # This will vary depending on the GATK3/GATK4.   
    if(grepl("^#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tNORMAL\\tTUMOR",line)) {
      # GATK3 found.
      tumourSlotNum = 2
      normalSlotNum = 1
    }
    
    
    # If the samples are named in the metadata, parse them out.
    if(grepl("^##normal_sample=",line)) {
      normalSampleName =  gsub("^^##normal_sample=", "", line)
    }
    if(grepl("^##tumor_sample=",line)) {
      tumourSampleName =  gsub("^^##tumor_sample=", "", line)
    }
    
    
    if(grepl(paste0("^#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t", normalSampleName, "\\t", tumourSampleName),line)) {
      # Looks like GATK4, normal slot 1 and tumour slot 2. 
      normalSlotNum = 1
      tumourSlotNum = 2
    }
    
    if(grepl(paste0("^#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t", tumourSampleName, "\\t", normalSampleName),line)) {
      # Looks like GATK4, normal slot 2 and tumour slot 1.
      
      normalSlotNum = 2
      tumourSlotNum = 1
      
    }
    
    
    
  } 
  close(con)
  
  # If GATK3 was not found in metadata and no tumour/normal sample names were
  # found in the metadata then assume we're dealing with GATK4 with tumour slot 1
  # normal slot 2. Not ideal I know but best we can do.
  # TODO!!!! Inform user with a warning here.
  if(tumourSlotNum == 0)
  {
    tumourSlotNum = 1
    normalSlotNum = 2
  }
  
  filterDescriptions = data.frame(ID,Description)
  
  colnames(filterDescriptions) <- c("ID","Description")
  
  
  # Load UCSC BSgenome and corresponding annotation.
  # TODO!!!! If these libraries are not installed ask user
  # if they want to install them now (move library() call from the
  # top to here).
  #  genome_str = gsub("GRCh","Hsapiens.UCSC.hg",genome_str)
  genome = paste0("BSgenome.Hsapiens.UCSC.hg", genome_str)
  txdb_str = paste0("TxDb.Hsapiens.UCSC.hg", genome_str,".knownGene")
  
  # If user has loaded hg38, then closes and loads a hg37 there may be
  # issues if hg38 is still loaded
  # Advise user too run (.packages()) and see what hg38 stuff is loaded,
  # Then unload them, ie.,
  #
  # detach("package:TxDb.Hsapiens.UCSC.hg38.knownGene", unload=TRUE)
  # detach("package:BSgenome.Hsapiens.UCSC.hg38", unload=TRUE)
  # 
  # TODO!!!! Are you sure this was the issue?
  # Perhaps it was just a memory issue / running out of RAM..
  if(genome == "BSgenome.Hsapiens.UCSC.hg19")
  {
    packagesAlreadyLoaded = (.packages())
    grepMatch = grepl("^BSgenome.Hsapiens.UCSC.hg3[8-9]$|^TxDb.Hsapiens.UCSC.hg3[8-9].knownGene$", packagesAlreadyLoaded)
    if(any(grepMatch))
    {
      problemLibs=paste0("     ","detach(\"package:",packagesAlreadyLoaded[grepMatch],"\", unload = TRUE)\n",collapse = "")     
      warningString=paste0("\nYou may need to unload the following libraries with\n\n",problemLibs,"\n before running this application with hg19 genome based UCSC annotation.\n")  
      warning(warningString)      
    }   
  }
  
  
  # Next 8 lines are very convoluted, however thay are they best I can come up with at the moment.
  cmd=paste0('if (!requireNamespace("',genome,'", quietly = TRUE))\n{\nprint("Loading BSgenome UCSC annotation.")\nBiocManager::install("',genome,'")\n}') # Install required library if not there already.
  source(textConnection(cmd))
  
  cmd=paste0('if (!requireNamespace("',txdb_str,'", quietly = TRUE))\n{\nprint("Loading BSgenome UCSC annotation.")\nBiocManager::install("',txdb_str,'")\n}') # Install required library if not there already.
  source(textConnection(cmd))
  
  cmd=paste0('library("',genome,'")')  # Equivalent to "library(genome)"
  source(textConnection(cmd))
  
  cmd=paste0('library("',txdb_str,'")')  # Equivalent to "library(txdb_str)"
  source(textConnection(cmd))
  
  # We now parse the required reference from the vcf metadata (above).
  
  ref_genome <- base::get(genome)
  ref_organism <- GenomeInfoDb::organism(ref_genome)
  ref_style <- seqlevelsStyle(ref_genome)
  genome_name <- genome(ref_genome)[[1]]
  if (!(class(ref_genome) == "BSgenome")) 
    stop("Name of a BSgenome object not found.")
  
  vcf = readVcf(csvPath, genome_name)
  
  # TODO!!!! sampleNames here is now redundant, remove it and references..
  sampleNames = rownames(colData(vcf))
  
  
  #    allVariant = ScanVcfParam(samples=sampleNames[1])
  allVariant = NULL
  
  
  # Just stick to the autosomal & XY info. from GenomeInfoDb,
  # ignore the rest..
  groups <- c(extractSeqlevelsByGroup(species = ref_organism, style = ref_style, group = "auto"), extractSeqlevelsByGroup(species = ref_organism, style = ref_style, group = "sex"))
  
  # Seqlevels in the VCF file often do not match those in the BSgenome/TxDb.
  # Hack to fix the age old problem of whether or not to prefix a chromosome number
  # with 'chr'. In some VCF's the chromosome contigs will be referenced by "chr1"
  # while in others it will be just "1". Prefix with a "chr" if its not there
  # already and move on.
  vcf <- renameSeqlevels(vcf, gsub("(^[XY]$)","chr\\1", gsub("(^[0-9]+$)","chr\\1", seqlevels(vcf))))
  #vcf <- renameSeqlevels(vcf, gsub("(^[0-9XY]+$)","chr\\1", seqlevels(vcf)))
  #vcf <- renameSeqlevels(vcf, paste0("chr", seqlevels(vcf)))  
  
  # Prune all the viral etc. contigs out of your VCF.
  # Just keep autosomal & XY. Then we will be able to annotate with TxDb etc..
  vcf <- keepSeqlevels(vcf, groups, pruning.mode = "tidy")
  
  
  ################ 
  
  
  # TODO!!!! update on success only...
  plot_data$vcfAnnotateObject <<- vcf 
  
  
  # Where there are multi allelic variants at a given loci,
  # we need to filter out the one with that gives the max TLOD value.
  # This code section replaces the logic that was provided previously in simple_parse.c
  
  start_time <- Sys.time()
  
  # Cache these slots to speed up looping through them.
  # Again the '[,1]' below subsets out the first sample (tumour sample) from geno slot.
  tlod_l = info(vcf)$TLOD
  alt_l = rowRanges(vcf)$ALT
  # (Take the allelic fraction from the first sample,, ie., tumour)
  af_l = geno(vcf[,tumourSlotNum])$AF
  af_norm_l = geno(vcf[,normalSlotNum])$AF
  
  
  
  # Format field differes between earlier versions (GATK3 & GATK4) of Mutect2.
  # This allows for that variation.
  
  if(exists("vcf[,tumourSlotNum])$DP"))
  {
    dp_l = as.integer(geno(vcf[,tumourSlotNum])$DP)
  }
  else
  {
    dp_l = lapply(geno(vcf[,tumourSlotNum])$AD, sum)
  }
  
  if(exists("vcf[,normalSlotNum])$DP"))
  {
    dp_norm_l = as.integer(geno(vcf[,normalSlotNum])$DP)
  }
  else
  {
    dp_norm_l = lapply(geno(vcf[,normalSlotNum])$AD, sum)
  }
  
  
  # Number of records in original VCF.
  oidx_length = length(info(vcf)$TLOD)
  
  # TODO!!!! put in a sanity check that all the lengths in next 3 lines are the same.
  tlod_t = integer(length(info(vcf)$TLOD))
  af_t = integer(length(geno(vcf[,tumourSlotNum])$AF))
  
  af_norm_t = integer(length(geno(vcf[,normalSlotNum])$AF))
  
  alts_t = character(length(rowRanges(vcf)$ALT))
  
  maLogicV = (lengths(info(vcf)$TLOD) > 1)
  
  # Are there any multiallelic entries in this VCF?
  
  # Does this VCF contain Multiallelic entries?
  # If so, for each multiallelic record, pick the one with the highest TLOD to use.  
  if(any(maLogicV))
  {
    # Indices to Multiallelic entries.
    maIdx = which(maLogicV)
    # Indices to single allelic entries.
    saIdx = which(!maLogicV)
    
    # Update the multiallelic entries to single allelic entries containing the one which gives the max TLOD.
    # Big, bad, time consuming loop.
    # TODO!!!! is there a better way of doing this? (btw., lapply is as bad/worse)
    for(i in maIdx) {
      # Cache this value to prevent having to perform the unlist twice.
      # (to speed things up)
      tmpTlod = unlist(tlod_l[i])
      
      maxIdx = which(tmpTlod == max(tmpTlod))
      # Check to see if there is a tie for max TLOD value.
      # If there is just pick the first one.
      # TODO!!!! should we do something else in this case?
      if (length(maxIdx) > 1)
      {
        maxIdx = maxIdx[1]
      }
      
      # Update the vectors with the element corresponding to the max TLOD value.
      tlod_t[i] = tlod_l[[i]][maxIdx]
      af_t[i] = af_l[[i]][maxIdx]
      af_norm_t[i] = af_norm_l[[i]][maxIdx]
      alts_t[i] = as.character(alt_l[[i]][maxIdx])
    }
    
    # Now populate single allelic entries in these vectors.
    # We know that all these 'saIdx' indices are single allelic entries so it's
    # safe to unlist them (pulling this stuff out of the loop saves time).
    tlod_t[saIdx] = unlist(info(vcf)$TLOD[saIdx])
    af_t[saIdx] = unlist(geno(vcf)$AF[,tumourSlotNum][saIdx])
    af_norm_t[saIdx] = unlist(geno(vcf)$AF[,normalSlotNum][saIdx])
    alts_t[saIdx] = as.character(unlist(rowRanges(vcf)$ALT[saIdx]))
  }
  else
  {
    tlod_t = as.numeric(info(vcf)$TLOD)
    af_t = as.numeric(geno(vcf)$AF[,tumourSlotNum])
    af_norm_t = as.numeric(geno(vcf)$AF[,normalSlotNum])
    alts_t = as.character(unlist(rowRanges(vcf)$ALT))
  }
  
  
  # print(head(data.frame(alts_t,tlod_t,af_t,af_norm_t)))
  
  # Now pull out the reference alleles.
  # We know these are one element lists already so no need to do anything other than the next line.
  refs_t = as.character(rowRanges(vcf)$REF)
  
  # Store alt_t alleles as we will need to modify SNVs among them for purposes of getting a mutational contexts
  ma_t = alts_t
  
  # Now get a mutational contexts vector for all variants.
  # Look at the trinucleotide sequence around each variant
  # (ie., one nucleotide before and one after it, start(vcf) - 1, end(vcf) + 1)).
  # TODO!!!! the genome name and version needs to match the VCF!!!!
  # Look in the VCF comments and alert the user.
  
  tnctx_t = unname(as.character(getSeq(get(genome), seqnames(vcf), start(vcf) - 1, end(vcf) + 1)))
  
  # Remember alt allele and context of all SNVs needs to be expressed in terms of C>X or T>X transitions
  # (for mutational context plot).
  # If it is not already in this format then we need to express the transition in terms of the reverse strand.
  # Find the mutations which need context adjustment.
  #
  # The subset below finds all SNVs with ref allele = "A" or "G"
  # (ie, stuff that needs to be rewritten  so it's expressed as a transition on the opposite strand)
  x = !is.na(match(alts_t, c("A", "C", "G", "T"))) & !is.na(match(refs_t, c("A", "G")))
  
  # Make a subset these trinucleotide mutational contexts
  y = tnctx_t[x]  
  
  # Change the context of these SNVs to reverse complement
  # of that context (ie. express in terms of opposite strand).
  y = reverse(chartr('ATGC', 'TACG', y))
  
  # replace the subset with reverse complement
  tnctx_t[x] = y
  
  # Finally change the corresponding ALT allele in this subset to its complement also.
  y = chartr('ATGC', 'TACG', alts_t[x])
  alts_t[x] = y
  
  end_time <- Sys.time()
  
  
  # Put in chr number with length(seqnames(vcf)) etc...
  # Again the '[,1]' below subsets out the first sample (tumour sample).
  
  # (BTW. leave POS = unname(as.character(ranges(vcf[,1]))) as an IRanges object for now. )
  # ( Later we can modify with POS = unname(as.character(ranges(vcf[,1]))) if required  )
  vcfTable = data.frame(CHROM = as.character(seqnames(vcf[,tumourSlotNum])), POS = ranges(vcf[,tumourSlotNum]), FILTER = as.character(rowRanges(vcf)$FILTER), AF = af_t, AF_IN_NORM = af_norm_t, DP = as.numeric(dp_l), DP_IN_NORM = as.numeric(dp_norm_l), T_LOD = as.numeric(tlod_t), REF = refs_t, ALT = ma_t, MA = alts_t, TNC = tnctx_t, OIDX = seq(1:oidx_length), CLASS = tnctx_t)
  
  
  ##################################################################################
  
  
  # How the hell did this get changed to a type integer?????????
  # TODO!!!! find out......
  # vcfTable$FILTER = as.character(vcfTable$FILTER)
  # print(head(rowRanges(vcf)$FILTER))
  
  # TODO!!!!
  # Why do we need to do this again here????
  # Why didn't it work in the line above vcfTable = data.frame.....
  vcfTable$FILTER = as.character(vcfTable$FILTER)
  
  # Find out what types of filters were flagged in this VCF.
  fTypes=unlist(strsplit(vcfTable$FILTER,';'))
  fTypes=sort(fTypes[!duplicated(fTypes)])
  
  # Remove the 'PASS' indications from the filter descriptions.
  # We are going to display that indication separately
  # (as it is not really a 'filter' in that sence.)
  fTypes=fTypes[!fTypes=="PASS"]
  
  
  # Now load the cancer signatures from COSMIC.
  # Store them for reference later (ie., in cosmicComboBreakdown() etc.)
  # Note: this has now been updated to v3 signatures as per below...
  # TODO!!!! should we cache this in case user doesn't have internet access??
  # cancer_signatures = read.table("https://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt", sep = "\t", header = TRUE)
  
  # if (dim(cancer_signatures)[1] != 96) 
  # stop("Signatures file does not have 96 tri-nucleotide enteries. Are you sure this is an SBS cosmic file?")
  
  
  ####### v3 changes ##############
  
  # v3signatures = read.csv("https://dcc.icgc.org/api/v1/download?fn=/PCAWG/mutational_signatures/Signatures/SP_Signatures/SigProfiler_reference_signatures/Sigprofiler_Exome_Signatures/sigProfiler_exome_SBS_signatures.csv")
  
  # TODO!!!! remove...
  # Load locally while internet down..
  v3signatures = read.csv("/home/sully/BIO_INFORMATICS/VCF_PARSE/Cancer_Informatics/COSMIC_v3/sigProfiler_exome_SBS_signatures.csv")
  
  
  if (nrow(v3signatures) != 96) 
    stop("Signatures file does not have 96 tri-nucleotide enteries. Are you sure this is an SBS cosmic file?")
  
  mTypes = paste0(substring(v3signatures[,2], 1, 1),"[",v3signatures[,1],"]",substring(v3signatures[,2], 3, 3))
  
  v3signatures = cbind(mTypes = mTypes, v3signatures)
  
  new_order = match(SIGNATURE_ORDER, v3signatures$mTypes)
  
  # Reorder cancer signatures dataframe
  v3signatures = v3signatures[as.vector(new_order),]
  
  # Add trinucletiode changes names as row.names
  row.names(v3signatures) = mTypes
  
  # Keep only 96 contributions of the signatures in matrix
  v3signatures = as.matrix(v3signatures[,4:ncol(v3signatures)])
  
  #################################
  
  ######### v2 legacy stuff ###########
  # Match the order of the mutation types to order they will be displayed on plot.
  # new_order = match(SIGNATURE_ORDER, cancer_signatures$Somatic.Mutation.Type)
  
  # Reorder cancer signatures dataframe
  # cancer_signatures = cancer_signatures[as.vector(new_order),]
  
  # Add trinucletiode changes names as row.names
  # row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
  
  # Keep only 96 contributions of the signatures in matrix
  # cancer_signatures = as.matrix(cancer_signatures[,4:33])
  
  # These signatures are stored in the vcfFile list (below) for future reference
  # (ie., later in cosmicComboBreakdown()).
  
  # Put all the relevant information together in the global vcfFile list.   
  # vcfFile <- list(pathname=csvPath, selectedSample=sampleNames[1], data=vcfTable, filterTypes=fTypes, filterDescriptions=filterDescriptions, cosmicSignatures=cancer_signatures, activeFilters=NULL, dp_breaks=NULL, xmin=NULL)
  
  
  ######################################################################################
  
  # Put all the relevant information together in the global vcfFile list.   
  vcfFile <- list(pathname=csvPath, selectedSample=sampleNames[1], tumourSampleName=tumourSampleName, normalSampleName=normalSampleName, data=vcfTable, filterTypes=fTypes, filterDescriptions=filterDescriptions, cosmicSignatures=v3signatures, activeFilters=NULL, dp_breaks=NULL, xmin=NULL)
  
  # Record it in the plot_data object. 
  plot_data$vcfFile <<- vcfFile 
  
  # Annotation stuff....
  
  # We parse this now from the vcf metadata and so can replace this line..
  # txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  
  # Next 2 lines are very convoluted, however thay are they best I can come up with at the moment.
  cmd=paste0("txdb <- ",txdb_str)
  source(textConnection(cmd))
  
  # Pull out annotation of coding sequences (takes about 10sec)
  coding <- predictCoding(vcf, txdb, seqSource=Hsapiens)
  
  # Pull out the gene names from geneid's in Entrez Gene identifiers annotation database.
  geneId = coding$GENEID
  geneNames = geneId
  symsIdx = !is.na(geneId)
  geneNames[symsIdx] = unlist(lookUp(geneId[symsIdx], 'org.Hs.eg', "SYMBOL"))
  
  
  # TODO!!!! FINISH HERE
  # Sometimes we get back NAs from Entrez Gene identifiers annotation database.
  # In that case set the name to 'UNKNOWN'
  if(sum(is.na(geneNames)) > 0)
  {
    
  }
  
  geneNames[is.na(geneNames)] = "UNKNOWN"
  
  # Update global list of all genes we have in this VCF.
  uniqueGeneList <<- as.character(unique(geneNames))
  
  # Pull out the transcript names using their ids from the transcript annotation package.
  k <- keys(txdb, keytype = "TXNAME")
  tx2txid <- select(txdb, k, "TXID", "TXNAME")
  
  # Pull out the location of these variants within the protein.
  # A mutation can affect a number of bases and therefore a number of amino acids.
  # The PROTEINLOC column of the coding slot is an IRanges object (list).
  # When only one amino acid is impacted it contains its location in the chain.
  # If a number of adjacent amino acids are impacted the list contains 2 elements,
  # a start and end location.
  # We only deal with the start location for plotting purposes.
  # The ref and variant amino acid fields (REFAA, VARAA) indicate how many
  # peptides the modification affects.
  #
  # Lists get messy in R.
  # Break this up into a 2 stage process, handling the single peptide changes first
  # and then those with multiple changes.
  
  p_loc_start <- vector(mode="integer", length=length(coding))
  tmpIdx = which(lengths(coding$PROTEINLOC) == 1)
  p_loc_start[tmpIdx] = unlist(coding$PROTEINLOC[tmpIdx])
  
  tmpIdx = which(lengths(coding$PROTEINLOC) != 1)
  x=(coding$PROTEINLOC[tmpIdx])
  p_loc_start[tmpIdx] = unname(unlist(lapply(x, `[[`, 1)))
  
  # Finally pull out the coding sequence lengths from the dm3 annotation.
  # We use this to work out the original protein length.
  # This will take about 10 secs to run.
  # TODO!!!! You know you could probably move this to the very begining when the shiny app is opened,
  # even before you load the VCF. That would split up the 25/30sec load time into two stages at least.
  dm3_txlens <- transcriptLengths(txdb, with.cds_len=TRUE,
                                  with.utr5_len=FALSE,
                                  with.utr3_len=FALSE)
  
  match_idx = match(coding$TXID, dm3_txlens$tx_id)
  
  # Protein length is deduced from the coding sequence length.
  # The '3' below is a codon length, the minus one corrects for the stop codon at the end.
  # So PROTEIN LENGTH = (CDSLENGTH/3)-1
  # as per next line, 
  
  # Pull all the required annotation together in a data.frame.
  # CONSEQUENCES below are frameshift, nonsense, nonsynonymous or synonymous.
  # We will refer back to this annotation dataframe for all annotation queries by the user.
  #
  # TODO!!!! Take out match(coding$QUERYID, vcfTable$OIDX) to a local variable from next line to speed it up.
  geneAnnotateDf <<- data.frame(CHROM = vcfTable$CHROM[match(coding$QUERYID, vcfTable$OIDX)],
                                VARIANT = ranges(vcf[match(coding$QUERYID, vcfTable$OIDX),tumourSlotNum]),
                                AF = vcfTable$AF[match(coding$QUERYID, vcfTable$OIDX)],
                                DP = vcfTable$DP[match(coding$QUERYID, vcfTable$OIDX)],
                                GENE_SYMBOL=geneNames,
                                QUERYID = coding$QUERYID,
                                CONSEQUENCE = coding$CONSEQUENCE,
                                TXID = coding$TXID,
                                TXNAME = tx2txid[coding$TXID,"TXNAME"],
                                PROTEIN_MOD_START = p_loc_start,
                                PROTEIN_ORIG_LENGTH = (((dm3_txlens$cds_len[match_idx])/3)-1),
                                REFAA = coding$REFAA,
                                VARAA = coding$VARAA)
  
  
  # To fix bug with some MuTect2 GATK3 files.  
  # In GATK3 some variants are rejected without listing a TLOD value.
  # To be pragmatic, we will just remove these variants from our analysis so NA's are
  # not introduced down the line.
  vcfTable = vcfTable[!is.na(vcfTable$T_LOD),]
  vcfFile$data <- vcfTable
  
  
  ##############################################################################################################
  
  return(vcfFile)
  
}

# 
# Returns a tri-nucleotide context mutational matrix (Sanger style).
# Matrix is compatable with MutationalPatterns library.
getSnvMutationalMatrix <- function(result, lower_frange = NULL, upper_frange = NULL)
{
  
  # Subset out all the SNVs within required allele freq. range.
  snvIdx = which(result$AF > lower_frange &
                   result$AF <= upper_frange &
                   !is.na(match(result$MA, c("A", "C", "G", "T"))) & 
                   !is.na(match(result$REF, c("A", "C", "G", "T"))))
  
  # TODO!!!!
  # We need to record this snvIdx in a global.
  # We will referr to this if the users requests us to save this SNV subset as a new VCF.
  plot_data$snvIdx <<- snvIdx
  
  # Create a context matrix to pass to ggplot and MutationalPatterns calls etc.
  
  # Give it one of every mutational context to start off.
  # This will prevent NA's being introduced down the line if a freq range that
  # does not contain any mutations of a given context is selected.
  # We will correct for every context having an additional variant below before returning.
  
  ctx = c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", 
          "C[C>A]T", "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", 
          "T[C>A]G", "T[C>A]T", "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "C[C>G]A", 
          "C[C>G]C", "C[C>G]G", "C[C>G]T", "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", 
          "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", 
          "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "G[C>T]A", "G[C>T]C", 
          "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T", "A[T>A]A", 
          "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", 
          "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", 
          "T[T>A]T", "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", 
          "C[T>C]G", "C[T>C]T", "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", 
          "T[T>C]C", "T[T>C]G", "T[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", 
          "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", 
          "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")
  
  
  for (i in snvIdx) {  
    x=strsplit(as.character(result$TNC[i]), "")
    # Rem., result$ALT is the ALT allele as it was listed in the VCF,
    # result$MA is the allele specified in terms of the cosmic mutational context etc.
    type = paste0(x[[1]][1],"[",x[[1]][2],">",as.character(unlist(result$MA[i])),"]",x[[1]][3])  
    ctx = c(ctx, type)
  }
  
  # TODO!!!! URHERE!!!!
  # We need to store the snvIdx somewhere so we can save the variants in the VCF within this range in a separate VCF.
  # We also need to record them, together with their annotation and tri-allelic context in a csv file.
  
  # Get a summary breakdown.
  tmp=summary(as.factor(ctx))
  
  # Match the order of the mutation types to order they will be displayed on plot.
  new_order = match(SIGNATURE_ORDER, names(tmp))
  
  # Reorder signatures vector.
  tmp = tmp[new_order]
  names(tmp) <- SIGNATURE_ORDER
  
  # Correct for the extra mutation of each context we added above.
  tmp = tmp-1
  
  # Create the mutation matrix.
  new_mat = matrix(tmp, 96, byrow=TRUE)
  rownames(new_mat) <- names(tmp)
  colnames(new_mat) <- "Tumour Sample"
  
  
  # Return the mutation matrix.
  # It will be used to plot the 96 context barplot
  # or to fit this sample to a combination of the v3 cosmic known signatures.
  return(new_mat)
  
}


# Neutral evolution test plot..
# 
nevolution <- function(result, lower_frange = 0.45, upper_frange = 0.55)
{
  
  
  # Subset out all the SNVs within required allele freq. range.
  snvIdx = which(result$AF > lower_frange &
                   result$AF <= upper_frange)
  
  
  
  n <- neutralitytest(result$AF[snvIdx], fmin = lower_frange, fmax = upper_frange)
  
  return(lsq_plot(n))
}

# Lollypop peptide plot..
# We also need to pass the session id to this as it updates the gene list (under 'select gene').
lpoppep <- function(session, result, lower_frange = 0.45, upper_frange = 0.55, geneName)
{
  
  # Avoid this event triggering while user is typing each character of gene name in.
  if(!any(geneName == uniqueGeneList) && geneName != "None Found")
  {
    return()
  }
  
  # Get a consistent legend colour map going here.
  #
  # Colours maftools uses..
  # Nonsence = 230:89:90 = #E6595A
  # Missence = 107:183:103 = #6BB767
  # Frame_shift_del = 92:155:197 = #5C9BC5
  # Splice site = 251:160:71
  
  # Leveles we have from predictCoding() call with transcript DB TxDb.Hsapiens.UCSC.hg38.knownGene
  #
  #> levels(geneAnnotateDf$CONSEQUENCE)
  #[1] "frameshift"    "nonsense"      "nonsynonymous" "synonymous"
  # 
  # So our colours are.
  cseqColours = c("#5C9BC5","#E6595A","#E6595A","#F0DADA")
  # Name them to keep the levels consistent.
  names(cseqColours) <- levels(geneAnnotateDf$CONSEQUENCE)
  colScale <- scale_colour_manual(name = "CONSEQUENCE",values = cseqColours)
  
  
  
  # Subset out all the SNVs within required allele freq. range.
  snvIdx = which(result$AF > lower_frange &
                   result$AF <= upper_frange)
  
  
  # pull out the annotation for this this subset using the oidx. 
  annSubsetIdx = !is.na(match(geneAnnotateDf$QUERYID, result$OIDX[snvIdx]))
  
  # This is the subset we're interested in.
  annaSubset = geneAnnotateDf[annSubsetIdx,]
  
  # Check if the user has just selected a new range having failed to find a protein impacted gene in last range.
  if(geneName=="None Found")
  {
    # If this is the case just show the first gene in the list.
    geneName=annaSubset$GENE_SYMBOL[1]
  }
  
  
  
  gChoices = as.character(unique(annaSubset$GENE_SYMBOL))
  
  # if gChoices != NA or has length zero..
  
  if(!is.na(gChoices) && (length(gChoices) > 0))
  {
    
    # Leave the current gene of interest selected if possible,
    # Pick first one if not.
    if(!any(geneName == gChoices))
    {
      # If the selected gene isn't in the list of impacted genes in this range
      # then we need to reset that list in the selectizeInput.
      updateGeneSelectionList(session, gChoices, gChoices[1])
    }
    else
    {
      updateGeneSelectionList(session, gChoices, geneName)
    }
    
  }
  else
  {
    updateGeneSelectionList(session, NA, NA)
  }
  
  
  hitList = annaSubset[grepl(paste("^",geneName,"$", sep = ""), annaSubset$GENE_SYMBOL),]
  
  
  # Remove duplicates.  ????
  #  hitList <- unique(hitList)
  
  
  # Move on if there's nothing to see here..
  if(!(nrow(hitList) > 0))
  {
    # OK, rather than returning NULL, printout an effective 'NULL plot'
    # so the user knows nothing was found.
    hit = data.frame(PROTEIN_MOD_START=c(0), AF=c(0))
    p=ggplot(data=hit, aes(x=PROTEIN_MOD_START, y=AF)) +
      #  scale_y_continuous(expand = expand_scale(mult = c(1, 2, 3, 4), add = c(0, 4,4,0))) +
      xlim(0, 100) +
      labs(y= "Allele frequency", x = "AA location within peptide") +
      geom_segment( aes(x=0, xend=100, y=0, yend=0)) +
      geom_segment( aes(x=0, xend=0, y=0, yend=-0.03*0.3)) +
      geom_segment( aes(x=100, xend=100, y=0, yend=-0.03*0.3)) 
    
    if(!is.na(geneName) && geneName != "None Found")
    {
      geneString = paste0(" for ",geneName)
    }
    else
    {
      geneString = "."
    }
    
    p = p + annotate("text",x=0,y=.3,hjust=0,label=paste0(hit[1,]$TXNAME,"\n","No protein impacting variants found",geneString,".\n(within this allele frequency range)"), size=3, fontface = "bold") +
      theme(
        legend.title = element_blank(),
        plot.title = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 9, face = "bold"),
        axis.title.y = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)
      )
    # Clear plot data, there is nothing to save to SVG here.    
    TopOneGgplot <<- NULL
    
    return(p)
    #    return(NULL)
  }  
  
  
  # Get the canonical (longest) sequence here and run with that for the plot..
  # Will usually be first in list but take no chances..
  # TODO!!!! all user to specify other sequencies to plot based on tx name.
  # TODO!!!! inform user if on the off chance, all variants are not listed in the plot displayed.
  hit = hitList[which(hitList$PROTEIN_ORIG_LENGTH==max(hitList$PROTEIN_ORIG_LENGTH))[1],]
  
  # Gather all variants in this transcript.
  hit = hitList[hitList$TXID == hit$TXID,]
  
  # Order this data.frame so the variants listed in the title
  # correspond to the order they are listed in the title.
  # Skip this, we don't need to do this if we don't re-order the title now.
  # TODO!!!! check this works as expected.
  # hit <- hit[order(hit$PROTEIN_MOD_START),]
  
  # For variants whose impact straddle more than one protein (ie. some double base substitutions),
  # all annotation entried will share the same 'PROTEIN_MOD_START'
  # We need to update this value for such entries by incrementing each subsequent protein
  # modification by one as we move past the modification start so we can display it
  # correctly on the plot.
  # Start in reverse order.
  # Remember loci and AA position in protein move in reverse order.
  # (locus coding first amino acid in peptide > locus of last).
  previousVarRecrd = NA
  for(i in nrow(hit):1) {
    # Check if a second protein is impacted by the same variant record.
    if( !is.na(previousVarRecrd) && (previousVarRecrd$VARIANT.names == hit[i,]$VARIANT.names) )
    {
      hit[i,]$PROTEIN_MOD_START = previousVarRecrd$PROTEIN_MOD_START + 1
    }
    previousVarRecrd = hit[i,]
  }
  
  # Lollypop plot of hit.
  # Split up this for the sake of readability.
  
  # TODO!!!! SPLICE_SITE variants....
  # This plot will not include splice site variants (like a maftools plot)!!
  # You have the information available here.
  # For example,
  #
  
  # allv <- locateVariants(vcf, txdb, AllVariants())
  
  # > levels(allv$LOCATION)
  # [1] "spliceSite" "intron"     "fiveUTR"    "threeUTR"   "coding"     "intergenic"
  # [7] "promoter
  
  # This will also contain a set of indices into the vcfTable, ie.,
  
  #> head(allv$QUERYID)
  #[1] 1 1 1 2 2 2
  
  # It could be pulled in quite easily here.
  
  # However I am not sure how maftools get away with displaying a splice_site mutation on a peptide plot ??
  # (the splice site is in the intron).
  # I guess they label the adjacent amino acid in the exon (?).
  
  # Also, (TODO!!!!) filter out the synonymous variants, doesn't make much sence to display them on a peptide plot(?maybe?)
  
  
  tLabel = paste0(hit$REFAA,hit$PROTEIN_MOD_START,hit$VARAA, sep="")
  tLabel[hit$CONSEQUENCE == "frameshift"] = paste0(hit[hit$CONSEQUENCE == "frameshift","PROTEIN_MOD_START"], "_f.shift:",hit[hit$CONSEQUENCE == "frameshift","VARIANT.width"])
  
  # Get pointers to insertions and deletions so we can better annotate them.
  # (as per http://atlasgeneticsoncology.org/Educ/NomMutID30067ES.html)
  deletionIdx = nchar(hit$REFAA) > nchar(hit$VARAA)
  insertionIdx = nchar(hit$REFAA) < nchar(hit$VARAA)
  
  
  if(is.logical(deletionIdx) && sum(deletionIdx) > 0)
  {
    tLabel[deletionIdx] = paste0(sub(paste0("^",hit[deletionIdx,]$VARAA),"",hit[deletionIdx,]$REFAA),hit[deletionIdx,"PROTEIN_MOD_START"]+1,"del")
  }
  
  if(is.logical(insertionIdx) && sum(insertionIdx) > 0)
  {
    # TODO!!!!
    # Improve annotation of deletions.
  }
  
  
  # Figure out where to place the legend.
  # where is the biggest gap, on the left or right?
  
  XlimRatio = 1.2
  
  if(XlimRatio*(hit$PROTEIN_ORIG_LENGTH[1] - max(hit$PROTEIN_MOD_START)) > min(hit$PROTEIN_MOD_START))
  {
    # Put it to the right.
    expandRight = 0.05
    legX = 0.8
    legY = 0.8
    XlimRatio = 1.2
    atateVjust = .7
    
    # Put it to the bottom or leave it on the top?
    if(hit[which(hit$PROTEIN_MOD_START==max(hit$PROTEIN_MOD_START))[1],"AF"] > max(hit$AF)/2)
    {
      legY = 0.4
    }
    
  }
  else
  {
    # Put it to the left.
    # Remember left top is also where the title will be.
    expandRight = 0
    # legX = 0.2
    legX = 0.08
    # legY = 0.3
    legY = 0.71
    XlimRatio = 1
    atateVjust = 1.9
  }
  
  # TODO!!!! you need to add a label to the title of this plot indicating the frequency
  # range from which it was drawn, just like the other TopOne plots.
  # paste(strwrap("this is a, fairly loooooong string but there are longer..", width = 30),collapse="\n")  
  
  # Plot the axes etc. first.
  p=ggplot(data=hit, aes(x=PROTEIN_MOD_START, y=AF)) +
    #  scale_y_continuous(expand = expand_scale(mult = c(1, 2, 3, 4), add = c(0, 4,4,0))) +
    scale_y_continuous(expand = c(expandRight, 0)) +
    xlim(0, XlimRatio*hit[1,]$PROTEIN_ORIG_LENGTH) +
    labs(y= "Allele frequency", x = "AA location within peptide") +
    geom_segment( aes(x=0, xend=hit$PROTEIN_ORIG_LENGTH, y=0, yend=0)) +
    geom_segment( aes(x=0, xend=0, y=0, yend=-0.03*max(hit$AF))) +
    geom_segment( aes(x=hit$PROTEIN_ORIG_LENGTH, xend=hit$PROTEIN_ORIG_LENGTH, y=0, yend=-0.03*max(hit$AF))) 
  
  
  
  # Now the long pin segments.
  p = p + ggtitle(paste(hit[1,]$GENE_SYMBOL," (",paste(strwrap(paste0(hit[order(hit$PROTEIN_MOD_START),]$VARIANT.names, collapse = ", "), width = 110),collapse="\n"),")", sep="")) + theme_bw() + theme(plot.title = element_text(hjust=0)) +
    geom_segment( aes(x=hit$PROTEIN_MOD_START, xend=hit$PROTEIN_MOD_START, y=0, yend=hit$AF), size = 0.3, color="#808080") 
  
  # Fill in the pinhead points, legend and fonts.
  p = p + geom_point(aes(x=hit$PROTEIN_MOD_START, y=hit$AF, color=hit$CONSEQUENCE), size=3) + colScale +
    geom_text_repel(data = hit, aes(x = hit$PROTEIN_MOD_START, y = hit$AF, label=tLabel), size=3.5, colour="blue", box.padding = .4, segment.size = .25) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    annotate("text",x=0,y=(max(hit$AF)+max(hit$AF)/5),hjust=0,vjust=atateVjust,label=hit[1,]$TXNAME, size=3, fontface = "bold") +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 7),
      legend.position = c(legX,legY),
      legend.background = element_rect(fill=alpha(NA, 0.6)),
      plot.title = element_text(size = 9, face = "bold"),
      axis.title.x = element_text(size = 9, face = "bold"),
      axis.title.y = element_text(size = 9, face = "bold"),
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9)
    )
  # Record this plot in TopOneGgplot incase the user wants to save it as svg.
  TopOneGgplot <<- p
  p
}

# TODO!!!! update comment..
# Possible the most exciting function name, ever..
cosmicComboBreakdown <- function(result, lower_frange = 0.45, upper_frange = 0.55)
{   
  # Get the mutational matrix for the SNVs within this allele frequence range.
  new_mat = getSnvMutationalMatrix(result, lower_frange, upper_frange)
  
  
  # Sum new_max to see how many SNVs we got.
  total_snvs = sum(new_mat)
  
  # Normalize.
  # If we want a realistic resid.norm then we need to do this.
  # First make sure the selected range returned some single base substitutions.
  # If it did not, we will print a null graph to let user know.
  numSbs = sum(new_mat)
  if(numSbs > 0)
    new_mat = new_mat/sum(new_mat)
  
  
  # How do the cosmic signatures breakdown within the observed mutational contexts
  # within this range? Lets see..
  
  # Fit mutation matrix to the COSMIC mutational signatures:
  # (using nonnegative linear least-squares fit method from pracma library)
  nonNegLsqFit = lsqnonneg(plot_data$vcfFile$cosmicSignatures, new_mat[,1])
  
  # Fix the order they will appear on the barplot.
  # (in ascending etc.. you need to do this or ggplot will re-order alphabetically)
  # df <- data.frame(signatures=as.numeric(gsub("Signature.", "", colnames(plot_data$vcfFile$cosmicSignatures))),
  #                 contribution=nonNegLsqFit$x)
  
  df <- data.frame(signatures=gsub("SBS", "", colnames(plot_data$vcfFile$cosmicSignatures)),
                   contribution=nonNegLsqFit$x)
  
  
  df$signatures <- factor(df$signatures, levels = df$signatures)
  
  # Plot it.
  p = ggplot(data=df, aes(x=signatures, y=contribution)) +
    geom_bar(stat="identity", color="blue", fill="white") +
    ggtitle(paste0("Signature Contribution\n[", signif(lower_frange,6), "," , signif(upper_frange,6), "), ", total_snvs, " SNVs. Norm of residuals = ", round(nonNegLsqFit$resid.norm,digits=6), ".")) +
    labs(y= "contribution (SNVs)", x = "v3 cosmic signature") + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_text(aes(label=signatures), vjust=-.5, color="blue", size=rel(2.5))
  
  # Record this plot in TopOneGgplot incase the user wants to save it as svg.
  TopOneGgplot <<- p
  return(p)
}

frange <- function(result, lower_frange = 0.45, upper_frange = 0.55)
{     
  # Get the mutational matrix for the SNVs within this allele frequence range.
  new_mat = getSnvMutationalMatrix(result, lower_frange, upper_frange)
  
  # Sum new_max to see how many SNVs we got..(make sure to do this before normalising!!)
  total_snvs = sum(new_mat)
  
  # Normalize.
  if(total_snvs > 0)
    new_mat = new_mat/sum(new_mat)
  
  # Format the data and scale for the plot.
  max_new_mat = max(new_mat)
  
  # Work out y-axis scaling (roughly)
  if(max_new_mat < 0.1)
    ymax = 0.1
  else if(max_new_mat < 0.2)
    ymax = 0.2
  else
    ymax = 0.3
  
  colors = c(
    "#2EBAED", "#000000", "#DE1C14",
    "#D4D2D2", "#ADCC54", "#F0D0CE")
  
  context = c(rep(c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT"), 3),
    rep(c(
      "ATA", "ATC", "ATG", "ATT",
      "CTA", "CTC", "CTG", "CTT",
      "GTA", "GTC", "GTG", "GTT",
      "TTA", "TTC", "TTG", "TTT"), 3))
  
  substitution = rep(c('C>A','C>G','C>T','T>A','T>C','T>G'), each = 16)
  substring(context, 2, 2) = "."
  df = data.frame(substitution = substitution, context = context)
  rownames(new_mat) = NULL
  df2 = cbind(df, as.data.frame(new_mat))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  value = NULL
  
  # Plot it.
  p = ggplot(data = df3, aes(x = context, y = value, 
                             fill = substitution, width = 1)) + geom_bar(stat = "identity", 
                                                                         colour = "black", size = 0.2) +
    ggtitle(paste0("[", signif(lower_frange,6), "," , signif(upper_frange,6), "), ", total_snvs, " SNVs.")) +
    scale_fill_manual(values = colors) + 
    facet_grid(. ~ substitution) + ylab("Relative contribution") +
    coord_cartesian(ylim = c(0, ymax)) + scale_y_continuous(breaks = seq(0, 
                                                                         ymax, 0.1)) + guides(fill = FALSE) + theme_bw() + 
    theme(axis.title.y = element_text(size = 12, vjust = 1), 
          axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
          axis.text.x = element_text(size = 5, angle = 90, 
                                     vjust = 0.4), strip.text.x = element_text(size = 9), 
          strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank(), 
          panel.spacing.x = unit(0, "lines"))
  
  # Record this plot in TopOneGgplot incase the user wants to save it as svg.
  TopOneGgplot <<- p
  return(p)
  
  
}




ui <- shinyUI(
  fluidPage(
    
    tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "app.css")),
    
    tags$head(tags$script(src = "app.js")),
    
    sidebarLayout(position = "right",
                  
                  sidebarPanel(
                    HTML(paste0('<div class="InvisableAtStart">')),
                    HTML(paste0('<div style="padding:0px; align: center; text-align: center;font-size:24px"><img style="width: 70%; margin-bottom:3px; margin-top:3px;" src="dhelix_logo_small.png"></div>')),
                    # The 'PASS' button and filters list.
                    HTML(paste0('<div align="center">')),
                    uiOutput("listFilters"),
                    
                    HTML(paste0('<div class="control-label" style="font-size: 12px; align: center; padding-top: 2%;">Threshold filters:</div>')),
                    
                    htmlOutput("tumourNameTag"),
                    
                    HTML(paste0('<div id="tlodThresh_slider" class="InvisableAtStart" style="display:inline-block; padding-right: 50px; margin-top: 4px;">')),
                    
                    numericInput("tlodThresh", "TLOD:", value = 0, min = 0, max = 1000, width = "80px"),
                    
                    HTML(paste0('</div>')),
                    
                    HTML(paste0('<div id="depthThresh_slider" class="InvisableAtStart" style="display:inline-block;">')),
                    numericInput("depthThresh", "Depth:", value = 0, min = 0, max = 1000, width = "80px"),
                    HTML(paste0('</div>')),
                    
                    htmlOutput("normalNameTag"),
                    
                    HTML(paste0('<div id="tlodThresh_slider" class="InvisableAtStart" style="display:inline-block; padding-right: 50px; margin-top: 4px;">')),
                    numericInput("afInNormal", "AF ceiling:", value = 0, min = 0, max = 1, step = 0.005, width = "65px"),
                    HTML(paste0('</div>')),
                    HTML(paste0('<div id="depthThresh_slider" class="InvisableAtStart" style="display:inline-block;">')),
                    numericInput("depthThreshInNormal", "Depth:", value = 0, min = 0, max = 1000, width = "80px"),
                    HTML(paste0('</div>')),
                    
                    
                    HTML(paste0('</div>')),
                    width = 3
                  ),
                  
                  
                  
                  mainPanel(   
                    # Background Logo.
                    HTML(paste0('<div style="width: 100% ; height: 100%; background-image: url(\'dhelix_logo.png\'); background-repeat: no-repeat; background-position: 50% 30%; background-size: 50%;">')),
                    
                    # The "TopOne" plot is a plot with a transparent background.
                    # It is positioned on top of the VAF plot and displays content when the user
                    # selects an allele frequency 'window' range with select and drag from the VAF plot.
                    # We don't need to set width here, it's set by args to output$TopOne <- renderPlot(...)
                    plotOutput("TopOne", dblclick = "topOnePlotDblclick", width = "100%", height = "7%"),
                    
                    # VAF plot.
                    # We don't use a double click event here nay more so remove..
                    #plotOutput("targetVcf", click = "targetVcfPlotDblclick", brush = "targetVcfPlotBrush", dblclick = "targetVcfPlotDblclick"),
                    plotOutput("targetVcf", dblclick = "targetVcfPlotDblclick", brush = "targetVcfPlotBrush"),
                    
                    
                    # A pre text area where the depth at a given VAF (specified by a user mouse click) is displayed.
                    # Directly below if is displayed the 'Select VCF' button to select the VCF file in question.
                    # For alignment purposes all page elements are arranged as part of a table from here on.
                    HTML(paste0('<div id="selectDiv" align="center" class="outer" style="background-image:url(\'selectVcf_at_start.png\'); background-repeat: no-repeat; background-position: 50% 20%; background-size: 50%;">')),
                    
                    HTML(paste0(
                      tags$div(id="info", 
                               class="shiny-html-output",
                               align="center",
                               style="font-family:monospace !important; font-size:80%; padding: 8px; background-color:transparent; colour: #333333")), ""),
                    
                    HTML(
                      paste0(
                        '<div align="center">',
                        '<div style="display:inline-block; padding-right: 10px; bottom: 15px">',
                        tags$button(id="evReactiveButton", 
                                    type="button", 
                                    class="btn btn-default action-button shiny-bound-input",
                                    style="border: 1px solid silver;padding:1px; font-size:95%;",
                                    HTML('Select VCF')),
                        '</div>',
                        '<div style="display:inline-block; padding-right: 10px; bottom: 15px">',
                        tags$button(id="writeSubsettedVcf", 
                                    type="button", 
                                    class="btn btn-default action-button InvisableAtStart",
                                    style="border: 1px solid silver;padding:1px; font-size:95%",
                                    HTML('Write VCF')),
                        '</div>',
                        '<div style="display:inline-block; padding-right: 86px; bottom: 15px">',
                        tags$button(id="saveAsSvg", 
                                    type="button", 
                                    class="btn btn-default action-button InvisableAtStart",
                                    style="border: 1px solid silver;padding:1px; font-size:95%",
                                    HTML('Save high res. PNG')),
                        '</div>')
                    ),
                    
                    HTML("<div class='InvisableAtStart' style='display:inline-block; padding-right: 10px;'>"),
                    numericInput("minF", "min freq", value = 0, min = 0, max = 1, step = 0.01, width = "80px"),
                    HTML("</div>"),
                    HTML("<div class='InvisableAtStart' style='display:inline-block; padding-left: 10px;'>"),
                    numericInput("maxF", "max freq", value = 0, min = 0, max = 1, step = 0.01, width = "80px"),
                    HTML("</div>"),
                    HTML("<div class='InvisableAtStart' style='display:inline-block; padding-left: 10px;'>"),
                    selectizeInput("geneList", "select gene",choices = c("Please load a VCF first"), selected = NULL, multiple = FALSE, options = list(), width = "120px"),                    
                    HTML("</div>"),            
                    HTML(paste0('<div class="InvisableAtStart" style="padding:2px" align="center">')),                      
                    radioButtons("afWindow", "Inset window function:",
                                 c("Signatures" = "SigC",
                                   "Mutational Processes" = "MutP",
                                   "Neutrality" = "Nevol",
                                   "Protein" = "PepS",
                                   "Candidate Filters" = "Filt"),
                                 inline=TRUE),
                    HTML(paste0('</div>')),
                    HTML("</div></div>"),          
                    # Close off div with background.          
                    HTML(paste0('</div>')),
                    width = 9
                  )
    )
  )
)


server <- shinyServer(function(input, output, session) {
  
  # Observe a brush stroke event.
  # Update the max and min freq. numeric input values accordingly.
  observeEvent(input$targetVcfPlotBrush, {    
    updateNumericInput(session, "minF", value = round(input$targetVcfPlotBrush$xmin,digits=2))
    updateNumericInput(session, "maxF", value = round(input$targetVcfPlotBrush$xmax,digits=2))   
  })
  
  observeEvent(input$topOnePlotDblclick$x, {
    # Check first to make sure TopOne window is active and signatures secondary analysis set.  
    # Return NULL otherwise.
    if(!(input$minF<input$maxF) || (input$afWindow != "SigC"))
    {
      return(NULL)
    }
    
    # TODO!!!!
    # We should really get the start of these bars on the plot from the ggplot object..
    # We could access this here fro the global, TopOneGgplot.
    # This approximation should do for now however..
    totalLength = 65
    sigNames = colnames(plot_data$vcfFile$cosmicSignatures)
    totalNumCosmicV3sigs = length(sigNames)  
    barLength = totalLength/totalNumCosmicV3sigs
    v3barNum = round(input$topOnePlotDblclick$x/barLength)
    httpLink = paste0("https://cancer.sanger.ac.uk/cosmic/signatures/SBS/", sigNames[v3barNum], ".tt")
    
    # Open the link to this signature description at Sanger.
    browseURL(httpLink) 
  })
  
  observeEvent(input$saveAsSvg, {
    
    timeStamp = format(Sys.time(), '%b_%d_%Y_%Hh%Mm%Ss')
    
    outputFile = paste0(gsub(pattern = "\\.vcf$", ".", plot_data$vcfFile$pathname), timeStamp, ".vcfView.png")
    
    # TODO!!!! URHERE..
    # Can't save as SVG
    # Sort out install issues with svglite and fix up below to save as SVG..
    
    if(!is.null(TopOneGgplot))
    {
      
      showModal(modalDialog(
        title = "Saving inset as PNG.",
        paste0("The inset window has been saved as '", outputFile,"'\n","(please close the inset window and resave if you wish to save the contents of the main window)."),
        easyClose = TRUE
      ))
      
      #  ggsave(file=outputFile, plot=TopOneGgplot, height = 12.5, width = 25, dpi = 64, type = "cairo")
      # note... can't get ggsave to work properly, giving up..
      
      # High resolution png.
      png(filename=outputFile, height = 800, width = 2400, res = 300)
      plot(TopOneGgplot)
      dev.off()
      
      
      # Unfortunately ggsave doesn't talk in terms of pixels so we need to convert here.
      #
      # Our main plot is height = 200, width = 400 pixels.
      # To convert from dpi, (https://community.rstudio.com/t/ggsave-aspect-ratio-whitespace-use-case-favicon-for-blogdown/6793)
      #
      # height = 2, width = 2, dpi = 16 to get a 32x32 pixel image.
      # So for our image,
      # height = 12.5, width = 25
      # note... can't get ggsave to work properly, giving up..
      # This section of comments now redundant..
      
      
    }
    else if(!is.null(targetVcfGgplot))
    {
      showModal(modalDialog(
        title = "Saving main window as PNG.",
        paste0("The main window has been saved as '", outputFile,"'."),
        easyClose = TRUE
      ))
      
      # In pixels, main plot is, height = 400, width = 930
      # Converting this gives us,
      #   ggsave(file=outputFile, plot=targetVcfGgplot,  height = 25, width = 58.125, dpi = 64)
      # ggsave(file=outputFile, plot=targetVcfGgplot,  height = 25, width = 49.9, dpi = 64)
      # note... can't get ggsave to work properly, giving up..
      
      # High resolution png.
      png(filename=outputFile, height = 1632, width = 3720, res = 300)
      plot(targetVcfGgplot)
      dev.off()
      
    }
    else
    {
      showModal(modalDialog(
        title = "No plot available to save.",
        paste0("Please select a VCF and create a plot before saving it."),
        easyClose = TRUE
      ))
      return(NULL)
    }
  })    
  
  observeEvent(input$writeSubsettedVcf, {
    
    # Whats left post filtering of vcf.
    result = filterResult()  
    
    if(is.null(plot_data))
      return(NULL)
    
    if(is.null(plot_data$vcfAnnotateObject))
      return(NULL)
    
    
    if(is.null(result))
      return(NULL)
    
    timeStamp = format(Sys.time(), '%b_%d_%Y_%Hh%Mm%Ss')
    
    outputFile = paste0(gsub(pattern = "\\.vcf$", ".", plot_data$vcfFile$pathname), timeStamp, ".vcfView.vcf")
    
    # If the user has selected a region for the 'TopOne' window then further subset
    # the file to be saved according to this selection.
    
    afFilterString = ""
    
    if(!is.null(input$targetVcfPlotBrush) && !is.null(plot_data) && !is.null(plot_data$snvIdx))
    {
      
      writeVcf(plot_data$vcfAnnotateObject[result$OIDX[plot_data$snvIdx]], outputFile)
      
      afFilterString= paste0("# Candidates variants were subsetted between allele frequency range, [", input$targetVcfPlotBrush$xmin, ", ",  input$targetVcfPlotBrush$xmax, ")")
    }
    
    else
    {
      writeVcf(plot_data$vcfAnnotateObject[result$OIDX], outputFile)
    }
    
    line=paste0("\n# vcfView post processing, ", timeStamp, ":\n#\n")
    
    if(length(input$activeFilters) > 0)
    {
      line=c(line, "# vcfView has removed candidate variants which failed the following filters,")
      line = c(line, gsub("^","# ",as.character(input$activeFilters)))
    }
    
    line = c(line, "\n# vcfView set the following thresholds,\n",
             paste0("# TLOD = ", as.character(input$tlodThresh)),
             paste0("# Tumour depth at variant locus = ", as.character(input$depthThresh)),
             paste0("# Normal depth at variant locus = ", as.character(input$depthThreshInNormal)),
             afFilterString)
    
    
    # Append the comments recording the modifications made to the VCF.
    write(line,file=outputFile,append=TRUE)
    
    
    showModal(modalDialog(
      title = "Saving subsetted VCF.",
      paste0("The subsetted VCF has been saved as '", outputFile,"'."),
      easyClose = TRUE
    ))
  })
  
  
  targetvcf <- eventReactive(input$evReactiveButton, {
    
    if(input$evReactiveButton<1)
      return(NULL)
    
    csvPath = file.choose()
    
    # TODO!!!! csvPath is a bit of a misnomer. Its the path to the VCF file. Update this variable name.
    if(is.null(csvPath))
      return(NULL)
    
    vcfFile = processVcf(csvPath, session)
    
    return(vcfFile);
    
  })
  
  # Set the tumour and normal sample names on the right panel.
  output$tumourNameTag <- renderText({ 
    vcfFile = targetvcf()
    return(paste0("<div style='font-size: 70%'>Tumour: ",vcfFile$tumourSampleName,"</div>")) })
  
  output$normalNameTag <- renderText({ 
    vcfFile = targetvcf()
    return(paste0("<div style='font-size: 70%'>Normal: ",vcfFile$normalSampleName,"</div>")) })  
  
  # Whats left after we have filtered the VCF according to option user has chosen from the UI.
  filterResult <- eventReactive(
    {
      input$activeFilters
      input$depthThresh
      input$depthThreshInNormal
      input$showPassedVariants
      input$tlodThresh
      input$afInNormal
      1
    },
    {      
      vcfFile = targetvcf()
      
      # Now we have out VCF grep out any filters we are not interested in.
      # If no filters have been set on the UI that's OK, we'll just end up with everything in the histogram.
      # Remember, we only want entries which contain all the contents of filter_vect,
      # nothing more, nothing less.
      # (don't be tempted to do something like result[grepl(filter_vect[1]|filter_vect[2] etc.., result$FILTER),])
      # filter_vect=c("PASS")
      # vcfFile$filters = ""
      result=vcfFile$data
      
      if(length(input$showPassedVariants) == 0)
      {
        return(NULL)
      }     
      
      # TODO URHERE!!!!
      # Change this to "display passed/failed variants"
      # Have we been asked to include 'PASS'ed variants in our plot?
      # (ie., not just those that have been filtered out)
      #    if(!is.null(input$includePassedVariants)) {
      #      # Filter them out if we've been asked not to display PASS variants.
      #      if(!input$includePassedVariants)
      #      {
      #        result = result[result$FILTER != "PASS",]
      #      }
      #    }
      
      
      secondPassLogicalVector = rep(FALSE,nrow(result))
      
      # Have we been asked remove all variants that failed a filter?
      # (ie. has the user clicked 'select all', setting the FILTER to 'SET_ALL'?)
      if(!is.null(input$activeFilters) && input$activeFilters == "SET_ALL" && input$showPassedVariants != "FailFOnly") {
        # Only include variants that PASSed (did not fail any filters).
        # If PASS variants have already been filtered out by the previous filter so be it.
        # We will end up with an empty result in that case.
        secondPassLogicalVector = (result$FILTER != "PASS")
        
      }
      
      # Have we been asked to display variants that failed any/all filters?
      # (ie. has the user clicked 'select all', setting the FILTER to 'SET_ALL' and chosen 'failed selected filters but passed all others'?)
      if(!is.null(input$activeFilters) && input$activeFilters == "SET_ALL" && input$showPassedVariants == "FailFOnly") {
        # Only include variants that failed.
        # (ie., exclude secondPassLogicalVector to exclude variants annotated as PASS in vcf.)
        
        secondPassLogicalVector = (result$FILTER == "PASS")
      }
      
      # Have we been asked to display variants that failed any/all filters
      # but yet the user has cleared all filters?
      # That is not possible.. exclude all entries as a result.
      
      if((is.null(input$activeFilters) || input$activeFilters == "CLEAR_ALL") && input$showPassedVariants == "FailFOnly") {
        # Only include variants that failed.
        # (ie., exclude secondPassLogicalVector to exclude variants annotated as PASS in vcf.)
        secondPassLogicalVector = rep(TRUE,nrow(vcfFile$data))
      }
      
      # Have we been requested to filter out some candidate variants that have failed selected filters?
      # Skip this filtering step if user has cleared all filters (either individually, or by selecting "CLEAR_ALL".
      # Also skip it if the user has selected SET_ALL as we dealt with this above in the
      # if(...input$activeFilters == "SET_ALL") statement.
      if(!is.null(input$activeFilters) && input$activeFilters != "CLEAR_ALL" && input$activeFilters != "SET_ALL") {
        # Have we been asked to show the variants that only failed selected filters?
        if(input$showPassedVariants == "FailFOnly")
        {
          grepString=paste("^",paste(input$activeFilters,collapse="$|^"),"$", sep = "")
          
          firstPass=lapply(strsplit(result$FILTER,';'), grepl, pattern=grepString)
          
          # Exclude entry unless it has failed all of the selected filteres.
          secondPassLogicalVector = !(lapply(firstPass,all) == TRUE)
        }
        else
        {
          grepString=paste("^",paste(input$activeFilters,collapse="$|^"),"$", sep = "")
          
          firstPass=lapply(strsplit(result$FILTER,';'), grepl, pattern=grepString)
          
          # If the entry has failed any of the selected filteres, exclude it.
          secondPassLogicalVector = (lapply(firstPass,any) == TRUE)
        }
      }
      
      # Have we been asked to show the variants that passed or failed the filters?
      if(input$showPassedVariants == "FailF")
      {
        secondPassLogicalVector = !secondPassLogicalVector
      }
      
      # Start your depth and tlod filtering.
      
      if(!is.null(input$tlodThresh)) {
        # Record anything that has failed either filter list or the T_LOD threshold.
        secondPassLogicalVector = secondPassLogicalVector | (result$T_LOD <= input$tlodThresh)     
      }    
      
      if(!is.null(input$depthThresh)) {
        # Record anything that has failed either filter list or tumour depth threshold.        
        secondPassLogicalVector = secondPassLogicalVector | (result$DP <= input$depthThresh)
      }
      
      if(!is.null(input$depthThreshInNormal)) {
        # Record anything that has failed either filter list or normal depth threshold.        
        secondPassLogicalVector = secondPassLogicalVector | (result$DP_IN_NORM <= input$depthThreshInNormal)
      }
      
      # If there is an AF ceiling set in the normal sample
      # then subset accordingly.
      if(!is.null(input$afInNormal) & input$afInNormal > 0) {
        # Record anything that has failed either filter list or the depth threshold.        
        secondPassLogicalVector = secondPassLogicalVector | (result$AF_IN_NORM >= input$afInNormal)
      }
      
      
      # Is anything is left after all filtered variants removed?
      result = result[!secondPassLogicalVector,]
      
      
      # If nothing is returned with this filter combination, ignore & return null.
      if(!nrow(result)) {
        return(NULL)
      }
      
      return(result)
    }) 
  
  output$TopOne <- renderPlot({
    
    vcfFile = targetvcf()
    
    # Whats left post filtering of vcf.
    result = filterResult()
    
    # The following events trigger this function.. 
    input$minF
    input$maxF
    
    
    if(is.na(input$minF) || is.na(input$maxF))
    {
      return(NULL)
    }
    
    xMin = input$minF
    xMax = input$maxF
    
    # Return NULL if no vaild range has been specified.
    if(!(xMin<xMax))
    {
      return(NULL)
    }
    
    if(is.null(vcfFile))
    {
      return(NULL) 
    }  
    
    # Update previous log and move on with processing.
    previousMinF <<- input$minF
    previousMaxF <<- input$maxF 
    
    
    # Now with the max and min frequency ranges are specified, handle the event as
    # specified by the afWindow input.
    if(input$afWindow == "SigC")
    {
      # Plot Cosmic signatures contributions.
      return(cosmicComboBreakdown(result,xMin,xMax))
    }
    else if(input$afWindow == "MutP")
    {
      # Plot mutational Processes contributions.
      return(frange(result,xMin,xMax))
    }
    else if(input$afWindow == "Nevol")
    {
      # Check for neutral evolution within this range..
      # (TODO!!!! put in checks to ensure range make sence with plody assumptions etc..)
      return(nevolution(result, xMin,xMax))
    }
    else if(input$afWindow == "PepS")
    {
      
      # Plot mutational Processes contributions.
      if(is.null(input$geneList))
      {
        # Alert user to select a gene of interest before calling this.
        showModal(modalDialog(
          title = "No gene of interest selected.",
          paste0("Please select a gene of interest from one of the options in the 'select gene' input."),
          easyClose = TRUE
        ))
        return(NULL)
      } 
      return(lpoppep(session, result, xMin,xMax, input$geneList))
    }   
    else
    {
      # TODO!!!! put remaining logic in this function into a separate function and call it from here.
      # Tidy this up!
    }
    
    # Otherwise just plot the filters..
    
    candidateVariantsFilterField = result[result$AF>xMin & result$AF<=xMax,"FILTER"]
    
    numCandidateVariants = length(candidateVariantsFilterField)
    
    percentOfPassedCandidates = 100 * (sum(candidateVariantsFilterField == "PASS") / numCandidateVariants)
    
    filterAggregation = result[result$AF>xMin & result$AF<=xMax,"FILTER"]
    
    failedCandidateVariantsFilterField = candidateVariantsFilterField[candidateVariantsFilterField != "PASS"]
    
    numFailedCandidateVariants = length(failedCandidateVariantsFilterField)
    
    filterAggregation=failedCandidateVariantsFilterField # TODO!!!! update 
    
    filterAggregationUl = unlist(strsplit(filterAggregation,';'))
    
    if(length(filterAggregationUl) > 0)
    {
      filterAggregationDf = data.frame(filter_types=filterAggregationUl)
      p = ggplot(filterAggregationDf, aes(x = filter_types, stat(count)/numFailedCandidateVariants)) + geom_bar() + scale_y_continuous(limits = c(0,1)) + ggtitle(paste0("[", signif(xMin,6), "," , signif(xMax,6), "), ", numCandidateVariants, " candidates, ", signif(percentOfPassedCandidates,6), "% PASS")) +
                labs(y="Failed fraction", x = "Filter types") +  theme(axis.text.x=element_text(family="serif", angle=45,hjust=1,vjust=0.95, colour="black", size = 12))
  # Record this plot in TopOneGgplot incase the user wants to save it as svg.
  TopOneGgplot <<- p
  return(p)

    }
    
    else
      return(NULL)
    
  }, height = 200, width = 600, bg="transparent")
  
  output$targetVcf <- renderPlot({ 
    
    vcfFile = targetvcf()
    
    # Whats left post filtering of vcf.
    result = filterResult()
    
    
    # Check to see if any data was returned in 'result'.
    # If not inform the user by constructing a null plot.
    if(is.null(result))
    {
      nullPlot = data.frame(density=c(0), AF=c(0))
      p=ggplot(data=nullPlot, aes(x=AF, y=density)) +
        ggtitle(vcfFile$pathname) +
        theme_bw() +
        theme(axis.line = element_line(size=1, colour = "black"),
              panel.grid.major = element_line(colour = "#d3d3d3"),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(), panel.background = element_blank(),
              plot.title = element_text(size = 8, family = "Tahoma", face = "italic"),
              text=element_text(family="Tahoma"),
              axis.text.x=element_text(colour="black", size = 9),
              axis.text.y=element_text(colour="black", size = 9)) +
        xlim(0, 1) +
        ylim(0, 1)
      p = p + annotate("text",x=0,y=.8,hjust=0,label="No variants found.", size=10, fontface = "bold")
      
      # Clear plot data, there is nothing to save to SVG here.
      targetVcfGgplot <<- NULL
      
      return(p)
    }
    
    # Btw., if you leave out the 'y = ..density..' below it will default to y=counts..
    
    # TODO!!!! make this user configurable by adding something to the UI.
    num_of_bins = 500;
    
    # Create a interpolation pallet to reflect this.
    colours = maxMedianBinDepth(num_of_bins, result)
    
    # Now we've got a plot update gene list annotation.
    # Make sure the gene list provided to the user reflects
    # only genes that are displayed in the plot so as not to clutter.
    
    # See if the max, min freq limits are set.
    # If they are only show genes within these limits.
    
    # Ignore if no valid range previously specified.
    if((previousMinF<previousMaxF) && previousMaxF != 0)
    {
      #
      # Subset out all the SNVs within required allele freq. range.
      snvIdx = which(result$AF > previousMinF &
                       result$AF <= previousMaxF)
      # pull out the annotation for this this subset using the oidx. 
      annSubsetIdx = !is.na(match(geneAnnotateDf$QUERYID, result$OIDX[snvIdx]))
      
      # However, if there is nothing found within this range,
      # pull out the annotation for this result subset using the oidx
      # without any freq range.
      # A 'not found' plot will be displayed to the user.
      if(!any(annSubsetIdx))
      {
        annSubsetIdx = !is.na(match(geneAnnotateDf$QUERYID, result$OIDX))
      }
      
    }
    else
    {
      # pull out the annotation for this result subset using the oidx
      # without any freq range. 
      annSubsetIdx = !is.na(match(geneAnnotateDf$QUERYID, result$OIDX))
    }
    
    # This is the subset we're interested in.
    annaSubset = geneAnnotateDf[annSubsetIdx,]
    
    # Update selectize options with genes that are referenced in this VCF.
    gChoices = as.character(unique(geneAnnotateDf[annSubsetIdx,]$GENE_SYMBOL))
    # currentGeneSelected = input$geneList
    
    if(!is.na(gChoices) && (length(gChoices) > 0))
    {
      updateGeneSelectionList(session, gChoices, gChoices[1])
    }
    else
    {
      updateGeneSelectionList(session, NA, NA)
    }
    
    g=ggplot(result, aes(x = AF, y = ..density..)) + geom_histogram(bins=num_of_bins, fill = colours)
    
    # plot_data <- list(pathname=NULL, sampleNames=NULL, selectedSample=NULL, gg_b_object=NULL, result=NULL, dp_breaks=NULL)
    
    # Plot.
    bla = g + geom_density(alpha=0.4) +
      ggtitle(vcfFile$pathname) +
      labs(y= "Density", x = "Allele frequency") +
      theme_bw() +
      theme(axis.line = element_line(size=1, colour = "black"),
            panel.grid.major = element_line(colour = "#d3d3d3"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(), panel.background = element_blank(),
            plot.title = element_text(size = 8, family = "Tahoma", face = "italic"),
            text=element_text(family="Tahoma"),
            axis.text.x=element_text(colour="black", size = 9),
            axis.text.y=element_text(colour="black", size = 9))
    
    # Uodate targetVcfGgplot in case user wants to prprint this out in svg format.
    targetVcfGgplot <<- bla
    
    return(bla)
    
  }, height = 408, width = 930)
  
  output$info <- renderText({
    
    vcfFile = targetvcf()
    
    if(is.null(vcfFile))
      return("Double click on plot for read depth at a given VAF...")
    
    if(is.null(input$targetVcfPlotDblclick$x))
      return("Double click on plot for read depth at a given VAF...")
    
    # Clear the Allele frequency range
    # (when the TopOne window is cleared by a user single click).
    # updateNumericInput(session, inputId="minF", value = 0)
    # updateNumericInput(session, inputId="maxF", value = 0)
    
    updateNumericInput(session, "minF", value = 0)
    updateNumericInput(session, "maxF", value = 0)
    
    
    # Clear plot data, there is nothing to save to SVG here.
    TopOneGgplot <<- NULL
    
    # Run through the breaks finding out which one this mouse click has landed in.
    
    # You might think we could find this out with something like, (since x ranges here from 0-1)
    # input$targetVcfPlotDblclick$x * <number of bins>
    # However with gg_plot the first bin is generally half this size of all the other bins (to center bins etc.).
    # Also, you may request a number of bins up front, but there is no guarantee you will get it (I think..)
    # I can find no guarantee of where gg_plot puts its breaks, other than looking directly at the gg_plot object.
    # Therefore I do it this way...
    
    i = 1
    while (i < 500) {
      if(input$targetVcfPlotDblclick$x < plot_data$gg_b_object$xmin[i])
      {
        i = i-1
        break
      }
      i = i+1
    }
    
    # TODO!!!!
    #
    # A better way to find out which bin this AF falls in (above) may be as follows,
    #
    # Start by putting it in last bin to start with and reset it if we find it's located elsewhere.
    #    i=length(gg_bd$xmin)
    #    
    #    if(aF < gg_bd$xmin[length(gg_bd$xmin)])
    #    {
    #      i = max(which(aF > gg_bd$xmin))
    #    }
    
    
    if(i<1)
      return("Click on plot for read depth at a given VAF...")
    
    num_bins_in = i
    
    # Pull the avg. depth out of the dp_breaks vector.
    # We worked this out before when creating the colour gradient.
    # return(paste0("VAF = ", input$targetVcfPlotDblclick$x, ", avg. depth = ", plot_data$dp_breaks[num_bins_in]))
    
    # TODO!!!! check the limits on the breaks here, is it xmin >= B < xmax or is it xmin > B =< xmax ???
    # I think it's the former but check it out....
    
    rtnString = paste0("VAF = [", signif(plot_data$gg_b_object$xmin[i],6), "," , signif(plot_data$gg_b_object$xmax[i],6), "), median depth = ", signif(plot_data$dp_breaks[num_bins_in],4), ", ", plot_data$gg_b_object$count[i], " variant")
    
    # Get our grammer right..
    if(plot_data$gg_b_object$count[i]>1)
    {
      rtnString = paste0(rtnString,"s ")
    }
    else
    {
      rtnString = paste0(rtnString," ")
    }
    
    ##################################################
    # TODO!!!!/ Done???? Do we want to keep this?
    #  Do we want to annotate these variants?
    # We could do something like,
    
    # Subset out all the SNVs within required allele freq. range.
    result = filterResult()  
    varIdx = which(result$AF > plot_data$gg_b_object$xmin[i] &
                     result$AF <= plot_data$gg_b_object$xmax[i])
    
    # pull out the annotation for this this subset using the oidx. 
    annSubsetIdx = !is.na(match(geneAnnotateDf$QUERYID, result$OIDX[varIdx]))
    
    # This is the subset we're interested in.
    annaSubset = geneAnnotateDf[annSubsetIdx,]
    
    # Only pull out what we've got gene symbol annotation for..
    gsIdx=nchar(as.character(annaSubset$GENE_SYMBOL)) > 0
    gsIdx[is.na(gsIdx)] = FALSE
    
    
    glist = unique(paste0("<option>(",annaSubset$GENE_SYMBOL[gsIdx]," [DP=",annaSubset$DP[gsIdx],"])</option>"))
    
    
    # rtnString = paste0(rtnString,'<div style="cursor: pointer; display:inline-block; color: blue; text-decoration: underline"><select style="cursor: pointer; border:0px; outline:0px; -webkit-appearance: none !important;">',paste(glist,collapse = ''),'</select>',length(as.character(annaSubset$GENE_SYMBOL)),',',nchar(as.character(annaSubset$GENE_SYMBOL)),',',lengths(as.character(annaSubset$GENE_SYMBOL)),'</div>')
    
    # Append the annotation.
    if(length(as.character(annaSubset$GENE_SYMBOL)) > 0)
    {
      rtnString = paste0(rtnString,'<div style="font-family:monospace !important; cursor: pointer; display:inline-block; color: blue; text-decoration: underline"><select style="background-color:transparent; font-family:monospace !important; font-size: 80%; cursor: pointer; border:0px; outline:0px; -webkit-appearance: none !important;">',paste(glist,collapse = ''),'</select></div>')
      
    }
    ##################################################
    
    return(rtnString)
    
  })
  
  output$listFilters <- renderUI({
    vcfFile = targetvcf()
    
    if(is.null(vcfFile))
      return(NULL) 
    
    req(targetvcf)
    
    htmlString="<div align='center' class='outer InvisableAtStart' style='padding:5px'>"
    
    
    htmlString=paste0(htmlString,radioButtons("showPassedVariants", "Show variants which:",
                                              c("Passed selected filters." = "PassF",
                                                "Failed one or more filters" = "FailF",
                                                "Failed selected filters but passed all others." = "FailFOnly"),
                                              inline=FALSE))
    
    htmlString=paste0(htmlString,"<div align='left' style='padding:0px;'>")
    
    
    
    htmlString=paste0(htmlString,"<div id='listFilters' class='shiny-html-output shiny-bound-output'><label class='control-label' for='listFilters'>Filter settings:</label><div id='activeFilters' style='margin-bottom:3px; margin-top:-11px;' class='form-group shiny-input-checkboxgroup shiny-input-container shiny-bound-input'><div class='shiny-options-group'>")
    
    
    for (i in vcfFile$filterTypes) {
      descr=vcfFile$filterDescriptions[vcfFile$filterDescriptions$ID == i,]$Description
      htmlString=paste0(htmlString,filterCheckboxHtml('activeFilters',i,descr))
    }
    
    htmlString=paste0(htmlString,"</div>")
    
    
    # Put in a select / clear all divs to shortcut checkbox selection.
    htmlString=paste0(htmlString,"<div style='cursor: pointer; display:inline-block; font-size: x-small; color: blue; text-decoration: underline' onclick='shinyToggleGroup(\"activeFilters\",true)'>select all</div><div style='display:inline-block; font-size: x-small; color: blue;'>&nbsp;|&nbsp;</div><div style='cursor: pointer; display:inline-block; font-size: x-small; color: blue; text-decoration: underline' onclick='shinyToggleGroup(\"activeFilters\",false);'>clear all</div>")
    
    
    # Leave it to the UI section to close off this div.
    # This is because (for some reason) we can't put a slider in her by pasting it onto the html string..
    htmlString=paste0(htmlString,"</div>")
    
    # Now we have the data we know the limits we want to set on our sliders..
    # Nevermind all this messing about with javascript..
    # updateSliderInput() to the rescue...
    #document.getElementById("tlodThresh_slider")..children[0].children[1].children[0].getAttribute("...");
    #document.getElementById("tlodThresh_slider").getElementsByTagName("span").getAttribute("...");
    
    minDepth = min(vcfFile$data$DP)
    maxDepth = max(vcfFile$data$DP)
    
    # midDepth = minDepth + (maxDepth - minDepth)/2
    # Actually better to leave the mid and min depth the same for now.
    # The user can push the bar up as they want.
    midDepth = minDepth
    
    
    # TODO!!!! Take out slider inputs, they're not fit for purpose here.
    # updateSliderInput(session, "depthThresh", "Depth threshold:", value = midDepth, min = minDepth, max = maxDepth)
    
    # We need to round the min and max values to the nearest 2 places of decimal to it fits with our step below.
    minDepth = round(minDepth,digits=2)
    maxDepth = round(maxDepth,digits=2)
    
    updateNumericInput(session, "depthThresh", value = minDepth, min = 3, max = maxDepth, step = 1)  
    minTlod = min(vcfFile$data$T_LOD)
    maxTlod = max(vcfFile$data$T_LOD)   
    
    # We need to round the min and max values to the nearest 2 places of decimal to it fits with our step below.
    minTlod = round(minTlod,digits=2)
    maxTlod = round(maxTlod,digits=2)
    
    updateNumericInput(session, "tlodThresh", value = 6.3, min = minTlod, max = maxTlod, step = 0.01)
    
    # Now we've got a plot update gene list the annotation.
    
    # pull out the annotation for this this subset using the oidx. 
    annSubsetIdx = !is.na(match(geneAnnotateDf$QUERYID, vcfFile$data$OIDX))
    
    # Update selectize options with genes that are referenced in this VCF.
    updateSelectizeInput(session,"geneList", choices = as.character(unique(geneAnnotateDf[annSubsetIdx,]$GENE_SYMBOL)))
    
    
    # When the div is loaded expose the sliders.
    # No point in having sliders etc. hanging around unless they are connected to something....
    htmlString=paste0(htmlString,'<script type="text/javascript">
                                  revealElements("InvisableAtStart");
                                  document.getElementById("selectDiv").style.backgroundImage="none";
                                  </script>')
    
    return(tags$html(HTML(htmlString)))
  })
  
  
})



shinyApp(ui=ui,server=server) 
