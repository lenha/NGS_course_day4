---
title: "Day 4 practicals"
author: "Lena Hansson"
date: "1 November 2016"
output: 
  html_document:
    toc: yes
    toc_depth: 5
    toc_float: true
---

<style type="text/css">

div.quest {
  background: #ECF8FF;
  border-left: 10px solid #3989CB;
  margin: 1.5em 10px;
  padding: 0.5em 10px;
  font-size: 14px;
}

div.diy {
  background: #FF7F50;
  border-left: 10px solid #7f0000;
  margin: 1.5em 10px;
  padding: 0.5em 10px;
  font-size: 14px;
}

div.assignment {
  background: #32CD32;
  border-left: 10px solid #556B2F;
  margin: 1.5em 10px;
  padding: 0.5em 10px;
  font-size: 14px;
}

h1 { font-size: 18px }
h2 { font-size: 16px }
h3 { font-size: 15px }


</style>

```{r setupDoc, echo=FALSE}
knitr::opts_chunk$set(echo = TRUE)
questNr <- diyNr <- assNr <- 0 # set all 3 counters to 0
```

Today we are going to  
1) be investigating the results from the alignment and the differential expression analysis that you went through on the previous occassion - do these make sense? do we trust them?  
2) learn how to use RMarkdown to more easily generate reproducible results  
3) learn how to use ggplot2 to make, and annotate, publication ready plots (PDF)  
 - PCA (incl gene loadings)  
 - Volcano plot  
4) Use IGV and Gviz to visualise coverage of genes


The idea is that you will be given most of the R code, but unlike the previous day you will be fored to "copy and paste", as well as think about what the code is doing.

During the course of the day we will go through the exercises in this document, and there will be three types of "exercises":  
<div class="diy">*Try it yourself*: Copy and paste the example code, but make a minor change, in order to test how it works</div>
<div class="quest">*Question*: These are discussion points which require you to think about a point and argue for/against a specific viewpoint</div>
<div class="assignment">*Assignment*: These specifies the plots that you should hand in as todays assignment</div>


# Data description

(taken from previous day!)

We will continue looking at the data from the previous day, eg *YAP and TAZ control peripheral myelination and the expression of laminin receptors in Schwann cells*. Poitelon et al. Nature Neurosci. 2016. (http://www.nature.com/neuro/journal/v19/n7/abs/nn.4316.html) In the experiments performed in this study, YAP and TAZ were knocked-down in Schwann cells to study myelination, using the sciatic nerve in mice as a model. The material for RNA-seq was collected from 2 conditions (wt and YAP(kd)TAZ(kd)), each in 3 biological replicates:

Accession |	Condition |	Replicate
----------|-----------|--------
SRR3222409 |	KO |	1
SRR3222410 |	KO |	2
SRR3222411 |	KO |	3
SRR3222412 |	WT |	1
SRR3222413 |	WT |	2
SRR3222414 |	WT |	3



# Reproducible research

Use RMarkdown for reproducible results as a way to connect your data with your "results" (eg plots), and discussion points.

## Useful links

RStudio + RMarkdown  
1) https://www.rstudio.com/wp-content/uploads/2015/02/rmarkdown-cheatsheet.pdf  
2) http://rmarkdown.rstudio.com/authoring_basics.html  
3) http://www.rstudio.com/wp-content/uploads/2015/03/rmarkdown-reference.pdf  
R  
1) Quick-R (for example http://www.statmethods.net/stats/index.html and http://www.statmethods.net/management/sorting.html)  
2) http://www.cyclismo.org/tutorial/R/


## An empty HTML document

<div class="diy">*Try it yourself* `r diyNr<-diyNr+1;diyNr`: Generate a R Markdown and take a look at the example (see below for how to)</div>
- Open RStudio
- Make a new R Markdown document (for example by clicking the top left arrow and choosing R Markdown).
- Name the document something use (eg change the title from Untitled to for example Day 4 practical) 
- Choose the HTML format.
Everything you do today will be collected in this document, which means you have have a "log" of all commands necessary to achieve the plots, along with the plots and your comments.  
- Then click on the Knit HTML button on the middle command line above your document, and see what it looks like.

<div class="diy">*Try it yourself* `r diyNr<-diyNr+1;diyNr`: Change the chunk option so that you can see the plot command in the generated html (hint echo)</div>  
<div class="diy">*Try it yourself* `r diyNr<-diyNr+1;diyNr`: Add a header above the plot by adding a # at the beginning of the line (see the useful links for other nice formatting options)</div>  
<div class="diy">*Try it yourself* `r diyNr<-diyNr+1;diyNr`: Re-knit the HTML file and ensure your two changes are included.</div>  

When you have a small simple file like this re-knitting every time is useful, but as you progress to more and more complex, and longer, R Markdown files, you could run your commands directly in RStudio and only re-knit once you want to see how your HTML page looks like. So go to the plot command of your RMarkdown file and press
```
Command and Enter
```
you will now see the plot in the bottom right corner of your RStudio


## Set up your environment

The first you should do is to 'setup' your environment, and this includes to load all required packages, and determine how the output should look like.

In general I find it easier to view "dependencies" by loading all R packages at the top, no matter how far down your document you will use them... Also I tend to add a comment (# to the right) of why I need a specific package, and also the link to the web-page containing the examples on how to use the package. The first time you load a package, it might not exists on your computer, so the below is a good code for testing if the package exists, and if not install it.

**Note** there are, at least, three major sources, or R packages  
* cran  
* bioconductor  
* github  
and they have slighly different installation instructions (for example library vs biocLite)!

```{r, echo=TRUE, warning=FALSE, message=FALSE }
if (!require("RColorBrewer")) {                       # test if the package can be found
  install.packages("RColorBrewer", dependencies = TRUE);  # install the package, if needed
}
library(RColorBrewer);                                # load the package
```

<div class="diy">*Try it yourself* `r diyNr<-diyNr+1;diyNr`: Load all the packages we will be using today (see below)</div>  

The full list of packages:
```{r, echo=TRUE, warning=FALSE, message=FALSE}
library(dplyr);                          # summarise your data
library(ggplot2);                        # ggplot
library(knitr);                          # tables with more options on what they will look like... 
library(limma);                          # voom
library(VennDiagram);                    # https://rstudio-pubs-static.s3.amazonaws.com/13301_6641d73cfac741a59c0a851feb99e98b.html
library(grid);                           # required by make3wayVenn
#library(RColorBrewer);
source("https://bioconductor.org/biocLite.R");   # a "collection" of many bioinformatics packages
library(Gviz);
library(GenomicRanges);                  # used by Gviz
library(BSgenome.Hsapiens.UCSC.hg19);    # used by Gviz
library(BSgenome.Mmusculus.UCSC.mm10);   # used by Gviz
```

Define where RStudio should be looking for additional files
```{r, echo=TRUE}
setwd("/Users/halena/Documents/WABI/Teaching/2016_GU_Seqcourse");
```
<div class="diy">*Try it yourself* `r diyNr<-diyNr+1;diyNr`: Copy the following files from yesterday into a single directory and then set the current working directory to this</div>
1) the tableCounts file (tableCounts_with_location.tab)  
2) results_DE.txt
3) the six *sorted* .bam files
4) the six .bai files for the above .bam files

(they can also be found on ruby on */home/teacher5/rna-seq-course*)

List *all* the files in the current working directory, there are two commands that means the same thing *list.files()* and *dir*.
```{r, echo=TRUE, eval=FALSE}
list.files()            # list ALL files
dir(pattern="*bam")     # list all BAM files
```
You will see again and again that there are multiple ways of achieving the same results in R, you have to figure out which one makes the most sense to you.


# QC your results

Its important to quality control (QC) your results, and this is very easily done in R by "viewing" your data.

## Load in the feature counts

```{r, echo=TRUE}
tableCounts <- read.table(file="tableCounts", sep="\t", header=TRUE, skip=1);
```
<div class="quest">*Question* `r questNr<-questNr+1;questNr`: Why do you think we wrote skip=1 in the above command?</div>

Look at the values in the table
```{r, echo=TRUE}
head(tableCounts);
```

If you are unsure about how a specific command is written you can always call *help* in either R or RStudio
```{r, echo=TRUE}
help(read.table)
```

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: What other options do the read.table function allow us?</div>

In order to look at the data in a better way inside RStudio, although it *wont* work inside your knitr
```{r, echo=TRUE, eval=FALSE}
View(tableCounts)
```

## Look at the feature counts

### Which gene is the most expressed?

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: why wont *max(tableCounts)* work?</div>

Different ways of using max to get the answer you need.
```{r, echo=TRUE}
max(tableCounts[,c("SRR3222409_Aligned.out.sam", "SRR3222410_Aligned.out.sam", "SRR3222411_Aligned.out.sam", "SRR3222412_Aligned.out.sam", "SRR3222413_Aligned.out.sam", "SRR3222414_Aligned.out.sam")]) # will work but is awkward
max(tableCounts[,7:12]) # will work but if you add a column it will break!
# for each gene get the max count
maxCountPerGene <- summarise(group_by(tableCounts, Geneid), maxCount=max(SRR3222409_Aligned.out.sam, SRR3222410_Aligned.out.sam, SRR3222411_Aligned.out.sam, SRR3222412_Aligned.out.sam, SRR3222413_Aligned.out.sam, SRR3222414_Aligned.out.sam));
max(maxCountPerGene$maxCount); # get the max count from this directly...
# use head and order your data
head(maxCountPerGene[order(maxCountPerGene$maxCount, decreasing=TRUE),], n=1);
```

check that it did what you expected it to
```{r}
head(maxCountPerGene);
```

how is the gene expressed across conditions?
```{r}
subset(tableCounts, Geneid %in% "ENSMUSG00000000001"); # test that the max count is indeed right
```

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: which way did you prefer?</div>
<div class="quest">*Question* `r questNr<-questNr+1;questNr`: What is this gene? Have a look at the Ensembl database and see if you recognise it</div>
http://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000056569;r=1:171150711-171161130

### For which sample (control or case) did we see it?
```{r, echo=TRUE}
subset(tableCounts, Geneid %in% "ENSMUSG00000056569");
```

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: Is this an outlier? or is this real?</div>

### Make a data frame of the counts
This will make it easier to access the alignment data later on
```{r}
onlyCounts = tableCounts[,7:12];
rownames(onlyCounts) = tableCounts$Geneid;
```

and if so you could have done
```{r}
onlyCounts[apply(onlyCounts[,-1],1,max)==max(onlyCounts[,-1]),]
```

### How many genes have a zero count across all six samples?
```{r, echo=TRUE}
# how many genes do we have at all?
nrow(maxCountPerGene);
# how many genes have a count of zero?
nrow(subset(maxCountPerGene, maxCount==0))
```

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: Would you include these in your analysis?</div>
<div class="quest">*Question* `r questNr<-questNr+1;questNr`: Where they included in the differential expression analysis?</div>


# Summarise the QC

Now you actually have quite a few plots, so in order to find your way around the HTML page, add a table of content by modifying your header section.

Have a look at how to specify your header section here:
http://rmarkdown.rstudio.com/html_document_format.html 

Once you added the TOC, go ahead and re-knit your document

Did you like the layout of the tableCounts table? No, then lets do something about it
```{r, echo=TRUE}
kable(head(tableCounts, row.names=FALSE));
```

Lets rename the counts columns to something more useful.
```{r, echo=TRUE}
names(tableCounts)[7:12] <- c("KO_1","KO_2", "KO_3", "WT_1", "WT_2", "WT_3");
```

Lets look at 10 genes with actual KO expression levels...
```{r, echo=TRUE}
kable(head(tableCounts[order(tableCounts$KO_1, decreasing=TRUE),], n=10), row.names=FALSE);
```

Now you have both the results (in the form of plots), along with the commands use to generate the plots, and even, in some ways, the data that was the basis for the plots, in the same document. In order to make your results even more re-producable, start adding comments, like figure legends to your plots.



# Volcano plots

In order to make a volcano plot, you need to read in the differential expression results from the previous day, and then you could technically simply plot it, but in general you are interested in visually highlighting the differentially regulated genes.

## Read in the differential results

Simply use read.table
```{r, echo=TRUE}
diffExp <- read.table("results_DE.txt",sep="\t",header=TRUE);
nrow(diffExp);
```

Then on the command line outside RStudio check the number of lines in the file against the above count
```{r, eval=FALSE}
wc -l results_DE.txt
```

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: did these return the same number of lines?</div>

The above wont work as there are lots of ' in the file! this is a general problem with reading files into R, wheras the below will work.
```{r, echo=TRUE}
diffExp <- read.table("results_DE.txt", sep="\t", header=TRUE, quote=""); 
```

### Look at the two knock-outs
```{r, echo=TRUE}
subset(diffExp, mgi_symbol %in% c("Taz","Yap1"))
```

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: Was the knock-out experiment successful?</div>


## Highlight significant genes
```{r, echo=TRUE}
diffExp$threshold = as.factor(diffExp$FDR<=0.01 & abs(diffExp$logFC)>=2)
```

## Construct the plot object
```{r, echo=TRUE}
g <- ggplot(data=diffExp, aes(x=logFC, y=-log10(PValue), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value")
print(g);
```

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: what do you think the error message from the above command means? Is this something we have to worry about?</div>

## Add the mouse symbols

Here we will add mouse symbols on top of the circles for the significantly expressed genes only, as adding it for all 13.000 genes will become unreadable.

```{r, echo=TRUE}
genesSignificant <- diffExp[diffExp$threshold==TRUE,]; # first pull out the significant genes
g + geom_text(data=genesSignificant, aes(label= mgi_symbol), size=1.2, colour="black");
```

The above works as genesSignificant looks exactly the same (names of columns) as diffExp, which is because its a subset of the rows in it.

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: is Taz and Yap among the significant genes?</div>
<div class="diy">*Try it yourself* `r diyNr<-diyNr+1;diyNr`: Change the above plot so that the text size is much larger</div>
<div class="assignment">*Assignment* `r assNr<-assNr+1;assNr`: Subset the volcano plot to only include -log10 pvalues above 5, how many of the significant genes are included in this plot?</div>
<div class="quest">*Question* `r questNr<-questNr+1;questNr`: Why do you think we subset the significant genes to only those with an *abs(logFC)*>2?</div>



# Venn diagrams

## Make a 3-way Venn diagram
```{r, echo=TRUE}
make3wayVenn <- function(A, B, C, categories){
  AB=intersect(x=A, y=B);
  AC=intersect(x=A, y=C);
  BC=intersect(x=B, y=C);
  overlap = intersect(x=AB, y=C);
  
  grid.newpage();
  v <- draw.triple.venn(area1=length(A),area2=length(B),area3=length(C),
                        n12=length(AB), n13=length(AC), n23=length(BC), 
                        n123=length(overlap),
                        category = categories, scaled=FALSE, 
                        print.mode=c('raw','percent'), sigdigs=0,
                        lty = 'blank', 
                        fill = c('skyblue', 'pink1', 'mediumorchid'));
  grid.draw(v);
}
```

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: What do you think the above code does?</div>
<div class="quest">*Question* `r questNr<-questNr+1;questNr`: Do you like it?</div>

```{r}
expressedInKO1 <- subset(tableCounts, KO_1>10);
expressedInKO2 <- subset(tableCounts, KO_2>10);
expressedInKO3 <- subset(tableCounts, KO_3>10);
make3wayVenn(expressedInKO1$Geneid, expressedInKO2$Geneid, expressedInKO3$Geneid,
             categories=c("KO_1","KO_2","KO_3"));
```

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: Is the overlap in expressed genes good?</div>
<div class="quest">*Question* `r questNr<-questNr+1;questNr`: Where would you like to set the cutoff for when a gene is expressed?</div>
<div class="assignment">*Assignment* `r assNr<-assNr+1;assNr`: Make a 3-way venn diagram for the control genes</div>

## More Venn Diagrams

<div class="diy">*Try it yourself* `r diyNr<-diyNr+1;diyNr`: type help("draw") to see a list of all the options from the VennDiagram package</div>

<div class="assignment">*Assignment* `r assNr<-assNr+1;assNr`: Make a 2-way venn diagram for the genes expressed in control vs the genes expressed in KO using the diffExp data frame</div>

Other Venn packages

In my opinion Venn diagrams is one of the plots you are asked for again and again that R does not provide a really good answer for. If you Google for *venn r* you will get a long list of possible packages. I like the one we used here because it gives me the possibility to easily set the layout of it, and to show the percentage as well.



# PCA plots

Take the normalised data (here we use log2), decide the annotations you want to add to the plot
```{r, echo=TRUE}
pseudoCount = log2(onlyCounts + 1);
annotations <- data.frame(Condition=c("KO","KO", "KO", "WT", "WT", "WT"), 
                          Names=names(pseudoCount));
```
Actually run the PCA analsysis
```{r}
pseudoCount.t <- t(pseudoCount);      # transpose the data
pcaObject <- prcomp(pseudoCount.t);   # perform the PCA analysis
```
Add the annotations to the pca results, calculate how much each principal component contributes
```{r}
X <- merge(x = pcaObject$x, y = annotations, 
           by.x="row.names", by.y = "Names", all.x=TRUE);
percentage <- round(pcaObject$sdev / sum(pcaObject$sdev) * 100, 2);
df <- as.data.frame(pcaObject$x);
percentageLabs <- paste( colnames(df), "-", paste( as.character(percentage), "%", sep="") );
```
Show the results as an annotated plot
```{r}
ggplot(X, aes_string(x="PC1", y="PC2", color="Condition"))+
  geom_point() +
  labs(title = "PCA") +
  xlab(percentageLabs[1]) + ylab(percentageLabs[2]);
```

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: do you think the KO and the WT are separated from each other?</div>
<div class="quest">*Question* `r questNr<-questNr+1;questNr`: is the percentage of variance explained by PC1 and PC2 high or low?</div>

## Make a gene loading plot

```{r, echo=TRUE}
showGeneLoadings <- function(pcaObject, N=30, numberOfPCToShow=4){
  gl <- as.data.frame(pcaObject$rotation);
  gl <- merge(x=gl, y=diffExp[,c("ensembl_gene_id","mgi_symbol")], by.x="row.names", by.y="ensembl_gene_id", all.x=TRUE);
  geneLoadings = array(data=0, dim=numberOfPCToShow);
  for(g in 1:numberOfPCToShow){
    p=paste("PC",g,sep="");
    if(p %in% colnames(gl)){
      gl.cur <- gl[, c("Row.names", "mgi_symbol", p)];
      names(gl.cur) = c("Ensembl", "Symbol", "Loading");
      gl.pc <- head(gl.cur[order(abs(gl.cur$Loading), decreasing=TRUE),], n=N); # pick the top N gene loadings
      sumInfluence=round(sum(t(abs(gl.pc$Loading))), digits=4);
      gl.pc$Symbol <- reorder(gl.pc$Symbol, gl.pc$Loading); # needed in order to make ggplot order the data in the plot
      gPlot <- ggplot(gl.pc, aes(x=Loading, y=factor(Symbol)))+
        geom_point(color="Blue", aes(order=Loading))+
        labs(x="Loading",y="Gene",title=sprintf("Top %i genes in %s", N, p));
      print(gPlot);
    }
  }
}
showGeneLoadings(pcaObject=pcaObject)
```

<div class="diy">*Try it yourself* `r diyNr<-diyNr+1;diyNr`: Show only the top 10 gene loadings for the 6 first principle components</div>
<div class="diy">*Try it yourself* `r diyNr<-diyNr+1;diyNr`: what would happen if you asked for the first 10 principal components?</div>
<div class="quest">*Question* `r questNr<-questNr+1;questNr`: how many principal components are there?</div>
<div class="quest">*Question* `r questNr<-questNr+1;questNr`: can you find Taz and Yap among the top gene loading genes? (does this make sense)</div>
<div class="quest">*Question* `r questNr<-questNr+1;questNr`: how does the top gene loading genes correspond to the differentially expressed genes?</div>
<div class="quest">*Question* `r questNr<-questNr+1;questNr`: PC2 seem to consist of genes mostly from a single family, which one and what does it do?</div>


# Visualise coverage of chromosomal neighborhood

## in IGV (standalone tool)

What we need:  
1) IGV - the actual program  
2) genome file (.fa) - "comes with" IGV  
3) sorted .bam file of interest(s)  
4) index of the sorted .bam file of interest(s), eg the .bai files

### Genome file
Inside IGV:
```
Genomes  
Load genomes from server  
Select 'Mouse mm10'  
```

### Load a .bam file
Inside IGV:  
```
File  
Load from file  
Select one of the .bam files - note that the .bam and the .bam.bai have to be located in the same folder
```

Do the above for all six .bam files

### View a gene
Inside IGV:  
```
In the textbox (top middle), type Taz and then click on the Go button  
```

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: Is Taz differentially expressed in your opinion?</div>
<div class="diy">*Try it yourself* `r diyNr<-diyNr+1;diyNr`: View one of the top genes found in the gene loadings of the principal components analysis, is this one differentially expressed?</div>


## Via Gviz in R

https://bioconductor.org/packages/release/bioc/html/Gviz.html  
https://bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.pdf  

### Example from the vignette
Plot CpG islands on chromosome 7 from the human genome hg19 - for a full descriptions of what the commands means, browse through the vignette on your own

```{r, echo=TRUE}
data(cpgIslands);
chr <- as.character(unique(seqnames(cpgIslands)))
gen <- genome(cpgIslands);
atrack <- AnnotationTrack(cpgIslands, name = "CpG");
plotTracks(atrack);
gtrack <- GenomeAxisTrack();
plotTracks(list(gtrack, atrack));
itrack <- IdeogramTrack(genome = gen, chromosome = chr);
plotTracks(list(itrack, gtrack, atrack));
```

#### Add gene information
```{r, echo=TRUE}
data(geneModels);
grtrack <- GeneRegionTrack(geneModels, genome = gen, chromosome = chr, name = "Gene Model")
plotTracks(list(itrack, gtrack, atrack, grtrack));
```

<div class="diy">*Try it yourself* `r diyNr<-diyNr+1;diyNr`: make the chromosome at the bottom instead of at the top</div>

```{r, echo=TRUE}
plotTracks(list(itrack, gtrack, atrack, grtrack), from = 26700000, to = 26750000);
plotTracks(list(itrack, gtrack, atrack, grtrack), extend.left = 0.5, extend.right = 1000000);
```

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: come up with an example in which it would be nice to visualise the up-stream region.</div>


#### Show the bases
```{r, echo=TRUE}
strack <- SequenceTrack(Hsapiens, chromosome = chr);
plotTracks(list(itrack, gtrack, atrack, grtrack, strack), 
           from = 26591822, to = 26591852, cex = 0.8);
```

#### Add gene symbols
```{r, echo=TRUE}
grtrack <- GeneRegionTrack(geneModels, genome = gen, 
                           chromosome = chr, name = "Gene Model", transcriptAnnotation = "symbol",
                           background.title = "brown");
plotTracks(list(itrack, gtrack, atrack, grtrack));
```

<div class="diy">*Try it yourself* `r diyNr<-diyNr+1;diyNr`: change the background color of the CpG track to yellow</div>

<div class="assignment">*Assignment* `r assNr<-assNr+1;assNr`: make a Gviz plot of the mouse genome highlighting your favorite gene, which tracks did you choose to include and why?</div>


### Show aligned reads

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: what types of data tracks can you visualise (hint page 19-20 in the vignette)</div>

#### Coverage (from bam)

Find the list of all BSgenomes - download the right one  
- http://bioconductor.org/packages/release/bioc/html/BSgenome.html
```
Find the mouse package name
install the package for mouse
````

Load the .bam file
```{r loadBam, echo=TRUE}
dTrack4 <- DataTrack(range = "SRR3222411_Aligned.out.bam.sorted.bam", genome = "mm10",
                     type = "h", name = "Coverage", window = NULL, 
                     options(ucscChromosomeNames=FALSE));
```

The ucscChromosomeNames=FALSE is important if the chromosome names in the .bam file does not contain the chr part, to get a complete list of chromosome names use samtools from the *command line*
```
samtools view SRR3222411_Aligned.out.bam.sorted.bam | cut -f 3 | sort | uniq -c
```

##### Coverage for Taz
```{r, echo=TRUE}
plotTracks(dTrack4, from = 74280718, to = 74292151, chromosome="X");
```

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: do the two representations agree?</div>

### Display parameters

Have a look at the settings sections of the Gviz manual (https://bioconductor.org/packages/devel/bioc/manuals/Gviz/man/Gviz.pdf)
```{r, echo=TRUE}
gtrack <- GenomeAxisTrack(col="#FF00FF", cex=1.2);
itrack <- IdeogramTrack(genome = gen, chromosome = chr, fontcolor="red");
plotTracks(list(itrack, gtrack, atrack, grtrack));
```

<div class="diy">*Try it yourself* `r diyNr<-diyNr+1;diyNr`: make the font for the GenomeAxisTrack larger and make the lines magenta, change the color of the chromosome name for the IdeogramTrack</div>
Re-plot the tracks to see that the above changes have taken effect</div>
<div class="diy">*Try it yourself* `r diyNr<-diyNr+1;diyNr`: Change another display option for one of the other visible tracks</div>

To see the available display parameters
```{r, echo=TRUE, eval=FALSE}
availableDisplayPars("GenomeAxisTrack");
availableDisplayPars("IdeogramTrack");
```

More Gviz examples
http://www.sthda.com/english/wiki/print.php?id=43



## IGV vs GViz

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: which representation (IGV or Gviz) did you prefer? and why?</div>

<div class="assignment">*Assignment* `r assNr<-assNr+1;assNr`: include the 6 .bam files (control and KO) as grouped line plots (page 21) to the previous assignment</div>


# Extra Assignments For The Querious

```{r, echo=FALSE}
questNr <- diyNr <- 0 # re-set the counter
```

# MA plots

## Microarray
An MA plot is an application of a Bland-Altman plot for visual representation of two channel DNA microarray gene expression data which has been transformed onto the M (log ratios) and A (mean average) scale. (from wikipedia)

## RNASeq
For weakly expressed genes, we have no chance of seeing differential expression, because the low read counts suffer from such high Poisson noise that any biological effect is drowned in the uncertainties from the sampling at a low rate. (http://www.bioconductor.org/help/workflows/rnaseqGene/)

An MA-plot is a plot of log-fold change (M-values, i.e. the log of the ratio of level counts for each  gene between two samples) against the log-average (A-values, i.e. the average level counts for each gene across the two samples).
The MA-plot is a useful to visualize reproducibility between samples of an experiment. From a MA-plot one can see if normalization is needed.
In MA plot, genes with similar expression levels in two samples will appear around the horizontal line
y = 0 (blue). A lowess t (in red) is plotted underlying a possible trend in the bias related to the mean expression.
(http://www.nathalievilla.org/doc/pdf/tutorial-rnaseq.pdf)

## Make a MA plot
````{r, echo=TRUE}
makeMAPlot <- function(dat, conditionColumns, caseColumns){
  x = rowMeans(as.matrix(dat[,conditionColumns]));
  y = rowMeans(as.matrix(dat[,caseColumns]));
  M = x - y;     ## M-values
  A = (x + y)/2; ## A-values
  df = data.frame(A, M);
  ggplot(df, aes(x = A, y = M)) + 
    geom_point(size = 1.5, alpha = 1/5) +
    geom_hline(color = 'blue3', yintercept=0) + 
    stat_smooth(se = FALSE, method = 'loess', color = 'red3');
}
makeMAPlot(pseudoCount, 1:3, 4:6);
```

<div class="quest">*Try it yourself* `r diyNr<-diyNr+1;diyNr`: what would it look like if we didn't log-transform the data?</div>
<div class="quest">*Question* `r questNr<-questNr+1;questNr`: do you think the data needs to be normalised?</div>

```{r, echo=TRUE}
v <- voom(onlyCounts);
makeMAPlot(v, 1:3, 4:6);
```

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: what did this change? was it for the better?</div>


# Is there a relation between read count and gene length?
```{r, echo=TRUE}
maxCountPerGeneAndLength <- summarise(group_by(tableCounts, Geneid, Length), maxCount=max(KO_1,KO_2, KO_3, WT_1, WT_2, WT_3));
ggplot(maxCountPerGeneAndLength, aes(x=Length, y=maxCount))+
  geom_point();
```

The above is not entirely useful, so log-transform your values
```{r, echo=TRUE}
ggplot(maxCountPerGeneAndLength, aes(x=Length, y=log2(maxCount)))+
  geom_point();
```

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: Is it "ok" to only log transform one axis, or should you do both?</div>

```{r, echo=TRUE}
ggplot(maxCountPerGeneAndLength, aes(x=log2(Length), y=log2(maxCount)))+
  geom_point();
```

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: Did that change your answer?</div>



# Generate a PDF file of a plot (or all plots from a R session)

## From inside RStudio
In the bottom right plot panel
```{r, echo=TRUE, eval=FALSE}
Export
PDF
```

## From 'normal' R
```{r, echo=TRUE, eval=FALSE}
pdf(file="volcanoPlot.pdf");
g <- ggplot(data=diffExp, aes(x=logFC, y=-log10(PValue), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value")
print(g);
dev.off() # not a good idea to run from inside your RMarkdown sesssion!
```
**Note** that if you run the above code (eg the pdf and dev commands) you will no longer see the plot from inside RStudio!

**Note** you can easily make a multipage PDF file containing all plots by calling *pdf()* at the start of your R session and *dev.off()* at the end...

## From RMarkdown
It is possible to "tell" RMarkdown to generate a "picture" of each plot
```{r, echo=TRUE, eval=FALSE}
knitr::opts_chunk$set(fig.path = 'plots/', dev = 'pdf'); # save all plots in the plots folder so I can access them from powerpoint as well..
```
The above will create a folder called plots and then make a .pdf file for each.

Of course these type of commands should be part of the 'setup' section at the top of your RMarkdown document.

<div class="diy">*Try it yourself* `r diyNr<-diyNr+1;diyNr`: Add the command to the top of your document, re-knit it, look at the "progress" at the bottom of your RStudio</div>

<div class="quest">*Question* `r questNr<-questNr+1;questNr`: Look at the list of generated pdf files, did you get as many as you expected? Why/Why not? Was the names useful?</div>

# Name your code chunks
Instead when creating a plot, try naming the code chunk eg to call the code chunk VolcanoPlot simply add it at the top of the current code chunk (after the r add a space followed by the name of the chunk).
```{r VolcanoPlot, echo=TRUE}
ggplot(data=diffExp, aes(x=logFC, y=-log10(PValue), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  theme(legend.position = "none") +
  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value")
```

now this plot will be called VolcanoPlot-1.pdf in the plots folder, and when you knit your document these names are shown at the bottom of your RStudio

