sinfa
=====

R code for Sequential INFormation Analysis

Sequential Information Analysis (SINFA) is a way to track which features in classes are transferred best accross 
a communication channel. 

The approach is used in the analysis of perception of diagnostic words by humans.  This code is modeled after the analysis 
published in JASA vol 54(5) by Wang and Bilger (1973).  One of the tables in the article is reproduced, so that you can
follow what is going on.  The features in the case are phonetic features of the phonemes of the words, the classes are 
phonemes. 

Start an R session
$ R
## read the code (once per session)
source("sinfa.R")
## read the data (this assigns global variables x (the confusion matric) and f (the phonetic features)
read.wang.bilger()
x
f
## perform the analysis
s <- sinfa(x, f, 7)
## show the analysis (equivalent to print(s))
s
## print the sequence of highest information carriers of the analysis
summary(s)
