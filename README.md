sinfa
=====

R code for Sequential INFormation Analysis

Sequential Information Analysis (SINFA) is a way to track which features in classes are transferred best accross 
a communication channel. 

The approach is used in the analysis of perception of diagnostic words by humans.  This code is modeled after the analysis published by Wang and Bilger in [JASA **54** _5_ (1973)](http://dx.doi.org/10.1121/1.1914417 "$30").  One of the tables in the article is reproduced, so that you can follow what is going on.  The features in the case are phonetic features of the phonemes of the words, the classes are phonemes. 

The analysis is in terms of information.  In the paper the terms "information" and "transferred information" are used, in the analysis tools we use the terms "entropy" and "mutual entropy". 

## Example ##

Start an R session (`$` is the shell prompt)

    $ R 

Then execute the sommands (`##` is a comment, the R prompt is not shown)

    ## read the code (once per session)
    source("sinfa.R")
    ## read the data (this assigns global variables x (the confusion matrix) and f (the phonetic features)
    read.wang.bilger()
    x 
    f 

This reads the phonetic features _wang-bilger-features.csv_ into global variable `f` and confusion matrix _wang-bilger-table6.csv_ into `x`, and prints them to the terminal 

    ## perform the analysis
    s <- sinfa(x, f, 7)

This performs the analysis (this may take a while) 7 levels deep, i.e., first looking for the feature with most relative information transfer, then taking that feature known and looking for the next feature with most relative information transfer, etc.  The results are stored in a structure of class "sinfa". 

    ## show the analysis (equivalent to print(s))
    s
    ## print the sequence of highest information carriers of the analysis
    summary(s)

There are two ways of displaying the information in a sinfa structure, `print(s)` and `summary(s)`.  `print(s)` prints the entire sinfa analysis, with full tables of entropy, mutual entropy and relative mutual entropy for all features, given the information of a number of features found in previous iterations.  

## Data ##

I entered the data from Wang and Bilger myself, you may find more of the tables in the article [here](http://people.cs.uchicago.edu/~dinoj/research/wangbilger.html "tables"). 

