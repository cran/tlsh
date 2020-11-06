## ---- echo=TRUE, message=FALSE, knitr::opts_chunk$set(cache=TRUE)-------------
library(blink)
library(plyr)
library(tlsh)
data(RLdata500)
head(RLdata500)
data.500 <- RLdata500[-c(2,4)]
head(data.500)

## -----------------------------------------------------------------------------
 blocks <- block_setup_v2(RLdata500, b=22, k=2)
 summary(blocks)

## -----------------------------------------------------------------------------
eval.blocksetup(RLdata500, b=26, key=identity.RLdata500)

## -----------------------------------------------------------------------------
(rr <- reduction.ratio.from.blocking(blocks)) 

