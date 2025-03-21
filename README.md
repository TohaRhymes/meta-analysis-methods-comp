## Code for "Biobanking as a tool for genomic research" article

* `plots.py` - Python script for drawing fig.2 from article.

*  All summary statistics can be found [here](https://drive.google.com/drive/folders/1xUMZRlYVHdfEbbSGNFaqDUH7UhHjAyxF?usp=sharing):
    * `se_not_scaled__RHEUMA_OTHER_WIDE___M061.TBL` - meta-analysis using METAL with `SCHEME=STDERR` and without effect scaling;
    * `se_scaled__RHEUMA_OTHER_WIDE___M061.TBL` - meta-analysis using METAL with `SCHEME=STDERR` and with effect scaling;
    * `ss__RHEUMA_OTHER_WIDE___M061.TBL` - meta-analysis using METAL with `SCHEME=SAMPLESIZE` and without effect scaling.

## Citation

When use, please cite:

```

@ARTICLE{Lazareva2022-gh,
  title    = "Biobanking as a Tool for Genomic Research: From Allele
              Frequencies to {Cross-Ancestry} Association Studies",
  author   = "Lazareva, Tatyana E and Barbitoff, Yury A and Changalidis, Anton
              I and Tkachenko, Alexander A and Maksiutenko, Evgeniia M and
              Nasykhova, Yulia A and Glotov, Andrey S",
  abstract = "In recent years, great advances have been made in the field of
              collection, storage, and analysis of biological samples. Large
              collections of samples, biobanks, have been established in many
              countries. Biobanks typically collect large amounts of biological
              samples and associated clinical information; the largest
              collections include over a million samples. In this review, we
              summarize the main directions in which biobanks aid medical
              genetics and genomic research, from providing reference allele
              frequency information to allowing large-scale cross-ancestry
              meta-analyses. The largest biobanks greatly vary in the size of
              the collection, and the amount of available phenotype and
              genotype data. Nevertheless, all of them are extensively used in
              genomics, providing a rich resource for genome-wide association
              analysis, genetic epidemiology, and statistical research into the
              structure, function, and evolution of the human genome. Recently,
              multiple research efforts were based on trans-biobank data
              integration, which increases sample size and allows for the
              identification of robust genetic associations. We provide
              prominent examples of such data integration and discuss important
              caveats which have to be taken into account in trans-biobank
              research.",
  journal  = "J Pers Med",
  volume   =  12,
  number   =  12,
  month    =  dec,
  year     =  2022,
  address  = "Switzerland",
  keywords = "GWAS; allele frequency; biobank; genomics; meta-analysis",
  language = "en"
}


```
