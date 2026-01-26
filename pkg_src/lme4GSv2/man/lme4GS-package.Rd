\name{lme4GS-package}
\alias{lme4GS}
\alias{lme4GS-package}
\title{lme4 for Genomic Selection}
\description{
Flexible functions that use 'lme4' as computational engine for
fitting models used in Genomic Selection (GS). GS is a technology used for genetic
improvement, and it has many advantages over phenotype-based selection. There
are several statistical models that adequately approach the statistical
challenges in GS, such as in linear mixed models (LMMs). The 'lme4' is the
standard package for fitting linear and generalized LMMs in the R-package,
but its use for genetic analysis is limited because it does not allow the
correlation between individuals or groups of individuals to be defined. The
'lme4GS' package is focused on fitting LMMs with covariance structures defined
by the user, bandwidth selection, and genomic prediction. The new package is
focused on genomic prediction of the models used in GS and can fit LMMs using
different variance-covariance matrices. Several examples of GS models are
presented using this package as well as the analysis using real data. For more
details see Caamal-Pat et.al. (2021) <doi:10.3389/fgene.2021.680569>.
}
