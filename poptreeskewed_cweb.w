\pdfoutput=1
\documentclass[a4paper,12pt]{cweb}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[lf]{Baskervaldx}
\usepackage[bigdelims,vvarbb]{newtxmath}
\usepackage{amsfonts, amsmath, amssymb}
\usepackage[super,round]{natbib}
\usepackage{fullpage}
\usepackage{marvosym}
\usepackage{graphicx}
\usepackage{bm}
\usepackage[round,numbers,super]{natbib}
\usepackage{color}
\usepackage{a4wide,fullpage}
\usepackage{setspace}
\usepackage{hyperref}
\usepackage{enumerate}
\usepackage{dsfont}
\usepackage[right]{lineno}
\usepackage{verbatim}
\usepackage{tabto}
\usepackage{lipsum}
\usepackage{refcheck}
\usepackage{environ}
\usepackage{orcidlink}
\setstretch{1.5}
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\}}}}
\newcommand{\EE}[1]{\ensuremath{\mathds{E}\left[ #1 \right]}}
\newcommand{\im}{\ensuremath{\imath} }
\newcommand{\jm}{\ensuremath{\jmath} }
\newcommand{\be}{\begin{equation}}
\newcommand{\ee}{\end{equation}}
\newcommand{\bd}{\begin{displaymath}}
\newcommand{\ed}{\end{displaymath}}
\newcommand{\prb}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }
\newcommand{\h}[1]{\ensuremath{\uptheta_{ #1 } } }
\newcommand{\VV}[1]{\ensuremath{ \mathbb{V}\left( #1 \right)}}
\newcommand{\hp}{\ensuremath{\theta_1}}
\newcommand{\hs}{\ensuremath{\theta_2}}
\newcommand{\D}{\ensuremath{\mathbb{D}}}
\newcommand{\F}{\ensuremath{\mathbb{F}}}
\newcommand{\G}{\ensuremath{\mathds{G}}}
\newcommand{\NN}{\ensuremath{\mathds{N}}}
\newcommand{\bt}[1]{\textcolor{blue}{\tt #1}}
\newcommand{\ao}{\ensuremath{{\alpha_1}}}
\newcommand{\at}{\ensuremath{{\alpha_2}}}
\newcommand{\OO}[1]{\ensuremath{\mathcal{O}\left( #1\right)}}
\makeatletter
\renewcommand{\maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{flushright}
  \textbf{\LARGE \@@title}

  \@@author
\end{flushright}\egroup
}
\makeatother
\title{ skewed : tree size under sweepstakes reproduction }
\author{  Bjarki Eldon\footnote{MfN Berlin, Germany\\ \today }\footnote{Supported by 
  Deutsche Forschungsgemeinschaft (DFG) - Projektnummer 273887127 
  through DFG SPP 1819: Rapid Evolutionary Adaptation grant STE 325/17-2 to Wolfgang Stephan; acknowledge  funding by Icelandic Centre of Research through an
Icelandic Research Fund Grant of Excellence no.\
185151-051 to  Einar \'Arnason, Katr\'in Halld\'orsd\'ottir, Wolfgang Stephan, AME, and BE.}\orcidlink{https://orcid.org/0000-0001-9354-2391} }

\begin{document}
\maketitle

\rule{\textwidth}{.8pt}


\begin{abstract}
 This code estimates the total tree size of a sample of a given size
 up to and including the total population.  The population is finite
 and haploid and evolves according to a model of sweepstakes
 reproduction.  The code is part of joint work with Jonathan A.\
 Chetwyn-Diggle and Alison Etheridge about trying to understand what
 happens to genealogies when sample size becomes `large'\cite{chetwyn-diggle_beta}. 
\end{abstract}


\tableofcontents


@* {\bf Copyright}. 

Copyright {\textcopyright} {\the\year}  Bjarki Eldon \newline

\smallskip

skewed : estimates tree size for a sample from a haploid population with sweepstakes reproduction

\smallskip

This document and any source code it contains  is distributed under the terms of the GNU General Public Licence (version $\ge 3$).  You
should have received a copy of the licence along with this file (see file COPYING).  


    The source codes  described in this document  are  free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This document and the code it contains   is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file (see COPYING).  If not, see \url{http://www.gnu.org/licenses/}.


@* {\bf compilation and execution}. 
\label{SEC:comp}

Use {\tt ctangle} to generate a C++ file ending in {\tt .c}, and {\tt
cweave} to generate a {\LaTeX} file ending in {\tt .tex} for compilin
the documentation.   Compiling  the C++ code:

{\tt c++ -Wall -Wextra -pedantic -std=c++20 -m64 -Og -march=native poptreeskewed.c -lm -lgsl -lgslcblas}

The code passes  {\tt valgrind} checks for memory leaks, and a basic {\tt cppcheck}.   


Call the compiled  code with

{\tt ./a.out  \$(shuf -i 1000-100000 -n1)}
@* {\bf intro}. 
\label{SEC:intro}


Consider a haploid population evolving according to sweepstakes
reproduction.  In any given generation the current $N$ individuals
independently produce juveniles (potential offspring) according to a
given law.  From the total pool of at least $N$ juveniles we sample
$N$ juveniles independently and uniformly at random without
replacement to form the next set of reproducing individuals.  The
random number of juveniles is given the distribution
\begin{equation}
\label{eq:1}
\prb{X_1 = k} = \one{1 \le k \le u(N)} \left( \frac{1}{k^\alpha} -  \frac{1}{(1+k)^\alpha} \right) \frac{1}{1 -  \tfrac{1}{(1 + u(N))^\alpha}},
\end{equation}
where $u(N) \in \{1,2,\ldots\}$ is an upper bound on the number of
juveniles.  Assuming an upper bound on the number of juveniles is
biologically reasonable.  We assume $1 < \alpha < 2$.



We are interested in the tree size $B^N(N)$ of the whole population
measured in generations.  We would like to see if we can identify a
scaling $f(N)$ so that $B^N(N)/f(N)$ converges (ideally almost surely,
but in $L_1$ would do), and we are hoping simulating this could give
us an idea of the scaling (although for sure can be misleading).


@* {\bf code}. 
\label{SEC:code}

The sampling is individual-based, meaning we consider a finite haploid
population, traverse the tree backwards-in-time generation by
generation, and measure the tree size in generations.  The code is
simple enough, unfortunately however, the sampling becomes damn slow
for $N\ge 10^6$.  The following sections describe each small part of
the code, each is self-explanatory.  


@*1 {\bf GSL random number generator }. 
\label{SEC:gslrng}

We will use the GSL random number generators for drawing a random
uniform from the unit interval, and in the hypergeometric sampling.

@<GSL rng@>=@#
gsl_rng * rngtype; 
static void setup_rng( unsigned long int s )
{
  const gsl_rng_type *T ; 
  gsl_rng_env_setup(); 
  T = gsl_rng_default ;
  rngtype = gsl_rng_alloc(T);
  gsl_rng_set( rngtype,  s) ;
}


@*1 {\bf constants}. 
\label{SEC:const}

For convenience we define here the key data, the $\alpha$ parameter
|CONST_ALPHA| and the cutoff $u(N)$ |CONST_CUTOFF| of the juvenile
distribution Eq~(\ref{eq:1}), and the population size $N$
|CONST_POP_SIZE|.  The code is configured to estimate the tree size of
the whole population, however  any sample size works.  

@<constants@>=@#
const double CONST_ALPHA = 1.05 ;
const size_t CONST_POP_SIZE = 1e1 ;
const double CONST_CUTOFF = 1.e1 ;

@*1 {\bf includes}. 
\label{SEC:includes}

Probably don't need all the libraries.  

@<includes@>=@#
#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <memory>
#include <utility>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <string>
#include <fstream>
#include <chrono>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


@*1 {\bf compute the cumulative density function }.
\label{SEC:cdf}

Compute the cumulative density function for sampling a number of juveniles using the mass function Eq~(\ref{eq:1}). 

@<cdf@>=@#
static void mass_function( std::vector<double>& vcdf )
{
/* set probability of zero juveniles to zero |vcdf[0] = 0.0|  \newline */
  vcdf.push_back(0.0) ;
  for( double j = 1; j <= CONST_CUTOFF ; ++j){
    /* set |vcdf[j]| to |vcdf[j-1]| $ + \prb{X_1 = j}$ Eq~(\ref{eq:1}) \newline */
    vcdf.push_back( vcdf.back() + (pow(1.0/j, CONST_ALPHA) -   pow(1.0/(1. + j), CONST_ALPHA))/( 1. -  pow( 1./(CONST_CUTOFF + 1.), CONST_ALPHA))  );}


}


@*1 {\bf sample one number of juveniles}. 
\label{SEC:samplejuvs}
Sample one realisation of number of juveniles by the inverse-CDF method using the CDF computed by |mass_function| \S~\ref{SEC:cdf}

@<sample number juveniles@>=@#
static size_t sample_juveniles( const std::vector<double>& vcdf, gsl_rng* r )
{
  size_t  j = 1;
  const double u = gsl_rng_uniform(r) ;
  while( vcdf[j] < u){
    ++j ; }
  return(j) ;
}


@*1 {\bf sample a pool of juveniles}. 
\label{SEC:pool}

Sample a pool of juveniles using |sample_juveniles| \S~\ref{SEC:samplejuvs}

@<pool@>=@#
static size_t sample_pool_juveniles( std::vector<size_t>& pool_juvs, const std::vector<double>& v_cdf, gsl_rng *r)
{
  pool_juvs.clear();
  size_t s = 0;
  for( size_t i = 0; i < CONST_POP_SIZE; ++i){
    pool_juvs.push_back( sample_juveniles(v_cdf, r) );
    s += pool_juvs.back() ; }
  /* return  the total number of juveniles; since $\prb{X_1 = 0} = 0$ we always have at least $N$ juveniles \newline */
 
  return (s);
}


@*1 {\bf compute new number of lines }. 
\label{SEC:newlines}

Compute new number of lines by assigning current lines to families
uniformly at random and without replacement; i.e.\ the joint
distribution of number of lines per family is multivariate
hypergeometric. We sample the marginals, and only need to record the
number of families with at least one line.
Conditional on $(X_1, \ldots, X_N) = (x_1, \ldots, x_N)$ and $m$ current number of lines
\begin{equation}
\label{eq:2}
\prb{(\nu_1,\ldots, \nu_N) = (m_1, \ldots, m_N) } =  \frac{ \binom{x_1}{m_1}\cdots  \binom{x_N}{m_N} }{ \binom{\sum_i x_i}{m} }.
\end{equation}
where $m_1 + \cdots + m_N = m$.  Then the new number of lines is
\begin{equation}
\label{eq:3}
m^\prime = \sum_i \one{m_i > 0}.
\end{equation}
We sample the marginals, meaning we sample consecutively $\nu_1$,
$\nu_2$, etc.\  updating the number of remaining lines  until all the current lines have been assigned a family,
or we have sampled $\nu_{N-1}$, in which case the remaining lines are assigned to $\nu_N$.  


@<new number lines@>=@#
static size_t r_new_number_lines( size_t k, const std::vector<size_t>& v_juvs, const size_t SN, gsl_rng *r )
{
  
  /* |k| is the current number of lines */
  size_t new_number_lines = 0 ;
  size_t n_others =  SN - v_juvs[0] ;
  /* sample $\nu_1$ from the marginal of Eq~(\ref{eq:2}) \newline */ 
  size_t x =  gsl_ran_hypergeometric( r,  v_juvs[0], n_others, k );
  /* update the new number of lines according  to Eq~(\ref{eq:3}) \newline */ 
  new_number_lines += (x > 0 ? 1 : 0);
  k -= x ;
  size_t i =0 ;
  while( (k > 0) && (i < CONST_POP_SIZE-1) ){
    ++i ;
    /* we can stop as soon as all lines 
       have been assigned to a family \newline */
    n_others -= v_juvs[i] ;
    x=  gsl_ran_hypergeometric( r,  v_juvs[i], n_others, k );
    new_number_lines += (x > 0 ? 1 : 0);
    k -= x ;
  }
 new_number_lines += (k > 0 ? 1 : 0);
 /* return $m^\prime$ in Eq~(\ref{eq:3}) \newline */
 return (new_number_lines) ;
}

@*1 {\bf estimate the tree size}.
\label{SEC:estimatetreesize}

Estimate the tree size  


@<estimatesize@>=@#
static void tree_size( gsl_rng *r)
{
 
  /* initialise the container for holding the CDF from Eq~(\ref{eq:1}) \newline */
  std::vector< double > v_cdf ;
  v_cdf.reserve( static_cast<size_t>(CONST_CUTOFF) + 1) ;

/* compute  the CDF from Eq~(\ref{eq:1}) \S~\ref{SEC:cdf} \newline */ 
  mass_function( v_cdf); @#


/* initialise the container for holding realisations of $X_1, \ldots, X_N$ the random number of juveniles \newline */
  std::vector<size_t> v_pjuvs ;
  v_pjuvs.reserve( CONST_POP_SIZE);

/* initialise the local variables \newline */
  size_t current_number_lines = CONST_POP_SIZE ;
  size_t SN = 0;
  size_t tree_size = 0; 
  while( current_number_lines > 1){
    tree_size += current_number_lines ;
    /* sample pool of juveniles \S~\ref{SEC:pool} \newline */
    SN = sample_pool_juveniles( v_pjuvs, v_cdf, r) ;
    /* compute the new number of lines \S~\ref{SEC:newlines} \newline */
    current_number_lines = r_new_number_lines( current_number_lines, v_pjuvs, SN, rngtype ) ;
  }
  printf("%lu\n", tree_size);
  
}




@*1 {\bf the |main| function}
\label{SEC:main}

@C

/* \S~\ref{SEC:includes} \newline */
@<includes@>@#
/* \S~\ref{SEC:gslrng} \newline */
@<GSL rng@>@#
@<constants@>@#
/* \S~\ref{SEC:cdf} \newline */
@<cdf@>@#
/* \S~\ref{SEC:samplejuvs} \newline */
@<sample number juveniles@>@#
/* \S~\ref{SEC:pool} \newline */
@<pool@>@#
/* \S~\ref{SEC:newlines} \newline */
@<new number lines@>@#
/* \S~\ref{SEC:estimatetreesize} \newline */
@<estimatesize@>@#
  int main(int argc, char *argv[])
{
/* estimate the tree size \S~\ref{SEC:estimatetreesize} \newline */

     /* initialise the GSL random number generator \S~\ref{SEC:gslrng} \newline */
  setup_rng( static_cast<unsigned long>( atol(argv[1]) ) );

         tree_size( rngtype ) ;
  gsl_rng_free(rngtype); 
  return (0);
  /* at end of |main| function \newline */
  /* See \S~\ref{SEC:comp} for compiling and running. \newline */
}


@* {\bf references}. 
\label{SEC:refs}

\bibliographystyle{plain}
\bibliography{newrefs}

@
\end{document}

