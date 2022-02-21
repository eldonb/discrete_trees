\pdfoutput=1
\documentclass[a4paper,12pt]{cweb}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[lf]{Baskervaldx}
\usepackage[bigdelims,vvarbb]{newtxmath}
\usepackage{amsfonts, amsmath, amssymb}
\usepackage{fullpage}
\usepackage{marvosym}
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
\usepackage{orcidlink}
\setstretch{1.5}
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\}}}}
\newcommand{\EE}[1]{\ensuremath{\mathds{E}\left[ #1 \right]}}%
\newcommand{\im}{\ensuremath{\imath} }%
\newcommand{\jm}{\ensuremath{\jmath} }%
\newcommand{\be}{\begin{equation}}%
\newcommand{\ee}{\end{equation}}%
\newcommand{\prb}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }%
\newcommand{\h}[1]{\ensuremath{\uptheta_{ #1 } } }%
\newcommand{\VV}[1]{\ensuremath{ \mathbb{V}\left( #1 \right)}}%
\newcommand{\hp}{\ensuremath{\theta_1}}
\newcommand{\hs}{\ensuremath{\theta_2}}
\newcommand{\D}{\ensuremath{\mathbb{D}}}
\newcommand{\F}{\ensuremath{\mathbb{F}} }
\newcommand{\G}{\ensuremath{\mathbb{G}} }
\newcommand{\bt}[1]{\textcolor{blue}{\tt #1}}
\makeatletter
\renewcommand{\maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{flushright}
  \textbf{\LARGE \@@title}

  \@@author
\end{flushright}\egroup
}
\makeatother
\title{Beta-coalescents when sample size is large}
\author{Bjarki Eldon\footnote{MfN Berlin, Germany} \footnote{Supported by 
  Deutsche Forschungsgemeinschaft (DFG) - Projektnummer 273887127 
%% „funded by the Deutsche Forschungsgemein-schaft (DFG, German Research Foundation) –Projektnummer(n)“.
%% 
  through DFG SPP 1819: Rapid Evolutionary Adaptation grant STE 325/17-2 to Wolfgang Stephan; acknowledge  funding by the Icelandic Centre of Research through an
Icelandic Research Fund Grant of Excellence no.\
185151-051 to  Einar \'Arnason, Katr\'in Halld\'orsd\'ottir, Alison M.\ Etheridge,   WS, and BE. BE also acknowledges Start-up module grants through SPP 1819  with Jere Koskela and Maite Wilke-Berenguer, and  with Iulia Dahmer. \\ \today} \orcidlink{https://orcid.org/0000-0001-9354-2391} }

\begin{document}
\maketitle

\rule{\textwidth}{.8pt}


\begin{abstract}
 This code samples discrete trees, i.e.\  traces gene genealogies of a
 sample from a finite population, records  the branch lengths, and
 estimates expected relative branch lengths $\EE{R_{i}^{N}(n)}$ for $1
 \le i \le n-1$.   We
 are interested in how $\EE{R_{i}^{N}(n)}$ behaves as a function of sample size $n$
when the sample is  from a haploid population of finite size $N$ and
the population evolves  according to random
 sweepstakes.   
\end{abstract}

\tableofcontents


@* {\bf Copyright}. 


Copyright {\textcopyright} {\the\year}  Bjarki Eldon \newline

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


@* {\bf Compilation,  output and execution}. 
\label{compile}

 This CWEB
      \citep{knuth1994cweb} document (the {\tt .w} file) can be
      compiled with {\tt cweave} to generate a {\tt .tex} file, and
      with {\tt ctangle} to generate a {\tt .c} \citep{kernighan1988c}
      file.

One can use {\tt cweave} to generate a {\tt .tex} file, and {\tt
ctangle} to generate a {\tt .c} file. To compile the C++ code (the {\tt
.c} file), one needs the GNU Scientific Library. 
Using a Makefile can be helpful, naming this file {\tt iguana.w}


 {\tt
iguana.pdf : iguana.tex \\
\tab\quad\quad\quad\quad cweave iguana.w \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        bibtex iguana \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        ctangle iguana \\
\tab\quad\quad\quad\quad        c++ -Wall -Wextra -pedantic -O3 -march=native -m64 iguana.c -lm -lgsl -lgslcblas \\
        
       
clean :  \\
\tab\quad\quad\quad\quad        rm -vf iguana.c iguana.tex \\
}


Use {\tt valgrind} to check for memory leaks:

{\tt valgrind -v --leak-check=full --show-leak-kinds=all <program call>}




@* {\bf introduction}. 
\label{intro}


We consider a haploid population of fixed  size $N$. Let $X^N, X_1^N, \ldots, X_N^N$ be i.i.d.\ discrete random variables taking values in $\{1, \ldots, \Psi_N\}$; the $X_1^N, \ldots, X_N^N$ denote the random number of juveniles independently   produced  in a given generation according to

\begin{equation}
\label{eq:1}
   \prb{X^N = k} = \frac{ (\Psi_N +1)^\alpha }{ (\Psi_N + 1)^\alpha -
1 } \left( \frac{1}{k^\alpha} - \frac{1}{(k+1)^{\alpha}} \right),
\quad 1 \leq k \leq \Psi_N.
\end{equation}
The mass in Eq \eqref{eq:1} is normalised so that $\prb{ 1 \leq
X^N \leq \Psi_N} =1 $, and $\prb{X^N = k} \ge \prb{X^N = k+1}$. Given
a pool of at least $N$ juveniles, we sample $N$ juveniles for the next
generation.  Leaving out an atom at zero gives $X_1^N + \cdots + X_N^N
\ge N$ almost surely, guaranteeing that we always have at least $N$
juveniles to choose from in each generation.  If $1 < \alpha < 2$ and
$\Psi_{N}/N \not\to 0$ the ancestral process tracing the random
ancestral relations of leaves  converges in finite-dimensional
distributions  to the Beta$(2-\alpha,\alpha)$-coalescent; if $\alpha
\ge 2$ or $\Psi_{N}/N \to 0$ converges to  Kingman. Thus, the model
described in Eq~\eqref{eq:1} is a flexible and realistic and mathematically tractable   model of
individual recruitment success. The Schweinsberg model, i.e.\ with
each individual (independently)  producing at most a finite number of juveniles, 
with number of
juveniles distributed according to Eq~\eqref{eq:1},  or some variant of it,  can and should replace the
Wright-Fisher model. 


Let $B_{i}^{N}(n)$ denote the random total  length of branches
supporting $i \in \{1, 2, \ldots, n\}$ leaves, with the length
measured in generations, and $n$ sample size.     Write $B^{N}(n) :=
B_{1}^{N}(n) + \cdots + B_{n-1}^{N}(n)$ for the random total tree length, and $R_{i}^{N}(n) :=
B_{i}^{N}(n)/B^{N}(n)$ the relative branch lengths.   We are interested in $\EE{R_{i}^{N}(n)}$, and
how   $\EE{R_{i}^{N}(n)}$ behaves as a function of $n$.  This has been
investigated in the case of evolution according to the Wright-Fisher
model. We are interested in this question in the case of evolution
according to random sweepstakes, with  juveniles produced according to
Eq~\eqref{eq:1} and then  sampled without replacement.  

@* {\bf Code}. 
\label{sec:code}

The included libraries are listed in \S~\ref{sec:includes}, and  the
global constants in \S~\ref{sec:constants}.  The random number
generators are defined in \S~\ref{SEC:rng}.  The mass function in
Eq~\eqref{eq:1} is computed in \S~\ref{sec:mass},  the corresponding
CDF for sampling juveniles is computed in \S~\ref{sec:cdf}.  A random
number of juveniles for a single individual is drawn in
\S~\ref{sec:randomjuvs}, and for all the $N$ individuals in
\S~\ref{sec:pool}. From the pool of juveniles (there are at least  $N$ of them
almost surely)   $N$ are sampled without replacement, and this is the
algorithm's 
Acciles's heel.  Sampling without replacement  is very inefficient.
To try to improve the efficiency in \S~\ref{sec:estimatecoalpr}  we  estimate the pairwise  coalescence
probability $c_{N}$ Eq~\eqref{eq:2}.  When there are two blocks
left, the random  time until the last two blocks merge is geometric with
success probability $c_{N}$.   When there are  three blocks left, we
can  use a similar approach,  i.e.\ estimate the corresponding
coalescence probabilities Eq~\eqref{eq:3} and Eq~\eqref{eq:4}, and
sample geometric times    In \S~\ref{sec:rmvhyper} we assign blocks
to  families given a realisation of number of juveniles per
individual.  If there are currently $n$ blocks, the  joint
distribution of number of blocks per family is multivariate
hypergeometric conditional on the number of juveniles.  Given a
realisation $x_{1}\ldots, x_{N}$ of $X_{1}, \ldots, X_{N}$  we approximate
the  number of blocks per family with
\begin{equation}
\label{eq:5}
\prb{\nu = (v_{1}, \ldots, v_{N})} =  \frac{ \binom{x_{1}}{v_{1}} \cdots \binom{x_{N}}{v_{N}}    }{ \binom{x_{1} + \cdots + x_{N}}{n} } 
\end{equation}
 This approximation well approximates  first sampling $N$ juveniles
 and using the surviving juveniles (offspring) in the hypergeometric
 Eq~\eqref{eq:5}.    The hypergeometric step in the algorithm is the
 bottleneck.  Given the number of blocks per family we can update the
 tree \S~\ref{sec:updatetree}. The tree is a vector of block sizes,
 and we merge blocks by first shuffling the order of the blocks, and
 then merging the rightmost blocks.  The new blocks are then
 appended to the  tree.  Thus, we only store the current configuration
 of the tree. Obviously if at most one block is assigned to each
 family then the tree is unchanged over the generation.    Given the
 current tree, the branch lengths for the current realisation  are updated
 in \S~\ref{sec:updateb}, and after each realisation the estimate of
 $\EE{R_{i}^{N}(n)}$ is updated in \S~\ref{sec:estimateri}.  We use
 the coalescence probability estimates \S~\ref{sec:estimatecoalpr} in
 \S~\ref{sec:threetwo}, in the case when there are at most three
 blocks left in the tree.  



@*1 {\bf includes}.
\label{sec:includes}

the included libraries


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
#include <list>
#include <string>
#include <fstream>
#include <forward_list>
#include <chrono>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_gamma.h>


@*1 {\bf constants}.
\label{sec:constants}


the global constants

@<constants@>=@#

/* the $\alpha$ parameter in Eq~\eqref{eq:1} \newline */
const double CONST_ALPHA = 1.05 ;
/* population size $N$ \newline */
const size_t CONST_POP_SIZE = 1e6 ;
/* the cutoff $\Psi_{N}$ in Eq~\eqref{eq:1} \newline */
const double CONST_CUTOFF = 1.0e6 ;
/* sample size \newline */
const size_t CONST_SAMPLE_SIZE = 10000 ;
/* number of experiments \newline */
const double CONST_NUMBER_EXPERIMENTS = 25. ;



@*1 {\bf the random number generator}. 
\label{SEC:rng}


@<gslrng@>=@#
/* the GSL random number engine \newline */
gsl_rng * rngtype ;
 /* obtain a seed out of thin air for the random number engine \newline */
 std::random_device randomseed;
  /* Standard mersenne twister  random number engine seeded with |randomseed()|
  \newline */
  std::mt19937_64 rng(randomseed());

/* set up and initialise the GSL random number generator \newline */
static void setup_rng( unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}



@*1 {\bf the mass function}. 
\label{sec:mass}


The mass function in Eq~\eqref{eq:1}

@<mass@>=@#
static double massfunction( const double j)
{

return ( (pow(1.0/j, CONST_ALPHA) -   pow(1.0/(1. + j),
CONST_ALPHA))/( 1. -  pow( 1./(CONST_CUTOFF + 1.), CONST_ALPHA)) ) ;
}


@*1 {\bf cdf}.
\label{sec:cdf}

 define the function for computing the 
   CDF  for the juvenile distribution in Eq~\eqref{eq:1}.

@<cdf@>=@#
static void mass_function( std::vector<double>& vcdf )
{
  vcdf.clear() ;
  /* we take $\prb{X^{N} = 0} = 0$ \newline */
  vcdf.push_back(0.0) ;
  for( double j = 1; j <= CONST_CUTOFF ; ++j){
  /* adding upp the mass function \S~\ref{sec:mass} \newline */
    vcdf.push_back( vcdf.back() + massfunction(j) );}

}


@*1 {\bf random number of juveniles}.
\label{sec:randomjuvs}

 the function for sampling  a random
   number of juveniles  returning $\min(j \in \mathds{N} : F(j) \ge  u
   )$ for $u$ a random uniform, and $F$ the CDF computed in \S~\ref{sec:cdf}

@<randomjuvs@>=@#
static size_t sample_juveniles( const std::vector<double>& vcdf  )
{
  size_t  j = 1;
  const double u = gsl_rng_uniform( rngtype ); 
  while( vcdf[j] < u){
    ++j ; }
  return(j) ;
}


@*1 {\bf pool of juveniles}.
\label{sec:pool}

sample a random pool of juveniles, i.e.\ each of the $N$ individuals
independently contributes a random number of juveniles according to
Eq~\eqref{eq:1}. 

@<pool@>=@#
static size_t sample_pool_juveniles( std::vector<size_t>& pool_juvs, const std::vector<double>& v_cdf)
{
  pool_juvs.clear();
  size_t s = 0;
  for( size_t i = 0; i < CONST_POP_SIZE; ++i){
  /* record the random number of juveniles for each individual
  \S~\ref{sec:randomjuvs} \newline */
    pool_juvs.push_back( sample_juveniles(v_cdf));
    s += pool_juvs.back() ; }
  /* return  the total number of juveniles $X_{1}^{N} + \cdots + X_{N}^{N}$ \newline */
  return (s);
}



@*1 {\bf estimate coalescence probabilities}. 
\label{sec:estimatecoalpr}


estimate coalescence probabilities for  speeding up reaching the most
recent common ancestor.  When only two blocks left we can sample a
geometric with success probability the pairwise coalescence
probability.  When only three blocks left can sample between a
pairwise merger and  a triple merger. Given a realisation $x_{1},
\ldots, x_{N}$ of $X_{1},
\ldots, X_{N}$ with $s_{N} := x_{1} + \cdots + x_{N}$  the pairwise coalescence
probability is
\begin{equation}
\label{eq:2}
c_{N} =  \sum_{j=1}^{N} \frac{x_{j}(x_{j} - 1)}{s_{N}(s_{N}-1)},
\end{equation}
a 3-merger when three blocks is
\begin{equation}
\label{eq:3}
c_{N}(3;3) = \sum_{j=1}^{N} \frac{(x_{j})_{3}}{(s_{N})_{3}},
\end{equation}
a 2-merger when three blocks is
\begin{equation}
\label{eq:4}
c_{N}(3;2) =  \sum_{j=1}^{N} \frac{3(x_{j})_{2}(s_{N} - x_{j})}{(s_{N})_{3}}.
\end{equation}

@<estimatecoalpr@>=@#
static void estimate_coalescence_probabilities( std::vector<double>& v_cN, const std::vector<double>& v_cdf,   std::vector<size_t>& v_pool_jvs)
{
  size_t SN {} ;
  /* estimate the coalescence probabilites from $10^{3}$ experiments \newline */
  for( size_t i = 0 ; i < 1000 ; ++i){
  /* sample a pool of juveniles and record the total number of
  juveniles \S~\ref{sec:pool} \newline */
    SN = sample_pool_juveniles( v_pool_jvs,  v_cdf) ;
    for( size_t j = 0 ; j < CONST_POP_SIZE; ++j){
      /* the pairwise probability Eq~\eqref{eq:2} \newline */
      v_cN[0] += (static_cast<double>( v_pool_jvs[j])/static_cast<double>(SN)) * (static_cast<double>(v_pool_jvs[j] - 1) / static_cast<double>(SN - 1)) ;
      /* a 3-merger when three blocks Eq~\eqref{eq:3}  \newline  */
      v_cN[1] += (static_cast<double>( v_pool_jvs[j])/static_cast<double>(SN)) *  (static_cast<double>( v_pool_jvs[j]-1)/static_cast<double>(SN-1)) *  (static_cast<double>( v_pool_jvs[j]-2)/static_cast<double>(SN-2)) ;
      /* a merger of two of three blocks Eq~\eqref{eq:4} \newline  */
      v_cN[2] += 3.*(static_cast<double>(SN - v_pool_jvs[j]) / static_cast<double>(SN)) * (static_cast<double>( v_pool_jvs[j])/static_cast<double>(SN - 1)) *(static_cast<double>(v_pool_jvs[j] - 1) / static_cast<double>(SN-2)) ;
    } }
  /* take the average \newline */
  v_cN[0] /= 1000. ;
  v_cN[1] /= 1000. ;
  v_cN[2] /= 1000. ;

  assert( v_cN[0] < 1.);
  assert( v_cN[1] < 1.);
  assert( v_cN[2] < 1.);
  assert( v_cN[2] >  v_cN[1] ); 
}


@*1 {\bf sample a multivariate hypergeometric }.
\label{sec:rmvhyper}

assign the ancestral blocks to families; given a realisation of random
numbers of juveniles the joint  number of blocks per family is a
multivariate hypergeometric.  We sample the marginals and update. 

@<rmvhyper@>=@#
static void rmvhyper( std::vector<size_t>& merger_sizes,  size_t k, const std::vector<size_t>& v_juvs, const size_t SN, gsl_rng *r )
{
  /* |k| is the current number of lines \newline */
  merger_sizes.clear();
  size_t number_of_new_lines = 0 ;
  size_t n_others =  SN - v_juvs[0] ;
  /* sample the number of blocks assigned to the first family  \newline  */
  size_t x = gsl_ran_hypergeometric( r, v_juvs[0], n_others, k);
  if( x > 1){
    /* only record  merger sizes \newline */
    merger_sizes.push_back(x ); }
    /* update the remaining number of blocks \newline */
  k -= x ;
  /* update new number of lines \newline */
  number_of_new_lines += ( x > 0 ? 1 : 0) ;
  size_t i =0 ;
    /* we can stop as soon as all lines 
       have been assigned to a family \newline */
  while( (k > 0) && (i < CONST_POP_SIZE-1) ){
  /* set the index to the one being sampled from \newline */
    ++i ;
    /* update |n_others|  \newline */
    n_others -= v_juvs[i] ;
    x =  gsl_ran_hypergeometric( r,  v_juvs[i], n_others, k );
    if( x > 1){
      merger_sizes.push_back( x) ; }
      /* update the remaining number of blocks \newline */
      k -= x ;
    /* update new number of lines \newline */
    number_of_new_lines += ( x > 0 ? 1 : 0) ;
  }
  /* check if at least two lines assigned to last individual \newline */
  if( k > 1){
    merger_sizes.push_back( k); }
  /* check if at least one line assigned to last individual \newline */
  if( k > 0){
    number_of_new_lines += 1;}
}


@*1 {\bf update the tree}. 
\label{sec:updatetree}

update the tree; the tree is a vector of block sizes. If there are
mergers we shuffle the tree and then   consequtively  merge blocks by summing and recoding  the size of the
merging blocks in each merger, removing the blocks that merge and
eventually adding the new blocks to the tree, the rightmost blocks
merging each time.  

@<updatetree@>=@#
static void update_tree( std::vector<size_t>& tree, const std::vector<size_t>& merger_sizes )
{
  std::vector<size_t> new_blocks {} ;
  if( merger_sizes.size() > 0){
    /* at least one merger \newline */
    new_blocks.clear() ;
    /* shuffle the tree \newline */
    std::ranges::shuffle( tree,  rng );
    /* loop over the mergers \newline */
    for( const auto &m: merger_sizes){
      /* |m| is number of blocks merging;  |m| is at least two; append
      new block to vector of  new blocks \newline */
      assert( m > 1) ;
      /* record the size of the new block by summing the sizes of the
      merging blocks \newline */
      new_blocks.push_back( std::accumulate( std::rbegin( tree), std::rbegin(tree) + m, 0) ) ;
      assert( new_blocks.back() > 1) ;
      /* remove the rightmost |m| merged  blocks from tree \newline */
      tree.resize( tree.size() - m) ;
    }
    /* append new blocks to tree \newline */
    tree.insert( tree.end(),  new_blocks.begin(),  new_blocks.end() ) ;
  }
  /* if no mergers then tree is unchanged */
}


@*1 {\bf update the branch lengths}. 
\label{sec:updateb}

update the branch lengths

@<updateb@>=@#
static void update_ebib( const std::vector< size_t>& tree,  std::vector<double>& vebib)
{
  for( const auto &b: tree){
    /* |b| is size of current block \newline */
    /* update the total tree size and then the branch length
    corresponding to the size of the block \newline */
    vebib[0] += 1.0 ;
    vebib[ b ]  += 1.0 ;}
}




@*1 {\bf update estimate of $\EE{R_{i}^{N}(n)}$}.
\label{sec:estimateri}

update estimate of  $\EE{R_{i}^{N}(n)}$ given branch lengths from one
realisation of a tree 

@<updateri@>=@#
static void update_estimate_ebib( const std::vector<double>& v_tmp,  std::vector<double>& v_ebib)
{
  for ( size_t i = 1 ; i < CONST_SAMPLE_SIZE ; ++i){
    v_ebib[i] += v_tmp[i]/v_tmp[0] ;}
}



@*1 {\bf three or two blocks left}. 
\label{sec:threetwo}

at most three blocks left, so sample times  using the estimates of the
coalescence probabilities \S~\ref{sec:estimatecoalpr}


@<threetwo@>=@#
static void three_or_two_blocks_left(  std::vector<double>& tmp_bib, const std::vector<double>& v_cN,  std::vector<size_t>& v_tree)
{

  double Tk = 0.;
  double Tkk = 0. ;
  size_t newblock {} ;
  switch( v_tree.size() ){
  case 3 : {
    /* three lines left so  sample the two waiting times for a 3-merger and a 2-merger \newline */
    Tk = static_cast<double>( gsl_ran_geometric(rngtype, v_cN[1] ) ) ;
    Tkk = static_cast<double>( gsl_ran_geometric(rngtype, v_cN[2]));
    if( Tk < Tkk){
      /* all three blocks merge;  update the branch lengths  \newline  */

      tmp_bib[0] += (3. * Tk) ;
      tmp_bib[ v_tree[0]] += Tk ;
      tmp_bib[ v_tree[1]] += Tk ;
      tmp_bib[ v_tree[2]] += Tk ;
      /* clear the tree */
      v_tree.clear() ;
      assert( v_tree.size() < 1);
    }
    else{
      /* a 2-merger occurs followed by a merger of the last two
      blocks \newline */
      tmp_bib[0] += (3. * Tkk) ;
      tmp_bib[ v_tree[0]] += Tkk ;
      tmp_bib[ v_tree[1]] += Tkk ;
      tmp_bib[ v_tree[2]] += Tkk ;
      /* shuffle the tree \newline */
      std::ranges::shuffle( v_tree,  rng );
      newblock = v_tree[1] + v_tree[2] ;
      v_tree.resize(1) ;
      v_tree.push_back( newblock);
      assert( v_tree.size() == 2 );
      /* sample waiting time until merger of last two blocks \newline */
      Tk = static_cast<double>( gsl_ran_geometric(rngtype, v_cN[0] ) ) ;
      tmp_bib[0] += (2. * Tk) ;
      tmp_bib[ v_tree[0]] += Tk ;
      tmp_bib[ v_tree[1]] += Tk ;
      v_tree.clear() ;
      assert( v_tree.size() < 1);
    }
    break ; }
  case 2 : {
    /* two blocks left \newline */
  Tk = static_cast<double>( gsl_ran_geometric(rngtype, v_cN[0] ) ) ;
    tmp_bib[0] += (2. * Tk) ;
    tmp_bib[ v_tree[0]] += Tk ;
    tmp_bib[ v_tree[1]] += Tk ;
    v_tree.clear() ;
    assert( v_tree.size() < 1);
    break ; }
  default : break ;
  }
}

@*1 {\bf estimate $\EE{R_{i}^{N}(n)}$ }. 
\label{sec:estimate}

estimate   $\EE{R_{i}^{N}(n)}$ from a given number of experiments. 

@<theestimator@>=@#
static void estimate_ebib( )
{
  std::vector< double > v_cdf ;
  v_cdf.reserve( static_cast<size_t>(CONST_CUTOFF) + 1) ;
  /* compute the CDF function for sampling juveniles \S~\ref{sec:cdf} \newline */
  mass_function( v_cdf);

  std::vector<size_t> v_number_juvs ;
  v_number_juvs.reserve( CONST_POP_SIZE);

/* the tree; initially all blocks are singletons \newline */
  std::vector<size_t> v_tree (CONST_SAMPLE_SIZE, 1) ;
  std::vector<size_t> v_merger_sizes {};
  v_merger_sizes.reserve( CONST_SAMPLE_SIZE );

  std::vector<double> v_tmp_ebib (CONST_SAMPLE_SIZE, 0.0) ;
  std::vector<double> v_ebib (CONST_SAMPLE_SIZE, 0.0) ;
  std::vector<double> v_coal_probs (3, 0.0) ;

  /* estimate the coalescence probs \S~\ref{sec:estimatecoalpr} \newline */
  estimate_coalescence_probabilities( v_coal_probs, v_cdf, v_number_juvs) ;
  
  size_t SN = 0;
  double number_experiments = CONST_NUMBER_EXPERIMENTS + 1.;
  while( --number_experiments  > 0.){
  /* initialise the tree as all singletons \newline */
  v_tree.clear();
  v_tree.assign( CONST_SAMPLE_SIZE, 1);
  /* initialise the container for the branch length for  the current
  realisation \newline */
  std::fill(  std::begin( v_tmp_ebib), std::end(  v_tmp_ebib ), 0.0 );
  assert( std::accumulate(  std::begin( v_tmp_ebib), std::end(
  v_tmp_ebib ), 0.0) == 0); 
  while( v_tree.size() > 1){
    
    /* record the branch lengths for the current tree configuration  \newline */
    update_ebib( v_tree, v_tmp_ebib) ;
    if( v_tree.size() > 3 ){
    /* sample pool of juveniles \newline */
    SN = sample_pool_juveniles( v_number_juvs, v_cdf) ;
    /* compute the merger sizes \S~\ref{sec:rmvhyper} \newline */
    rmvhyper( v_merger_sizes, v_tree.size(),  v_number_juvs, SN, rngtype) ;
    /* update the tree \S~\ref{sec:updatetree} \newline */
    update_tree( v_tree, v_merger_sizes);}
    else{
      /* at most three blocks left \S~\ref{sec:threetwo} \newline */
      three_or_two_blocks_left( v_tmp_ebib, v_coal_probs, v_tree) ;
    }
  }
  /* update estimate of   $\EE{R_{i}^{N}(n)}$ \newline  */
  update_estimate_ebib( v_tmp_ebib, v_ebib);
  }
  /* print the  estimate of $\EE{R_{i}^{N}}$ \newline */
  for( const auto&r: v_ebib){
    std::cout << r << '\n' ; }
}
 


@*1 {\bf the main module}.
\label{sec:main}

The |main| function

@C

@<includes@>@#
@<gslrng@>@#
@<constants@>@#
@<mass@>@#
@<cdf@>@#
@<randomjuvs@>@#
@<pool@>@#
@<estimatecoalpr@>@#
@<rmvhyper@>@#
@<updatetree@>@#
@<updateb@>@#
@<updateri@>@#
@<threetwo@>@#
@<theestimator@>@#

int main(int argc, char *argv[])
{

/* initialise the GSL random number generator |rngtype|
\S~\ref{SEC:rng} \newline */
setup_rng(  static_cast<unsigned long int>(atoi(argv[1])) );

/* estimate $\EE{R_{i}^{N}(n)}$ \S~\ref{sec:estimate} \newline */
estimate_ebib( ) ;

gsl_rng_free( rngtype ); 
return GSL_SUCCESS ; 
}


@* {\bf conclusion}.
\label{sec:concl}


sampling without versus with replacement  results in qualitatively
different behaviour of $\EE{R_{i}^{N}(n)}$ as a function of sample
size.  The Schweinsberg model\cite{schweinsberg03}, where  individuals produce at most
finite number of juveniles, from which one samples without replacement
to form a new set of reproducing individuals, is a natural, realistic,
and mathematically tractable way of  modeling  recruitment dynamics.


@* {\bf bibliography}.
\label{sec:refs}


\bibliographystyle{plain}
\bibliography{refs}

@
\end{document}