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
\usepackage{url}
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

  \@@date
\end{flushright}\egroup
}
\makeatother
\title{Exact expected branch lengths for Beta-coalescents}
\author{Bjarki Eldon\footnote{MfN Berlin, Germany\\ \today }\orcidlink{https://orcid.org/0000-0001-9354-2391}\footnote{Supported by 
  Deutsche Forschungsgemeinschaft (DFG) - Projektnummer 273887127 
  through DFG SPP 1819: Rapid Evolutionary Adaptation grant STE 325/17-2 to Wolfgang Stephan; acknowledge  funding by Icelandic Centre of Research through an
Icelandic Research Fund Grant of Excellence no.\
185151-051 to  Einar \'Arnason, Katr\'in Halld\'orsd\'ottir, Wolfgang Stephan, Alison Etheridge, and BE.} \  }

\begin{document}
\maketitle

\rule{\textwidth}{.8pt}



\begin{abstract}
This code computes exact expected branch lengths  for  Beta-coalescents using recursions \cite{BBE2013a}.  The incomplete Beta-coalescent is derived in \cite{chetwyn-diggle_beta}; if you use the results and/or the code  please cite  \cite{chetwyn-diggle_beta} when it comes out.   
\end{abstract}




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
.c} file), one needs the GNU Scientific Library (GSL).


Use {\tt splint} to check the code:

{\tt splint ebikplusbeta\_cweb.c}


Use {\tt valgrind} to check for memory leaks:

{\tt valgrind -v --leak-check=full --show-leak-kinds=all <program call> }



@* {\bf introduction}. 
\label{intro}


We consider $\Lambda$-coalescents with transition rates, for $1 <\alpha < 2$ and $2 \le k \le n$, and $C> 0$ is a constant of proportionality, 
\begin{equation}
\label{eq:lambdarates}
\lambda_{n,k} = \binom{n}{k}\left(   \one{k=2} +   \frac{C}{B(x;2-\alpha,\alpha)} \int_0^x t^{k-\alpha-1}(1-t)^{n-k + \alpha - 1}dt \right),
\end{equation}
where $B(x;2-\alpha,\alpha) =  \int_0^x t^{1-\alpha}(1-t)^{\alpha - 1}dt$, $0 < x \le 1$.   We restrict to the case $1 < \alpha < 2$ and \cite{chetwyn-diggle_beta}
\begin{equation}
\label{eq:cutoff}
x =   \frac{K }{K + m_\infty}
\end{equation}
where $K > 0$ is a constant and 
\begin{equation}
\label{eq:2}
m_\infty =  1 + \frac{2^{1-\alpha}}{\alpha - 1}.
\end{equation}

Let $B_i(n)$ denote the random length of branches supporting $i \in \{1, \ldots, n-1\}$ leaves. The code computes the exact expected values $\EE{B_i(n)}$ using recursions \cite{BBE2013a}.  



@* {\bf Code}. 
\label{SEC:code}

@*1 {\bf Includes}. 
\label{incl}

The included libraries. 

@<Includes@>=@#
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h> 
#include <gsl/gsl_vector.h> 
#include <gsl/gsl_matrix.h> 
#include <gsl/gsl_sf_pow_int.h> 
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_sf_elementary.h> 
#include <gsl/gsl_sf_gamma.h> 
#include <gsl/gsl_fit.h> 
#include <gsl/gsl_multifit_nlin.h> 
#include <gsl/gsl_sf_exp.h> 
#include <gsl/gsl_sf_log.h> 
#include <gsl/gsl_sort.h> 
#include <gsl/gsl_statistics_double.h> 
#include <gsl/gsl_integration.h> 
#include <gsl/gsl_errno.h> 


@*1 {\bf The beta function}. 
\label{SEC:betafunction}

Compute the (incomplete) beta function
\begin{equation}
\label{eq:1}
B(x; a,b) :=   \int_0^x t^{a-1}(1-t)^{b-1}dt
\end{equation}
for $0< x \le 1$ and $a,b > 0$.  

@<beta function@>=@#
/* return the (incomplete) beta function \newline  */
static double betafunc( const double x, const double a, const  double b )
{
  /* the GSL incomplete beta function is normalised by the complete beta function  */
  return( x < 1. ? gsl_sf_beta_inc( a, b, x) * gsl_sf_beta( a, b) : gsl_sf_beta(a,b)  ); 
}


@*1 {\bf The Beta-coalescent rate}.
\label{SEC:betacoalrate}


Compute the total Beta-coal rate of merging blocks as in Eq~(\ref{eq:lambdarates}).  

@<betarate@>=@#
static double  lbetank( const size_t n, const size_t j,  const  double a, const double K )
{
  /* compute beta rate ; |a| is $\alpha$ \newline */
  /* |n| is current number of blocks \newline */
  /* |j| is number of blocks to merge \newline */
  assert( n > 1);
  assert( j > 1);
  assert( j <= n);

/* compute the cutoff Eq~(\ref{eq:cutoff}) \newline */
  double x = ( K > 0. ?  K / (K +  (1 + pow(2, 1. -a)/(a-1.) ) ) : 1. ) ;
  /* compute the Beta-part of Eq~(\ref{eq:lambdarates}) \newline */
  return(  gsl_sf_choose( (unsigned)n, (unsigned)j ) *  betafunc( x, ((double)j) - a,  ((double)n) + a -  ((double)j) ) / betafunc(x,  2. - a, a) ) ;
}


@*1 {\bf The jump rate}.  
\label{SEC:jumprate}

Return the jump rate of jumping from $i$ to $j$ blocks using the rate in Eq~(\ref{eq:lambdarates}).

@<jump rate@>=@#
static double qij(const size_t i, const size_t j,  const double a, const double K, const double Cconst)
{
  /* compute jump rate of block counting process */
  /* |a| is $\alpha$ \newline  */
  /* |K| is the cutoff constant $K$ \newline */
  /* |Cconst| is the constant of proportionality $C$ \newline */
  assert( i > 1);
  assert( j < i);
  assert( j > 0);

/* take |Cconst = 1| in case of no Kingman part \newline */
  return(  (Cconst*lbetank( i, i-j+1, a, K)) +  (1.* (i-j == 1 ? gsl_sf_choose( (unsigned)i,2) : 0.) ));
}

@*1 {\bf The jump matrices}.
\label{SEC:jumpmatrices}

Compute the matrices of jump rates and probabilities needed for the recursions \cite{BBE2013a}.

@<compute jump matrices@>=@#
static void QP( const size_t n, const double a, const double K, const double Cconst,    gsl_matrix * Q, gsl_matrix *P  )
{

  /* compute matrices  qij and pij \newline */

  size_t i, j;
  double s = 0;
  double x = 0;
  for( i = 2; i <= n ; ++i){
    assert( i <= n);
    s = 0. ;
    for( j = 1; j < i ; ++j){
    /* compute the jump rate |qij| \ref{SEC:jumprate} \newline */
      x = qij( i, j, a, K, Cconst); 
      s = s +  x;
      gsl_matrix_set( Q, i, j, x);
      gsl_matrix_set( P, i, j, x); }
    gsl_matrix_set( Q, i, i, s);
    for( j = 1; j < i ; j++){
      assert( j < i);
      gsl_matrix_set( P, i, j,  gsl_matrix_get( P, i, j)/s ); } }
}


@*1 {\bf Compute the $g$ matrix}. 
\label{SEC:gmatrix}


@<compute g matrix@>=@#
static void gmatrix( const size_t n,  gsl_matrix * G , gsl_matrix * Q, gsl_matrix * P)
{

  size_t  i, k, m ;
  double s = 0.0 ;

  /* initialise the diagonal */
  for ( i = 2 ; i <= n ; ++i){
    gsl_matrix_set( G, i, i, 1./gsl_matrix_get( Q, i,i) );}

  for( i = 3; i <= n ; ++i){
    assert( i <= n );
    for( m = 2; m < i; ++m){
        s = 0. ;
        for( k = m; k < i ; ++k){
          assert( i <= n);
          assert( k <= n);
          assert( m <= n);
          s = s +  gsl_matrix_get( P, i, k) * gsl_matrix_get( G, k, m) ; }
      gsl_matrix_set( G, i, m, s); } }
}


@*1 {\bf Compute the matrix  $p^{(n)}[k,b]$}.
\label{SEC:matrixpnkb}

@<compute pnkb@>=@#
/* compute the matrix  $p^{(n)}[k,b]$ \newline */
static void pnb( const size_t n, gsl_matrix * lkb, gsl_matrix * G, gsl_matrix * P  )
{

  /* j is nprime from Prop A1 in paper \newline */
  /* lnb is the matrix  $p^{(nprime)}[k,b]$; used for each fixed k  \newline */
  gsl_matrix * lnb = gsl_matrix_calloc( n+1, n+1);
  
  size_t  k, b, j, i ;
  double s = 0.0 ;
  gsl_matrix_set( lkb, n, 1, 1.0 ); 
  for ( k = 2 ; k < n ; k++){
    for( i = k ; i <= n ; i++){
    for ( b = 1; b <= i - k + 1 ; b++){
      gsl_matrix_set( lnb, i, b, ( k == i ? ( b == 1 ? 1.0 : 0.0) : 0.0 ) );
      s = 0. ;
      for( j = k ; j < i ; j++){
        gsl_matrix_set( lnb,  i, b, gsl_matrix_get( lnb, i,b) +  (b > i - j ? (((double)(b - i + j)) * gsl_matrix_get( lnb, j, b - i + j) * (gsl_matrix_get( P, i, j) * gsl_matrix_get( G, j, k) / gsl_matrix_get( G, i, k)) / ((double)j)) : 0.0 ) +  (b < j ? ((((double)(j-b))*gsl_matrix_get( lnb, j, b) *  (gsl_matrix_get( P, i, j) * gsl_matrix_get( G, j, k) / gsl_matrix_get( G, i, k)) / ((double)j) )) : 0.0 ) ) ;
      } } }

    for( j = 1 ; j <= (n - k + 1) ; j++){
      gsl_matrix_set( lkb, k, j, gsl_matrix_get( lnb, n, j) ); }
    gsl_matrix_set_zero( lnb ); }

/* free the |lnb| matrix \newline */
  gsl_matrix_free( lnb );

}

@*1 {\bf compute the expected spectrum}. 
\label{SEC:ebi}

Compute the exact expected spectrum 


@<compute ebi@>=@#
static void compute_ebi( const size_t n, const  double a, const  double K, const double Cconst )
{
/* |n| is sample size; number of leaves \newline */
/* |a| is $\alpha$ \newline */
  /* |K| is the cutoff constant $K$ \newline */
  /* |Cconst| is the constant of proportionality $C$ \newline */

  gsl_matrix * P = gsl_matrix_calloc( n+1, n+1);
  gsl_matrix * Q = gsl_matrix_calloc( n+1, n+1);
  gsl_matrix * G = gsl_matrix_calloc( n+1, n+1);
  gsl_matrix * Pn = gsl_matrix_calloc( n+1, n+1);

  /* compute the Q and P matrices \S~\ref{SEC:jumpmatrices}  \newline */
  QP( n, a, K, Cconst,  Q, P);
  /* compute the $g$ matrix \S~\ref{SEC:gmatrix} \newline */
  gmatrix( n,  G, Q, P);
  /* compute $pnb$ matrix \S~\ref{SEC:matrixpnkb} \newline */
  pnb( n, Pn, G, P);
  
  size_t b, k;
  double s = 0.0 ;
  double eb = 0.0 ;
  double * v_ebi  = (double *)calloc( n, sizeof(double) ); 

  for ( b = 1; b < n ; b++){
    s = 0. ;
    for( k = 2; k <= n-b+1 ; k++){
      s = s +  (gsl_matrix_get( Pn, k, b) * ((double)k) * gsl_matrix_get( G, n, k)); }
    v_ebi[b] = s ; 
    eb = eb + s ;}

/* print out $\EE{B_i(n)}/\EE{B(n)}$ \newline */
  for ( b = 1; b < n ; ++b){
    printf("%g\n", v_ebi[b]/eb); }

  /* free memory \newline */
  gsl_matrix_free( P);
  gsl_matrix_free( Q);
  gsl_matrix_free( G);
  gsl_matrix_free( Pn);
  free( v_ebi);
}





@*1 {\bf the main module}. 
\label{main}

The |main| function for calling |compute_ebi| \S~\ref{SEC:ebi} .

@C

@<Includes@>@#
/* see \S\ref{incl} \newline */
@<beta function@>@#
/* \S\ref{SEC:betafunction} \newline */
@<betarate@>@#
/* \S\ref{SEC:betacoalrate} \newline */
@<jump rate@>@#
/* \S\ref{SEC:jumprate} \newline */
@<compute jump matrices@>@#
/* \S\ref{SEC:jumpmatrices} \newline */
@<compute g matrix@>@#
/* \S\ref{SEC:gmatrix} \newline */
@<compute pnkb@>@#
/* \S\ref{SEC:matrixpnkb} \newline */
@<compute ebi@>@#
/* \S\ref{SEC:ebi} \newline */

 int main( int argc, char *argv[])
{

/* \S~\ref{SEC:ebi} \newline */
  compute_ebi( (size_t)atoi(argv[1]),  atof(argv[2]),  atof(argv[3]),  atof(argv[4])  ) ;

return GSL_SUCCESS ;
}


@* {\bf references}. 


\bibliographystyle{plain}
\bibliography{refs}



@
\end{document}
