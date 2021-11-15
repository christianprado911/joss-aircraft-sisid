# Longitudinal and Latero-directional Identification
---
title: `Application of reverse derivative on aircraft system Identification`
tags:
  - matlab
  - reverse
  - identification
  - aircraft

authors:
  - name: Christian Prado S. Machado [co-first author] 
    orcid: 0000-0001-7254-8188
    affiliation: 1
  - name: Dimas Abreu Archanjo Dutra [co-first author]
    orcid: 0000-0001-8202-744X
    affiliation: 2

affiliations:
   - name: Master student, mechanical engineering graduate program, Universidade Federal de Minas Gerais - Brazil
     index: 1
   - name: Teacher and researcher, mechanical engineering department, Universidade Federal de Minas Gerais - Brzil
     index: 2
 
date: 15 november 2021
bibliography: bib.bib
---

#Summary

During the system identification process, the use of derivatives is part of the evaluation and optimization process. The use of reverse derivative could decrease the computational load during derivatives calculation if the number of outputs is greater than the number of inputs, which is our object of interest. A comparision was also made with other types of numerical derivation, such as forward, central and complex step. The output error method was the base of the identification process. 

# Statemend of need

At this work, we developed a program capable of receiving the input and output data of maneuvers of an aircraft, more especifically short period and roll provided by Ravindra at his book "Flight Vehicle System Identification - Second Edition" (Ref. 1 at bib.bib), so the aerodynamic derivatives could be evaluated in a simplified aircraft model. The plataform used was chosen dua to the speed of execution in its processes and by its presence in the engineering market. 

This software was developed in order to develop the student's knowledge about system identification at aircraft maneuvers and enable the assessment of the proposed problem. 

# Mathematics

The forward, central and complex step differentiantion are defined as:

\begin{align}
\displaystyle
    F'(x) &\approx \frac{F(x+h) - F(x)}{h},
    &
    F'(x) &\approx \frac{F(x+h) - F(x-h)}{2h},
    &
    F'(x) &\approx
    \frac{\operatorname{Im}[F(x+\imath h)]}{h};
\end{align}

where $F\colon\reals\to\reals$ its any function. The forward and reverse mode use the chain rule to differentiate the functions. In numerical applications it is necessary to use symbolic functions to implement both types of differentions as  well as done analytically. 

The computational cost can be calculated from any function that uses derivatives, so we can elaborate a table relative to a example function. 
