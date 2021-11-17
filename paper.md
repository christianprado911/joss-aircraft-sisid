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

# Summary

During the system identification process, the use of derivatives is part of the evaluation and optimization process. The reverse derivative could decrease the computational load during derivatives calculation if the number of outputs is greater than the number of inputs, which is our object of interest. A comparision was also made with other types of numerical derivation, such as forward, central and complex step. The output error method was the base of the identification process. 

# Statemend of need

At this work, we developed a program capable of receiving the input and output data of maneuvers of an aircraft, more especifically short period and roll provided by Ravindra at his book "Flight Vehicle System Identification - Second Edition" `@ravindra:2015`, so the aerodynamic derivatives could be evaluated in a simplified aircraft model. The plataform used was chosen dua to the speed of execution in its processes and by its presence in the engineering market. 

This software was developed in order to develop the student's knowledge about system identification at aircraft maneuvers and enable the assessment of the proposed problem. 

# Mathematics

The forward, central and complex step differentiantion are defined as:

\begin{align}
    F'(x) &\approx \frac{F(x+h) - F(x)}{h},
    &
    F'(x) &\approx \frac{F(x+h) - F(x-h)}{2h},
    &
    F'(x) &\approx
    \frac{\operatorname{Im}[F(x+\imath h)]}{h};
\end{align}

where  $F: \mathbb{R}\to\mathbb{R} $ its any function. These are types of finite-difference methods that are computationally expensive and their accuracy are dependent of the size of de perturbation $h$ in forward and central step, or need to be work in a complex context in the last case. The forrward and central differences leads to truncation and round-off error due to machine  discretization on numbers, while the complex step does not suffer from such problems, but neither does it lead to an exact result. See pictures on \autoref{fig:finitediff}, where these values were calculated in relation to the exact (analytical) differentiation.

![Erro vs size of perturbation h.\label{fig:finite}]![finite](https://user-images.githubusercontent.com/52748683/142127193-3601f748-6a69-4780-8cff-e46eff70b1e2.png)


In a different way, forward (direct) and reverse (also knows as adjoint method) method uses the chain rule to differentiate functions in any case and are analytically accurate. For numerical applications it is necessary to use symbolic functions to implement both types of differentials. The adjoint variable `@adjoint:2008` $\bar{v}_ i$ are introduced here to indicate the value of the sensitivity of that variable relative to all successors and it is multiplied by the adjoint of the successor.

$$
\bar{v}_j = \sum_{i<j} \frac{\partial \phi_i(v_j)_{j<i}}{\partial v_j} \bar{v}_i
$$

where $<$ indicates the successor function that the actual function are dependent. Here the forward method have to be done one time to, after that, we can compute the value of all adjoints needed, dependent on what output and input we are interest. The index value and the value of each intermediate function itselfe have to be save so we can evaluate every adjoint needed, they can not be overwriting. This lead to a particular forward method plus a returning swep that despite having additional steps for the calculation of each derivative, it allows to evaluate only the inputs of interest (such as successor variables), making the calculation of each output restricted only as inputs on which it is dependent.  

The computational cost can be calculated from any function that uses derivatives, for instance, the cost function $J$  of the output error method applied on aircraft system identification `@ravindra:2015` (sec 4.3). 

\begin{equation}
    J(\theta) := \frac12 \sum_{k=1}^N [z(t_k) - y(t_k; \theta)] \trans R\inv [z(t_k) - y(t_k; \theta)],
\end{equation}

 where $x\in\reals^{n_x}$ are the states, $y\in\reals^{n_y}$ are the outputs, $\theta\in\reals^{n_\theta}$ are the unknown parameters, $R\in\reals^{n_y\times n_y}$ is the covariance matrix  and $z$ are the measured value from data. The forward (direct) method leads to 
 
 \begin{equation}
 \nabla J(\theta) = -\sum_{k=1}^N \frac{\dd y(t_k;\theta)}{\dd \theta}\trans R\inv [z(t_k) - y(t_k; \theta)].
 \end{equation}
 
 and the adjoint method 
 
\begin{equation}
 \nabla J(\theta) = \sum_{k=1}^N [\nabla_\theta g(x(t_k;\theta), u(t_k), \theta) \trans \bar y(t_k;\theta) + T\nabla_\theta f(x(t_k;\theta), u(t_k), \theta) \trans \bar x(t_{k+1};\theta)].
\end{equation}

Evaluating the cost of each derivative we can build a table comparing the sum and multiplication processes

\begin{table}[!htb]
\begin{tabular}{c|l|l|}
\cline{2-3}
\multicolumn{1}{l|}{}                & \multicolumn{1}{c|}{\textbf{Forward}} & \multicolumn{1}{c|}{\textbf{Adjoint}} \\ \hline
\multicolumn{1}{|c|}{\textbf{mult.}} & $n_{\theta} n_y (1 + n_y)$            & $n_{\theta}(1 + n_y + 2n_x)$          \\ \hline
\multicolumn{1}{|c|}{\textbf{sum}}   & $n_y(1 + n_{\theta})$                 & $n_{\theta}(1 + n_x + n_y)$           \\ \hline
\end{tabular}
\end{table}

# Acknowledgments

I thank CAPES for the scholarship granted and my professor Dimas for all his support and dedication. 

# References

