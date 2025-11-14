# Paper repository

This repository contains the official code accompanying the paper:

**Knockoffs for low dimensions: changing the nominal level post-hoc to gain power while controlling the FDR**

## Overview

Knockoffs provide a powerful framework for controlled variable selection with false discovery rate (FDR) guarantees. While widely used in high-dimensional settings, their performance can deteriorate in low-dimensional scenarios with sparse signals. A key limitation is that the knockoff filter requires a minimum number of selections determined by the nominal FDR level.

In this work, we introduce a **post-hoc adjustment using e-values** that allows the nominal FDR level to be adapted *after* running the knockoff procedure. This leads to two benefits:

-   **If the standard knockoff filter makes no discoveries:**\
    We can safely increase the nominal level to potentially recover meaningful signals.

-   **If the standard knockoff filter makes discoveries:**\
    We can often reduce the nominal level to improve precision.

These improvements come at no cost: the post-hoc results are always at least as informative as those from the original knockoff procedure. We also extend the method to recently proposed **derandomized knockoff** techniques.

## Repository Structure

-   `util.R` — Implementation of all methods used in the paper, including:
    -   original knockoffs\
    -   derandomized knockoffs\
    -   post-hoc knockoffs\
    -   post-hoc derandomized knockoffs
-   `benchtm_example.R` — Code for reproducing the case study from **Section 6** of the paper.
