# [MEEG007] Computational Fluid Dynamics Programming Problems(ssodam distribution ver.)
## Warning
본 문서는 기계공학과 학석사 공용과목 MEEG007 전산유체역학 과제 중 프로그래밍/수치해석 과제를 돕기 위하여 업로드되었습니다. 프로그램 코드를 참고하여도 좋으나 내용을 충분히 숙지하는 것을 권고드리고 과제물 표절 등의 이유로 인한 불이익에 대해서 본 필자는 책임을 지지 않겠습니다.

## About 2021-2 Coursework
Author|Contact
---|---
Heesung Kim (Senior in Department of Mechanical Engineering, Sogang University)|mech9917@gmail.com

Documentation Date : 22. June, 2022

## Programming Language/version/envs
- Python 100%
- Python ver 3.9
- Anaconda3-Pycharm

## Imports
- NumPy (for linear algebraic problems)
- Pandas (for data and tables)
- matplotlib (for plotting)
- math (for exponentials)

## Contents
### Chapter 4 The Finite Volume Method for Diffusion Problems
  - 1D Steady State Diffusion
    - Example 4.1 Source-free heat conduction in an insulated rod
    - Example 4.2 1-D steady conduction with heat source
    - Example 4.3 Cooling of a circular fin by means of convective heat transfer along its length
### Chapter 5 The Finite Volume Method for Convection-Diffusion Problems
  - The Central Differencing Scheme
    - Example 5.1 1-D convection-diffusion problem with T(0)=1, T(L)=0
      - Case 1 0.1m/s, 5 nodes
      - Case 2 2.5m/s, 5 nodes
      - Case 3 2.5m/s, 20 nodes
  - The Upwind Differencing Scheme
    - Example 5.2 Upwind scheme for 0.1m/s, 2.5m/s with 5-point grid
  - The Hybrid Differencing Scheme
    - Example 5.3 Hybrid scheme for 2.5m/s, 5-point and 25-point grid
  - Higher order Differencing Schemes
    - Example 5.4 QUICK scheme for 0.2m/s, 5-point grid
### Chapter 7 Solution of Discretized Equations
  - TDMA(Tri-Diagonal Matrix Algorithm)
    - Example 7.1 TDMA for example 4.3
  - Application of TDMA to 2-D problems
    - Example 7.2 Two-dimensional diffusion problem
