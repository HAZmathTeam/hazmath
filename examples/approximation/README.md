# Testing the rational approximation of functions in HAZMATH and HAZNICS

---

The examples below are for approximating the function: f(x) = x^{-0.5} + 0.25*x^{0.5} on the interval [a,1], a<<1.

- to make the hazmath library from this directory:
    
 ---
 
 make -C ../.. distclean ; make -C ../.. config haznics=yes suitesparse=yes lapack=yes shared=yes ; make -C ../.. install
 
 ---

- to test AAA algorithm:
    - in the makefile in this directory ("examples/approximation")  set SRCFILE = test_aaa
    - type "make clean ; make"
    - to test the AAA algorithm run "./test_aaa" as follows:  
    
  ---
  
    ./test_aaa.ex &lt;&lt;EOF_FRAC &gt;frac_aaa.m  
           -0.50 0.50 1.00 0.25 0.00001 1.00  
     EOF_FRAC  
  
  ---
 
- to test the BRASIL algorithm:
    - in the makefile in this directory ("examples/approximation")  set SRC_SRCFILE = test_brasil
    - type "make clean ; make"
    - to test the BRASIL algorithm run "./test_aaa" as follows:  
    
  ---
  
    ./test_brasil.ex &lt;&lt;EOF_FRAC &gt;frac_brasil.m  
           -0.50 0.50 1.00 0.25 0.00001 1.00  
     EOF_FRAC  
  
  ---

- In both cases use octave or matlab to check the approximation. 

# References

1. AAA algorithm: See: Yuji Nakatsukasa, Olivier Sete, and Lloyd N. Trefethen. The AAA 
algorithm for rational approximation. SIAM J. Sci.Comput., 40(3):A1494–A1522, 2018.  
[https://doi.org/10.1137/16M1106122](https://doi.org/10.1137/16M1106122)

2. BRASIL algorithm: C. Hofreither, "An algorithm for best rational approximation based on   
barycentric rational interpolation." Numerical Algorithms, 2021  
[https://doi.org/10.1007/s11075-020-01042-0](https://doi.org/10.1007/s11075-020-01042-0)

---

END of README for Rational Approximation 

(c) HAZMATH 2009-
