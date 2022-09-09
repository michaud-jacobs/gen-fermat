// Magma code to support the calculations in the paper: On some Generalized Fermat Equations of the form x^2+y^2l=z^p.

// This code carries out the computations needed for the proof of Theorem 1.

p:=7;
L<zet>:=CyclotomicField(p);
K:=sub<L | zet+1/zet>;
OK:=Integers(K);

q:=3*OK; // The unique prime of K above 3.
nq:=Norm(q);

H:=HilbertCuspForms(K,2^3*OK); 
Hnew:=NewSubspace(H);
Dimension(Hnew);   // 5
NewformDecomposition(Hnew); // 5 forms

