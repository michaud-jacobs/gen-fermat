// Magma code to support the calculations in the paper: On some Generalized Fermat Equations of the form x^2+y^2l=z^p.

// This code computes the quantity B defined in Theorem 1 of [15] for Proposition 5.2 of the paper and verifies the dimensions of the Hilbert cusp forms at the necessary levels.

// Computing B
// The code works for any prime p.

p:=11; 
L<zet>:=CyclotomicField(p);
K<t>:=sub<L | zet+1/zet>; 
OK:=Integers(K); 

Aut,_,nu:=AutomorphismGroup(K);

Gal:=[nu(Aut.1)^i : i in [0..((p-1)/2)-1]]; // elements of Gal(K/Q)

UK,psi:=UnitGroup(K);
r:=#Generators(UK);
es:=[K ! (psi(UK.i)) : i in [1..r]];    // basis for unit group

S:=[];                      // we compute the possible non-constant isogeny characters
for i in [1..2^r-2] do
    seq1:=Intseq(i,2);
    seq2:=seq1 cat [0*j : j in [1..r-#seq1]];
    seq:=[12*ss : ss in seq2];
    S:=S cat [seq];
end for;

Ns:=[];          // we compute the twisted norms
for j in [1..r] do
    Nsjs:=[];
    for s in S do
        Nsjs := Nsjs cat [  &* [ (Gal[i](es[j]))^s[i] : i in [1..r]]  ];
    end for;
    Ns:=Ns cat [Nsjs];
end for;


As_s:=[];                // we compute the A_s values
for j in [1..2^r-2] do
    I:=&+[(Ns[i][j]-1)*OK : i in [1..r]];
    if I eq 0*OK then
          As:=0;
          As_s:=As_s cat [As];
    else  As:=Integers() ! Norm(I);
          As_s:=As_s cat [As];
    end if;
end for;

B:=LCM(As_s);        // the required value B.
Factorisation(B);

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

// Dimensions of space of Hilbert cusp forms new at level N_l

p:=11;
L<zet>:=CyclotomicField(p);
K:=sub<L | zet+1/zet>; 
OK:=Integers(K);
H:=HilbertCuspForms(K,2^3*OK);
Hnew:=NewSubspace(H);
Dimension(Hnew); // a few minutes for p=13 and p=17.
