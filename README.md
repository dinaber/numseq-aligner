# numseq-aligner
# my_aligner.m
Aligns multiple number sequences similarly to "multialign" in Matlab, but for numbers (and not NT or AA) 

 Input : matrix with numbers. Each row is treated as seperate number sequence. 
 Output: aligned matrix with numbers. Each coloumn represents a single
         number. zeros stand for space holders.

 Method: 
   The matrix is transffered into adjacency direct matrix. 
   A topological ordering algorithm for directed graphs is
   applied to produce the ordered template including all numbers. 
   The sequences are then mapped to the tamplate.

 Important: 
   ###The zeros are disregarded. 
   ###Any duplication are treated but not triplications!
   ###The algorithem must received acyclic graph - > self loops must be removed.
 
 EXAMPLE:
   Input: [1 2 3 4 5 0;
           3 9 5 7 0 0;
           1 5 8 0 0 0]

   Output:[1 2 3 4 0 5 0 0;
           0 0 3 0 9 5 7 0;
           1 0 0 0 0 5 0 8]
 
 
