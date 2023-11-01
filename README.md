# AlternativeSplicing

 Data: Data.csv contains structured data in tab delimited columns about genome variables for different species spanning the whole tree of life:
 
   Column 1: GFF file name of whole-genome assembly annotated files downloaded from the NCBI database
   
   Column 2: Name of organism
   
   Column 3,4,5: Taxonomy classification at different levels. Classification done by the NCBI Taxonomy Database
   
   Column 6: Genome Size. Output from 'genomevariables.cpp'
   
   Column 7: Coding Size. Output from 'genomevariables.cpp'
   
   Column 8: Gene Size. Output from 'genomevariables.cpp'
   
   Column 9: Alternative Splicing Ratio. Output from 'genomevariables.cpp'


 Test: 
 
   > g++ genomevariables.cpp -o o
   > ./o GCF_004115265.1_mRhiFer1_v1.p_genomic.gff 


 Output: Creates an Output file with the following information:
 
         > #GCF_004115265.1_mRhiFer1_v1.p_genomic.gff
         > Genome Size: 2075768562
         > Coding Size: 33885542
         > Gene Size: 984758131
         > ALternative Splicing Ratio: 2.9527092114979303
