/*********************************************************************/  
//  Estimates genome variables from GFF genome annotation files.
//  Genome variables: 
//  -> GenomeSize
//  -> Coding
//  -> GeneSize
//  -> AlternativeSplicingRatio
/*********************************************************************/
//
// Input: GFF or GFF3 file
// Output: Genome variables
//
// g++ genomevariables.cpp -o o
// -o gff_file_name.gff
//
//It assumes:
//All rows not starting with # correspond to a data type
//All rows have 9 columns, with a fixed structure type data
/*********************************************************************/
//  Programado por: Rebeca de la Fuente :  4/04/2022                                 
/*********************************************************************/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <iterator>
#include <map>
#include <limits>

typedef std::numeric_limits< double > dbl;

using namespace std;

#define  SIZE_LINE  5000


  struct  Data{
  string type;
  long int start;
  long int end;
  long int size;
  string ID;
  string Parent;
  long int region;
  long int id; 
  long int idm; 
  long int num_mRNA;
  long int num_CDS_portions;
  long int num_CDS;
  long int size_mRNA;
  long int coding;
  long int sumcds;
  long int size_CDS;
  long double gp;   
  int duplicates;
  };
  




int main( int argc, char *argv[])
{

//====================== INPUT GFF FILES   ================== //

 FILE    *fin;
 
 if (argc<2) {
     printf("\nUsa:\n> GRADO_PA  File_IN.GFF\n");
     exit(0);
     }

 if ( (fin=fopen(argv[1], "r"))==NULL) {
     printf ("\n Problems with %s ??.\n", argv[1]);
     fflush(stdin); getchar(); exit(0);
     }

//=====================  OUTPUT FILE    ================ //
 
 char  ficha2[200];
 ofstream fout;
 strcpy(ficha2, argv[1]);  strcat(ficha2, "_genome_variables.txt");
 fout.open(ficha2);

//=================== READ GFF or GFF3 FILES   ==================== //


 //Read GFF files and store data in structures 
 // Node types:
 // -> genes (nodes_genes)
 // -> mRNA (nodes_mRNA)
 // -> CDS (nodes_CDS)
 // Node attributes:
 // -> type
 // -> ID
 // -> Parent
 // -> start
 // -> end
 // -> region
 
 
 char *linea;
 linea = (char *) calloc (SIZE_LINE, sizeof(char));
 if (linea == NULL)    exit(EXIT_FAILURE);

 string word;
 string type;
 long int start,end;
 
 Data aux;
 vector<Data> nodes_genes;
 vector<Data> nodes_mRNA;
 vector<Data> nodes_CDS;
 vector<Data> region;
 
 long int genome_size;
 long int num_regions;
 long int num_genes;
 long int num_fcds;
 long int num_cds;


 string finalid;
 string finalidm;
 long int contadorid;
 
 num_genes=0;
 genome_size=0;
 num_regions=-1;
 int k=0;
 while (!feof(fin))
 { 
	
      fgets(linea, SIZE_LINE, fin);
      //---> skipping comment lines
      if (linea[0]!='#')
      { 
	
          stringstream str(linea);
          int column=1;
          while(getline(str, word, '\t'))
          {
                if(column==3)
                {
                   aux.type=word;
                }
                if(column==4)
                {
                   aux.start=stoi(word);
                }
                if(column==5)
                {
                   aux.end=stoi(word);
                }
                if(column==9)
                {
                   unsigned first = word.find("ID=");
                   unsigned last = word.find(";");
                   int ff=last;
                   string strNew = word.substr (first,last-first);
                   strNew.erase(0,3);
                   aux.ID=strNew;
                                         
                   if (word.find("Parent=") != string::npos )
                   {
                       first = word.find("Parent=");
                       last = word.find_first_of(";",first);
                       strNew = word.substr (first,last-first);
                       strNew.erase(0,7);
                       aux.Parent=strNew;
                       //cout<<"Yes, string contains the Parent<<endl;
                   }
                   else
                   {
                      aux.Parent="NoParent";
                      // cout<<"No, string do not contains the Parent<<endl;
                   }       
                }
                
          column++;
          }
          if(aux.type=="region" || aux.type=="chromosome" || aux.type=="supercontig")
          {
             genome_size=genome_size+aux.end-aux.start+1;
             num_regions++;
             region.push_back(aux);
          }          
          aux.region=num_regions;
          
        
          if(aux.type=="gene")
          {        
             aux.size=aux.end-aux.start+1;  
             finalid=aux.ID;
             contadorid=k;
             nodes_genes.push_back(aux);
             num_genes++;
             k++;
          }
          if(aux.type=="mRNA")
          {          
             if(finalid==aux.Parent)
             {
                aux.id=contadorid;
                finalidm=aux.ID;
                nodes_mRNA.push_back(aux);
             }
          }
          if(aux.type=="CDS")
          {
             if(finalid==aux.Parent || finalidm==aux.Parent)
             {
                aux.id=contadorid;
                nodes_CDS.push_back(aux);
             }          
          }
	
      }				  
 }
 num_regions++;


//============================== COMPUTATION GENOME VARIABLES   ============================ //

 
 for(int i=0;i<nodes_genes.size();i++)
 {
     nodes_genes[i].num_mRNA=0;
     nodes_genes[i].num_CDS=0;
     nodes_genes[i].num_CDS_portions=0;
     nodes_genes[i].size_mRNA=0;
     nodes_genes[i].size_CDS=0;
 }
 


 vector<string> libraryCDS[num_genes];
 vector<Data> libraryCDS_genes[num_genes];


 int esta;
 int id_aux;
 for(int i=0;i<nodes_CDS.size();i++)
 {

     id_aux=nodes_CDS[i].id;
     nodes_genes[id_aux].num_CDS_portions++;
     nodes_genes[id_aux].size_CDS=nodes_genes[id_aux].size_CDS+nodes_CDS[i].end-nodes_CDS[i].start+1;
     libraryCDS_genes[id_aux].push_back(nodes_CDS[i]); 
     esta=0;
     for(int k=0;k<libraryCDS[id_aux].size();k++)
     {
         if(libraryCDS[id_aux][k]==nodes_CDS[i].ID)
         {
            esta=1;
         }
     }
     if(esta==0)
     {
        libraryCDS[id_aux].push_back(nodes_CDS[i].ID);  
        nodes_genes[id_aux].num_CDS++;
        
     } 
     
 }
 
 for(int i=0;i<nodes_mRNA.size();i++)
 {
     id_aux=nodes_mRNA[i].id;
     nodes_genes[id_aux].num_mRNA++;
     nodes_genes[id_aux].size_mRNA=nodes_genes[id_aux].size_mRNA+nodes_mRNA[i].end-nodes_mRNA[i].start+1;     
 }
 
 

 int corrupto;  
 for(int i=0;i<nodes_genes.size();i++)
 {
     int aux_cds=0;
     long int i0=nodes_genes[i].start;
     long int i1=nodes_genes[i].end;
     long int s=nodes_genes[i].size;
     vector<int> nucleotidos(s);

     for(int k=0;k<s;k++)
     {
         nucleotidos[k]=0;
     }
     for(int j=0;j<libraryCDS_genes[i].size();j++)
     {
            corrupto=0;
            int ik0=libraryCDS_genes[i][j].start-i0;
            int ik1=libraryCDS_genes[i][j].end-i0;
            if(ik0<0 || ik1>s)
            {
               corrupto=1;
            }
            if(corrupto==0)
            {
               for(int k=ik0;k<ik1;k++)
               {
                   nucleotidos[k]++;
               } 
            }             
     }
     
     nodes_genes[i].coding=0;
     nodes_genes[i].sumcds=0;
     double aux_gp;
     for(int k=0;k<s;k++)
     {
         nodes_genes[i].sumcds=nodes_genes[i].sumcds+nucleotidos[k];
         if(nucleotidos[k]!=0)
         {
            nodes_genes[i].coding++;
         }
     }
     nodes_genes[i].gp=0;
     if(nodes_genes[i].coding>0)
     {
        aux_gp=nodes_genes[i].sumcds;
        aux_gp=aux_gp/nodes_genes[i].coding;
        nodes_genes[i].gp=aux_gp;
     }       
 }
 

 int num_genes_cds=0;
 double grado_paralelismo=0;
 int sum_coding=0;
 int sum_sumcds=0;
 double gp_media=0;
 double numnodes=0;
 long int size_genes=0;
 long int size_genes_cds=0;
 for(int i=0;i<nodes_genes.size();i++)
 {
     if(nodes_genes[i].num_CDS!=0)
     {
        num_genes_cds++;
        size_genes_cds=size_genes_cds+nodes_genes[i].size;
     }
     size_genes=size_genes+nodes_genes[i].size;
     sum_coding=sum_coding+nodes_genes[i].coding;
     sum_sumcds=sum_sumcds+nodes_genes[i].sumcds;
     if(nodes_genes[i].gp>0)
     {
        gp_media=gp_media+nodes_genes[i].gp;
        numnodes++;
     }
 }
 grado_paralelismo=sum_sumcds;
 grado_paralelismo=grado_paralelismo/sum_coding;
 gp_media=gp_media/numnodes;
 
 
 

 //Alternative splicing ratio
 
  vector<Data> CDSregion[num_regions];
  for(int i=0;i<nodes_CDS.size();i++)
  {
      int k=nodes_CDS[i].region;
      CDSregion[k].push_back(nodes_CDS[i]);
  }
  
  

  int codingr=0;
  int sumcdsr=0;
  double gpr=0;
  for(int i=0;i<num_regions;i++)
  {
  
      long int s=region[i].end-region[i].start+1;
      vector<int> nucleotidos(s);
      for(int k=0;k<s;k++)
      {
          nucleotidos[k]=0;
      }
      for(int j=0;j<CDSregion[i].size();j++)
      {
          corrupto=0; 
          long int i0=CDSregion[i][j].start;
          long int i1=CDSregion[i][j].end;
          if(i0<0 || i1>s)
          {
             corrupto=1;
          }
          if(corrupto==0)
          {
             for(int k=i0;k<i1;k++)
             {
                 nucleotidos[k]++;
             } 
          }         
      }
      
      for(int k=0;k<s;k++)
      { 
          sumcdsr=sumcdsr+nucleotidos[k];
          if(nucleotidos[k]!=0)
          {
             codingr++;
          }    
      }
  }
  gpr=sumcdsr;
  gpr=gpr/codingr;
 
//============================== DISPLAY RESULTS   ============================ //

 long double GenomeSize=genome_size;
 long double Coding=codingr;
 long double AlternativeSplicingRatio=gpr;
 long double GeneSize=size_genes;
  
 fout.precision(dbl::max_digits10);
 
 fout<<"#"<<argv[1]<<endl;
 fout<<"Genome Size: "<<GenomeSize<<endl;
 fout<<"Coding Size: "<<Coding<<endl;
 fout<<"Gene Size: "<<GeneSize<<endl;
 fout<<"ALternative Splicing Ratio: "<<AlternativeSplicingRatio<<endl;
 
 







 

}



