/*
   Copyright 2017 University of Washington                          3-clause BSD license

   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

   3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <sys/stat.h>
#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <fstream>
#include <set>

using namespace std;

#include "mango-hash.h"

int phd_read_cmdline_enzyme_cut_params(char *argv[], enzyme_cut_params &params)
{
   string prot_file(argv[1]);
   int semi_tryptic = atoi(argv[2]);
   string precut_amino(argv[3]);
   string postcut_amino(argv[4]);
   string prenocut_amino(argv[5]);
   string postnocut_amino(argv[6]);
   int internal_lysine = atoi(argv[7]);
   string hash_file(argv[8]);

   // Validate the parameters: no intersection of pre and post and break and no-break

   cout << endl 
         << "Parameters of the run: " << endl
         << "Protein database file name is: " << prot_file << endl
         << "Tryptic/Semi-tryptic run: " << semi_tryptic << endl
         << "Pre break amino acids: " << precut_amino << endl
         << "Post break amino acids: " << postcut_amino << endl
         << "Pre no break amino acids: " << prenocut_amino << endl
         << "Post no break amino acids: " << postnocut_amino << endl
         << "Internal lysine (missed cleavage): " << internal_lysine << endl 
         << "Hash file given is: " << hash_file << endl 
         <<endl;

   params.precut_amino = precut_amino;
   params.postcut_amino = postcut_amino;
   params.prenocut_amino = prenocut_amino;
   params.postnocut_amino = postnocut_amino;
   params.missed_cleavage = internal_lysine;
   params.semi_tryptic = semi_tryptic;

   return 0;
}

int main(int argc, char *argv[]) 
{
   /*
    Usage: <protein_database> <missed_cleavage> <cut_aminoacid> <ignore_prolin> <peptide_hash_file>
   */

   if (argc != 9) {
      cout << "./a.out <protein_database> <semi/tryptic> <precut_aa> " 
            << "<postcut_aa> <prenocut_aa> <postnocut_aa> <internal_lysine> <saved_hash_file>" << endl;
      exit(1);
   }

   enzyme_cut_params enz_params;
   protein_hash_db_t phdp;
   if (!phd_read_cmdline_enzyme_cut_params(argv, enz_params)) {
      phdp = phd_retrieve_hash_db(argv[1], enz_params, argv[8]);   
   } else {
      cout << "Parameters are not proper" << endl;
      exit(1);
   }
   
   //retrieving peptides of mass 256
   int mass;
   for (int mass = 0; mass < 5000; mass++) {

      vector<peptide_hash_database::phd_peptide> *peptides = phdp->phd_get_peptides_ofmass(mass);

      if ((*peptides).size() <= 0) continue;
      for (peptide_hash_database::phd_peptide peptide : *peptides) {
         cout << "Peptide sequence of mass " << mass << ": " << peptide.phdpep_sequence() << " is contained in proteins ";

         for (int i = 0; i < peptide.phdpep_protein_list_size(); i++) {
             cout << " Name: " << peptide.phdpep_protein_list(i).phdpro_name() << " Id: " << peptide.phdpep_protein_list(i).phdpro_id();
         }

         cout << endl;
      }
   }

}
