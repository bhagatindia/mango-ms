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
#include <cmath>
#include <algorithm>

#include "protein_pep_hash.pb.h"

using namespace std;

#include "mango-hash.h"
#include "../mango_DataInternal.h"

#define MIN_PEPTIDE_MASS 1
#define MAX_PEPTIDE_MASS 5000

#define MAX_EDGES 30

float pp_amino_acid_mass[MAX_EDGES] = {  
               71.037113805, //A 
               99999, //B not there
               160.03064805, //103.009184505 + 57.021464, //C
               115.026943065, //D
               129.042593135, //E
               147.068413945, //F
               57.021463735, //G
               137.058911875, //H
               113.084064015, //I
               99999, //J
               128.094963050, //K
               113.084064015, //L
               131.040484645, //M
               114.042927470, //N
               132.089877680, //O
               97.052763875, //P
               128.058577540, //Q
               156.101111050, //R
               87.032028435, //S
               101.047678505, //T
               150.95363, //U
               99.068413945, //V
               186.079312980, //W
               99999, //X
               163.063328575, //Y
               99999 //Z 
            };

float protein_hash_db_::phd_calculate_mass_peptide(const string peptide)
{
   float mass = 0;
   for (const char &c : peptide) {
      mass += pp_amino_acid_mass[c - 'A'];
   }
   return mass;
}

float phd_calculate_mass_peptide(const string peptide)
{
   float mass = 0;
   for (const char &c : peptide) {
      mass += pp_amino_acid_mass[c - 'A'];
   }
   return mass;
}

vector<peptide_hash_database::phd_peptide>* protein_hash_db_::phd_get_peptides_ofmass(int mass)
{
   // another memory leak
   vector<peptide_hash_database::phd_peptide> *ret = new vector<peptide_hash_database::phd_peptide>;

   peptide_hash_database::phd_peptide_mass pepm = phd_file_entry.phdpepm(mass);
   if (pepm.phdpmass_mass() == mass)
   {
      for (int i = 0; i < pepm.phdpmass_peptide_list_size(); i++)
      {
         ret->push_back(pepm.phdpmass_peptide_list(i));
      }
   }
   return ret;
}

#define PEP_WITHIN_TOLERANCE(mass_given, tolerance, mass_computed) \
  ((mass_computed <= (mass_given + tolerance)) && (mass_computed >= (mass_given - tolerance)))

vector<peptide_hash_database::phd_peptide>* protein_hash_db_::phd_get_peptides_ofmass_tolerance(float mass_given, float tolerance)
{
   // another memory leak
   vector<peptide_hash_database::phd_peptide> *ret = new vector<peptide_hash_database::phd_peptide>;

   int mass_min = floor(mass_given - tolerance);
   int mass_max = ceil(mass_given + tolerance);

   for (int mass = mass_min; mass <= mass_max; mass++)
   {
      peptide_hash_database::phd_peptide_mass pepm = phd_file_entry.phdpepm(mass);
      if (pepm.phdpmass_mass() == mass) {
         for (int i = 0; i < pepm.phdpmass_peptide_list_size(); i++) {
            float mass_computed = phd_calculate_mass_peptide(pepm.phdpmass_peptide_list(i).phdpep_sequence());
            if (PEP_WITHIN_TOLERANCE(mass_given, tolerance, mass_computed)) {
               ret->push_back(pepm.phdpmass_peptide_list(i));
            }
         }
      }
   }
   return ret;
}

// Define a free function in the library to free the memory

void phd_split_string(std::string str, std::string splitBy, std::vector<std::string>& tokens)
{
   /* Store the original string in the array, so we can loop the rest
    * of the algorithm. */
   tokens.push_back(str);

   // Store the split index in a 'size_t' (unsigned integer) type.
   size_t splitAt;
   // Store the size of what we're splicing out. One character is good enough.
   size_t splitLen = 1;
   // Create a string for temporarily storing the fragment we're processing.
   std::string frag;
   // Loop infinitely - break is internal.
   while(true)
   {
      /* Store the last string in the vector, which is the only logical
       * candidate for processing. */
      frag = tokens.back();
      /* The index where the split is. */
      splitAt = frag.find_first_of(splitBy);
      // If we didn't find a new split point...
      if(splitAt == string::npos)
      {
         // Break the loop and (implicitly) return.
         break;
      }
      /* Put everything from the left side of the split where the string
       * being processed used to be. */
      tokens.back() = frag.substr(0, splitAt + 1);
      /* Push everything from the right side of the split to the next empty
       * index in the vector. */
      tokens.push_back(frag.substr(splitAt+splitLen, frag.size()-(splitAt+splitLen)));
   }
}

void phd_read_protein_database(string file, peptide_hash_database::phd_file& pfile)
{
   // For every record we can read from the file, append it to our resulting data
   peptide_hash_database::phd_protein *record;
   int pcount = 0;

   std::ifstream fstream(file.c_str());
   std::string line, prd_id, pr_seq;   

   string split_char = " ";
   vector<string> split_tokens;

   while (std::getline(fstream, line)) {
      if (line.size() > 0 && line.at(0) == '>') {
         line.erase(0,1);                 // don't include '>' in protein name
         record = pfile.add_phdpro();
         pcount++;
         phd_split_string(line, split_char, split_tokens);
         record->set_phdpro_name(split_tokens[0]);
         record->set_phdpro_id(pcount);
      } else {
         // append or start a new one
         line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
         (record->mutable_phdpro_pepseq())->append(line);
      }
      split_tokens.clear();
   }

   return;
}

struct range {
   int start, length, missed, left, right;
};

void phd_basic_cut_pre_post(enzyme_cut_params params, const string protein_seq,
                              vector<range *> &splits)
{
   const string spres = params.precut_amino;
   const string sposts = params.postcut_amino;
   size_t pos = 0, length = protein_seq.length(), npres, nposts, cpos, start = pos;

   //cout << "Protein length: " << length << " and pre split string " << spres << " and post split string is " << sposts << " and position is " << pos << endl;
   //cout << "Protein sequence is: " << protein_seq << endl;

   while ((pos != std::string::npos) && (start < length)) {
      range *r = new range; r->start = pos;

      cpos = protein_seq.find_first_of(spres + sposts, start);

      npres = protein_seq.find_first_of(spres, start);
      nposts = protein_seq.find_first_of(sposts, start);

      if (cpos == std::string::npos) {
         r->length = length - r->start;
         splits.push_back(r);
         break;
      } else if (cpos == npres) {
         pos = cpos;
         r->length = pos - r->start;
         start = pos + 1;
      } else if (cpos == nposts) {
         // Check for K followed by Proline
         pos = cpos + 1;
         r->length = pos - r->start;
         start = pos;
      }
      //cout << "Splits are " << npres << " " << nposts << " " << " final position is " << pos << endl;

      splits.push_back(r);
   }

   // Now we get the basic peptides
/* 
   cout << "We got the basic peptides based on the pre and post conditions" << endl;

   for (range *r: splits) {
      cout << "Range: Start is " << r->start << " and length is " << r->length << " " << protein_seq.substr(r->start, r->length) << endl;
   }
*/

}

void phd_handle_post_merge(enzyme_cut_params params, const string protein_seq,
                           vector<range *> &splits, vector<range *> &nocut_splits)
{
   // We have to take care of pre nocut amino acids and postcut amino acids and merge
   const string nocut_spres = params.prenocut_amino;
   const string nocut_spost = params.postnocut_amino;

   vector<range *>::iterator split_iterator;

   //cout << "Post no cut " << nocut_spost << endl;

   for (split_iterator = splits.begin(); split_iterator != splits.end(); split_iterator++) {
      range *r = *split_iterator;

      range *r_new = new range;
      r_new->start = r->start;
      r_new->length = r->length;
      r_new->missed = r_new->left = r_new->right = 0;

      while (split_iterator + 1 != splits.end()) {
         range *r_next = *(split_iterator + 1);
         //cout << "Split is " << r_next->start << " length is " << r_next->length << " the character is " << (protein_seq.c_str())[r_next->start] << endl;
         if (nocut_spost.find_first_of((protein_seq.c_str())[r_next->start]) != std::string::npos) {
            //cout << "Found a split" << endl;
            r_new->length += r_next->length;
            split_iterator++;
         } else break;
      }
      
      nocut_splits.push_back(r_new);
   }
   
/*
   cout << "No cut splits" << endl;

   for (range *r: nocut_splits) {
      cout << "Range: Start is " << r->start << " and length is " << r->length << " " << protein_seq.substr(r->start, r->length) << endl;
   }
*/

}

void phd_handle_missed_cleavage(enzyme_cut_params params, const string protein_seq,
                                 vector<range *> &nocut_splits, vector<range *> &final_splits)
{
   // We have to handle missing cleavage
   int missed_cleavage = params.missed_cleavage;
   //cout << "Missed cleavage " << missed_cleavage << endl;

   if (missed_cleavage) {
      int num_cuts = nocut_splits.size();
      for (int i = 0; i < num_cuts; i++) {
         range *r = nocut_splits.at(i);
	 // Before merging, check for the internal K: if so, dont merge, just emit and continue
	 string peptide = protein_seq.substr(r->start, r->length);
	 if (((int)peptide.find_first_of("K") != r->length - 1) && (peptide.find_first_of("K") != std::string::npos)) {
            range *r_new = new range;
            r_new->start = r->start;
            r_new->length = r->length;
            r_new->missed = peptide.find_first_of("K");
            r_new->left = r_new->right = 0;
            final_splits.push_back(r_new);
	    continue;
	 }

	 if (i+1 >= num_cuts) continue;

         range *r_next = nocut_splits.at(i+1);
         //cout << "Peeking at the end " << (protein_seq.c_str())[r->start + r->length -1] << endl;
         if ((protein_seq.c_str())[r->start + r->length -1] == 'K') {
            range *r_new = new range;
            r_new->start = r->start;
            r_new->length = r->length + r_next->length;
            r_new->missed = r->length;
            r_new->left = r_new->right = 0;
            final_splits.push_back(r_new);
         }
      }
      // Add the last one to the final splits
      /*
      range *r = nocut_splits.at(num_cuts - 1);
      range *r_new = new range;
      r_new->start = r->start;
      r_new->length = r->length;
      r_new->missed = r_new->left = r_new->right = 0;
      final_splits.push_back(r_new);
      */
   } else {
      for (range *r : nocut_splits) {
         range *r_new = new range;
         r_new->start = r->start;
         r_new->length = r->length;
         r_new->missed = r_new->left = r_new->right = 0;
         final_splits.push_back(r_new);
      }
   }

/*
   for (range *r: final_splits) {
      cout << "Range: Start is " << r->start << " and length is " << r->length << 
               " and missed cleavage " << r->missed << " Left is " << r->left << 
               " Right is " << r->right << " and the peptide is " << protein_seq.substr(r->start, r->length) << endl;
   }
*/

}

void phd_add_missed_semi_tryptic(range* r, vector<range *> &nocut_splits)
{
   // Left tryptic
   //cout << "Left trypic peptides for missed cleavages " << r->start << " " << r->length << " " << r->missed << endl;
   for (int pos = r->start + r->missed + 1; pos < (r->start + r->length); pos++) {
      range *r_new = new range;
      r_new->start = r->start;
      r_new->length= pos - r->start;
      r_new->missed = 1;
      r_new->right = 0;
      r_new->left = 1;
      nocut_splits.push_back(r_new);
   }

   // Right tryptic
   for (int start = r->start + 1; start < (r->start + r->missed); start++) {
      range *r_new = new range;
      r_new->start = start;
      r_new->length= r->start + r->length - start;
      r_new->missed = 1;
      r_new->right = 1;
      r_new->left = 0;
      nocut_splits.push_back(r_new);
   } 
}

void phd_add_semi_tryptic( range* r, vector<range *> &nocut_splits)
{
   // Left tryptic
   int start = r->start;
   for (int len = 1; len < r->length; len++) {
      range *r_new = new range;
      r_new->start = r->start;
      r_new->length= len;
      r_new->missed = 0;
      r_new->right = 0;
      r_new->left = 1;
      nocut_splits.push_back(r_new);
   }

   // right tryptic
   for (start = r->start + 1; start < (r->start + r->length -1); start++) {
      range *r_new = new range;
      r_new->start = start;
      r_new->length= r->start + r->length - start;
      r_new->missed = 0;
      r_new->right = 1;
      r_new->left = 0;
      nocut_splits.push_back(r_new);
   } 
}

void phd_handle_semi_tryptic(enzyme_cut_params params, vector<range *> &nocut_splits)
{
   // We have to add sem-tryptic peptides
   int semi_tryptic = params.semi_tryptic;
   //cout << "Semi-tryptic " << semi_tryptic << endl;

   if (semi_tryptic) {
      int num_count = nocut_splits.size();
      for (int i = 0; i < num_count; i++) {
         if ((nocut_splits.at(i))->missed) {
            phd_add_missed_semi_tryptic(nocut_splits.at(i), nocut_splits);
         } else {
            phd_add_semi_tryptic(nocut_splits.at(i), nocut_splits);
         }
      }
   }

}

void phd_add_peptide_into_hash (string peptide, 
                                 const peptide_hash_database::phd_protein pro_seq,
                                 peptide_hash_database::phd_peptide_mass *pfile_pepm)
{
   int found = 0;
   peptide_hash_database::phd_peptide *found_peptide;
   // Loop over the list of peptides to find whether this peptide is always first
   for (int i = 0; i < pfile_pepm->phdpmass_peptide_list_size(); i++)
   {
       if (!strcmp(pfile_pepm->phdpmass_peptide_list(i).phdpep_sequence().c_str(), 
                    peptide.c_str())) {
           if (found) {
              cout << "A duplicate found when it is not supposed to be. So, panicking" << endl;
              exit(1);
           }
           found = 1;
           found_peptide = pfile_pepm->mutable_phdpmass_peptide_list(i);
       }
   }

   if (!found) {
       found_peptide = pfile_pepm->add_phdpmass_peptide_list();
       found_peptide->set_phdpep_sequence(peptide);
   }

   // Add the protein to the found peptide list
   peptide_hash_database::phd_protein* fp = found_peptide->add_phdpep_protein_list();
   fp->set_phdpro_name(pro_seq.phdpro_name());
   fp->set_phdpro_id(pro_seq.phdpro_id());
}

void phd_split_protein_sequence_peptides(enzyme_cut_params params, 
                                          const peptide_hash_database::phd_protein pro_seq,
                                          peptide_hash_database::phd_file &pfile)
{
   //FIXME: Clean up memory. There are memory leaks
   const string protein_seq = pro_seq.phdpro_pepseq();
   vector<range *> splits, nocut_splits, final_splits;

   phd_basic_cut_pre_post(params, protein_seq, splits);

   phd_handle_post_merge(params, protein_seq, splits, nocut_splits);

   phd_handle_missed_cleavage(params, protein_seq, nocut_splits, final_splits);

   phd_handle_semi_tryptic(params, final_splits);

   for (range *r: final_splits) {
      string peptide = protein_seq.substr(r->start, r->length);
      int mass = phd_calculate_mass_peptide(peptide);
      /*
      cout << "Range: Start is " << r->start << " and length is " << r->length << 
               " and missed cleavage " << r->missed << " Left is " << r->left << 
               " Right is " << r->right << " and the peptide is " << peptide << 
               " and its mass is " << mass << endl;
      */
      if (MIN_PEPTIDE_MASS < mass && mass < MAX_PEPTIDE_MASS) {
         phd_add_peptide_into_hash(peptide, pro_seq, pfile.mutable_phdpepm(mass));
      }
   }

   // Cleanup splits
   for (range *r: splits) delete r;

   // Clean up nocut splits
   for (range *r: nocut_splits) delete r;

   // Clean up final splits
   for (range *r: final_splits) delete r;
}

void phd_add_peptide_hash_database (peptide_hash_database::phd_file &pfile, 
                                    enzyme_cut_params cut_params)
{
   const peptide_hash_database::phd_header hdr = pfile.phdhdr();
// char szTmp[24];

   for (int i = 0; i < MAX_PEPTIDE_MASS; i++) {
      peptide_hash_database::phd_peptide_mass *pepm = pfile.add_phdpepm();
      pepm->set_phdpmass_mass(i);
   }

// cout << "Splitting peptides for protein id: ";
   for (int i = 0; i < pfile.phdpro_size(); i++) {
//    sprintf(szTmp, "%d", i);
//    cout << szTmp;
//    cout.flush();
//    for (int ii=0; ii<strlen(szTmp);ii++)
//       cout << '\b';
      phd_split_protein_sequence_peptides(cut_params, pfile.phdpro(i), pfile);
   }
// cout << endl;

   // If the parameter is semi-tryptic, add all left and right semi-tryptic peptides
}

void phd_save_hash_db(peptide_hash_database::phd_file &pfile, const char *hash_file)
{
   /*
   peptide_hash_database::phd_header hdr;

   // Create header which is mandatory
   hdr = pfile.phdhdr();
   */
// cout << "Creating file" << endl;

   std::ofstream ofs;
   ofs.open (hash_file, ios::out | ios::trunc | ios::binary);
   if (!pfile.SerializeToOstream(&ofs)) {
      cout << "Cannot write to hash file" << endl;
      exit(1);
   }

// cout << "Wrote file" << endl;
   ofs.close();
}

void phd_print_hash_file_params(peptide_hash_database::phd_header hdr)
{
   cout << endl 
         << "Parameters of the file: " << endl
         << "Protein database file name is: " << hdr.phdhdr_precut_amino()<< endl
         << "Tryptic/Semi-tryptic run: " << hdr.phdhdr_semi_tryptic() << endl
         << "Pre break amino acids: " << hdr.phdhdr_precut_amino() << endl
         << "Post break amino acids: " << hdr.phdhdr_postcut_amino() << endl
         << "Pre no break amino acids: " << hdr.phdhdr_prenocut_amino() << endl
         << "Post no break amino acids: " << hdr.phdhdr_postnocut_amino() << endl
         << "Internal lysine (missed cleavage): " << hdr.phdhdr_missed_cleavage() << endl 
         <<endl;

}

void phd_populate_hdr_params(enzyme_cut_params params, 
                              peptide_hash_database::phd_header *phdr)
{
   phdr->set_phdhdr_version(1);
   phdr->set_phdhdr_protein_source_filename("DEADBEEF");
   phdr->set_phdhdr_protein_source_file_digest("DEADBEEF");
   phdr->set_phdhdr_num_proteins(10);
   phdr->set_phdhdr_hash_file_name("DEADBEEF");
   phdr->set_phdhdr_hash_file_digest("DEADBEEF");

   phdr->set_phdhdr_precut_amino(params.precut_amino);
   phdr->set_phdhdr_postcut_amino(params.postcut_amino);
   phdr->set_phdhdr_prenocut_amino(params.prenocut_amino);
   phdr->set_phdhdr_postnocut_amino(params.postnocut_amino);
   phdr->set_phdhdr_missed_cleavage(params.missed_cleavage);
   phdr->set_phdhdr_semi_tryptic(params.semi_tryptic);
}

void phd_create_hash_file (const char *protein_file, enzyme_cut_params params, 
                              const char *phd_file)
{
   peptide_hash_database::phd_file pfile;

// cout << "Read protein database" << endl;
   phd_read_protein_database(protein_file, pfile);

// cout << "Creating peptides based on the enzyme digestion criteria" << endl;
   phd_add_peptide_hash_database(pfile, params);

// cout << "Populating header parameters" << endl;
   phd_populate_hdr_params(params, pfile.mutable_phdhdr());

   phd_save_hash_db(pfile, phd_file);
}

int phd_load_hash_file (const char *hash_file, peptide_hash_database::phd_file &pfile)
{
// cout << "File name is " << hash_file << endl;
   fstream input(hash_file, ios::in | ios::binary);
   if (!input) {
      cout << hash_file << ": File not found.  Creating a new file." << endl;
      return 1;
   }

   if (pfile.ParseFromIstream(&input)) {
//    cout << "Read the hash file: " << hash_file << endl;
      peptide_hash_database::phd_header phdr = pfile.phdhdr();
//    phd_print_hash_file_params(phdr);
      return 0;
   } else {
      cout << "File read has problems" << endl;
   }

   return 1;
}

int phd_compare_enzyme_cut_params_hash_file (enzyme_cut_params params,
                                    peptide_hash_database::phd_header phdr)
{
   int ret_value = 1;

// cout << "Precut amino acid: " << params.precut_amino << " file " << phdr.phdhdr_precut_amino() << endl;
   if (params.precut_amino.compare(phdr.phdhdr_precut_amino())) ret_value = 0;

// cout << "Postcut amino acid: " << params.postcut_amino << " file " << phdr.phdhdr_postcut_amino() << endl;
   if (params.postcut_amino.compare(phdr.phdhdr_postcut_amino())) ret_value = 0;

// cout << "Pre nocut amino acid: " << params.prenocut_amino << " file " << phdr.phdhdr_prenocut_amino() << endl;
   if (params.prenocut_amino.compare(phdr.phdhdr_prenocut_amino())) ret_value = 0;

// cout << "Post nocut amino acid: " << params.postnocut_amino << " file " << phdr.phdhdr_postnocut_amino() << endl;
   if (params.postnocut_amino.compare(phdr.phdhdr_postnocut_amino())) ret_value = 0;

// cout << "Missed cleavage: " << params.missed_cleavage << " file " << phdr.phdhdr_missed_cleavage() << endl;
   if ((params.missed_cleavage != phdr.phdhdr_missed_cleavage())) ret_value = 0;

// cout << "Semi tryptic : " << params.semi_tryptic << " file " << phdr.phdhdr_semi_tryptic() << endl;
   if ((params.semi_tryptic != phdr.phdhdr_semi_tryptic())) ret_value = 0;

   return ret_value;
}

inline int phd_hash_file_match_criteria (enzyme_cut_params params, const char *phd_file)
{
   peptide_hash_database::phd_file pfile;
   int ret_value = 0;

// cout << "File name is " << phd_file << endl;
   fstream input(phd_file, ios::in | ios::binary);
   if (!input) {
      cout << " " << phd_file << ": File not found.  Creating a new file." << endl;
      return 0;
   }

   if (pfile.ParseFromIstream(&input)) {
//    cout << "Reading the hash file and comparing the header param: " << phd_file << endl;
      peptide_hash_database::phd_header phdr = pfile.phdhdr();
      ret_value = phd_compare_enzyme_cut_params_hash_file(params, phdr);
//    if (ret_value)ucout << "Comparing returned equal parameters" << endl;
//    else cout << "Comparing returned non-equal parameters" << endl;
   }

   input.close();
   return ret_value;
}
      
protein_hash_db_t phd_retrieve_hash_db (const char *protein_file, 
                              enzyme_cut_params params, const char *phd_file)
{
   GOOGLE_PROTOBUF_VERIFY_VERSION;

//   cout << "In phd_retrieve_hash_db" << endl;

   if (!phd_hash_file_match_criteria(params, phd_file)) {
      // Check the existance of both files
      // Parse the file and compare with the enzyme cut params and the digests
      // If something doesnt match, create one and return entry
//    cout << "Creating hash file" << endl;
      phd_create_hash_file(protein_file, params, phd_file);
//    cout << "Created hash database" << endl;
   }

   protein_hash_db_t ret_entry = new protein_hash_db_;
// cout << "Loading hash database" << endl;
   phd_load_hash_file(phd_file, ret_entry->phd_file_entry);

   return ret_entry;
}

