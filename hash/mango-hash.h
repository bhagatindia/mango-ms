/*
   Copyright 2017 University of Washington                          3-clause BSD license

   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

   3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <string>

#include "protein_pep_hash.pb.h"

struct enzyme_cut_params {
   int         missed_cleavage;
   string      precut_amino;
   string      postcut_amino;
   string      prenocut_amino;
   string      postnocut_amino;
   int         semi_tryptic;
};

struct protein_hash_db_ {
   peptide_hash_database::phd_file phd_file_entry;
   vector<peptide_hash_database::phd_peptide>* phd_get_peptides_ofmass(int mass);
   vector<peptide_hash_database::phd_peptide>* phd_get_peptides_ofmass_tolerance(float mass_given, float tolerance);
   float phd_calculate_mass_peptide(const string peptide);
};

typedef protein_hash_db_* protein_hash_db_t;

protein_hash_db_t phd_retrieve_hash_db (const char *protein_file,
                                        enzyme_cut_params params,
                                        const char *phd_file);

