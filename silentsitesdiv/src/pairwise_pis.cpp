// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_map>
#include <string>
using namespace Rcpp;

// --------------------------------------------------------------------------------
// (1) “num” and “aa” table from kaks.c, to map a 3‐letter codon → 0..63 and then → amino acid
// --------------------------------------------------------------------------------

static int num_codon(const char *cod) {
  // Return 0..63 for a valid codon (ACGT only), or 64 if invalid.
  static const char bases[] = "ACGT";
  if (strchr(bases, cod[0]) == NULL ||
      strchr(bases, cod[1]) == NULL ||
      strchr(bases, cod[2]) == NULL)
    return 64;  // indicates “skip this codon”

  int n1 = 0, n2 = 0, n3 = 0;
  if (cod[0] == 'C') n1 = 1;
  if (cod[0] == 'G') n1 = 2;
  if (cod[0] == 'T') n1 = 3;
  if (cod[1] == 'C') n2 = 1;
  if (cod[1] == 'G') n2 = 2;
  if (cod[1] == 'T') n2 = 3;
  if (cod[2] == 'C') n3 = 1;
  if (cod[2] == 'G') n3 = 2;
  if (cod[2] == 'T') n3 = 3;
  return 16 * n1 + 4 * n2 + n3;
}

// We will fill this array once at startup:
static int aa_of_codon[64] = {0};
static bool aa_table_initialized = false;

static void init_aa_table() {
  if (aa_table_initialized) return;
  aa_table_initialized = true;

  // Copy‐paste of “aa[0]=… aa[63]=…” from kaks.c’s prefastlwl() (standard code)
  aa_of_codon[0]  = 17;  // AAA → Lys
  aa_of_codon[1]  = 16;  // AAC → Asn
  aa_of_codon[2]  = 17;  // AAG → Lys
  aa_of_codon[3]  = 16;  // AAT → Asn
  aa_of_codon[4]  = 13;  // ACA → Thr
  aa_of_codon[5]  = 13;  // ACC → Thr
  aa_of_codon[6]  = 13;  // ACG → Thr
  aa_of_codon[7]  = 13;  // ACT → Thr
  // Note: in kaks.c, code_mt can switch some values; here we assume standard (not mitochondrial).
  aa_of_codon[8]  = 18;  // AGA → Arg
  aa_of_codon[9]  = 14;  // AGC → Ser
  aa_of_codon[10] = 18;  // AGG → Arg
  aa_of_codon[11] = 14;  // AGT → Ser
  aa_of_codon[12] = 7;   // ATA → Ile
  aa_of_codon[13] = 7;   // ATC → Ile
  aa_of_codon[14] = 5;   // ATG → Met
  aa_of_codon[15] = 7;   // ATT → Ile
  aa_of_codon[16] = 15;  // CAA → Gln
  aa_of_codon[17] = 4;   // CAC → His
  aa_of_codon[18] = 15;  // CAG → Gln
  aa_of_codon[19] = 4;   // CAT → His
  aa_of_codon[20] = 9;   // CCA → Pro
  aa_of_codon[21] = 9;   // CCC → Pro
  aa_of_codon[22] = 9;   // CCG → Pro
  aa_of_codon[23] = 9;   // CCT → Pro
  aa_of_codon[24] = 18;  // CGA → Arg
  aa_of_codon[25] = 18;  // CGC → Arg
  aa_of_codon[26] = 18;  // CGG → Arg
  aa_of_codon[27] = 18;  // CGT → Arg
  aa_of_codon[28] = 6;   // CTA → Leu
  aa_of_codon[29] = 6;   // CTC → Leu
  aa_of_codon[30] = 6;   // CTG → Leu
  aa_of_codon[31] = 6;   // CTT → Leu
  aa_of_codon[32] = 19;  // GAA → Glu
  aa_of_codon[33] = 20;  // GAC → Asp
  aa_of_codon[34] = 19;  // GAG → Glu
  aa_of_codon[35] = 20;  // GAT → Asp
  aa_of_codon[36] = 11;  // GCA → Ala
  aa_of_codon[37] = 11;  // GCC → Ala
  aa_of_codon[38] = 11;  // GCG → Ala
  aa_of_codon[39] = 11;  // GCT → Ala
  aa_of_codon[40] = 12;  // GGA → Gly
  aa_of_codon[41] = 12;  // GGC → Gly
  aa_of_codon[42] = 12;  // GGG → Gly
  aa_of_codon[43] = 12;  // GGT → Gly
  aa_of_codon[44] = 8;   // GTA → Val
  aa_of_codon[45] = 8;   // GTC → Val
  aa_of_codon[46] = 8;   // GTG → Val
  aa_of_codon[47] = 8;   // GTT → Val
  aa_of_codon[48] = 0;   // TAA → STOP
  aa_of_codon[49] = 3;   // TAC → Tyr
  aa_of_codon[50] = 0;   // TAG → STOP
  aa_of_codon[51] = 3;   // TAT → Tyr
  aa_of_codon[52] = 14;  // TCA → Ser
  aa_of_codon[53] = 14;  // TCC → Ser
  aa_of_codon[54] = 14;  // TCG → Ser
  aa_of_codon[55] = 14;  // TCT → Ser
  aa_of_codon[56] = 2;   // TGA → Trp  (standard code)
  aa_of_codon[57] = 10;  // TGC → Cys
  aa_of_codon[58] = 2;   // TGG → Trp
  aa_of_codon[59] = 10;  // TGT → Cys
  aa_of_codon[60] = 6;   // TTA → Leu
  aa_of_codon[61] = 1;   // TTC → Phe
  aa_of_codon[62] = 6;   // TTG → Leu
  aa_of_codon[63] = 1;   // TTT → Phe
}

// --------------------------------------------------------------------------------
// (2) Precompute “synonymous sites” for each codon (Li 1985 approx):
//     For each codon idx, for each of its 3 positions, see how many single‐base
//     changes leave the amino acid unchanged.  Then divide total by 3.
// --------------------------------------------------------------------------------

static double syn_sites_per_codon[64] = {0.0};
static bool syn_table_initialized = false;

static void init_syn_sites() {
  if (syn_table_initialized) return;
  syn_table_initialized = true;

  init_aa_table();  // make sure aa_of_codon[] is populated

  // For each of the 64 codons, build its 3‐letter string and count how many
  // single‐nucleotide substitutions remain in the same amino acid class.
  static const char *bases = "ACGT";
  for (int idx = 0; idx < 64; idx++) {
    // Reconstruct the 3‐letter codon string from idx:
    int n1 = idx / 16;
    int n2 = (idx % 16) / 4;
    int n3 = idx % 4;
    char codon[4];
    codon[0] = "ACGT"[n1];
    codon[1] = "ACGT"[n2];
    codon[2] = "ACGT"[n3];
    codon[3] = '\0';
    int aa1 = aa_of_codon[idx];
    double count_syn = 0.0;

    // For each of the 3 positions in this codon, try substituting each of the
    // other 3 bases and see if the AA stays the same:
    for (int pos = 0; pos < 3; pos++) {
      for (int b = 0; b < 4; b++) {
        char alt = bases[b];
        if (alt == codon[pos]) continue;
        char alt_codon[4] = { codon[0], codon[1], codon[2], '\0' };
        alt_codon[pos] = alt;
        int idx2 = num_codon(alt_codon);
        if (idx2 < 64 && aa_of_codon[idx2] == aa1) {
          count_syn += 1.0;
        }
      }
    }
    // Divide by 3 to get the “synonymous site” count for this codon (Li 1985)
    syn_sites_per_codon[idx] = count_syn / 3.0;
  }
}

// --------------------------------------------------------------------------------
// (3) Helper: count nucleotide mismatches between two codons (assumed length 3)
// --------------------------------------------------------------------------------

static inline int codon_nuc_diff(const char *c1, const char *c2) {
  int diff = 0;
  if (c1[0] != c2[0]) diff++;
  if (c1[1] != c2[1]) diff++;
  if (c1[2] != c2[2]) diff++;
  return diff;
}

// --------------------------------------------------------------------------------
// (4) The R‐callable function: takes a CharacterVector of codon‐aligned DNA,
//     returns an R “dist” of pairwise πₛ.  (Lower‐triangle ordering.)
// --------------------------------------------------------------------------------

// [[Rcpp::export]]
SEXP pairwise_piS(CharacterVector r_seqs) {
  int n = r_seqs.size();
  if (n < 2) {
    stop("Need at least two sequences.");
  }

  // Convert to std::vector<std::string> (uppercase)
  std::vector<std::string> seqs(n);
  for (int i = 0; i < n; i++) {
    seqs[i] = Rcpp::as<std::string>(r_seqs[i]);
    // Convert to uppercase to match our num_codon logic
    for (auto &ch : seqs[i]) {
      if (ch >= 'a' && ch <= 'z') ch = ch - 'a' + 'A';
    }
  }

  // Check all sequences are same length and multiple of 3
  int L = seqs[0].size();
  if (L % 3 != 0) {
    stop("Sequences must be codon‐aligned (length multiple of 3).");
  }
  for (int i = 1; i < n; i++) {
    if ((int)seqs[i].size() != L) {
      stop("All sequences must have the same length.");
    }
  }

  // Ensure our lookup tables are initialized:
  init_aa_table();
  init_syn_sites();

  int nCodons = L / 3;
  int out_len = (n * (n - 1)) / 2;
  NumericVector out(out_len);  // holds lower‐triangle πₛ
  int pos_out = 0;

  // Temporary buffers for codons:
  char codon1[4], codon2[4];
  codon1[3] = codon2[3] = '\0';

  // Main pairwise loops (i < j), in “dist” order: for i=0..n-2, for j=i+1..n-1
  for (int i = 0; i < n - 1; i++) {
    const std::string &s1 = seqs[i];
    for (int j = i + 1; j < n; j++) {
      const std::string &s2 = seqs[j];

      double syn_site_sum = 0.0;
      int syn_diff_sum = 0;

      // Walk codon by codon
      for (int c = 0; c < nCodons; c++) {
        // Extract 3‐letter codon from each sequence:
        codon1[0] = s1[3*c + 0];
        codon1[1] = s1[3*c + 1];
        codon1[2] = s1[3*c + 2];
        codon2[0] = s2[3*c + 0];
        codon2[1] = s2[3*c + 1];
        codon2[2] = s2[3*c + 2];

        // Map to indices 0..63 (or 64 if invalid)
        int idx1 = num_codon(codon1);
        int idx2 = num_codon(codon2);

        // Skip if either codon contains non-ACGT or is a stop/invalid
        if (idx1 == 64 || idx2 == 64) {
          continue;
        }

        // Add up synonymous “sites” from codon1
        syn_site_sum += syn_sites_per_codon[idx1];

        // If codons differ AND encode the same amino acid → synonymous differences
        if (idx1 != idx2 && aa_of_codon[idx1] == aa_of_codon[idx2]) {
          // Count how many nucleotides differ in this codon
          syn_diff_sum += codon_nuc_diff(codon1, codon2);
        }
      }

      // Compute πₛ = (observed syn_differences)/(synonymous_sites)
      double piS = R_NaReal;
      if (syn_site_sum > 0.0) {
        piS = ((double) syn_diff_sum) / syn_site_sum;
      }
      out[pos_out++] = piS;
    }
  }

  // Now wrap “out” into a dist object
  SEXP res = PROTECT(out);
  Rf_setAttrib(res, Rf_install("Size"), wrap(n));
  Rf_setAttrib(res, Rf_install("Labels"), r_seqs.names());
  Rf_setAttrib(res, Rf_install("Diag"), wrap(false));
  Rf_setAttrib(res, Rf_install("Upper"), wrap(false));
  Rf_setAttrib(res, Rf_install("method"), Rf_mkString("pairwise_piS"));

  // Make class “dist”
  Rf_setAttrib(res, Rf_install("class"), Rf_mkString("dist"));
  UNPROTECT(1);
  return res;
}
