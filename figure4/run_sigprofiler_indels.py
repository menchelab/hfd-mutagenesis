from SigProfilerExtractor import sigpro as sig
   # to get input from vcf files

data = "/users/mathilde.meyenberg/sigprofiler/INDEL2"

def main_function():
    sig.sigProfilerExtractor("vcf", "output_indel", data, reference_genome="mm10", opportunity_genome='mm10', context_type = "ID", exome = False,
    minimum_signatures=1, maximum_signatures=4, nmf_replicates=100, resample = True, gpu=False,
    nmf_init="random", precision= "single", matrix_normalization= "gmm", seeds= "random",
    min_nmf_iterations= 10000, max_nmf_iterations=1000000, nmf_test_conv= 10000, nmf_tolerance= 1e-15,
    cosmic_version=3.3, make_decomposition_plots=True)

if __name__=="__main__":
   main_function()
