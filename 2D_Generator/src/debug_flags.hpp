#ifndef PHS_DEBUG_FLAGS
#define PHS_DEBUG_FLAGS 1

class debug_flags{
    public:
    bool save_component_matrices_petsc;
    bool save_component_matrices_ascii;
    bool print_eigenvalue_information;
    bool run_performance_metrics;
    std::string filename;
    bool perform_eigenvalue_check;
    double eigenvalue_tolerance;
    bool verbose_mode;
    double eigenvector_output_tolerance;


    debug_flags(){
        save_component_matrices_petsc = false;
        save_component_matrices_ascii = false;
        print_eigenvalue_information = false;
        run_performance_metrics = false;
        filename = "results.eigenset";
        perform_eigenvalue_check = false;
        eigenvalue_tolerance = 0.000001;
        verbose_mode = false;
        eigenvector_output_tolerance = 0.00000001;
        
    }
};
#endif