
typedef struct __ising_parm_struct
{
    double j_scalar;
    double bx_scalar;
    double bz_scalar;
}  IsingParmStruct;

enum Argument_Mode
{
    FROM_CSV,
    RANGE
};


template <class T>
T get_number_from_string(std::string input_string, std::string context = "no context provided."){
    double return_value;
    try
    {
        return_value = std::stod(input_string);
    }
    catch (const std::invalid_argument &e)
    {
        throw std::runtime_error("Invalid argument exception in get_number_from_string: " + context);
    }
    catch (const std::out_of_range &e)
    {
        throw std::runtime_error("Out of range exception in get_number_from_string: " + context);
    }
    return (T) return_value;

}

/* Reads the solve from a CSV*/
std::vector<IsingParmStruct> get_parms_from_csv(std::string file_name)
{
    std::vector<IsingParmStruct> output_parameters;
    std::ifstream file_descriptor(file_name);
    if (file_descriptor.good())
    {
        std::string line;
        int line_counter = 1;
        while (std::getline(file_descriptor, line))
        {
            double temp_j, temp_bx, temp_bz;
            std::string builder;
            int current_value = 0;

            for (int i = 0; i < line.length(); i++)
            {
                if (line.length() != 0)
                {
                    if (line[i] != ',' && i != line.length() - 1)
                    {
                        builder += line[i];
                    }
                    else
                    {
                        if (line.length() - 1 == i)
                            builder += line[i];
                        double temp = get_number_from_string<double>(builder);
                        builder = "";
                        if (current_value == 0)
                            temp_j = temp;
                        else if (current_value == 1)
                            temp_bx = temp;
                        else
                            temp_bz = temp;
                        current_value++;
                    }
                }
            }
            output_parameters.push_back((IsingParmStruct){.j_scalar = temp_j, .bx_scalar = temp_bx, .bz_scalar = temp_bz});
            line_counter++;
        }
    }
    else
    {
        throw std::runtime_error("Could not open file");
    }
    return output_parameters;
}

std::vector<IsingParmStruct> get_parms_from_range(double j, int num_bx, double init_bx, double stop_bx, int num_bz, double init_bz, double stop_bz)
{
    std::vector<IsingParmStruct> output_parameters;
    double current_row = init_bx;
    for (int bx_idx = 0; bx_idx < num_bx; bx_idx++)
    {
        double current_col = init_bz;
        for (int bz_idx = 0; bz_idx < num_bz; bz_idx++)
        {
            output_parameters.push_back((IsingParmStruct){.j_scalar = j, .bx_scalar = current_row, .bz_scalar = current_col});
            current_col += ((stop_bz - init_bz) / (num_bz - 1));
        }
        current_row += ((stop_bx - init_bx) / (num_bx - 1));
    }
    return output_parameters;
}


void process_debug_args(int debug_arg_start, int arg_count, char ** args, debug_flags * flags){

    for(int i = debug_arg_start; i < arg_count; i++ ){
        std::string current_arg(args[i]);
        if(current_arg == std::string("--save-components-petsc")){
            flags->save_component_matrices_petsc = true;
        }
        else if(current_arg == std::string("--save-components-ascii")){
            flags->save_component_matrices_ascii = true;
        }
        else if(current_arg == std::string("--print-eigenvalue-information")){
            flags->print_eigenvalue_information = true;
        }
        else if(current_arg == std::string("--run-performance-metrics")){
            flags->run_performance_metrics = true;
        }
        else if(current_arg  == std::string("--perform-eigenvalue-check")){
            flags->perform_eigenvalue_check = true;
        }
        else if(current_arg  == std::string("--verbose")){
            flags->verbose_mode = true;
        }
        else if(current_arg.substr(0,7) == "--file="){
            flags->filename = current_arg.substr(5);
        }
        else{
            std::cout << "Warning: Invalid argument \" " << current_arg << " \"." << std::endl;
        }
    }

}