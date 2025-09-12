#  Created by Sam Champer, 2020.
#  A product of the Messer Lab, http://messerlab.org/slim/

#  Sam Champer, Ben Haller and Philipp Messer, the authors of this code, hereby
#  place the code in this file into the public domain without restriction.
#  If you use this code, please credit SLiM-Extras and provide a link to
#  the SLiM-Extras repository at https://github.com/MesserLab/SLiM-Extras.
#  Thank you.

# This is an example of how to use Python as a driver for SLiM.
# Output is formated as a csv, but just printed to stdout by default.

# In order to reconfigure this file for your research project, the
# run_slim() and configure_slim_command_line() functions do not need to be modified.
# Changes you would likely want to make are to the argument parser in main(),
# in order to pass your desired variables to SLiM, and to the parse_slim()
# function, where you could do your desired operations on the output of SLiM.



#try: 10 times each parameter

from argparse import ArgumentParser
import subprocess


def parse_slim(slim_string):
    """
    Parse the output of SLiM to extract whatever data we're looking for.
    If we want to do a more complex analysis on the output of the SLiM file,
    this is where we do it.
    Args:
        slim_string: the entire output of a run of SLiM.
    Return
        output: the desired output we want from the SLiM simulation.
    """
    # The example SLiM file has been configured such that all the
    # output we want is printed on lines that start with "OUT:"
    # so we'll discard all other output lines."
    output = ""
    lines = slim_string.split('OUT')
    del lines[0]
    for line in lines:
        if line.startswith(":"):
            output += line.split(":")[1]
    return output


def run_slim(command_line_args):
    """
    Runs SLiM using subprocess.
    Args:
        command_line_args: list; a list of command line arguments.
    return: The entire SLiM output as a string.
    """
    slim = subprocess.Popen(command_line_args,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            universal_newlines=True)
    out, err = slim.communicate()
    # For debugging purposes:
    # std.out from the subprocess is in slim.communicate()[0]
    # std.error from the subprocess is in slim.communicate()[1]
    # Errors from the process can be printed with:
    # print(err)
    return out


def configure_slim_command_line(args_dict):
    """
    Sets up a list of command line arguments for running SLiM.
    Args:
        args_dict: a dictionary of arg parser arguments.
    Return
        clargs: A formated list of the arguments.
    """
    # We're running SLiM, so the first arg is simple:
    clargs = "slim "
    # The filename of the source file must be the last argument:
    source = args_dict.pop("source")
    # Add each argument from arg parser to the command line arguemnts for SLiM:
    for arg in args_dict:
        if isinstance(args_dict[arg], bool):
            clargs += f"-d {arg}={'T' if args_dict[arg] else 'F'} "
        else:
            clargs += f"-d {arg}={args_dict[arg]} "
    # Add the source file, and return the string split into a list.
    clargs += source
    return clargs.split()


def main():
    """
    1. Configure using argparse.
    2. Generate the command line list to pass to subprocess through the run_slim() function.
    3. Run SLiM.
    4. Process the output of SLiM to extract the information we want.
    5. Print the results.
    """
    # Get args from arg parser:
    parser = ArgumentParser()
    
    parser.add_argument('-src', '--source', default="pan_1d_daisy_TRUE_som0823.slim", type=str,
                        help=r"SLiM file to be run. Default 'pan_1d_daisy_TRUE_som0823.slim'")
    
    
    parser.add_argument('-header', '--print_header', action='store_true', default=False,
                        help='If this is set, python prints a header for a csv file.')

    parser.add_argument('-drop_size', '--DROP_SIZE', default=0.5, type=float,
                        help='The drop size of daisy drive . Default 0.5.')

    parser.add_argument('-somatic_fitness_f', '--SOMATIC_FITNESS_MUTLIPLIER_F', default=1.0, type=float,
                        help='The SOMATIC_FITNESS_MUTLIPLIER_F of daisy drive . Default 1.0.')
    #parser.add_argument('-drop_radius', '--DROP_RADIUS', default=0.16, type=float,
     #                   help='The drop radius of daisy drive . Default 0.16.')

    #parser.add_argument('-density_interaction_distance', '--DENSITY_INTERACTION_DISTANCE', default=0.02, type=float,
     #                   help='The density interaction distance of daisy drive. Default 0.02.')


    args_dict = vars(parser.parse_args())

    #11*50*100
    if args_dict.pop("print_header", None):
        print("result,size_all,generation,Drop_Size,SOMATIC_FITNESS_MUTLIPLIER_F,slice,drive1_freq,drive2_freq,drive3_freq,pop_size,drop_proportion")


    # The '-header' argument prints a header for the output. This can
    # help generate a nice CSV by adding this argument to the first SLiM run:
    # Next, assemble the command line arguments in the way we want to for SLiM:
    clargs = configure_slim_command_line(args_dict)

#initial dataframe: first run
    slim_result = run_slim(clargs)
    
    parsed_result = parse_slim(slim_result)

    print(parsed_result)


    #print(parsed_result)#
if __name__ == "__main__":
    main()



