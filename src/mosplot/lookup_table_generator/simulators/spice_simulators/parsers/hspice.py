"""
Author: Raphael Gonzalez (RAffA), Mathew Spencer
github: RaffaGonzo


Copyright 2021 Raphael Gonzalez, Mathew Spencer

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
import struct
from os.path import basename, split, join
from os import mkdir
from shutil import rmtree
import pickle
import sys

COUNT = 0

def read_data_block(data_block, vrsn):
    data_format = {'9601': ('<f', 4),  # float
                   '2001': ('<d', 8)}  # double precision
    assert vrsn in data_format.keys(), 'only the 2001 and 9601 formats are supported'
    data_bytes = (data_block[0 + i: data_format[vrsn][1] + i] for i in range(0, len(data_block), data_format[vrsn][1]))  # break the block into chunks of 4 (or 8) bytes each depending on the encoding..?
    block_lst = []
    for byte_set in data_bytes:
        val = struct.unpack(data_format[vrsn][0], byte_set)[0]
        if val == 1.0000000150474662e+30 or val == 1e+30:  # I think there may be a precision error | 32 bit doesn't have the same level of precision as 64 bit... but neither of these is large for float
            global COUNT
            COUNT += 1
        block_lst.append(val)
    return block_lst

def read_binary_signal_file(file_path):
    """
    reads the binary file (Hspice output for POST = 1)
    :param file_path: path to file, taken as a string
    :return: file_header_block which contains information about the test-bench settings, data_list which is a 2D list containing the floating point
    """
    with open(file_path, 'rb') as file:
        decode_info = file.read(16)  # this is the block header, first 4 bytes is an endian indicator, next 4 is a blank, then another endian indicator, then last 4 are the block size
        block_size = decode_info[-4:]
        block_size_int = struct.unpack('<i', block_size)[0]  # the struct will return a tuple no matter what, since we know there is one element in the tuple we can extract it on the spot
        file_header_block = file.read(block_size_int).decode('UTF-8')
        version = file_header_block[20:24]
        version = version.strip(' ')
        if not version:
            version = file_header_block[16:20]
        footer = file.read(4)
        data_list = []
        while True:
            block_header = file.read(16)
            if not block_header:  # this allows for the assertion to not get invoked when it reaches the end of the file
                break
            block_size = block_header[-4:]
            block_size_int = struct.unpack('<i', block_size)[0]
            data = file.read(block_size_int)
            data_out = read_data_block(data, version)
            data_list.append(data_out)
            footer = file.read(4)

        assert struct.unpack('<i', footer)[0] == block_size_int, 'footer size: ' + str(struct.unpack('<i', footer)[0]) + ' header size: ' + str(block_size_int)  # a check against any errors or corruption

        return file_header_block, data_list

def get_var_name_idx(var_lst):
    for i, item in enumerate(var_lst):
        try:
            int(item)
        except ValueError:
            return i
    return None

def parse_var_name(name):
    name_parts = name.split('(')
    new_name_parts = []
    for nm in name_parts:
        nm = nm.replace('.', '_')
        nm = nm.replace(':', '_')
        new_name_parts.append(nm)
    try:
        if new_name_parts[0] == new_name_parts[1]:
            return new_name_parts[0]
        return new_name_parts[0] + '_' + new_name_parts[1]
    except IndexError:
        return new_name_parts[0]

def parse_header(header_str):
    header_lst = header_str.split(' ')
    header_lst = list(filter(lambda x: x != '', header_lst))
    header_lst = header_lst[:2] + header_lst[15:]  # this chops out the copy right stuff
    analysis_type = header_lst[2]
    var_info = header_lst[3:-1]
    first_var_units = {
        '1': 'sec',
        '2': 'freq',
        '3': 'volt'
    }
    gnrl_var_units = {
        '1': 'volt',  # 1	voltage as dependent var in sweep or transient
        '2': 'volt',  # 2	voltage as dependent var in AC analysis
        '8': 'Amp',  # 8	current
        '9': 'Amp',  # 9	magnitude of current (not verified)
        '10': 'Amp',  # 10	real part of current (not verified)
        '11': 'iAmp',  # 11	imaginary part of current (not verified)
        '15': 'Amp',  # 15	current
    }
    # from https://github.com/krichter722/gwave/blob/master/doc/hspice-output.txt
    # First comes an integer for each variable, indicating the variable type.
    #
    # The first variable is the independent variable.  The integral types for this
    # also indicates the analysis type.
    # 1	Transient analysis; type is time.
    # 2	AC analysis; type is frequency
    # 3	DC sweep; type is voltage
    #
    # The type numbers for the dependent variables are:
    #
    # 1	voltage as dependent var in sweep or transient
    # 2	voltage as dependent var in AC analysis
    # 8	current
    # 9	magnitude of current (not verified)
    # 10	real part of current (not verified)
    # 11	imaginary part of current (not verified)
    # 15	current
    # 16 the units stuff is not used but perhaps it would be useful to other people.
    var_idx = get_var_name_idx(var_info)
    var_nums, var_names = var_info[:var_idx], var_info[var_idx:]  # these two variables should have identical length
    var_names = [parse_var_name(name) for name in var_names]
    return len(var_names), var_names

def dict_to_matlab_str(data_dict):
    rtn_str = ''
    for key, arr in data_dict.items():
        rtn_str += key + '= {'
        arr_str = str(arr)
        arr_str = arr_str[1:-1]
        arr_str = arr_str.strip(',')
        arr_str = arr_str.replace('],', '];\n')
        rtn_str += arr_str + '};'
    return rtn_str

def dict_to_csv_str(data_dict):
    key_lst = list(data_dict.keys())
    csv_str = ''
    for key in key_lst: csv_str += key + ', '
    csv_str = csv_str[:-2] + '\n'  # init the csv str

    for i in range(len(data_dict[key_lst[0]])):  # essentially for each sweep
        for j in range(len(data_dict[key_lst[0]][i])):  # for each item index inside the 2d arr
            for key in key_lst:  # select the right column in a defined order
                try:
                    csv_str += str(data_dict[key][i][j]) + ', '
                except IndexError:  # only should happen when there is a sweep
                    try:
                        csv_str += ' , '
                    except NameError:
                        pass

            csv_str = csv_str[:-2] + '\n'
    csv_str = csv_str[:-1]  # take out final return char
    return csv_str

def dict_to_mat_obj(data_dict):
    if sio in sys.modules and np in sys.modules:
        for var, sweep_lst in data_dict.items():
            data_dict[var] = np.array(sweep_lst)
        return data_dict

def break_by_sweep(data_lst, vrsn):
    """
    breaks all the data pulled from .tr* file into sweeps, if there are none it will just returns a 1 x n 2d list
    :param data_lst: a list with the values read from the binary file
    :param vrsn: the version of binary file being parsed
    :return: 2d list of data separated by sweep
    """
    vrsn_dict = {
        '2001': 1e30,
        '9601': 1.0000000150474662e+30
    }
    if not data_lst:
        return []
    i = data_lst.index(vrsn_dict[vrsn])
    return [data_lst[:i]] + break_by_sweep(data_lst[i + 1:], vrsn)

def general_make_dict(var_lst, data_lst, header_str, MULTI_SWEEP_FLAG):
    var_dict = {var: [] for var in var_lst}
    if MULTI_SWEEP_FLAG:
        var_lst = var_lst[:-1]
        sweep_val_lst = []

    else:
        sweep_val_lst = []  # won't have anything added to it in this case though
    var_count = len(var_lst)

    for i, sweep in enumerate(data_lst):
        if MULTI_SWEEP_FLAG:
            if int(sweep[0]) == sweep[0]:
                sweep_val_lst.append(int(sweep[0]))
            else:
                sweep_val_lst.append(float(sweep[0]))
            sweep = sweep[1:]
        for var in var_dict.keys():
            var_dict[var].append([])
        for j, val in enumerate(sweep):
            var = var_lst[j % var_count]
            var_dict[var][i].append(val)
    header_lst = header_str.split(' ')
    header_lst = list(filter(lambda x: x != '', header_lst))
    header_lst = header_lst[17:]
    dict_key = ''
    for var in header_lst[:-1]:  # by the end of this loop if there is a sweep dict_key will have the sweep_var
        try:
            int(var)
        except ValueError:
            dict_key = var

    for i, lst in enumerate(var_dict[parse_var_name(dict_key)]): var_dict[parse_var_name(dict_key)][i] = lst

    if MULTI_SWEEP_FLAG:
        return var_dict, MULTI_SWEEP_FLAG, dict_key, sweep_val_lst
    return var_dict, MULTI_SWEEP_FLAG, None, None

def ac_make_dict(var_lst, data_lst, header_str, MULTI_SWEEP_FLAG):
    var_lst = [['HERTZ']] + [[f'{parse_var_name(var)}_Mag', f'{parse_var_name(var)}_Phase'] for var in var_lst if var != 'HERTZ']  # increase the # of vars and prep for Re and Im parts to each
    var_lst = sum(var_lst, [])  # flatten
    var_dict = {var: [[] for i in range(len(data_lst))] for var in var_lst}
    var_count = len(var_lst)
    for sweep_idx, lst in enumerate(data_lst):
        i = 0
        while i < len(lst):
            if i % var_count == 0:
                var_dict['HERTZ'][sweep_idx].append(lst[i])
                i += 1
            else:
                var_dict[var_lst[i % var_count]][sweep_idx].append(lst[i])
                var_dict[var_lst[(i % var_count) + 1]][sweep_idx].append(lst[i + 1])
                i += 2

    return var_dict, MULTI_SWEEP_FLAG, None, None  # the Nones make this function play nice with the way I used it later.

def write_to_dict(data_lst, header_str, ext):
    # this function puts it all together so that we get a dictionary with variable names for keys and lists of the respective values. If there are multiple sweeps then there will be multiple lists, the lists are always 2D
    var_count, var_lst = parse_header(header_str)
    # var_dict = {}
    version = header_str[20:24]
    version = version.strip(' ')
    if not version:
        version = header_str[16:20]

    data_lst = sum(data_lst, [])
    data_lst = break_by_sweep(data_lst, version)
    MULTI_SWEEP_FLAG = len(data_lst) > 1
    assert COUNT > 0  # for debugging

    func_dict = {'tr': general_make_dict,
                 'sw': general_make_dict,
                 'ac': ac_make_dict}

    return func_dict[ext](var_lst, data_lst, header_str, MULTI_SWEEP_FLAG)

def write_to_file(write_path, file_content, ext=None):
    if type(file_content) == str:
        with open(write_path, 'w+') as file:
            file.write(file_content)
    elif ext == 'mat':
        sio.savemat(write_path, file_content)
    else:
        with open(write_path, 'wb') as file:
            pickle.dump(file_content, file, protocol=pickle.HIGHEST_PROTOCOL)

def get_outfile_name(path, ext=None):
    file_name = basename(path)
    outfile_name = file_name.replace('.', '_')
    if ext is not None:
        outfile_name += f'.{ext}'
    return outfile_name

def import_export_binary(path, ext, from_ext):
    """
    this function will import *.tr*, *.sw*, and *.ac* binary and put a similarly named csv, matlab *.m, or python.pickle file/ folder of files in the same directory
    :param path: (str) path to the *.tr* binary to be imported
    :param ext: (str) the extension of the output file [only options *.m, *.csv, and *.pickle]
    :return: None
    """
    dict_to_pickle_dict = lambda dict_data: dict_data

    ext_dict = {'m': dict_to_matlab_str,
                'csv': dict_to_csv_str,
                'pickle': dict_to_pickle_dict,
                'mat': dict_to_mat_obj}
    try:
        header, data = read_binary_signal_file(path)
        main_dict, sweep_flag, sweep_var_name, sweep_values = write_to_dict(data, header, from_ext)
        # use sweep_flag to tell if it should be a folder of csvs or single file
        if sweep_flag and ext == 'csv':
            content = ''
            fpath = path + '_' + sweep_var_name + '('
            outfile_name = get_outfile_name(path, ext)
            folder_name = get_outfile_name(outfile_name)  # gets rid of dot
            og_path = str(split(fpath)[0])
            path = join(og_path, folder_name)
            try:
                mkdir(path)
            except FileExistsError:  # this will overwrite an existing folder with the same name
                rmtree(path)
                mkdir(path)

            for i in range(len(list(main_dict.values())[0])):
                temp_dict = {}
                for key in main_dict.keys():
                    temp_dict[key] = [main_dict[key][i]]
                file_name = get_outfile_name(fpath) + str(sweep_values[i]) + ').csv'
                full_path = join(path, file_name)
                content = dict_to_csv_str(temp_dict)
                write_to_file(full_path, content)
        else:
            outfile_name = get_outfile_name(path, ext)
            dir_path, _ = split(path)
            file_path = join(str(dir_path), outfile_name)
            content = ext_dict[ext](main_dict)
            write_to_file(file_path, content, ext)
        return main_dict, content

    except KeyError:
        raise ValueError('the extension must NOT have the "." in it ex: "csv" | the only extension options are m, csv, and pickle')
    except FileNotFoundError as err:
        print(err)

def import_export_ascii(path, ext):
    dict_to_pickle_dict = lambda dict_data: dict_data

    ext_dict = {'m': dict_to_matlab_str,
                'csv': dict_to_csv_str,
                'pickle': dict_to_pickle_dict,
                'mat': dict_to_mat_obj}
    try:
        main_dict = signal_file_ascii_read(path)
        content_obj = ext_dict[ext](main_dict)
        sweep_flag = len(list(main_dict.values())[0]) > 1
        sweep_var_name = ''
        sweep_values = []
        for key, val in main_dict.items():
            if len(val[0]) == 1 and sweep_flag:
                sweep_var_name = key
                sweep_values = sum(val, [])
                break
        if sweep_flag and ext == 'csv':
            fpath = f'{path}_{sweep_var_name}('
            outfile_name = get_outfile_name(path, ext)
            folder_name = get_outfile_name(outfile_name)  # gets rid of dot
            og_path = str(split(fpath)[0])
            path = join(og_path, folder_name)
            try:
                mkdir(path)
            except FileExistsError:  # this will overwrite an existing folder with the same name
                rmtree(path)
                mkdir(path)

            for i in range(len(list(main_dict.values())[0])):
                temp_dict = {}
                for key in main_dict.keys():
                    temp_dict[key] = [main_dict[key][i]]
                file_name = get_outfile_name(fpath) + str(sweep_values[i]) + ').csv'
                full_path = join(path, file_name)
                content = dict_to_csv_str(temp_dict)
                write_to_file(full_path, content)
        else:
            outfile_name = get_outfile_name(path, ext)
            dir_path, _ = split(path)
            file_path = join(str(dir_path), outfile_name)
            write_to_file(file_path, content_obj, ext)
        return main_dict, content_obj
    except KeyError:
        raise ValueError('the only extension options are "m," "csv," or "pickle"')
    except FileNotFoundError as err:
        print(err)

def usage():
    usage_str = """Usage:
    python Hspice_parse <input file path> <extension of output file(s)>

    *** note: This parser only supports the 2001, 9601 formats
              designated by ".option post_version=9601" for example

    supported input file types are:
        - AC simulations *.ac*
        - DC sweeps *.sw*
        - Transient simulations *.tr*

    supported output file types are:
        - CSV note that if there are multiple sweeps in the input
          file this option will output a folder of CSVs
        - pickle
        - matlab's *.m

    -h or --help to get this message
                """
    print(usage_str)

def get_from_ext(path):  # this function is for auto detecting the input file extention
    file_name = basename(path)
    name_components = file_name.split('.')
    entire_ext = name_components[-1]
    return entire_ext[:2]  # returns the first two chars of the extension i.e. "ac" "sw" "tr" and ignores the numbers after

def reformat(sweep_lst):
    """
    This function will change the original output of signal_file_ascii_read() to the format that
    the write to dict function makes.
    :param sweep_lst: list of dictionaries
    :return: dictionary of 2D lists
    """
    rtn_dict = {key: [] for key in sweep_lst[0].keys()}
    for sweep in sweep_lst:
        for key, value in sweep.items():
            rtn_dict[key].append(value)
    return rtn_dict

###################### ascii parsing ###################
# author: Mathew Spencer

# class simulation_result():
#     """All the results germane to a single hSpice simulation including
#     results, measures, parameters, temperature and alter number.
#
#     Contains dicts called 'results' 'measures' 'parameters' and
#     'misc' each of which is keyed with named values and has values
#     corresponding to the results.
#     """
#     pass


# class sweep_result():
#     """Aggregates the results of individual simulations
#     """
#     pass

def signal_file_ascii_read(post2file):
    """Accepts a raw sweep or measure file from spice in ASCII format
    (.tr*, .sw*, .mt*, etc.) and returns a list of dicts containing the
    signal and parameter names as keys.  Signal name keys have an array
    containing the simulation results as values.  Parameter name keys have
    the parameter as a value.  There is one dict in the array for each
    simulation sweep point."""

    ## Preamble part of header
    f = open(post2file)
    l = f.readline()  # 1st line contains string w/ # of variables
    nauto = int(l[0:4])
    nprobe = int(l[4:8])
    nsweepparam = int(l[8:12])
    l = f.readline()  # 2nd and 3rd lines are copyright and stuff
    l = f.readline()
    ndataset = int(l.split()[-1])  # 3rd line ends with the number of data sets

    ## Number and name part of header
    l = f.readline()  # 4th line+ useful, but may be wrapped
    while l.find('$&%#') == -1:
        l = l + f.readline()
    l = l.replace('\n', '')
    simparams = l.split()[:-1]  # Throw away terminator string
    datatypes = simparams[0:nauto + nprobe]  # NOTUSED
    varnames = simparams[nauto + nprobe:2 * (nauto + nprobe)]
    paramnames = simparams[2 * (nauto + nprobe):]

    # Transform varnames and paramnames for Matlab
    varnames = [x.partition('(')[0] if x.startswith('x') else x for x in varnames]
    varnames = [x.replace('(', '_') for x in varnames]
    varnames = [x.replace('.', '_') for x in varnames]
    varnames = [x.replace(':', '_') for x in varnames]
    paramnames = [x.replace(':', '_') for x in paramnames]
    paramnames = ['param_' + x for x in paramnames]

    # Read data block
    all_sweep_results = []
    l = f.readline().strip()
    while l:
        # Set up data storage structure
        this_sweep_result = {}
        for name in varnames:
            this_sweep_result[name] = []
        for name in paramnames:
            this_sweep_result[name] = 0

        ## Data block for one sweep
        numbers = []
        fieldwidth = l.find('E') + 4
        while True:
            lastrow = l.find('0.1000000E+31') != -1
            while l:
                numbers.append(float(l[0:fieldwidth]))
                l = l[fieldwidth:]
            if lastrow:
                break
            else:
                l = f.readline().strip()
        numbers = numbers[:-1]  # throw away terminator
        params = numbers[:nsweepparam]
        for index in range(len(varnames)):
            this_sweep_result[varnames[index]] = numbers[nsweepparam + index::nauto + nprobe]
        this_sweep_result.update(zip(paramnames, [[x] for x in params]))
        all_sweep_results.append(this_sweep_result)
        l = f.readline().strip()

    f.close()
    return reformat(all_sweep_results)

def signal_array_to_matlab_string(signal_array, accum=True):  # Todo: This should be obsolete now!
    """ Takes in a signal array produced by signal_file_read and returns a string that
    is suitable for inclusion in a .m file to be loaded into matlab.  Makes 2D arrays
    for DC and AC sims.  1D arrays numbered by sweep for transients b/c length uneven.
    """
    matlab_string = ''
    if accum:
        # first, slice across the signal array to  gather like signals into 2D arrays
        signal_accumulation = signal_array[0]
        for name, signal in signal_accumulation.items():
            signal_accumulation[name] = ['%.5e' % x for x in signal]
        for simulation_result in signal_array[1:]:
            for name, signal in simulation_result.items():
                signal_accumulation[name].extend([';'] + ['%.5e' % x for x in signal])
        # Second write out to a .m file string
        for name, signal in signal_accumulation.items():
            matlab_string = matlab_string + name + '= [' + ' '.join(signal) + '];\n'
    else:
        # some cases, esp transient sims, can't be accumulated b/c different sizes
        matlab_string = ''
        signal_accumulation = signal_array[0]
        for name, signal in signal_accumulation.items():
            signal_accumulation[name] = [['%.5e' % x for x in signal]]
        for simulation_result in signal_array[1:]:
            for name, signal in simulation_result.items():
                signal_accumulation[name].append(['%.5e' % x for x in signal])
        for name, signal_list in signal_accumulation.items():
            matlab_string = matlab_string + name + '={'
            for signal in signal_list:
                matlab_string = matlab_string + '[ ' + ' '.join(signal) + '];\n'
            matlab_string = matlab_string + '};'
    return matlab_string

def measure_file_read(measure_file):
    """Reads a measure file and resturns a dict keyed on the names of the columns
    with values that are a list of results"""
    f = open(measure_file)
    l = f.readline().strip()  # 1st line auto generated
    param_count = int(l.split('PARAM_COUNT=')[-1].strip())
    l = f.readline().strip()  # 2nd line is the comment
    l = f.readline().strip()  # 3rd line starts var names, last is alter#
    while (l.find('alter#') == -1):
        l = l + ' ' + f.readline().strip()
    l = l.replace('#', '')
    varnames = l.split()
    varnum = len(varnames)
    measure_result = {name: [] for name in varnames}
    l = f.readline().strip()  # 4th line starts data clumps to EOF

    while l:
        vals = l.split()
        if len(vals) != varnum:  # accumulate results over multiple lines
            l = l + ' ' + f.readline().strip()
        else:  # then write to measure dict
            for name, val in zip(varnames, vals):
                if val == 'failed':  # skip over failed measures
                    val = 0  # silent failures are the best kind!
                measure_result[name].append([float(val)])  # todo that I added [] to make the lists 2D not sure if necessary
            l = f.readline().strip()
    f.close()
    return measure_result

def measure_result_to_matlab_string(measure_result):  # this should be obsolete now
    matlab_string = ''
    for name, signal in measure_result.items():
        matlab_string = matlab_string + name + '= [' + ' '.join(['%.5e' % x for x in signal]) + '];\n'
    return matlab_string
# end  of Spencer's (lightly modified) code

def is_binary(file_path):
    with open(file_path) as file:
        try:
            file.readline()
            return False  # this means that it is ascii
        except UnicodeDecodeError:
            return True  # this means that it is binary

def import_export(path, ext):
    if is_binary(path):
        from_ext = get_from_ext(path)
        import_export_binary(path, ext, from_ext)
    else:
        import_export_ascii(path, ext)

if __name__ == '__main__':
    import sys
    try:
        if sys.argv[1] == '-h' or sys.argv[1] == '--help' or sys.argv[1] == '-help':
            usage()
        else:
            path = sys.argv[1]
            extension = sys.argv[2]
            if extension == 'mat':
                try:
                    import scipy.io as sio
                    import numpy as np
                except ImportError:
                    print("To use Matlab's .mat extension, you must have Scipy and Numpy installed on your machine. See https://www.scipy.org/install.html and https://numpy.org/install for details.")
                    exit(1)
            import_export(path, extension)
    except:
        usage()
        resp_str = '\n\n This command needs the path to the file for import and the extension of the file being output.\n The extension can be either m, csv, or pickle. No other extensions are supported with this file at this time. If you are seeing this message, then something went wrong.'
        print(resp_str)
