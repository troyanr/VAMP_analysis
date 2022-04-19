
import os
import re
import yamale

_FILE_DIR = os.path.dirname(os.path.abspath(__file__))
SCHEMA_NAME = '{}/params_schema.yaml'.format(_FILE_DIR)


def read_params(fname, schema_name=SCHEMA_NAME):
    '''
    Read parameters from params YAML file.

    Parameters
    ----------
    fname: str
        Path to params.yaml file
    schame_name: str
        Path to params.yaml schema.
        Default is python/params_schema.yaml

    Returns
    -------
    params: dict
        Dictionary of params values
    '''
    
    # read files
    schema = yamale.make_schema(schema_name)
    data = yamale.make_data(fname)
    
    # validate params with schema
    try:
        yamale.validate(schema, data)
    except yamale.YamaleError as e:
        print('Params validation failed!')
        print(e)
        return None
    data = data[0][0]

    # expand file regexes if necissary
    data_files = os.listdir(data['path_to_data_dir'])
    def expand_matches(file_pattern):
        pattern = re.sub(r'\\(.)', r'\1', file_pattern)
        pattern = re.compile(pattern)
        matches = [f for f in data_files if pattern.search(f)]
        return matches

    for input_output, sample in data['file_names'].items():
        if input_output == 'input':
            if isinstance(sample, str):
                matches = expand_matches(sample)
                data['file_names'][input_output] = matches
        else:
            for sample_name, replicate in sample.items():
                if isinstance(replicate, str):
                    matches = expand_matches(replicate)
                    data['file_names'][input_output][sample_name] = matches

    # fix sequence logo colors (convert lists to tuple)
    for feature, color in data['logo_features'].items():
        if isinstance(color, list):
            data['logo_features'][feature] = tuple(color)

    return data


if __name__ == '__main__':
    assert(read_params('{}/../params.yaml'.format(_FILE_DIR)))

