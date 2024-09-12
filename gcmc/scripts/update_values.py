import re
import os

def cpp_type_to_numpy(cpp_type):
    type_mapping = {
        'int': 'np.int32',
        'float': 'np.float32',
        'double': 'np.float64',
        'char': 'np.dtype("S1")',
        'unsigned int': 'np.uint32',
        'uint': 'np.uint32',
        'long': 'np.int64',
        'unsigned long': 'np.uint64',
        'short': 'np.int16',
        'unsigned short': 'np.uint16',
    }
    return type_mapping.get(cpp_type, cpp_type)

def parse_cpp_header(header_file):
    with open(header_file, 'r') as f:
        content = f.read()

    structs = {}
    constants = {}

    # Parse structs
    struct_pattern = r'struct\s+(\w+)\s*\{([\s\S]*?)\};'
    for match in re.finditer(struct_pattern, content):
        struct_name = match.group(1)
        struct_content = match.group(2)
        fields = []
        for line in struct_content.split('\n'):
            line = line.strip()
            if line and not line.startswith('//'):
                field_match = re.match(r'(\w+(?:\s+\w+)?)\s+(\w+)(?:\[(\d+)\])?;', line)
                if field_match:
                    field_type, field_name, array_size = field_match.groups()
                    numpy_type = cpp_type_to_numpy(field_type)
                    if array_size:
                        if field_type == 'char':
                            fields.append((field_name, f"np.dtype('S{array_size}')"))
                        else:
                            fields.append((field_name, numpy_type, int(array_size)))
                    else:
                        fields.append((field_name, numpy_type))
        structs[struct_name] = fields

    # Parse constants, including both #define and const declarations
    constant_patterns = [
        r'#define\s+(\w+)(?:\s+(.+))?',  # For #define, content is now optional
        r'const\s+\w+\s+(\w+)\s*=\s*(.+?);'  # For const declarations
    ]
    
    for pattern in constant_patterns:
        for match in re.finditer(pattern, content):
            constant_name = match.group(1)
            constant_value = match.group(2)
            
            if constant_value is None or constant_value.strip() == '':
                # Skip empty macros or macros without values
                continue
            
            constant_value = constant_value.strip()
            
            # Remove trailing comments (if any)
            constant_value = re.split(r'\s*//.*', constant_value)[0].strip()
            
            # Skip if the constant value is just another macro or complex preprocessor directive
            if constant_value.startswith('#') or constant_value.startswith('('):
                continue
            
            # Add the constant
            constants[constant_name] = constant_value

    return structs, constants

def generate_python_code(structs, constants):
    code = "import numpy as np\n\n"

    # Generate GcmcConstants class
    code += "class GcmcConstants:\n"
    for name, value in constants.items():
        code += f"    {name} = {value}\n"
    code += "\n"

    # Generate classes for struct types
    for struct_name, fields in structs.items():
        code += f"{struct_name}_dtype = np.dtype([\n"
        for field in fields:
            if len(field) == 2:
                field_name, field_type = field
                if field_type.startswith('np.dtype'):
                    code += f"    ('{field_name}', {field_type}),\n"
                else:
                    code += f"    ('{field_name}', {field_type}),\n"
            elif len(field) == 3:
                field_name, field_type, array_size = field
                if field_type in structs:
                    code += f"    ('{field_name}', {field_type}_dtype, ({array_size},)),\n"
                else:
                    code += f"    ('{field_name}', {field_type}, ({array_size},)),\n"
        code += "])\n\n"

    return code

def update_values_py(header_file, values_file):
    structs, constants = parse_cpp_header(header_file)
    generated_code = generate_python_code(structs, constants)

    # Create new content to append
    new_content = (
        '\n\n"""Automatically generated values from gcmc.h"""\n\n'
        "# Auto-generated from gcmc.h\n" +
        generated_code
    )

    # Append the new content to the file
    with open(values_file, 'a') as f:
        f.write(new_content)

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    header_file = os.path.join(script_dir, '..', 'cpp', 'gcmc.h')
    values_file = os.path.join(script_dir, '..', 'python', 'gcmcH.py')
    update_values_py(header_file, values_file)
    print(f"Created {values_file}")
