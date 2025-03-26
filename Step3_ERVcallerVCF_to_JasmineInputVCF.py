#!/usr/bin/env python3
"""
Convert ERVcaller VCF output to Jasmine-compatible VCF format.
"""

import re
import argparse
import os

def parse_infor(infor_field):
    """Parse the INFOR field to extract relevant information."""
    # INFOR=NAME,START,END,LEN,DIRECTION,STATUS
    parts = infor_field.split(';')
    for part in parts:
        if part.startswith('INFOR='):
            infor_data = part[6:].split(',')
            if len(infor_data) >= 4:  # Ensure we have enough data
                variant_len = infor_data[3]
                direction = infor_data[4] if len(infor_data) > 4 else None
                return variant_len, direction
    return None, None

def convert_vcf_for_jasmine(input_file, output_file):
    """
    Convert ERVcaller VCF to Jasmine-compatible format.
    
    Adds:
    1. Unique variant ID (CHROM_POS_ALT)
    2. SVLEN in INFO field
    3. Ensures SVTYPE is present
    4. Adds STRANDS info based on direction in INFOR
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            line = line.strip()
            
            # Write header lines unchanged
            if line.startswith('#'):
                # Add SVLEN and STRANDS INFO field definitions if they don't exist
                if line.startswith('##INFO') and 'SVLEN' not in line and line == list(filter(lambda x: x.startswith('##INFO'), open(input_file))).pop():
                    outfile.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">\n')
                    outfile.write('##INFO=<ID=STRANDS,Number=1,Type=String,Description="Strand orientation of the adjacency">\n')
                outfile.write(line + '\n')
                continue
            
            # Process data lines
            fields = line.split('\t')
            if len(fields) < 8:
                print(f"Warning: Skipping malformed line: {line}")
                continue
            
            chrom, pos, id_field, ref, alt, qual, filter_field, info = fields[:8]
            format_field = fields[8] if len(fields) > 8 else ''
            sample_data = fields[9:] if len(fields) > 9 else []
            
            # Extract variant type from ALT field
            svtype = None
            if alt.startswith('<'):
                # Extract the variant type from inside the angle brackets
                match = re.search(r'<([^>]+)>', alt)
                if match:
                    svtype_full = match.group(1)
                    # Further extract just the type (after colon if present)
                    if ':' in svtype_full:
                        svtype = svtype_full.split(':')[-1]
                    else:
                        svtype = svtype_full
            
            # Create unique ID combining chromosome and alternative
            unique_id = f"{chrom}_{pos}_{alt.replace('<', '').replace('>', '')}"
            
            # Parse INFO field to get variant length and direction
            variant_len, direction = parse_infor(info)
            
            # Prepare new INFO field
            info_parts = info.split(';')
            new_info_parts = []
            
            # Add SVTYPE if not already in INFO
            svtype_present = False
            for part in info_parts:
                if part.startswith('SVTYPE='):
                    svtype_present = True
                new_info_parts.append(part)
            
            if not svtype_present and svtype:
                new_info_parts.append(f"SVTYPE={svtype}")
            
            # Add SVLEN if we have it
            if variant_len and variant_len != 'NULL':
                new_info_parts.append(f"SVLEN={variant_len}")
            
            # Add STRANDS based on direction
            if direction:
                if direction == '+':
                    strands = '+-'
                elif direction == '-':
                    strands = '-+'
                else:
                    strands = '??' # Unknown
                new_info_parts.append(f"STRANDS={strands}")
            
            # Create new INFO field
            new_info = ';'.join(new_info_parts)
            
            # Write the modified line
            new_fields = [chrom, pos, unique_id, ref, alt, qual, filter_field, new_info]
            if format_field:
                new_fields.append(format_field)
            new_fields.extend(sample_data)
            
            outfile.write('\t'.join(new_fields) + '\n')

def process_multiple_files(file_list, output_dir):
    """
    Process multiple VCF files and create a file list for Jasmine.
    
    Args:
        file_list: List of input VCF files
        output_dir: Directory to save converted files
    
    Returns:
        Path to the created file list
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    jasmine_file_list = os.path.join(output_dir, "jasmine_file_list.txt")
    with open(jasmine_file_list, 'w') as list_file:
        for input_file in file_list:
            filename = os.path.basename(input_file)
            output_file = os.path.join(output_dir, f"jasmine_{filename}")
            
            convert_vcf_for_jasmine(input_file, output_file)
            list_file.write(f"{output_file}\n")
    
    return jasmine_file_list

def main():
    parser = argparse.ArgumentParser(description='Convert ERVcaller VCF to Jasmine-compatible format')
    parser.add_argument('-i', '--input', required=True, help='Input VCF file or directory of VCF files')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory for converted files')
    parser.add_argument('-j', '--jasmine_output', help='Output file for Jasmine merged results')
    
    args = parser.parse_args()
    
    # Determine if input is a single file or directory
    if os.path.isdir(args.input):
        vcf_files = [os.path.join(args.input, f) for f in os.listdir(args.input) if f.endswith('.vcf')]
        if not vcf_files:
            print(f"No VCF files found in {args.input}")
            return
    else:
        vcf_files = [args.input]
    
    # Process the files and create a file list for Jasmine
    file_list_path = process_multiple_files(vcf_files, args.output_dir)
    
    print(f"Conversion complete. Jasmine file list created at: {file_list_path}")
    print("\nTo run Jasmine, use the following command:")
    jasmine_output = args.jasmine_output or os.path.join(args.output_dir, "jasmine_merged.vcf")
    print(f"jasmine file_list={file_list_path} out_file={jasmine_output}")

if __name__ == "__main__":
    main()
