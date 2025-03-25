#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CF2Insights - A tool for analysis and visualization of mutation data

This script extracts and visualizes mutation information from BAM files,
and generates region-organized HTML reports for easy viewing and sharing of results.
"""

import os
import sys
import argparse
import yaml
import pandas as pd
import concurrent.futures
import traceback
import base64
import jinja2
import logging
from collections import defaultdict, Counter
from typing import List, Tuple, Dict, Any
import time
import re

# Import in_sight modules
from ..vis_flow import base_visualization

# Default configuration file path
DEFAULT_CONFIG_PATH = 'cf2insights_config.yaml'

# Setup logging
def setup_logging(log_file="cf2insights.log"):
    """Set up logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger("CF2Insights")

logger = setup_logging()

def load_config(config_path):
    """Load configuration from YAML file"""
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        logger.info(f"Configuration loaded from {config_path}")
        return config
    except Exception as e:
        logger.error(f"Error loading configuration from {config_path}: {e}")
        logger.debug(traceback.format_exc())
        sys.exit(1)

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Analyze and visualize mutation data')
    
    # Configuration file option
    parser.add_argument('--config', '-c', default=DEFAULT_CONFIG_PATH,
                        help=f'Path to configuration file (default: {DEFAULT_CONFIG_PATH})')
    
    # Allow overriding config file with command line options
    parser.add_argument('--mutation_file', '-m',
                        help='Mutation information file (TSV format, cells as rows, mutations as columns)')
    parser.add_argument('--bam_path', '-b',
                        help='BAM file path')
    parser.add_argument('--reference_fasta', '-r',
                        help='Reference genome FASTA file path')
    parser.add_argument('--r_script',
                        help='R script path for plotting')
    parser.add_argument('--template_file',
                        help='HTML report template file path')
    parser.add_argument('--output_dir', '-o',
                        help='Base output directory path')
    parser.add_argument('--min_base_quality', type=int,
                        help='Minimum base quality value')
    parser.add_argument('--target_count', type=int,
                        help='Target read count')
    parser.add_argument('--max_workers', type=int,
                        help='Maximum number of worker threads for parallel processing')
    parser.add_argument('--align_samples', action='store_true',
                        help='Align samples')
    parser.add_argument('--mutation_parser', 
                        help='Rule for parsing mutation IDs (default, suffix, or regex:pattern)')
    
    return parser.parse_args()

def get_config():
    """Get configuration from file and command line arguments"""
    args = parse_args()
    
    # Load configuration from file
    if os.path.exists(args.config):
        config = load_config(args.config)
    else:
        # If no config file specified or doesn't exist, create empty config
        logger.warning(f"Configuration file {args.config} not found")
        config = {}
    
    # Override config with any command line arguments
    for key, value in vars(args).items():
        # Skip 'config' key and None values
        if key == 'config' or value is None:
            continue
        # For boolean flags, only override if explicitly provided on command line
        if key == 'align_samples' and not sys.argv.count('--align_samples'):
            continue
        config[key] = value
            
    # Check for required parameters
    required_params = ['mutation_file', 'bam_path', 'reference_fasta', 'r_script']
    missing_params = [param for param in required_params if param not in config]
    
    if missing_params:
        logger.error(f"Missing required parameters: {', '.join(missing_params)}")
        logger.error("Please provide them in the configuration file or as command line arguments")
        sys.exit(1)
        
    # Set defaults for optional parameters if not present
    if 'output_dir' not in config:
        config['output_dir'] = './output_base_visualization'
    if 'template_file' not in config:
        config['template_file'] = 'template_simple.html'
    if 'min_base_quality' not in config:
        config['min_base_quality'] = 10
    if 'target_count' not in config:
        config['target_count'] = 200
    if 'max_workers' not in config:
        config['max_workers'] = 10
    if 'align_samples' not in config:
        config['align_samples'] = True
    if 'mutation_parser' not in config:
        config['mutation_parser'] = 'default'
    
    return config

def parse_mutation_file(mut_file: str, parse_rule: str = "default") -> Dict[str, Any]:
    """
    Parse mutation file and return structured mutation information
    
    Args:
        mut_file: Path to mutation information file (TSV format)
        parse_rule: Rule for parsing mutation IDs. Options:
                    - "default": chromosome_position_ref_alt (e.g. chr1_39034563_T_A)
                    - "suffix": any_prefix_chromosome_position_ref_alt (e.g. P4_cSCC_chr1_1234_A_G)
                    - "regex:pattern": Use custom regex pattern with named groups
        
    Returns:
        Dict containing cell IDs and mutation information
        Format: {
            'cb_labels': [cell_id1, cell_id2,...],
            'mutations': {
                'chr22:22895504-22895504': {
                    'chrom': 'chr22',
                    'pos': 22895504,
                    'ref': 'C',
                    'alt': 'A',
                    'positive_cells': ['CAAGGCCAGACGACGT-1', ...]
                },
                ...
            }
        }
    """
    logger.info(f"Parsing mutation file: {mut_file} with rule: {parse_rule}")
    # Read TSV file
    df = pd.read_csv(mut_file, sep='\t', index_col=0)
    
    result = {
        'cb_labels': df.index.tolist(),
        'mutations': {}
    }
    
    # Parse each mutation site
    for mut_id in df.columns:
        # Apply the specified parsing rule
        if parse_rule == "default":
            # Parse chromosome_position_ref_alt format
            parts = mut_id.split('_')
            if len(parts) != 4:
                raise ValueError(f"Invalid mutation ID format: {mut_id}")
            
            chrom, pos, ref, alt = parts
            
        elif parse_rule == "suffix":
            # Parse any_prefix_chromosome_position_ref_alt format
            # For example: P4_cSCC_chr1_1234_A_G -> last 4 parts are what we need
            parts = mut_id.split('_')
            if len(parts) < 4:
                raise ValueError(f"Invalid mutation ID format: {mut_id}")
            
            # Extract the last 4 parts
            chrom, pos, ref, alt = parts[-4], parts[-3], parts[-2], parts[-1]
            
            # Verify the chromosome format (should start with 'chr')
            if not chrom.startswith('chr'):
                raise ValueError(f"Invalid chromosome format in mutation ID: {mut_id}")
            
        elif parse_rule.startswith("regex:"):
            # Use custom regex pattern
            pattern = parse_rule[6:]  # Remove 'regex:' prefix
            match = re.match(pattern, mut_id)
            if not match:
                raise ValueError(f"Mutation ID {mut_id} doesn't match pattern {pattern}")
            
            try:
                chrom = match.group('chrom')
                pos = match.group('pos')
                ref = match.group('ref')
                alt = match.group('alt')
            except IndexError:
                raise ValueError(f"Regex pattern must contain named groups: chrom, pos, ref, alt")
            
        else:
            raise ValueError(f"Unsupported parse_rule: {parse_rule}")
        
        region = f"{chrom}:{pos}-{pos}"  # Single base position
        
        # Record mutation information
        result['mutations'][region] = {
            'chrom': chrom,
            'pos': int(pos),
            'ref': ref,
            'alt': alt,
            'positive_cells': df.index[df[mut_id] == 1].tolist()
        }
    
    return result

# # 默认规则处理 chr1_39034563_T_A 格式
# result = parse_mutation_file("mutations.tsv")

# # 处理前缀格式如 P4_cSCC_chr1_1234_A_G
# result = parse_mutation_file("mutations.tsv", parse_rule="suffix")

# # 使用自定义正则表达式处理复杂格式
# pattern = r"(?P<sample>\w+)_(?P<type>\w+)_(?P<chrom>chr\w+)_(?P<pos>\d+)_(?P<ref>[ACGTN])_(?P<alt>[ACGTN])"
# result = parse_mutation_file("mutations.tsv", parse_rule=f"regex:{pattern}")

def generate_report_by_region(region_data: Dict, images_dir: str, reports_dir: str, 
                             region_id: str, template_file: str) -> None:
    """
    Generate an independent HTML report for each region using structured directory images
    
    Args:
        region_data: Dictionary containing region processing results
        images_dir: Directory path for this region's images
        reports_dir: Report output directory
        region_id: Region's unique identifier
        template_file: HTML template file path
    """
    # Get basic region information
    first_item = region_data['items'][0]
    
    # Prepare template data
    template_data = {
        'region_id': region_id,
        'chromosome': first_item['chrom'],
        'position': first_item['pos'],
        'reference': first_item['ref'],
        'alternate': first_item['alt'],
        'samples': []
    }
    
    # Add data for each cell
    for item in region_data['items']:
        plot_path = os.path.join(images_dir, f"{item['cb_label']}.png")
        
        # Read image and convert to base64
        image_base64 = ''
        if os.path.exists(plot_path):
            with open(plot_path, 'rb') as img_file:
                image_bytes = img_file.read()
                image_base64 = base64.b64encode(image_bytes).decode('utf-8')
        
        # 这里需要添加sampler_info
        # 假设当前item中有结果数据，或者从另一个地方获取
        # 如果没有可以创建一个默认值
        sampler_info = item.get('sampler_info', {'total': 0, 'sampled': 0})
        
        sample_data = {
            'cell_barcode': item['cb_label'],
            'image_base64': image_base64,
            'is_empty': not os.path.exists(plot_path),
            'note_message': 'Plot not generated' if not os.path.exists(plot_path) else '',
            'sampler_info': sampler_info  # 添加采样器信息
        }
        template_data['samples'].append(sample_data)
    
    # Generate report
    try:
        with open(template_file, 'r') as f:
            template_content = f.read()
        
        template = jinja2.Template(template_content)
        rendered_content = template.render(**template_data)
        
        # Write HTML file to reports directory
        output_html = os.path.join(reports_dir, f"{region_id}_report.html")
        with open(output_html, 'w', encoding='utf-8') as f:
            f.write(rendered_content)
            
        logger.info(f"Generated HTML report for region {region_id}")
        
    except Exception as e:
        logger.error(f"Error generating report for region {region_id}: {str(e)}")
        logger.debug(traceback.format_exc())

def process_visualization_tasks(config, processing_list):
    """
    Process all visualization tasks in parallel, without grouping by region
    
    Args:
        config: Configuration dictionary
        processing_list: List of tasks to process
    """
    # Create layered output directory structure
    reports_dir = os.path.join(config['output_dir'], 'reports', 'region_reports')
    images_base_dir = os.path.join(config['output_dir'], 'images', 'by_region')
    
    # Create necessary directories
    os.makedirs(reports_dir, exist_ok=True)
    os.makedirs(images_base_dir, exist_ok=True)
    
    # Pre-create region directories to avoid conflicts during concurrent creation
    all_regions = set()
    for item in processing_list:
        region_id = f"{item['chrom']}_{item['pos']}_{item['ref']}_{item['alt']}"
        all_regions.add(region_id)
    
    for region_id in all_regions:
        region_images_dir = os.path.join(images_base_dir, region_id)
        os.makedirs(region_images_dir, exist_ok=True)
    
    # Create list of all tasks
    all_tasks = []
    
    # Analyze all processing tasks for parallel execution
    for item in processing_list:
        region_id = f"{item['chrom']}_{item['pos']}_{item['ref']}_{item['alt']}"
        region_images_dir = os.path.join(images_base_dir, region_id)
        
        # Expected output file path
        output_filename = f"{item['cb_label']}.png"
        expected_output = os.path.join(region_images_dir, output_filename)
        
        # If file doesn't exist, add to task list
        if not os.path.exists(expected_output):
            task_args = {
                'bam_path': config['bam_path'],
                'region': item['region'],
                'reference_fasta': config['reference_fasta'],
                'ref_base': item['ref'],
                'r_script_path': config['r_script'],
                'output_dir': region_images_dir,
                'prefix': item['cb_label'],
                'run_visualization': True,
                'min_base_quality': config['min_base_quality'],
                'target_count': config['target_count'],
                'filter_tags': {'CB': item['cb_label']},
                'align_samples': config['align_samples']
            }
            print(f"task_args: {task_args}")
            all_tasks.append((task_args, item))
            logger.info(f"Prepared task: {item['cb_label']} at {item['region']}")
        else:
            logger.info(f"Skipped existing visualization: {item['cb_label']} at {item['region']}")
    
    logger.info(f"Total of {len(all_tasks)} tasks to execute with {config['max_workers']} workers")
    
    # Create process pool and submit all tasks
    with concurrent.futures.ProcessPoolExecutor(max_workers=config['max_workers']) as executor:
        # Create a dictionary to map futures to their task info
        future_to_task = {}
        
        # Submit all tasks
        for task_args, item in all_tasks:
            logger.info(f"Submitting task: {item['cb_label']} at {time.strftime('%H:%M:%S')}")
            future = executor.submit(base_visualization, **task_args)
            future_to_task[future] = item
        
        # 处理并存储结果
        results = {}
        
        # Process results as they complete
        for future in concurrent.futures.as_completed(future_to_task):
            item = future_to_task[future]
            try:
                result = future.result()
                # 保存采样器信息
                region_id = f"{item['chrom']}_{item['pos']}_{item['ref']}_{item['alt']}"
                if region_id not in results:
                    results[region_id] = {}
                
                # 保存该任务的sampler_info
                results[region_id][item['cb_label']] = {
                    'sampler_info': result.get('sampler_info', {'total': 0, 'sampled': 0})
                }
                
                logger.info(f"Completed task: {item['cb_label']} at {item['region']} at {time.strftime('%H:%M:%S')}")
            except Exception as e:
                logger.error(f"Visualization task error: {item['cb_label']} at {item['region']}, error: {e}")
                logger.debug(traceback.format_exc())
        
        # 将结果传递给报告生成函数
        generate_all_reports(processing_list, images_base_dir, reports_dir, config['template_file'], results)

def generate_all_reports(processing_list, images_base_dir, reports_dir, template_file, results=None):
    """Generate reports for all regions"""
    # Group by region
    region_groups = defaultdict(lambda: {'items': []})
    for item in processing_list:
        region_id = f"{item['chrom']}_{item['pos']}_{item['ref']}_{item['alt']}"
        
        # 克隆item以避免修改原始数据
        item_copy = item.copy()
        
        # 如果有结果数据，添加到item
        if results and region_id in results and item['cb_label'] in results[region_id]:
            item_copy.update(results[region_id][item['cb_label']])
            
        region_groups[region_id]['items'].append(item_copy)
    
    logger.info(f"Starting to generate reports for {len(region_groups)} regions")
    
    # Generate report for each region
    for region_id, group_data in region_groups.items():
        region_images_dir = os.path.join(images_base_dir, region_id)
        generate_report_by_region(group_data, region_images_dir, reports_dir, region_id, template_file)

def print_statistics(processing_list, config):
    """Print processing data statistics"""
    # Initialize statistics containers
    mutation_counter = defaultdict(int)  # Count processing times for each mutation site
    cell_counter = Counter()             # Count processing times for each cell
    unique_mutations = set()             # Unique mutation site set
    unique_cells = set()                 # Unique cell set

    # Traverse processing list for statistics
    for item in processing_list:
        # Count mutation sites (using chromosome+position+ref+alt as unique identifier)
        mutation_id = f"{item['chrom']}_{item['pos']}_{item['ref']}_{item['alt']}"
        mutation_counter[mutation_id] += 1
        unique_mutations.add(mutation_id)
        
        # Count cells
        cell_counter[item['cb_label']] += 1
        unique_cells.add(item['cb_label'])

    # Output statistics
    logger.info(f"Total task count: {len(processing_list):,}")
    logger.info(f"Unique mutation sites: {len(unique_mutations):,}")
    logger.info(f"Involved cells: {len(unique_cells):,}\n")

    logger.info("Mutation processing count TOP5:")
    for mut, count in Counter(mutation_counter).most_common(5):
        logger.info(f"· {mut}: {count} times")

    logger.info("\nCell processing task count TOP5:")
    for cell, count in cell_counter.most_common(5):
        logger.info(f"· {cell}: {count} times")
    
    # Save statistics to CSV
    stats_dir = os.path.join(config['output_dir'], 'data', 'stats')
    os.makedirs(stats_dir, exist_ok=True)
    
    # Save mutation statistics
    mut_stats = pd.DataFrame(list(mutation_counter.items()), 
                            columns=['mutation_id', 'count'])
    mut_stats.to_csv(os.path.join(stats_dir, 'mutation_stats.csv'), index=False)
    
    # Save cell statistics
    cell_stats = pd.DataFrame(list(cell_counter.items()), 
                             columns=['cell_barcode', 'count'])
    cell_stats.to_csv(os.path.join(stats_dir, 'cell_stats.csv'), index=False)

def main():
    """Main function"""
    # Get configuration
    config = get_config()
    
    # debug 
    print(f"config: {config}")
    
    # Create base output directory
    os.makedirs(config['output_dir'], exist_ok=True)
    
    # Parse mutation file
    mutation_data = parse_mutation_file(config['mutation_file'], config['mutation_parser'])
    
    # Create data directory and save raw input
    data_dir = os.path.join(config['output_dir'], 'data', 'raw')
    os.makedirs(data_dir, exist_ok=True)
    
    # Get all region and cell combinations to process
    processing_list = []
    for region, info in mutation_data['mutations'].items():
        for cb_label in info['positive_cells']:
            processing_list.append({
                'region': region,
                'cb_label': cb_label,
                'chrom': info['chrom'],
                'pos': info['pos'],
                'ref': info['ref'],
                'alt': info['alt']
            })
    
    # Print statistics
    print_statistics(processing_list, config)
    
    # Process data - use new parallel processing function
    process_visualization_tasks(config, processing_list)
    
    logger.info("Processing complete!")

if __name__ == "__main__":
    # Only execute when run as a script
    main()

