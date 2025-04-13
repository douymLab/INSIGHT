import os
from in_sight.reads_manager import create_combined_read_manager
from in_sight.utils import adjust_chromosome_name
from in_sight.add_attr_info import hsnp_processor, str_data_processor
from in_sight.analysis.hmm_model import StrReads
from in_sight.utils import get_flanking_sequences, determine_file_type, get_file_mode
from importlib import resources
import in_sight
import os
import subprocess
from typing import List, Tuple, Union, Dict
import time
import io
from PIL import Image
import cairosvg
from dataclasses import dataclass, field
from typing import Any
import pandas as pd
import psutil
import logging
import traceback
from typing import Optional
import pysam
@dataclass
class STRRegion:
    chrom: str
    start: int
    end: int
    str_unit: str
    str_id: str
    str_unit_length: int
    average_length: float
    additional_info: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        # Convert types if needed
        self.start = int(self.start)
        self.end = int(self.end)
        self.str_unit_length = int(self.str_unit_length)
        self.average_length = float(self.average_length)

    def __getitem__(self, key):
        if key in self.__dict__:
            return getattr(self, key)
        return self.additional_info.get(key)

    def __setitem__(self, key, value):
        if key in self.__dict__:
            setattr(self, key, value)
        else:
            self.additional_info[key] = value

    def to_dict(self) -> Dict[str, Any]:
        result = {k: v for k, v in self.__dict__.items() if k != 'additional_info' and v is not None}
        result.update(self.additional_info)
        return result

def create_str_region_from_tsv(tsv_path: str, str_id: str) -> STRRegion:
    """Create a STRRegion object from a TSV file"""
    df = pd.read_csv(tsv_path, sep='\t')
    str_region_info = df[df['str_id'] == str_id].to_dict(orient='records')[0]
    
    # Get the expected field names from STRRegion (excluding additional_info)
    expected_fields = [f for f in STRRegion.__dataclass_fields__.keys() if f != 'additional_info']
    
    # Split the dictionary into expected fields and additional fields
    main_fields = {k: v for k, v in str_region_info.items() if k in expected_fields}
    additional_fields = {k: v for k, v in str_region_info.items() if k not in expected_fields}
    
    # Create the STRRegion with additional_info containing extra fields
    return STRRegion(**main_fields, additional_info=additional_fields)

def parse_region(region_str: str) -> Tuple[str, int, int]:
    """Parse region string into chrom, start, end components"""
    if ':' not in region_str:
        raise ValueError(f"Invalid region format: {region_str}")
    
    chrom, pos_part = region_str.split(':', 1)
    pos_part = pos_part.replace(',', '')  # 移除可能存在的千位分隔符
    
    if '-' in pos_part:
        start_str, end_str = pos_part.split('-', 1)
    else:
        start_str = end_str = pos_part
    
    try:
        start = int(start_str)
        end = int(end_str)
    except ValueError:
        raise ValueError(f"Invalid position in region: {region_str}")
    
    # 新增处理：当start=end时自动将end+1
    if start == end:
        end += 1
    
    if start > end:
        raise ValueError(f"Start position ({start}) cannot be greater than end position ({end})")
    
    return chrom, start, end

def run_r_visualization(base_df_path, read_df_path, reference_df_path, output_file="output.svg", 
                       mut_start=None, mut_end=None, mut_type=None, start_axis=None,
                       r_script_path=resources.files(in_sight).joinpath('r_script/plot_plain_with_full_features.R')):
    """
    Call Rscript for data visualization.

    Args:
        base_df_path (str): Path to base dataframe CSV file.
        read_df_path (str): Path to read dataframe CSV file.
        reference_df_path (str): Path to reference dataframe CSV file.
        output_file (str): Path to output image file (PNG).
        mut_start (int, optional): Mutation start position.
        mut_end (int, optional): Mutation end position.
        mut_type (str, optional): Mutation type.
        start_axis (int, optional): Start position of the reference sequence.

    Returns:
        None
    """
    # Ensure Rscript exists
    if not os.path.exists(r_script_path):
        raise FileNotFoundError(f"R script not found at {r_script_path}")

    # Prepare command arguments
    cmd = ["Rscript", r_script_path, 
           base_df_path, 
           read_df_path, 
           reference_df_path,
           output_file]
    
    # Add optional arguments if provided
    if mut_start is not None:
        cmd.append(str(mut_start))
        if mut_end is not None:
            cmd.append(str(mut_end))
            if mut_type is not None:
                cmd.append(str(mut_type))
                if start_axis is not None:
                    cmd.append(str(start_axis))

    # Call Rscript
    try:
        print(cmd)
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Rscript execution failed: {e}")

    print(f"Visualization saved to: {output_file}")

def str_visulization(bam_path: str,
                    region: Union[str, List[str], List[Tuple[str, str]]],
                    reference_fasta: str = None,
                    str_id: str = None,
                    str_region_info_path: str = None,
                    r_script_path: str = resources.files(in_sight).joinpath('r_script/plot_base.R'),
                    output_dir: str = "./",
                    prefix: str = "pileup_analysis",
                    run_visualization: bool = True,
                    mutation_region_label: str = None,  # 指定哪个区域标签用于突变范围
                    hsnp_region_label: str = None,     # 指定哪个区域应用hsnp处理器
                    hsnp_info: Dict[str, Dict[int, Tuple[str, str]]] = None,  # 添加hsnp信息参数
                    align_samples: bool = False,
                    start_axis: int = 1,
                    mut_type: str = "STR",
                    **kwargs):
    """
    Perform single base region analysis and generate visualizations
    
    Parameters:
        bam_path: Path to BAM/CRAM file
        region: Analysis region(s) in one of these formats:
               - Single string: "chrom:start-end" or "chrom:pos"
               - List of strings: ["chrom:start-end", "chrom:pos", ...]
               - List of tuples: [(region_string, label), ...] where each tuple contains a region string and a custom label
        reference_fasta: Path to reference genome
        str_id: STR ID
        str_region_info_path: Path to STR region TSV file
        r_script_path: Path to R script for visualization
        output_dir: Output directory (default current directory)
        prefix: Output file prefix (default 'pileup_analysis')
        run_visualization: Whether to execute R visualization (default True)
        mutation_region_label: Label of the region to use for mutation start/end coordinates (default None, which uses all regions)
        hsnp_region_label: Label of the region to apply hsnp processor (default None)
        hsnp_info: HSNP information dictionary for the hsnp processor (required if hsnp_region_label is specified)
        align_samples: Whether to align samples (default False)
        start_axis: Start position of the reference sequence (default 1)
        **kwargs: Additional parameters passed to sampler or read manager
    Returns:
        Dictionary containing generated DataFrames and file paths
    """
    # 处理不同类型的区域输入
    if isinstance(region, str):
        # 单个区域字符串
        regions_list = [(region, "analysis_region")]
    elif isinstance(region, list):
        if not region:
            raise ValueError("Region list cannot be empty")
        
        # 检查是否包含自定义标签（每个元素是元组）
        if isinstance(region[0], tuple):
            # 列表中的元组: [(region_string, label), ...]
            regions_list = region
        else:
            # 列表中的字符串: [region1, region2, ...]
            regions_list = [(reg, f"region_{i+1}") for i, reg in enumerate(region)]
    else:
        raise ValueError("Region must be a string, list of strings, or list of (string, label) tuples")
    
    # 解析所有区域
    parsed_regions = []
    chrom = None
    
    for reg_info in regions_list:
        reg_str, label = reg_info
        curr_chrom, start, end = parse_region(reg_str)
        
        # 确保所有区域都在同一条染色体上
        if chrom is None:
            chrom = curr_chrom
        elif chrom != curr_chrom:
            raise ValueError(f"All regions must be on the same chromosome. Found {chrom} and {curr_chrom}")
        
        # 添加到解析区域列表，使用提供的标签
        parsed_regions.append((start, end, label))
    
    # 检查hsnp相关参数
    if hsnp_region_label is not None:
        if hsnp_info is None:
            raise ValueError("hsnp_info must be provided when hsnp_region_label is specified")
            
        # 验证指定的hsnp区域标签存在
        if not any(label == hsnp_region_label for _, _, label in parsed_regions):
            raise ValueError(f"Region with label '{hsnp_region_label}' not found")
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 创建分析管理器，使用解析后的多个区域
    manager = create_combined_read_manager(
        bam_path,
        chrom,
        parsed_regions,  # 现在传递多个区域及其自定义标签
        reference_fasta
    )

    # 处理数据
    manager.clear_processors()
    
    # 如果指定了hsnp_region_label，为对应区域添加hsnp处理器
    if hsnp_region_label is not None and hsnp_info is not None:
        # 找到指定标签的区域
        hsnp_regions = [r for r in parsed_regions if r[2] == hsnp_region_label]
        if hsnp_regions:
            hsnp_start, hsnp_end, _ = hsnp_regions[0]
            print(f"Adding HSNP processor for region '{hsnp_region_label}' ({hsnp_start}-{hsnp_end})")
            # 为该区域添加hsnp处理器
            manager.add_processor(hsnp_processor(hsnp_info, (hsnp_start, hsnp_end)))
            
    # 增加str_data
    if str_region_info_path is not None and str_id is not None:
        str_region = create_str_region_from_tsv(str_region_info_path, str_id)
        str_id, str_reads = process_single_str_id(str_id, bam_path, reference_fasta, str_region)
        manager.add_processor(str_data_processor(str_reads.str_data))
    
    manager.process_reads()
    
    manager.sort_reads('start')
    
    if align_samples:
        print(f"align_samples: {align_samples}")
        ref_start = min(r[0] for r in parsed_regions) - 150
        ref_end = max(r[1] for r in parsed_regions) + 150
        print(f"ref_start: {ref_start}, ref_end: {ref_end}")
    else:
        ref_start = None
        ref_end = None
    read_df, base_df, reference_df = manager.generate_dataframes(ref_start=ref_start, ref_end=ref_end)
    
    sampler_info = manager.sampler_info
    
    # 生成输出路径
    file_paths = None
    
    # 生成参考序列
    if not read_df.empty:  
        # 生成输出路径
        file_paths = {
            'base_df': os.path.join(output_dir, f"{prefix}_base.csv"),
            'read_df': os.path.join(output_dir, f"{prefix}_read.csv"),
            'reference_df': os.path.join(output_dir, f"{prefix}_reference.csv"),
            'plot': os.path.join(output_dir, f"{prefix}.pdf")
        }
        
        if align_samples:
            ## cut all dataframes to same reference position range
            base_df = base_df[(base_df['ref_pos'] >= ref_start) & (base_df['ref_pos'] <= ref_end)]
            reference_df = reference_df[(reference_df['ref_pos'] >= ref_start) & (reference_df['ref_pos'] <= ref_end)]

        # 保存文件
        base_df.to_csv(file_paths['base_df'], index=False)
        read_df.to_csv(file_paths['read_df'], index=False)
        reference_df.to_csv(file_paths['reference_df'], index=False)
    else:
        run_visualization = False
    
    # 执行可视化
    if run_visualization and file_paths:
        # 根据参数决定使用哪个区域的坐标
        if mutation_region_label is not None:
            # 查找指定标签的区域
            target_regions = [r for r in parsed_regions if r[2] == mutation_region_label]
            if not target_regions:
                raise ValueError(f"Region with label '{mutation_region_label}' not found")
            # 使用指定标签的区域坐标
            mut_start = target_regions[0][0]
            mut_end = target_regions[0][1]
            print(f"Using mutation region with label '{mutation_region_label}': {mut_start}-{mut_end}")
        else:
            # 使用所有区域的范围（默认行为）
            mut_start = min(region[0] for region in parsed_regions)
            mut_end = max(region[1] for region in parsed_regions)
            print(f"Using full range for mutation region: {mut_start}-{mut_end}")
        
        # 准备可视化参数
        vis_kwargs = {
            'base_df_path': file_paths['base_df'],
            'read_df_path': file_paths['read_df'],
            'reference_df_path': file_paths['reference_df'],
            'output_file': file_paths['plot'],
            'mut_start': mut_start,
            'mut_end': mut_end,
            'mut_type': mut_type,
            'r_script_path': r_script_path,
            'start_axis': start_axis
        }
        
        # 从 kwargs 中提取可能的可视化参数
        for key in ['mut_type', 'start_axis']:
            if key in kwargs:
                vis_kwargs[key] = kwargs[key]
                
        run_r_visualization(**vis_kwargs)

        # 将pdf转换为png
        pdf_path = file_paths['plot']
        png_path = pdf_path.replace('.pdf', '.png')
        convert_pdf_to_png(pdf_path, png_path)
    
    return {
        'file_paths': file_paths,
        'sampler_info': sampler_info
    }

def convert_pdf_to_png(pdf_path: str, png_path: str, dpi: int = 300, retries: int = 3,
                      background_color: str = 'white', max_dpi: int = 1200) -> bool:
    """Convert PDF to PNG using ImageMagick with automatic DPI adjustment and error handling.
    
    Args:
        pdf_path (str): Path to input PDF file
        png_path (str): Path to output PNG file
        dpi (int): Initial DPI value for conversion (default: 300)
        retries (int): Number of retry attempts with reduced DPI (default: 3)
        background_color (str): Background color for the PNG (default: 'white')
        max_dpi (int): Maximum allowed DPI value (default: 1200)
        
    Returns:
        bool: True if conversion successful, False otherwise
        
    Raises:
        FileNotFoundError: If input PDF file does not exist
        ValueError: If DPI value is invalid
    """
    start_time = time.time()
    current_dpi = min(max_dpi, max(1, dpi))  # 确保DPI在有效范围内
    
    # 验证输入文件
    if not os.path.exists(pdf_path):
        raise FileNotFoundError(f"Input PDF file not found: {pdf_path}")
    
    # 确保输出目录存在
    os.makedirs(os.path.dirname(png_path), exist_ok=True)
    
    for attempt in range(retries + 1):
        try:
            # 构建ImageMagick命令
            cmd = [
                'convert',
                '-density', str(current_dpi),
                '-background', background_color,
                '-alpha', 'remove',
                '-quality', '95',
                '-strip',  # 移除元数据以减小文件大小
                pdf_path,
                png_path
            ]
            
            # 执行转换命令
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )
            
            # 验证输出文件
            if os.path.exists(png_path) and os.path.getsize(png_path) > 0:
                print(f"Successfully converted {pdf_path} to {png_path} (DPI: {current_dpi})")
                return True
            else:
                raise RuntimeError("Output file is empty or was not created")
                
        except subprocess.CalledProcessError as e:
            print(f"ERROR: Conversion failed at DPI {current_dpi}: {e.stderr}")
            if attempt < retries:
                # 指数退避策略
                current_dpi = max(1, current_dpi // 2)
                print(f"Retrying with DPI {current_dpi}...")
                time.sleep(1)  # 添加短暂延迟
            else:
                print("ERROR: Conversion failed after all attempts")
                return False
        except Exception as e:
            print(f"Unexpected error during PDF to PNG conversion: {str(e)}")
            return False
    
    return False

def convert_svg_to_png(svg_path: str, png_path: str, dpi: int = 300, retries: int = 3, 
                      background_color: str = 'white', max_dpi: int = 1200) -> bool:
    """Convert SVG to PNG with automatic DPI adjustment and error handling.
    
    Args:
        svg_path (str): Path to input SVG file
        png_path (str): Path to output PNG file
        dpi (int): Initial DPI value for conversion (default: 300)
        retries (int): Number of retry attempts with reduced DPI (default: 3)
        background_color (str): Background color for the PNG (default: 'white')
        max_dpi (int): Maximum allowed DPI value (default: 1200)
        
    Returns:
        bool: True if conversion successful, False otherwise
        
    Raises:
        FileNotFoundError: If input SVG file does not exist
        ValueError: If DPI value is invalid
    """
    start_time = time.time()
    current_dpi = min(max_dpi, max(1, dpi))  # 确保DPI在有效范围内
    
    # 验证输入文件
    if not os.path.exists(svg_path):
        raise FileNotFoundError(f"Input SVG file not found: {svg_path}")
    
    # 确保输出目录存在
    os.makedirs(os.path.dirname(png_path), exist_ok=True)
    
    for attempt in range(retries + 1):
        try:
            # 禁用PIL的最大图像像素限制
            Image.MAX_IMAGE_PIXELS = None
            
            # 读取SVG文件
            with open(svg_path, 'rb') as svg_file:
                svg_data = svg_file.read()
            
            # 使用cairosvg进行转换
            png_data = cairosvg.svg2png(
                bytestring=svg_data,
                dpi=current_dpi,
                background_color=background_color
            )
            
            # 使用PIL处理图像
            with Image.open(io.BytesIO(png_data)) as img:
                # 确保图像是RGB模式
                if img.mode in ('RGBA', 'LA'):
                    background = Image.new('RGB', img.size, background_color)
                    background.paste(img, mask=img.split()[-1])
                    img = background
                elif img.mode != 'RGB':
                    img = img.convert('RGB')
                
                # 保存图像
                img.save(
                    png_path,
                    "PNG",
                    dpi=(current_dpi, current_dpi),
                    optimize=True,
                    quality=95
                )
            
            # 验证输出文件
            if os.path.exists(png_path) and os.path.getsize(png_path) > 0:
                print(f"Successfully converted to PNG: {png_path} (DPI: {current_dpi})")
                return True
            else:
                raise RuntimeError("Output file is empty or was not created")
                
        except Exception as e:
            print(f"ERROR: Conversion failed at DPI {current_dpi}: {str(e)}")
            if attempt < retries:
                # 指数退避策略
                current_dpi = max(1, current_dpi // 2)
                print(f"Retrying with DPI {current_dpi}...")
                time.sleep(1)  # 添加短暂延迟
            else:
                print("ERROR: Conversion failed after all attempts")
                return False
    
    return False

def process_str_reads(bam_path: str, str_id: str, reference_fasta: str) -> Dict[str, StrReads]:
    """Process STR reads for a single STR site
    
    Args:
        bam_path: Path to BAM/CRAM file
        str_id: Single STR ID to process
        reference_fasta: Path to reference genome
        
    Returns:
        Dictionary containing processed STR reads for the given ID
    """
    str_reads_dict = {}
    
    # 直接处理单个str_id
    result = process_single_str_id(str_id, bam_path=bam_path, reference_fasta=reference_fasta)
    
    # 处理结果
    if result is not None:
        str_id, str_reads = result
        if str_reads is not None:
            str_reads_dict[str_id] = str_reads
    
    return str_reads_dict

def get_memory_usage() -> str:
    """Get current memory usage in a human-readable format"""
    process = psutil.Process(os.getpid())
    memory_mb = process.memory_info().rss / 1024 / 1024
    return f"{memory_mb:.2f} MB"

def process_single_str(
    bam_path: str,
    mode: str,
    ref_filename: str,
    single_region: Dict,
    flanking_seqs: Tuple[str, str]
) -> StrReads:
    """Process a single STR region"""
    left_flanking_three_bp, right_flanking_three_bp = flanking_seqs
    
    with pysam.AlignmentFile(bam_path, mode, reference_filename=ref_filename) as alignment_file:
        str_reads = StrReads()
        str_reads.compute_str_data(
            alignment_file,
            single_region['chrom'],
            single_region['start'],
            single_region['end']-1,
            single_region['str_id'],
            single_region['str_unit'],
            left_flanking_three_bp,
            right_flanking_three_bp
        )
        return str_reads

def process_single_str_id(str_id: str, bam_path: str, reference_fasta: str, str_region: STRRegion) -> Tuple[str, Optional[StrReads]]:
    """Process a single STR ID"""
    start_time = time.time()
    logging.info(f"Processing site: {str_id}")
    
    try:
        single_region = str_region
        
        flanking_seqs = get_flanking_sequences(reference_fasta, single_region)
        file_type = determine_file_type(bam_path)
        
        mode = get_file_mode(file_type)
        str_reads = process_single_str(
            bam_path, mode, reference_fasta, single_region, flanking_seqs
        )
        
        elapsed_time = time.time() - start_time
        logging.info(f"Completed processing {str_id} in {elapsed_time:.2f}s")
        logging.info(f"Memory usage: {get_memory_usage()}")
        
        return str_id, str_reads
        
    except Exception as e:
        logging.error(f"Error processing {str_id}: {str(e)}")
        logging.error(traceback.format_exc())
        return str_id, None

