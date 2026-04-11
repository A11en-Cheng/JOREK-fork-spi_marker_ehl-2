#!/usr/bin/env python3
"""
JOREK H5 to CSV Converter with Multi-threading
Converts jorek*.h5 sequences to fragments*.csv using jorek2_diagno_spi binary.
"""

import os
import sys
import re
import shutil
import subprocess
import threading
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Tuple, Set, Dict
import argparse

class JorekConverter:
    def __init__(self, binary_path: str, namelist_path: str, 
                 source_dir: str = ".", target_dir: str = "csv_output",
                 temp_base: str = "/tmp/jorek_convert", timeout: int = 300):
        self.binary_path = Path(binary_path).resolve()
        self.namelist_path = Path(namelist_path).resolve()
        self.source_dir = Path(source_dir).resolve()
        self.target_dir = Path(target_dir).resolve()
        self.temp_base = Path(temp_base)
        self.timeout = timeout
        
        # Create directories
        self.target_dir.mkdir(parents=True, exist_ok=True)
        self.temp_base.mkdir(parents=True, exist_ok=True)
        
        # Thread safety
        self.lock = threading.Lock()
        self.active_threads: Dict[int, str] = {}

    def extract_namelist_file_refs(self) -> Set[str]:
        """Extract likely filename references from the namelist text."""
        refs: Set[str] = set()

        try:
            text = self.namelist_path.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            return refs

        for raw_line in text.splitlines():
            line = raw_line.split("!", 1)[0].strip()
            if not line or "=" not in line:
                continue

            _, rhs = line.split("=", 1)

            # Quoted values in namelists are typical file paths.
            for quoted in re.findall(r"['\"]([^'\"]+)['\"]", rhs):
                val = quoted.strip()
                if val and val.lower() not in {"none", "__no_filename__"}:
                    refs.add(val)

            # Some namelists use unquoted tokens.
            for token in re.findall(r"\b[\w./-]+\b", rhs):
                lower = token.lower()
                if lower in {"t", "f", "true", "false", "none", "__no_filename__"}:
                    continue
                if token.isdigit():
                    continue
                if any(ch in token for ch in ("/", ".")) or token.startswith("jorek_"):
                    refs.add(token)

        return refs

    def link_auxiliary_inputs(self, thread_dir: Path, file_refs: Set[str]):
        """Link namelist-referenced input files into a thread workspace."""
        search_dirs = [self.namelist_path.parent, self.source_dir]

        for ref in file_refs:
            ref_path = Path(ref)
            target = thread_dir / ref_path.name
            if target.exists():
                continue

            candidates: List[Path] = []
            if ref_path.is_absolute():
                candidates.append(ref_path)
            else:
                for base in search_dirs:
                    candidates.append(base / ref)
                    candidates.append(base / ref_path.name)

            src = next((c for c in candidates if c.is_file()), None)
            if src is None:
                continue

            try:
                os.symlink(src, target)
            except OSError:
                # Symlink can fail on some filesystems; fallback to copy.
                shutil.copy2(src, target)
        
    def parse_step_range(self, range_str: str) -> Set[int]:
        """Parse step range specification like '1500-1700,1900' or '1500,1700,1900'"""
        steps = set()
        ranges = range_str.split(',')
        
        for r in ranges:
            if '-' in r:
                # Handle range like 1500-1900
                parts = r.split('-')
                if len(parts) == 2:
                    start, end = map(int, parts)
                    steps.update(range(start, end + 1))
                elif len(parts) == 3:
                    # Handle range with step like 1500-100-1900
                    start, step, end = map(int, parts)
                    steps.update(range(start, end + 1, step))
            else:
                # Single step
                steps.add(int(r))
        
        return steps

    def find_h5_files(self, steps: Set[int]) -> List[Tuple[int, Path]]:
        """Find H5 files matching the specified steps"""
        files = []
        pattern = re.compile(r'jorek(\d{6})\.h5')
        
        for h5_file in self.source_dir.glob("jorek*.h5"):
            match = pattern.match(h5_file.name)
            if match:
                step_num = int(match.group(1))
                if step_num in steps:
                    files.append((step_num, h5_file))
        
        # Sort by step number
        files.sort(key=lambda x: x[0])
        return files

    def setup_thread_workspace(self, thread_id: int) -> Path:
        """Create isolated workspace for a thread"""
        thread_dir = self.temp_base / f"thread_{thread_id}"
        thread_dir.mkdir(exist_ok=True)
        return thread_dir

    def cleanup_thread_workspace(self, thread_id: int):
        """Clean up thread workspace"""
        thread_dir = self.temp_base / f"thread_{thread_id}"
        if thread_dir.exists():
            shutil.rmtree(thread_dir)

    def convert_single_step(self, step_num: int, h5_file: Path, thread_id: int,
                            file_refs: Set[str]) -> bool:
        """Convert single H5 file to CSV"""
        thread_dir = self.setup_thread_workspace(thread_id)
        
        try:
            # Mark thread as active
            with self.lock:
                self.active_threads[thread_id] = f"step_{step_num:06d}"
            
            print(f"Thread {thread_id}: Processing step {step_num:06d}")
            
            # Copy H5 file to thread workspace as jorek_restart.h5
            target_h5 = thread_dir / "jorek_restart.h5"
            shutil.copy2(h5_file, target_h5)
            
            # Copy namelist file
            target_namelist = thread_dir / self.namelist_path.name
            shutil.copy2(self.namelist_path, target_namelist)

            # Link profile and other auxiliary inputs referenced by namelist.
            self.link_auxiliary_inputs(thread_dir, file_refs)
            
            try:
                # jorek2_diagno_spi reads namelist from stdin in this mode.
                with open(target_namelist, "r", encoding="utf-8") as namelist_fp:
                    result = subprocess.run(
                        [str(self.binary_path)],
                        cwd=thread_dir,
                        stdin=namelist_fp,
                        capture_output=True,
                        text=True,
                        timeout=self.timeout
                    )
                
                if result.returncode != 0:
                    print(f"Thread {thread_id}: Error converting step {step_num:06d}")
                    if result.stdout:
                        print(f"STDOUT: {result.stdout.strip()}")
                    print(f"STDERR: {result.stderr}")
                    return False
                
                # Find and rename output .dat files
                moved_count = self.process_output_files(step_num, thread_dir)
                if moved_count == 0:
                    print(f"Thread {thread_id}: No fragments*.dat produced for step {step_num:06d}")
                    if result.stdout:
                        print(f"STDOUT: {result.stdout.strip()}")
                    if result.stderr:
                        print(f"STDERR: {result.stderr.strip()}")
                    return False
                
                print(f"Thread {thread_id}: Successfully converted step {step_num:06d}")
                return True
                
            except subprocess.TimeoutExpired:
                print(f"Thread {thread_id}: Timeout converting step {step_num:06d} after {self.timeout}s")
                return False
            except Exception as e:
                print(f"Thread {thread_id}: Exception converting step {step_num:06d}: {e}")
                return False
                
        finally:
            # Clean up and mark thread as inactive
            self.cleanup_thread_workspace(thread_id)
            with self.lock:
                if thread_id in self.active_threads:
                    del self.active_threads[thread_id]

    def process_output_files(self, step_num: int, thread_dir: Path) -> int:
        """Rename and move output .dat files to .csv with step number"""
        moved_count = 0
        for dat_file in thread_dir.glob("fragments*.dat"):
            # Create new filename with step number
            stem = dat_file.stem
            new_name = f"{stem}_{step_num:06d}.csv"
            target_path = self.target_dir / new_name
            
            # Rename and move file
            shutil.move(str(dat_file), str(target_path))

            moved_count += 1

        return moved_count

    def convert_parallel(self, step_range: str, max_workers: int = None):
        """Main conversion function with parallel processing"""
        file_refs = self.extract_namelist_file_refs()

        # Parse step range
        steps = self.parse_step_range(step_range)
        print(f"Found {len(steps)} steps to process: {sorted(steps)}")
        
        # Find matching H5 files
        h5_files = self.find_h5_files(steps)
        if not h5_files:
            print("No matching H5 files found!")
            return
        
        print(f"Found {len(h5_files)} matching H5 files")
        
        # Determine number of workers
        if max_workers is None:
            max_workers = min(len(h5_files), os.cpu_count() or 4)
        
        print(f"Using {max_workers} threads for conversion")
        print(f"Target directory: {self.target_dir}")
        
        # Thread counter
        thread_counter = 0
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all conversion tasks
            future_to_step = {}
            for step_num, h5_file in h5_files:
                thread_counter += 1
                thread_id = thread_counter
                future = executor.submit(
                    self.convert_single_step, 
                    step_num, 
                    h5_file, 
                    thread_id,
                    file_refs
                )
                future_to_step[future] = (step_num, thread_id)
            
            # Monitor progress
            completed = 0
            failed = 0
            total = len(future_to_step)
            
            for future in as_completed(future_to_step):
                step_num, thread_id = future_to_step[future]
                try:
                    success = future.result()
                    if success:
                        completed += 1
                    else:
                        failed += 1
                except Exception as e:
                    print(f"Unexpected error processing step {step_num:06d}: {e}")
                    failed += 1
                
                # Print progress
                print(f"Progress: {completed + failed}/{total} completed "
                      f"({completed} successful, {failed} failed)")
                
                # Print active threads
                with self.lock:
                    if self.active_threads:
                        active_str = ", ".join([f"thread {tid}({step})" 
                                              for tid, step in self.active_threads.items()])
                        print(f"Active: {active_str}")

        print(f"\nConversion completed: {completed} successful, {failed} failed")
        
        if failed > 0:
            print("Some conversions failed. Check the output above for details.")

    def cleanup_all(self):
        """Clean up all temporary files"""
        if self.temp_base.exists():
            shutil.rmtree(self.temp_base)
            print(f"Cleaned up temporary directory: {self.temp_base}")

def main():
    parser = argparse.ArgumentParser(
        description="Convert JOREK H5 sequences to CSV using jorek2_diagno_spi",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Step range examples:
  1500,1700,1900        - Specific steps
  1500-1900             - Range from 1500 to 1900
  1500-1900,2000-3000   - Multiple ranges
  1500-100-1900         - Range with step size (1500, 1600, 1700, 1800, 1900)
  1500-1900,3100,3250   - Mixed specification
        """
    )
    
    parser.add_argument("binary", help="Path to jorek2_diagno_spi binary")
    parser.add_argument("namelist", help="Path to namelist file")
    parser.add_argument("steps", help="Step range specification")
    
    parser.add_argument("-j", "--threads", type=int, 
                       help="Number of threads (default: CPU count)")
    parser.add_argument("-s", "--source-dir", default=".", 
                       help="Source directory with H5 files (default: current)")
    parser.add_argument("-t", "--target-dir", default="csv_output",
                       help="Target directory for CSV files (default: csv_output)")
    parser.add_argument("--temp-dir", default="/tmp/jorek_convert",
                       help="Temporary directory (default: /tmp/jorek_convert)")
    parser.add_argument("--timeout", type=int, default=300,
                       help="Timeout per step in seconds (default: 300)")
    parser.add_argument("--keep-temp", action="store_true",
                       help="Keep temporary files for debugging")
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.isfile(args.binary):
        print(f"Error: Binary file not found: {args.binary}")
        sys.exit(1)
    
    if not os.path.isfile(args.namelist):
        print(f"Error: Namelist file not found: {args.namelist}")
        sys.exit(1)
    
    if not os.path.isdir(args.source_dir):
        print(f"Error: Source directory not found: {args.source_dir}")
        sys.exit(1)
    
    # Create converter
    converter = JorekConverter(
        binary_path=args.binary,
        namelist_path=args.namelist,
        source_dir=args.source_dir,
        target_dir=args.target_dir,
        temp_base=args.temp_dir,
        timeout=args.timeout
    )
    
    try:
        # Run conversion
        converter.convert_parallel(args.steps, args.threads)
        
    except KeyboardInterrupt:
        print("\nConversion interrupted by user")
    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)
    finally:
        # Cleanup
        if not args.keep_temp:
            converter.cleanup_all()

if __name__ == "__main__":
    main()