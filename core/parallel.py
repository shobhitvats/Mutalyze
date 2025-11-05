"""
Parallel Execution Module
Handles multiprocessing for mutation analysis
"""

import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from typing import List, Dict, Callable, Any, Optional
from pathlib import Path
import logging
from functools import partial

try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False
    logging.warning("tqdm not available. Progress bars will be disabled.")

logger = logging.getLogger(__name__)


class ParallelExecutor:
    """Manages parallel execution of mutation calculations"""
    
    def __init__(self, max_workers: Optional[int] = None, use_threads: bool = False):
        """
        Initialize parallel executor
        
        Args:
            max_workers: Maximum number of parallel workers (None = CPU count)
            use_threads: If True, use threads instead of processes
        """
        if max_workers is None:
            max_workers = mp.cpu_count()
        
        self.max_workers = max_workers
        self.use_threads = use_threads
        
        logger.info(f"Initialized parallel executor with {max_workers} workers "
                   f"({'threads' if use_threads else 'processes'})")
    
    def map_mutations(self, func: Callable, mutations: List[Dict], 
                     show_progress: bool = True,
                     chunk_size: int = 1) -> List[Any]:
        """
        Apply a function to mutations in parallel
        
        Args:
            func: Function to apply to each mutation
            mutations: List of mutation dictionaries
            show_progress: Show progress bar
            chunk_size: Chunk size for batching
            
        Returns:
            List of results
        """
        executor_class = ThreadPoolExecutor if self.use_threads else ProcessPoolExecutor
        
        results = []
        
        with executor_class(max_workers=self.max_workers) as executor:
            # Submit all tasks
            futures = {executor.submit(func, mutation): i 
                      for i, mutation in enumerate(mutations)}
            
            # Collect results with progress bar
            if show_progress and TQDM_AVAILABLE:
                progress = tqdm(total=len(mutations), 
                              desc="Processing mutations",
                              unit="mutation")
            else:
                progress = None
            
            for future in as_completed(futures):
                idx = futures[future]
                try:
                    result = future.result()
                    results.append((idx, result))
                except Exception as e:
                    logger.error(f"Mutation {idx} failed: {e}")
                    results.append((idx, None))
                
                if progress:
                    progress.update(1)
            
            if progress:
                progress.close()
        
        # Sort by original index
        results.sort(key=lambda x: x[0])
        return [r[1] for r in results]
    
    def map_with_context(self, func: Callable, mutations: List[Dict],
                        context: Dict, show_progress: bool = True) -> List[Any]:
        """
        Apply a function with shared context to mutations in parallel
        
        Args:
            func: Function that takes (mutation, context)
            mutations: List of mutations
            context: Shared context dictionary
            show_progress: Show progress bar
            
        Returns:
            List of results
        """
        # Create partial function with context
        func_with_context = partial(func, context=context)
        
        return self.map_mutations(func_with_context, mutations, show_progress)
    
    def batch_process(self, func: Callable, items: List[Any],
                     batch_size: int = 10, 
                     show_progress: bool = True) -> List[Any]:
        """
        Process items in batches
        
        Args:
            func: Function to apply to each batch
            items: List of items to process
            batch_size: Number of items per batch
            show_progress: Show progress bar
            
        Returns:
            Flattened list of results
        """
        # Create batches
        batches = [items[i:i+batch_size] for i in range(0, len(items), batch_size)]
        
        executor_class = ThreadPoolExecutor if self.use_threads else ProcessPoolExecutor
        
        all_results = []
        
        with executor_class(max_workers=self.max_workers) as executor:
            futures = {executor.submit(func, batch): i 
                      for i, batch in enumerate(batches)}
            
            if show_progress and TQDM_AVAILABLE:
                progress = tqdm(total=len(batches), 
                              desc="Processing batches",
                              unit="batch")
            else:
                progress = None
            
            for future in as_completed(futures):
                idx = futures[future]
                try:
                    batch_results = future.result()
                    all_results.append((idx, batch_results))
                except Exception as e:
                    logger.error(f"Batch {idx} failed: {e}")
                    all_results.append((idx, []))
                
                if progress:
                    progress.update(1)
            
            if progress:
                progress.close()
        
        # Sort and flatten
        all_results.sort(key=lambda x: x[0])
        flattened = []
        for _, batch_results in all_results:
            if isinstance(batch_results, list):
                flattened.extend(batch_results)
            else:
                flattened.append(batch_results)
        
        return flattened


class MutationTaskManager:
    """Manages mutation analysis tasks with caching and recovery"""
    
    def __init__(self, cache_dir: str = "data/mutation_cache"):
        """
        Initialize task manager
        
        Args:
            cache_dir: Directory for caching results
        """
        import pickle
        from pathlib import Path
        
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.pickle = pickle
        
    def get_cache_path(self, task_id: str) -> Path:
        """Get cache file path for a task"""
        return self.cache_dir / f"{task_id}.pkl"
    
    def save_results(self, task_id: str, results: Any):
        """Save results to cache"""
        cache_path = self.get_cache_path(task_id)
        try:
            with open(cache_path, 'wb') as f:
                self.pickle.dump(results, f)
            logger.debug(f"Saved results to cache: {task_id}")
        except Exception as e:
            logger.warning(f"Failed to cache results: {e}")
    
    def load_results(self, task_id: str) -> Optional[Any]:
        """Load results from cache"""
        cache_path = self.get_cache_path(task_id)
        if cache_path.exists():
            try:
                with open(cache_path, 'rb') as f:
                    results = self.pickle.load(f)
                logger.info(f"Loaded cached results: {task_id}")
                return results
            except Exception as e:
                logger.warning(f"Failed to load cached results: {e}")
        return None
    
    def clear_cache(self, task_id: Optional[str] = None):
        """Clear cache for a specific task or all tasks"""
        if task_id:
            cache_path = self.get_cache_path(task_id)
            if cache_path.exists():
                cache_path.unlink()
                logger.info(f"Cleared cache: {task_id}")
        else:
            for cache_file in self.cache_dir.glob("*.pkl"):
                cache_file.unlink()
            logger.info("Cleared all cache")


def parallel_map(func: Callable, items: List[Any], 
                max_workers: Optional[int] = None,
                show_progress: bool = True) -> List[Any]:
    """
    Convenience function for parallel mapping
    
    Args:
        func: Function to apply
        items: List of items
        max_workers: Maximum parallel workers
        show_progress: Show progress bar
        
    Returns:
        List of results
    """
    executor = ParallelExecutor(max_workers=max_workers)
    return executor.map_mutations(func, items, show_progress)


def calculate_optimal_workers(num_tasks: int, 
                             max_workers: Optional[int] = None) -> int:
    """
    Calculate optimal number of workers for a given number of tasks
    
    Args:
        num_tasks: Number of tasks to execute
        max_workers: Maximum allowed workers
        
    Returns:
        Optimal number of workers
    """
    cpu_count = mp.cpu_count()
    
    if max_workers is None:
        max_workers = cpu_count
    
    # Use fewer workers for small tasks
    if num_tasks < cpu_count:
        optimal = max(1, num_tasks)
    else:
        optimal = min(cpu_count, max_workers)
    
    return optimal
