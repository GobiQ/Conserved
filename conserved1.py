import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
try:
    from Bio.SeqUtils import GC
except ImportError:
    # For newer versions of Biopython
    def GC(seq):
        seq = seq.upper()
        gc_count = seq.count('G') + seq.count('C')
        total_count = len(seq)
        if total_count == 0:
            return 0.0
        return (gc_count / total_count) * 100

import requests
import io
import time
from typing import List, Dict, Tuple, Optional
import re
from collections import defaultdict

# Set page config
st.set_page_config(
    page_title="Genome Conservation Analysis Tool",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

class NCBIGenomeAnalyzer:
    """Class to handle NCBI data retrieval and conservation analysis"""
    
    def __init__(self, email: str, organism_type: str = "Bacteria"):
        Entrez.email = email
        self.conservation_threshold = 0.8
        self.organism_type = organism_type
        
    def search_genomes(self, species: str, max_results: int = 50) -> List[Dict]:
        """Search for genome sequences directly in nucleotide database"""
        try:
            st.info(f"Searching nucleotide database for {species} genomes...")
            
            # Organism-specific search strategies
            if self.organism_type == "Viroid":
                search_terms = self._get_viroid_search_terms(species)
                size_filter = lambda x: 50 <= x <= 2000  # Very broad viroid range
            elif self.organism_type == "Virus":
                search_terms = self._get_virus_search_terms(species)
                size_filter = lambda x: 1000 <= x <= 500000
            else:
                search_terms = self._get_standard_search_terms(species)
                size_filter = lambda x: x >= 10000
            
            sequences = []
            
            for i, search_term in enumerate(search_terms):
                try:
                    st.info(f"Search strategy {i+1}/{len(search_terms)}: {search_term}")
                    
                    handle = Entrez.esearch(
                        db="nucleotide", 
                        term=search_term, 
                        retmax=max_results * 2,
                        sort="relevance"
                    )
                    search_results = Entrez.read(handle)
                    handle.close()
                    
                    if not search_results.get('IdList'):
                        st.warning(f"No results for strategy {i+1}")
                        continue
                    
                    st.success(f"Found {len(search_results['IdList'])} sequences")
                    
                    # Process sequences with better error handling
                    batch_sequences = self._process_sequence_batch(
                        search_results['IdList'], species, size_filter, max_results
                    )
                    sequences.extend(batch_sequences)
                    
                    if len(sequences) >= max_results:
                        break
                        
                except Exception as search_error:
                    st.warning(f"Search strategy {i+1} failed: {search_error}")
                    continue
            
            # Remove duplicates
            unique_sequences = self._deduplicate_sequences(sequences, max_results)
            
            st.success(f"Found {len(unique_sequences)} unique genome sequences")
            return unique_sequences
            
        except Exception as e:
            st.error(f"Error in genome search: {str(e)}")
            # Create test entry to ensure UI works
            return [{
                'sequence_id': 'TEST_001',
                'title': f'Test {species} Sequence',
                'organism': species,
                'length': 300,
                'source': 'test'
            }]
    
    def _get_viroid_search_terms(self, species: str) -> List[str]:
        """Generate comprehensive viroid-specific search terms"""
        # Extract key terms from species name
        terms = species.lower().split()
        base_terms = []
        
        # Handle different viroid naming patterns
        if "viroid" in species.lower():
            # Already contains viroid
            base_name = species.replace(" viroid", "").replace(" Viroid", "")
            base_terms.append(base_name)
        else:
            base_terms.append(species)
        
        # Add common abbreviations
        if "hop latent" in species.lower():
            base_terms.extend(["HpLVd", "HLVd", "hop latent"])
        elif "potato spindle tuber" in species.lower():
            base_terms.extend(["PSTVd", "potato spindle tuber"])
        elif "citrus exocortis" in species.lower():
            base_terms.extend(["CEVd", "citrus exocortis"])
        elif "chrysanthemum stunt" in species.lower():
            base_terms.extend(["CSVd", "chrysanthemum stunt"])
        
        search_terms = []
        
        # Build comprehensive search strategy
        for base_term in base_terms:
            search_terms.extend([
                f'"{base_term}"[Organism] AND viroid',
                f'{base_term}[Organism] AND viroid',
                f'"{base_term} viroid"[Title]',
                f'{base_term} viroid[Title]',
                f'"{base_term}" AND viroid AND complete',
                f'{base_term} viroid complete',
                f'{base_term} viroid sequence',
                f'"{base_term}"[Organism]',
                f'{base_term}[Organism]'
            ])
        
        # Add generic viroid searches
        if len(terms) > 1:
            genus = terms[0]
            search_terms.extend([
                f'{genus}[Organism] AND viroid',
                f'viroid AND {genus}',
                f'{genus} viroid'
            ])
        
        # Remove duplicates while preserving order
        unique_terms = []
        for term in search_terms:
            if term not in unique_terms:
                unique_terms.append(term)
        
        return unique_terms
    
    def _get_virus_search_terms(self, species: str) -> List[str]:
        """Generate comprehensive virus-specific search terms"""
        base_terms = [species]
        
        # Add common virus abbreviations
        if "tobacco mosaic" in species.lower():
            base_terms.extend(["TMV", "tobacco mosaic"])
        elif "sars" in species.lower():
            base_terms.extend(["SARS-CoV-2", "SARS-CoV", "coronavirus"])
        elif "influenza" in species.lower():
            base_terms.extend(["flu", "influenza"])
        
        search_terms = []
        
        for base_term in base_terms:
            search_terms.extend([
                f'"{base_term}"[Organism] AND "complete genome"',
                f'"{base_term}"[Organism] AND genome',
                f'{base_term}[Organism] AND "complete genome"',
                f'{base_term}[Organism] AND genome',
                f'"{base_term}" AND "complete genome"',
                f'{base_term} complete genome',
                f'{base_term} genome',
                f'"{base_term}"[Organism]',
                f'{base_term}[Organism]'
            ])
        
        # Add genus-level searches
        terms = species.split()
        if len(terms) > 1:
            genus = terms[0]
            search_terms.extend([
                f'{genus}[Organism] AND "complete genome"',
                f'{genus}[Organism] AND virus',
                f'{genus} virus'
            ])
        
        return list(dict.fromkeys(search_terms))  # Remove duplicates
    
    def _get_standard_search_terms(self, species: str) -> List[str]:
        """Generate comprehensive search terms for bacteria/eukaryotes"""
        base_terms = [species]
        
        search_terms = []
        
        for base_term in base_terms:
            search_terms.extend([
                f'"{base_term}"[Organism] AND "complete genome"',
                f'"{base_term}"[Organism] AND chromosome',
                f'"{base_term}"[Organism] AND "reference genome"',
                f'{base_term}[Organism] AND "complete genome"',
                f'{base_term}[Organism] AND chromosome',
                f'{base_term}[Organism] AND "reference genome"',
                f'"{base_term}" AND "complete genome"',
                f'{base_term} complete genome',
                f'{base_term} chromosome',
                f'{base_term} genome',
                f'"{base_term}"[Organism]',
                f'{base_term}[Organism]'
            ])
        
        # Add genus and species-level searches
        terms = species.split()
        if len(terms) >= 2:
            genus = terms[0]
            species_epithet = terms[1]
            search_terms.extend([
                f'{genus}[Organism] AND "complete genome"',
                f'{genus}[Organism] AND chromosome',
                f'{genus} {species_epithet}[Organism]',
                f'{genus} {species_epithet} genome'
            ])
        
        return list(dict.fromkeys(search_terms))  # Remove duplicates
    
    def _process_sequence_batch(self, id_list: List[str], species: str, size_filter, max_results: int) -> List[Dict]:
        """Process sequence IDs with robust error handling"""
        sequences = []
        
        try:
            # Process in smaller batches
            batch_size = 10
            for i in range(0, min(len(id_list), 50), batch_size):  # Limit to first 50
                batch_ids = id_list[i:i+batch_size]
                
                try:
                    handle = Entrez.esummary(db="nucleotide", id=','.join(batch_ids))
                    summaries = Entrez.read(handle)
                    handle.close()
                    
                    for summary in summaries:
                        try:
                            seq_id = str(summary.get('AccessionVersion', summary.get('Id', 'Unknown')))
                            title = str(summary.get('Title', 'Unknown'))
                            length = int(summary.get('Length', 0))
                            organism = str(summary.get('Organism', species))
                            
                            if size_filter(length) and seq_id != 'Unknown':
                                sequences.append({
                                    'sequence_id': seq_id,
                                    'title': title,
                                    'organism': organism,
                                    'length': length,
                                    'source': 'nucleotide'
                                })
                            
                            if len(sequences) >= max_results:
                                break
                                
                        except (ValueError, TypeError, KeyError):
                            continue
                
                except Exception as batch_error:
                    st.warning(f"Batch processing error: {batch_error}")
                    continue
                
                if len(sequences) >= max_results:
                    break
        
        except Exception as e:
            st.warning(f"Processing failed, creating fallback entries: {e}")
            # Create fallback entries
            for i, seq_id in enumerate(id_list[:10]):
                sequences.append({
                    'sequence_id': seq_id,
                    'title': f'{species} sequence {i+1}',
                    'organism': species,
                    'length': 300,
                    'source': 'fallback'
                })
        
        return sequences
    
    def _deduplicate_sequences(self, sequences: List[Dict], max_results: int) -> List[Dict]:
        """Remove duplicates and sort sequences"""
        seen_ids = set()
        unique_sequences = []
        
        # Sort appropriately for organism type
        if self.organism_type == "Viroid":
            sequences.sort(key=lambda x: abs(x.get('length', 300) - 300))
        else:
            sequences.sort(key=lambda x: x.get('length', 0), reverse=True)
        
        for seq in sequences:
            seq_id = seq.get('sequence_id', '').strip()
            if seq_id and seq_id != 'Unknown' and seq_id not in seen_ids:
                seen_ids.add(seq_id)
                unique_sequences.append(seq)
                if len(unique_sequences) >= max_results:
                    break
        
        return unique_sequences
    
    def fetch_all_sequences_simultaneously(self, sequence_ids: List[str]) -> Dict[str, str]:
        """Fetch ALL sequences at once for true comparative analysis"""
        sequences = {}
        
        st.info(f"Fetching {len(sequence_ids)} sequences for comparative analysis...")
        progress_bar = st.progress(0)
        
        for i, sequence_id in enumerate(sequence_ids):
            try:
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=sequence_id,
                    rettype="fasta",
                    retmode="text"
                )
                sequence_data = handle.read()
                handle.close()
                
                parsed_sequences = list(SeqIO.parse(io.StringIO(sequence_data), "fasta"))
                if parsed_sequences:
                    sequence = str(parsed_sequences[0].seq).upper()
                    sequences[sequence_id] = sequence
                    
                    if len(sequences) <= 3:
                        st.success(f"Fetched {sequence_id}: {len(sequence)} bp")
                
                progress_bar.progress((i + 1) / len(sequence_ids))
                
            except Exception as e:
                st.warning(f"Failed to fetch {sequence_id}: {str(e)}")
                continue
        
        st.success(f"Successfully fetched {len(sequences)}/{len(sequence_ids)} sequences")
        return sequences
    
    def generate_test_sequences(self, species: str, num_sequences: int = 5) -> Dict[str, str]:
        """Generate synthetic test sequences with known conservation patterns"""
        st.info("Generating synthetic test sequences with known conservation patterns")
        
        if self.organism_type == "Viroid":
            base_length = 300
            # Create a conserved core with variable regions
            conserved_regions = [
                "GGAACCTTGTCGGATCCGAGGATCCGTCGAAC",  # Highly conserved
                "CTCGAGCTTGGACTCGAGCT",               # Moderately conserved
                "GGCCTTAACCGGTT"                     # Highly conserved
            ]
        else:
            base_length = 1000
            conserved_regions = [
                "ATGCGATCGGATCCGAGGATCCGTCGAACCTGG",
                "TTCGAGCTTGGACTCGAGCTATGCC",
                "GGCCTTAACCGGTTAAGCTT",
                "CGATCGATCGTAGCTAGCTA"
            ]
        
        sequences = {}
        
        for i in range(num_sequences):
            sequence = ""
            pos = 0
            
            for j, conserved_seq in enumerate(conserved_regions):
                # Add variable region before conserved region
                if j > 0:
                    var_length = 50 if self.organism_type == "Viroid" else 200
                    variable_region = self._generate_random_sequence(var_length, f"variable_{j}_{i}")
                    sequence += variable_region
                    pos += var_length
                
                # Add conserved region with some mutations
                if i == 0:
                    # Reference sequence - no mutations
                    sequence += conserved_seq
                else:
                    # Introduce some mutations (10-20% depending on region)
                    mutation_rate = 0.1 if j % 2 == 0 else 0.2  # Alternate high/low conservation
                    mutated_seq = self._introduce_mutations(conserved_seq, mutation_rate, i)
                    sequence += mutated_seq
                
                pos += len(conserved_seq)
            
            # Fill to target length
            if len(sequence) < base_length:
                remaining = base_length - len(sequence)
                sequence += self._generate_random_sequence(remaining, f"tail_{i}")
            
            sequences[f"TEST_SEQ_{i+1:02d}"] = sequence[:base_length]
        
        st.success(f"Generated {len(sequences)} test sequences with known conservation patterns")
        return sequences

    def _generate_random_sequence(self, length: int, seed_str: str) -> str:
        """Generate random DNA sequence with consistent seed"""
        import random
        random.seed(hash(seed_str) % (2**32))
        bases = 'ATGC'
        return ''.join(random.choice(bases) for _ in range(length))

    def _introduce_mutations(self, sequence: str, mutation_rate: float, seq_index: int) -> str:
        """Introduce random mutations at specified rate"""
        import random
        random.seed((hash(sequence) + seq_index) % (2**32))
        bases = 'ATGC'
        mutated = list(sequence)
        
        for i in range(len(mutated)):
            if random.random() < mutation_rate:
                current_base = mutated[i]
                available_bases = [b for b in bases if b != current_base]
                mutated[i] = random.choice(available_bases)
        
        return ''.join(mutated)
    
    def true_comparative_analysis(self, sequences: Dict[str, str], window_size: int = 25, step_size: int = 5) -> pd.DataFrame:
        """Perform comparative analysis with improved alignment strategy"""
        if len(sequences) < 2:
            st.error("Need at least 2 sequences for comparative analysis")
            return pd.DataFrame()
        
        st.info(f"Performing comparative analysis on {len(sequences)} sequences...")
        
        seq_ids = list(sequences.keys())
        seq_list = list(sequences.values())
        
        # Better length handling - use overlapping regions only
        lengths = [len(seq) for seq in seq_list]
        min_length = min(lengths)
        max_length = max(lengths)
        
        st.info(f"Sequence lengths: {min_length} - {max_length} bp")
        
        # Choose analysis strategy based on length variation
        length_variation = (max_length - min_length) / min_length
        
        if length_variation > 0.1:  # >10% variation
            # Use minimum length for fair comparison
            analysis_length = min_length
            st.info(f"High length variation detected. Using {analysis_length} bp for comparison")
        else:
            # Use common length
            analysis_length = min_length
        
        # Ensure window size is appropriate
        if window_size > analysis_length:
            window_size = max(5, analysis_length // 10)
            step_size = max(1, window_size // 5)
            st.warning(f"Adjusted window size to {window_size}bp for short sequences")
        
        # Prepare sequences (trim to common length, no padding)
        aligned_sequences = []
        for seq_id, sequence in sequences.items():
            trimmed_seq = sequence[:analysis_length].upper()
            aligned_sequences.append(trimmed_seq)
            
        # Sliding window analysis
        results = []
        total_windows = max(1, (analysis_length - window_size) // step_size + 1)
        progress_bar = st.progress(0)
        
        for i, pos in enumerate(range(0, analysis_length - window_size + 1, step_size)):
            windows = [seq[pos:pos + window_size] for seq in aligned_sequences]
            
            # Skip windows with too many gaps or Ns
            valid_windows = []
            for window in windows:
                n_count = window.count('N') + window.count('-')
                if n_count / len(window) < 0.3:  # Less than 30% gaps/Ns
                    valid_windows.append(window)
            
            if len(valid_windows) < 2:
                continue
                
            conservation_score = self._calculate_conservation_score(valid_windows)
            identity_percentage = self._calculate_identity_percentage(valid_windows)
            consensus_sequence = self._generate_consensus(valid_windows)
            
            # Calculate additional metrics
            gc_content = self._calculate_gc_content(consensus_sequence)
            complexity_score = self._calculate_sequence_complexity(valid_windows)
            
            results.append({
                'start': pos + 1,
                'end': pos + window_size,
                'conservation_score': conservation_score,
                'identity_percentage': identity_percentage,
                'consensus_sequence': consensus_sequence,
                'gc_content': gc_content,
                'complexity_score': complexity_score,
                'num_sequences': len(valid_windows),
                'length': window_size,
                'individual_windows': valid_windows
            })
            
            if i % 50 == 0:
                progress_bar.progress(min(i / total_windows, 1.0))
        
        progress_bar.progress(1.0)
        df = pd.DataFrame(results)
        st.success(f"Analysis complete: {len(df)} windows analyzed from {analysis_length} bp")
        return df
    
    def _calculate_conservation_score(self, windows: List[str]) -> float:
        """Calculate conservation score with weighted scoring instead of binary threshold"""
        if not windows or len(windows) < 2:
            return 0.0
        
        window_length = len(windows[0])
        total_conservation = 0.0
        
        for pos in range(window_length):
            bases = [window[pos] for window in windows if pos < len(window)]
            if not bases:
                continue
            
            # Count each base type
            base_counts = {}
            for base in bases:
                base_counts[base] = base_counts.get(base, 0) + 1
            
            # Calculate conservation as the fraction of the most common base
            max_count = max(base_counts.values())
            conservation_ratio = max_count / len(bases)
            
            # Use weighted scoring instead of binary threshold
            if conservation_ratio >= 0.9:
                position_score = 1.0
            elif conservation_ratio >= 0.8:
                position_score = 0.8
            elif conservation_ratio >= 0.7:
                position_score = 0.6
            elif conservation_ratio >= 0.6:
                position_score = 0.4
            elif conservation_ratio >= 0.5:
                position_score = 0.2
            else:
                position_score = 0.0
            
            total_conservation += position_score
        
        return total_conservation / window_length
    
    def _calculate_identity_percentage(self, windows: List[str]) -> float:
        """Calculate percentage where sequences are highly similar (not requiring 100% identity)"""
        if not windows or len(windows) < 2:
            return 0.0
        
        window_length = len(windows[0])
        similar_positions = 0
        
        for pos in range(window_length):
            bases = [window[pos] for window in windows if pos < len(window)]
            if not bases:
                continue
                
            # Count base frequencies
            base_counts = {}
            for base in bases:
                base_counts[base] = base_counts.get(base, 0) + 1
            
            # Consider position "similar" if dominant base is >70%
            max_count = max(base_counts.values())
            if max_count / len(bases) >= 0.7:
                similar_positions += 1
        
        return (similar_positions / window_length) * 100
    
    def _generate_consensus(self, windows: List[str]) -> str:
        """Generate consensus sequence"""
        if not windows:
            return ""
        
        consensus = ""
        window_length = max(len(window) for window in windows)
        
        for pos in range(window_length):
            bases = [window[pos] for window in windows if pos < len(window)]
            if not bases:
                consensus += "N"
                continue
            
            base_counts = {}
            for base in bases:
                base_counts[base] = base_counts.get(base, 0) + 1
            
            most_common_base = max(base_counts, key=base_counts.get)
            consensus += most_common_base
        
        return consensus
    
    def _calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC content safely"""
        if not sequence:
            return 0.0
        valid_bases = [b for b in sequence.upper() if b in 'ATGCU']
        if not valid_bases:
            return 0.0
        gc_count = sum(1 for b in valid_bases if b in 'GC')
        return (gc_count / len(valid_bases)) * 100

    def _calculate_sequence_complexity(self, windows: List[str]) -> float:
        """Calculate sequence complexity (Shannon entropy)"""
        if not windows:
            return 0.0
        
        # Combine all windows to calculate overall complexity
        combined_seq = ''.join(windows)
        
        # Count base frequencies
        base_counts = {}
        for base in combined_seq:
            if base in 'ATGCU':
                base_counts[base] = base_counts.get(base, 0) + 1
        
        if not base_counts:
            return 0.0
        
        total_bases = sum(base_counts.values())
        entropy = 0.0
        
        for count in base_counts.values():
            if count > 0:
                p = count / total_bases
                entropy -= p * np.log2(p)
        
        return entropy

    def compare_custom_sequence(self, custom_sequence: str, loaded_sequences: Dict[str, str], 
                              min_match_length: int = 20, max_mismatches: int = 2) -> pd.DataFrame:
        """Compare a custom sequence against all loaded genome sequences"""
        if not custom_sequence or len(custom_sequence) < min_match_length:
            return pd.DataFrame()
        
        custom_sequence = custom_sequence.upper().replace(' ', '').replace('\n', '')
        
        results = []
        
        for seq_id, genome_sequence in loaded_sequences.items():
            # Find all matches of the custom sequence in the genome
            matches = self._find_sequence_matches(custom_sequence, genome_sequence, 
                                                min_match_length, max_mismatches)
            
            for match in matches:
                results.append({
                    'sequence_id': seq_id,
                    'match_start': match['start'],
                    'match_end': match['end'],
                    'match_length': match['length'],
                    'mismatches': match['mismatches'],
                    'identity_percentage': match['identity'],
                    'match_sequence': match['matched_seq'],
                    'genome_region': genome_sequence[match['start']:match['end']],
                    'context_start': max(0, match['start'] - 50),
                    'context_end': min(len(genome_sequence), match['end'] + 50),
                    'context_sequence': genome_sequence[max(0, match['start'] - 50):min(len(genome_sequence), match['end'] + 50)]
                })
        
        return pd.DataFrame(results)
    
    def _find_sequence_matches(self, query_seq: str, target_seq: str, 
                             min_length: int, max_mismatches: int) -> List[Dict]:
        """Find approximate matches of query sequence in target sequence"""
        matches = []
        query_len = len(query_seq)
        
        # Use sliding window approach for approximate matching
        for i in range(len(target_seq) - min_length + 1):
            for window_len in range(min_length, min(query_len + 1, len(target_seq) - i + 1)):
                target_window = target_seq[i:i + window_len]
                
                # Calculate mismatches for this window
                mismatches = self._calculate_mismatches(query_seq[:window_len], target_window)
                
                if mismatches <= max_mismatches:
                    identity = ((window_len - mismatches) / window_len) * 100
                    
                    matches.append({
                        'start': i + 1,  # 1-based
                        'end': i + window_len,
                        'length': window_len,
                        'mismatches': mismatches,
                        'identity': identity,
                        'matched_seq': target_window
                    })
        
        # Sort by identity percentage (highest first), then by length
        matches.sort(key=lambda x: (-x['identity'], -x['length']))
        
        # Remove overlapping matches (keep the best one)
        filtered_matches = []
        for match in matches:
            is_overlapping = False
            for existing in filtered_matches:
                if (match['start'] <= existing['end'] and match['end'] >= existing['start']):
                    is_overlapping = True
                    break
            
            if not is_overlapping:
                filtered_matches.append(match)
        
        return filtered_matches[:10]  # Limit to top 10 matches per sequence
    
    def _calculate_mismatches(self, seq1: str, seq2: str) -> int:
        """Calculate number of mismatches between two sequences"""
        if len(seq1) != len(seq2):
            return max(len(seq1), len(seq2))
        
        mismatches = 0
        for i in range(len(seq1)):
            if seq1[i] != seq2[i]:
                mismatches += 1
        
        return mismatches

def create_conservation_map(sequences: Dict[str, str], results_df: pd.DataFrame) -> go.Figure:
    """Create a clear graphical conservation map with color-coded regions"""
    
    # Create conservation regions data
    genome_length = min(len(seq) for seq in sequences.values())
    
    # Define conservation levels with clear colors
    conservation_levels = []
    colors = []
    
    for _, row in results_df.iterrows():
        start = row['start']
        end = row['end']
        conservation = row['conservation_score']
        
        # Color code based on conservation level
        if conservation >= 0.9:
            color = 'darkgreen'
            level = 'Highly Conserved (90%+)'
        elif conservation >= 0.8:
            color = 'green'
            level = 'Conserved (80-90%)'
        elif conservation >= 0.6:
            color = 'yellow'
            level = 'Moderately Conserved (60-80%)'
        elif conservation >= 0.4:
            color = 'orange'
            level = 'Poorly Conserved (40-60%)'
        else:
            color = 'red'
            level = 'Variable (<40%)'
        
        conservation_levels.append({
            'start': start,
            'end': end,
            'conservation': conservation,
            'level': level,
            'color': color
        })
    
    # Create the conservation map visualization
    fig = go.Figure()
    
    # Add conservation regions as colored bars
    for region in conservation_levels:
        fig.add_trace(go.Scatter(
            x=[region['start'], region['end']],
            y=[1, 1],
            mode='lines',
            line=dict(color=region['color'], width=8),
            name=region['level'],
            showlegend=False,
            hovertemplate=f"Position: {region['start']}-{region['end']}<br>" +
                         f"Conservation: {region['conservation']:.1%}<br>" +
                         f"Level: {region['level']}<extra></extra>"
        ))
    
    # Add conservation score as a line plot above
    fig.add_trace(go.Scatter(
        x=results_df['start'],
        y=results_df['conservation_score'] + 1.5,
        mode='lines',
        name='Conservation Score',
        line=dict(color='purple', width=2),
        hovertemplate='Position: %{x}<br>Conservation: %{customdata:.1%}<extra></extra>',
        customdata=results_df['conservation_score']
    ))
    
    # Add legend manually
    legend_traces = [
        go.Scatter(x=[None], y=[None], mode='lines', line=dict(color='darkgreen', width=8), name='Highly Conserved (90%+)'),
        go.Scatter(x=[None], y=[None], mode='lines', line=dict(color='green', width=8), name='Conserved (80-90%)'),
        go.Scatter(x=[None], y=[None], mode='lines', line=dict(color='yellow', width=8), name='Moderately Conserved (60-80%)'),
        go.Scatter(x=[None], y=[None], mode='lines', line=dict(color='orange', width=8), name='Poorly Conserved (40-60%)'),
        go.Scatter(x=[None], y=[None], mode='lines', line=dict(color='red', width=8), name='Variable (<40%)')
    ]
    
    for trace in legend_traces:
        fig.add_trace(trace)
    
    fig.update_layout(
        title="Genome Conservation Map - Color-Coded by Conservation Level (Hover for details)",
        xaxis_title="Genomic Position (bp)",
        yaxis_title="Conservation Level",
        height=400,
        yaxis=dict(
            tickvals=[1, 2],
            ticktext=['Conservation Regions', 'Conservation Score'],
            range=[0.5, 2.5]
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="center",
            x=0.5
        )
    )
    
    return fig

def create_sequence_alignment_view(sequences: Dict[str, str], results_df: pd.DataFrame) -> go.Figure:
    """Create a sequence alignment view showing all sequences with conservation coloring"""
    
    seq_ids = list(sequences.keys())
    seq_list = list(sequences.values())
    
    # Limit to reasonable display length
    display_length = min(1000, min(len(seq) for seq in seq_list))
    
    fig = go.Figure()
    
    # Color map for bases
    base_colors = {'A': 'red', 'T': 'blue', 'G': 'green', 'C': 'orange', 'U': 'blue', 'N': 'gray', '-': 'white'}
    
    # Create a heatmap-style visualization
    for i, (seq_id, sequence) in enumerate(sequences.items()):
        # Get conservation scores for positions
        conservation_scores = []
        for pos in range(1, display_length + 1):
            # Find conservation score for this position
            conservation = 0
            for _, row in results_df.iterrows():
                if row['start'] <= pos <= row['end']:
                    conservation = row['conservation_score']
                    break
            conservation_scores.append(conservation)
        
        # Create sequence display with conservation background
        fig.add_trace(go.Scatter(
            x=list(range(1, display_length + 1)),
            y=[i] * display_length,
            mode='markers',
            marker=dict(
                size=8,
                color=conservation_scores,
                colorscale='RdYlGn',  # Red-Yellow-Green scale
                cmin=0,
                cmax=1,
                colorbar=dict(title="Conservation Score") if i == 0 else None,
                line=dict(width=0.5, color='black')
            ),
            text=[sequence[j] for j in range(display_length)],
            textposition="middle center",
            textfont=dict(size=8, color='black'),
            name=seq_id,
            hovertemplate=f'{seq_id}<br>Position: %{{x}}<br>Base: %{{text}}<br>Conservation: %{{marker.color:.2f}}<extra></extra>'
        ))
    
    fig.update_layout(
        title=f"Multi-Sequence Alignment View (First {display_length} bp)",
        xaxis_title="Genomic Position (bp)",
        yaxis_title="Sequences",
        yaxis=dict(
            tickvals=list(range(len(seq_ids))),
            ticktext=seq_ids,
            autorange='reversed'
        ),
        height=max(300, len(seq_ids) * 40),
        showlegend=False
    )
    
    return fig

def create_detailed_conservation_heatmap(sequences: Dict[str, str], results_df: pd.DataFrame) -> go.Figure:
    """Create a detailed heatmap showing conservation at each position"""
    
    seq_ids = list(sequences.keys())
    seq_list = list(sequences.values())
    
    # Calculate conservation at each position
    min_length = min(len(seq) for seq in seq_list)
    display_length = min(500, min_length)  # Limit for performance
    
    position_conservation = []
    
    for pos in range(display_length):
        # Get bases at this position from all sequences
        bases_at_pos = [seq[pos] for seq in seq_list]
        
        # Calculate conservation (fraction of sequences with most common base)
        base_counts = {}
        for base in bases_at_pos:
            base_counts[base] = base_counts.get(base, 0) + 1
        
        max_count = max(base_counts.values())
        conservation = max_count / len(bases_at_pos)
        position_conservation.append(conservation)
    
    # Create heatmap data
    heatmap_data = []
    for seq in seq_list:
        row = []
        for pos in range(display_length):
            # Use conservation score as the value, but color by base type
            base = seq[pos]
            conservation = position_conservation[pos]
            
            # Encode base type and conservation
            if base == 'A':
                value = conservation + 0.0
            elif base in ['T', 'U']:
                value = conservation + 1.0
            elif base == 'G':
                value = conservation + 2.0
            elif base == 'C':
                value = conservation + 3.0
            else:
                value = conservation + 4.0  # N or gap
            
            row.append(value)
        heatmap_data.append(row)
    
    # Create custom colorscale
    colorscale = [
        [0.0, 'lightcoral'],    # Low conservation A
        [0.2, 'red'],          # High conservation A
        [0.25, 'lightblue'],   # Low conservation T/U
        [0.45, 'blue'],        # High conservation T/U
        [0.5, 'lightgreen'],   # Low conservation G
        [0.7, 'green'],        # High conservation G
        [0.75, 'lightyellow'], # Low conservation C
        [0.95, 'orange'],      # High conservation C
        [1.0, 'gray']          # N/gaps
    ]
    
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_data,
        x=list(range(1, display_length + 1)),
        y=seq_ids,
        colorscale=colorscale,
        showscale=True,
        colorbar=dict(
            title="Base Type & Conservation",
            tickvals=[0.1, 1.1, 2.1, 3.1, 4.1],
            ticktext=["A", "T/U", "G", "C", "N/Gap"]
        ),
        hovertemplate='Sequence: %{y}<br>Position: %{x}<br>Conservation: %{customdata:.2f}<extra></extra>',
        customdata=[[position_conservation[j] for j in range(display_length)] for _ in range(len(seq_list))]
    ))
    
    fig.update_layout(
        title=f"Conservation Heatmap by Position (First {display_length} bp)",
        xaxis_title="Genomic Position (bp)",
        yaxis_title="Sequences",
        height=max(400, len(seq_ids) * 30)
    )
    
    return fig

def create_alignment_heatmap(sequences: Dict[str, str], results_df: pd.DataFrame) -> go.Figure:
    """Create a heatmap showing sequence alignment across all positions"""
    seq_ids = list(sequences.keys())
    seq_list = list(sequences.values())
    
    # Get the length to analyze (use shortest sequence)
    min_length = min(len(seq) for seq in seq_list)
    analysis_length = min(min_length, 1000)  # Limit for visualization performance
    
    # Create base identity matrix
    base_to_num = {'A': 1, 'T': 2, 'G': 3, 'C': 4, 'U': 2, 'N': 0, '-': 0}
    
    # Build matrix for heatmap
    matrix = []
    for seq_id, sequence in sequences.items():
        row = [base_to_num.get(base, 0) for base in sequence[:analysis_length]]
        matrix.append(row)
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=matrix,
        x=list(range(1, analysis_length + 1)),
        y=seq_ids,
        colorscale=[
            [0, 'white'],      # Gaps/Unknown
            [0.25, 'red'],     # A
            [0.5, 'blue'],     # T/U  
            [0.75, 'green'],   # G
            [1, 'orange']      # C
        ],
        showscale=True,
        colorbar=dict(
            title="Base Type",
            tickvals=[0.5, 1.5, 2.5, 3.5],
            ticktext=["Gap/N", "A", "T/U", "G", "C"]
        )
    ))
    
    fig.update_layout(
        title="Sequence Alignment Heatmap (Colored by Base Type)",
        xaxis_title="Genomic Position (bp)",
        yaxis_title="Sequences",
        height=max(300, len(seq_ids) * 50),
        xaxis=dict(showgrid=True, gridwidth=1, gridcolor='lightgray'),
        yaxis=dict(showgrid=True, gridwidth=1, gridcolor='lightgray')
    )
    
    return fig

def create_conservation_track(results_df: pd.DataFrame, sequences: Dict[str, str]) -> go.Figure:
    """Create a detailed conservation track showing conservation at each position"""
    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        subplot_titles=('Conservation Score by Position', 'Base Composition at Each Position'),
        vertical_spacing=0.1,
        row_heights=[0.6, 0.4]
    )
    
    # Top plot: Conservation score
    fig.add_trace(
        go.Scatter(
            x=results_df['start'], 
            y=results_df['conservation_score'],
            mode='lines+markers',
            name='Conservation Score',
            line=dict(color='purple', width=2),
            marker=dict(size=4),
            hovertemplate='Position: %{x}<br>Conservation: %{y:.2f}<extra></extra>'
        ),
        row=1, col=1
    )
    
    # Add conservation threshold line
    fig.add_hline(y=0.8, line_dash="dash", line_color="red", 
                  annotation_text="High Conservation (80%)", row=1, col=1)
    
    # Bottom plot: Base composition diversity
    # Calculate base diversity for each window
    diversity_scores = []
    for _, row in results_df.iterrows():
        windows = row.get('individual_windows', [])
        if windows:
            # Calculate Shannon diversity for each position in the window
            window_diversity = []
            for pos in range(len(windows[0])):
                bases = [window[pos] for window in windows if pos < len(window)]
                base_counts = {}
                for base in bases:
                    base_counts[base] = base_counts.get(base, 0) + 1
                
                # Shannon diversity
                total = sum(base_counts.values())
                diversity = 0
                for count in base_counts.values():
                    if count > 0:
                        p = count / total
                        diversity -= p * np.log2(p)
                window_diversity.append(diversity)
            
            avg_diversity = np.mean(window_diversity) if window_diversity else 0
        else:
            avg_diversity = 0
        
        diversity_scores.append(avg_diversity)
    
    fig.add_trace(
        go.Scatter(
            x=results_df['start'],
            y=diversity_scores,
            mode='lines',
            name='Base Diversity',
            line=dict(color='green', width=2),
            hovertemplate='Position: %{x}<br>Diversity: %{y:.2f}<extra></extra>'
        ),
        row=2, col=1
    )
    
    fig.update_layout(
        title="Conservation Track Analysis",
        height=600,
        showlegend=True
    )
    
    fig.update_xaxes(title_text="Genomic Position (bp)", row=2, col=1)
    fig.update_yaxes(title_text="Conservation Score", row=1, col=1)
    fig.update_yaxes(title_text="Base Diversity", row=2, col=1)
    
    return fig

def create_sequence_match_visualization(matches_df: pd.DataFrame, custom_sequence: str) -> go.Figure:
    """Create visualization showing where custom sequence matches occur in genomes"""
    if matches_df.empty:
        return go.Figure()
    
    fig = go.Figure()
    
    # Color code by identity percentage
    colors = []
    for identity in matches_df['identity_percentage']:
        if identity >= 95:
            colors.append('darkgreen')
        elif identity >= 90:
            colors.append('green')
        elif identity >= 80:
            colors.append('yellow')
        elif identity >= 70:
            colors.append('orange')
        else:
            colors.append('red')
    
    # Create scatter plot
    fig.add_trace(go.Scatter(
        x=matches_df['match_start'],
        y=matches_df['sequence_id'],
        mode='markers',
        marker=dict(
            size=matches_df['match_length'] * 2,  # Size based on match length
            color=colors,
            line=dict(width=1, color='black'),
            opacity=0.7
        ),
        text=[f"Pos: {row['match_start']}-{row['match_end']}<br>"
              f"Length: {row['match_length']} bp<br>"
              f"Identity: {row['identity_percentage']:.1f}%<br>"
              f"Mismatches: {row['mismatches']}" 
              for _, row in matches_df.iterrows()],
        hovertemplate='%{text}<extra></extra>',
        name='Sequence Matches'
    ))
    
    fig.update_layout(
        title=f"Custom Sequence Matches in Loaded Genomes<br><sub>Query: {custom_sequence[:50]}{'...' if len(custom_sequence) > 50 else ''}</sub>",
        xaxis_title="Genomic Position (bp)",
        yaxis_title="Genome Sequences",
        height=max(400, len(matches_df['sequence_id'].unique()) * 50),
        showlegend=False
    )
    
    return fig

def create_sequence_alignment_display(matches_df: pd.DataFrame, custom_sequence: str) -> None:
    """Display detailed sequence alignments for matches"""
    if matches_df.empty:
        st.warning("No matches found for the custom sequence.")
        return
    
    st.subheader("Detailed Sequence Alignments")
    
    # Group matches by sequence
    for seq_id in matches_df['sequence_id'].unique():
        seq_matches = matches_df[matches_df['sequence_id'] == seq_id].sort_values('identity_percentage', ascending=False)
        
        with st.expander(f"Matches in {seq_id} ({len(seq_matches)} matches)", expanded=len(seq_matches) <= 3):
            for i, (_, match) in enumerate(seq_matches.iterrows()):
                st.write(f"**Match {i+1}:** Position {match['match_start']}-{match['match_end']} "
                        f"(Identity: {match['identity_percentage']:.1f}%, Length: {match['match_length']} bp)")
                
                # Show alignment
                query_seq = custom_sequence[:match['match_length']]
                target_seq = match['match_sequence']
                
                # Create alignment display
                alignment_display = ""
                for j, (q, t) in enumerate(zip(query_seq, target_seq)):
                    if q == t:
                        alignment_display += "|"
                    else:
                        alignment_display += " "
                
                col1, col2, col3 = st.columns([1, 1, 1])
                with col1:
                    st.write("**Query:**")
                    st.code(query_seq)
                with col2:
                    st.write("**Alignment:**")
                    st.code(alignment_display)
                with col3:
                    st.write("**Target:**")
                    st.code(target_seq)
                
                # Show context
                st.write("**Context (±50 bp):**")
                context = match['context_sequence']
                context_start = match['context_start']
                context_end = match['context_end']
                
                # Highlight the match in context
                match_rel_start = match['match_start'] - context_start - 1
                match_rel_end = match['match_end'] - context_start
                
                highlighted_context = (context[:match_rel_start] + 
                                    "[" + context[match_rel_start:match_rel_end] + "]" + 
                                    context[match_rel_end:])
                
                st.code(f"Position {context_start+1}-{context_end}: {highlighted_context}")
                st.write("---")

def create_sequence_browser(sequences: Dict[str, str], results_df: pd.DataFrame):
    """Create an interactive sequence browser"""
    seq_ids = list(sequences.keys())
    
    # Position selector
    max_position = min(len(seq) for seq in sequences.values())
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        start_pos = st.number_input(
            "Start position:", 
            min_value=1, 
            max_value=max_position-50, 
            value=1,
            help="Select the starting position (1-based) for the sequence view. This determines where the sequence browser begins displaying sequences."
        )
    
    with col2:
        window_length = st.selectbox(
            "View window size:",
            [50, 100, 200, 500],
            index=0,
            help="Number of bases to display in the sequence browser. Larger windows show more context but may be harder to read. Smaller windows focus on specific regions."
        )
    
    end_pos = min(start_pos + window_length - 1, max_position)
    
    # Extract sequences for the selected region
    st.write(f"**Showing positions {start_pos}-{end_pos}:**")
    st.info("💡 **Sequence Browser Help**: Use the controls above to navigate through sequences. Green circles (🟢) indicate high conservation, yellow (🟡) moderate conservation, and red (🔴) low conservation. The position ruler helps you locate specific bases.")
    
    # Create alignment display
    alignment_data = []
    for seq_id, sequence in sequences.items():
        seq_segment = sequence[start_pos-1:end_pos]
        alignment_data.append({
            'Sequence ID': seq_id,
            'Sequence': seq_segment,
            'Length': len(seq_segment)
        })
    
    # Display as formatted text with position markers
    st.write("**Position markers:**")
    
    # Create position ruler
    ruler_top = ""
    ruler_bottom = ""
    for i in range(len(alignment_data[0]['Sequence'])):
        pos = start_pos + i
        if pos % 10 == 0:
            ruler_top += str(pos)[-2] if pos >= 10 else " "
            ruler_bottom += str(pos)[-1]
        else:
            ruler_top += " "
            ruler_bottom += str(pos)[-1] if pos % 5 == 0 else "."
    
    st.code(f"     {ruler_top}\n     {ruler_bottom}")
    
    # Display sequences with conservation highlighting
    for data in alignment_data:
        seq_id = data['Sequence ID']
        sequence = data['Sequence']
        
        # Color code based on conservation
        colored_sequence = ""
        for i, base in enumerate(sequence):
            pos = start_pos + i
            
            # Find conservation score for this position
            conservation_score = 0
            for _, row in results_df.iterrows():
                if row['start'] <= pos <= row['end']:
                    conservation_score = row['conservation_score']
                    break
            
            # Color based on conservation level
            if conservation_score >= 0.8:
                colored_sequence += f"🟢{base}"  # High conservation
            elif conservation_score >= 0.6:
                colored_sequence += f"🟡{base}"  # Medium conservation
            else:
                colored_sequence += f"🔴{base}"  # Low conservation
        
        st.write(f"**{seq_id[:20]}:** {colored_sequence}")
    
    # Legend
    st.write("**Legend:** 🟢 High conservation (≥80%) | 🟡 Medium conservation (60-80%) | 🔴 Low conservation (<60%)")
    
    # Show conservation statistics for this region
    region_stats = results_df[
        (results_df['start'] >= start_pos) & 
        (results_df['end'] <= end_pos)
    ]
    
    if not region_stats.empty:
        avg_conservation = region_stats['conservation_score'].mean()
        max_conservation = region_stats['conservation_score'].max()
        min_conservation = region_stats['conservation_score'].min()
        
        st.write(f"**Region Statistics (positions {start_pos}-{end_pos}):**")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Average Conservation", f"{avg_conservation:.2f}")
        with col2:
            st.metric("Maximum Conservation", f"{max_conservation:.2f}")
        with col3:
            st.metric("Minimum Conservation", f"{min_conservation:.2f}")

def test_ncbi_connection(email: str) -> bool:
    """Test NCBI connection"""
    try:
        Entrez.email = email
        handle = Entrez.esearch(db="nucleotide", term="Escherichia coli", retmax=1)
        test_results = Entrez.read(handle)
        handle.close()
        return bool(test_results.get('IdList'))
    except Exception as e:
        st.error(f"NCBI connection failed: {str(e)}")
        return False

def main():
    st.title("Genome Conservation Analysis Tool")
    st.markdown("Identifies conserved regions by comparing multiple genome sequences from NCBI.")
    
    # Add helpful information about the tool
    with st.expander("ℹ️ How this tool works", expanded=False):
        st.markdown("""
        **This tool helps you find conserved regions in genomic sequences by:**
        
        1. **Searching NCBI**: Finds multiple genome sequences for your species of interest
        2. **Comparative Analysis**: Compares sequences using sliding window analysis
        3. **Conservation Scoring**: Calculates conservation scores for each region
        4. **Results Display**: Shows the most conserved regions with detailed information
        
        **Key Features:**
        - Works with any organism type (viral, bacterial, eukaryotic, viroid)
        - Auto-optimizes parameters based on sequence length and organism type
        - Shows consensus sequences and individual sequence alignments
        - Provides downloadable results for further analysis
        - Includes test mode for learning and troubleshooting
        
        **Tips for best results:**
        - Use the scientific name of your organism
        - Select the correct organism type for optimized parameters
        - Try different window sizes if conservation is low
        - Include more sequences for better statistical power
        """)
    
    # Sidebar configuration
    with st.sidebar:
        st.header("Configuration")
        
        # Email for NCBI
        email = st.text_input(
            "Email (required for NCBI):",
            placeholder="your.email@example.com",
            help="NCBI requires an email address for API access. This is used to identify your requests and is required by NCBI's usage policies. Your email is not stored or shared."
        )
        
        if not email or '@' not in email:
            st.warning("Please provide a valid email address")
            st.stop()
        
        # Species selection
        st.subheader("Species Selection")
        example_species = st.selectbox(
            "Choose example or enter custom:",
            ["Custom", "Homo sapiens", "Escherichia coli", "Saccharomyces cerevisiae", 
             "Hop latent viroid", "Potato spindle tuber viroid", "Tobacco mosaic virus"],
            help="Select a pre-configured species or choose 'Custom' to enter your own. Pre-configured species have optimized search strategies for better results."
        )
        
        if example_species == "Custom":
            species = st.text_input(
                "Species name:",
                placeholder="Enter scientific name",
                help="Enter the scientific name using standard nomenclature (e.g., 'Hop latent viroid', 'Escherichia coli', 'Homo sapiens'). Use the full scientific name for best results."
            )
        else:
            species = example_species
        
        if not species:
            st.warning("Please enter a species name")
            st.stop()
        
        # Test mode option
        st.subheader("Testing & Debugging")
        test_mode = st.checkbox("Use synthetic test data", 
                               help="Generate synthetic test sequences with known conservation patterns. Useful for testing the tool and understanding how conservation analysis works. No NCBI access required.")
        
        if test_mode:
            st.info("Test mode: Will generate synthetic sequences with predefined conservation patterns")
        
        # Organism type selection (removed auto-detect)
        organism_type = st.selectbox(
            "Organism type:",
            ["Viroid", "Virus", "Bacteria", "Fungi", "Plant", "Animal"],
            help="Select the type of organism you're analyzing. This affects search strategies, sequence length expectations, and analysis parameters. Viroids are very small (200-400 bp), viruses are small to medium (1-500 kb), while bacteria and eukaryotes are much larger."
        )
        
        # Auto-adjust parameters based on organism type AND sequence data
        st.subheader("Analysis Parameters")
        
        if 'sequences' in st.session_state and st.session_state['sequences']:
            # Get actual sequence lengths for better parameter adjustment
            seq_lengths = [seq['length'] for seq in st.session_state['sequences'] if seq['length'] > 0]
            if seq_lengths:
                avg_length = sum(seq_lengths) / len(seq_lengths)
                max_length = max(seq_lengths)
                min_length = min(seq_lengths)
                
                st.info(f"Detected sequences: {min_length}-{max_length} bp (avg: {avg_length:.0f} bp)")
                
                # Adjust parameters based on actual sequence data
                if organism_type == "Viroid" or avg_length <= 500:
                    default_window = max(10, min_length // 10)
                    default_step = max(2, default_window // 5)
                    default_max_sequences = min(20, len(st.session_state['sequences']))
                    param_info = f"Viroid Mode: {default_window}bp windows, {default_step}bp steps"
                    
                elif organism_type == "Virus" or avg_length <= 50000:
                    default_window = max(50, int(avg_length * 0.01))
                    default_step = default_window // 2
                    default_max_sequences = 15
                    param_info = f"Virus Mode: {default_window}bp windows, {default_step}bp steps"
                    
                else:
                    default_window = 1000
                    default_step = 500  
                    default_max_sequences = 10
                    param_info = f"Large genome mode: {default_window}bp windows"
            else:
                if organism_type == "Viroid":
                    default_window = 25
                    default_step = 5
                    default_max_sequences = 20
                    param_info = "Viroid Mode: 25bp windows, 5bp steps"
                else:
                    default_window = 1000
                    default_step = 500
                    default_max_sequences = 10  
                    param_info = "Standard mode: 1000bp windows"
        else:
            if organism_type == "Viroid":
                default_window = 25
                default_step = 5
                default_max_sequences = 20
                param_info = "Viroid Mode: 25bp windows, 5bp steps"
            elif organism_type == "Virus":
                default_window = 200
                default_step = 50
                default_max_sequences = 15
                param_info = "Virus Mode: 200bp windows, 50bp steps"
            else:
                default_window = 1000
                default_step = 500
                default_max_sequences = 10
                param_info = "Standard Mode: 1000bp windows, 500bp steps"
        
        st.success(f"Auto-optimized: {param_info}")
        
        # Advanced parameters
        with st.expander("Advanced: Custom Parameters"):
            custom_params = st.checkbox("Override auto-optimized parameters", 
                                      value=False,
                                      help="Check this to manually adjust analysis parameters instead of using the auto-optimized settings")
            
            if custom_params:
                window_size = st.slider("Window size (bp):", 5, 2000, default_window,
                                      help="Size of the sliding window for conservation analysis. Smaller windows detect shorter conserved regions but may be noisy. Larger windows smooth out noise but may miss short conserved elements.")
                step_size = st.slider("Step size (bp):", 1, 500, default_step,
                                    help="How many base pairs to move the window for each analysis step. Smaller steps give higher resolution but take longer. Larger steps are faster but may miss conserved regions between steps.")
                max_sequences = st.slider("Max sequences:", 2, 50, default_max_sequences,
                                        help="Maximum number of sequences to analyze. More sequences give better conservation statistics but take longer to process. NCBI has rate limits, so very high numbers may cause delays.")
            else:
                window_size = default_window
                step_size = default_step
                max_sequences = default_max_sequences
    
    # Main interface
    col1, col2 = st.columns([1, 1])
    
    with col1:
        if test_mode:
            if st.button("Generate Test Sequences", type="primary", 
                        help="Generate synthetic test sequences with known conservation patterns. This is useful for testing the tool and understanding how conservation analysis works."):
                analyzer = NCBIGenomeAnalyzer(email, organism_type)
                
                with st.spinner("Generating test sequences..."):
                    test_sequences = analyzer.generate_test_sequences(species, max_sequences)
                    
                    # Store as if they were real sequences
                    test_seq_list = []
                    for seq_id, sequence in test_sequences.items():
                        test_seq_list.append({
                            'sequence_id': seq_id,
                            'title': f'Synthetic test sequence {seq_id}',
                            'organism': species,
                            'length': len(sequence),
                            'source': 'synthetic'
                        })
                    
                    st.session_state['sequences'] = test_seq_list
                    st.session_state['test_sequences'] = test_sequences
                    st.session_state['species'] = species
                    st.session_state['organism_type'] = organism_type
                    
                st.success("Test sequences generated! You can now run analysis to verify the tool works.")
        else:
            if st.button("Search Genome Sequences", type="primary",
                        help="Search NCBI's nucleotide database for genome sequences of the specified species. This will find and retrieve sequences for comparative analysis."):
                analyzer = NCBIGenomeAnalyzer(email, organism_type)
                
                with st.spinner(f"Searching for {species} sequences..."):
                    sequences = analyzer.search_genomes(species, max_sequences)
                
                if sequences:
                    st.session_state['sequences'] = sequences
                    st.session_state['species'] = species
                    st.session_state['organism_type'] = organism_type
                    st.success(f"Found {len(sequences)} sequences for {species}")
                else:
                    st.error(f"No sequences found for {species}")
    
    with col2:
        if st.button("Test NCBI Connection",
                    help="Test your connection to NCBI's database to ensure you can retrieve sequences. This helps diagnose connection issues before running the full analysis."):
            if test_ncbi_connection(email):
                st.success("NCBI connection successful!")
            else:
                st.error("NCBI connection failed")
    
    # Display sequences and analysis interface
    if 'sequences' in st.session_state:
        st.subheader(f"Found Sequences for {st.session_state['species']}")
        
        sequence_data = []
        for seq in st.session_state['sequences']:
            sequence_data.append({
                'Sequence ID': seq['sequence_id'],
                'Title': seq['title'][:60] + '...' if len(seq['title']) > 60 else seq['title'],
                'Organism': seq['organism'],
                'Length (bp)': f"{seq['length']:,}"
            })
        
        df_sequences = pd.DataFrame(sequence_data)
        st.dataframe(df_sequences, use_container_width=True)
        
        current_organism_type = st.session_state.get('organism_type', organism_type)
        
        # ALWAYS use comparative analysis for ALL organisms
        st.success(f"{current_organism_type.upper()} MODE: Comparative Conservation Analysis")
        st.info("All available sequences will be compared to find conserved regions across multiple genomes/isolates")
        st.info(f"Using optimized parameters: {window_size}bp windows, {step_size}bp steps")
        
        # Show how many sequences will be analyzed
        num_available = len(st.session_state['sequences'])
        num_to_analyze = min(num_available, max_sequences)
        
        st.write(f"**Available sequences:** {num_available}")
        st.write(f"**Will analyze:** {num_to_analyze} (limited by max_sequences parameter)")
        
        if st.button("Analyze ALL Sequences for Conservation", type="primary",
                    help="Perform comparative conservation analysis on all available sequences. This will download sequences from NCBI (if not using test data) and identify the most conserved regions across all sequences."):
            analyzer = NCBIGenomeAnalyzer(email, current_organism_type)
            sequence_ids = [seq['sequence_id'] for seq in st.session_state['sequences']]
            
            with st.spinner(f"Performing comparative conservation analysis on {num_to_analyze} sequences..."):
                try:
                    # Use test sequences if available, otherwise fetch from NCBI
                    if 'test_sequences' in st.session_state:
                        sequences = st.session_state['test_sequences']
                        st.info("Using synthetic test sequences for analysis")
                    else:
                        # Fetch sequences
                        sequences = analyzer.fetch_all_sequences_simultaneously(sequence_ids[:num_to_analyze])
                    
                    if len(sequences) < 1:
                        st.error("Could not fetch any sequences")
                        return
                    
                    if len(sequences) == 1:
                        st.warning("Only 1 sequence available - showing composition analysis")
                        st.info("Note: True conservation analysis requires multiple sequences for comparison")
                        
                        # Single sequence composition analysis
                        seq_id = list(sequences.keys())[0]
                        sequence = sequences[seq_id]
                        
                        results = []
                        for pos in range(0, len(sequence) - window_size + 1, step_size):
                            window_seq = sequence[pos:pos + window_size]
                            gc_content = (window_seq.count('G') + window_seq.count('C')) / len(window_seq) * 100
                            
                            results.append({
                                'start': pos + 1,
                                'end': pos + window_size,
                                'gc_content': gc_content,
                                'sequence': window_seq
                            })
                        
                        df_results = pd.DataFrame(results)
                        st.session_state['single_results'] = df_results
                        st.session_state['analyzed_sequence'] = sequence
                        st.rerun()
                    
                    else:
                        # Multi-sequence comparative analysis
                        st.success(f"Fetched {len(sequences)} sequences - performing comparative analysis")
                        
                        df_analysis = analyzer.true_comparative_analysis(sequences, window_size, step_size)
                        
                        if not df_analysis.empty:
                            st.session_state['comparative_results'] = df_analysis
                            st.session_state['compared_sequences'] = sequences
                            st.session_state['num_sequences_compared'] = len(sequences)
                            st.success(f"Comparative analysis complete! {len(df_analysis)} windows analyzed")
                            st.rerun()
                        else:
                            st.error("Analysis produced no results")
                            
                except Exception as e:
                    st.error(f"Analysis failed: {e}")
                    import traceback
                    st.code(traceback.format_exc())
    
    # Display single sequence results
    if 'single_results' in st.session_state:
        st.header("Single Viroid Sequence Analysis")
        
        df_results = st.session_state['single_results']
        sequence = st.session_state['analyzed_sequence']
        
        st.info(f"Analyzed sequence: {len(sequence)} bp")
        
        fig = px.line(df_results, x='start', y='gc_content', title='GC Content Along Sequence')
        st.plotly_chart(fig, use_container_width=True)
        
        st.subheader("Sequence Windows")
        display_results = df_results.head(10)
        st.dataframe(display_results)
        
        with st.expander("Full Sequence"):
            st.code(sequence)
    
    # Display comparative results
    if 'comparative_results' in st.session_state:
        st.header("Comparative Conservation Analysis Results")
        
        df_results = st.session_state['comparative_results']
        compared_sequences = st.session_state['compared_sequences']
        num_sequences = st.session_state['num_sequences_compared']
        
        # Summary statistics
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Sequences Compared", num_sequences)
        with col2:
            st.metric("Windows Analyzed", len(df_results))
        with col3:
            # Show the top conservation score instead of counting regions above threshold
            max_conservation = df_results['conservation_score'].max()
            st.metric("Highest Conservation Score", f"{max_conservation:.1%}")
        with col4:
            avg_identity = df_results['identity_percentage'].mean()
            st.metric("Average Identity", f"{avg_identity:.1f}%")
        
        # Conservation visualization
        st.write("**Conservation Analysis Charts** - Hover over points for detailed information")
        fig = make_subplots(
            rows=3, cols=1,
            shared_xaxes=True,
            subplot_titles=(
                'Conservation Score (fraction conserved) - Higher values indicate more conserved regions', 
                'Identity Percentage (all identical) - Percentage of positions where sequences are identical', 
                'GC Content of Consensus - GC content of the consensus sequence at each position'
            ),
            vertical_spacing=0.08
        )
        
        fig.add_trace(
            go.Scatter(x=df_results['start'], y=df_results['conservation_score'], 
                      mode='lines', name='Conservation Score', line=dict(color='purple', width=2)),
            row=1, col=1
        )
        
        fig.add_trace(
            go.Scatter(x=df_results['start'], y=df_results['identity_percentage'], 
                      mode='lines', name='Identity %', line=dict(color='red', width=2)),
            row=2, col=1
        )
        
        fig.add_trace(
            go.Scatter(x=df_results['start'], y=df_results['gc_content'], 
                      mode='lines', name='GC Content', line=dict(color='blue', width=2)),
            row=3, col=1
        )
        
        fig.update_layout(height=700, showlegend=False, 
                         title_text=f"Comparative Analysis: {st.session_state.get('species', 'Unknown')} ({num_sequences} sequences)")
        fig.update_xaxes(title_text="Genomic Position (bp)", row=3, col=1)
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Most conserved regions with sequences
        st.subheader("Most Conserved Regions with Sequences")
        
        # Always show the most conserved regions, regardless of absolute conservation level
        all_regions = df_results.copy()
        all_regions = all_regions.sort_values('conservation_score', ascending=False)
        
        # Show top regions (up to 20, or all if fewer than 20)
        num_top_regions = min(20, len(all_regions))
        top_conserved_regions = all_regions.head(num_top_regions)
        
        if not top_conserved_regions.empty:
            # Calculate statistics for context
            max_conservation = top_conserved_regions['conservation_score'].max()
            min_conservation = top_conserved_regions['conservation_score'].min()
            avg_conservation = top_conserved_regions['conservation_score'].mean()
            
            st.write(f"Showing top {len(top_conserved_regions)} most conserved regions")
            st.write(f"Conservation range: {min_conservation:.1%} - {max_conservation:.1%} (avg: {avg_conservation:.1%})")
            
            # Add context about conservation levels
            if max_conservation < 0.3:
                st.info("⚠️ Low overall conservation detected. These are the most conserved regions available, even though conservation is generally low.")
                st.write("💡 **Tip:** Consider using smaller window sizes or including more sequences to improve conservation detection.")
            elif max_conservation < 0.6:
                st.info("ℹ️ Moderate conservation detected. These regions show the highest conservation levels in your dataset.")
                st.write("💡 **Tip:** These regions may still be biologically significant despite moderate conservation scores.")
            else:
                st.success("✅ High conservation detected. These regions show strong conservation across sequences.")
            
            # Show conservation distribution
            st.write("**Conservation Score Distribution:**")
            conservation_stats = top_conserved_regions['conservation_score'].describe()
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Max", f"{conservation_stats['max']:.1%}")
            with col2:
                st.metric("75th percentile", f"{conservation_stats['75%']:.1%}")
            with col3:
                st.metric("Median", f"{conservation_stats['50%']:.1%}")
            with col4:
                st.metric("Min", f"{conservation_stats['min']:.1%}")
            
            for i, (idx, region) in enumerate(top_conserved_regions.iterrows()):
                with st.expander(f"Region {i+1}: Position {region['start']}-{region['end']} (Conservation: {region['conservation_score']:.1%})"):
                    col1, col2 = st.columns([1, 2])
                    
                    with col1:
                        st.metric("Conservation Score", f"{region['conservation_score']:.1%}")
                        st.metric("Identity Percentage", f"{region['identity_percentage']:.1f}%")
                        st.metric("GC Content", f"{region['gc_content']:.1f}%")
                        st.metric("Length", f"{region['length']} bp")
                    
                    with col2:
                        st.write("**Consensus Sequence:**")
                        consensus_seq = region['consensus_sequence']
                        
                        # Format sequence nicely (60 characters per line)
                        formatted_consensus = ""
                        for j in range(0, len(consensus_seq), 60):
                            line_num = j // 60 + 1
                            line_seq = consensus_seq[j:j+60]
                            formatted_consensus += f"{line_num:3d}: {line_seq}\n"
                        
                        st.code(formatted_consensus, language=None)
                        
                        # Show individual sequences for comparison
                        st.write("**Individual Sequences in this Region:**")
                        region_start = region['start'] - 1  # Convert to 0-based
                        region_end = region['end']
                        
                        for seq_id, full_sequence in compared_sequences.items():
                            if region_end <= len(full_sequence):
                                region_sequence = full_sequence[region_start:region_end]
                                st.code(f"{seq_id}: {region_sequence}", language=None)
            
            # Download option
            csv = top_conserved_regions.to_csv(index=False)
            st.download_button(
                label="Download Most Conserved Regions (CSV)",
                data=csv,
                file_name=f"most_conserved_regions_{st.session_state.get('species', 'unknown').replace(' ', '_')}.csv",
                mime="text/csv",
                help="Download the most conserved regions as a CSV file. This includes position coordinates, conservation scores, consensus sequences, and other metrics for further analysis in Excel or other tools."
            )
            
        else:
            st.warning("No regions found for analysis. This might indicate an issue with the sequence data or analysis parameters.")
    
    # Custom Sequence Comparison Section
    if 'compared_sequences' in st.session_state and st.session_state['compared_sequences']:
        st.header("🔍 Custom Sequence Comparison")
        st.write("Compare your own sequence against the loaded genome sequences to find matches and similarities.")
        
        # Custom sequence input
        col1, col2 = st.columns([2, 1])
        
        with col1:
            custom_sequence = st.text_area(
                "Enter your sequence:",
                placeholder="Enter DNA/RNA sequence (ATGC or AUGC)...",
                height=100,
                help="Enter the sequence you want to search for in the loaded genomes. Can be DNA (ATGC) or RNA (AUGC). Spaces and newlines will be automatically removed."
            )
        
        with col2:
            st.write("**Search Parameters:**")
            min_match_length = st.number_input(
                "Minimum match length:",
                min_value=10,
                max_value=1000,
                value=20,
                help="Minimum length of sequence matches to report"
            )
            max_mismatches = st.number_input(
                "Maximum mismatches:",
                min_value=0,
                max_value=10,
                value=2,
                help="Maximum number of mismatches allowed in a match"
            )
        
        if custom_sequence:
            # Clean the sequence
            clean_sequence = custom_sequence.upper().replace(' ', '').replace('\n', '').replace('\r', '')
            
            # Validate sequence
            valid_bases = set('ATGCU')
            if not all(base in valid_bases for base in clean_sequence):
                st.error("❌ Invalid sequence. Please use only A, T, G, C (DNA) or A, U, G, C (RNA) bases.")
            elif len(clean_sequence) < min_match_length:
                st.warning(f"⚠️ Sequence too short. Minimum length is {min_match_length} bases.")
            else:
                st.success(f"✅ Valid sequence: {len(clean_sequence)} bases")
                
                if st.button("🔍 Search for Sequence Matches", type="primary",
                           help="Search for your custom sequence in all loaded genome sequences"):
                    
                    with st.spinner("Searching for sequence matches..."):
                        analyzer = NCBIGenomeAnalyzer(email, st.session_state.get('organism_type', 'Bacteria'))
                        matches_df = analyzer.compare_custom_sequence(
                            clean_sequence, 
                            st.session_state['compared_sequences'],
                            min_match_length, 
                            max_mismatches
                        )
                        
                        if not matches_df.empty:
                            st.session_state['custom_sequence_matches'] = matches_df
                            st.session_state['custom_sequence'] = clean_sequence
                            st.success(f"Found {len(matches_df)} matches across {matches_df['sequence_id'].nunique()} sequences!")
                        else:
                            st.warning("No matches found. Try reducing the minimum match length or increasing the maximum mismatches allowed.")
                            st.session_state['custom_sequence_matches'] = pd.DataFrame()
                            st.session_state['custom_sequence'] = clean_sequence
    
    # Display custom sequence results
    if 'custom_sequence_matches' in st.session_state and not st.session_state['custom_sequence_matches'].empty:
        st.header("🎯 Custom Sequence Match Results")
        
        matches_df = st.session_state['custom_sequence_matches']
        custom_seq = st.session_state['custom_sequence']
        
        # Summary statistics
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Matches", len(matches_df))
        with col2:
            st.metric("Sequences with Matches", matches_df['sequence_id'].nunique())
        with col3:
            best_identity = matches_df['identity_percentage'].max()
            st.metric("Best Match Identity", f"{best_identity:.1f}%")
        with col4:
            avg_identity = matches_df['identity_percentage'].mean()
            st.metric("Average Identity", f"{avg_identity:.1f}%")
        
        # Match visualization
        st.subheader("Match Locations")
        fig = create_sequence_match_visualization(matches_df, custom_seq)
        st.plotly_chart(fig, use_container_width=True)
        
        # Detailed results table
        st.subheader("Match Details")
        
        # Filter options
        col1, col2 = st.columns(2)
        with col1:
            min_identity_filter = st.slider("Filter by minimum identity:", 0, 100, 0)
        with col2:
            selected_sequences = st.multiselect(
                "Filter by sequence:",
                options=matches_df['sequence_id'].unique(),
                default=matches_df['sequence_id'].unique()
            )
        
        # Apply filters
        filtered_matches = matches_df[
            (matches_df['identity_percentage'] >= min_identity_filter) &
            (matches_df['sequence_id'].isin(selected_sequences))
        ]
        
        if not filtered_matches.empty:
            # Display summary table
            display_cols = ['sequence_id', 'match_start', 'match_end', 'match_length', 
                           'identity_percentage', 'mismatches']
            st.dataframe(
                filtered_matches[display_cols].rename(columns={
                    'sequence_id': 'Sequence ID',
                    'match_start': 'Start Position',
                    'match_end': 'End Position', 
                    'match_length': 'Length (bp)',
                    'identity_percentage': 'Identity (%)',
                    'mismatches': 'Mismatches'
                }),
                use_container_width=True
            )
            
            # Detailed alignments
            create_sequence_alignment_display(filtered_matches, custom_seq)
            
            # Download option
            csv = filtered_matches.to_csv(index=False)
            st.download_button(
                label="📥 Download Match Results (CSV)",
                data=csv,
                file_name=f"sequence_matches_{st.session_state.get('species', 'unknown').replace(' ', '_')}.csv",
                mime="text/csv",
                help="Download all sequence match results as a CSV file for further analysis"
            )
        else:
            st.info("No matches meet the current filter criteria. Try adjusting the filters.")

if __name__ == "__main__":
    main()
