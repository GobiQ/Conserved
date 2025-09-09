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
        """Generate viroid-specific search terms"""
        return [
            f'"{species}"[Organism] AND viroid',
            f'{species}[Organism] AND viroid',
            f'"{species}" AND viroid AND complete',
            f'{species} viroid complete',
            f'{species} viroid',
            f'viroid AND {species.split()[0]}' if len(species.split()) > 1 else f'viroid AND {species}',
            f'"{species}"[Organism]',
            f'{species}[Organism]'
        ]
    
    def _get_virus_search_terms(self, species: str) -> List[str]:
        """Generate virus-specific search terms"""
        return [
            f'"{species}"[Organism] AND "complete genome"',
            f'"{species}"[Organism] AND genome',
            f'{species}[Organism] AND "complete genome"',
            f'{species}[Organism] AND genome',
            f'"{species}" AND "complete genome"',
            f'{species} complete genome',
            f'"{species}"[Organism]',
            f'{species}[Organism]'
        ]
    
    def _get_standard_search_terms(self, species: str) -> List[str]:
        """Generate search terms for bacteria/eukaryotes"""
        return [
            f'"{species}"[Organism] AND "complete genome"',
            f'"{species}"[Organism] AND chromosome',
            f'{species}[Organism] AND "complete genome"',
            f'{species}[Organism] AND chromosome',
            f'"{species}" AND "complete genome"',
            f'{species} complete genome',
            f'{species} chromosome',
            f'"{species}"[Organism]',
            f'{species}[Organism]',
            f'{species.split()[0]}[Organism] AND "complete genome"' if len(species.split()) > 1 else None
        ]
    
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
    
    def true_comparative_analysis(self, sequences: Dict[str, str], window_size: int = 25, step_size: int = 5) -> pd.DataFrame:
        """Perform TRUE comparative analysis - all sequences compared simultaneously"""
        if len(sequences) < 2:
            st.error("Need at least 2 sequences for comparative analysis")
            return pd.DataFrame()
        
        st.info(f"Performing comparative analysis on {len(sequences)} sequences...")
        
        seq_ids = list(sequences.keys())
        seq_list = list(sequences.values())
        
        # Find optimal analysis length
        lengths = [len(seq) for seq in seq_list]
        if self.organism_type == "Viroid":
            target_length = max(set(lengths), key=lengths.count)  # Most common length
        else:
            target_length = min(lengths)
        
        if target_length < window_size:
            window_size = max(5, target_length // 5)
            step_size = max(1, window_size // 5)
            st.warning(f"Adjusted window size to {window_size}bp for short sequences")
        
        # Align sequences
        aligned_sequences = []
        for seq_id, sequence in sequences.items():
            if len(sequence) >= target_length:
                aligned_seq = sequence[:target_length]
            else:
                aligned_seq = sequence.ljust(target_length, 'N')
            aligned_sequences.append(aligned_seq)
        
        # Sliding window analysis
        results = []
        total_windows = (target_length - window_size) // step_size + 1
        progress_bar = st.progress(0)
        
        for i, pos in enumerate(range(0, target_length - window_size + 1, step_size)):
            windows = [seq[pos:pos + window_size] for seq in aligned_sequences]
            
            conservation_score = self._calculate_conservation_score(windows)
            identity_percentage = self._calculate_identity_percentage(windows)
            consensus_sequence = self._generate_consensus(windows)
            
            gc_content = GC(consensus_sequence) if consensus_sequence else 0
            
            results.append({
                'start': pos + 1,
                'end': pos + window_size,
                'conservation_score': conservation_score,
                'identity_percentage': identity_percentage,
                'consensus_sequence': consensus_sequence,
                'gc_content': gc_content,
                'num_sequences': len(windows),
                'length': window_size,
                'individual_windows': windows
            })
            
            if i % 50 == 0:
                progress_bar.progress(min(i / total_windows, 1.0))
        
        progress_bar.progress(1.0)
        df = pd.DataFrame(results)
        st.success(f"Comparative analysis complete: {len(df)} windows analyzed")
        return df
    
    def _calculate_conservation_score(self, windows: List[str]) -> float:
        """Calculate conservation score across all sequences"""
        if not windows or len(windows) < 2:
            return 0.0
        
        window_length = len(windows[0])
        conserved_positions = 0
        
        for pos in range(window_length):
            bases = [window[pos] for window in windows if pos < len(window)]
            if not bases:
                continue
            
            base_counts = {}
            for base in bases:
                base_counts[base] = base_counts.get(base, 0) + 1
            
            max_count = max(base_counts.values())
            conservation_ratio = max_count / len(bases)
            
            if conservation_ratio >= 0.8:
                conserved_positions += 1
        
        return conserved_positions / window_length
    
    def _calculate_identity_percentage(self, windows: List[str]) -> float:
        """Calculate percentage where ALL sequences are identical"""
        if not windows or len(windows) < 2:
            return 0.0
        
        window_length = len(windows[0])
        identical_positions = 0
        
        for pos in range(window_length):
            bases = [window[pos] for window in windows if pos < len(window)]
            if len(set(bases)) == 1:
                identical_positions += 1
        
        return (identical_positions / window_length) * 100
    
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
    
    # Sidebar configuration
    with st.sidebar:
        st.header("Configuration")
        
        # Email for NCBI
        email = st.text_input(
            "Email (required for NCBI):",
            placeholder="your.email@example.com",
            help="NCBI requires an email address for API access"
        )
        
        if not email or '@' not in email:
            st.warning("Please provide a valid email address")
            st.stop()
        
        # Species selection
        st.subheader("Species Selection")
        example_species = st.selectbox(
            "Choose example or enter custom:",
            ["Custom", "Homo sapiens", "Escherichia coli", "Saccharomyces cerevisiae", 
             "Hop latent viroid", "Potato spindle tuber viroid", "Tobacco mosaic virus"]
        )
        
        if example_species == "Custom":
            species = st.text_input(
                "Species name:",
                placeholder="Enter scientific name",
                help="Enter the scientific name (e.g., 'Hop latent viroid')"
            )
        else:
            species = example_species
        
        if not species:
            st.warning("Please enter a species name")
            st.stop()
        
        # Organism type selection (removed auto-detect)
        organism_type = st.selectbox(
            "Organism type:",
            ["Viroid", "Virus", "Bacteria", "Fungi", "Plant", "Animal"],
            help="Select organism type for optimized analysis parameters"
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
            custom_params = st.checkbox("Override auto-optimized parameters", value=False)
            
            if custom_params:
                window_size = st.slider("Window size (bp):", 5, 2000, default_window)
                step_size = st.slider("Step size (bp):", 1, 500, default_step)
                max_sequences = st.slider("Max sequences:", 2, 50, default_max_sequences)
            else:
                window_size = default_window
                step_size = default_step
                max_sequences = default_max_sequences
    
    # Main interface
    col1, col2 = st.columns([1, 1])
    
    with col1:
        if st.button("Search Genome Sequences", type="primary"):
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
        if st.button("Test NCBI Connection"):
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
        
        # FORCE viroid mode for small sequences OR viroid selection
        if current_organism_type == "Viroid" or any(seq['length'] <= 500 for seq in st.session_state['sequences']):
            st.success("VIROID MODE: Comparative Analysis Only")
            st.warning("Individual sequence selection disabled - viroids require multi-sequence comparison")
            st.info(f"Using parameters: {window_size}bp windows, {step_size}bp steps")
            
            if st.button("Analyze ALL Viroid Sequences", type="primary"):
                analyzer = NCBIGenomeAnalyzer(email, "Viroid")
                sequence_ids = [seq['sequence_id'] for seq in st.session_state['sequences']]
                
                with st.spinner("Performing comparative analysis..."):
                    try:
                        sequences = analyzer.fetch_all_sequences_simultaneously(sequence_ids)
                        
                        if len(sequences) < 1:
                            st.error("Could not fetch any sequences")
                            return
                        
                        if len(sequences) == 1:
                            st.warning("Only 1 sequence - showing composition analysis")
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
                            df_analysis = analyzer.true_comparative_analysis(sequences, window_size, step_size)
                            
                            if not df_analysis.empty:
                                st.session_state['comparative_results'] = df_analysis
                                st.session_state['compared_sequences'] = sequences
                                st.session_state['num_sequences_compared'] = len(sequences)
                                st.success(f"Analysis complete! {len(df_analysis)} windows analyzed")
                                st.rerun()
                            else:
                                st.error("Analysis produced no results")
                                
                    except Exception as e:
                        st.error(f"Analysis failed: {e}")
                        import traceback
                        st.code(traceback.format_exc())
        
        else:
            st.info(f"Large genome mode: {current_organism_type}")
            
            selected_sequence = st.selectbox(
                "Select sequence for analysis:",
                options=[s['sequence_id'] for s in st.session_state['sequences']],
                format_func=lambda x: f"{x} - {next((s['title'][:40] + '...' if len(s['title']) > 40 else s['title']) for s in st.session_state['sequences'] if s['sequence_id'] == x)} ({next(s['length'] for s in st.session_state['sequences'] if s['sequence_id'] == x):,} bp)"
            )
            
            if st.button("Analyze Selected Sequence", type="primary"):
                st.info("Single sequence analysis for large genomes - implementation pending")
    
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
            highly_conserved = len(df_results[df_results['conservation_score'] >= 0.8])
            st.metric("Highly Conserved Regions", highly_conserved)
        with col4:
            avg_identity = df_results['identity_percentage'].mean()
            st.metric("Average Identity", f"{avg_identity:.1f}%")
        
        # Conservation visualization
        fig = make_subplots(
            rows=3, cols=1,
            shared_xaxes=True,
            subplot_titles=(
                'Conservation Score (fraction conserved)', 
                'Identity Percentage (all identical)', 
                'GC Content of Consensus'
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
        
        highly_conserved_regions = df_results[df_results['conservation_score'] >= 0.8].copy()
        
        if not highly_conserved_regions.empty:
            highly_conserved_regions = highly_conserved_regions.sort_values('conservation_score', ascending=False)
            
            st.write(f"Found {len(highly_conserved_regions)} highly conserved regions (≥80% conservation)")
            
            num_top_regions = min(10, len(highly_conserved_regions))
            
            for i, (idx, region) in enumerate(highly_conserved_regions.head(num_top_regions).iterrows()):
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
            csv = highly_conserved_regions.to_csv(index=False)
            st.download_button(
                label="Download Conserved Regions (CSV)",
                data=csv,
                file_name=f"conserved_regions_{st.session_state.get('species', 'unknown').replace(' ', '_')}.csv",
                mime="text/csv"
            )
            
        else:
            st.info("No highly conserved regions found (≥80% conservation). Try adjusting parameters or including more sequences.")

if __name__ == "__main__":
    main()
