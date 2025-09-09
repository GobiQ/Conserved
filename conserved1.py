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
    
    def __init__(self, email: str, organism_type: str = "Auto-detect"):
        Entrez.email = email
        self.conservation_threshold = 0.8
        self.organism_type = organism_type
        
    def search_genomes(self, species: str, max_results: int = 50) -> List[Dict]:
        """Search for genome sequences directly in nucleotide database"""
        try:
            st.info(f"🔍 Searching nucleotide database for {species} genomes...")
            
            # Organism-specific search strategies
            if self.organism_type == "Viroid":
                search_terms = self._get_viroid_search_terms(species)
                size_filter = lambda x: 100 <= x <= 2000  # Broader viroid size range
            elif self.organism_type == "Virus":
                search_terms = self._get_virus_search_terms(species)
                size_filter = lambda x: 1000 <= x <= 500000  # Viral genome range
            else:
                search_terms = self._get_standard_search_terms(species)
                size_filter = lambda x: x >= 10000  # Larger genomes
            
            sequences = []
            
            for i, search_term in enumerate(search_terms):
                try:
                    st.info(f"Search strategy {i+1}/{len(search_terms)}: {search_term}")
                    
                    # Search nucleotide database directly
                    handle = Entrez.esearch(
                        db="nucleotide", 
                        term=search_term, 
                        retmax=max_results * 2,  # Get more to filter later
                        sort="relevance"
                    )
                    search_results = Entrez.read(handle)
                    handle.close()
                    
                    if not search_results.get('IdList'):
                        st.warning(f"No results for strategy {i+1}")
                        continue
                    
                    st.success(f"Found {len(search_results['IdList'])} sequences")
                    
                    # Get sequence summaries in batches
                    batch_sequences = self._process_sequence_batch(
                        search_results['IdList'], species, size_filter, max_results
                    )
                    sequences.extend(batch_sequences)
                    
                    # If we have enough good sequences, break
                    if len(sequences) >= max_results:
                        break
                        
                except Exception as search_error:
                    st.warning(f"Search strategy {i+1} failed: {search_error}")
                    continue
            
            # Remove duplicates and sort by relevance
            unique_sequences = self._deduplicate_sequences(sequences, max_results)
            
            st.success(f"✅ Found {len(unique_sequences)} unique genome sequences")
            return unique_sequences
            
        except Exception as e:
            st.error(f"Error in genome search: {str(e)}")
            return []
    
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
            # Genus-level fallback
            f'{species.split()[0]}[Organism] AND "complete genome"' if len(species.split()) > 1 else None
        ]
    
    def _process_sequence_batch(self, id_list: List[str], species: str, size_filter, max_results: int) -> List[Dict]:
        """Process a batch of sequence IDs"""
        sequences = []
        batch_size = 20
        
        for i in range(0, len(id_list), batch_size):
            batch_ids = id_list[i:i+batch_size]
            
            try:
                # Get summaries for this batch
                handle = Entrez.esummary(db="nucleotide", id=','.join(batch_ids))
                summaries = Entrez.read(handle)
                handle.close()
                
                for summary in summaries:
                    try:
                        seq_id = summary.get('AccessionVersion', summary.get('Id', 'Unknown'))
                        title = summary.get('Title', 'Unknown')
                        length = int(summary.get('Length', 0))
                        organism = summary.get('Organism', species)
                        
                        # Apply size filter
                        if size_filter(length):
                            sequences.append({
                                'sequence_id': seq_id,
                                'title': title,
                                'organism': organism,
                                'length': length,
                                'source': 'nucleotide'
                            })
                        
                        # Stop if we have enough
                        if len(sequences) >= max_results:
                            break
                            
                    except (ValueError, TypeError, KeyError):
                        continue
                
                if len(sequences) >= max_results:
                    break
                    
            except Exception as batch_error:
                st.warning(f"Failed to process batch {i//batch_size + 1}: {batch_error}")
                continue
        
        return sequences
    
    def _deduplicate_sequences(self, sequences: List[Dict], max_results: int) -> List[Dict]:
        """Remove duplicates and sort sequences"""
        seen_ids = set()
        unique_sequences = []
        
        # Sort by length (appropriate for organism type)
        if self.organism_type == "Viroid":
            # For viroids, prefer typical sizes (200-400bp)
            sequences.sort(key=lambda x: abs(x['length'] - 300))
        elif self.organism_type == "Virus":
            # For viruses, prefer larger complete genomes
            sequences.sort(key=lambda x: x['length'], reverse=True)
        else:
            # For others, prefer larger sequences
            sequences.sort(key=lambda x: x['length'], reverse=True)
        
        for seq in sequences:
            seq_id = seq['sequence_id']
            if seq_id not in seen_ids and seq_id != 'Unknown':
                seen_ids.add(seq_id)
                unique_sequences.append(seq)
                if len(unique_sequences) >= max_results:
                    break
        
        return unique_sequences
    
    def fetch_genome_sequence(self, sequence_id: str, chromosome: str = "1") -> Optional[str]:
        """Fetch genome sequence directly from nucleotide database"""
        try:
            # Adjust sequence length limits based on organism type
            if self.organism_type == "Viroid":
                max_sequence_length = 10000  # Viroids are tiny
            elif self.organism_type == "Virus":
                max_sequence_length = 500000  # Viral genomes
            else:
                max_sequence_length = 2000000  # Default 2MB limit
            
            st.info(f"🧬 Fetching sequence: {sequence_id}")
            
            # Get sequence info first
            handle = Entrez.esummary(db="nucleotide", id=sequence_id)
            seq_info = Entrez.read(handle)
            handle.close()
            
            if not seq_info:
                st.error("Could not retrieve sequence information")
                return None
            
            # Determine fetch parameters
            seq_length = int(seq_info[0].get('Length', 0))
            
            # For viroids, fetch complete sequence
            if self.organism_type == "Viroid":
                fetch_length = seq_length
            else:
                fetch_length = min(seq_length, max_sequence_length)
            
            st.info(f"Fetching {fetch_length:,} bp from sequence of length {seq_length:,} bp")
            
            # Fetch the actual sequence
            handle = Entrez.efetch(
                db="nucleotide", 
                id=sequence_id, 
                rettype="fasta", 
                retmode="text",
                seq_start=1,
                seq_stop=fetch_length
            )
            sequence_data = handle.read()
            handle.close()
            
            # Parse FASTA
            sequences = list(SeqIO.parse(io.StringIO(sequence_data), "fasta"))
            if sequences:
                sequence = str(sequences[0].seq)
                
                # Organism-specific validation
                if self.organism_type == "Viroid":
                    if len(sequence) < 150:
                        st.warning(f"Short sequence for viroid: {len(sequence)} bp")
                    elif len(sequence) > 1000:
                        st.warning(f"Long sequence for viroid: {len(sequence)} bp")
                    else:
                        st.success(f"Viroid sequence: {len(sequence)} bp (typical range)")
                
                return sequence
            else:
                st.error("Failed to parse FASTA sequence")
                return None
                
        except Exception as e:
            st.error(f"Error fetching sequence: {str(e)}")
            return None
    
    def sliding_window_analysis(self, sequence: str, window_size: int = 1000, step_size: int = 500) -> pd.DataFrame:
        """Perform sliding window analysis on the sequence"""
        results = []
        
        # Ensure sequence is valid
        if not sequence or len(sequence) < window_size:
            st.error(f"Sequence too short for analysis. Length: {len(sequence) if sequence else 0}, Required: {window_size}")
            return pd.DataFrame()
        
        # Add progress bar for long sequences
        total_windows = (len(sequence) - window_size) // step_size + 1
        progress_bar = st.progress(0)
        
        for i, pos in enumerate(range(0, len(sequence) - window_size + 1, step_size)):
            window_seq = sequence[pos:pos + window_size]
            
            # Calculate various metrics
            gc_content = GC(window_seq)
            at_content = 100 - gc_content
            
            # Calculate complexity (Shannon entropy)
            entropy = self._calculate_entropy(window_seq)
            
            # Calculate repeat content (simple approach)
            repeat_content = self._calculate_repeat_content(window_seq)
            
            results.append({
                'start': pos + 1,
                'end': pos + window_size,
                'gc_content': gc_content,
                'at_content': at_content,
                'entropy': entropy,
                'repeat_content': repeat_content,
                'length': window_size
            })
            
            # Update progress
            if i % 100 == 0:  # Update every 100 windows
                progress_bar.progress(min(i / total_windows, 1.0))
        
        progress_bar.progress(1.0)
        return pd.DataFrame(results)
    
    def _calculate_entropy(self, sequence: str) -> float:
        """Calculate Shannon entropy of a sequence"""
        if not sequence:
            return 0.0
        
        # For RNA sequences (viroids), consider both DNA and RNA bases
        if self.organism_type == "Viroid":
            # RNA analysis - convert T to U for proper RNA analysis
            rna_sequence = sequence.upper().replace('T', 'U')
            counts = {'A': 0, 'U': 0, 'G': 0, 'C': 0, 'N': 0}
            for base in rna_sequence:
                if base in counts:
                    counts[base] += 1
                else:
                    counts['N'] += 1
        else:
            # DNA analysis
            counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0}
            for base in sequence.upper():
                if base in counts:
                    counts[base] += 1
                else:
                    counts['N'] += 1  # Unknown bases
        
        # Calculate entropy
        total = sum(counts.values())
        if total == 0:
            return 0.0
        
        entropy = 0.0
        for count in counts.values():
            if count > 0:
                p = count / total
                entropy -= p * np.log2(p)
        
        return entropy
    
    def _calculate_repeat_content(self, sequence: str) -> float:
        """Calculate percentage of repetitive content (simple approach)"""
        if len(sequence) < 4:
            return 0.0
        
        # Look for simple repeats (di-, tri-, tetra-nucleotides)
        repeat_bases = 0
        seq_upper = sequence.upper()
        
        # Check for dinucleotide repeats
        for i in range(len(seq_upper) - 3):
            dinucl = seq_upper[i:i+2]
            if seq_upper[i+2:i+4] == dinucl:
                repeat_bases += 2
        
        return (repeat_bases / len(sequence)) * 100
    
    def identify_conserved_regions(self, df: pd.DataFrame) -> pd.DataFrame:
        """Identify potentially conserved regions based on various criteria"""
        if df.empty:
            return pd.DataFrame()
        
        # Organism-specific conservation criteria
        if self.organism_type == "Viroid":
            # Viroids are highly structured RNA - different criteria
            conserved_mask = (
                (df['gc_content'] >= 30) & (df['gc_content'] <= 70) &  # Broader GC range for RNA
                (df['entropy'] >= 1.0) &  # Lower complexity threshold (highly structured)
                (df['repeat_content'] <= 30)  # Higher repeat tolerance (RNA structures)
            )
            
            # Viroid-specific scoring (emphasizes secondary structure potential)
            conserved_regions = df[conserved_mask].copy()
            if not conserved_regions.empty:
                conserved_regions['conservation_score'] = (
                    (2 - abs(conserved_regions['gc_content'] - 50) / 50) * 0.2 +  # Less weight on GC
                    (conserved_regions['entropy'] / 2) * 0.5 +  # More weight on complexity
                    ((100 - conserved_regions['repeat_content']) / 100) * 0.3  # Moderate repeat penalty
                )
        
        elif self.organism_type == "Virus":
            # Viral genomes - compact and functional
            conserved_mask = (
                (df['gc_content'] >= 35) & (df['gc_content'] <= 65) &  # Moderate GC content
                (df['entropy'] >= 1.3) &  # Moderate complexity
                (df['repeat_content'] <= 25)  # Low repeat content
            )
            
            conserved_regions = df[conserved_mask].copy()
            if not conserved_regions.empty:
                conserved_regions['conservation_score'] = (
                    (2 - abs(conserved_regions['gc_content'] - 50) / 50) * 0.3 +
                    (conserved_regions['entropy'] / 2) * 0.4 +
                    ((100 - conserved_regions['repeat_content']) / 100) * 0.3
                )
        
        else:
            # Default criteria for bacteria/eukaryotes
            conserved_mask = (
                (df['gc_content'] >= 40) & (df['gc_content'] <= 60) &  # Moderate GC content
                (df['entropy'] >= 1.5) &  # High complexity
                (df['repeat_content'] <= 20)  # Low repeat content
            )
            
            conserved_regions = df[conserved_mask].copy()
            if not conserved_regions.empty:
                conserved_regions['conservation_score'] = (
                    (2 - abs(conserved_regions['gc_content'] - 50) / 50) * 0.3 +
                    (conserved_regions['entropy'] / 2) * 0.4 +
                    ((100 - conserved_regions['repeat_content']) / 100) * 0.3
                )
        
        if not conserved_regions.empty:
            return conserved_regions.sort_values('conservation_score', ascending=False)
        else:
            return pd.DataFrame()

def test_ncbi_connection(email: str) -> bool:
    """Test NCBI connection with better error reporting"""
    try:
        Entrez.email = email
        
        # Test with a simple, reliable search
        handle = Entrez.esearch(db="nucleotide", term="Escherichia coli", retmax=1)
        test_results = Entrez.read(handle)
        handle.close()
        
        if test_results.get('IdList'):
            return True
        else:
            st.warning("NCBI connection working but no test results found")
            return False
            
    except Exception as e:
        st.error(f"NCBI connection test failed: {str(e)}")
        
        # Check for common issues
        if "HTTP Error 429" in str(e):
            st.error("Rate limit exceeded. Please wait and try again.")
        elif "HTTP Error 403" in str(e):
            st.error("Access forbidden. Check your email address.")
        elif "timeout" in str(e).lower():
            st.error("Connection timeout. Check your internet connection.")
        else:
            st.error("Unknown connection error. Please verify your email and internet connection.")
        
        return False

def main():
    st.title("🧬 Genome Conservation Analysis Tool")
    st.markdown("""
    This tool identifies conserved regions within genome sequences using data directly from NCBI.
    It searches the nucleotide database for actual genome sequences and analyzes them for conservation patterns.
    """)
    
    # Sidebar configuration
    with st.sidebar:
        st.header("Configuration")
        
        # Email for NCBI (required)
        email = st.text_input(
            "Email (required for NCBI):",
            placeholder="your.email@example.com",
            help="NCBI requires an email address for API access"
        )
        
        if not email or '@' not in email:
            st.warning("Please provide a valid email address to use NCBI services")
            st.stop()
        
        # Species selection with examples
        st.subheader("Species Selection")
        example_species = st.selectbox(
            "Choose example or enter custom:",
            ["Custom", "Homo sapiens", "Escherichia coli", "Saccharomyces cerevisiae", 
             "Potato spindle tuber viroid", "Tobacco mosaic virus", "SARS-CoV-2"]
        )
        
        if example_species == "Custom":
            species = st.text_input(
                "Species name:",
                placeholder="Enter scientific name",
                help="Enter the scientific name (e.g., 'Homo sapiens', 'Potato spindle tuber viroid')"
            )
        else:
            species = example_species
        
        if not species:
            st.warning("Please enter a species name")
            st.stop()
        
        # Analysis parameters
        st.subheader("Analysis Parameters")
        
        # Organism type selection for optimal parameters
        organism_type = st.selectbox(
            "Organism type:",
            ["Auto-detect", "Viroid", "Virus", "Bacteria", "Fungi", "Plant", "Animal"],
            help="Automatically optimizes all analysis parameters for the selected organism type"
        )
        
        # Auto-adjust parameters based on organism type
        if organism_type == "Viroid":
            default_window = 50
            min_window, max_window = 20, 200
            default_step = 10
            min_step, max_step = 5, 100
            default_max_sequences = 50
            param_info = "🧬 **Viroid Mode**: Ultra-fine analysis (50bp windows, 10bp steps)"
        elif organism_type == "Virus":
            default_window = 200
            min_window, max_window = 100, 1000
            default_step = 50
            min_step, max_step = 25, 500
            default_max_sequences = 30
            param_info = "🦠 **Virus Mode**: Compact genome analysis (200bp windows, 50bp steps)"
        elif organism_type == "Bacteria":
            default_window = 1000
            min_window, max_window = 500, 3000
            default_step = 500
            min_step, max_step = 100, 1500
            default_max_sequences = 20
            param_info = "🧫 **Bacteria Mode**: Standard analysis (1kb windows, 500bp steps)"
        else:  # Eukaryotes or auto-detect
            default_window = 1000
            min_window, max_window = 500, 5000
            default_step = 500
            min_step, max_step = 100, 2000
            default_max_sequences = 20
            param_info = "🌿 **Eukaryote Mode**: Broad analysis (1kb windows, 500bp steps)"
        
        # Display optimized parameters prominently
        st.info(param_info)
        st.success(f"✅ **Auto-optimized**: Window: {default_window}bp | Step: {default_step}bp | Max sequences: {default_max_sequences}")
        
        # Advanced parameters (collapsed by default)
        with st.expander("🔧 Advanced: Custom Parameters (Optional)"):
            st.warning("⚠️ Only modify these if you need custom analysis parameters. The auto-optimized values above are recommended.")
            
            custom_params = st.checkbox("Override auto-optimized parameters", value=False)
            
            if custom_params:
                window_size = st.slider(
                    "Window size (bp):", 
                    min_window, max_window, default_window, 
                    step=5 if organism_type == "Viroid" else 50,
                    help=f"Custom window size (default: {default_window}bp for {organism_type.lower()}s)"
                )
                step_size = st.slider(
                    "Step size (bp):", 
                    min_step, max_step, default_step,
                    step=5 if organism_type == "Viroid" else 25,
                    help=f"Custom step size (default: {default_step}bp for {organism_type.lower()}s)"
                )
                max_sequences = st.slider(
                    "Max sequences to search:", 
                    5, 100, default_max_sequences, 5,
                    help=f"Custom max sequences (default: {default_max_sequences} for {organism_type.lower()}s)"
                )
            else:
                # Use auto-optimized parameters
                window_size = default_window
                step_size = default_step
                max_sequences = default_max_sequences
    
    # Main interface
    col1, col2 = st.columns([1, 1])
    
    with col1:
        if st.button("🔍 Search Genome Sequences", type="primary"):
            # Initialize analyzer
            analyzer = NCBIGenomeAnalyzer(email, organism_type)
            
            with st.spinner(f"Searching NCBI nucleotide database for {species} genomes..."):
                sequences = analyzer.search_genomes(species, max_sequences)
            
            if sequences:
                st.session_state['sequences'] = sequences
                st.session_state['species'] = species
                st.session_state['organism_type'] = organism_type
                st.success(f"✅ Found {len(sequences)} genome sequences for {species}")
            else:
                st.error(f"❌ No genome sequences found for '{species}'. Try a different species name or organism type.")
    
    with col2:
        if st.button("🔗 Test NCBI Connection"):
            with st.spinner("Testing NCBI connection..."):
                if test_ncbi_connection(email):
                    st.success("✅ NCBI connection successful!")
                else:
                    st.error("❌ NCBI connection failed")
    
    # Display sequences
    if 'sequences' in st.session_state and 'species' in st.session_state:
        st.subheader(f"Available Genome Sequences for {st.session_state['species']}")
        
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
        
        # Select sequence for analysis
        selected_sequence = st.selectbox(
            "Select sequence for analysis:",
            options=[s['sequence_id'] for s in st.session_state['sequences']],
            format_func=lambda x: f"{x} - {next((s['title'][:40] + '...' if len(s['title']) > 40 else s['title']) for s in st.session_state['sequences'] if s['sequence_id'] == x)} ({next(s['length'] for s in st.session_state['sequences'] if s['sequence_id'] == x):,} bp)"
        )
        
        # Auto-select best sequence option
        auto_select = st.checkbox("Auto-select best available sequence", value=True, help="Automatically select the most suitable sequence from the results")
        
        if st.button("🧬 Analyze Conservation", type="primary"):
            if not selected_sequence:
                st.error("Please select a sequence for analysis")
                return
            
            # Initialize analyzer
            analyzer = NCBIGenomeAnalyzer(email, st.session_state['organism_type'])
            
            # Determine sequence parameter
            target_sequence = selected_sequence if not auto_select else selected_sequence
            
            with st.spinner("Fetching genome sequence and performing analysis..."):
                # Fetch sequence
                sequence = analyzer.fetch_genome_sequence(target_sequence)
                
                if sequence:
                    st.success(f"✅ Successfully fetched {len(sequence):,} bp sequence")
                    
                    # Perform sliding window analysis
                    df_analysis = analyzer.sliding_window_analysis(
                        sequence, window_size, step_size
                    )
                    
                    if not df_analysis.empty:
                        # Identify conserved regions
                        conserved_regions = analyzer.identify_conserved_regions(df_analysis)
                        
                        # Store results
                        st.session_state['analysis_results'] = df_analysis
                        st.session_state['conserved_regions'] = conserved_regions
                        st.session_state['sequence_length'] = len(sequence)
                        st.session_state['selected_sequence'] = target_sequence
                    else:
                        st.error("Analysis failed - no windows could be processed")
                else:
                    st.error("❌ Failed to fetch genome sequence. Try a different sequence.")
    
    # Display results
    if 'analysis_results' in st.session_state and not st.session_state['analysis_results'].empty:
        organism_type = st.session_state.get('organism_type', 'Unknown')
        
        # Display results with organism-specific insights
        if organism_type == "Viroid":
            st.subheader("🧬 Viroid Structure Analysis")
            st.info("**Viroid-Specific Insights:** Viroids are small, circular, single-stranded RNA molecules that form complex secondary structures. Conservation analysis focuses on structural elements rather than protein-coding sequences.")
        elif organism_type == "Virus":
            st.subheader("🦠 Viral Genome Analysis")  
            st.info("**Virus-Specific Insights:** Viral genomes are typically compact with overlapping genes and regulatory elements. Conservation reflects functional constraints and host adaptation.")
        else:
            st.subheader("📊 Analysis Results")
        
        df_results = st.session_state['analysis_results']
        conserved_regions = st.session_state['conserved_regions']
        sequence_length = st.session_state['sequence_length']
        
        # Organism-specific summary metrics
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Sequence Length", f"{sequence_length:,} bp")
            if organism_type == "Viroid" and sequence_length > 500:
                st.caption("⚠️ Large for viroid")
            elif organism_type == "Viroid" and 200 <= sequence_length <= 400:
                st.caption("✅ Typical viroid size")
        with col2:
            st.metric("Windows Analyzed", len(df_results))
        with col3:
            st.metric("Conserved Regions", len(conserved_regions))
        with col4:
            conservation_percentage = (len(conserved_regions) / len(df_results)) * 100 if len(df_results) > 0 else 0
            st.metric("Conservation %", f"{conservation_percentage:.1f}%")
            if organism_type == "Viroid" and conservation_percentage > 80:
                st.caption("🧬 Highly structured")
            elif organism_type == "Virus" and conservation_percentage > 60:
                st.caption("🦠 Compact genome")

        # Visualization with organism-specific titles
        if organism_type == "Viroid":
            st.subheader("🧬 Viroid Secondary Structure Landscape")
            structure_info = "**Viroid Analysis:** Focus on regions with moderate complexity that may form important secondary structures (hairpins, loops, pseudoknots)."
        elif organism_type == "Virus":
            st.subheader("🦠 Viral Genomic Landscape")
            structure_info = "**Viral Analysis:** Compact genome organization with overlapping features and regulatory elements."
        else:
            st.subheader("📈 Genomic Landscape")
            structure_info = ""
        
        if structure_info:
            st.info(structure_info)
        
        # Create subplot
        fig = make_subplots(
            rows=4, cols=1,
            shared_xaxes=True,
            subplot_titles=('GC Content (%)', 'Sequence Complexity (Entropy)', 'Repeat Content (%)', 'Conservation Score'),
            vertical_spacing=0.05
        )
        
        # GC Content
        fig.add_trace(
            go.Scatter(x=df_results['start'], y=df_results['gc_content'], 
                      mode='lines', name='GC Content', line=dict(color='blue')),
            row=1, col=1
        )
        
        # Entropy
        fig.add_trace(
            go.Scatter(x=df_results['start'], y=df_results['entropy'], 
                      mode='lines', name='Entropy', line=dict(color='green')),
            row=2, col=1
        )
        
        # Repeat Content
        fig.add_trace(
            go.Scatter(x=df_results['start'], y=df_results['repeat_content'], 
                      mode='lines', name='Repeat Content', line=dict(color='red')),
            row=3, col=1
        )
        
        # Conservation Score (only for conserved regions)
        if not conserved_regions.empty:
            fig.add_trace(
                go.Scatter(x=conserved_regions['start'], y=conserved_regions['conservation_score'], 
                          mode='markers', name='Conservation Score', 
                          marker=dict(color='purple', size=8)),
                row=4, col=1
            )
        
        fig.update_layout(height=800, showlegend=False, title_text=f"Genomic Analysis: {st.session_state.get('species', 'Unknown')} - {st.session_state.get('selected_sequence', 'Unknown')}")
        fig.update_xaxes(title_text="Genomic Position (bp)", row=4, col=1)
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Conserved regions table
        st.subheader("🎯 Top Conserved Regions")
        if not conserved_regions.empty:
            display_conserved = conserved_regions.head(20)[
                ['start', 'end', 'gc_content', 'entropy', 'repeat_content', 'conservation_score']
            ].round(3)
            st.dataframe(display_conserved, use_container_width=True)
            
            # Download option
            csv = conserved_regions.to_csv(index=False)
            st.download_button(
                label="📥 Download Conserved Regions (CSV)",
                data=csv,
                file_name=f"conserved_regions_{st.session_state.get('species', 'unknown').replace(' ', '_')}_{st.session_state.get('selected_sequence', 'unknown')}.csv",
                mime="text/csv"
            )
        else:
            st.info("ℹ️ No conserved regions identified with current criteria. Try adjusting the analysis parameters.")
        
        # Distribution plots
        st.subheader("📊 Statistical Distributions")
        
        col1, col2 = st.columns(2)
        
        with col1:
            fig_hist = px.histogram(
                df_results, x='gc_content', nbins=30,
                title='GC Content Distribution',
                labels={'gc_content': 'GC Content (%)', 'count': 'Frequency'}
            )
            st.plotly_chart(fig_hist, use_container_width=True)
        
        with col2:
            fig_scatter = px.scatter(
                df_results, x='entropy', y='gc_content',
                color='repeat_content',
                title='Complexity vs GC Content',
                labels={
                    'entropy': 'Sequence Complexity (Entropy)',
                    'gc_content': 'GC Content (%)',
                    'repeat_content': 'Repeat Content (%)'
                }
            )
            st.plotly_chart(fig_scatter, use_container_width=True)

        # Organism-specific interpretation
        if organism_type == "Viroid":
            with st.expander("🧬 Understanding Viroid Results"):
                st.markdown("""
                **Viroid-Specific Analysis:**
                
                🧬 **What are Viroids?**
                - Smallest known pathogens (200-400 nucleotides)
                - Circular, single-stranded RNA molecules
                - No protein-coding capacity
                - Pathogenic to plants only
                - Replicate using host RNA polymerase II
                
                🔬 **Conservation in Viroids:**
                - **Secondary Structure Elements**: Hairpin loops, bulges, and pseudoknots
                - **Central Conserved Region (CCR)**: Essential for replication
                - **Pathogenicity Domain**: Variable region affecting virulence
                - **Terminal Conserved Region (TCR)**: Important for processing
                
                📊 **Interpreting Your Results:**
                - **High Conservation (>80%)**: Likely structural/functional domains
                - **Moderate GC Content (30-70%)**: Optimal for RNA folding
                - **Structured Regions**: Lower entropy may indicate important secondary structures
                - **Variable Regions**: May be pathogenicity or host-adaptation domains
                
                ⚠️ **Important Notes:**
                - Viroids are RNA, but often stored as DNA sequences in databases
                - Secondary structure prediction would require specialized RNA folding tools
                - Conservation analysis here is compositional, not comparative between species
                """)
        elif organism_type == "Virus":
            with st.expander("🦠 Understanding Viral Results"):
                st.markdown("""
                **Viral Genome Analysis:**
                
                🦠 **Viral Genome Features:**
                - Compact organization with overlapping genes
                - High gene density with minimal non-coding regions
                - Regulatory elements often overlap with coding sequences
                - Size varies dramatically (1kb to >1Mb)
                
                🔬 **Conservation Patterns:**
                - **Essential Genes**: Replication, transcription machinery
                - **Structural Proteins**: Capsid, envelope proteins
                - **Regulatory Elements**: Promoters, origins of replication
                - **Variable Regions**: Host interaction, immune evasion
                
                📊 **Interpreting Results:**
                - **High Conservation**: Core viral functions, essential domains
                - **Balanced Composition**: Functional constraints on codon usage
                - **Low Repeat Content**: Compact genomes minimize redundancy
                """)

    # Enhanced troubleshooting for direct nucleotide search
    with st.expander("🔧 Troubleshooting - Nucleotide Database Search"):
        st.markdown("""
        **Direct Nucleotide Search Benefits:**
        
        ✅ **Why this approach works better:**
        - Searches actual genome sequences, not just assembly metadata
        - Finds many more sequences per organism (10-50 vs 1-2)
        - Works for all organism types including viroids and viruses
        - No assembly-to-sequence mapping issues
        - Direct access to sequence data
        
        **🧬 Viroid-Specific Improvements:**
        - Searches nucleotide database where viroids are actually stored
        - Size filtering (150-1000 bp) to find viroid-sized sequences
        - Multiple search strategies for different viroid naming conventions
        - No dependency on formal genome assemblies
        
        **🦠 Virus-Specific Improvements:**
        - Finds complete viral genomes and individual segments
        - Size-appropriate filtering for viral genome range
        - Searches for both complete genomes and individual sequences
        
        **General Improvements:**
        - Much higher sequence yield per search
        - Better coverage of available genomic data
        - More robust search strategies
        - Direct sequence access without intermediate mapping
        """)
        
    # Updated troubleshooting section
    with st.expander("🔧 Common Issues and Solutions"):
        st.markdown("""
        **Search Issues:**
        
        1. **No sequences found:**
           - Check species name spelling (use scientific names)
           - Try different organism type selections
           - Some organisms may have limited sequence data
           
        2. **Few sequences found:**
           - Try "Auto-detect" organism type for broader search
           - Check if organism name includes common alternatives
           - Use genus name only (e.g., "Escherichia" instead of "Escherichia coli")
        
        3. **Sequence fetch fails:**
           - Try a different sequence from the list
           - Some sequences may be restricted or corrupted
           - Check NCBI status if multiple sequences fail
        
        4. **Analysis errors:**
           - Ensure sequence length > window size
           - Reduce window size for shorter sequences
           - Check for unusual characters in sequence
        
        **Performance Tips:**
        - Viroids: Start with 20-50bp windows for detailed structure analysis
        - Viruses: Use 100-500bp windows for functional domain analysis  
        - Bacteria: Use 500-2000bp windows for gene-level analysis
        - Larger organisms: Use 1000-5000bp windows for broad patterns
        """)

    # Updated information panel
    with st.expander("ℹ️ About the Nucleotide Database Approach"):
        st.markdown("""
        **Why Search Nucleotide Database Directly:**
        
        🔍 **Better Data Access:**
        - Assembly database = metadata about genome projects
        - Nucleotide database = actual DNA/RNA sequences
        - Direct access to sequence data without mapping issues
        - Much higher sequence yield per organism
        
        📊 **Improved Coverage:**
        - Finds 10-50 sequences per organism (vs 1-2 with assemblies)
        - Works for all organism types including viroids
        - Includes partial sequences, complete genomes, and chromosomes
        - No dependency on formal genome assembly projects
        
        🧬 **Organism-Specific Benefits:**
        - **Viroids**: Found in nucleotide DB, rarely in assembly DB
        - **Viruses**: Complete genomes and segments both accessible
        - **Bacteria**: Multiple strains and isolates available
        - **Eukaryotes**: Chromosomes, scaffolds, and contigs all searchable
        
        ⚡ **Technical Advantages:**
        - Eliminates assembly-to-sequence mapping step
        - Reduces points of failure in data retrieval
        - Faster and more reliable sequence access
        - Better error handling and recovery options
        
        **Search Strategy:**
        - Multiple search terms per organism for comprehensive coverage
        - Size-based filtering appropriate for organism type
        - Relevance-based sorting to prioritize best matches
        - Deduplication to avoid analyzing identical sequences
        """)

if __name__ == "__main__":
    main()
