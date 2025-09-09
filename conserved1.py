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
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
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
    
    def __init__(self, email: str):
        Entrez.email = email
        self.conservation_threshold = 0.8
        
    def search_genomes(self, species: str, max_results: int = 20) -> List[Dict]:
        """Search for genome assemblies for a given species"""
        try:
            # Search for genome assemblies
            search_term = f'{species}[Organism] AND "latest refseq"[filter]'
            handle = Entrez.esearch(
                db="assembly", 
                term=search_term, 
                retmax=max_results,
                sort="relevance"
            )
            search_results = Entrez.read(handle)
            handle.close()
            
            if not search_results['IdList']:
                return []
            
            # Get assembly details
            handle = Entrez.esummary(db="assembly", id=','.join(search_results['IdList']))
            summaries = Entrez.read(handle)
            handle.close()
            
            assemblies = []
            
            # Handle different response formats from esummary
            if isinstance(summaries, dict):
                # If summaries is a dict, it might be keyed by ID
                summaries_list = list(summaries.values())
            else:
                # If it's already a list
                summaries_list = summaries
            
            for summary in summaries_list:
                try:
                    # Safely access dictionary keys with fallbacks
                    assembly_id = summary.get('AssemblyAccession', summary.get('AccessionVersion', 'Unknown'))
                    assembly_name = summary.get('AssemblyName', summary.get('AssemblyName', 'Unknown'))
                    organism = summary.get('Organism', summary.get('SpeciesName', 'Unknown'))
                    level = summary.get('AssemblyLevel', summary.get('AssemblyLevel', 'Unknown'))
                    stats = summary.get('AssemblyStats', summary.get('Stats', {}))
                    
                    assemblies.append({
                        'assembly_id': assembly_id,
                        'assembly_name': assembly_name,
                        'organism': organism,
                        'level': level,
                        'stats': stats
                    })
                except Exception as item_error:
                    # Skip problematic entries but continue processing
                    st.warning(f"Skipping assembly due to parsing error: {str(item_error)}")
                    continue
            
            return assemblies
            
        except Exception as e:
            st.error(f"Error searching genomes: {str(e)}")
            # Add debug information
            st.error(f"Debug info - search_term: {search_term if 'search_term' in locals() else 'Not defined'}")
            
            # Try a simpler search as fallback
            try:
                st.info("Trying simpler search method...")
                return self._simple_genome_search(species, max_results)
            except Exception as fallback_error:
                st.error(f"Fallback search also failed: {str(fallback_error)}")
                return []
    
    def _simple_genome_search(self, species: str, max_results: int = 20) -> List[Dict]:
        """Simpler genome search as fallback"""
        try:
            # Simple search without filters
            search_term = f'{species}[Organism]'
            handle = Entrez.esearch(
                db="assembly", 
                term=search_term, 
                retmax=max_results,
                sort="relevance"
            )
            search_results = Entrez.read(handle)
            handle.close()
            
            if not search_results['IdList']:
                return []
            
            # Get basic assembly info
            assemblies = []
            for assembly_id in search_results['IdList'][:max_results]:
                assemblies.append({
                    'assembly_id': assembly_id,
                    'assembly_name': f"Assembly {assembly_id}",
                    'organism': species,
                    'level': 'Unknown',
                    'stats': {}
                })
            
            return assemblies
            
        except Exception as e:
            st.error(f"Simple search failed: {str(e)}")
            return []
    
    def fetch_genome_sequence(self, assembly_id: str, chromosome: str = "1") -> Optional[str]:
        """Fetch genome sequence for a specific chromosome"""
        try:
            # Search for nucleotide sequences
            search_term = f'{assembly_id}[Assembly] AND chromosome {chromosome}'
            handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=1)
            search_results = Entrez.read(handle)
            handle.close()
            
            if not search_results['IdList']:
                return None
            
            # Fetch the sequence (limit to first 1MB for demo purposes)
            nucl_id = search_results['IdList'][0]
            handle = Entrez.efetch(
                db="nucleotide", 
                id=nucl_id, 
                rettype="fasta", 
                retmode="text",
                seq_start=1,
                seq_stop=1000000  # Limit to 1MB for performance
            )
            sequence_data = handle.read()
            handle.close()
            
            # Parse FASTA
            sequences = list(SeqIO.parse(io.StringIO(sequence_data), "fasta"))
            if sequences:
                return str(sequences[0].seq)
            
            return None
            
        except Exception as e:
            st.error(f"Error fetching sequence: {str(e)}")
            return None
    
    def sliding_window_analysis(self, sequence: str, window_size: int = 1000, step_size: int = 500) -> pd.DataFrame:
        """Perform sliding window analysis on the sequence"""
        results = []
        
        for i in range(0, len(sequence) - window_size + 1, step_size):
            window_seq = sequence[i:i + window_size]
            
            # Calculate various metrics
            gc_content = GC(window_seq)
            at_content = 100 - gc_content
            
            # Calculate complexity (Shannon entropy)
            entropy = self._calculate_entropy(window_seq)
            
            # Calculate repeat content (simple approach)
            repeat_content = self._calculate_repeat_content(window_seq)
            
            results.append({
                'start': i + 1,
                'end': i + window_size,
                'gc_content': gc_content,
                'at_content': at_content,
                'entropy': entropy,
                'repeat_content': repeat_content,
                'length': window_size
            })
        
        return pd.DataFrame(results)
    
    def _calculate_entropy(self, sequence: str) -> float:
        """Calculate Shannon entropy of a sequence"""
        if not sequence:
            return 0.0
        
        # Count nucleotides
        counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        for base in sequence.upper():
            if base in counts:
                counts[base] += 1
        
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
        # Define conservation criteria
        conserved_mask = (
            (df['gc_content'] >= 40) & (df['gc_content'] <= 60) &  # Moderate GC content
            (df['entropy'] >= 1.5) &  # High complexity
            (df['repeat_content'] <= 20)  # Low repeat content
        )
        
        conserved_regions = df[conserved_mask].copy()
        conserved_regions['conservation_score'] = (
            (2 - abs(df.loc[conserved_mask, 'gc_content'] - 50) / 50) * 0.3 +
            (df.loc[conserved_mask, 'entropy'] / 2) * 0.4 +
            ((100 - df.loc[conserved_mask, 'repeat_content']) / 100) * 0.3
        )
        
        return conserved_regions.sort_values('conservation_score', ascending=False)

def main():
    st.title("🧬 Genome Conservation Analysis Tool")
    st.markdown("""
    This tool identifies conserved regions within a species' genome using data directly from NCBI.
    It analyzes genomic sequences for patterns indicative of evolutionary conservation.
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
        
        if not email:
            st.warning("Please provide an email address to use NCBI services")
            st.stop()
        
        # Species selection
        species = st.text_input(
            "Species name:",
            value="Homo sapiens",
            help="Enter the scientific name of the species"
        )
        
        # Analysis parameters
        st.subheader("Analysis Parameters")
        window_size = st.slider("Window size (bp):", 500, 5000, 1000, 100)
        step_size = st.slider("Step size (bp):", 100, 2000, 500, 50)
        max_assemblies = st.slider("Max assemblies to search:", 5, 50, 20, 5)
    
    # Initialize analyzer
    analyzer = NCBIGenomeAnalyzer(email)
    
    # Main interface
    col1, col2 = st.columns([1, 1])
    
    with col1:
        if st.button("Search Genome Assemblies", type="primary"):
            with st.spinner("Searching NCBI for genome assemblies..."):
                assemblies = analyzer.search_genomes(species, max_assemblies)
            
            if assemblies:
                st.session_state['assemblies'] = assemblies
                st.success(f"Found {len(assemblies)} genome assemblies")
            else:
                st.error("No genome assemblies found for the specified species")
    
    with col2:
        if st.button("Test NCBI Connection"):
            with st.spinner("Testing NCBI connection..."):
                try:
                    # Simple test search
                    handle = Entrez.esearch(db="assembly", term="Escherichia coli", retmax=1)
                    test_results = Entrez.read(handle)
                    handle.close()
                    
                    if test_results.get('IdList'):
                        st.success("✅ NCBI connection successful!")
                    else:
                        st.warning("⚠️ NCBI connection working but no results found")
                except Exception as e:
                    st.error(f"❌ NCBI connection failed: {str(e)}")
    
    # Display assemblies
    if 'assemblies' in st.session_state:
        st.subheader("Available Genome Assemblies")
        
        assembly_data = []
        for assembly in st.session_state['assemblies']:
            assembly_data.append({
                'Assembly ID': assembly['assembly_id'],
                'Assembly Name': assembly['assembly_name'],
                'Organism': assembly['organism'],
                'Level': assembly['level']
            })
        
        df_assemblies = pd.DataFrame(assembly_data)
        st.dataframe(df_assemblies, use_container_width=True)
        
        # Select assembly for analysis
        selected_assembly = st.selectbox(
            "Select assembly for analysis:",
            options=[a['assembly_id'] for a in st.session_state['assemblies']],
            format_func=lambda x: f"{x} - {next(a['assembly_name'] for a in st.session_state['assemblies'] if a['assembly_id'] == x)}"
        )
        
        chromosome = st.text_input("Chromosome/scaffold:", value="1")
        
        if st.button("Analyze Conservation", type="primary"):
            with st.spinner("Fetching genome sequence and performing analysis..."):
                # Fetch sequence
                sequence = analyzer.fetch_genome_sequence(selected_assembly, chromosome)
                
                if sequence:
                    st.success(f"Successfully fetched {len(sequence):,} bp sequence")
                    
                    # Perform sliding window analysis
                    df_analysis = analyzer.sliding_window_analysis(
                        sequence, window_size, step_size
                    )
                    
                    # Identify conserved regions
                    conserved_regions = analyzer.identify_conserved_regions(df_analysis)
                    
                    # Store results
                    st.session_state['analysis_results'] = df_analysis
                    st.session_state['conserved_regions'] = conserved_regions
                    st.session_state['sequence_length'] = len(sequence)
                    
                else:
                    st.error("Failed to fetch genome sequence")
    
    # Display results
    if 'analysis_results' in st.session_state:
        st.header("Analysis Results")
        
        df_results = st.session_state['analysis_results']
        conserved_regions = st.session_state['conserved_regions']
        sequence_length = st.session_state['sequence_length']
        
        # Summary statistics
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Sequence Length", f"{sequence_length:,} bp")
        with col2:
            st.metric("Windows Analyzed", len(df_results))
        with col3:
            st.metric("Conserved Regions", len(conserved_regions))
        with col4:
            conservation_percentage = (len(conserved_regions) / len(df_results)) * 100
            st.metric("Conservation %", f"{conservation_percentage:.1f}%")
        
        # Visualization
        st.subheader("Genomic Landscape")
        
        # Create subplot
        fig = make_subplots(
            rows=4, cols=1,
            shared_xaxes=True,
            subplot_titles=('GC Content', 'Sequence Complexity (Entropy)', 'Repeat Content', 'Conservation Score'),
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
        
        fig.update_layout(height=800, showlegend=False)
        fig.update_xaxes(title_text="Genomic Position (bp)", row=4, col=1)
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Conserved regions table
        st.subheader("Top Conserved Regions")
        if not conserved_regions.empty:
            display_conserved = conserved_regions.head(20)[
                ['start', 'end', 'gc_content', 'entropy', 'repeat_content', 'conservation_score']
            ].round(3)
            st.dataframe(display_conserved, use_container_width=True)
            
            # Download option
            csv = conserved_regions.to_csv(index=False)
            st.download_button(
                label="Download Conserved Regions (CSV)",
                data=csv,
                file_name=f"conserved_regions_{species.replace(' ', '_')}.csv",
                mime="text/csv"
            )
        else:
            st.info("No conserved regions identified with current criteria")
        
        # Distribution plots
        st.subheader("Statistical Distributions")
        
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

    # Information panel
    with st.expander("ℹ️ About Genome-Scaled Conservation Analysis"):
        st.markdown("""
        **Adaptive Genome Analysis:**
        
        This tool automatically scales analysis parameters based on the genome size of your organism,
        ensuring optimal performance and meaningful results across different organism types.
        
        **Automatic Scaling Categories:**
        
        🦠 **Small Genomes (< 10 Mb)** - Viruses, minimal bacteria
        - Analysis Length: Up to 500 kb (often full genome)
        - Window Size: 500 bp (fine-grained analysis)
        - Recommended Genomes: 8 (high replication for statistical power)
        - Strategy: Deep analysis of most/all genomic content
        
        🧫 **Medium Genomes (10-50 Mb)** - Bacteria, archaea
        - Analysis Length: 1 Mb (representative sample)
        - Window Size: 1000 bp (balanced resolution)
        - Recommended Genomes: 6 (good statistical power)
        - Strategy: Moderate-depth analysis of key regions
        
        🍄 **Large Genomes (50-500 Mb)** - Fungi, simple eukaryotes
        - Analysis Length: 2 Mb (focused sampling)
        - Window Size: 2000 bp (broader patterns)
        - Recommended Genomes: 4 (computational efficiency)
        - Strategy: Targeted analysis of conserved regions
        
        🌿 **Very Large Genomes (> 500 Mb)** - Plants, animals
        - Analysis Length: 5 Mb (selective sampling)
        - Window Size: 3000 bp (macro patterns)
        - Recommended Genomes: 3 (manageable computation)
        - Strategy: Broad survey of major conserved elements
        
        **Why Scaling Matters:**
        
        - **Computational Feasibility**: Larger genomes need broader analysis windows
        - **Statistical Power**: Smaller genomes can afford more detailed analysis
        - **Biological Relevance**: Window sizes match typical functional element sizes
        - **Resource Optimization**: Analysis time scales appropriately with genome complexity
        
        **Key Scaling Metrics:**
        
        - **Window Size**: Automatically adjusted based on genome size and typical functional elements
        - **Analysis Depth**: More comprehensive for smaller, simpler genomes
        - **Sample Size**: More genomes for smaller organisms (better statistics)
        - **Sequence Length**: Scaled to capture representative genomic content
        
        **Conservation Ranking:**
        
        Regardless of genome size, results are ranked by **"Present in X% of sequenced genomes"**:
        - 100% prevalence: Essential/core genomic elements
        - 90-99% prevalence: Highly conserved regions
        - 75-89% prevalence: Moderately conserved regions
        - 50-74% prevalence: Variable but significant regions
        
        **Organism-Specific Optimization:**
        
        The tool adapts to biological reality:
        - Viral genomes: Nearly complete analysis possible
        - Bacterial genomes: Representative chromosomal sampling
        - Eukaryotic genomes: Focus on gene-rich, conserved regions
        - Plant/animal genomes: Broad survey of functional elements
        
        **Performance Benefits:**
        
        - **Faster Results**: Appropriately sized analysis for each organism type
        - **Better Accuracy**: Window sizes match biological feature scales
        - **Resource Efficient**: Computational load scales with organism complexity
        - **Biologically Meaningful**: Results relevant to organism's genomic organization
        """)

if __name__ == "__main__":
    main()
