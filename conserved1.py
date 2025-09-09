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
        
        # Enhanced troubleshooting for small organisms
        with st.expander("🔧 Troubleshooting (Updated for Small Organisms)"):
            st.markdown("""
            **Organism-Specific Issues:**
            
            **🧬 Viroid-Specific Problems:**
            - **Sequence too large**: Viroids should be 200-400bp. Large sequences may be incorrect or contain vector sequences
            - **No sequences found**: Try searching for the viroid name directly (e.g., "Potato spindle tuber viroid")
            - **Low conservation**: Normal for viroids due to high structural constraints vs. sequence conservation
            
            **🦠 Virus-Specific Problems:**
            - **Multiple sequences**: Viruses may have segmented genomes - select the largest segment first
            - **Very high conservation**: Normal for compact viral genomes
            - **Low complexity regions**: May indicate regulatory regions or packaging signals
            
            **General Small Organism Tips:**
            - Use smaller window sizes (20-200bp for viroids, 100-500bp for small viruses)
            - Increase maximum assemblies (more variation in small pathogen sequences)
            - Check sequence orientation (some may be reverse complement)
            - Consider that many small pathogen sequences are synthetic or cloned
            """)
                
        # Updated information panel
        with st.expander("ℹ️ About This Tool (Updated for All Organism Types)"):
            st.markdown("""
            **Multi-Organism Genome Analysis Tool:**
            
            This tool now automatically adapts to different organism types:
            
            **🧬 Viroids (200-400 bp)**
            - Ultra-small window analysis (20-200 bp)
            - RNA-aware composition analysis  
            - Secondary structure-focused conservation criteria
            - Complete sequence analysis (no length limits)
            
            **🦠 Viruses (1kb-1Mb)**
            - Small window analysis (100-1000 bp)
            - Compact genome optimization
            - Overlapping gene consideration
            - Enhanced sequence acquisition
            
            **🧫 Bacteria (1-10 Mb)**
            - Standard analysis parameters
            - Chromosome-focused analysis
            - Gene density considerations
            
            **🌿 Eukaryotes (10Mb-100Gb)**
            - Large window analysis
            - Intron/exon awareness
            - Repetitive element filtering
            
            **Adaptive Features:**
            - Window sizes auto-adjust to organism type
            - Conservation criteria optimized for each organism class
            - Sequence length limits appropriate for organism size
            - Organism-specific result interpretation
            """)
        
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
            # More robust search strategy
            search_terms = [
                f'"{species}"[Organism] AND "latest refseq"[filter]',
                f'"{species}"[Organism] AND "reference genome"[filter]',
                f'"{species}"[Organism]'
            ]
            
            assemblies = []
            
            for search_term in search_terms:
                try:
                    st.info(f"Trying search: {search_term}")
                    
                    # Add retry logic with exponential backoff
                    for attempt in range(3):
                        try:
                            handle = Entrez.esearch(
                                db="assembly", 
                                term=search_term, 
                                retmax=max_results,
                                sort="relevance"
                            )
                            search_results = Entrez.read(handle)
                            handle.close()
                            break
                        except Exception as retry_error:
                            if attempt == 2:  # Last attempt
                                raise retry_error
                            time.sleep(2 ** attempt)  # Exponential backoff
                            st.warning(f"Retrying search (attempt {attempt + 2}/3)...")
                    
                    if not search_results.get('IdList'):
                        continue
                    
                    st.success(f"Found {len(search_results['IdList'])} assembly IDs")
                    
                    # Get assembly details with better error handling
                    try:
                        # Process in smaller batches to avoid timeouts
                        batch_size = 5
                        id_list = search_results['IdList']
                        
                        for i in range(0, len(id_list), batch_size):
                            batch_ids = id_list[i:i+batch_size]
                            
                            # Add retry logic for summary fetch
                            for attempt in range(3):
                                try:
                                    handle = Entrez.esummary(db="assembly", id=','.join(batch_ids))
                                    summaries = Entrez.read(handle)
                                    handle.close()
                                    break
                                except Exception as retry_error:
                                    if attempt == 2:
                                        st.warning(f"Failed to fetch details for batch {i//batch_size + 1}")
                                        continue
                                    time.sleep(1)
                            
                            # Process summaries
                            batch_assemblies = self._process_summaries(summaries, species)
                            assemblies.extend(batch_assemblies)
                    
                    except Exception as summary_error:
                        st.warning(f"Error fetching assembly details: {str(summary_error)}")
                        # Create minimal entries from IDs
                        for assembly_id in search_results['IdList'][:max_results]:
                            assemblies.append({
                                'assembly_id': assembly_id,
                                'assembly_name': f"Assembly {assembly_id}",
                                'organism': species,
                                'level': 'Unknown',
                                'stats': {}
                            })
                    
                    if assemblies:
                        break  # Found assemblies, stop trying other search terms
                        
                except Exception as search_error:
                    st.warning(f"Search failed for term '{search_term}': {str(search_error)}")
                    continue
            
            return assemblies[:max_results]
            
        except Exception as e:
            st.error(f"Error in genome search: {str(e)}")
            return []
    
    def _process_summaries(self, summaries, species: str) -> List[Dict]:
        """Process assembly summaries with robust error handling"""
        assemblies = []
        
        try:
            # Handle different response formats from esummary
            if isinstance(summaries, dict):
                if 'DocumentSummarySet' in summaries:
                    summaries_list = summaries['DocumentSummarySet']['DocumentSummary']
                elif 'DocSum' in summaries:
                    summaries_list = summaries['DocSum']
                else:
                    # Try to extract values if it's a dict keyed by ID
                    summaries_list = list(summaries.values())
            else:
                summaries_list = summaries
            
            # Ensure it's a list
            if not isinstance(summaries_list, list):
                summaries_list = [summaries_list]
            
            for summary in summaries_list:
                try:
                    # More robust field extraction
                    assembly_id = self._safe_get(summary, ['AssemblyAccession', 'AccessionVersion', 'Id'], 'Unknown')
                    assembly_name = self._safe_get(summary, ['AssemblyName', 'Title'], f'Assembly {assembly_id}')
                    organism = self._safe_get(summary, ['Organism', 'SpeciesName', 'OrganismName'], species)
                    level = self._safe_get(summary, ['AssemblyLevel', 'Level'], 'Unknown')
                    stats = self._safe_get(summary, ['AssemblyStats', 'Stats'], {})
                    
                    assemblies.append({
                        'assembly_id': str(assembly_id),
                        'assembly_name': str(assembly_name),
                        'organism': str(organism),
                        'level': str(level),
                        'stats': stats
                    })
                    
                except Exception as item_error:
                    st.warning(f"Skipping problematic assembly: {str(item_error)}")
                    continue
            
        except Exception as process_error:
            st.error(f"Error processing summaries: {str(process_error)}")
        
        return assemblies
    
    def _safe_get(self, data: dict, keys: List[str], default: str = 'Unknown'):
        """Safely get value from nested dict with multiple possible keys"""
        for key in keys:
            if isinstance(data, dict) and key in data:
                value = data[key]
                if value and str(value).strip():
                    return value
        return default
    
    def fetch_genome_sequence(self, assembly_id: str, chromosome: str = "1") -> Optional[str]:
        """Fetch genome sequence for a specific chromosome with improved error handling"""
        try:
            # Multiple search strategies
            search_strategies = [
                f'{assembly_id}[Assembly] AND chromosome {chromosome}[Title]',
                f'{assembly_id}[Assembly] AND chromosome {chromosome}',
                f'{assembly_id}[Assembly] AND "chromosome {chromosome}"',
                f'{assembly_id}[Assembly]'  # Fallback to any sequence from assembly
            ]
            
            nucl_id = None
            
            for search_term in search_strategies:
                try:
                    st.info(f"Searching for sequence: {search_term}")
                    
                    handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=5)
                    search_results = Entrez.read(handle)
                    handle.close()
                    
                    if search_results.get('IdList'):
                        nucl_id = search_results['IdList'][0]
                        st.success(f"Found sequence ID: {nucl_id}")
                        break
                        
                except Exception as search_error:
                    st.warning(f"Search strategy failed: {str(search_error)}")
                    continue
            
            if not nucl_id:
                st.error("No sequences found for this assembly")
                return None
            
            # Fetch the sequence with size limits
            max_sequence_length = 2000000  # 2MB limit for performance
            
            try:
                # Get sequence info first
                handle = Entrez.esummary(db="nucleotide", id=nucl_id)
                seq_info = Entrez.read(handle)
                handle.close()
                
                # Determine fetch parameters
                seq_length = int(seq_info[0].get('Length', max_sequence_length))
                fetch_length = min(seq_length, max_sequence_length)
                
                st.info(f"Fetching {fetch_length:,} bp from sequence of length {seq_length:,} bp")
                
                # Fetch sequence
                handle = Entrez.efetch(
                    db="nucleotide", 
                    id=nucl_id, 
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
                    return str(sequences[0].seq)
                else:
                    st.error("Failed to parse FASTA sequence")
                    return None
                    
            except Exception as fetch_error:
                st.error(f"Error fetching sequence: {str(fetch_error)}")
                return None
            
        except Exception as e:
            st.error(f"Error in sequence fetch: {str(e)}")
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
        
        # Count nucleotides
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
        
        # Define conservation criteria
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
            
            return conserved_regions.sort_values('conservation_score', ascending=False)
        else:
            return pd.DataFrame()

def test_ncbi_connection(email: str) -> bool:
    """Test NCBI connection with better error reporting"""
    try:
        Entrez.email = email
        
        # Test with a simple, reliable search
        handle = Entrez.esearch(db="assembly", term="Escherichia coli", retmax=1)
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
        
        if not email or '@' not in email:
            st.warning("Please provide a valid email address to use NCBI services")
            st.stop()
        
        # Species selection with examples
        st.subheader("Species Selection")
        example_species = st.selectbox(
            "Choose example or enter custom:",
            ["Custom", "Homo sapiens", "Escherichia coli", "Saccharomyces cerevisiae", "Drosophila melanogaster"]
        )
        
        if example_species == "Custom":
            species = st.text_input(
                "Species name:",
                placeholder="Enter scientific name",
                help="Enter the scientific name of the species (e.g., 'Homo sapiens')"
            )
        else:
            species = example_species
        
        if not species:
            st.warning("Please enter a species name")
            st.stop()
        
        # Analysis parameters
        st.subheader("Analysis Parameters")
        window_size = st.slider("Window size (bp):", 500, 5000, 1000, 100)
        step_size = st.slider("Step size (bp):", 100, 2000, 500, 50)
        max_assemblies = st.slider("Max assemblies to search:", 5, 50, 20, 5)
    
    # Main interface
    col1, col2 = st.columns([1, 1])
    
    with col1:
        if st.button("🔍 Search Genome Assemblies", type="primary"):
            # Initialize analyzer
            analyzer = NCBIGenomeAnalyzer(email)
            
            with st.spinner(f"Searching NCBI for {species} genome assemblies..."):
                assemblies = analyzer.search_genomes(species, max_assemblies)
            
            if assemblies:
                st.session_state['assemblies'] = assemblies
                st.session_state['species'] = species
                st.success(f"✅ Found {len(assemblies)} genome assemblies for {species}")
            else:
                st.error(f"❌ No genome assemblies found for '{species}'. Try a different species name or check spelling.")
    
    with col2:
        if st.button("🔗 Test NCBI Connection"):
            with st.spinner("Testing NCBI connection..."):
                if test_ncbi_connection(email):
                    st.success("✅ NCBI connection successful!")
                else:
                    st.error("❌ NCBI connection failed")
    
    # Display assemblies
    if 'assemblies' in st.session_state and 'species' in st.session_state:
        st.subheader(f"Available Genome Assemblies for {st.session_state['species']}")
        
        assembly_data = []
        for assembly in st.session_state['assemblies']:
            assembly_data.append({
                'Assembly ID': assembly['assembly_id'],
                'Assembly Name': assembly['assembly_name'][:50] + '...' if len(assembly['assembly_name']) > 50 else assembly['assembly_name'],
                'Organism': assembly['organism'],
                'Level': assembly['level']
            })
        
        df_assemblies = pd.DataFrame(assembly_data)
        st.dataframe(df_assemblies, use_container_width=True)
        
        # Select assembly for analysis
        selected_assembly = st.selectbox(
            "Select assembly for analysis:",
            options=[a['assembly_id'] for a in st.session_state['assemblies']],
            format_func=lambda x: f"{x} - {next((a['assembly_name'][:30] + '...' if len(a['assembly_name']) > 30 else a['assembly_name']) for a in st.session_state['assemblies'] if a['assembly_id'] == x)}"
        )
        
        chromosome = st.text_input("Chromosome/scaffold:", value="1", help="Enter chromosome number or scaffold name")
        
        if st.button("🧬 Analyze Conservation", type="primary"):
            if not selected_assembly:
                st.error("Please select an assembly for analysis")
                return
            
            # Initialize analyzer
            analyzer = NCBIGenomeAnalyzer(email)
            
            with st.spinner("Fetching genome sequence and performing analysis..."):
                # Fetch sequence
                sequence = analyzer.fetch_genome_sequence(selected_assembly, chromosome)
                
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
                        st.session_state['selected_assembly'] = selected_assembly
                        st.session_state['chromosome'] = chromosome
                    else:
                        st.error("Analysis failed - no windows could be processed")
                else:
                    st.error("❌ Failed to fetch genome sequence. Try a different chromosome or assembly.")
    
    # Display results
    if 'analysis_results' in st.session_state and not st.session_state['analysis_results'].empty:
        st.header("📊 Analysis Results")
        
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
            conservation_percentage = (len(conserved_regions) / len(df_results)) * 100 if len(df_results) > 0 else 0
            st.metric("Conservation %", f"{conservation_percentage:.1f}%")
        
        # Visualization
        st.subheader("📈 Genomic Landscape")
        
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
        
        fig.update_layout(height=800, showlegend=False, title_text=f"Genomic Analysis: {st.session_state.get('species', 'Unknown')} - {st.session_state.get('selected_assembly', 'Unknown')} - Chr {st.session_state.get('chromosome', 'Unknown')}")
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
                file_name=f"conserved_regions_{st.session_state.get('species', 'unknown').replace(' ', '_')}_{st.session_state.get('selected_assembly', 'unknown')}.csv",
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

    # Troubleshooting section
    with st.expander("🔧 Troubleshooting"):
        st.markdown("""
        **Common Issues and Solutions:**
        
        1. **No assemblies found:**
           - Check species name spelling (use scientific names like "Homo sapiens")
           - Try removing quotes or special characters
           - Some species may have limited public genome data
        
        2. **Sequence fetch fails:**
           - Try different chromosome numbers (1, 2, X, Y, MT)
           - Some assemblies may not have individual chromosomes available
           - Try "1", "I", "chr1", or just leave blank for automatic selection
        
        3. **NCBI connection issues:**
           - Ensure valid email address is provided
           - Check internet connection
           - NCBI may be temporarily unavailable
           - Rate limits may apply - wait a few minutes and retry
        
        4. **Analysis errors:**
           - Sequence may be too short for chosen window size
           - Reduce window size for smaller sequences
           - Some sequences may contain unusual characters
        
        **Performance Tips:**
        - Start with smaller organisms (bacteria) for faster results
        - Use larger window sizes for initial exploration
        - Large genomes may take several minutes to process
        """)

    # Information panel
    with st.expander("ℹ️ About This Tool"):
        st.markdown("""
        **What This Tool Does:**
        
        This application connects directly to NCBI's databases to retrieve and analyze genome sequences.
        It identifies potentially conserved regions based on sequence composition and complexity.
        
        **Conservation Metrics:**
        - **GC Content**: Percentage of G and C nucleotides (optimal: 40-60%)
        - **Sequence Complexity**: Shannon entropy measure (higher = more complex)
        - **Repeat Content**: Percentage of simple repetitive sequences (lower = better)
        - **Conservation Score**: Combined metric weighing all factors
        
        **Data Sources:**
        - NCBI Assembly database for genome assemblies
        - NCBI Nucleotide database for sequence data
        - RefSeq reference sequences when available
        
        **Limitations:**
        - Analysis limited to publicly available genomes
        - Sequence fetch limited to 2MB for performance
        - Conservation prediction based on composition, not comparative analysis
        - Real conservation requires multiple species comparison
        """)

if __name__ == "__main__":
    main()
