import io
import panel as pn
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from bokeh.models import ColumnDataSource, Rect, Text, TapTool, CustomJS
from bokeh.models import HoverTool
from bokeh.layouts import gridplot
from bokeh.plotting import figure, show
import subprocess
import numpy as np
import tkinter as tk
from sklearn.linear_model import LinearRegression

# Define functions for sequence alignment and visualization
def muscle_alignment(seqs):
    filename = 'temp.faa'
    with open(filename, 'w') as f:
        SeqIO.write(seqs, f, "fasta")
    muscle_executable = r"muscle3.8.31_i86win32.exe"
    cmd = f"{muscle_executable} -in {filename} -out temp.txt"
    subprocess.run(cmd, shell=True, check=True)
    align = list(SeqIO.parse('temp.txt', 'fasta'))
    return align

def get_colors(sequence):
    # Define a color mapping for nucleotide characters
    color_map = {'A': 'red', 'T': 'green', 'G': 'orange', 'C': 'blue', '-': 'white'}
    
    # Map each character in the sequence to its corresponding color
    colors = [color_map[char] for char in sequence]
    
    return colors
def view_alignment(aln, fontsize="9pt", plot_width=2000, max_positions=2000):
    seqs = [rec.seq for rec in aln]
    ids = [rec.id for rec in aln]
    N = len(seqs[0])
    S = len(seqs)

    # Limit the number of positions to display
    display_range = min(N, max_positions)

    # Prepare data for rendering
    text = []
    colors = []
    x = []
    y = []

    for i, seq in enumerate(seqs):
        # Simplify sequence for display and limit to display range
        simplified_seq = simplify_sequence(seq[:display_range])
        text.extend(simplified_seq)
        colors.extend(get_colors(simplified_seq))

        # Position data for glyphs
        x.extend(list(range(1, len(simplified_seq) + 1)))
        y.extend([i] * len(simplified_seq))

    # Create data source for plotting
    source = ColumnDataSource(dict(x=x, y=y, text=text, colors=colors))

    plot_height = S * 15 + 50

    # Create main figure (overview)
    p = figure(title=None, width=plot_width, height=50, x_range=(0, display_range + 1),
               y_range=(-1, S), tools="xpan, xwheel_zoom, reset, save", min_border=0, toolbar_location='below')
    rects = Rect(x="x", y="y", width=1, height=1, fill_color="colors", line_color=None, fill_alpha=0.6)
    p.add_glyph(source, rects)
    p.yaxis.visible = True
    p.grid.visible = True

    # Create detailed figure (zoomed-in)
    p1 = figure(title=None, width=plot_width, height=plot_height, x_range=(0, min(display_range, 100)),
                y_range=ids, tools="xpan, reset", min_border=0, toolbar_location='below')
    glyph = Text(x="x", y="y", text="text", text_align='center', text_color="black",
                 text_font="monospace", text_font_size=fontsize)
    p1.add_glyph(source, glyph)

    # Add rectangles as glyphs for the detailed figure
    rects = Rect(x="x", y="y", width=1, height=1, fill_color="colors", line_color=None, fill_alpha=0.4)
    p1.add_glyph(source, rects)

    p1.grid.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0

    # Combine plots into a grid layout
    p = gridplot([[p], [p1]], toolbar_location='below')
    return p

def simplify_sequence(seq):
    """Simplify a sequence by grouping consecutive identical nucleotides."""
    simplified_seq = []
    current_char = ''
    for char in seq:
        if char != current_char:
            simplified_seq.append(char)
            current_char = char
    return simplified_seq

def detect_sequence_type(sequence):
    # Convert the sequence to uppercase for case-insensitive comparisons
    sequence = sequence.upper()

    # Define valid nucleotide and amino acid characters
    nucleotide_bases = {'A', 'C', 'G', 'T', 'U'}
    amino_acids = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}

    # Count occurrences of each nucleotide and amino acid
    count_nucleotides = sum(1 for base in sequence if base in nucleotide_bases)
    count_amino_acids = sum(1 for aa in sequence if aa in amino_acids)

    # Calculate total length of the sequence
    total_length = len(sequence)

    # Calculate percentage of nucleotides and amino acids in the sequence
    nucleotide_percentage = count_nucleotides / total_length
    amino_acid_percentage = count_amino_acids / total_length

    # Determine sequence type based on composition and additional checks
    if total_length < 10:
        return 'Too Short'

    # Determine sequence type based on composition
    if nucleotide_percentage > 0.50:  # If >50% of characters are nucleotides
        if 'U' in sequence:
            return 'RNA'  # Likely RNA sequence
        else:
            return 'DNA'  # Likely DNA sequence
    elif amino_acid_percentage > 0.90:  # If >90% of characters are amino acids
        # Perform additional checks for protein-like characteristics
        has_start_codon = sequence.startswith('ATG')
        has_stop_codon = any(sequence.endswith(stop) for stop in ['TAA', 'TAG', 'TGA'])
        contains_valid_amino_acids = all(aa in amino_acids for aa in sequence)

        if has_start_codon and has_stop_codon and contains_valid_amino_acids:
            return 'Protein'  # Likely protein sequence with start and stop codons
        else:
            return 'Peptide'  # Likely short amino acid sequence (peptide) or non-standard protein
    else:
        return 'Unknown'  # If neither RNA, DNA, nor Protein criteria are met, classify as Unknown

# Define mutation analysis function
def plot_mutation_analysis(seq1, seq2):
    mutation_positions = []
    mutation_types = []

    for i, (char1, char2) in enumerate(zip(seq1, seq2), start=1):
        if char1 != char2:
            mutation_positions.append(i)
            mutation_types.append('Mismatch')
        else:
            mutation_types.append('Match')

    source = ColumnDataSource(data=dict(x=mutation_positions, y=[0]*len(mutation_positions), 
                                        mutation_types=mutation_types))

    p = figure(title="Mutation Analysis", plot_width=800, plot_height=300,
               x_range=(0, len(seq1) + 1), y_range=(-1, 1),
               toolbar_location=None, tools="")

    p.circle(x='x', y='y', size=10, source=source, color='red', legend_field='mutation_types')

    p.title.text_font_size = '14pt'
    p.xaxis.axis_label = 'Position'
    p.yaxis.visible = False
    p.grid.visible = False

    p.legend.title = "Mutation Type"
    p.legend.location = "top_left"
    p.legend.label_text_font_size = "10pt"

    return p

    
# Define event handler for upload button
def upload_and_display(event):
    try:
        all_sequences = []
        ids = []

        for file_contents in file_input.value:
            file_text = file_contents.decode("utf-8")
            sequences = list(SeqIO.parse(io.StringIO(file_text), "fasta"))
            all_sequences.extend(sequences)
            ids.extend([rec.id for rec in sequences])

        if all_sequences:
            seq_pane.object = f"Loaded {len(all_sequences)} sequences from {len(file_input.value)} files."

            aln = muscle_alignment(all_sequences)
            bokeh_plot = view_alignment(aln, fontsize="7pt", plot_width=1000)
            bokeh_pane.object = bokeh_plot

            seq_types = [detect_sequence_type(str(rec.seq)) for rec in all_sequences]
            seq_types_summary = {}
            for seq_id, seq_type in zip(ids, seq_types):
                if seq_type in seq_types_summary:
                    seq_types_summary[seq_type].append(seq_id)
                else:
                    seq_types_summary[seq_type] = [seq_id]

            seq_types_str = "\n".join(f"- {seq_type}: {', '.join(seq_ids)}" for seq_type, seq_ids in seq_types_summary.items())
            seq_types_pane.object = f"Sequence Types Summary:\n{seq_types_str}"

    except Exception as e:
        seq_pane.object = f"Error loading sequences: {str(e)}"

def compare_sequences_callback(event):
    if len(file_input.value) < 2:
        comparison_pane.object = "Please upload at least two sequences for comparison."
        return

    try:
        seq1 = str(SeqIO.parse(io.StringIO(file_input.value[0].decode("utf-8")), "fasta").__next__().seq)
        seq2 = str(SeqIO.parse(io.StringIO(file_input.value[1].decode("utf-8")), "fasta").__next__().seq)
        
        differences, same_positions, similarity_percentage = compare_sequences_helper(seq1, seq2)
        
        diff_str = "\n".join(differences)
        same_str = "\n".join(same_positions)
        
        comparison_pane.object = f"Sequence Differences:\n{diff_str}\n\nSame Positions:\n{same_str}\n\nSimilarity Percentage: {similarity_percentage:.2f}%"

    except Exception as e:
        comparison_pane.object = f"Error comparing sequences: {str(e)}"




def compare_sequences_helper(seq1, seq2):
    differences = []
    same_positions = []
    identical_count = 0
    total_positions = 0

    for i, (char1, char2) in enumerate(zip(seq1, seq2), start=1):
        total_positions += 1
        if char1 != char2:
            differences.append(f"At position {i}: '{char1}' (Sequence 1) != '{char2}' (Sequence 2)")
        else:
            same_positions.append(f"At position {i}: '{char1}'")
            identical_count += 1

    # Calculate percentage similarity
    similarity_percentage = (identical_count / total_positions) * 100

    return differences, same_positions, similarity_percentage

def plot_mutation_analysis(seq1, seq2):
    mutation_positions = []
    mutation_types = []

    for i, (char1, char2) in enumerate(zip(seq1, seq2), start=1):
        if char1 != char2:
            mutation_positions.append(i)
            mutation_types.append('Mismatch')

    source = ColumnDataSource(data=dict(x=mutation_positions, y=[0]*len(mutation_positions), 
                                        mutation_types=mutation_types))

    p = figure(title="Mutation Analysis", width=800, height=600,
               x_axis_label='Amino Acid Substitutions per Site',
               y_axis_label='Disease Mutations per Replacement Site',
               toolbar_location=None, tools="")

    # Calculate observed values (for demonstration only)
    observed_x = [0.5, 1.2, 2.0, 3.1, 4.5]
    observed_y = [2, 3, 4, 5, 6]

    # Plot observed data points using scatter() instead of circle()
    p.scatter(x=observed_x, y=observed_y, size=8, color='navy', alpha=0.6, legend_label='Observed Data')

    # Plot expected regression lines
    model = LinearRegression()
    model.fit(np.array(observed_x).reshape(-1, 1), observed_y)
    regression_line = model.predict(np.array(observed_x).reshape(-1, 1))

    p.line(observed_x, regression_line, line_width=2, color='red', legend_label='Regression Line')

    # Plot expected lines (uniform and evolutionarily-influenced distributions)
    uniform_expected = np.linspace(min(observed_x), max(observed_x), 100)
    evolution_expected = 2 * uniform_expected

    p.line(uniform_expected, uniform_expected, line_dash='dashed', line_width=2, color='green', legend_label='Uniform Expected')
    p.line(uniform_expected, evolution_expected, line_dash='dotted', line_width=2, color='orange', legend_label='Evolution Expected')


    p.title.text_font_size = '14pt'
    p.legend.location = "top_left"
    p.legend.label_text_font_size = "10pt"

    return p


# Define event handler for mutation button
def perform_mutation_analysis(event):
    if len(file_input.value) < 2:
        mutation_pane.object = "Please upload at least two sequences for mutation analysis."
        return

    try:
        seq1 = str(SeqIO.parse(io.StringIO(file_input.value[0].decode("utf-8")), "fasta").__next__().seq)
        seq2 = str(SeqIO.parse(io.StringIO(file_input.value[1].decode("utf-8")), "fasta").__next__().seq)
        mutation_plot = plot_mutation_analysis(seq1, seq2)

        # Update the mutation_pane with the Bokeh plot object
        mutation_pane.object = mutation_plot

    except Exception as e:
        # Handle errors by displaying an error message
        mutation_pane.object = f"Error performing mutation analysis: {str(e)}"

# Define Panel app components
file_input = pn.widgets.FileInput(accept=".fasta", multiple=True)
upload_btn = pn.widgets.Button(name='Upload & Display', button_type='primary')
seq_pane = pn.pane.Str()
seq_types_pane = pn.pane.Str()
bokeh_pane = pn.pane.Bokeh()
comparison_pane = pn.pane.Str()
mutation_pane = pn.pane.Bokeh()
compare_btn = pn.widgets.Button(name='Compare Sequences', button_type='primary')
mutation_btn = pn.widgets.Button(name='Perform Mutation Analysis', button_type='primary')

# Assign event handlers to buttons
upload_btn.on_click(upload_and_display)
compare_btn.on_click(compare_sequences_callback)
mutation_btn.on_click(perform_mutation_analysis)

# Create Panel app layout
upload_section = pn.Column(file_input, upload_btn, seq_pane, seq_types_pane, bokeh_pane)
compare_section = pn.Column(compare_btn, comparison_pane)
mutation_section = pn.Column(file_input, mutation_btn, mutation_pane)

# Display the app layout
app = pn.Tabs(
    ('Upload & Display', upload_section),
    ('Sequence Comparison', compare_section),
    ('Mutation Analysis', mutation_section)
)


# Display the app
app.show()
