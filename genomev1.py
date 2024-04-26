import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog, ttk, scrolledtext
from Bio import SeqIO, Seq, AlignIO, Align
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import pickle
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from Bio.Align.Applications import MuscleCommandline
from io import StringIO
import os
import shutil
from threading import Thread
import string
import panel as pn
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Range1d,Rect, Text
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot

class GenomicApp:
    def __init__(self, master):
        self.master = master
        self.master.title("Genomic Application")

        self.sequence_data = []
        self.loaded_files = []  #loaded_files attribute
        

        # Bioinformatic-style colors
        self.bio_GREEN = "#4CAF50"  # Green
        self.bio_BLUE = "#2196F3"  # Blue
        self.bio_YELLOW = "#FFC107"  # Yellow
        self.bio_WHITE = "#FFFFFF"  # White
        self.bio_GRAY = "#A9A9A9"  # Gray

        # Style
        style = ttk.Style()
        style.configure("TButton", padding=5, font=('Arial', 12))
        style.configure("TFrame", background=self.bio_WHITE)
        style.configure("TLabel", font=('Arial', 12, 'bold'))
        style.configure("TNotebook", background=self.bio_GRAY)

        # Notebook (Tabbed interface)
        self.notebook = ttk.Notebook(self.master)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Load Sequences tab
        self.load_frame = ttk.Frame(self.notebook, style="TFrame")
        self.notebook.add(self.load_frame, text="Load Sequences")

        self.create_load_tab()

        # Display loaded files on the left
        self.loaded_files_listbox = tk.Listbox(self.master, height=5, width=30)
        self.loaded_files_listbox.pack(side=tk.LEFT, fill=tk.BOTH, padx=10, pady=10)
        self.update_loaded_files_listbox()

        proceed_button = tk.Button(self.master, text="Preview", command=self.proceed_with_selected_file)
        proceed_button.pack(side="left")

        # Alignment tab
        self.alignment_frame = ttk.Frame(self.notebook, style="TFrame")
        self.notebook.add(self.alignment_frame, text="Alignment")

        self.create_alignment_tab()

        # Compare Sequences tab
        self.compare_frame = ttk.Frame(self.notebook, style="TFrame")
        self.notebook.add(self.compare_frame, text="Compare Sequences")

        self.create_compare_tab()

        # Plot Data tab
        self.plot_frame = ttk.Frame(self.notebook, style="TFrame")
        self.notebook.add(self.plot_frame, text="Plot Data")

        self.create_plot_tab()

    def create_load_tab(self):
        button_style = {'style': 'TButton', 'width': 20}
        label_style = {'style': 'TLabel'}
        frame_style = {'style': 'TFrame'}

        # Frame
        load_frame = ttk.Frame(self.load_frame, **frame_style)
        load_frame.pack(pady=10)

        # Label
        ttk.Label(load_frame, text="Please follow these tips when loading sequences:", **label_style).pack(pady=5)

        tips_text = (
            "1. File Format: Make sure the file is in FASTA format.\n"
            "2. Minimum and Maximum Capacity: There is no strict limit, but very large files may take longer to load.\n"
            "3. Multiple Files: You can load multiple files by holding down the Ctrl key while selecting files."
        )

        ttk.Label(load_frame, text=tips_text, wraplength=900, **label_style).pack(pady=10)

        # Load Button
        ttk.Button(load_frame, text="Load Sequences", command=self.load_sequences_from_file, **button_style).pack(pady=10)

        # Convert to FASTA Button
        ttk.Button(load_frame, text="Convert to FASTA", command=self.convert_to_fasta, **button_style).pack(pady=10)

        # Tooltip
        ttk.Label(load_frame, text="Tip: Convert Binary text to FASTA", foreground='gray').pack(pady=5)

    def load_sequences(self):
        # Call the correct method to load sequences
        self.load_sequences_from_file()
        
    def convert_to_fasta(self, is_nucleotide=False):
        # Ask the user to enter binary text
        binary_text = simpledialog.askstring("Convert to FASTA", "Enter binary text:")
        if not binary_text:
            return  # User canceled or closed the dialog

        try:
            if is_nucleotide:
                fasta_sequence = self.binary_to_nucleotide_sequence(binary_text)
            else:
                # Convert binary text to ASCII text
                fasta_sequence = ''.join([chr(int(binary_text[i:i+8], 2)) for i in range(0, len(binary_text), 8)])

            # Add the sequence to the sequence data
            self.sequence_data.append(fasta_sequence)

            # Display a message
            messagebox.showinfo("Conversion Complete", "Binary text converted and added to sequences.")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during conversion:\n{str(e)}")

    def binary_to_nucleotide_sequence(self, binary_text):
        nucleotide_mapping = {
            '00': 'A',
            '01': 'C',
            '10': 'G',
            '11': 'T'
        }

        nucleotides = [nucleotide_mapping.get(binary_text[i:i+2], 'N') for i in range(0, len(binary_text), 2)]
        return ''.join(nucleotides)

    def update_loaded_files_listbox(self):
        # Clear existing items in the listbox
        self.loaded_files_listbox.delete(0, tk.END)

        # Insert loaded file names into the listbox
        for file_path in self.loaded_files:
            file_name = os.path.basename(file_path)  # Extract only the file name
            self.loaded_files_listbox.insert(tk.END, file_name)

    def load_sequences_from_file(self):
        # Ask the user to select one or more FASTA files
        file_paths = filedialog.askopenfilenames(title="Select FASTA Files", filetypes=[("FASTA Files", "*.fasta")])
        if not file_paths:
            return  # User canceled or closed the dialog

        try:
            # Load sequences from the selected files
            for file_path in file_paths:
                self.loaded_files.append(file_path)  # Store file path
                with open(file_path, 'r') as file:
                    for record in SeqIO.parse(file, "fasta"):
                        sequence = str(record.seq)
                        # Store sequence
                        self.sequence_data.append(sequence)

            # Display the number of loaded sequences
            messagebox.showinfo("Sequences Loaded", f"{len(self.sequence_data)} sequences loaded successfully.")

            # Update the loaded files listbox
            self.update_loaded_files_listbox()

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during loading sequences:\n{str(e)}")

    def proceed_with_selected_file(self):
        # Get the selected file from the listbox
        selected_index = self.loaded_files_listbox.curselection()
        if not selected_index:
            messagebox.showinfo("No File Selected", "Please select a file from the loaded list.")
            return

        # Retrieve the selected file path
        selected_file_index = selected_index[0]
        selected_file_path = self.loaded_files[selected_file_index]

        # Perform further analysis or processing with the selected file
        self.process_selected_file(selected_file_path)

    def process_selected_file(self, file_path):
        # Here you can implement additional processing based on the selected file
        try:
            # For demonstration, let's display the content of the selected file
            with open(file_path, 'r') as file:
                file_content = file.read()
            messagebox.showinfo("File Content", f"Content of {file_path}:\n{file_content}")

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during processing:\n{str(e)}")

    def create_alignment_tab(self):
        button_style = {'style': 'TButton', 'width': 20}
        label_style = {'style': 'TLabel'}
        frame_style = {'style': 'TFrame'}

        # Frame
        alignment_frame = ttk.Frame(self.alignment_frame, **frame_style)
        alignment_frame.pack(pady=5)

        # Label
        ttk.Label(alignment_frame, text="Sequence Alignment", **label_style).pack(pady=5)

        alignment_button = ttk.Button(alignment_frame, text="Perform Muscle Alignment", command=self.perform_muscle_alignment, **button_style)
        alignment_button.pack(pady=10)


        # Display Aligned Sequences with both vertical and horizontal scrolling
        self.aligned_text = scrolledtext.ScrolledText(alignment_frame, wrap=tk.WORD, font=('Verdana', 12), width=225, height=40)
        self.aligned_text.pack(pady=5)

        # Create horizontal scrollbar
        x_scrollbar = ttk.Scrollbar(alignment_frame, orient=tk.HORIZONTAL, command=self.aligned_text.xview)
        x_scrollbar.pack(side=tk.BOTTOM, fill=tk.X)

        # Create vertical scrollbar
        y_scrollbar = ttk.Scrollbar(alignment_frame, orient=tk.VERTICAL, command=self.aligned_text.yview)
        y_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Attach the scrollbars to the ScrolledText widget
        self.aligned_text.configure(xscrollcommand=x_scrollbar.set, yscrollcommand=y_scrollbar.set)

    def detect_sequence_type(self, sequence):
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

    def perform_muscle_alignment(self):
        if len(self.sequence_data) < 1:
            messagebox.showinfo("Insufficient Sequences", "Please load at least one sequences for alignment.")
            return
        try:
            # Create a temporary input file for MUSCLE
            input_file_path = "temp_input.fasta"
            with open(input_file_path, "w") as temp_input_file:
                temp_input_file.write("\n".join([f">{i}\n{seq}" for i, seq in enumerate(self.sequence_data)]))

            # Specify the MUSCLE command
            self.muscle_executable = r"E:\Genome\muscle3.8.31_i86win32.exe"
            muscle_output_file = "aligned_sequences.fasta"
            muscle_command = [self.muscle_executable, "-in", input_file_path, "-out", muscle_output_file]

            # Run the MUSCLE command using subprocess
            process = subprocess.run(muscle_command, capture_output=True, text=True)

            # Check for errors
            if process.returncode != 0:
                raise Exception(f"Error during MUSCLE alignment:\n{process.stderr}")

            # Read the alignment results from the output file
            with open(muscle_output_file, "r") as output_file:
                aligned_results = output_file.read()

             # Display the alignment with coloring sequence
            self.display_aligned_sequences(aligned_results)

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during MUSCLE alignment:\n{str(e)}")
    
    def display_aligned_sequences(self, aligned_results):
        # Display aligned sequences in the scrolled text widget
        self.aligned_text.delete(1.0, tk.END)  # Clear previous content

        # Parse the MUSCLE alignment results
        alignment = AlignIO.read(StringIO(aligned_results), "fasta")

        # Display detailed information about each aligned sequence
        for record in alignment:
            seq_name = record.id
            start_position = record.annotations.get("start", 1)
            end_position = record.annotations.get("end", len(record.seq))
            sequence = str(record.seq)

            # Determine sequence type
            seq_type = self.detect_sequence_type(sequence)

            # Create a formatted string with sequence details and type
            seq_detail_str = f"Name: {seq_name}\n"
            seq_detail_str += f"Start Position: {start_position}\n"
            seq_detail_str += f"End Position: {end_position}\n"
            seq_detail_str += f"Sequence Type: {seq_type}\n"

            # Insert sequence details into the scrolled text widget
            self.aligned_text.insert(tk.END, seq_detail_str)

            # Insert a new line after the sequence details
            self.aligned_text.insert(tk.END, "\n")

            # Insert the colored sequence separately
            self.insert_colored_sequence(sequence)

            # Insert another new line after the colored sequence
            self.aligned_text.insert(tk.END, "\n\n")

        # Scroll to the end of the text widget
        self.aligned_text.see(tk.END)

    def insert_colored_sequence(self, sequence):
        # Iterate over each character in the sequence
        for char in sequence:
            # Determine the color based on the nucleotide character
            if char in 'Aa':
                color = "green"
            elif char in 'Tt':
                color = "red"
            elif char in 'Gg':
                color = "blue"
            elif char in 'Cc':
                color = "orange"
            elif char in "Uu":
                color = "purple"
            else:
                color = "black"  # Default color for unknown characters

            # Insert the character with the specified color and apply tag
            self.aligned_text.insert(tk.END, char, color)
            self.aligned_text.tag_configure(color, foreground=color)

    def create_compare_tab(self):
        button_style = {'style': 'TButton', 'width': 20}
        label_style = {'style': 'TLabel'}
        frame_style = {'style': 'TFrame'}

        # Frame
        compare_frame = ttk.Frame(self.compare_frame, **frame_style)
        compare_frame.pack(pady=10)

        # Label
        ttk.Label(compare_frame, text="Compare Sequences", **label_style).pack(pady=5)

        # Compare Button
        ttk.Button(compare_frame, text="Compare Sequences", command=self.compare_sequences, **button_style).pack(pady=10)

    def compare_sequences(self):
        if len(self.sequence_data) < 1:
            messagebox.showinfo("Insufficient Sequences", "Please load at least one sequences for comparison.")
            return

        # Select a FASTA file for comparison
        file_path = filedialog.askopenfilename(title="Select FASTA File for Comparison", filetypes=[("FASTA Files", "*.fasta")])
        if not file_path:
            return  # User canceled or closed the dialog

        # Load sequences from the selected file
        try:
            with open(file_path, 'r') as file:
                sequences_to_compare = [str(record.seq) for record in SeqIO.parse(file, "fasta")]
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during loading sequences:\n{str(e)}")
            return

        # Check if there are at least two sequences for comparison
        if len(sequences_to_compare) < 2:
            messagebox.showinfo("Insufficient Sequences", "Please provide a FASTA file with at least two sequences for comparison.")
            return

        # Get the first two sequences for detailed comparison
        seq1 = self.sequence_data[0]
        seq2 = sequences_to_compare[0]

        # Compare sequences and identify differences
        differences = []
        for i, (nuc1, nuc2) in enumerate(zip(seq1, seq2), start=1):
            if nuc1 != nuc2:
                differences.append(f"At position {i}: {nuc1} (Reference) != {nuc2} (Comparison)")

        # Display the differences with colored labels
        if differences:
            self.display_sequence_differences(seq1, seq2, differences)
        else:
            messagebox.showinfo("No Differences", "The sequences are identical.")

    def display_sequence_differences(self, ref_seq, comp_seq, differences):
        # Create a new window to display the differences
        diff_window = tk.Toplevel(self.master)
        diff_window.title("Sequence Differences")

        # Create a scrolled text widget for displaying the differences
        scrolled_text = scrolledtext.ScrolledText(diff_window, wrap=tk.WORD, font=('Verdana', 10))
        scrolled_text.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

        # Define color mappings for each nucleotide
        color_mappings = {'A': 'green', 'T': 'red', 'C': 'blue', 'G': 'orange', 'U': 'purple'}

        # Insert headers with addresses
        #scrolled_text.insert(tk.END, "Address   |   ")
        #for i in range(1, len(ref_seq) + 1):
        #    scrolled_text.insert(tk.END, f"{i}    ")
        #scrolled_text.insert(tk.END, "\n" + "-" * 70 + "\n")

        # Insert nucleotides with colored text into the scrolled text widget
        for i, (nuc_ref, nuc_comp) in enumerate(zip(ref_seq, comp_seq), start=1):
            label_text = nuc_ref
            color = color_mappings.get(nuc_ref, 'black')  # Default to black for unknown nucleotides
            tag_name = f"tag_{i}"

            # Highlight differences with a background color
            if nuc_ref != nuc_comp:
                label_text = f"[{label_text}]"  # Highlight differences in square brackets
                color = 'Purple'  # Background color for differences

            # Insert the nucleotide with the specified color
            scrolled_text.insert(tk.END, f"{i}   |   {label_text}    ", tag_name)
            scrolled_text.tag_configure(tag_name, foreground=color, background='lightgray')  # Set background color

        # Append the differences below the nucleotides
        scrolled_text.insert(tk.END, "\n\n" + "\n".join(differences), 'tag_differences')
        scrolled_text.tag_configure('tag_differences', foreground='black', background='white')  # Set background color for differences

        # Make the scrolled text widget expand to fill the available space
        diff_window.grid_rowconfigure(0, weight=1)
        diff_window.grid_columnconfigure(0, weight=1)


    def create_plot_tab(self):
        button_style = {'style': 'TButton', 'width': 20}
        label_style = {'style': 'TLabel'}
        frame_style = {'style': 'TFrame'}

        # Frame
        plot_frame = ttk.Frame(self.plot_frame, **frame_style)
        plot_frame.pack(pady=10)

        # Label
        ttk.Label(plot_frame, text="Plot Data", **label_style).pack(pady=5)

        # Plot Button
        ttk.Button(plot_frame, text="Plot Data", command=self.plot_data, **button_style).pack(pady=10)

    def plot_data(self):
        if not self.sequence_data:
            messagebox.showinfo("No Sequences", "Please load sequences before plotting.")
            return

        try:
            # Ask the user for plot options
            plot_options = ["Bar Plot", "Histogram", "Custom Sequence Plot"]
            selected_option = simpledialog.askstring("Plot Data", "Select plot option:", options=plot_options)
            if selected_option is None:
                return  # User canceled or closed the dialog

            if selected_option == "Custom Sequence Plot":
                selected_sequences = self.select_sequences_for_plot()
                if not selected_sequences:
                    return  # User canceled or closed the dialog
            else:
                selected_sequences = self.sequence_data

            # Ask the user for plot customization options
            plot_color = simpledialog.askstring("Plot Data", "Enter plot color (e.g., 'red', '#00FF00'):")
            if plot_color is None:
                return  # User canceled or closed the dialog

            # Example: Plot the length of sequences
            sequence_lengths = [len(seq) for seq in selected_sequences]

            # Plotting
            fig, ax = plt.subplots()

            if selected_option == "Bar Plot":
                ax.bar(range(1, len(selected_sequences) + 1), sequence_lengths, color=plot_color)
            elif selected_option == "Histogram":
                ax.hist(sequence_lengths, bins=20, color=plot_color)
                ax.set_xlabel("Sequence Length")
                ax.set_ylabel("Frequency")
                ax.set_title("Sequence Length Distribution")

            ax.set_xlabel("Sequence")
            ax.set_ylabel("Length")
            ax.set_title("Sequence Lengths")

            # Embed the plot in Tkinter
            canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
            canvas.draw()
            canvas.get_tk_widget().pack()

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during plotting:\n{str(e)}")

    def select_sequences_for_plot(self):
        # Allow the user to select specific sequences for plotting
        selected_sequences = []
        try:
            selected_indices = simpledialog.askstring("Select Sequences", "Enter indices of sequences to plot (comma-separated):")
            if selected_indices is None:
                return None  # User canceled or closed the dialog

            indices = [int(index.strip()) for index in selected_indices.split(",")]
            selected_sequences = [self.sequence_data[index - 1] for index in indices]
        except Exception as e:
            messagebox.showerror("Error", f"Invalid input for sequence indices:\n{str(e)}")

        return selected_sequences



    def save_sequences(self):
        if not self.sequence_data:
            messagebox.showinfo("No Sequences", "No sequences to save.")
            return

        # Ask the user for the output file
        file_path = filedialog.asksaveasfilename(title="Save Sequences", filetypes=[("Text Files", "*.txt")])
        if not file_path:
            return  # User canceled or closed the dialog

        try:
            # Save sequences to the selected file
            with open(file_path, 'w') as file:
                for i, sequence in enumerate(self.sequence_data, start=1):
                    file.write(f'>Sequence{i}\n{sequence}\n')

            messagebox.showinfo("Save Complete", f"Sequences saved to: {file_path}")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during saving:\n{str(e)}")

    def load_pickle(self):
        # Ask the user to select a pickle file
        file_path = filedialog.askopenfilename(title="Load Pickle File", filetypes=[("Pickle Files", "*.pkl")])
        if not file_path:
            return  # User canceled or closed the dialog

        try:
            # Load sequences from the pickle file
            with open(file_path, 'rb') as file:
                self.sequence_data = pickle.load(file)

            # Update the listbox with the loaded sequences
            self.listbox.delete(0, tk.END)
            for sequence in self.sequence_data:
                self.listbox.insert(tk.END, sequence)

            print("Sequences loaded from pickle file.")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during loading pickle file:\n{str(e)}")

    def search_sequences(self):
        try:
            # Ask the user for the search term and search options
            search_term = simpledialog.askstring("Search Sequences", "Enter a search term:")
            if search_term is None:
                return  # User canceled or closed the dialog

            search_options = ["Search Anywhere", "Search at Beginning", "Search at End", "Search at Specific Address", "Search for Specific Pattern"]
            selected_option = simpledialog.askitem("Search Sequences", "Select search option:", search_options)
            if selected_option is None:
                return  # User canceled or closed the dialog

            # Perform the search based on the selected option
            search_results = []

            for index, sequence in enumerate(self.sequence_data):
                if selected_option == "Search Anywhere" and search_term.lower() in sequence.lower():
                    search_results.append((index, sequence))
                elif selected_option == "Search at Beginning" and sequence.lower().startswith(search_term.lower()):
                    search_results.append((index, sequence))
                elif selected_option == "Search at End" and sequence.lower().endswith(search_term.lower()):
                    search_results.append((index, sequence))
                elif selected_option == "Search at Specific Address":
                    address = simpledialog.askinteger("Search Sequences", "Enter the specific address (1-based index):")
                    if address is None:
                        return  # User canceled or closed the dialog
                    if 1 <= address <= len(sequence) and sequence[address - 1:address - 1 + len(search_term)].lower() == search_term.lower():
                        search_results.append((index, sequence))
                elif selected_option == "Search for Specific Pattern":
                    pattern_options = ["Exact Match", "Regular Expression"]
                    selected_pattern_option = simpledialog.askitem("Search Sequences", "Select pattern matching option:", pattern_options)
                    if selected_pattern_option is None:
                        return  # User canceled or closed the dialog

                    if selected_pattern_option == "Exact Match" and search_term.lower() == sequence.lower():
                        search_results.append((index, sequence))
                    elif selected_pattern_option == "Regular Expression":
                        import re
                        if re.search(search_term, sequence, flags=re.IGNORECASE):
                            search_results.append((index, sequence))

            # Display the search results
            if search_results:
                result_message = f"Found {len(search_results)} sequence(s) matching the search term."
                for index, sequence in search_results:
                    result_message += f"\n\nIndex: {index + 1}\nSequence: {sequence}"
                messagebox.showinfo("Search Results", result_message)
            else:
                messagebox.showinfo("No Matches", "No sequences found matching the search term.")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during searching:\n{str(e)}")

    def display_sequences(self, sequences):
        try:
            # Ask the user for display options
            display_options = ["Highlight Specific Pattern", "Navigate to Specific Address"]
            selected_option = simpledialog.askitem("Display Sequences", "Select display option:", display_options)
            if selected_option is None:
                return  # User canceled or closed the dialog

            # Create a scrolled text widget to display sequences with rich text formatting
            text_widget = scrolledtext.ScrolledText(self.master, wrap=tk.WORD, width=80, height=20, font=('Arial', 12))
            text_widget.pack(pady=10)

            for sequence_index, sequence in enumerate(sequences, start=1):
                text_widget.insert(tk.END, f"Sequence {sequence_index}:\n{sequence}\n\n")

            if selected_option == "Highlight Specific Pattern":
                pattern_to_highlight = simpledialog.askstring("Display Sequences", "Enter pattern to highlight:")
                if pattern_to_highlight is None:
                    return  # User canceled or closed the dialog

                for sequence_index, sequence in enumerate(sequences, start=1):
                    self.highlight_pattern(text_widget, sequence, pattern_to_highlight, sequence_index)
            elif selected_option == "Navigate to Specific Address":
                address_to_navigate = simpledialog.askinteger("Display Sequences", "Enter address to navigate to (1-based index):")
                if address_to_navigate is None:
                    return  # User canceled or closed the dialog

                for sequence_index, sequence in enumerate(sequences, start=1):
                    self.display_sequence_around_address(text_widget, sequence, address_to_navigate, sequence_index)

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during displaying sequences:\n{str(e)}")

    def highlight_pattern(self, text_widget, sequence, pattern, sequence_index):
        import re

        # Use regular expression to find and highlight the pattern
        matches = re.finditer(pattern, sequence, flags=re.IGNORECASE)

        for match in matches:
            start, end = match.span()
            text_widget.tag_add(f"sequence_{sequence_index}", f"{text_widget.index(tk.END)}-{len(sequence)}c")
            text_widget.tag_config(f"sequence_{sequence_index}", background="yellow")
            text_widget.tag_add(f"sequence_{sequence_index}", f"{text_widget.index(tk.END)}-{len(sequence)}c")

    def display_sequence_around_address(self, text_widget, sequence, address, sequence_index):
        # Display a portion of the sequence around the specified address
        start_index = max(0, address - 10)
        end_index = min(len(sequence), address + 10)

        displayed_sequence = sequence[start_index:end_index]
        text_widget.insert(tk.END, f"Sequence {sequence_index}, Navigated to address {address}:\n{displayed_sequence}\n\n")


# Create the main window
root = tk.Tk()
app = GenomicApp(root)
root.mainloop()
