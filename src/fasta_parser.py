class FastaParser:
    def __init__(self, file_path, limit=None):
        """
        Initializes the FastaParser with the given file path and limit.
        Reads the entire file into memory upon initialization.
        
        Parameters:
        - file_path: The path to the FASTA file.
        - limit: The maximum number of sequences to read from the file. 
                 If None, all sequences are read.
        """
        self.file_path = file_path
        self.limit = limit
        self.sequences = []
        
        self._parse_fasta_file()

    def _parse_fasta_file(self):
        """
        Parses the FASTA file, reading all sequences into memory.
        """
        header = None
        seq_lines = []
        iter = 1
        
        with open(self.file_path, 'r') as f:
            for line in f:
                if self.limit and iter >= self.limit:
                    break  # Stop if the limit is reached
                line = line.strip()
                if not line:  # Skip empty lines
                    continue
                if line.startswith('>'):
                    if header is not None:
                        # Add the sequence to the list of sequences
                        iter += 1
                        self.sequences.append(''.join(seq_lines).upper())
                    header = line[1:]  # Update header (unused here)
                    seq_lines = []     # Reset sequence accumulator
                else:
                    seq_lines.append(line)
            if header is not None:
                # Add the last sequence after finishing the file
                iter += 1
                self.sequences.append(''.join(seq_lines).upper())
