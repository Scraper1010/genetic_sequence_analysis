from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                            QHBoxLayout, QTextEdit, QLabel, QPushButton, 
                            QScrollArea, QFrame, QFileDialog, QSplitter,
                            QProgressBar)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt5.QtGui import QFont, QPalette, QColor, QIcon
from Bio import SeqIO
from collections import Counter
import sys
import time
from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import QMimeData

class DarkTheme:
    @staticmethod
    def apply(app):
        app.setStyle("Fusion")
        palette = QPalette()
        palette.setColor(QPalette.Window, QColor(53, 53, 53))
        palette.setColor(QPalette.WindowText, Qt.white)
        palette.setColor(QPalette.Base, QColor(25, 25, 25))
        palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
        palette.setColor(QPalette.ToolTipBase, Qt.white)
        palette.setColor(QPalette.ToolTipText, Qt.white)
        palette.setColor(QPalette.Text, Qt.white)
        palette.setColor(QPalette.Button, QColor(53, 53, 53))
        palette.setColor(QPalette.ButtonText, Qt.white)
        palette.setColor(QPalette.BrightText, Qt.red)
        palette.setColor(QPalette.Link, QColor(42, 130, 218))
        palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
        palette.setColor(QPalette.HighlightedText, Qt.black)
        app.setPalette(palette)

class AnalysisWorker(QThread):
    finished = pyqtSignal(dict)
    progress = pyqtSignal(int)
    
    def __init__(self, sequence):
        super().__init__()
        self.sequence = sequence
    
    def run(self):
        try:
            # Basic sequence analysis
            self.progress.emit(10)
            length = len(self.sequence)
            
            # Detect if sequence is DNA or RNA
            is_rna = 'U' in self.sequence
            sequence_type = "RNA" if is_rna else "DNA"
            
            self.progress.emit(20)
            gc_content = (self.sequence.count('G') + self.sequence.count('C')) / length * 100
            
            self.progress.emit(30)
            # Nucleotide frequencies
            nucleotide_counts = Counter(self.sequence)
            nucleotide_percentages = {nt: (count / length) * 100 for nt, count in nucleotide_counts.items()}
            
            self.progress.emit(50)
            # Codon frequencies
            codons = [self.sequence[i:i+3] for i in range(0, len(self.sequence) - 2, 3)]
            codon_counts = Counter(codons)
            
            self.progress.emit(70)
            # AT/GC Skew
            at_skew = (self.sequence.count('A') - self.sequence.count('T')) / (self.sequence.count('A') + self.sequence.count('T') + 1e-6)
            gc_skew = (self.sequence.count('G') - self.sequence.count('C')) / (self.sequence.count('G') + self.sequence.count('C') + 1e-6)
            
            self.progress.emit(90)
            # Count N's in sequence
            n_count = self.sequence.count('N')
            
            # Filter out N's from sequence
            filtered_sequence = ''.join(base for base in self.sequence if base != 'N')
            
            # Generate transformations
            complement = self.sequence.translate(str.maketrans('ATGCU', 'TACGA' if not is_rna else 'UACGA'))
            reverse_complement = complement[::-1]
            transcription = self.sequence.replace('T', 'U') if not is_rna else self.sequence
            
            results = {
                'length': length,
                'gc_content': gc_content,
                'at_skew': at_skew,
                'gc_skew': gc_skew,
                'nucleotide_counts': nucleotide_counts,
                'nucleotide_percentages': nucleotide_percentages,
                'codon_counts': codon_counts,
                'sequence': self.sequence,
                'filtered_sequence': filtered_sequence,
                'n_count': n_count,
                'sequence_type': sequence_type,
                'complement': complement,
                'reverse_complement': reverse_complement,
                'transcription': transcription
            }
            
            self.progress.emit(100)
            self.finished.emit(results)
            
        except Exception as e:
            self.finished.emit({'error': str(e)})

class ResultFrame(QFrame):
    def __init__(self, title, content, show_copy_button=True):
        super().__init__()
        self.setFrameShape(QFrame.StyledPanel)
        self.setStyleSheet("""
            QFrame {
                background-color: #2d2d2d;
                border-radius: 5px;
                padding: 15px;
                margin: 5px;
            }
        """)
        layout = QVBoxLayout()
        
        # Title and copy button layout
        title_layout = QHBoxLayout()
        title_label = QLabel(title)
        title_label.setStyleSheet("""
            color: #42b0f5;
            font-weight: bold;
            font-size: 18px;
        """)
        
        title_layout.addWidget(title_label)
        
        if show_copy_button:
            # Create a container for the copy button and feedback
            copy_container = QWidget()
            copy_container.setStyleSheet("""
                QWidget {
                    background-color: #3d3d3d;
                    border-radius: 15px;
                    padding: 2px;
                }
            """)
            copy_container_layout = QHBoxLayout(copy_container)
            copy_container_layout.setContentsMargins(0, 0, 0, 0)
            
            self.copy_feedback = QLabel("")
            self.copy_feedback.setStyleSheet("""
                color: #4CAF50;
                font-size: 16px;
                font-weight: bold;
            """)
            self.copy_feedback.setFixedWidth(100)
            
            copy_button = QPushButton("📋")
            copy_button.setStyleSheet("""
                QPushButton {
                    background-color: #3d3d3d;
                    color: #42b0f5;
                    border: none;
                    font-size: 24px;
                    padding: 4px;
                    border-radius: 13px;
                }
                QPushButton:hover {
                    background-color: #4d4d4d;
                    color: #3a9ee0;
                }
                QPushButton:pressed {
                    background-color: #5d5d5d;
                    color: #2a8ed0;
                }
            """)
            copy_button.setFixedSize(40, 40)
            copy_button.setCursor(Qt.PointingHandCursor)
            copy_button.clicked.connect(lambda: self.copy_to_clipboard(content))
            
            copy_container_layout.addWidget(self.copy_feedback)
            copy_container_layout.addWidget(copy_button)
            copy_container.setFixedSize(150, 40)
            
            title_layout.addStretch()
            title_layout.addWidget(copy_container)
        
        content_label = QLabel(content)
        content_label.setStyleSheet("""
            color: white;
            font-size: 16px;
        """)
        content_label.setWordWrap(True)
        
        layout.addLayout(title_layout)
        layout.addWidget(content_label)
        self.setLayout(layout)
    
    def copy_to_clipboard(self, text):
        clipboard = QApplication.clipboard()
        clipboard.setText(text)
        self.copy_feedback.setText("Copied!")
        QTimer.singleShot(2000, lambda: self.copy_feedback.setText(""))

class FileLoader(QThread):
    finished = pyqtSignal(str)
    error = pyqtSignal(str)
    
    def __init__(self, file_name):
        super().__init__()
        self.file_name = file_name
    
    def run(self):
        try:
            record = SeqIO.read(self.file_name, "fasta")
            self.finished.emit(str(record.seq))
        except Exception as e:
            self.error.emit(str(e))

class GeneticAnalyzerGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("🧬 Genetic Sequence Analyzer")
        self.setMinimumSize(1200, 900)
        self.current_file = None
        self.sequence = None
        
        # Main widget and layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        layout = QVBoxLayout(main_widget)
        
        # Create splitter for input and results
        splitter = QSplitter(Qt.Vertical)
        
        # Input section
        input_frame = QFrame()
        input_frame.setStyleSheet("""
            QFrame {
                background-color: #2d2d2d;
                border-radius: 5px;
                padding: 15px;
            }
        """)
        input_layout = QVBoxLayout()
        
        # File input section
        file_layout = QHBoxLayout()
        input_label = QLabel("Enter DNA Sequence:")
        input_label.setStyleSheet("""
            color: #42b0f5;
            font-weight: bold;
            font-size: 18px;
        """)
        
        self.file_path_label = QLabel("No file selected")
        self.file_path_label.setStyleSheet("""
            color: #888888;
            font-size: 14px;
        """)
        self.file_path_label.setWordWrap(True)
        
        self.browse_button = QPushButton("Browse FASTA File")
        self.browse_button.setStyleSheet("""
            QPushButton {
                background-color: #42b0f5;
                color: white;
                border: none;
                border-radius: 3px;
                padding: 8px;
                font-weight: bold;
                font-size: 14px;
                min-height: 40px
            }
            QPushButton:hover {
                background-color: #3a9ee0;
            }
            QPushButton:disabled {
                background-color: #2d2d3d;
                color: #888888;
            }
        """)
        self.browse_button.clicked.connect(self.browse_file)
        
        self.reset_file_button = QPushButton("Reset File")
        self.reset_file_button.setStyleSheet("""
            QPushButton {
                background-color: #ff4444;
                color: white;
                border: none;
                border-radius: 3px;
                padding: 8px;
                font-weight: bold;
                font-size: 14px;
                min-height: 40px
            }
            QPushButton:hover {
                background-color: #ff3333;
            }
            QPushButton:disabled {
                background-color: #2d2d3d;
                color: #888888;
            }
        """)
        self.reset_file_button.clicked.connect(self.reset_file)
        self.reset_file_button.setVisible(False)
        
        file_layout.addWidget(input_label)
        file_layout.addWidget(self.file_path_label)
        file_layout.addWidget(self.browse_button)
        file_layout.addWidget(self.reset_file_button)
        
        self.sequence_input = QTextEdit()
        self.sequence_input.setPlaceholderText("Enter your DNA sequence here or load a FASTA file...")
        self.sequence_input.setStyleSheet("""
            QTextEdit {
                background-color: #1e1e1e;
                color: white;
                border: 1px solid #3d3d3d;
                border-radius: 3px;
                padding: 8px;
                font-size: 16px;
            }
            QTextEdit:focus {
                border: 1px solid #42b0f5;
                color: white;
            }
            QTextEdit:focus {
                border: 1px solid #42b0f5;
                color: white;
            }
        """)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setStyleSheet("""
            QProgressBar {
                border: 2px solid #3d3d3d;
                border-radius: 10px;
                text-align: center;
                background-color: #1e1e1e;
                min-height: 30px;
            }
            QProgressBar::chunk {
                background-color: #42b0f5;
                border-radius: 8px;
                margin: 2px;
            }
        """)
        self.progress_bar.setVisible(False)
        
        self.analyze_button = QPushButton("Analyze Sequence")
        self.analyze_button.setStyleSheet("""
            QPushButton {
                background-color: #42b0f5;
                color: white;
                border: none;
                border-radius: 3px;
                padding: 10px;
                font-weight: bold;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: #3a9ee0;
            }
        """)
        self.analyze_button.clicked.connect(self.analyze_sequence)
        
        input_layout.addLayout(file_layout)
        input_layout.addWidget(self.sequence_input)
        input_layout.addWidget(self.progress_bar)
        input_layout.addWidget(self.analyze_button)
        input_frame.setLayout(input_layout)
        
        # Results section
        self.results_scroll = QScrollArea()
        self.results_scroll.setWidgetResizable(True)
        self.results_scroll.setStyleSheet("""
            QScrollArea {
                border: none;
                background-color: #1e1e1e;
            }
        """)
        
        self.results_widget = QWidget()
        self.results_layout = QVBoxLayout(self.results_widget)
        self.results_scroll.setWidget(self.results_widget)
        
        # Add widgets to splitter
        splitter.addWidget(input_frame)
        splitter.addWidget(self.results_scroll)
        splitter.setSizes([300, 600])
        
        # Add splitter to main layout
        layout.addWidget(splitter)
        
        # Apply dark theme
        DarkTheme.apply(QApplication.instance())
    
    def browse_file(self):
        file_name, _ = QFileDialog.getOpenFileName(
            self,
            "Open FASTA File",
            "",
            "FASTA Files (*.fasta *.fa);;All Files (*)"
        )
        if file_name:
            self.current_file = file_name
            self.file_path_label.setText(f"Selected: {file_name}")
            self.file_path_label.setStyleSheet("""
                color: #42b0f5;
                font-size: 14px;
            """)
            self.browse_button.setText("Change File")
            self.reset_file_button.setVisible(True)
            self.sequence_input.clear()
            self.sequence_input.setPlaceholderText("File selected. Click 'Analyze Sequence' to process.")
            self.sequence_input.setEnabled(False)
            self.sequence_input.setStyleSheet("""
                QTextEdit {
                    background-color: #1e1e1e;
                    color: #888888;
                    border: 1px solid #3d3d3d;
                    border-radius: 3px;
                    padding: 8px;
                    font-size: 16px;
                }
            """)
    
    def reset_file(self):
        self.current_file = None
        self.file_path_label.setText("No file selected")
        self.file_path_label.setStyleSheet("""
            color: #888888;
            font-size: 14px;
        """)
        self.browse_button.setText("Browse FASTA File")
        self.reset_file_button.setVisible(False)
        self.sequence_input.clear()
        self.sequence_input.setEnabled(True)
        self.sequence_input.setPlaceholderText("Enter your DNA sequence here or load a FASTA file...")
        self.sequence_input.setStyleSheet("""
            QTextEdit {
                background-color: #1e1e1e;
                color: white;
                border: 1px solid #3d3d3d;
                border-radius: 3px;
                padding: 8px;
                font-size: 16px;
            }
            QTextEdit:focus {
                border: 1px solid #42b0f5;
                color: white;
            }
            QTextEdit:focus {
                border: 1px solid #42b0f5;
                color: white;
            }
        """)
    
    def analyze_sequence(self):
        # Clear previous results
        for i in reversed(range(self.results_layout.count())): 
            self.results_layout.itemAt(i).widget().setParent(None)
        
        # Get sequence from either file or input
        if self.current_file:
            # Disable all inputs
            self.sequence_input.setEnabled(False)
            self.browse_button.setEnabled(False)
            self.analyze_button.setEnabled(False)
            self.reset_file_button.setEnabled(False)
            self.progress_bar.setVisible(True)
            self.progress_bar.setValue(0)
            
            # Load file in background
            self.file_loader = FileLoader(self.current_file)
            self.file_loader.finished.connect(self.start_analysis)
            self.file_loader.error.connect(self.file_load_error)
            self.file_loader.start()
        else:
            sequence = self.sequence_input.toPlainText().strip()
            if not sequence:
                self.show_error("Please enter a sequence")
                return
            
            # Check for spaces
            if ' ' in sequence:
                self.show_error("Sequence cannot contain spaces")
                return
            
            # Check minimum length
            if len(sequence) < 10:
                self.show_error("Sequence must be at least 10 characters long")
                return
                
            self.start_analysis(sequence.upper())
    
    def show_error(self, message):
        """Show error message in results area"""
        self.results_layout.addWidget(ResultFrame("❌ Error", message, show_copy_button=False))
    
    def start_analysis(self, sequence):
        # Create and start worker thread
        self.worker = AnalysisWorker(sequence)
        self.worker.progress.connect(self.update_progress)
        self.worker.finished.connect(self.show_results)
        self.worker.start()
    
    def file_load_error(self, error_message):
        self.sequence_input.setEnabled(True)
        self.browse_button.setEnabled(True)
        self.analyze_button.setEnabled(True)
        self.progress_bar.setVisible(False)
        
        self.file_path_label.setText("Error loading file")
        self.file_path_label.setStyleSheet("""
            color: #ff4444;
            font-size: 14px;
        """)
        self.browse_button.setText("Browse FASTA File")
        self.current_file = None
    
    def update_progress(self, value):
        # Smooth progress bar animation
        current_value = self.progress_bar.value()
        if value > current_value:
            self.progress_bar.setValue(value)
    
    def show_results(self, results):
        # Re-enable only the reset button
        self.sequence_input.setEnabled(False)
        self.browse_button.setEnabled(False)
        self.reset_file_button.setEnabled(False)
        self.analyze_button.setEnabled(True)
        
        # Update button to reset state
        self.analyze_button.setText("Reset Analysis")
        self.analyze_button.setStyleSheet("""
            QPushButton {
                background-color: #ff4444;
                color: white;
                border: none;
                border-radius: 3px;
                padding: 10px;
                font-weight: bold;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: #ff3333;
            }
        """)
        self.analyze_button.clicked.disconnect()
        self.analyze_button.clicked.connect(self.reset_analysis)
        
        self.progress_bar.setVisible(False)
        
        if 'error' in results:
            self.results_layout.addWidget(ResultFrame("❌ Error", results['error']))
            return
        
        # Store results for report generation
        self.analysis_results = results
        
        # Create button container
        button_container = QWidget()
        button_layout = QHBoxLayout(button_container)
        button_layout.setContentsMargins(0, 0, 0, 15)
        
        # Add save report button
        self.save_button = QPushButton("💾 Save Analysis Report")
        self.save_button.setStyleSheet("""
            QPushButton {
                background-color: #4CAF50;
                color: white;
                border: none;
                border-radius: 3px;
                padding: 10px;
                font-weight: bold;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
        """)
        self.save_button.clicked.connect(self.save_report)
        button_layout.addWidget(self.save_button)
        
        # Add save RNA transcription button
        self.save_rna_button = QPushButton("🧬 Save RNA Transcription")
        self.save_rna_button.setStyleSheet("""
            QPushButton {
                background-color: #2196F3;
                color: white;
                border: none;
                border-radius: 3px;
                padding: 10px;
                font-weight: bold;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: #1976D2;
            }
        """)
        self.save_rna_button.clicked.connect(lambda: self.save_sequence(results['transcription'], "RNA_transcription"))
        button_layout.addWidget(self.save_rna_button)
        
        # Add save complement button
        self.save_complement_button = QPushButton("🔄 Save Complement")
        self.save_complement_button.setStyleSheet("""
            QPushButton {
                background-color: #9C27B0;
                color: white;
                border: none;
                border-radius: 3px;
                padding: 10px;
                font-weight: bold;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: #7B1FA2;
            }
        """)
        self.save_complement_button.clicked.connect(lambda: self.save_sequence(results['complement'], "complement"))
        button_layout.addWidget(self.save_complement_button)
        
        # Add save reverse complement button
        self.save_reverse_complement_button = QPushButton("🔄 Save Reverse Complement")
        self.save_reverse_complement_button.setStyleSheet("""
            QPushButton {
                background-color: #FF9800;
                color: white;
                border: none;
                border-radius: 3px;
                padding: 10px;
                font-weight: bold;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: #F57C00;
            }
        """)
        self.save_reverse_complement_button.clicked.connect(lambda: self.save_sequence(results['reverse_complement'], "reverse_complement"))
        button_layout.addWidget(self.save_reverse_complement_button)
        
        # Add button container to results layout
        self.results_layout.addWidget(button_container)
        
        # Create result frames
        self.results_layout.addWidget(ResultFrame("🧬 Sequence Type", results['sequence_type']))
        self.results_layout.addWidget(ResultFrame("📄 Sequence Length", f"{format(results['length'], ',')} bases"))
        self.results_layout.addWidget(ResultFrame("🧪 GC Content", f"{results['gc_content']:.2f}%"))
        self.results_layout.addWidget(ResultFrame("🔬 AT Skew", f"{results['at_skew']:.4f}"))
        self.results_layout.addWidget(ResultFrame("🔬 GC Skew", f"{results['gc_skew']:.4f}"))

        # Nucleotide composition
        comp_text = "\n".join([f"{nt}: {format(count, ',')} ({pct:.2f}%)" 
                              for nt, (count, pct) in zip(results['nucleotide_counts'].keys(), 
                                                         zip(results['nucleotide_counts'].values(), 
                                                             results['nucleotide_percentages'].values()))])
        
        self.results_layout.addWidget(ResultFrame("🔢 Nucleotide Composition", comp_text))
        
        # Most common codons
        codon_text = "\n".join([f"{codon}: {format(count, ',')} times" 
                              for codon, count in results['codon_counts'].most_common(5)])
        self.results_layout.addWidget(ResultFrame("🔗 Most Common Codons", codon_text))
        
        # Missing genes (N count)
        if results['n_count'] > 0:
            self.results_layout.addWidget(ResultFrame("🔍 Missing Genes (N Count)", f"{format(results['n_count'], ',')} positions"))
            # Show filtered sequence previews only if there are N's
            filtered_preview_first = results['filtered_sequence'][:50]
            filtered_first_blocks = [filtered_preview_first[i:i+10] for i in range(0, len(filtered_preview_first), 10)]
            filtered_first_display = " ".join(filtered_first_blocks)
            self.results_layout.addWidget(ResultFrame("🧬 First 50 Bases (N's removed)", filtered_first_display))
            
            filtered_preview_last = results['filtered_sequence'][-50:]
            filtered_last_blocks = [filtered_preview_last[i:i+10] for i in range(0, len(filtered_preview_last), 10)]
            filtered_last_display = " ".join(filtered_last_blocks)
            self.results_layout.addWidget(ResultFrame("🧬 Last 50 Bases (N's removed)", filtered_last_display))

        # Original sequence display with 10-base deviations
        first_50 = results['sequence'][:50]
        first_blocks = [first_50[i:i+10] for i in range(0, len(first_50), 10)]
        first_display = " ".join(first_blocks)
        self.results_layout.addWidget(ResultFrame("🔍 First 50 Bases", first_display))
    
        last_50 = results['sequence'][-50:]
        last_blocks = [last_50[i:i+10] for i in range(0, len(last_50), 10)]
        last_display = " ".join(last_blocks)
        self.results_layout.addWidget(ResultFrame("🔍 Last 50 Bases", last_display))
    
    def generate_report(self):
        """Generate a formatted HTML report from the analysis results."""
        results = self.analysis_results
        
        # Format the values before inserting them into the template
        gc_skew = f"{results['gc_skew']:.4f}"
        at_skew = f"{results['at_skew']:.4f}"
        gc_content = f"{results['gc_content']:.2f}"
        n_count_percentage = f"{(results['n_count']/results['length']*100):.2f}"
        
        report = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{results['sequence_type']} Sequence Analysis Report</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            line-height: 1.6;
            color: #e0e0e0;
            background-color: #1e1e2f;
            max-width: 900px;
            margin: 0 auto;
            padding: 20px;
        }}
        h1, h2, h3 {{
            color: #42b0f5;
        }}
        h1 {{
            border-bottom: 2px solid #42b0f5;
            padding-bottom: 10px;
        }}
        h2 {{
            border-bottom: 1px solid #3d3d5f;
            padding-bottom: 5px;
            margin-top: 30px;
        }}
        h3 {{
            margin-top: 20px;
            color: #3a9ee0;
        }}
        .timestamp {{
            color: #888888;
            font-style: italic;
            margin-bottom: 30px;
        }}
        ul {{
            list-style-type: none;
            padding-left: 10px;
        }}
        li {{
            margin-bottom: 8px;
            padding-left: 20px;
            position: relative;
        }}
        li:before {{
            content: "•";
            color: #42b0f5;
            font-weight: bold;
            position: absolute;
            left: 0;
        }}
        .sequence-preview {{
            font-family: monospace;
            background-color: #2d2d3f;
            color: #e0e0e0;
            padding: 15px;
            border-radius: 5px;
            overflow-x: auto;
            white-space: pre-wrap;
            border: 1px solid #3d3d5f;
            margin-bottom: 15px;
            font-size: 16px;
            letter-spacing: 1px;
        }}
        .footer {{
            margin-top: 40px;
            padding-top: 10px;
            border-top: 1px solid #3d3d5f;
            color: #888888;
            font-style: italic;
        }}
        .key {{
            font-weight: bold;
            color: #42b0f5;
        }}
        .legend {{
            background-color: #2d2d3f;
            padding: 15px;
            border-radius: 5px;
            margin-top: 15px;
            border: 1px solid #3d3d5f;
        }}
        .legend h4 {{
            color: #42b0f5;
            margin-top: 0;
            margin-bottom: 10px;
        }}
        .legend-item {{
            display: inline-block;
            margin-right: 20px;
            margin-bottom: 5px;
        }}
        .legend-key {{
            font-weight: bold;
            color: #42b0f5;
        }}
    </style>
</head>
<body>
    <h1>🧬 {results['sequence_type']} Sequence Analysis Report</h1>
    <div class="timestamp">Generated on: {time.strftime('%Y-%m-%d %H:%M:%S')}</div>
    
    <h2>📊 Basic Statistics</h2>
    <ul>
        <li><span class="key">Sequence Type:</span> {results['sequence_type']}</li>
        <li><span class="key">Sequence Length:</span> {format(results['length'], ',')} bases</li>
        <li><span class="key">GC Content:</span> {gc_content}%</li>
        <li><span class="key">Missing Bases (N):</span> {format(results['n_count'], ',')} positions ({n_count_percentage}%)</li>
    </ul>
    
    <h2>🔬 Nucleotide Analysis</h2>
    <ul>"""
        
        # Nucleotide composition
        for nt, (count, pct) in zip(results['nucleotide_counts'].keys(), 
                                   zip(results['nucleotide_counts'].values(), 
                                       results['nucleotide_percentages'].values())):
            report += f'\n        <li><span class="key">{nt}:</span> {format(count, ",")} ({pct:.2f}%)</li>'
        
        # Add nucleotide legend
        report += """
    </ul>
    <div class="legend">
        <h4>Nucleotide Codes:</h4>
        <div class="legend-item"><span class="legend-key">A</span> = Adenine</div>
        <div class="legend-item"><span class="legend-key">C</span> = Cytosine</div>
        <div class="legend-item"><span class="legend-key">G</span> = Guanine</div>
        <div class="legend-item"><span class="legend-key">T</span> = Thymine</div>
        <div class="legend-item"><span class="legend-key">U</span> = Uracil (RNA only)</div>
        <div class="legend-item"><span class="legend-key">N</span> = Any nucleotide (A/C/G/T)</div>
        <div class="legend-item"><span class="legend-key">R</span> = Purine (A/G)</div>
        <div class="legend-item"><span class="legend-key">Y</span> = Pyrimidine (C/T)</div>
        <div class="legend-item"><span class="legend-key">S</span> = Strong bonds (G/C)</div>
        <div class="legend-item"><span class="legend-key">W</span> = Weak bonds (A/T)</div>
        <div class="legend-item"><span class="legend-key">K</span> = Keto (G/T)</div>
        <div class="legend-item"><span class="legend-key">M</span> = Amino (A/C)</div>
        <div class="legend-item"><span class="legend-key">B</span> = Not A (C/G/T)</div>
        <div class="legend-item"><span class="legend-key">D</span> = Not C (A/G/T)</div>
        <div class="legend-item"><span class="legend-key">H</span> = Not G (A/C/T)</div>
        <div class="legend-item"><span class="legend-key">V</span> = Not T (A/C/G)</div>
    </div>
    
    <h2>⚖️ Sequence Bias Analysis</h2>
    <ul>
        <li><span class="key">AT Skew:</span> """ + at_skew + """</li>
        <li><span class="key">GC Skew:</span> """ + gc_skew + """</li>
    </ul>"""
        
        # Codon analysis with formatted numbers
        report += "\n    <h2>🧬 Codon Analysis</h2>\n    <ul>"
        for codon, count in results['codon_counts'].most_common(5):
            report += f'\n        <li><span class="key">{codon}:</span> {format(count, ",")} times</li>'
        report += "\n    </ul>"
        
        # Sequence previews with 10-base deviations
        report += """
    <h2>🔍 Sequence Preview</h2>
    <div class="preview-section">
        <h3>First 50 Bases (N's removed)</h3>
        <div class="sequence-preview">"""
        
        # First 50 filtered bases with 10-base deviations
        filtered_preview_first = results['filtered_sequence'][:50]
        filtered_first_blocks = [filtered_preview_first[i:i+10] for i in range(0, len(filtered_preview_first), 10)]
        report += " ".join(filtered_first_blocks)
        
        report += """</div>
        <h3>Last 50 Bases (N's removed)</h3>
        <div class="sequence-preview">"""
        
        # Last 50 filtered bases with 10-base deviations
        filtered_preview_last = results['filtered_sequence'][-50:]
        filtered_last_blocks = [filtered_preview_last[i:i+10] for i in range(0, len(filtered_preview_last), 10)]
        report += " ".join(filtered_last_blocks)
        
        report += """</div>
        <h3>Original Sequence (First 50)</h3>
        <div class="sequence-preview">"""
        
        # First 50 original bases with 10-base deviations
        first_50 = results['sequence'][:50]
        first_blocks = [first_50[i:i+10] for i in range(0, len(first_50), 10)]
        report += " ".join(first_blocks)
        
        report += """</div>
        <h3>Original Sequence (Last 50)</h3>
        <div class="sequence-preview">"""
        
        # Last 50 original bases with 10-base deviations
        last_50 = results['sequence'][-50:]
        last_blocks = [last_50[i:i+10] for i in range(0, len(last_50), 10)]
        report += " ".join(last_blocks)
        
        report += """</div>
    </div>
    
    <div class="footer">
        <p>Generated by Genetic Sequence Analyzer</p>
        <p>Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
</body>
</html>"""
        
        return report
    
    def save_report(self):
        """Save the analysis report to a file."""
        if not hasattr(self, 'analysis_results'):
            return
        
        # Get the script directory to use as default save location
        import os
        script_dir = os.path.dirname(os.path.abspath(__file__))
        
        # Create default filename with current date and time
        default_filename = f"DNA_Analysis_{time.strftime('%Y-%m-%d_%H-%M-%S')}.html"
        default_path = os.path.join(script_dir, default_filename)
            
        file_name, _ = QFileDialog.getSaveFileName(
            self,
            "Save Analysis Report",
            default_path,
            "HTML Files (*.html);;All Files (*)"
        )
        
        if file_name:
            report = self.generate_report()
            try:
                with open(file_name, 'w', encoding='utf-8') as f:
                    f.write(report)
                
                # Show success message without copy button
                success_frame = QFrame()
                success_frame.setStyleSheet("""
                    QFrame {
                        background-color: #4CAF50;
                        border-radius: 5px;
                        padding: 15px;
                        margin: 5px;
                    }
                """)
                success_layout = QVBoxLayout(success_frame)
                
                success_label = QLabel(f"✅ Report successfully saved to {file_name}")
                success_label.setStyleSheet("""
                    color: white;
                    font-weight: bold;
                    font-size: 16px;
                """)
                success_label.setWordWrap(True)
                
                success_layout.addWidget(success_label)
                self.results_layout.insertWidget(1, success_frame)  # Insert after save button
                
                # Auto-remove success message after 5 seconds
                QTimer.singleShot(5000, lambda: success_frame.setParent(None))
                
            except Exception as e:
                # Show error message without copy button
                error_frame = QFrame()
                error_frame.setStyleSheet("""
                    QFrame {
                        background-color: #ff4444;
                        border-radius: 5px;
                        padding: 15px;
                        margin: 5px;
                    }
                """)
                error_layout = QVBoxLayout(error_frame)
                
                error_label = QLabel(f"❌ Error Saving Report: {str(e)}")
                error_label.setStyleSheet("""
                    color: white;
                    font-weight: bold;
                    font-size: 16px;
                """)
                error_label.setWordWrap(True)
                
                error_layout.addWidget(error_label)
                self.results_layout.insertWidget(1, error_frame)  # Insert after save button
                
                # Auto-remove error message after 5 seconds
                QTimer.singleShot(5000, lambda: error_frame.setParent(None))
    
    def reset_analysis(self):
        # Clear results
        for i in reversed(range(self.results_layout.count())): 
            self.results_layout.itemAt(i).widget().setParent(None)
        
        # Re-enable all inputs
        self.sequence_input.setEnabled(True)
        self.browse_button.setEnabled(True)
        self.reset_file_button.setEnabled(True)
        self.analyze_button.setEnabled(True)
        
        # Clear input and file selection
        self.sequence_input.clear()
        self.current_file = None
        self.file_path_label.setText("No file selected")
        self.file_path_label.setStyleSheet("""
            color: #888888;
            font-size: 14px;
        """)
        self.browse_button.setText("Browse FASTA File")
        self.reset_file_button.setVisible(False)
        
        # Reset button to analyze state
        self.analyze_button.setText("Analyze Sequence")
        self.analyze_button.setStyleSheet("""
            QPushButton {
                background-color: #42b0f5;
                color: white;
                border: none;
                border-radius: 3px;
                padding: 10px;
                font-weight: bold;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: #3a9ee0;
            }
        """)
        self.analyze_button.clicked.disconnect()
        self.analyze_button.clicked.connect(self.analyze_sequence)

    def save_sequence(self, sequence, sequence_type):
        """Save a sequence as a FASTA file."""
        if not sequence:
            return
        
        # Get the script directory to use as default save location
        import os
        script_dir = os.path.dirname(os.path.abspath(__file__))
        
        # Create default filename with current date and time
        default_filename = f"{sequence_type}_{time.strftime('%Y-%m-%d_%H-%M-%S')}.fasta"
        default_path = os.path.join(script_dir, default_filename)
            
        file_name, _ = QFileDialog.getSaveFileName(
            self,
            f"Save {sequence_type.replace('_', ' ').title()}",
            default_path,
            "FASTA Files (*.fasta);;All Files (*)"
        )
        
        if file_name:
            try:
                with open(file_name, 'w', encoding='utf-8') as f:
                    f.write(f">{sequence_type}\n")
                    # Write sequence in chunks of 80 characters
                    for i in range(0, len(sequence), 80):
                        f.write(sequence[i:i+80] + "\n")
                
                # Show success message
                success_frame = QFrame()
                success_frame.setStyleSheet("""
                    QFrame {
                        background-color: #4CAF50;
                        border-radius: 5px;
                        padding: 15px;
                        margin: 5px;
                    }
                """)
                success_layout = QVBoxLayout(success_frame)
                
                success_label = QLabel(f"✅ {sequence_type.replace('_', ' ').title()} successfully saved to {file_name}")
                success_label.setStyleSheet("""
                    color: white;
                    font-weight: bold;
                    font-size: 16px;
                """)
                success_label.setWordWrap(True)
                
                success_layout.addWidget(success_label)
                self.results_layout.insertWidget(1, success_frame)  # Insert after button container
                
                # Auto-remove success message after 5 seconds
                QTimer.singleShot(5000, lambda: success_frame.setParent(None))
                
            except Exception as e:
                # Show error message
                error_frame = QFrame()
                error_frame.setStyleSheet("""
                    QFrame {
                        background-color: #ff4444;
                        border-radius: 5px;
                        padding: 15px;
                        margin: 5px;
                    }
                """)
                error_layout = QVBoxLayout(error_frame)
                
                error_label = QLabel(f"❌ Error Saving {sequence_type.replace('_', ' ').title()}: {str(e)}")
                error_label.setStyleSheet("""
                    color: white;
                    font-weight: bold;
                    font-size: 16px;
                """)
                error_label.setWordWrap(True)
                
                error_layout.addWidget(error_label)
                self.results_layout.insertWidget(1, error_frame)  # Insert after button container
                
                # Auto-remove error message after 5 seconds
                QTimer.singleShot(5000, lambda: error_frame.setParent(None))

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = GeneticAnalyzerGUI()
    window.show()
    sys.exit(app.exec_())