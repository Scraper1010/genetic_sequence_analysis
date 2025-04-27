# 🧬 Genetic Sequence Analyzer 🌟

A modern, user-friendly desktop application for analyzing DNA and RNA sequences. Built with PyQt5, this tool provides comprehensive sequence analysis with a beautiful dark-themed interface. ✨

## ✨ Features

- **Sequence Analysis** 🔍
  - DNA/RNA sequence type detection 🧪
  - GC content calculation 📊
  - AT/GC skew analysis 📈
  - Nucleotide composition analysis 🧬
  - Codon frequency analysis 🔢
  - Missing gene (N) detection and filtering 🔎

- **Sequence Transformations** 🔄
  - RNA transcription 🧬
  - DNA/RNA complement generation 🔁
  - Reverse complement generation 🔂

- **File Operations** 💾
  - Load sequences from FASTA files 📂
  - Save analysis reports in HTML format 📄
  - Export transformed sequences in FASTA format 📤

- **User Interface** 🎨
  - Modern dark theme 🌙
  - Real-time progress tracking 📊
  - Interactive sequence previews 👀
  - Copy-to-clipboard functionality 📋
  - Responsive layout with split views 📱

## 📋 Requirements

- Python 3.6 or higher 🐍
- PyQt5 🖥️
- Biopython 🧬

## 🚀 Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/genetic-sequence-analyzer.git
cd genetic-sequence-analyzer
```

2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

## 💻 Usage

1. Run the application:
```bash
python main.py
```

2. Input your sequence:
   - Type or paste a DNA/RNA sequence directly into the text box ✍️
   - OR load a sequence from a FASTA file using the "Browse FASTA File" button 📂

3. Click "Analyze Sequence" to process the sequence 🔍

4. View the analysis results:
   - Basic statistics (length, GC content, etc.) 📊
   - Nucleotide composition 🧬
   - Sequence bias analysis 📈
   - Codon analysis 🔢
   - Sequence previews 👀

5. Save your results:
   - Generate and save an HTML report 📄
   - Export transformed sequences (RNA, complement, reverse complement) 📤

## 📊 Analysis Features

### Basic Statistics 📈
- Sequence type detection (DNA/RNA) 🧬
- Sequence length 📏
- GC content percentage 📊
- Missing bases (N) count and percentage 🔍

### Nucleotide Analysis 🧬
- Detailed nucleotide composition 📝
- Percentage distribution 📊
- Comprehensive nucleotide code legend 📚

### Sequence Bias 📈
- AT skew calculation 📊
- GC skew calculation 📊

### Codon Analysis 🔢
- Top 5 most frequent codons 🏆
- Codon frequency distribution 📊

### Sequence Previews 👀
- First and last 50 bases 🔍
- Filtered sequences (N's removed) 🧹
- 10-base block formatting for readability 📝

## 💾 File Formats

### Input 📥
- FASTA files (*.fasta, *.fa) 📂
- Raw sequence text ✍️

### Output 📤
- HTML reports (*.html) 📄
- FASTA files (*.fasta) 📂

## 🎨 Interface Features

- **Dark Theme** 🌙
  - Modern, eye-friendly color scheme 🎨
  - High contrast for better readability 👀
  - Consistent styling across all elements ✨

- **Interactive Elements** 🖱️
  - Progress bar for analysis status 📊
  - Copy buttons for quick sequence copying 📋
  - Auto-dismissing success/error messages 💫
  - Responsive button states 🔄

## 🔧 Technical Details

- Built with PyQt5 for the GUI 🖥️
- Uses Biopython for sequence handling 🧬
- Multi-threaded analysis for better performance ⚡
- HTML report generation with CSS styling 🎨
- FASTA file format support 📂

## 📝 License

This project is licensed under the MIT License - see the LICENSE file for details. 📜

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request. 🌟

## ⚠️ Known Limitations

- Maximum sequence length is limited by system memory 💾
- Some complex analyses may take longer for very large sequences ⏳
- FASTA files must be properly formatted 📝

## 📞 Support

For support, please open an issue in the GitHub repository or contact the maintainers. 💌 