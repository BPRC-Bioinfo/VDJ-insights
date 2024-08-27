from flask import Flask, send_from_directory, redirect, request, jsonify
from pathlib import Path
import webbrowser
import threading
import sys

app = Flask(__name__)

# Global variable for the base directory
base_dir = None

# Function to dynamically set base_dir


def set_base_dir(input_dir):
    global base_dir
    base_dir = Path(input_dir).resolve()

# Automatically redirect from the root URL ("/") to "/html/index.html"


@app.route('/')
def redirect_to_index():
    return redirect('/html/index.html')

# Serve any HTML file from the dynamically set base directory


@app.route('/html/<path:filename>')
def serve_html(filename):
    if base_dir is None:
        return "Base directory not set", 500
    return send_from_directory(base_dir / 'html', filename)

# Serve CSS files from the dynamically set base directory


@app.route('/css/<path:filename>')
def serve_css(filename):
    if base_dir is None:
        return "Base directory not set", 500
    return send_from_directory(base_dir / 'css', filename)

# Serve JS files from the dynamically set base directory


@app.route('/js/<path:filename>')
def serve_js(filename):
    if base_dir is None:
        return "Base directory not set", 500
    return send_from_directory(base_dir / 'js', filename)

# Serve image files from the dynamically set base directory


@app.route('/images/<path:filename>')
def serve_images(filename):
    if base_dir is None:
        return "Base directory not set", 500
    return send_from_directory(base_dir / 'images', filename)

# New endpoint to handle sequence saving


@app.route('/save_sequences', methods=['POST'])
def save_sequences():
    new_sequences = request.json.get('sequences', {})
    fasta_file_path = base_dir.parent / ".tool/library/library.fasta"

    if not new_sequences:
        return jsonify({"error": "No sequences provided"}), 400

    total_new = append_to_fasta(fasta_file_path, new_sequences)
    if total_new > 0:
        message = f"Added a total of {total_new} new sequences to the library"
    else:
        message = "No new sequences found to add"
    return jsonify({"message": message}), 200

# Read existing FASTA file and return a dictionary of {header: sequence}


def read_fasta(file_path):
    sequences = {}
    header = None
    sequence = []

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    sequences[header] = ''.join(sequence)
                header = line
                sequence = []
            else:
                sequence.append(line)
        if header:
            sequences[header] = ''.join(sequence)  # Add the last sequence

    return sequences

# Append new sequences to FASTA file only if not already present


def append_to_fasta(file_path, new_sequences):
    existing_sequences = read_fasta(file_path)
    total_new = 0
    with open(file_path, 'a') as file:
        for header, sequence in new_sequences.items():
            if header not in existing_sequences and sequence not in existing_sequences.values():
                total_new += 1
                file.write(f"{header}\n{sequence}\n")
    return total_new


def open_browser():
    # Open the default web browser to the given URL after the Flask app starts
    webbrowser.open_new('http://127.0.0.1:8000/html/index.html')


if __name__ == '__main__':
    # Parse the base directory from command-line arguments
    if len(sys.argv) > 1:
        set_base_dir(sys.argv[1])
    else:
        print("Error: No base directory provided")
        sys.exit(1)

    # Open the browser only after the first request, ensuring the server is ready
    threading.Timer(1, open_browser).start()

    # Start the Flask app
    app.run(debug=True, host='0.0.0.0', port=8000)
