from flask import Flask, send_from_directory, redirect
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


def open_browser():
    # Open the default web browser to the given URL after the Flask app starts
    webbrowser.open_new('http://127.0.0.1:8000/html/index.html')


if __name__ == '__main__':
    pass
