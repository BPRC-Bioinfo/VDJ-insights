from flask import Flask, send_from_directory, redirect
from pathlib import Path
import webbrowser
import threading

app = Flask(__name__)

# Set the base directory to point to the parent directory of 'scripts', where 'source' is located
base_dir = Path(__file__).resolve().parent.parent

# Automatically redirect from the root URL ("/") to "/html/index.html"


@app.route('/')
def redirect_to_index():
    return redirect('/html/index.html')

# Serve any HTML file from the source/html folder


@app.route('/html/<path:filename>')
def serve_html(filename):
    return send_from_directory(base_dir / 'source' / 'html', filename)

# Serve CSS files from the source/css folder


@app.route('/css/<path:filename>')
def serve_css(filename):
    return send_from_directory(base_dir / 'source' / 'css', filename)

# Serve JS files from the source/js folder


@app.route('/js/<path:filename>')
def serve_js(filename):
    return send_from_directory(base_dir / 'source' / 'js', filename)

# Serve image files from the source/images folder


@app.route('/images/<path:filename>')
def serve_images(filename):
    return send_from_directory(base_dir / 'source' / 'images', filename)


def open_browser():
    # Open the default web browser to the given URL after the Flask app starts
    webbrowser.open_new('http://127.0.0.1:8000/html/index.html')


if __name__ == '__main__':
    # Start a thread to open the browser automatically
    threading.Timer(1, open_browser).start()

    # Run the Flask app and listen on 127.0.0.1:8000
    app.run(debug=True, port=8000)
