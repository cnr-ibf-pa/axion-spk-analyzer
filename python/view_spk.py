"""
view_spk.py - open an Axion .spk file in the HTML viewer.

    python view_spk.py path/to/recording.spk [--port 8765] [--no-browser]

Starts a small local web server (127.0.0.1 only) that serves spk_viewer.html
plus the requested .spk file, then opens the browser on it.
You can also just open spk_viewer.html directly and drag-and-drop a file.
"""
import argparse
import http.server
import os
import urllib.parse
import webbrowser

APP_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # folder with the HTML pages


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("spk", nargs="?", help=".spk file to open (optional)")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--no-browser", action="store_true")
    a = ap.parse_args()

    spk_path = os.path.abspath(a.spk) if a.spk else None
    if spk_path and not os.path.isfile(spk_path):
        ap.error(f"file not found: {spk_path}")
    spk_name = os.path.basename(spk_path) if spk_path else None

    class Handler(http.server.SimpleHTTPRequestHandler):
        def __init__(self, *args, **kw):
            super().__init__(*args, directory=APP_DIR, **kw)

        def translate_path(self, path):
            p = urllib.parse.unquote(urllib.parse.urlparse(path).path)
            if spk_path and p == "/data/" + spk_name:
                return spk_path
            if p in ("/", ""):
                return os.path.join(APP_DIR, "spk_viewer.html")
            return super().translate_path(path)

        def log_message(self, *args):
            pass

    url = f"http://127.0.0.1:{a.port}/spk_viewer.html"
    if spk_path:
        url += "?src=" + urllib.parse.quote("/data/" + spk_name)
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", a.port), Handler)
    print(f"Serving viewer at {url}  (Ctrl+C to stop)")
    if not a.no_browser:
        webbrowser.open(url)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
