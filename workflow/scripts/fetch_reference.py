"""Fetch/copy a reference FASTA or GTF into place for the index rules.

Used by the prepare_genome / prepare_gtf rules (see rules/index.smk). The
source -- taken from snakemake.params.src -- is either a local path (plain, or
gzip/bzip2-compressed; the compression is sniffed from the magic bytes, not
assumed from the file extension) or an http(s)/ftp(s)/file:// URL. The result
is always written uncompressed to snakemake.output[0].

Runs under snakemake's `script:` directive, so `snakemake` is injected; the
source may be a URL, which is why the rules deliberately declare no `input`.
"""
import bz2
import gzip
import os
import shutil
import sys
import urllib.error
import urllib.request
from urllib.parse import urlparse

BLOCK = 1024 * 1024


def _log(*parts):
    print(" ".join(str(p) for p in parts), file=sys.stderr, flush=True)


def _is_url(path):
    return urlparse(str(path)).scheme in ("http", "https", "ftp", "file")


def _fail(msg):
    print(msg, file=sys.stderr, flush=True)
    sys.exit(1)


def _sanity_check(path, src):
    """Cheap sanity checks so a wrong/broken source fails here with a clear
    message instead of deep inside a tool (e.g. STAR) hours later."""
    size = os.path.getsize(path)
    if size == 0:
        _fail(f"fetch_reference: {src!r} produced an empty file at {path}")
    with open(path) as fh:
        first = fh.read(1)
    if path.endswith(".fa") and first != ">":
        _fail(
            f"fetch_reference: {path} does not look like a FASTA (first byte "
            f"is {first!r}, expected '>'). Check the --fasta source {src!r}."
        )


def main():
    # snakemake's `script:` directive does not capture the process's output
    # into the rule log on its own -- redirect stderr to it so progress and
    # any failure messages are recorded where the user is told to look.
    sys.stderr = open(snakemake.log[0], "w")  # noqa: F821
    src = str(snakemake.params.src)  # noqa: F821
    dest = str(snakemake.output[0])  # noqa: F821
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)

    tmp = dest + ".part"
    if os.path.exists(tmp):
        os.remove(tmp)

    try:
        if _is_url(src):
            _log(f"fetch_reference: downloading {src}")
            try:
                with urllib.request.urlopen(src) as resp, open(tmp, "wb") as out:
                    shutil.copyfileobj(resp, out, BLOCK)
            except urllib.error.URLError as exc:
                _fail(f"fetch_reference: failed to download {src}: {exc}")
            _log(f"fetch_reference: downloaded to {tmp}")
        else:
            if not os.path.exists(src):
                _fail(f"fetch_reference: source not found: {src!r}")
            _log(f"fetch_reference: copying {src}")
            shutil.copyfile(src, tmp)

        with open(tmp, "rb") as fh:
            magic = fh.read(4)

        if magic[:2] == b"\x1f\x8b":
            opener, mode = gzip.open, "gzip"
        elif magic[:3] == b"BZh":
            opener, mode = bz2.open, "bzip2"
        else:
            opener, mode = None, None

        if opener is not None:
            _log(f"fetch_reference: decompressing ({mode}) -> {dest}")
            with opener(tmp, "rb") as fin, open(dest, "wb") as fout:
                shutil.copyfileobj(fin, fout, BLOCK)
            os.remove(tmp)
        else:
            os.replace(tmp, dest)

        _sanity_check(dest, src)
        _log(f"fetch_reference: wrote {dest}")
    except OSError as exc:
        _fail(f"fetch_reference: failed writing {dest}: {exc}")
    finally:
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass


if __name__ == "__main__":
    main()
