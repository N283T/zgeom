#!/usr/bin/env python3
"""Install a SHA-256-pinned official Zig archive on a GitHub runner."""

import argparse
import hashlib
import os
from pathlib import Path
import shutil
import tarfile
import tempfile
import urllib.request
import zipfile


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--directory", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    runner_temp = Path(os.environ["RUNNER_TEMP"])
    install_root = runner_temp / "zig-install"
    if install_root.exists():
        shutil.rmtree(install_root)
    install_root.mkdir()

    suffix = ".zip" if args.url.endswith(".zip") else ".tar.xz"
    with tempfile.NamedTemporaryFile(dir=runner_temp, suffix=suffix, delete=False) as output:
        archive_path = Path(output.name)
        digest = hashlib.sha256()
        request = urllib.request.Request(args.url, headers={"User-Agent": "zgeom-ci"})
        with urllib.request.urlopen(request, timeout=120) as response:
            while chunk := response.read(1024 * 1024):
                output.write(chunk)
                digest.update(chunk)

    actual = digest.hexdigest()
    expected = args.sha256.lower()
    if actual != expected:
        archive_path.unlink(missing_ok=True)
        raise SystemExit(f"SHA-256 mismatch: expected {expected}, got {actual}")

    if suffix == ".zip":
        with zipfile.ZipFile(archive_path) as archive:
            archive.extractall(install_root)
    else:
        with tarfile.open(archive_path, "r:xz") as archive:
            archive.extractall(install_root, filter="data")
    archive_path.unlink()

    zig_dir = install_root / args.directory
    executable = zig_dir / ("zig.exe" if os.name == "nt" else "zig")
    if not executable.is_file():
        raise SystemExit(f"Zig executable not found after extraction: {executable}")

    github_path = os.environ.get("GITHUB_PATH")
    if not github_path:
        raise SystemExit("GITHUB_PATH is not set")
    with Path(github_path).open("a", encoding="utf-8") as output:
        output.write(f"{zig_dir}{os.linesep}")
    print(f"Installed and verified {args.url} ({actual})")


if __name__ == "__main__":
    main()
