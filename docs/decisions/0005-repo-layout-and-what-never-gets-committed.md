# 0005: Use a src-layout package, keep data out of git, and keep candidate coordinates private

**Decision:** Code lives in `src/minefinder/` and gets installed in editable mode. Everything under `data/` and `outputs/` is gitignored. The repo stays **private** until there's a cleaned public version.

## In plain words

- **`src/` layout**: the package sits one folder down, so Python can only import it if it's properly installed (`pip install -e .`). That sounds like a hassle, but it catches the "works on my machine because I happened to run from the repo folder" bug early. Editable mode (`-e`) means your edits take effect immediately, with no reinstall.
- **No data in git**: LiDAR tiles are gigabytes, and git is built for text. The raw inventories are public and re-downloadable, so the *download scripts* are what we commit, not the files.
- **No candidate coordinates in git**: our outputs are locations of *possibly unrecorded mines*. That's a safety issue (open shafts), a land-access issue, and a looting issue (artifacts on federal land are protected). They stay in `outputs/`, which never gets committed. A private repo is a second layer of protection on top of that.

## Formally

- PyPA's **src layout** prevents the implicit `sys.path` import of the working directory. `pytest` is configured with `pythonpath = ["src"]`, and the package is installed with setuptools' `packages.find(where="src")`.
- `.gitignore` excludes `data/*/*` but keeps `.gitkeep` placeholders, because git doesn't track empty folders.
- Project non-negotiable: **never include exact coordinates of candidate sites in anything meant for publication.** For a public portfolio version later, publish metrics, aggregate maps (e.g. hex bins at 1 km or coarser) and code, never point lists.

## Where it lives

`pyproject.toml`, `.gitignore`, the `data/` and `outputs/` folders.

## Revisit when

You want the repo public. Then make a separate clean repo or orphan branch, check that the git history contains no coordinates, and add a `SECURITY.md`-style note on responsible use.
