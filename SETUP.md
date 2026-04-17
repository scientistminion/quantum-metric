# GitHub & Docs Setup Guide

Follow these steps to get your code on GitHub, set up automatic testing, and publish docs to GitHub Pages.

## 1. Drop these files into your existing `quantum-metric/` directory

After unzipping, you should have this added to your project:

```
quantum-metric/
├── .github/
│   └── workflows/
│       ├── docs.yml         ← builds & deploys docs on every push to main
│       └── tests.yml        ← runs pytest on py3.9–3.12
└── docs/
    ├── Makefile
    ├── requirements.txt
    └── source/
        ├── conf.py
        ├── index.md
        ├── installation.md
        ├── quickstart.md
        ├── cli.md
        ├── python_api.md
        ├── theory.md
        ├── changelog.md
        ├── api/index.md
        └── examples/
            ├── index.md
            └── ag_fcc.md
```

## 2. Edit docs to use your real GitHub username

In `docs/source/conf.py`, find and replace:

- `"yourusername"` → your actual GitHub username
- `"Your Name"` → your real name

Also in `docs/source/index.md`:

- Replace every `yourusername` with your actual GitHub username

## 3. Build docs locally (test before pushing)

```bash
cd docs
pip install -r requirements.txt
make html

# Open build/html/index.html in your browser
firefox build/html/index.html   # or `xdg-open` on Linux, `open` on Mac
```

Fix any warnings, then push.

## 4. Create the GitHub repo and push

```bash
cd quantum-metric/       # the project root

git init
git add .
git commit -m "Initial commit"

# Create the empty repo on GitHub first, then:
git remote add origin https://github.com/YOUR_USERNAME/quantum-metric.git
git branch -M main
git push -u origin main
```

## 5. Enable GitHub Pages

1. Go to your repo on GitHub
2. Click **Settings** → **Pages**
3. Under **Build and deployment**, set **Source** to **GitHub Actions**
4. Save

## 6. Trigger the docs build

Either push another commit or go to **Actions** → **Build and deploy docs** → **Run workflow**.

After ~2 minutes, your docs will be live at:

```
https://YOUR_USERNAME.github.io/quantum-metric/
```

Add a link to this URL in your project's README.md so visitors find the docs.

## 7. (Optional) Protect the main branch

In **Settings** → **Branches** → **Add rule**:

- Branch name pattern: `main`
- Require a pull request before merging
- Require status checks to pass (`Tests` and `Build and deploy docs`)

This prevents you from accidentally pushing broken code to main.

## Local docs development loop

While editing docs:

```bash
cd docs
make clean && make html
# Then reload build/html/index.html
```

Or for live reload (optional):

```bash
pip install sphinx-autobuild
sphinx-autobuild source build/html
# Now edit any .md file; the browser auto-reloads
```

## Adding a new example

1. Create `docs/source/examples/my_material.md`
2. Add it to the toctree in `docs/source/examples/index.md`
3. Rebuild

## Publishing to PyPI later

Once your docs are live and you're ready to release to PyPI:

```bash
pip install build twine
python -m build
twine upload --repository testpypi dist/*    # test first
twine upload dist/*                           # real
```

After that, the PyPI badge in your docs will light up green and anyone can `pip install quantum-metric`.
