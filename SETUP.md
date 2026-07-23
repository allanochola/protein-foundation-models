# SETUP

This document describes the complete workflow used in this repository.

The workflow is intentionally simple.

```
GitHub
    ↓
Google Colab
    ↓
Run experiment
    ↓
Review results
    ↓
Update documentation
    ↓
Commit
    ↓
Push to GitHub
```

GitHub is the permanent record.

Google Colab is temporary compute.

---

# One-Time Setup

These steps only need to be completed once.

## 1. Create the GitHub repository

Create a new public repository.

Do **not** allow GitHub to create:

- README
- .gitignore
- LICENSE

The repository already contains them.

Upload the repository contents.

---

## 2. Create a GitHub Personal Access Token

Navigate to

```
Settings
Developer Settings
Personal Access Tokens
Fine-grained Tokens
```

Generate a new token.

Repository access

```
Only select repositories
```

Permissions

```
Contents

Read and write
```

Nothing else is required.

Copy the token immediately.

GitHub only displays it once.

---

## 3. Add the token to Colab

Open Google Colab.

Click the key icon.

Create

```
GH_TOKEN
```

Paste the token.

Enable

```
Notebook access
```

Without this toggle

```python
userdata.get("GH_TOKEN")
```

will fail.

---

## 4. Install Jupytext

Local machine

```bash
pip install jupytext
```

Jupytext is used to generate notebooks from the canonical Python scripts.

The repository tracks

```
.py
```

not

```
.ipynb
```

Generated notebooks are disposable.

---

# Every Colab Session

Always follow these steps in order.

---

## Step 1

Create a new Colab notebook.

For Experiment 01

```
Runtime

↓

Change runtime type

↓

CPU
```

GPU is unnecessary.

---

## Step 2

Load the GitHub token.

```python
from google.colab import userdata

token = userdata.get("GH_TOKEN")

assert token is not None
```

---

## Step 3

Configure Git credentials.

```python
from pathlib import Path

credentials = Path("/root/.git-credentials")

credentials.write_text(
    f"https://x-access-token:{token}@github.com\n"
)

credentials.chmod(0o600)
```

Enable the credential helper.

```python
!git config --global credential.helper store
```

Configure Git identity.

```python
!git config --global user.name "Your Name"

!git config --global user.email "your@email.com"
```

---

## Step 4

Clone the repository.

```bash
git clone https://github.com/USERNAME/protein-foundation-models.git
```

Move into the repository.

```bash
cd protein-foundation-models
```

---

## Step 5

Always pull first.

```bash
git pull --rebase origin main
```

Never begin work on an outdated copy.

---

## Step 6

Install dependencies.

```bash
pip install -r requirements.txt
```

Dependencies belong in

```
requirements.txt
```

Never install packages inside experiment scripts.

---

## Step 7

Run the experiment.

Example

```bash
python notebooks/01a_msa_depth_published_scores.py
```

Run scripts from the repository root.

Scripts resolve output paths relative to themselves.

---

## Step 8

Review outputs.

Typical outputs are

```
results/

figures/

experiments/
```

Inspect

- CSV files
- figures
- provenance
- logs

before committing.

---

## Step 9

Review Git status.

```bash
git status
```

Then

```bash
git diff --stat
```

Stage only the files belonging to the experiment.

Example

```bash
git add \
notebooks/01a_msa_depth_published_scores.py \
results/esm2_ladder_by_depth.csv \
figures/esm2_scaling_by_depth.png \
models/esm2.md
```

Avoid

```bash
git add .
```

Large experiments often generate files that should never enter version control.

---

## Step 10

Commit.

Commit messages should describe findings.

Good

```
01a: interaction negative above 650M
```

Bad

```
update

changes

new version
```

---

## Step 11

Push.

```bash
git push origin main
```

Verify GitHub reflects the commit.

---

## Step 12

Delete credentials.

```python
credentials.unlink(missing_ok=True)
```

Colab machines are temporary.

Credentials should be too.

---

# Long Experiments

Colab sessions terminate unexpectedly.

Every experiment that runs for more than approximately thirty minutes should checkpoint.

The rule is

```
checkpoint

↓

commit

↓

push
```

A checkpoint committed only inside the Colab VM disappears when the runtime dies.

---

# Reproducibility Rules

Every experiment must

- begin with a protocol
- record software versions
- record input dataset version
- pin external repositories by commit SHA
- produce deterministic output paths
- save provenance
- update documentation if conclusions change

---

# Repository Philosophy

Documentation is not separate from experiments.

Every experiment may modify

```
models/

models.csv

README
```

if the evidence changes.

Claims are living objects.

Experiments determine whether they survive.

---

# Canonical Workflow

```
Write protocol

↓

Run experiment

↓

Review outputs

↓

Interpret

↓

Update documentation

↓

Commit

↓

Push
```

The protocol always comes before the experiment.

Documentation always comes after it.
